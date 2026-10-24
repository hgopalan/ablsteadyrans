How-to guides
=============

Short recipes for common tasks. Each assumes you have run :doc:`tutorial_column` and know how an input file looks.
The full list of inputs is in :doc:`inputs`.

.. contents:: On this page
   :local:
   :depth: 1

Change an input from the command line
-------------------------------------

Add ``key=value`` after the input file:

.. code-block:: bash

   mpirun -n 4 build-rel/gust inputs steady.max_iterations=500 io.plot_file=test_plt

Use this for keys that the input file does not set: a key given twice is refused. To change a key the file sets,
edit the file (or a copy of it).

Restart a run
-------------

Every run writes a checkpoint at the end, ``<io.chk_file>NNNNN``. To continue from it, for example with more
iterations:

#. Copy the input file. In the copy, remove every ``init.*`` key (a restart takes its state from the checkpoint)
   and raise ``steady.max_iterations`` if needed.
#. Give the checkpoint directory on the command line:

   .. code-block:: bash

      mpirun -n 4 build-rel/gust inputs_restart io.restart=steady_chk01000

- The mesh and terrain inputs and the floating-point precision must match the checkpoint. The number of ranks and
  ``mesh.max_grid_size`` may change.
- ``steady.max_iterations`` counts the iterations of the whole run, restarts included.
- The run stops if it would overwrite the checkpoint it read. Change ``io.chk_file`` if needed.
- A restarted run reproduces the uninterrupted run bit for bit on the same ranks and box layout, over flat ground and
  over terrain. A checkpoint written before audit A2 (October 2026) restarts too; over terrain the continued run then
  follows the current wall model (:doc:`boundary_conditions`), so it is not the run the checkpoint came from.

Use your own terrain
--------------------

**From a point cloud.** Write the terrain as ``x y z`` lines (metres, one point per line) and set

.. code-block:: text

   terrain.type = xyz
   terrain.file = my_terrain.xyz
   terrain.idw_neighbors = 6     # points used for each mesh node
   terrain.idw_power     = 2     # inverse-distance weight exponent

Each mesh node gets the inverse-distance-weighted height of its nearest points. The points should cover the whole
domain, including a margin beyond its edges.

**From an ERF terrain file.** Set ``terrain.type = erf_grid`` and ``terrain.file`` to ERF's structured terrain file,
as in :doc:`tutorial_askervein`.

**Check the mesh first.** Run with ``steady.stage = none`` (the default): GUST builds the mesh, prints its
quality report, writes a plotfile and a checkpoint, and stops. Look at the largest slope and the smallest layer thickness. If STF
smoothing reports capped levels, lower ``terrain.gamma_m`` or raise the domain top.

Set the wind direction with inflow and outflow
----------------------------------------------

``bc.inflow_direction`` is the direction the wind blows *towards*, in degrees from the x axis, counter-clockwise.
Make the faces the wind enters through ``inflow`` and the faces it leaves through ``outflow``:

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Wind from
     - ``bc.inflow_direction``
     - Faces
   * - West
     - 0
     - ``bc.xlo = inflow``, ``bc.xhi = outflow``, ``ylo``/``yhi`` slip or periodic
   * - South-west
     - 45
     - ``bc.xlo = bc.ylo = inflow``, ``bc.xhi = bc.yhi = outflow``
   * - North
     - 270
     - ``bc.yhi = inflow``, ``bc.ylo = outflow``, ``xlo``/``xhi`` slip or periodic

A periodic direction needs ``mesh.is_periodic = 1`` in that direction and ``periodic`` on both of its faces. A
non-periodic direction needs ``mesh.is_periodic = 0``.

Choose the vertical grid
------------------------

- **Uniform layers:** set nothing; the layer thickness is the domain height divided by ``mesh.n_cell`` z.
- **Stretched layers to a fixed top:** set ``mesh.dz0`` (the first layer); the stretching ratio is computed so the
  layers reach ``mesh.prob_hi`` z.
- **Stretched layers with a fixed ratio:** also set ``mesh.stretch_ratio``; the top is wherever the layers end.
- **Your own levels:** list the ``mesh.n_cell`` z + 1 heights, starting at 0, in a file and set
  ``mesh.zlevels_file``.

Put a few cells inside the surface layer: a first cell of a few metres for typical terrain studies.

Add stratification and a damping layer
--------------------------------------

Give the potential-temperature sounding as heights and values:

.. code-block:: text

   bc.theta_profile_z = 0.   700.  1000.
   bc.theta_profile   = 300. 300.  309.
   bc.zhi_theta_grad  = 0.03          # K/m through the top
   forcing.gravity    = 9.81          # the default; 0 switches buoyancy off

In stable air, hills launch gravity waves. Absorb them with a Rayleigh layer near the top, on :math:`w` and
:math:`\theta` together:

.. code-block:: text

   bc.rayleigh_zdamp    = 700.        # start height (m)
   bc.rayleigh_dampcoef = 0.2         # strength (1/s)
   bc.rayleigh_fields   = w theta

With inflow and outflow faces, a sponge layer (``bc.sponge_width``, ``bc.sponge_strength``) damps reflections at the
sides.

Add Coriolis forcing and a geostrophic wind
-------------------------------------------

.. code-block:: text

   forcing.coriolis    = 1
   forcing.latitude    = 45.
   forcing.geostrophic = 1
   forcing.ug          = 10.
   forcing.vg          = 0.

Start such cases from the single-column profile, ``init.first_guess = column`` (:doc:`initial_conditions`). It
removes the slow inertial adjustment and typically halves the iterations.

.. _converge-faster:

Converge faster
---------------

On a CPU build the defaults are line LU-SGS with two sweeps and a CFL cap of 8, the setting with the fewest
iterations on Askervein (:ref:`askervein-erf-mesh`). A GPU build keeps z-line
Jacobi; give the same settings there explicitly. The hill cases were measured at a cap of 16:

.. code-block:: text

   steady.sweep        = lusgs        # the CPU default
   steady.zline_sweeps = 2            # the default
   steady.cfl_max      = 8.           # the default with lusgs; 16 on the hill cases
   steady.residual_ref_iterations = 5 # the default; a longer reference window holds the CFL at cfl0
   init.first_guess    = column       # with Coriolis forcing; otherwise profile

Then, depending on the case:

- **Larger terrain meshes:** add grid sequencing. It needs ``mesh.n_cell`` and ``mesh.max_grid_size`` in x and y
  divisible by :math:`2^{L-1}`:

  .. code-block:: text

     mesh.hierarchy_levels = 3
     steady.sequencing     = 1

- **The geostrophic wind has u and v of opposite sign:** add ``steady.line_alternate = 1``.
- **Thin first cells (aspect ratio far above 10) on flat or gentle terrain:** add ``simplec.w_vnn_cap = 16``. Over
  steep terrain the default slope margin (``simplec.w_vnn_cap_slope_margin = 4``) lowers the cap where needed; see
  :ref:`Slope margin <slope-margin>` in :doc:`results`.
- **Many small boxes:** the LU-SGS sweeps lag the coupling across box faces, so very small boxes
  (``mesh.max_grid_size``) cost extra iterations.

Fix a run that does not converge
--------------------------------

Look at the log first: which residual stays high, and does the CFL number keep dropping back?

- **Residuals grow after the CFL reaches its cap:** lower ``steady.cfl_max`` (for example 16 → 8, or 4 → 2 with
  z-line Jacobi).
- **Oscillation in stratified flow:** add a Rayleigh layer on ``w theta``; keep ``simplec.pressure_reltol`` at
  1e-3 or tighter.
- **The θ residual levels off well below its start:** a boundary layer under an inversion warms slowly without end.
  Lower ``steady.residual_drop_orders`` (for example to 2) and rely on ``steady.residual_tol``, as in
  :doc:`tutorial_flat_abl`.
- **Very slow decay of w on thin cells:** see the vertical-velocity cap above.

A run that stops at ``steady.max_iterations`` still writes its plotfile and checkpoint, so you can inspect it and
restart.

View and check the results
--------------------------

- **Visualization.** Open the plotfile directory in ParaView ("AMReX/BoxLib Grid Reader"), VisIt or yt.
- **Profiles in Python.** ``extern_erf/erf_plotfile.py`` reads plotfiles with the Python standard library only; see
  :doc:`tutorial_column`.
- **ERF's checks.** Set ``io.turb_plot_file``. ERF's scripts ``extern_erf/check_neutral.py``, ``check_hill.py`` and
  ``check_hill3d.py`` then read the turbulence plotfile directly.
- **Compare two runs.** AMReX's ``amrex_fcompare`` compares two plotfiles field by field.

Run on GPUs and on Kestrel
--------------------------

Build with CUDA for an NVIDIA H100:

.. code-block:: bash

   cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release -DGUST_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=9.0
   cmake --build build-gpu -j 8

Run one MPI rank per GPU. On a GPU, ``steady.line_order = diagonal`` lets the LU-SGS sweeps run one column per thread.

On NLR Kestrel, ``tools/kestrel/build_kestrel.sh cpu`` or ``gpu`` builds with the right modules (:doc:`building`).
