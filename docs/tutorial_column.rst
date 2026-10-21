Tutorial 1: A neutral surface layer
===================================

This first tutorial runs the simplest case: wind over flat, uniform ground with no temperature effects. The
answer is known exactly, so you can check that ablsteady reproduces it. The run takes a few seconds on a laptop.

**You will learn to:** run ablsteady, read its log, find its output and compare a wind profile with theory.

**Before you start:** build ablsteady (:doc:`building`). The commands below assume you are in the repository
directory and the executable is ``build-rel/ablsteady``.

The physical problem
--------------------

A constant pressure gradient drives the wind over flat ground with roughness length :math:`z_0 = 0.1` m. The domain
is 1000 m high. At steady state the ground stress balances the driving force, so the friction velocity is known:

.. math::

   u_*^2 = F\,H \quad\Rightarrow\quad u_* = \sqrt{1.225\times10^{-4} \times 1000} = 0.35\ \mathrm{m/s}.

Near the ground the wind follows the logarithmic law and the TKE is constant:

.. math::

   u(z) = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad k = \frac{u_*^2}{C_{\mu 0}^2} = 0.396\ \mathrm{m^2/s^2}.

Nothing changes in the horizontal, so a small 4 × 4 × 64 mesh is enough.

The input file
--------------

Open ``Tests/Canonical_RANS/Neutral_Column/inputs_steady_column``. Without its comments it reads:

.. code-block:: text

   # Mesh: 4 x 4 x 64 cells, 2560 m x 2560 m x 1000 m, periodic sides
   mesh.n_cell        = 4 4 64
   mesh.prob_hi       = 2560. 2560. 1000.
   mesh.is_periodic   = 1 1 0
   mesh.max_grid_size = 4 4
   terrain.type = flat

   # Boundaries: rough ground with a wall model, slip top
   bc.zlo = surface_layer
   bc.zhi = slip_wall
   bc.surface_z0 = 0.1
   bc.theta_profile_z = 0.
   bc.theta_profile   = 300.

   # Initial field: a log-law profile (deliberately with the wrong u* = 0.25 m/s)
   bc.inflow_profile   = log_law
   bc.inflow_ustar     = 0.25
   bc.inflow_z0        = 0.1
   bc.inflow_direction = 0.
   init.first_guess = profile

   # Turbulence model and driving force
   turb.max_geom_lscale = 1.0e12
   turb.face_mean       = logarithmic
   forcing.pressure_gradient = 1.225e-4 0. 0.

   # Solver and convergence criteria
   steady.stage                = simplec
   steady.sweep                = zjacobi
   steady.cfl0                 = 1.
   steady.cfl_max              = 50.
   steady.max_iterations       = 20000
   steady.residual_drop_orders = 6
   steady.residual_tol         = 1.e-8

   # Output
   io.plot_file      = steady_column_plt
   io.turb_plot_file = steady_column_turb
   io.chk_file       = steady_column_chk

Every input is ``namespace.key = value``. The namespaces group the inputs by topic: ``mesh``, ``bc`` (boundary
conditions), ``init`` (initial conditions), ``steady`` (solver), and so on. :doc:`inputs` lists them all.

Three settings are worth noting:

- The initial profile uses :math:`u_* = 0.25` m/s, not 0.35 m/s. The solver has to find the right answer itself.
- ``steady.sweep = zjacobi`` keeps the z-line Jacobi sweeps of the iteration counts below. Without it a CPU build
  uses line LU-SGS (the default since the Phase 2 checkpoint, :doc:`results`).
- ``steady.cfl_max = 50`` is large. It works here because every column is identical; 3-D cases need smaller values
  (:doc:`howto`).

Run the case
------------

Run it in a new directory, on two MPI ranks:

.. code-block:: bash

   mkdir -p runs/column && cd runs/column
   mpirun -n 2 ../../build-rel/ablsteady ../../Tests/Canonical_RANS/Neutral_Column/inputs_steady_column

Read the log
------------

The log starts with the mesh and the initial field:

.. code-block:: text

   ablsteady: (4,4,64) cells, 1 boxes (each spans the full z column)
   ablsteady: initializer: 0 Poisson iterations (gmres_mlmg, reltol 1.000e-10), divergence 0.000000000e+00 ...
   ablsteady: mesh quality
     min cell thickness / flat layer thickness  = 1.000000000e+00
     ...
   ablsteady: SIMPLEC stage: 2 zjacobi sweeps, pressure gmres_mlmg reltol 0.001, alpha_p 1, w_vnn_cap off, cfl0 1, cfl_max 50

Over flat ground the initial profile already conserves mass, so the initializer needs no Poisson iterations. Then
the solver prints one line per iteration (the first ten, then every tenth):

.. code-block:: text

   ablsteady: it 300 cfl 5.000e+01 r 1.935e-02 res 9.669e-05 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.352e-03 mon 6.937e+00 2.165e-01 3.362e-01 poisson 0 ...

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Meaning
   * - ``it``
     - Iteration number.
   * - ``cfl``
     - Current CFL number. It grows from ``steady.cfl0`` to ``steady.cfl_max`` as the residual falls.
   * - ``r``
     - Largest residual relative to its reference value.
   * - ``res``
     - Scaled residuals of :math:`u, v, w, p, \theta, k`. The column has no :math:`v, w, p, \theta` imbalance, so
       those are 0.
   * - ``mon``
     - Monitors: the largest wind speed, the mean TKE and the mean friction velocity :math:`u_*`.
   * - ``poisson``, ``div``
     - Iterations of the pressure Poisson solve, and the divergence before and after the correction.
   * - ``clipped``
     - Cells where :math:`k` was raised to its floor.

The run ends with the verdict:

.. code-block:: text

   ablsteady: it 2231 cfl 5.000e+01 r 1.373e-07 res 6.863e-10 ... mon 7.228e+00 2.375e-01 3.500e-01 ...
   ablsteady: steady after 2231 iterations (2231 over the run)
   ablsteady: stage account: 2231 iterations, 3.459247747 s, 0 AU
   ablsteady: wrote steady_column_plt02231
   ablsteady: wrote steady_column_turb02231
   ablsteady: wrote steady_column_chk02231

The mean friction velocity (the last monitor) is 0.3500 m/s, the exact value.

The output
----------

The run directory now holds three outputs (:doc:`output`):

- ``steady_column_plt02231``: the solution (:math:`u, v, w, p, \theta, k` and the height), an AMReX plotfile;
- ``steady_column_turb02231``: turbulence fields with ERF's variable names;
- ``steady_column_chk02231``: a checkpoint for restarting.

Plotfiles open in ParaView, VisIt or yt.

Compare with the log law
------------------------

The repository includes a small plotfile reader that needs only the Python standard library. Save this script as
``profile.py`` in the run directory:

.. code-block:: python

   import math
   import sys
   sys.path.insert(0, "../../extern_erf")   # ERF's plotfile reader
   import erf_plotfile

   z, avg, _ = erf_plotfile.planar_averages("steady_column_plt02231", ["u", "tke"])
   ustar, kappa, z0 = 0.35, 0.41, 0.1
   for k in range(6):
       loglaw = ustar / kappa * math.log((z[k] + z0) / z0)
       print(f"z = {z[k]:6.2f} m   u = {avg['u'][k]:.3f}   log law = {loglaw:.3f}   k = {avg['tke'][k]:.4f}")

Run it with ``python3 profile.py``:

.. code-block:: text

   z =   7.81 m   u = 3.731   log law = 3.731   k = 0.3960
   z =  23.44 m   u = 4.620   log law = 4.662   k = 0.4382
   z =  39.06 m   u = 5.029   log law = 5.097   k = 0.4067
   z =  54.69 m   u = 5.302   log law = 5.383   k = 0.3879
   z =  70.31 m   u = 5.507   log law = 5.597   k = 0.3745
   z =  85.94 m   u = 5.670   log law = 5.768   k = 0.3636

The first cell matches the log law exactly, because the wall model imposes it there, and :math:`k` equals
:math:`u_*^2/C_{\mu 0}^2 = 0.396`. Higher up the wind falls slightly below the log law (by about 1–2 %), because the
stress decreases with height in this driven column.

Next steps
----------

- :doc:`tutorial_flat_abl` adds Coriolis forcing and a temperature inversion.
- Try ``mpirun -n 1`` instead of 2: the result is the same.
