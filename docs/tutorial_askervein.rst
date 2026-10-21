Tutorial 4: Real terrain — Askervein hill
=========================================

Askervein is a 116 m high hill on the island of South Uist, Scotland. A field campaign in 1982–83 measured the wind
on masts along lines across the hill, and it has been a standard test for terrain flow models since. This tutorial
runs it on a coarse mesh in a few seconds.

**You will learn to:** read terrain from a file, set up inflow and outflow boundaries, use a stretched vertical
grid and choose the mass-consistent transmissivities.

The setup
---------

The wind blows from the west (+x) with a log-law profile, :math:`u_* = 0.6` m/s and :math:`z_0 = 0.03` m. The air
is neutral and there is no Coriolis force, as in ERF's input file. The input file is
``Tests/Real_Terrain/Askervein/inputs_steady_askervein``:

.. code-block:: text

   # 6 km x 6 km, 30 x 30 x 18 cells; 20 m first layer stretched by 1.08
   mesh.n_cell        = 30 30 18
   mesh.prob_lo       = -3000. -3000. 0.
   mesh.prob_hi       =  3000.  3000. 440.
   mesh.max_grid_size = 16
   mesh.dz0           = 20.
   mesh.stretch_ratio = 1.08

   # Terrain from ERF's terrain file (given on the command line)
   terrain.type      = erf_grid
   terrain.smoothing = 0

   # Inflow on the west face, outflow on the east, slip walls on the sides
   bc.xlo = inflow
   bc.xhi = outflow
   bc.ylo = slip_wall
   bc.yhi = slip_wall
   bc.zlo = surface_layer
   bc.zhi = slip_wall
   bc.surface_z0       = 0.03
   bc.inflow_profile   = log_law
   bc.inflow_ustar     = 0.6
   bc.inflow_z0        = 0.03
   bc.inflow_direction = 0.

   # Neutral air: uniform theta and no gravity
   bc.theta_profile_z  = 0.
   bc.theta_profile    = 300.
   forcing.gravity     = 0.

   # Mass-consistent initial field with alpha_v / alpha_h = 0.5
   init.first_guess = profile
   init.alpha_h     = 1.0
   init.alpha_v     = 0.5

- **The mesh is not periodic.** ``mesh.is_periodic`` is left at its default ``0 0 0``, and every side face gets a
  boundary type.
- **Stretched layers.** The first layer is 20 m thick; each layer above is 1.08 times thicker
  (:doc:`terrain_coordinates`).
- **Transmissivities.** :math:`\alpha_v/\alpha_h = 0.5` lets the initial field go partly around the hill rather than
  only over it (:doc:`initial_conditions`).

Run the case
------------

The terrain file ``extern_erf/askervein.txt`` is passed on the command line, which adds the key
``terrain.file`` to the inputs:

.. code-block:: bash

   mkdir -p runs/askervein && cd runs/askervein
   mpirun -n 2 ../../build-rel/ablsteady ../../Tests/Real_Terrain/Askervein/inputs_steady_askervein \
       terrain.file=../../extern_erf/askervein.txt

The log reports the mass-consistent adjustment, and that it used the ``gmres_fft`` Poisson back end because
:math:`\alpha_v \ne \alpha_h`:

.. code-block:: text

   ablsteady: initializer: 13 Poisson iterations (gmres_fft, reltol 1.000e-10), divergence 4.048451025e-01 before, 6.095532275e-12 after
   ablsteady: initializer: init.alpha_v / init.alpha_h = 0.5 needs the split z-face coefficient, which only solver.poisson_backend = gmres_fft has; ...

The run converges in 175 iterations:

.. code-block:: text

   ablsteady: steady after 175 iterations (175 over the run)

Look at the result
------------------

Open ``steady_askervein_plt00156`` in ParaView or VisIt and color a horizontal slice near the ground by ``u``. The
wind speeds up over the hilltop and slows in the lee.
The variable ``z`` holds the height of each cell centre.

Toward a validation run
-----------------------

This mesh has 200 m horizontal cells, too coarse to compare with the mast data. ERF's input file uses 300 × 300 cells.
To refine, increase ``mesh.n_cell`` and reduce ``mesh.dz0``, and see :ref:`converge-faster` for the solver
settings that suit larger meshes.

What you have seen
------------------

- Real terrain comes from a file (``terrain.type = erf_grid`` or ``xyz``).
- Inflow and outflow faces replace periodic sides.
- The mass-consistent transmissivities shape the initial field over complex terrain.
