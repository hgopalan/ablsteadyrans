Tutorial 3: Flow over a hill
============================

This tutorial adds terrain: a 100 m high, two-dimensional ridge under the boundary layer of :doc:`tutorial_flat_abl`.
It is ERF's ``Neutral_Hill_2D`` case. You will run it twice, first with the standard settings and then with a
faster initial field. Each run takes 10 to 30 seconds.

**You will learn to:** define analytic terrain, use the faster line LU-SGS solver, damp waves near the top, and
start from a better initial field.

The physical problem
--------------------

The ridge is a Witch of Agnesi profile,

.. math::

   h(x) = \frac{h_0}{1 + (x/a)^2},

with height :math:`h_0 = 100` m and half-width :math:`a = 500` m. The wind speeds up over the crest and slows down in
the lee. The domain is 5120 m long and periodic, with one cell in y, so the flow is two-dimensional.

The input file
--------------

``Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d_lusgs`` keeps the sounding and forcing of tutorial 2 and
adds:

.. code-block:: text

   # A 2-D domain: 128 x 1 x 64 cells, centred on the hill
   mesh.n_cell        = 128 1 64
   mesh.prob_lo       = -2560. 0. 0.
   mesh.prob_hi       =  2560. 40. 1000.
   mesh.max_grid_size = 64 1

   # The terrain: a Witch of Agnesi ridge along y
   terrain.type            = agnesi
   terrain.hill_shape      = ridge_x
   terrain.hill_height     = 100.
   terrain.hill_half_width = 500.
   terrain.smoothing       = 0

   # A Rayleigh damping layer above 700 m, on w and theta
   bc.rayleigh_zdamp    = 700.
   bc.rayleigh_dampcoef = 0.2
   bc.rayleigh_fields   = w theta

   # Line LU-SGS with two sweeps and a CFL cap of 16
   steady.sweep           = lusgs
   steady.zline_sweeps    = 2
   steady.cfl_max         = 16.
   simplec.rotation_relax = 1.

- ``terrain.smoothing = 0`` uses the basic terrain-following mapping, as ERF's input file does
  (:doc:`terrain_coordinates`).
- The hill launches gravity waves in the stable air aloft. The Rayleigh layer absorbs them before they reflect
  from the top (:doc:`boundary_conditions`).
- Line LU-SGS (:doc:`solution_algorithm`) allows a CFL cap of 16, against 4 for z-line Jacobi.

Run the case
------------

.. code-block:: bash

   mkdir -p runs/hill && cd runs/hill
   mpirun -n 2 ../../build-rel/ablsteady ../../Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d_lusgs

Over terrain the initial field needs a Poisson solve to conserve mass. The log reports it:

.. code-block:: text

   ablsteady: initializer: 12 Poisson iterations (gmres_mlmg, reltol 1.000e-10), divergence 1.815778083e-01 before, 4.144210521e-12 after

and the run ends after about 3080 iterations:

.. code-block:: text

   ablsteady: steady after 3048 iterations (3048 over the run)

Check the speed-up
------------------

ERF's check script compares the crest speed-up :math:`U_\mathrm{top}/U_\mathrm{upstream} - 1` in the lowest cells:

.. code-block:: bash

   python3 ../../extern_erf/check_hill.py --physics steady_hill2d_lusgs_turb03077

.. code-block:: text

   crest speed-up (U_top/U_up - 1), k=0               0.594267            0.5        0.3 range yes
   crest speed-up (U_top/U_up - 1), k=1               0.499086            0.5        0.3 range yes
   crest speed-up (U_top/U_up - 1), k=2               0.432988            0.5        0.3 range yes
   0 check(s) failed

The wind near the ground is about 50 % faster at the crest than upstream.

A faster start: the single-column initial field
-----------------------------------------------

The log-law initial profile has no wind turning with height, so the solver spends many iterations on the slow
inertial adjustment. ``init.first_guess = column`` first solves the flat-ground column (as in tutorial 2, but in one
dimension, in a fraction of a second) and uses its Ekman spiral as the initial profile (:doc:`initial_conditions`).
The input file ``inputs_steady_hill2d_column`` differs only in this line:

.. code-block:: text

   init.first_guess = column

.. code-block:: bash

   mpirun -n 2 ../../build-rel/ablsteady ../../Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d_column

The log reports the column solve, then the 3-D run:

.. code-block:: text

   ablsteady: column guess: converged after 129 steps (0.018 s), 64 levels, u* 0.3785 m/s, u 10.0000 v -0.0000 m/s at the top
   ablsteady: steady after 1231 iterations (1231 over the run)

The run converges in 1231 iterations instead of 3048. The column profile is written to ``column_guess.dat``.

What you have seen
------------------

- Terrain is set in the ``terrain`` namespace; the mesh follows it automatically.
- Line LU-SGS (``steady.sweep = lusgs``) converges in fewer iterations than z-line Jacobi.
- For cases with Coriolis forcing, ``init.first_guess = column`` halves the iterations.

Next: :doc:`tutorial_askervein` uses real terrain and an inflow boundary.
