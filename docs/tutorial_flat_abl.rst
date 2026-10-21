Tutorial 2: An atmospheric boundary layer
=========================================

This tutorial adds two features of a real atmosphere to :doc:`tutorial_column`: the Earth's rotation (Coriolis
force) and a temperature inversion that caps the boundary layer. It is ERF's ``Neutral_ABL_Flat`` case, solved
directly for the steady state. The run takes about 15 seconds.

**You will learn to:** drive the flow with a geostrophic wind, prescribe a temperature sounding and check the result
with ERF's own check script.

The physical problem
--------------------

Above the boundary layer, the wind is in geostrophic balance: the Coriolis force balances the large-scale pressure
gradient, and the wind is :math:`(u_g, v_g) = (10, 0)` m/s. Near the ground, friction slows the wind and turns it
to the left (in the northern hemisphere). The result is the Ekman spiral.

The potential temperature is 300 K up to 700 m and rises above, 0.03 K/m. This stable inversion limits the depth
of the boundary layer.

The input file
--------------

``Tests/Canonical_RANS/Neutral_ABL_Flat/inputs_steady_neutral_flat`` differs from tutorial 1 in these lines:

.. code-block:: text

   # A larger horizontal domain, split into four boxes
   mesh.n_cell        = 8 8 64
   mesh.max_grid_size = 4 4

   # Temperature sounding: 300 K to 700 m, 309 K at 1000 m, 0.03 K/m through the top
   bc.theta_profile_z = 0. 700. 1000.
   bc.theta_profile   = 300. 300. 309.
   bc.zhi_theta_grad  = 0.03

   # Coriolis forcing with a geostrophic wind of 10 m/s (f = 1e-4 1/s)
   forcing.coriolis        = 1
   forcing.coriolis_3d     = 0
   forcing.latitude        = 90.
   forcing.rotation_period = 125663.7061435917
   forcing.geostrophic     = 1
   forcing.ug              = 10.
   forcing.vg              = 0.

   # Initial profile: the log law with u* = 0.35 m/s (tutorial 1: 0.25 m/s)
   bc.inflow_ustar     = 0.35

   # Solver: four z-line sweeps, CFL cap 2
   steady.sweep                   = zjacobi
   steady.cfl_max                 = 2.
   steady.zline_sweeps            = 4
   steady.residual_ref_iterations = 100
   steady.residual_drop_orders    = 2
   turb.tke_floor = 1.e-6

The pressure-gradient force of tutorial 1 is replaced by the geostrophic forcing. The lines
``turb.max_geom_lscale`` and ``turb.face_mean`` are removed, so ERF's defaults apply: a length-scale cap of 30 m and
arithmetic face averages.

**Why only two orders of residual drop?** The boundary layer slowly entrains the warm air of the inversion and warms
(by about 0.01 K per day) without ever stopping. The :math:`\theta` and :math:`w` residuals therefore level off
about 2.5 orders below their starting value. A drop of 2 orders, together with ``steady.residual_tol``, still gives a
converged wind field.

Run the case
------------

.. code-block:: bash

   mkdir -p runs/flat && cd runs/flat
   mpirun -n 2 ../../build-rel/ablsteady ../../Tests/Canonical_RANS/Neutral_ABL_Flat/inputs_steady_neutral_flat

The run ends with

.. code-block:: text

   ablsteady: it 2840 cfl 2.000e+00 r 6.816e-03 res 1.793e-05 6.202e-06 5.720e-06 6.324e-08 7.754e-06 1.488e-05 mon 1.114e+01 1.346e-01 3.774e-01 ...
   ablsteady: steady after 2840 iterations (2840 over the run)

Unlike tutorial 1, all six residuals are non-zero: the Coriolis force couples :math:`u` and :math:`v`, and the
inversion couples :math:`w` and :math:`\theta`. The largest wind speed is 11.1 m/s, slightly above the geostrophic
wind (a low-level jet below the inversion), and :math:`u_* = 0.377` m/s.

Check the result with ERF's script
----------------------------------

ERF provides a check script for this case. It reads the turbulence plotfile, which uses ERF's variable names:

.. code-block:: bash

   python3 ../../extern_erf/check_neutral.py --physics steady_neutral_flat_turb02840 steady_neutral_flat_turb02840/surf_hist.dat

It prints one line per check and ends with ``0 check(s) failed``. Among the checks:

.. code-block:: text

   u_star [m/s] plausible                         0.377402           0.35       0.15   abs yes
   KE(k=0)/u*^2 vs 1/Cmu0^2                        3.23251         3.2325       0.05   rel yes
   |U|(k=0) vs log law                              4.0235         4.0235        0.1   rel yes

What you have seen
------------------

- Coriolis forcing and a sounding are set in the ``forcing`` and ``bc`` namespaces.
- With rotation and stratification the solver needs a smaller CFL cap than the column (2 instead of 50).
- ERF's check scripts work directly on ablsteady's turbulence plotfile.

Next: :doc:`tutorial_hill` adds terrain.
