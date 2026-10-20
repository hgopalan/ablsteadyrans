Results
========

This page collects the measured accuracy and cost of each part. The cases are those of :doc:`cases`; the tests
are in :doc:`testing`.

Operator accuracy
-----------------

Manufactured flow over a 120 m cosine hill (slopes up to 0.38), BTF mapping, geometric layers, 16 × 12 × 8 to
64 × 48 × 32 cells:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Operator
     - L2 order
     - Max-norm order
   * - Convection, unlimited MUSCL
     - 1.99–2.10
     - 1.94–2.10
   * - Convection, ``van_leer`` / ``van_albada``
     - 1.5–2.2
     - about 1 (extrema)
   * - Convection, ``upwind``
     - 1.00
     - 0.96
   * - Scalar and momentum diffusion
     - 1.99–2.00
     - 1.96–1.99

The Poisson layer is second order away from sloped ground. Next to sloped ground and sloped Neumann sides it is first
order, as ERF's operator is.

Wall distance
-------------

Tucker's distance against the exact distance on ERF's decks:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Case
     - Max relative error
     - Max error, d < 100 m
   * - Flat ground (128 × 1 × 64)
     - 0.13 %
     - 0.03 m
   * - 2-D ridge, 100 m (128 × 1 × 64)
     - 5.3 %
     - 0.16 cells
   * - 3-D hill, 100 m (64 × 64 × 20)
     - 9.2 %
     - 0.13 cells

Initial field against FastWindTerrain
-------------------------------------

FastWindTerrain's bump case (100 m Gaussian hill on a 0.1 slope, power-law inflow): hilltop speed-up
:math:`|u_h|/S_0(h)`, ablsteady / FastWindTerrain.

.. list-table::
   :header-rows: 1
   :widths: 12 22 22 22 22

   * - :math:`\alpha_v`
     - 20 m
     - 40 m
     - 80 m
     - 150 m
   * - 1.0
     - 1.346 / 1.335
     - 1.253 / 1.247
     - 1.157 / 1.155
     - 1.084 / 1.083
   * - 0.5
     - 1.398 / 1.394
     - 1.240 / 1.236
     - 1.118 / 1.117
     - 1.052 / 1.052
   * - 0.1
     - 1.154 / 1.158
     - 1.088 / 1.092
     - 1.052 / 1.052
     - 1.032 / 1.032

The two codes agree to within 0.3–0.8 % of the speed-up at every height.

SIMPLEC runs to steady
----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Case
     - Iterations
     - Check
   * - 1-D neutral column (RT-1.15.1)
     - 2231
     - :math:`u_*` 0.3499999 against 0.35; first-cell :math:`|U|` and :math:`k` exact to 1e-6
   * - Lid-driven cavity, Re 100 (UT-1.15.3)
     - 1760
     - within 0.0083 of Ghia, Ghia and Shin (1982)
   * - ``Neutral_ABL_Flat`` (RT-1.15.2)
     - about 2850
     - :math:`u_*` 0.377 m/s; ERF's ``check_neutral.py --physics`` passes
   * - **2-D hill (RT-1.16.1)**
     - **4590**
     - crest speed-up 0.59 / 0.50 / 0.43; ERF's ``check_hill.py --physics`` passes
   * - **3-D hill (RT-1.16.2)**
     - **2092**
     - crest speed-up 0.39 / 0.27 / 0.21; ERF's ``check_hill3d.py --physics`` passes
   * - Smooth plateau (beyond the plan)
     - 8635
     - within 0.37 m/s of ERF's 24 h run in :math:`u`
   * - Askervein, coarse (RT-1.16.3)
     - 156
     - convergence only; the line A/AA comparison needs Kestrel

**Cost against ERF.** On the 2-D hill one iteration costs 1.3 ERF steps. The steady state takes 49 s on two ranks,
against 116 s for ERF's 6 h run and about 460 s for 24 h. On the 3-D hill the steady state takes 148 s on four
ranks, against 242 s for ERF's 4 h run, which is still a quarter of an inertial period from steady.

**Sweeps and CFL cap on the 2-D hill:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Sweeps / cap
     - Iterations
     - Wall (s)
     - Verdict
   * - 4 / 2
     - 20 000
     - 119
     - not steady
   * - 8 / 4
     - 5173
     - 71
     - steady
   * - **16 / 4**
     - **4590**
     - **49 (idle)**
     - **steady** (the deck)

**Line LU-SGS on the 2-D hill (sub-phase 2.1).** With two symmetric sweeps and a CFL cap of 16, the 2-D hill reaches
the verdict in 3062 iterations (23 s on two ranks) against 4590 iterations and 49 s for sixteen Jacobi sweeps at a
cap of 4. Every ERF check passes; the steady fields are within the verdict tolerance of the Jacobi result (:math:`u`
within 0.017 m/s and :math:`w` within 2e-4 m/s). On the 3-D hill, two LU-SGS sweeps at a cap of 16 take 823
iterations (51 s on four ranks) against 2092 iterations (148 s) for sixteen Jacobi sweeps at a cap of 4; ERF's
checks pass in both cases.

.. list-table:: Line LU-SGS on the 2-D hill (two ranks)
   :header-rows: 1
   :widths: 22 18 18

   * - Sweeps / cap
     - Iterations
     - Wall (s)
   * - 1 / 4
     - 6306
     - 41
   * - 2 / 4
     - 4937
     - 35
   * - 1 / 8
     - 4962
     - 34
   * - 2 / 8
     - 3619
     - 27
   * - 1 / 16
     - 4369
     - 31
   * - **2 / 16**
     - **3062**
     - **23**
   * - 2 / 32
     - 2886
     - 22
   * - 4 / 32
     - 5278
     - 84
   * - 2 / 64
     - 8431
     - 106
   * - 4 / 64
     - 20 000
     - 481 (not steady)
   * - 2 / 128
     - 20 000
     - 405 (not steady)

**The converged-column first guess (sub-phase 2.1b).** With ``init.first_guess = column`` (:doc:`initial_conditions`: the
deck's column solved as a 1-D problem in 0.02–0.06 s, its u, v, k broadcast by height above the ground) the
inertial spin-up of the profile guess is gone:

.. list-table:: The first guess and the verdict (LU-SGS decks; two ranks for the 2-D hill and the flat case, four
   for the 3-D hill; wall times with the machine shared)
   :header-rows: 1
   :widths: 30 26 16 16

   * - case
     - first guess
     - iterations
     - wall (s)
   * - 2-D hill
     - profile + projection
     - 3062
     - 23 (idle)
   * - 2-D hill
     - column + projection
     - 1237
     - 98 (shared)
   * - 2-D hill
     - column alone
     - 1202
     - 95 (shared)
   * - 3-D hill
     - profile + projection
     - 823
     - 51 (idle)
   * - 3-D hill
     - column + projection
     - 402
     - 1441 (shared, starved)
   * - Neutral_ABL_Flat
     - profile + projection
     - 2839
     - 13 (idle)
   * - Neutral_ABL_Flat
     - column + projection
     - not steady by 20 000
     - (:doc:`initial_conditions`)

On the 2-D hill the cross-wind residual's plateau (0.1–0.2 of its reference until iteration 1500 from the profile)
is gone: at iteration 300 the column run's u and v residuals are an order below the profile run's, and the verdict
comes at 1237 iterations against 3062; column alone and column with the projection differ by 3 %. The column
seed's *first* residual is larger than the profile's (u 5.8e-2 against 3.1e-2, v 3.5e-2 against 3.1e-3): the
lifted spiral carries an advective imbalance :math:`U\,\partial_z v\,h'` on the flanks that the profile guess (v
= 0) does not have and the projection does not remove; the balance shows after the CFL ramp. On the 3-D hill the verdict comes at 402 iterations against 823. On the flat case the column seed is the 1-D
steady state of u and v to round-off (first-iteration u residual 3e-11 against 3e-2), but its k at the inversion
base differs from the 3-D closure's and the run falls into a slow k oscillation there that does not reach the
verdict (:doc:`initial_conditions`); the flat deck keeps the profile.

**Askervein (RT-1.16.3, the convergence half).** The coarse mesh of RT-1.12.1 (30 × 30 × 18) with the log-law
inflow, the outflow and slip sides, no Coriolis, four sweeps, cap 2: the verdict after 156 iterations (2 s on two
ranks). The residual is largest in the two columns next to the outflow (100 times the interior, a checkerboard in
:math:`R_w`); with ``steady.residual_tol`` 1e-10 it keeps falling there (1e-9 after 706 iterations, the interior at
1e-13), so the outflow needs neither a sponge nor an outlet drag. The 5 % comparison with ERF's URANS on lines A
and AA needs the full mesh and the reference run (Kestrel) and is not in this version.

**Time per step** (idle machine). On the 2-D hill one SIMPLEC iteration with 16 sweeps costs 1.3 ERF steps (10.8
against 8.1 ms on two ranks); the steady state takes 4590 iterations against the 57 600 steps of ERF's 24 h (3.1
against ERF's 6 h physics run, which is not steady aloft): 49 s against 116 s (6 h) and about 460 s (24 h). On the
3-D hill one iteration costs 2.8 ERF steps (71 against 25 ms on four ranks) and the steady state 2092 iterations
against ERF's 9600 steps of 4 h: 148 s against 242 s, with ERF still a quarter of an inertial period from steady.

Newton–Krylov
-------------

On the Ekman column of UT-3.1.2, plain Newton reduces the scaled residual by 0.403, 0.171, 0.0131, 2.09e-4 and
8.26e-8 in its first five steps, with :math:`r_{n+1}/r_n^2` between 0.25 and 0.57 (quadratic convergence). At the
default Krylov tolerance 1e-2 every restart length from 10 up takes 3, 5 and 8 GMRES iterations in the first three
steps; restart 5 reaches the 2000-iteration cap on the third.

Ten plain Newton steps on the column at rtol 1e-2 (GMRES iterations in all):

.. list-table::
   :header-rows: 1

   * - Restart
     - 10
     - 30 (default)
     - 60
     - **200**
   * - Iterations
     - 4942 (cap reached on step 4)
     - 1469
     - 647
     - **346**

On a 30 m hill under the wall law (9 × 6 × 8, inflow and outflow, :math:`\Delta\tau` = 20 s, the defaults)
four steps take :math:`|R|` from 11.3 to 0.211, every solve meeting its tolerance.

Pseudo-transient continuation
-----------------------------

The continuation of sub-phase 3.2 (:doc:`solver`) on the same hill (six uneven boxes), on the hill with
:math:`u` three times larger, and on the Ekman column, GMRES(200), up to :math:`|R| \le 10^{-8}|R_0|`
(``Continuation.DISABLED_ParameterSweep``; steps accepted, GMRES iterations in all):

.. list-table::
   :header-rows: 1

   * - Settings
     - Hill
     - Hill, :math:`u \times 3`
     - Column
   * - **Defaults** (:math:`c` 1e-4, :math:`\beta` 0.5, :math:`p` 1, CFL cap 10, infinite :math:`\Delta\tau` at the cap)
     - **8, 976**
     - **8, 967**
     - **17, 510** (11 backtracks)
   * - No line search (``steady.newton_lambda_min = 0.99``)
     - 8, 976
     - 8, 967
     - fails at 3.0e-2
   * - ``steady.newton_dtau_inf_at_cap = 0``
     - 17, 1790
     - 16, 1395
     - 9.1e-4 after 60
   * - :math:`c` = 0.1
     - 8, 976
     - 8, 967
     - 15, 318
   * - :math:`c` = 0.4
     - 8, 976
     - 8, 967
     - fails on the first step
   * - :math:`\beta` = 0.3
     - 8, 976
     - 8, 967
     - 17, 437
   * - :math:`\beta` = 0.8
     - 8, 976
     - 8, 967
     - 19, 628
   * - :math:`p` = 0.5
     - 9, 1126
     - 9, 1088
     - fails at 3.5e-2
   * - :math:`p` = 2
     - 8, 976
     - 8, 967
     - 15, 355
   * - CFL cap 100
     - 9, 1119
     - 9, 1036
     - 22, 700
   * - CFL cap 1000
     - 10, 1282
     - 10, 1168
     - 20, 337
   * - ``steady.cfl0 = 10`` (Newton from the start)
     - 6, 940
     - 13, 2585
     - 11, 370
   * - Plain Newton (``NewtonKrylov::step``, no globalization)
     - 6 steps
     - diverges (:math:`|R|` 99.5, 377, 12 971, ...)
     - 9 steps, :math:`|R|` rises once

:math:`|R|` with the defaults: hill 11.3, 2.18, 0.842, 0.282, 0.124, 0.0403, 4.0e-4, 4.0e-6, 3.9e-8 (CFL 1, 2, 4,
8, then infinite); hill :math:`u \times 3` 99.5, 42.1, 12.3, 2.42, 0.928, 0.0734, 2.0e-3, 1.9e-5, 1.9e-7; column
0.760, 0.568, 0.387, 0.223, 0.0912, 0.0241, 0.0232, 0.0178, 0.0167, 1.1e-3, 1.06e-3, 9.95e-4, 9.61e-4,
9.38e-4, 4.5e-4, 4.5e-6, 5.2e-7, 5.0e-9 (eleven backtracks, the slow steps at :math:`\lambda < 1`). No setting showed limit cycling (a CFL alternating between cut and growth); the
failures are retries running out at a small CFL. GMRES(30) on the hill: 23 steps to 6.6e-9 without the switch to
infinite :math:`\Delta\tau`; with it, twelve steps at the 2000-iteration limit near 6e-3 and a failure at 7.3e-8.

In UT-3.2.1 (one box, GMRES(30), :math:`u \times 3`, plain Newton throughout) ``NewtonKrylov::step`` raises
:math:`|R|` by 3.53 on the first step, where :math:`k` goes negative. The continuation keeps :math:`k` at its floor,
takes three full steps (reductions 0.670, 0.783, 0.826) and then three steps at :math:`\lambda = 1/2`, where the full
step would raise :math:`|R|` by 2.77.
