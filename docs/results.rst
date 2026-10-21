Results
========

This page collects the measured accuracy and cost of ablsteady. Each study states its setup, shows the numbers in a
table and lists the main findings. Recommended settings are in **bold**. The cases are described in
:doc:`test_cases`; the tests that check these numbers are listed in :doc:`testing`.

Summary of recommended settings
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 30 35

   * - Setting
     - Recommended
     - Section
   * - Line solver
     - ``steady.sweep = lusgs``, 2 sweeps
     - `Line LU-SGS`_
   * - CFL cap
     - ``steady.cfl_max = 16``
     - `Line LU-SGS`_
   * - Initial field, terrain with Coriolis
     - ``init.first_guess = column``
     - `Initial field: profile or column`_
   * - Initial field, flat with an inversion
     - ``init.first_guess = profile``
     - `Initial field: profile or column`_
   * - Thin first cells
     - ``simplec.w_vnn_cap = 16`` (slope margin 4, the default)
     - `Vertical-velocity cap`_
   * - Grid sequencing (terrain)
     - ``steady.sequencing = 1``, 2–3 levels
     - `Grid sequencing`_
   * - Poisson back end
     - ``gmres_mlmg`` (default)
     - `Poisson back end`_

Verification
------------

Spatial accuracy
^^^^^^^^^^^^^^^^

*Setup.* Method of manufactured solutions on a 120 m cosine hill (slopes up to 0.38), BTF mapping, geometric layers,
meshes from 16 × 12 × 8 to 64 × 48 × 32 cells.

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Operator
     - Order, L2 norm
     - Order, max norm
   * - Convection, unlimited MUSCL
     - 1.99–2.10
     - 1.94–2.10
   * - Convection, ``van_leer`` / ``van_albada``
     - 1.5–2.2
     - about 1 (at extrema)
   * - Convection, ``upwind``
     - 1.00
     - 0.96
   * - Scalar and momentum diffusion
     - 1.99–2.00
     - 1.96–1.99

*Findings.*

- The operators reach their design order.
- The Poisson operator is second order away from sloped ground and first order next to sloped ground and sloped
  Neumann sides, as ERF's operator is.
- A uniform wind has zero discrete divergence on all three terrain mappings, to 1e-12.

Wall distance
^^^^^^^^^^^^^

*Setup.* Tucker's wall distance compared with the exact distance on ERF's meshes.

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Case
     - Max relative error
     - Max error for d < 100 m
   * - Flat ground (128 × 1 × 64)
     - 0.13 %
     - 0.03 m
   * - 2-D ridge, 100 m (128 × 1 × 64)
     - 5.3 %
     - 0.16 cells
   * - 3-D hill, 100 m (64 × 64 × 20)
     - 9.2 %
     - 0.13 cells

Mass-consistent initial field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* FastWindTerrain's bump case: a 100 m Gaussian hill on a 0.1 slope with a power-law inflow. The table gives
the hilltop speed-up :math:`|u_h|/S_0(h)` at four heights, as ablsteady / FastWindTerrain.

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

*Finding.* The two codes agree to within 0.3–0.8 % at every height.

Steady solutions
----------------

Convergence on the test cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* Each case runs to the steady verdict (:ref:`convergence-criteria`) with the settings of its input file,
on two ranks (four for the 3-D hill).

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Case
     - Iterations
     - Check against reference
   * - 1-D neutral column
     - 2231
     - :math:`u_*` = 0.3499999 against 0.35; first-cell speed and :math:`k` exact to 1e-6
   * - Lid-driven cavity, Re = 100
     - 1760
     - within 0.0083 of Ghia, Ghia and Shin (1982)
   * - ``Neutral_ABL_Flat``
     - 2840
     - :math:`u_*` = 0.377 m/s; ERF's ``check_neutral.py --physics`` passes
   * - 2-D hill
     - 4590 (Jacobi); 3077 (LU-SGS)
     - crest speed-up 0.59 / 0.50 / 0.43; ERF's ``check_hill.py --physics`` passes
   * - 3-D hill
     - 2092 (Jacobi); 823 (LU-SGS)
     - crest speed-up 0.39 / 0.27 / 0.21; ERF's ``check_hill3d.py --physics`` passes
   * - Smooth plateau
     - 8635
     - within 0.37 m/s of ERF's 24 h run in :math:`u`
   * - Askervein, coarse mesh
     - 156
     - convergence only (below)

The column, flat, 2-D hill (LU-SGS) and Askervein counts were re-measured on ``main`` at commit ``7bf696c``
(October 2026). The other numbers on this page were measured when each feature was added; later changes to the
code shift iteration counts by up to about 1 % (for example 3062 → 3077 on the 2-D hill).

The crest speed-ups are those of the z-line Jacobi runs (16 sweeps, CFL cap 4). The LU-SGS runs reach the same
steady state to the convergence tolerance, and ERF's checks pass for both (`Line LU-SGS`_).

**Askervein.** The coarse mesh (30 × 30 × 18) with a log-law inflow, outflow and slip sides and no Coriolis reaches
the verdict after 156 iterations (2 s on two ranks). The residual is largest in the two columns next to the outflow,
but keeps falling there with a tighter tolerance (1e-9 after 706 iterations), so the outflow needs no sponge. The
comparison with ERF's unsteady RANS on lines A and AA uses ERF's own setup on the full mesh (`Phase 2 checkpoint`_).

Cost compared with ERF
^^^^^^^^^^^^^^^^^^^^^^

*Setup.* Idle machine; ERF run at the pinned commit with its own input files.

.. list-table::
   :header-rows: 1
   :widths: 24 22 18 36

   * - Case
     - ablsteady
     - Cost of one iteration
     - ERF
   * - 2-D hill, 2 ranks
     - **23 s** (LU-SGS, timed at 3062 it.); 49 s (Jacobi, 4590 it.)
     - 1.3 ERF steps (Jacobi)
     - 116 s for 6 h; about 460 s for 24 h
   * - 3-D hill, 4 ranks
     - **51 s** (LU-SGS, 823 it.); 148 s (Jacobi, 2092 it.)
     - 2.8 ERF steps (Jacobi)
     - 242 s for 4 h, still a quarter of an inertial period from steady

*Finding.* ablsteady reaches a converged steady state in a fraction of the time ERF needs for a run that is not yet
steady.

Solver settings
---------------

z-line Jacobi: sweeps and CFL cap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* Largest stable CFL number from a Fourier analysis of the z-line Jacobi iteration
(``tools/sweeps_cfl_scan.py``), and the 2-D hill with three settings.

.. list-table::
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Sweeps per iteration
     - 2
     - 4
     - 8
     - 16
   * - Largest stable CFL
     - 0.9
     - 2.4
     - 3.7
     - 6.1

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Sweeps / CFL cap
     - Iterations
     - Wall time (s)
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
     - **49**
     - **steady**

Line LU-SGS
^^^^^^^^^^^

*Setup.* 2-D hill, two ranks, ``steady.sweep = lusgs``.

.. list-table::
   :header-rows: 1
   :widths: 30 20 20

   * - Sweeps / CFL cap
     - Iterations
     - Wall time (s)
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
     - 20 000 (not steady)
     - 481
   * - 2 / 128
     - 20 000 (not steady)
     - 405

*Findings.*

- Two sweeps at a CFL cap of 16 converge in 3062 iterations and 23 s, against 4590 iterations and 49 s for 16 Jacobi
  sweeps at a cap of 4.
- On the 3-D hill, LU-SGS takes 823 iterations (51 s, four ranks) against 2092 (148 s) for Jacobi.
- The steady fields agree with the Jacobi result to the convergence tolerance (:math:`u` within 0.017 m/s), and all
  ERF checks pass.
- Caps above 32 slow the iteration down.

Alternating sweep direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* ``steady.line_alternate`` on four ranks, LU-SGS with ``line_ghosts = zero``, an overlap of 2, diagonal order
and two sweeps. The first four rows are the truly neutral 3-D hill, run to round-off; the last two run to their input
file's convergence criteria (two orders).

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - Case
     - Off
     - On
     - z-line Jacobi
   * - :math:`(u_g, v_g) = (10, 0)`, CFL cap 16
     - 1792
     - 1813
     - 3495
   * - :math:`(u_g, v_g) = (10, 0)`, CFL cap 64
     - 2997
     - 2991
     - —
   * - :math:`(u_g, v_g) = (10, -10)`, CFL cap 16
     - not converged at 20 000
     - **3052**
     - 6051
   * - :math:`(u_g, v_g) = (10, -10)`, CFL cap 64
     - not converged at 20 000
     - 3585
     - —
   * - Conventionally neutral 3-D hill, cap 16
     - 1100
     - 1146
     - —
   * - 2-D hill
     - 3062
     - 3061
     - —

*Findings.*

- When :math:`u` and :math:`v` have opposite signs, the plain LU-SGS sweeps stall; alternating sweeps converge in
  half the Jacobi iterations.
- When the wind keeps one sign, alternating costs 0–4 % more iterations. The default is off.

Box count and overlap
^^^^^^^^^^^^^^^^^^^^^

*Setup.* Conventionally neutral 3-D hill (64 × 64 × 20, periodic), LU-SGS with two sweeps at a cap of 16. One box on
one rank; 4 and 16 boxes on four ranks. Iterations to the verdict.

.. list-table::
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - Boxes
     - ``lagged``
     - ``zero``, overlap 0
     - ``zero``, overlap 2
     - ``zero``, overlap 4
   * - 1
     - **860**
     - 998
     - 948
     - 926
   * - 4
     - **823**
     - 1205
     - 1100
     - 1020
   * - 16
     - **974**
     - 1552
     - 1310
     - 1270

*Findings.*

- On CPUs, ``zero`` needs 16–59 % more iterations than ``lagged`` without overlap, and 8–30 % more with an overlap
  of 4. The communication it saves is expected to matter on GPUs and many ranks; that is not yet measured.
- More boxes cost more iterations in both modes.
- ``line_order = diagonal`` gives the same iterates and, on CPUs, the same wall time as ``sequential``.

Initial field: profile or column
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* LU-SGS input files; two ranks for the 2-D hill and the flat case, four for the 3-D hill. The column itself
takes 0.02–0.06 s. Wall times marked "shared" were measured on a busy machine and are not comparable.

.. list-table::
   :header-rows: 1
   :widths: 30 26 16 16

   * - Case
     - Initial field
     - Iterations
     - Wall time (s)
   * - 2-D hill
     - profile + adjustment
     - 3062
     - 23
   * - 2-D hill
     - **column + adjustment**
     - **1237**
     - 98 (shared)
   * - 2-D hill
     - column only
     - 1202
     - 95 (shared)
   * - 3-D hill
     - profile + adjustment
     - 823
     - 51
   * - 3-D hill
     - **column + adjustment**
     - **402**
     - 1441 (shared)
   * - ``Neutral_ABL_Flat``
     - **profile + adjustment**
     - **2839**
     - 13
   * - ``Neutral_ABL_Flat``
     - column + adjustment
     - not steady at 20 000
     - —

*Findings.*

- On the hills the column start removes the slow inertial adjustment and halves the iterations or better.
- The first residual of the column start is larger (an advective imbalance on the hill flanks), but it falls faster
  once the CFL ramp is over.
- On the flat case the column's :math:`k` at the inversion base differs from the 3-D solution, and the run
  oscillates there. Use ``profile`` for that case.

Vertical-velocity cap
^^^^^^^^^^^^^^^^^^^^^

*Setup.* Askervein mesh (140 × 140 cells at 43 m, layers stretched by 1.08), LU-SGS with two sweeps at a CFL cap of
16, ``gmres_mlmg``. Iterations to the verdict for several caps :math:`V` (``simplec.w_vnn_cap``). "Diverged" means
the residuals rose to 0.2–0.8 after the CFL reached its cap; "—" means not run. :math:`V = 4` and 8 were measured
before the lower limit of 16 was introduced; they can no longer be set. This table was measured without the slope
margin (``simplec.w_vnn_cap_slope_margin = 0``) on earlier builds; `Slope margin`_ below has the default, measured
later (the Askervein 1 m run with :math:`V = 16` and no margin now takes 219 iterations, and the 0.4 m run with
:math:`V = 64` 338).

.. list-table::
   :header-rows: 1
   :widths: 28 12 12 12 12 12 12

   * - Case
     - Off
     - :math:`V = 4`
     - :math:`V = 8`
     - :math:`V = 16`
     - :math:`V = 32`
     - :math:`V = 64`
   * - Flat, 0.4 m first cell
     - 1620
     - 194
     - —
     - **196**
     - —
     - —
   * - Flat, 1 m first cell
     - 485
     - —
     - —
     - **216**
     - —
     - —
   * - Flat, 1 m, CFL cap 32
     - 855
     - —
     - —
     - **285**
     - —
     - —
   * - Hill, 1 m first cell
     - 324
     - diverged
     - diverged
     - **200**
     - 200
     - 201
   * - Hill, 0.4 m first cell
     - about 2800
     - —
     - —
     - diverged
     - —
     - 278

*Findings.*

- On flat ground the cap reduces the iterations 2.2 to 8.3 times.
- The capped and uncapped steady states agree to the convergence tolerance: 2.7e-7 m/s in the velocity at 0.4 m, and
  3e-6 in the 10 m speed-ups on lines A, AA and B at 1 m.
- Over steep terrain with very thin cells a small cap without the slope margin can diverge. The margin
  (`Slope margin`_) removes this.
- With z-line Jacobi (four sweeps, CFL cap 4) the cap has little effect: at that CFL the VNN number is small.

.. _slope-margin:

Slope margin
""""""""""""

*Setup.* As above, on this change's branch (after #38), with the cap :math:`V = 16` and several slope margins :math:`s`
(``simplec.w_vnn_cap_slope_margin``, :doc:`solution_algorithm`). "Off" is the run without the
cap. The steep ridge is ``Tests/Canonical_RANS/Steep_Ridge_Thin_Cell`` (32.5° maximum slope, 0.4 m first cell); its
verdict is 3 orders, the Askervein verdicts are the default. "Stalled" means the residual ratio stayed between 0.1
and 0.7 and the run was stopped.

.. list-table::
   :header-rows: 1
   :widths: 30 12 14 12 12 12 12

   * - Case
     - Off
     - :math:`s = 0`
     - :math:`s = 1`
     - :math:`s = 2`
     - :math:`s = 4`
     - :math:`s = 8`
   * - Askervein, 1 m first cell
     - 324
     - 219
     - —
     - 219
     - **219**
     - —
   * - Askervein, 0.4 m first cell
     - about 2800
     - diverged
     - —
     - stalled
     - **225**
     - 225
   * - Steep ridge, 0.4 m first cell
     - 3983
     - not steady at 3500
     - stalled
     - 2714
     - **2568**
     - 3111

*Findings.*

- The margin changes only the iteration path. On flat ground the cap is unchanged and the results are
  bit-identical. On the ridge (6-order verdict) the steady states with the cap off, with :math:`s = 2` and with
  :math:`V = 64` and no margin agree to 3e-8 m/s.
- :math:`s = 1` (the bound itself) is not enough. The bound assumes frozen coefficients and the horizontal and
  vertical coefficients both at :math:`1/\Delta\tau` (:doc:`solution_algorithm`). :math:`s = 2` stalls on
  Askervein at 0.4 m. :math:`s = 4` is the smallest margin tried that converges on all three cases.
- On Askervein at 0.4 m, :math:`s = 4` needs 225 iterations against 338 with the larger cap :math:`V = 64` and no
  margin.
- On the ridge about 85 % of the iterations are a transient that is the same in every configuration: the lee
  separation forms while the SER law lowers the CFL to about 0.08. The margin shortens only the tail.

Poisson back end
^^^^^^^^^^^^^^^^

*Setup.* Askervein terrain at 43 m horizontal spacing (140 × 140 cells, layers stretched by 1.08 to about 750 m),
LU-SGS with two sweeps at a CFL cap of 16, ``init.alpha_v = init.alpha_h = 1``. Mean Poisson iterations per
pressure correction.

.. list-table::
   :header-rows: 1
   :widths: 20 15 20 20 25

   * - First layer (m)
     - dx/dz
     - SIMPLEC iterations
     - ``gmres_fft``
     - ``gmres_mlmg``
   * - 8
     - 5
     - 196
     - 6.6
     - **2.6**
   * - 4
     - 11
     - 203
     - 7.7
     - **2.6**
   * - 2
     - 21
     - 230
     - 9.6
     - **2.5**
   * - 1
     - 43
     - 324
     - 12.8
     - **2.1**

*Findings.*

- Both back ends give the same SIMPLEC iterations.
- ``gmres_mlmg`` keeps the Poisson work per correction flat as the first layer thins; ``gmres_fft`` needs more.
- The back ends stop on different norms (volume-weighted for ``gmres_mlmg``). On these runs ``gmres_mlmg`` left 10–30
  times less divergence at the same ``simplec.pressure_reltol``.

Grid sequencing
^^^^^^^^^^^^^^^

*Setup.* ``steady.sequencing = 1`` (:doc:`solution_algorithm`). Level-0 iterations to the normal convergence criteria,
with the coarse-level iterations in brackets (coarsest first). Default settings unless noted.

.. list-table::
   :header-rows: 1
   :widths: 40 15 25 20

   * - Case
     - 1 level
     - With sequencing
     - Levels
   * - ``Neutral_Hill_3D``
     - 2091
     - **706** (383, 51)
     - 3
   * - Truly neutral 3-D hill (to round-off)
     - 3495
     - **2688** (239, 38)
     - 3
   * - Askervein, neutral
     - 156
     - **130** (100)
     - 2
   * - Askervein, conventionally neutral
     - 220
     - **124** (100)
     - 2
   * - Askervein, 1 m first cell, ``w_vnn_cap = 16``, no slope margin
     - 200
     - **157** (100)
     - 2
   * - Askervein, 1 m first cell, ``w_vnn_cap = 16``, no slope margin
     - 200
     - **151** (100, 40)
     - 3

Wall time on four ranks, ``Neutral_Hill_3D``: 145 s without sequencing, 57 s with three levels and the reference
hold on every level (``steady.sequencing_hold = 1``), 54 s with the default.

.. list-table:: Effect of the reference hold and the starting CFL (level-0 iterations)
   :header-rows: 1
   :widths: 40 20 20 20

   * - Case
     - ``sequencing_hold = 1``
     - **Default**
     - ``cfl_start = previous``
   * - ``Neutral_Hill_3D``, 3 levels
     - 774 (383, 127)
     - **706 (383, 51)**
     - 704 (383, 51)
   * - Truly neutral 3-D hill, 3 levels
     - 2773 (239, 101)
     - **2688 (239, 38)**
     - 2709 (239, 28)
   * - Askervein, neutral, 2 levels
     - 164 (100)
     - **130 (100)**
     - 130 (100)
   * - Askervein, conventionally neutral, 2 levels
     - 147 (100)
     - **124 (100)**
     - 124 (100)

*Findings.*

- Three levels cut the level-0 iterations of the 3-D hill by a factor of three and the wall time by more than half.
- Skipping the reference hold on the seeded levels (the default) saves most on the coarse levels (127 → 51,
  101 → 38).
- The gain is the transient the coarse levels remove. A deep stop (round-off) gains less in proportion than a
  two-order stop, since the final convergence rate is that of level 0.
- The starting field still matters: on ``Neutral_Hill_3D`` with three levels (and the hold on every level), the
  default projection needs 774 level-0 iterations, no projection 1442, and :math:`\alpha_v = 0.1` 1007.
- The converged state matches the run without sequencing to 1e-8 (``rt_steady_hill3d_neutral_sequencing``).
- Equal verdicts do not mean equal end states: each run stops at the first iteration that meets the criteria, so
  final residuals can differ by up to an order of magnitude within the same criteria.

Phase 2 checkpoint
^^^^^^^^^^^^^^^^^^

*Setup.* Askervein on ERF's mesh (300 × 300 × 18 cells, 20 m, a 20 m first layer stretched by 1.08), the neutral
deck of RT-1.16.3 (``inputs_steady_askervein``, wind along x) unless noted, 8 ranks, 25 boxes of 64 × 64 columns.
Iterations to the steady verdict, measured on one CPU (Apple M2 Pro) with other jobs running; wall times and
allocation units are left for Kestrel. "Not steady" means no verdict within the iteration limit (1000 or 2000).

.. list-table:: Sweep and CFL cap (PT-2.6.1; the deck's 100-iteration reference)
   :header-rows: 1
   :widths: 30 10 10 10 10 10 10

   * - Sweep
     - Cap 2
     - Cap 4
     - Cap 8
     - Cap 16
     - Cap 32
     - Cap 64
   * - ``zjacobi``, 2 sweeps
     - —
     - not steady
     - —
     - —
     - —
     - —
   * - ``zjacobi``, 4 sweeps
     - 332
     - 260
     - 246
     - 426
     - not steady
     - —
   * - ``zjacobi``, 8 sweeps
     - —
     - —
     - 205
     - —
     - —
     - —
   * - ``zjacobi``, 16 sweeps
     - —
     - 210
     - 190
     - 271
     - 649
     - —
   * - ``lusgs``, 1 sweep
     - —
     - —
     - 201
     - 265
     - —
     - —
   * - **lusgs, 2 sweeps**
     - 294
     - 211
     - **190**
     - 254
     - 377
     - 611
   * - ``lusgs``, 4 sweeps
     - —
     - —
     - 189
     - 253
     - —
     - —

``steady.line_order = diagonal`` and ``steady.line_alternate = 1`` give 190 iterations at a cap of 8, as the
sequential order. On the conventionally neutral deck on its own mesh (30 × 30 × 26, two ranks) the counts are
close: ``zjacobi`` with 4 sweeps 220 (cap 2) and 209 (cap 8), with 16 sweeps 208, ``lusgs`` 209 (cap 8) and 221
(cap 16).

*Findings.*

- ``lusgs`` with two sweeps at a cap of 8 needs the fewest iterations of the configurations with two sweeps or
  fewer per iteration. ``zjacobi`` matches it only with 16 sweeps per iteration.
- Every sweep has a cap above which iterations rise before the run becomes unstable: the best cap is below the
  stability limit. The plan's measure, the highest stable CFL, would choose a slower setting.
- ``zjacobi`` with its default two sweeps does not reach steady on Askervein at a cap of 4.

.. list-table:: Reference window :math:`K` (``steady.residual_ref_iterations``; cap 8; iterations in the hold,
   the CFL ramp and at the cap)
   :header-rows: 1
   :widths: 26 14 12 12 12 24

   * - Run
     - Iterations
     - Hold
     - Ramp
     - Cap
     - :math:`\sum` CFL to :math:`r < 10^{-2}`
   * - ``lusgs``, :math:`K` = 100
     - 190
     - 100
     - 13
     - 77
     - 357
   * - ``lusgs``, :math:`K` = 50
     - 156
     - 50
     - 26
     - 80
     - 358
   * - ``lusgs``, :math:`K` = 20
     - 141
     - 20
     - 41
     - 80
     - 356
   * - **lusgs, K = 5 (default)**
     - **137**
     - 5
     - 52
     - 80
     - 356
   * - ``zjacobi``, 4 sweeps, :math:`K` = 100
     - 246
     - 100
     - 49
     - 97
     - 725
   * - ``zjacobi``, 4 sweeps, :math:`K` = 5
     - 204
     - 5
     - 102
     - 97
     - 754

*Findings.*

- The residual falls with the accumulated pseudo-time :math:`\sum_n \mathrm{CFL}_n`, not with the iteration count:
  for every :math:`K` the relative residual reaches :math:`10^{-2}` at the same :math:`\sum` CFL within 2 %
  (``lusgs``) and 4 % (``zjacobi``). The hold runs at ``steady.cfl0`` and adds pseudo-time slowly.
- A shorter hold leaves a larger residual at its end, so the ramp (:math:`\mathrm{CFL} = \mathrm{cfl0}\,r^{-p}`)
  takes longer; the iterations at the cap do not change. From 100 to 50 the run saves 34 iterations, not 50.
- ``steady.cfl0`` scales the whole ramp: with :math:`K` = 5, ``lusgs`` needs 118 iterations with ``cfl0 = 2``
  and 192 with ``cfl0 = 0.5``.
- ``zjacobi`` needs twice the pseudo-time of ``lusgs`` to reach the same residual.

.. list-table:: Grid sequencing (cap 8; level-0 iterations, the coarse levels in brackets)
   :header-rows: 1
   :widths: 40 20 20 20

   * - Sweep
     - 1 level
     - 2 levels
     - 3 levels
   * - ``zjacobi``, 4 sweeps
     - 246
     - 211 (127)
     - 213 (100, 68)
   * - ``zjacobi``, 16 sweeps
     - 190
     - 111 (111)
     - 112 (100, 42)
   * - ``lusgs``, 2 sweeps
     - 190
     - 110 (111)
     - 112 (100, 42)
   * - **lusgs, 2 sweeps, K = 5**
     - **137**
     - **110 (44)**
     - —

*Findings.*

- Two levels cut the level-0 iterations of ``lusgs`` by 42 %; each coarse iteration costs about a quarter of a
  fine one. A third level adds nothing on level 0.
- With :math:`K` = 5 the coarse level needs 44 iterations instead of 111; the level-0 count is the same.

.. list-table:: Box count (PT-2.6.3; ``lusgs``, two sweeps, cap 8; iterations)
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Boxes
     - ``lagged``
     - ``zero``, overlap 0
     - ``zero``, overlap 2
     - ``zero``, overlap 4
   * - 1
     - 189
     - 189
     - —
     - —
   * - 4
     - 190
     - 191
     - 190
     - 190
   * - 16
     - 190
     - 196
     - 193
     - 191
   * - 64
     - 191
     - 204
     - 199
     - 195

*Findings.*

- With the lagged ghosts (the default) the box count changes the iterations by at most 1 %.
- Without communication (``zero``) 64 boxes cost 8 % more iterations; an overlap of 4 columns recovers half of it.
  The overlap default stays 0 with lagged ghosts; on a GPU the choice depends on the cost of the exchange (PR-32).

**Comparison with ERF (RT-2.6.1).** ERF's Askervein deck at the pinned commit (``inputs_anel``: anelastic URANS,
the k-equation closure, fifth-order horizontal and third-order vertical upwind advection, the power-law inflow from
210° through the west and south faces) ran to 4800 s on the same mesh (10 378 steps, 8 ranks) and was averaged over
2400–4800 s (``ablsteady_refavg``; the two half-window averages differ by at most 7.6e-4 m/s). The 10 m speeds on
lines A and AA then changed by less than 1e-3 of their final values after about 1500 steps and 3e-4 after 1730.
ablsteady runs the same setup (``Tests/Real_Terrain/Askervein_210``):

.. list-table:: Askervein, wind from 210°, 10 m speeds against ERF's URANS average
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Configuration
     - Iterations
     - Hilltop A / AA
     - Upwind masts
     - Lee masts (d > 300 m)
   * - Phase 1 (``zjacobi``, 4 sweeps, cap 2, :math:`K` = 100)
     - 538
     - +3.0 % / +1.5 %
     - +0.9 to +1.5 %
     - −5.0 to −7.7 %
   * - **Phase 2 default (lusgs, cap 8, K = 5, 2 levels)**
     - **96 (42)**
     - +3.0 % / +1.5 %
     - +0.9 to +1.5 %
     - −5.0 to −7.7 %

*Findings.*

- The two configurations reach the same lines (to 0.005 %); the Phase 2 default needs 5.6 times fewer level-0
  iterations.
- Up to 300 m past the hilltop every mast is within 4.5 % of ERF; the four lee masts beyond are 5.0–7.7 % slower.
  ablsteady's second-order advection (MUSCL with van Leer) against ERF's fifth order sets the wake there; the
  same trend (second-order upwind gives lower lee speeds than a higher-order scheme) is reported for WindNinja by
  Wagenbrenner et al. (2019).
- A tighter verdict (six orders instead of four, 270 iterations) changes no line speed by more than 0.001 %.
- The upwind and hilltop excess comes from the inflow :math:`k`. ablsteady's inflow :math:`k` is the log-law
  equilibrium with the unbounded length :math:`\kappa(z + z_0)`, which grows with height on ERF's power-law
  profile (8 m²/s² at 316 m, against 0.31 in ERF). With the closure's capped length (``turb.max_geom_lscale``)
  in the inflow table the upwind masts agree within 0.6 % and the hilltops within 2 %; the lee masts move further
  (to −9.4 %). That change is not in this version.

Newton–Krylov (in development)
------------------------------

*Setup.* The Ekman column of UT-3.1.2 with plain Newton steps.

- The scaled residual falls by 0.403, 0.171, 0.0131, 2.09e-4 and 8.26e-8 in the first five steps, with
  :math:`r_{n+1}/r_n^2` between 0.25 and 0.57: quadratic convergence.
- On a 30 m hill under the wall model (9 × 6 × 8 cells, :math:`\Delta\tau` = 20 s), four steps reduce :math:`|R|`
  from 11.3 to 0.211.

GMRES iterations over ten Newton steps on the column (Krylov tolerance 1e-2):

.. list-table::
   :header-rows: 1
   :widths: 24 19 19 19 19

   * - Restart length
     - 10
     - 30 (default)
     - 60
     - 200
   * - GMRES iterations
     - 4942 (cap reached)
     - 1469
     - 647
     - **346**

The defaults are provisional and will be tuned together with the Newton stage.

Pseudo-transient continuation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The continuation (:doc:`solution_algorithm`) on the same hill (six uneven boxes), on the hill with
:math:`u` three times larger, and on the Ekman column, GMRES(200), up to :math:`|R| \le 10^{-8}|R_0|`
(``Continuation.DISABLED_ParameterSweep``; steps accepted, GMRES iterations in all):

.. list-table::
   :header-rows: 1

   * - Settings
     - Hill
     - Hill, :math:`u \times 3`
     - Column
   * - **Defaults, GMRES(200)** (:math:`c` 1e-4, :math:`\beta` 0.5, :math:`p` 1, CFL cap 10, infinite :math:`\Delta\tau` at the cap)
     - **8, 976**
     - **8, 967**
     - **17, 510** (11 backtracks)
   * - No line search (``steady.newton_lambda_min = 0.99``)
     - 8, 976
     - 8, 967
     - fails at :math:`|R|/|R_0|` = 3.0e-2
   * - ``steady.newton_dtau_inf_at_cap = 0``
     - 17, 1790
     - 16, 1395
     - :math:`|R|/|R_0|` = 9.1e-4 after 60
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
     - fails at :math:`|R|/|R_0|` = 3.5e-2
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
   * - ``steady.cfl0 = 10`` (CFL 10, then Newton)
     - 6, 940
     - 13, 2585
     - 11, 370
   * - ``steady.newton_ptc_growth = 1.5``
     - 8, 976
     - 8, 967
     - 17, 510
   * - **steady.newton_ptc_growth = 1.5, p = 0.5**
     - **9, 1126**
     - **9, 1088**
     - **26, 475**
   * - ``steady.newton_ptc_growth = 1.5``, no line search
     - 8, 976
     - 8, 967
     - :math:`|R|/|R_0|` = 2.8e-3 after 60
   * - Leave infinite :math:`\Delta\tau` after a damped or unconverged step (tried, not kept)
     - 8, 976
     - 8, 967
     - fails at :math:`|R|/|R_0|` = 8.8e-5
   * - Plain Newton (``NewtonKrylov::step``, no globalization)
     - 6 steps
     - diverges (:math:`|R|` 99.5, 377, 12 971, ...)
     - 9 steps, :math:`|R|` rises once

:math:`|R|` with the defaults: hill 11.3, 2.18, 0.842, 0.282, 0.124, 0.0403, 4.0e-4, 4.0e-6, 3.9e-8 (CFL 1, 2, 4,
8, then infinite); hill :math:`u \times 3` 99.5, 42.1, 12.3, 2.42, 0.928, 0.0734, 2.0e-3, 1.9e-5, 1.9e-7; column
0.760, 0.568, 0.387, 0.223, 0.0912, 0.0241, 0.0232, 0.0178, 0.0167, 1.1e-3, 1.06e-3, 9.95e-4, 9.61e-4,
9.38e-4, 4.5e-4, 4.5e-6, 5.2e-7, 5.0e-9 (eleven backtracks, the slow steps at :math:`\lambda < 1`). No setting showed limit cycling (a CFL alternating between cut and growth); the
failures are retries running out at a small CFL. GMRES(30) on the hill: 23 steps to 6.6e-9 without the switch to
infinite :math:`\Delta\tau`; with it (the default restart of 30 with the default switch), twelve steps at the 2000-iteration limit near 6e-3 and a failure at 7.3e-8;
leaving infinite :math:`\Delta\tau` after a damped or unconverged step fails earlier, at 2.6e-5 after 11 steps.
``steady.newton_ptc_growth = 1.5`` changes no result of the defaults above and turns the two column failures at a
finite CFL into convergence (:math:`p` = 0.5: 26 steps) or slow progress (no line search). With :math:`p` = 0.5 it
converges on all three cases; it is the recommended setting until the defaults are tuned.

In UT-3.2.1 (one box, GMRES(30) stopped at 200 iterations, :math:`u \times 3`, plain Newton throughout)
``NewtonKrylov::step`` raises :math:`|R|` 3.56-fold on the first step, where :math:`k` goes negative, although its
Krylov relative residual is 0.013: the rise is nonlinear. The continuation keeps :math:`k` at its floor; on macOS
arm64 it takes full steps 0–2 and 4 and backtracks on steps 3 (:math:`\lambda = 1/2`) and 5 (:math:`\lambda = 1/4`).
These solves stop at their iteration limit, so the sequence depends on the platform's rounding. With the former
2000-iteration limit, Linux GCC took three full steps (reductions 0.670, 0.783, 0.826) and then three at
:math:`\lambda = 1/2`, where the full step would raise :math:`|R|` by 2.77. The test checks the properties (Armijo,
monotone :math:`|R|`, the :math:`k` floor, at least one backtrack), not the sequence.

The Newton preconditioner
^^^^^^^^^^^^^^^^^^^^^^^^^

*Setup.* UT-3.3.3: the open hill (9 × 6 × 8 cells, six boxes, a 30 m hill, inflow and outflow, the wall law), four
Newton steps at a fixed CFL, Krylov tolerance 1e-3, GMRES(200). Mean GMRES iterations per Newton step:

.. list-table::
   :header-rows: 1

   * - Preconditioner
     - CFL 10
     - CFL 1e2
     - CFL 1e3
     - CFL 1e4
   * - ``identity``
     - 163.5
     - —
     - —
     - 227.5
   * - **simplec (default), zjacobi**
     - **10.25**
     - **10.5**
     - **10.5**
     - **10.5**
   * - simplec with the stage's denominator
     - 18
     - —
     - —
     - 360–381

UT-3.3.3 runs ``identity`` and the stage's denominator only at CFL 10 and 1e4 (—: not run). In the continuation
loop on the same hill (UT-3.3 (a)) both preconditioners converge in 8 steps; ``simplec`` needs 54 GMRES iterations
in all and ``identity`` 976. The stage's denominator varies between runs at large CFL (359 and 381 measured at 1e4)
because GMRES at the 200-vector restart is sensitive to rounding; the bound fails either way.
