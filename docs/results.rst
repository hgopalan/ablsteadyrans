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

Newton–Krylov
-------------

On the Ekman column of UT-3.1.2, plain Newton reduces the scaled residual by 0.403, 0.171, 0.0131, 2.09e-4 and
9.41e-8 in its first five steps, with :math:`r_{n+1}/r_n^2` between 0.25 and 0.65 (quadratic convergence). At the
default Krylov tolerance 1e-2 every restart length from 10 up takes 3, 5 and 8 GMRES iterations.
