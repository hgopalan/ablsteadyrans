SIMPLEC stage
=============

The SIMPLEC stage (sub-phase 1.15, ``Source/Simplec.{H,cpp}``) is the first solver the application runs: one
pseudo-time iteration of the segregated pressure–velocity–scalar solver, handed to the driver (:doc:`driver`) as a
``Stage``. ``steady.stage = simplec`` turns it on; the default ``none`` keeps the write-and-stop behaviour of the
earlier versions.

One iteration
-------------

With the cell-centred local step :math:`\Delta\tau` from the driver:

#. the state's ghost cells and the wall state at :math:`q` (:doc:`boundary`, :doc:`surface`);
#. the residual :math:`R(q)` with the closure's fields at :math:`q` (:doc:`residual`, :doc:`closure`); its scaled norms
   are what the stage reports;
#. the momentum predictor :math:`(1/\Delta\tau + A_1)\,\delta u = -R` (:doc:`predictor`): :math:`u^*`;
#. the pressure correction (below): :math:`p' = p' + \alpha_p\,\phi`, :math:`u = u^* - (1/d)\,\nabla\phi`;
#. :math:`\theta` and :math:`k`: one z-line Jacobi step each in delta form, :math:`(1/\Delta\tau + A_1)\,\delta q =
   -R_q`, with :math:`A_1` from ``implicit_rows`` (the wall-k rows and the point-implicit k sink
   :math:`\partial S_k/\partial k` on the diagonal, :doc:`closure`); :math:`k` is then clipped at ``turb.tke_floor``.
   The convection of :math:`\theta` is re-evaluated with the corrected velocity's transport (the second half of the
   point-implicit buoyancy below); its diffusion, the residual of :math:`k` and the closure's sources are those of the
   start of the iteration, so the iteration costs one residual evaluation;
#. the closure accepts its fields as the relaxation history (``turb.nut_relax``, :doc:`closure`).

Pressure correction
-------------------

The predicted velocity :math:`u^*` is not divergence-free. With the SIMPLEC denominator :math:`d = a_P - \sum a_{nb}`
the predictor stores on every face (:math:`1/\Delta\tau` plus the folds, :doc:`predictor`), the correction solves

.. math::
   :label: simplec_poisson

   \nabla\cdot\Bigl(\frac{1}{d}\,\nabla\phi\Bigr) = \nabla\cdot u^*, \qquad
   u = u^* - \frac{1}{d}\,\nabla\phi, \qquad p' = p' + \alpha_p\,\phi,

with the Poisson layer's face coefficients :math:`\beta = 1/d` (:doc:`poisson`; the flux :math:`-(1/d)\nabla\phi`
is the solver's coefficient-weighted flux, so the corrected divergence is exactly the solver's residual) and the
boundary types of the run's faces (``poisson_bc_from_faces``: Dirichlet on an outflow, Neumann elsewhere). The
stiffness of the point-implicit buoyancy (below) enters :math:`d` and makes :math:`\beta` anisotropic inside an
inversion (up to :math:`N^2\Delta\tau^2` between the vertical and horizontal faces: 4 at :math:`\Delta\tau` = 64 s);
the FFT preconditioner assumes an isotropic ratio, so the GMRES count grows with the CFL (5–10 iterations at the
regression deck's cap). Faces with
:math:`d = 0` (the Dirichlet rows of the open faces, the ground and top faces of :math:`w`) take :math:`\beta = 1`
and no correction. The correction is applied to :math:`(u, v, \Omega)` and :math:`w` follows (``w_from_omega``), as
the initializer does. In the singular case (no Dirichlet face) the solver removes the mean of the right-hand side: a
net source the boundaries cannot carry, refused above round-off (:math:`10^{-8}` of :math:`U/h`).

**Why SIMPLEC.** SIMPLE drops the neighbours' corrections, :math:`u' = -\nabla p'/a_P`, and must under-relax
:math:`p'`; SIMPLEC takes them equal to the cell's own, :math:`u' = -\nabla p'/(a_P - \sum a_{nb})`, so
:math:`\alpha_p = 1` (Van Doormaal and Raithby 1984). With the pseudo-time term :math:`d \to 1/\Delta\tau` away from
the folds and as :math:`\Delta\tau \to \infty` the correction tends to the plain projection (:math:`\phi \propto
1/\Delta\tau`, the velocity correction finite: UT-1.15.4).

**The tolerance.** The plan asked for a loose relative tolerance of about 0.1. With a stable stratification that is
unstable: the divergence the solver leaves advects :math:`\theta` vertically, the buoyancy of that error drives
:math:`w`, and the loop grows. On ERF's Neutral_ABL_Flat (an inversion from 700 m) the run diverges within ten
iterations at 0.1 and 1e-2 and converges at 1e-3 and 1e-4, so the default ``simplec.pressure_reltol`` is 1e-3
(five to ten GMRES iterations with the FFT preconditioner instead of two).

Point-implicit buoyancy and rotation
------------------------------------

Two conservative couplings make the explicit pseudo-time iteration oscillate or grow where turbulence does not damp
them:

- :math:`w` and :math:`\theta` form the oscillator :math:`\partial_\tau w = (g/\theta_0)\,\theta'`,
  :math:`\partial_\tau\theta' = -w\,\partial_z\theta` of frequency :math:`N`; an explicit step amplifies it by
  :math:`(1 + N^2\Delta\tau^2)^{1/2}` whatever the CFL, and in an inversion :math:`\nu_t \to 0`;
- :math:`u` and :math:`v` form the rotation :math:`\partial_\tau u = f(v - v_g)`, :math:`\partial_\tau v = -f(u - u_g)`,
  undamped above the boundary layer.

Both get the linearised backward Euler coupling in the predictor (:doc:`predictor`): the diagonal of the row gains
the stiffness (:math:`(g/\theta_0)\max(\partial_z\theta, 0)\,\tau_\theta` on :math:`w`; :math:`f^2\tau_c` on :math:`u`
and :math:`v`), into the SIMPLEC denominator :math:`d` as well, and the right-hand side gains the force of the change
the other equation's own residual predicts (:math:`-(g/\theta_0)\,\tau_\theta\,R_\theta`; :math:`\mp f\,\tau_c
R_{v,u}`). Here :math:`\tau_\theta = 1/(1/\Delta\tau + a_{P,\theta})` is the step :math:`\theta`'s own implicit rows
allow (:math:`\Delta\tau` alone predicts far too large a change where the diffusion of the inversion's kink is
stiff: the 2 × 2 analysis of the review gives :math:`|G| = 0.32` at :math:`N\tau = 3` against 0.91 with
:math:`\Delta\tau`, and with that loss of damping the first iterations of Neutral_ABL_Flat diverged), and
:math:`\tau_c = 1/\max(|f|, a_P - \sum a_{nb})` (the net implicit diagonal: the drag and damping folds; the
transport weights cancel against their neighbours): with :math:`\tau_c = \Delta\tau` the rotation damps by
:math:`(f\Delta\tau)^2/2 \approx 10^{-4}` per iteration at the horizontal CFL, and the inertial oscillation of
Neutral_ABL_Flat persisted for 20 000 iterations; one inertial time where friction is weak damps it by about
:math:`f\Delta\tau` per iteration, and the net diagonal where a drag dominates (the wall-drag fold of the first
cell) keeps the cross term from slowing that relaxation (the review found 0.9 per iteration instead of
:math:`1/(1 + a\Delta\tau)` with a plain :math:`1/|f|`; the full :math:`a_P` would not do, since its advection
weights are transport and brought the inertial oscillation back). Every added term vanishes at a fixed
point (:math:`R_\theta = R_u = R_v = 0` and the delta form), so the steady solution is unchanged; the terms are a
preconditioner of the pseudo-time path. UT-1.15.7 pins the formulas.

What the stage reports
----------------------

The six residual norms are RMS values per unit volume of :math:`R_u, R_v, R_w` (scaled by :math:`U^2/H`), of the
continuity residual :math:`\nabla\cdot u` at the start of the iteration (:math:`U/H`) and of :math:`R_\theta`,
:math:`R_k` (:math:`U\theta_s/H`, :math:`U k_s/H`), with :math:`U, \theta_s, k_s` the ``jvp_scales`` of the residual
inputs (10 m/s, 10 K, 0.1 m²/s²) and :math:`H` the domain height. The monitors are the largest face speed, the mean
:math:`k` and, with a surface layer, the mean :math:`u_*`. The application prints, for the first ten iterations and
then every tenth, the CFL, the relative residual, the six norms, the monitors, the Poisson iteration count, the
divergence before and after the correction and the number of cells clipped at the k floor.

The application loop
--------------------

With ``steady.stage = simplec`` the application reads ``bc.*``, ``solver.poisson_*``, ``turb.*``, ``forcing.*``,
``convection.*``, ``simplec.*`` and ``steady.*`` (:doc:`inputs`), builds the first guess (``init.first_guess =
profile``; the inflow profile is used on periodic domains too) or reads the checkpoint, runs the driver to its verdict
or ``steady.max_iterations``, prints the stage account, and writes the plotfile, the turbulence diagnostics plotfile
(``io.turb_plot_file``, :doc:`output`) and the checkpoint with the driver history and the closure's :math:`\nu_t`
history. A restart continues bit for bit (``rt_steady_column_restart_*``): the stage keeps its state in :math:`q` and
the closure's history, and its ``init`` restores the diffusivity the driver's step limit reads (one iteration behind
the state, as in the run it continues). A diverged run aborts after writing its output.

Results
-------

**The 1-D neutral column (RT-1.15.1).** A doubly periodic 4 × 4 × 64 column of 1000 m, surface layer (z0 = 0.1 m),
slip top, uniform :math:`\theta`, driven by the force :math:`F = u_*^2/H` with :math:`u_* = 0.35` m/s, from a log
law with :math:`u_* = 0.25` m/s; CFL cap 50 (nothing varies horizontally, so the lagged neighbours are exact).

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Quantity
     - Value
     - Check
   * - Iterations to the steady verdict (6 orders, 1e-8)
     - 2231 (2.8 s on one core)
     - —
   * - Mean :math:`u_*` against :math:`\sqrt{FH} = 0.35`
     - 0.3499999
     - 1e-4 (force balance)
   * - :math:`|U|` of the first cell against the wall law
     - exact
     - 1e-6
   * - :math:`|U|` of the second cell against the log law
     - 0.9 % low
     - 1 %
   * - Cells 2–7 against the log law
     - 1.3–1.9 % low
     - 2 %
   * - :math:`k` of the first cell against :math:`u_*^2/C_{\mu 0}^2`
     - exact
     - 1e-6

The profile sits below the log law above the second cell because the force-driven column has the stress
:math:`u_*^2(1 - z/H)`, not a constant, and because the discrete strain overshoots the log law's in the lowest cells
(RT-1.9.1); the local-equilibrium closed form for the linear stress is 1.0 % above the result in cells 2–7.

**The lid-driven cavity (UT-1.15.3).** Re = 100 on 64 × 64 cells (no-slip walls on x, the ground and the lid; two
periodic cells in y), run by the driver to its verdict in 1760 iterations (14 s): :math:`u` on the vertical centreline
within 0.0045 and :math:`w` on the horizontal one within 0.0083 of Ghia, Ghia and Shin's (1982) tabulated values
(the criterion is 0.01, 1 % of the lid speed).

**Neutral_ABL_Flat (RT-1.15.2).** ERF's deck (8 × 8 × 64, geostrophic wind 10 m/s at latitude 90°, :math:`f =
10^{-4}` 1/s, inversion from 700 m, slip top); four sweeps, CFL cap 2, pressure tolerance 1e-3. The findings above
(the tolerance, the point-implicit buoyancy and rotation) come from this case. The capped boundary layer entrains
the inversion and warms at 1.4e-7 K/s (0.012 K/day) without end (the mixed layer's :math:`R_\theta` is uniform; with
gravity off every residual but :math:`\theta`'s falls to round-off), so :math:`\theta` and :math:`w` keep a residual
2.4–2.6 orders below their reference while every residual is below ``steady.residual_tol``: the deck asks for a drop
of 2 orders. The regression run reaches the verdict in about 2850 iterations (13 s on two ranks); the mean
:math:`u_*` is 0.377 m/s and ERF's ``check_neutral.py --physics`` passes every check (:doc:`testing`).

**The bursts and their cause.** With the rotation's step at a full inertial time (:math:`\tau_c = 1/|f|`) some
runs burst: within ten iterations the u and v residuals grow fifteenfold, the SER ramp collapses the CFL to 0.01,
and the run crawls towards a state with :math:`k \approx 5` m²/s² between 500 and 900 m and :math:`\theta` mixed to
800 m (the shear across the inversion erodes it). Dumps of the fields around a burst show a horizontally varying
:math:`w` of alternating sign through the whole inviscid part of the inversion (:math:`k` at its floor,
:math:`\nu_t \approx 5\times 10^{-6}` m²/s) growing by 1.4 per iteration with :math:`\theta` still uniform: an
internal gravity wave fed by the rotation's right-hand-side term, which carries the *whole* residual of the other
component, :math:`f\tau_c R_v`, pressure and wave parts included; in the inviscid layer the net diagonal vanishes,
:math:`\tau_c` reaches :math:`10^4` s and the injection is at full strength. Whether the wave is seeded and the layer
thick enough was decided by round-off (the clang and gcc builds separate at iteration 90, where :math:`k` first
reaches its floor), the bursts came at CFL 0.8 as well as 2, more sweeps only delayed them, and the smooth limiter,
a tighter pressure tolerance, a larger k floor and the :math:`\nu_t` relaxation did not remove them. Bounding the step
to :math:`r` inertial times (``simplec.rotation_relax``, 0.3) removes them: with :math:`r` = 0.1, 0.2, 0.3, 0.5 the
bursting path reaches the verdict in 7262, 5655, 2705 and 2613 iterations, with :math:`r = 1` it bursts at
iteration 450, and with the default every path tried (two compilers, one to three ranks) converges in 2550–2860
iterations (0.2: 5600–5900). The deck keeps the k floor of 1e-6 so that round-off at the machine-epsilon floor does
not pick the path.

Terrain cases (sub-phase 1.16)
------------------------------

**The lagged sweeps over terrain.** ERF's 2-D hill (``Tests/Decks/inputs_steady_hill2d``) with the flat deck's
settings (four sweeps, CFL cap 2) brings the residual down to 0.1–0.2 of its reference in 2000–4000 iterations and
then grows slowly (an e-folding of about 400 iterations, in the :math:`w` residual first, a tilted pattern of three
wavelengths across the periodic domain). The growth is there with gravity off, with zero lapse rate, with Coriolis
off, at CFL 1 (slower), with ``rotation_relax`` 0.3, with a stronger or deeper Rayleigh layer, with all five fields
damped, with :math:`\alpha_p = 0.5`, with ``turb.nut_relax`` 0.3, and on one rank as on two (bit for bit); a 1 m and
a 1 cm hill grow the same way, so any horizontal variation seeds it (the flat deck never does: nothing varies
horizontally). It is gone with first-order upwind convection (no defect correction) and gone with 8, 16 or 64
sweeps at the same second-order residual: the lagged x neighbours of the z-line Jacobi sweeps are the carrier, with
the thresholds of :doc:`predictor` (CFL 2.4 at four sweeps, 3.7 at eight, 6.1 at sixteen), passed locally near the
crest at a nominal cap of 2. ERF's anelastic solver never splits the operator (an explicit RK3 at Courant 0.375
with fifth-order upwind, inside its stability region), which is why it never sees this. OpenFOAM's SIMPLEC converges
the predictor's linear system before the pressure step, which the lagged sweeps approximate; 16 sweeps and 64 give
the same path here. The terrain decks run 16 sweeps with a cap of 4 (the fastest of the converged variants) and ``simplec.rotation_relax``
= 1: with the planar-mean split (:doc:`predictor`) the default 0.3 converges too, in twice the iterations (table).

**ERF's 2-D hill (RT-1.16.1).** 128 × 1 × 64 cells, the Witch of Agnesi (100 m, half-width 500 m), the flat
deck's sounding and forcing, Rayleigh damping of :math:`w` and :math:`\theta` from 700 m (0.2 1/s). The verdict and
ERF's ``check_hill.py --physics``, every check passing: crest speed-up 0.59, 0.50 and 0.43 in the three lowest cells
(Jackson and Hunt's :math:`2h/L` = 0.4, the window 0.2–0.8), upstream :math:`u_*` 0.354 m/s, the upstream log law
within 0.1, 1.7 and 4.2 % in the three lowest cells. Against ERF 92410d57's own run of the deck on the same
machine: at 6 h (ERF's physics run) the crest speed-up profiles agree within 0.03 and :math:`w` within 0.05 m/s, and
:math:`u` differs by up to 1.2 m/s in the 400–700 m band, where ERF is a third of an inertial period into its spin-up
(its :math:`u` still changes by 0.65 m/s between 8 and 12 h); at 24 h ERF is within 0.33 m/s (RMS 0.12) of the steady
state while still changing by 0.28 m/s per 6 h.

.. list-table:: Sweeps and CFL cap on the 2-D hill (two ranks, M-series Mac, ``rotation_relax`` 1 unless noted; the
   machine shared with other runs except where marked idle)
   :header-rows: 1
   :widths: 22 18 18 18

   * - sweeps / cap
     - iterations
     - wall (s)
     - verdict
   * - 4 / 2
     - 20 000
     - 119
     - not steady (grows from 2400)
   * - 8 / 2
     - 8124
     - 97
     - steady
   * - 16 / 2
     - 7847
     - 112
     - steady
   * - 8 / 4
     - 5173
     - 71
     - steady
   * - **16 / 4**
     - **4590**
     - **77 (49 idle)**
     - **steady** (the deck)
   * - 16 / 4, ``rotation_relax`` 0.3
     - 8577
     - 93 (idle)
     - steady

**ERF's 3-D hill (RT-1.16.2).** 64 × 64 × 20 cells, the hill of revolution (100 m, 500 m), Rayleigh damping of
:math:`w` and :math:`\theta` from 550 m; 16 sweeps, cap 4, ``rotation_relax`` 1: the verdict after 2092 iterations
(148 s on four ranks, idle machine; 4538 iterations at ``rotation_relax`` 0.3, 8475 at a cap of 2) and ERF's
``check_hill3d.py --physics`` passing
every check (crest speed-up 0.39, 0.27, 0.21 against :math:`1.6h/L` = 0.32 within a factor of two, upstream
:math:`u_*` 0.367 m/s, the log law within 0.0, 6.6 and 12 %). ERF's 4 h run of its deck (9600 steps, 590 s on four
ranks) passes the same checks and differs from the steady state by up to 1.6 m/s in :math:`u` (RMS 1.2) and 0.1 m/s
in :math:`w`: four hours is a quarter of an inertial period.

**A smooth plateau (beyond the plan).** The same settings on a periodic plateau with tanh ramps
(``Tests/Decks/inputs_steady_ramp2d``, the terrain file from ``tools/make_ramp_terrain.py``: 100 m high, ramps of 250
m at ±750 m, the maximum slope 0.2 of the Witch of Agnesi, both codes reading the same file to 1e-10 m): the verdict
after 8635 iterations (95 s on two ranks, ``rotation_relax`` 0.3), ERF's 24 h run within 0.37 m/s of the steady
state in :math:`u` (RMS 0.14) and 0.13 m/s in :math:`w`, the crest speed-up profiles within 0.04 (0.52 against 0.51
in the first cell), as on the hill.

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

Variations from the plan
------------------------

- ``simplec.pressure_reltol`` defaults to 1e-3, not the plan's 0.1 (above).
- The point-implicit buoyancy and rotation of the predictor are not in the plan; without them the stratified and
  rotating cases do not converge.
- The scalars use the start-of-iteration residual with the convection of :math:`\theta` re-evaluated after the
  correction (the plan's "segregated sweeps" leave the ordering open).
- Lateral ``no_slip_wall`` faces and a wall velocity (``bc.<face>_wall_velocity``) are added for the cavity; the inflow
  profile serves as the first guess on periodic domains (``init.first_guess = profile`` without an inflow face).
- The constant pressure-gradient force ``forcing.pressure_gradient`` (ERF's ``abl_pressure_grad``) drives the 1-D
  column; the plan's column had no forcing specified.
- Sub-phase 1.16: the terrain decks run 16 z-line sweeps at a CFL cap of 4 (the plan's four sweeps grow over
  terrain, above); the ground and top Dirichlet value of :math:`w` is the no-flux value; the Rayleigh layer selects
  its fields and relaxes :math:`\theta` to the sounding at the absolute height; the rotation's point-implicit term
  splits the other residual into its horizontal mean and the fluctuation (:doc:`predictor`). The H100-against-CPU
  comparison (RT-1.16.4) and the Askervein line comparison are deferred to Kestrel.
