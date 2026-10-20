Solver
======

The steady equations :math:`R(q) = 0` are solved by a pseudo-time iteration. Each iteration takes a step

.. math::
   :label: pseudo-time

   \left(\frac{1}{\Delta\tau} + A\right)\delta q = -R(q), \qquad q \leftarrow q + \delta q,

where :math:`\Delta\tau` is a local pseudo-time step and :math:`A` approximates the Jacobian. A *stage* chooses
:math:`A` and how to solve for :math:`\delta q`. Phase 1 has one stage, SIMPLEC. The Newton–Krylov stage (Phase 3)
and FAS multigrid (Phase 4) plug into the same driver.

**Why pseudo-time.** The step is small where the flow changes fast and large where it is near steady. Because
:math:`R` is the full second-order residual, a converged state is the same whatever :math:`A` is.

When a run is steady
--------------------

A run is steady when the residuals and the monitored quantities pass together.

**Residuals.** With :math:`R_q^n` the scaled norm of equation :math:`q` at iteration :math:`n`,

.. math::

   R_q^n < \tau_R \quad\text{and}\quad R_q^n \le 10^{-N} R_q^\mathrm{ref}, \qquad
   R_q^\mathrm{ref} = \max\Bigl(\max_{m < K} R_q^m,\ \tau_R\Bigr) .

:math:`\tau_R`, :math:`N` and :math:`K` are ``steady.residual_tol`` (5e-5), ``steady.residual_drop_orders`` (4) and
``steady.residual_ref_iterations`` (5). The reference is the early peak, so the drop does not depend on how low the
first iteration happened to start.

**Monitors.** Each monitor :math:`M_j` (for example the hilltop speed-up or :math:`u_*`) must stay within
:math:`T_j = \tau_M\max(|M_j|, s_j)` over the last :math:`W` iterations, and Aitken's estimate of the remaining
distance to the limit must be within the same amount:

.. math::

   \max_{n-W \le m \le n}|M_j^m - M_j^n| \le T_j, \qquad
   \left|\frac{\rho\,d_2}{1 - \rho}\right| \le T_j, \qquad \rho = \frac{d_2}{d_1},\ -1 < \rho < 1,

where :math:`d_1, d_2` are the changes over the two halves of the window, :math:`\tau_M` is
``steady.monitor_rel_tol`` (1e-3) and :math:`W` is ``steady.monitor_window`` (50).

**Why the extrapolation.** A small change over the window does not bound the distance to the limit. For a slowly
converging monitor the window test alone can stop 10 % away from it; with Aitken the error at the stop is about
:math:`\tau_M`.

The verdict is ``Steady``, ``NotSteady`` (listing what fails) or ``Diverged`` (a NaN or infinity, final).
Oscillations with a period up to :math:`4W` are never declared steady.

Driver
------

**Local pseudo-time step.** In every cell,

.. math::

   \Delta\tau = \frac{\mathrm{CFL}}{\dfrac{|u_c|}{\Delta x} + \dfrac{|v_c|}{\Delta y} + \dfrac{|\Omega_c|}{J\,\Delta\zeta}
   + 2\nu\left(\dfrac{1}{\Delta x^2} + \dfrac{1}{\Delta y^2}\right)} .

Vertical diffusion is implicit in every stage, so it does not bound the step.

**SER ramp.** The CFL starts at ``steady.cfl0`` and grows as the residual falls (switched evolution relaxation):

.. math::

   \mathrm{CFL}_{n+1} = \min\bigl(\mathrm{CFL}_\mathrm{max},\ g\,\mathrm{CFL}_n,\ \mathrm{CFL}_0\,r_n^{-p}\bigr),

with :math:`r_n` the largest residual relative to its reference, :math:`p` = ``steady.ser_exponent`` and :math:`g` =
``steady.ser_max_growth``.

**Stage switches.** The driver moves to the next stage after ``steady.switch_window`` iterations below
``steady.switch_residual``, and back when the residual grows by ``steady.switch_back_growth`` over a window. No two
switches fall within one window.

**History and restart.** The driver records every iteration (stage, CFL, residuals, monitors, wall time,
allocation units) and writes it to the checkpoint. A restarted run reproduces the records and the verdict of an
uninterrupted one bit for bit.

SIMPLEC stage
-------------

``steady.stage = simplec``. One iteration:

#. fill the ghost cells and compute the wall state at :math:`q`;
#. evaluate :math:`R(q)` with the closure at :math:`q`;
#. momentum predictor: :math:`(1/\Delta\tau + A_1)\,\delta\mathbf u = -R_\mathbf u` gives :math:`\mathbf u^*`;
#. pressure correction: :math:`\mathbf u = \mathbf u^* - \nabla\phi/d`, :math:`p \leftarrow p + \alpha_p\phi`;
#. one z-line step each for :math:`\theta` and :math:`k` in the same delta form; clip :math:`k` at the floor.

:math:`A_1` is the first-order operator of :doc:`discretization`.

**Momentum predictor.** Each z column is solved exactly (Thomas algorithm), with the x and y neighbours lagged from
the previous sweep (z-line Jacobi). ``steady.zline_sweeps`` sets the sweeps per step. A Fourier scan
(``tools/sweeps_cfl_scan.py``) gives the largest stable CFL:

.. list-table::
   :header-rows: 1
   :widths: 25 15 15 15 15

   * - Sweeps
     - 2
     - 4
     - 8
     - **16**
   * - Largest stable CFL
     - 0.9
     - 2.4
     - 3.7
     - **6.1**

The flat ABL deck runs 4 sweeps at a cap of 2. The terrain decks run **16 sweeps at a cap of 4**.

**Line LU-SGS sweeps.** ``steady.sweep = lusgs`` (sub-phase 2.1; ``zjacobi`` is the default) keeps the exact z-line
solve but visits each box's columns in order: a forward pass uses updated lower x/y neighbours, then a backward pass
uses updated upper neighbours. One ``steady.zline_sweeps`` sweep is this symmetric pair, and one is allowed. For a
wind whose :math:`u` and :math:`v` have the same sign in the ordering, a pass resolves the first-order upwind
coupling exactly; the pure-advection defect-correction factor stays below one at every CFL, approaching one for the
odd–even mode as the CFL grows. With opposite signs, each pass leaves one upwind neighbour lagged. The first-order
line solve alone still contracts, by :math:`c_u c_v/((1+c_u)(1+c_v)) < 1` per sweep, close to one at a large CFL;
the defect-correction iteration (the second-order residual) is stable with one sweep only to a CFL of about 2.5,
and with two to about 33. Couplings across box faces use the previous pass (the
default ``steady.line_ghosts = lagged``, below), so the iteration path depends on ``mesh.max_grid_size`` and boxes
shorter than the CFL limit the solve. The hill results
are in :doc:`results`; ``tools/sweeps_cfl_scan.py --thresholds`` prints the scalar stability factors.

**Box-local sweeps, overlap and face ownership.** Sub-phase 2.2 adds two ways to treat the neighbours across a box
face (``Source/LineSweeps.H``):

- ``steady.line_ghosts = lagged`` (the default; sub-phase 2.1, bit for bit): ``delta.FillBoundary`` before every
  pass, so the neighbours across a box face or a periodic seam take the previous pass's values. Both copies of a
  shared face are solved; ``OverrideSync`` keeps the lowest box index's copy.
- ``steady.line_ghosts = zero``: no communication inside the sweeps. Each box sweeps the region it owns, grown by
  :math:`n_o` = ``steady.line_overlap`` columns in x and y, with :math:`\delta = 0` outside it, and keeps only the
  owned points. The region is clipped to the domain on a non-periodic axis and is not grown along an axis one cell
  wide. :math:`n_o` is at most 8 and at most the narrowest box width along each axis it applies to, so the grown
  region never reaches past the neighbouring box (or one period). The system (:math:`A`, :math:`a_P`, the
  right-hand side) on the grown region comes from one ``FillBoundary`` of copies with :math:`n_o` ghost cells. Along a periodic axis one cell wide (y in a 2-D deck)
  both neighbours are the column itself; their coefficients move to the diagonal, which solves that coupling
  exactly.

Each step is a restricted additive Schwarz step with the box-local LU-SGS as the subdomain solve. The outer
defect correction recomputes the residual, so the converged state does not depend on :math:`n_o` or the box layout;
the iterates do. Inside a box the error of a seam decays over about :math:`c` columns downstream (:math:`c` the
CFL); the overlap moves the zero boundary :math:`n_o` columns away from the owned points.

On the C grid (zero ghosts) each box owns the faces of its valid box in its nodal direction except its high face,
which is the next box's low face. The domain's high face belongs to the last box; on a periodic axis the domain's
low face is that face's image, so the seam also belongs to the last box. Each face-centred unknown is swept once per
pass. The non-owned copies take the owner's value through ``OverrideSync`` with this owner mask (AMReX's default
mask, the lowest box index, would keep a copy that was not swept). Cells and z faces belong to their box.

**GPU diagonal ordering.** ``steady.line_order = diagonal`` visits the planes :math:`(i - i_0) + (j - j_0) = s` in
increasing :math:`s` (forward) or decreasing (backward), one column per thread inside a plane. A column's lower x/y
neighbours lie on plane :math:`s - 1` and its upper neighbours on :math:`s + 1`, so each column sees the same
values as in the lexicographic order: the result is the sequential sweep's, bit for bit (UT-2.2.4), with either
``steady.line_ghosts``. Each plane is one kernel launch, so a box of :math:`n_x \times n_y` columns takes
:math:`n_x + n_y - 1` launches per pass instead of one. The GPU timings and the per-machine defaults come in a later
sub-phase.

**Alternating ordering** (beyond the plan of sub-phase 2.2). With ``steady.line_alternate = 1`` the odd sweeps of
each call (the second, fourth, ...) take the j-reversed pair, :math:`(+i, -j)` forward then :math:`(-i, +j)`
backward (with ``steady.line_order = diagonal`` the wavefronts are :math:`(i - i_0) + (j_1 - j) = s`, their own
kernels). Where :math:`u` and :math:`v` have opposite signs one of the two pairs follows the wind, so its forward
pass resolves the first-order upwind coupling exactly instead of lagging one neighbour (UT-2.2.6). It needs
``steady.zline_sweeps`` ≥ 2. Measured on four ranks (the truly neutral Hill_3D deck ``inputs_steady_hill3d_neutral_lusgs``:
zero ghosts, :math:`n_o = 2`, diagonal, two sweeps; iterations to round-off):

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - Case
     - off
     - on
     - ``zjacobi``
   * - :math:`(u_g, v_g) = (10, 0)`, cap 16
     - 1792
     - 1813
     - 3495
   * - :math:`(u_g, v_g) = (10, 0)`, cap 64
     - 2997
     - 2991
     -
   * - :math:`(u_g, v_g) = (10, -10)`, cap 16
     - not converged at 20 000
     - 3052
     - 6051
   * - :math:`(u_g, v_g) = (10, -10)`, cap 64
     - not converged at 20 000
     - 3585
     -
   * - Conventionally neutral Hill_3D (``inputs_steady_hill3d`` with the same LU-SGS settings, cap 16, two orders)
     - 1100
     - 1146
     -
   * - 2-D hill (``inputs_steady_hill2d_lusgs``, one cell in y)
     - 3062
     - 3061
     -

The ``zjacobi`` column is the deck's own 16 sweeps at a cap of 4. With an opposite-sign geostrophic wind the
ascending pairs stall (the TKE residual drops 0.04 orders in 20 000 iterations), while the alternating ordering
converges in half the ``zjacobi`` iterations. Where the wind keeps one sign it costs 0–4 % more iterations. The
default stays ``0`` (bit-identical); the per-machine defaults are chosen in a later sub-phase.

**Box count and overlap on CPU.** The conventionally neutral Hill_3D (``inputs_steady_hill3d``, 64 × 64 × 20,
periodic in x and y) with ``steady.sweep = lusgs``, two sweeps, a cap of 16 and the deck's two-order stop; 1 box on 1
rank, 4 and 16 boxes on 4 ranks; iterations to the verdict:

.. list-table::
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - Boxes
     - ``lagged``
     - ``zero``, :math:`n_o = 0`
     - ``zero``, :math:`n_o = 2`
     - ``zero``, :math:`n_o = 4`
   * - 1
     - 860
     - 998
     - 948
     - 926
   * - 4
     - 823
     - 1205
     - 1100
     - 1020
   * - 16
     - 974
     - 1552
     - 1310
     - 1270

On CPU the box-local sweeps take 16–59 % more iterations than the lagged ghosts at :math:`n_o = 0` and 8–30 % more
at :math:`n_o = 4`: δ = 0 at the seams (the periodic seam of a single box included) costs more than the overlap
recovers. What they save, the ``FillBoundary`` of every pass, matters on a GPU and on many ranks and is measured
with the per-machine defaults. More boxes cost iterations with either mode (823 on 4 boxes, 974 on 16, lagged). The
diagonal ordering gives the same iterates as the sequential one and the same wall time on CPU (16 boxes,
:math:`n_o = 2`: 1310 iterations, 57 s on four ranks).

**Pressure correction.** With the SIMPLEC denominator :math:`d = a_P - \sum a_{nb}`,

.. math::

   \nabla\cdot\Bigl(\frac{1}{d}\,\nabla\phi\Bigr) = \nabla\cdot\mathbf u^*, \qquad
   \mathbf u = \mathbf u^* - \frac{1}{d}\,\nabla\phi, \qquad p \leftarrow p + \alpha_p\,\phi .

**Why SIMPLEC.** SIMPLE drops the neighbours' corrections and must under-relax :math:`p`. SIMPLEC takes them equal
to the cell's own, so :math:`\alpha_p = 1` (Van Doormaal and Raithby 1984). ``simplec.pressure_reltol`` is 1e-3: a
looser tolerance lets the left-over divergence drive :math:`w` through buoyancy, and the run diverges.

**Point-implicit buoyancy and rotation.** Two couplings oscillate without damping in pseudo-time where :math:`\nu_t`
is small: :math:`w`–:math:`\theta` (frequency :math:`N`) and :math:`u`–:math:`v` (frequency :math:`f`). An explicit
step amplifies each by :math:`(1 + \omega^2\Delta\tau^2)^{1/2}`. The predictor treats both with a linearised backward
Euler step:

.. math::

   w:\ \ a_P \mathrel{+}= \frac{g}{\theta_0}\max(\partial_z\theta, 0)\,\tau_\theta,\quad
   \mathrm{rhs} \mathrel{-}= \frac{g}{\theta_0}\,\tau_\theta R_\theta ; \qquad
   u, v:\ \ a_P \mathrel{+}= f^2\max(\tau_c, \tau_\mathrm{own}),\quad
   \mathrm{rhs} \mathrel{\mp}= f\bigl(\tau_c\langle R\rangle + \tau_\mathrm{own}(\overline R - \langle R\rangle)\bigr).

- :math:`\tau_\theta = 1/(1/\Delta\tau + a_{P,\theta})` is the step :math:`\theta`'s own rows allow.
- :math:`\langle R\rangle` is the level mean of the other component's residual and :math:`\overline R` its local
  value. The mean takes :math:`\tau_c = 1/\max(|f|/r, a_P - \sum a_{nb})`, with :math:`r` =
  ``simplec.rotation_relax`` (0.3); the fluctuation takes the row's own step :math:`\tau_\mathrm{own}`.

**Why split the rotation.** The inertial oscillation is horizontally uniform, so only the level mean needs the
inertial time scale. Applying it to the local residual fed an internal gravity wave in the inviscid inversion.
Every added term vanishes when :math:`R = 0`, so the steady solution is unchanged.

**The w rows' VNN cap** (``simplec.w_vnn_cap``, off by default). On thin first cells the :math:`w` residual decays
slowly. The cause is the SIMPLEC denominator of the :math:`w` rows. The vertical neighbour weights
:math:`a_{ZM}, a_{ZP}` appear in both :math:`a_P` and :math:`\sum a_{nb}`, so they cancel and
:math:`d \approx 1/\Delta\tau`. The implicit z-line solve moves :math:`w` by about :math:`1/(1/\Delta\tau + K)` per
unit force, with :math:`K = a_{ZM} + a_{ZP}` (mostly the vertical diffusion). The pressure correction moves it by
:math:`\Delta\tau`. A correction that varies from layer to layer is therefore too large by the vertical VNN number

.. math::

   \mathrm{VNN} = \Delta\tau\,K,

and the next predictor takes most of it back. With a 0.4 m first cell under 43 m cells VNN reaches about 380, a few
metres above the ground. A linear model of one SIMPLEC iteration (one horizontal Fourier mode, a stretched log-law
column, an exact z-line solve) gives a slowest decay of :math:`1 - \rho \approx 1/(2\,\mathrm{VNN}_\max)`, which
matches the solver's late :math:`w` tail. With the cap :math:`V > 0` the :math:`w` rows take

.. math::

   \frac{1}{\Delta\tau} \;\rightarrow\; \frac{1}{\Delta\tau} + \frac{a_{ZM} + a_{ZP}}{V}
   \quad\text{in } a_P \text{ and in } d,
   \qquad \mathrm{VNN}_\mathrm{eff} = \frac{\mathrm{VNN}}{1 + \mathrm{VNN}/V} < V,

except on the Dirichlet ground and top rows. u, v, θ and k keep their rows.

- **Why the answer does not change.** The predictor is in delta form, so its fixed point is :math:`R = 0` whatever
  the diagonal. On the staggered grid there is no Rhie–Chow interpolation, so :math:`d` enters only the pressure
  correction: it changes the path to the fixed point, not the fixed point.
- **What it costs.** One kernel over the :math:`w` faces per iteration; nothing with the cap off (0 is
  bit-identical to the solver without the input).
- **Why both** :math:`a_P` **and** :math:`d`. In the same model, capping :math:`d` alone diverges for
  :math:`V \le 8`, and capping :math:`a_P` alone gains nothing.
- **Why the floor of 16.** On flat ground the model diverges for :math:`V \le 2` at CFL 16 and above, and for
  :math:`V = 4` at CFL 64. On the Askervein hill with a 1 m first cell :math:`V = 4` and 8 burst after the jump to
  CFL 16, while 16, 32 and 64 all converged in 200 iterations.
- **Terrain: a known limit.** The pressure correction's z-face flux scales the whole :math:`\mathbf n\cdot\nabla
  p'`, slope cross term included, by :math:`\beta_z = 1/d_w`. With the cap, :math:`\beta_z \ll \beta_x` near the
  ground, so the corrected :math:`w` gains :math:`-h(\beta_x - \beta_z)\,\partial_x p'` (:math:`h` the slope),
  and the operator becomes indefinite once :math:`\beta_z/\beta_x` falls below about :math:`h^2`. The model is
  stable roughly while :math:`1 + \mathrm{VNN}/V < 1/h^2`. On the Askervein hill with a 0.4 m first cell,
  :math:`V = 16` burst at the jump to CFL 16, which the uncapped run makes cleanly. With a 1 m first cell
  (max aspect ratio 50.78, slopes to 23°) :math:`V \ge 16` converged and :math:`V = 4` and 8 burst. The application
  therefore warns when the cap is on over terrain and the mesh's max aspect ratio, max(dx, dy)/thickness, is above
  the tested 50.78. The aspect ratio is a proxy: the condition depends on VNN, :math:`V` and the local slope.
  A larger cap is safer: on the 0.4 m Askervein hill (max aspect ratio 127) :math:`V = 64` converged in 278
  iterations where 16 burst. Without the cap the :math:`w` residual converges slowly on such meshes (the uncapped
  0.4 m hill needed about 2800 iterations). The warning does not stop the run. The slowest decay of a long horizontal wave (:math:`\lambda
  \gtrsim` 3 km) is not the :math:`w` mode and the cap barely changes it.

Measured on Askervein (140 × 140 cells at 43 m, the layers stretched by 1.08, line LU-SGS with two sweeps at a CFL
cap of 16, ``gmres_mlmg``), iterations to the steady verdict:

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
   * - Flat, 0.4 m
     - 1620
     - 194
     - not run
     - **196**
     - not run
     - not run
   * - Flat, 1 m
     - 485
     - not run
     - not run
     - **216**
     - not run
     - not run
   * - Flat, 1 m, CFL cap 32
     - 855
     - not run
     - not run
     - **285**
     - not run
     - not run
   * - Hill, 1 m
     - 324
     - burst
     - burst
     - **200**
     - 200
     - 201
   * - Hill, 0.4 m
     - about 2800
     - not run
     - not run
     - burst
     - not run
     - 278

"Burst": after the jump to the CFL cap the residuals rose to 0.2–0.8 and the CFL collapsed; those runs were stopped
unsteady (V = 4 at 870 iterations, V = 8 at 1470, V = 16 at 0.4 m at 460). With z-line Jacobi (four sweeps, CFL cap
4) the 1 m flat case needs 292 iterations either way: at that CFL the VNN is small and the cap has little to do. The
capped and uncapped steady states agree to the convergence tolerance: at 0.4 m with a tight verdict (6 orders,
5e-7; 353 iterations) u, v and w agree to 2.7e-7 m/s, and on the 1 m hill the 10 m speed-ups on lines A, AA and B
agree to 3e-6 (:doc:`testing`, ``rt_thin_cell_w_vnn_cap``). The pressure correction takes more Poisson iterations
with the cap (1 m hill: 2.83 per correction against 2.14), fewer in all (566 against 693).

Grid sequencing
---------------

With ``steady.sequencing = 1`` and ``mesh.hierarchy_levels`` = L ≥ 2 (sub-phase 2.5; ``Source/Sequencing.H``), a
fresh run solves on the coarse levels of the x/y hierarchy (:doc:`governing_curvilinear`) first and starts each finer
level from the coarser solution:

1. Level L − 1 (the coarsest): the initializer (``init.*``) on that level's mesh, then the SIMPLEC stage under its own
   driver to the *loose stop*.
2. Levels L − 2 … 1: ``prolong_state`` from the level above (bilinear per ζ-level, k clipped, the level's boundary
   conditions and a projection; :doc:`governing_curvilinear`), then SIMPLEC to the loose stop.
3. Level 0 (the run's mesh): ``prolong_state`` from level 1, then the run's driver to the strict stop
   (``steady.residual_*``, ``steady.monitor_*``), as without sequencing.

Every level has its own boundary conditions, closure (its :math:`\nu_t` history starts empty), wall distance,
pressure solver and driver, built from the run's inputs; the SER ramp starts at ``steady.cfl0`` on every level.
``init.*`` acts on the coarsest level only: every finer level starts from a prolonged, projected state, so
``init.projection = 0`` leaves only the coarsest start unprojected. ``steady.max_iterations`` caps each level
separately; the closing line ``steady after N iterations`` and the checkpoint count level 0's iterations, and a
summary table after it lists every level's iterations, wall time, allocation units and how it stopped. A coarse
level that reaches the cap is prolonged as it is, with a warning; one that diverges aborts the run.

**The residual reference.** The strict criterion asks each scaled residual :math:`R_q` to drop
``steady.residual_drop_orders`` orders below its reference :math:`R_q^{\mathrm{ref}}`, the largest of its first
``steady.residual_ref_iterations`` values. A prolonged start is already close to the solution, so a level's own
first residuals are small, and level 0 would be held to a lower absolute residual than a run without sequencing. With
sequencing every finer level starts its references at the coarsest level's (``Driver::seed_reference``), which, like
a run without sequencing, starts from the initializer; the first ``residual_ref_iterations`` values can still raise
them. So the drop is measured from the initial field's residual on every level. The relative residual
:math:`r_n` of the SER ramp uses the same reference: after the reference iterations it is already small on a finer
level, and the CFL then climbs to ``steady.cfl_max`` at ``steady.ser_max_growth`` per iteration (on
``Neutral_Hill_3D`` within 10 iterations of the hold, against about 200 without sequencing; part of the measured gain,
about 50–100 iterations, comes from this faster ramp). The coarsest level's references are not larger than
level 0's own in the cases measured (``Neutral_Hill_3D``: 0.66–1.0 of them, the coarse terrain being smoother), so the
seeding makes level 0's drop test slightly stricter there; a seed above level 0's own value would loosen the drop
test by at most :math:`\log_{10}(\text{seed}/\text{own})` orders, and the absolute and monitor tests not at all. The
summary table's last line prints both for every residual.

**The loose stop.** A coarse level stops after an iteration when

.. math::

   R_q^n \le \max\!\left(R_q^{\mathrm{ref}}, \tau\right) 10^{-d_\ell} \quad\text{or}\quad R_q^n < \tau
   \qquad \text{for every residual } q,

after at least ``steady.residual_ref_iterations`` iterations (both branches), with :math:`\tau` = ``steady.residual_tol`` and
:math:`d_\ell` = ``steady.sequencing_drop_orders`` (2), or when the strict criterion holds. With
``steady.sequencing_monitor_tol`` = :math:`\mu > 0` it also stops when every monitor :math:`m` changed by at most
:math:`\mu \max(|m|, s_m)` over the last ``steady.monitor_window`` iterations (the strict monitor test without its
extrapolation). That branch is off by default: right after a prolongation the monitors (largest face speed, mean k,
mean :math:`u_*`) are already near their limits and change by less than 1 % over a window while the residuals have
not dropped. On ``Neutral_Hill_3D`` with :math:`\mu = 0.01`, level 1 stopped after 51 iterations with no residual
drop and level 0 took 985 iterations, against 774 with the residual branch alone.

**Cost.** An iteration on level l costs about :math:`4^{-l}` of a level-0 iteration (the same z column, a quarter
of the columns per level), plus the Poisson solve. Measured on four ranks (Apple M-series, double precision):

.. list-table::
   :header-rows: 1

   * - Case (strict stop)
     - Levels
     - Coarse iterations (2, 1)
     - Level-0 iterations
     - Wall time (s)
   * - ``Neutral_Hill_3D`` (2 orders)
     - 1
     - —
     - 2091
     - 144.8
   * - ``Neutral_Hill_3D``
     - **3**
     - **383, 127**
     - **774**
     - **57.1**
   * - ``Neutral_Hill_3D``
     - 2
     - 692 (level 1)
     - 843
     - 69.4
   * - Truly neutral hill (round-off)
     - 1
     - —
     - 3495
     - 39.6
   * - Truly neutral hill
     - 3
     - 239, 101
     - 2773
     - 32.5
   * - Truly neutral hill, ``sequencing_drop_orders = 3``
     - 3
     - 323, 147
     - 2472
     - 28.8

The ``Neutral_Hill_3D`` deck stops at two orders (with its monitor test), the truly neutral deck
(``Truly_Neutral_Hill_3D/inputs_steady_hill3d_neutral``) at round-off (1e-11 and eight orders). The starting field
matters with sequencing (PT-2.5.2, ``Neutral_Hill_3D``, three levels): the plain projection (``init.projection = 1``,
:math:`\alpha_v = \alpha_h`, the default) needs 774 level-0 iterations, the first guess without projection
(``init.projection = 0``) 1442, and the mass-consistent field with ``init.alpha_v = 0.1`` 1007. Without sequencing
the three need 2091, 2042 and 2044. The default stays the plain projection.

The gain is the transient the coarse levels remove. Once the error is down to the slowest mode of the level-0
iteration, every run converges at the same rate, so a deep stop (round-off, 8 orders) gains less in proportion than
an engineering stop (2 orders). The converged state does not depend on sequencing: ``rt_steady_hill3d_neutral_sequencing``
matches the run without it to 1e-8 (:doc:`testing`).

**Restart and memory.** The coarse levels run on a fresh start only: a checkpoint holds the level-0 state and its
driver history (with the seeded references), and a restart with ``steady.sequencing = 1`` continues level 0 and says
so. The coarse levels' wall time and allocation units are in the log, not in the checkpoint. The level-0 solver
(stage, pressure solver, wall distance) is built before the coarse levels run and stays allocated while they do, and
the prolongation into level 0 builds a second level-0 pressure solver for its projection, so the peak memory is
above that of level 0 alone.

Poisson solver
--------------

The initializer, the pressure correction and the Tucker wall distance solve

.. math::

   \frac{1}{J}\Bigl[\delta_\xi\bigl(a_x\beta_x f_x\bigr) + \delta_\eta\bigl(a_y\beta_y f_y\bigr)
   + \delta_\zeta\bigl(a_z\beta_z f_z\bigr)\Bigr] = \mathrm{rhs},

where :math:`f` are ERF's face gradients with the terrain cross terms and :math:`\beta` are face coefficients
(:math:`1/d` for SIMPLEC, :math:`\Sigma` for the initializer). ``solver.poisson_backend`` selects:

- ``gmres_mlmg`` (default) and ``mlmg``: AMReX's ``MLTerrainPoisson`` with geometric multigrid, under GMRES or
  alone. It coarsens only in x and y and solves each column exactly, the usual treatment of cells much thinner
  than they are wide. It gives the same solution as back end A to 1e-9 and needs no FFT;
- ``gmres_fft``: ERF's ``TerrainPoisson`` with ``amrex::GMRES`` and ERF's FFT preconditioner (the flat
  stretched Laplacian, FFT in x and y, tridiagonal in z).

**Why gmres_mlmg is the default.** On the Askervein terrain at 43 m horizontal spacing, the two back ends give the
same SIMPLEC iterations to the steady verdict when started from the same field, but the column solves keep back end
B's work per pressure correction nearly flat as the first layer thins (mean Poisson iterations per correction):

.. list-table::
   :header-rows: 1
   :widths: 25 15 30 30

   * - First layer (m)
     - dx/dz
     - ``gmres_fft``
     - ``gmres_mlmg``
   * - 8
     - 5
     - 196 iterations, 6.6
     - 196 iterations, 2.6
   * - 4
     - 11
     - 203 iterations, 7.7
     - 203 iterations, 2.6
   * - 2
     - 21
     - 230 iterations, 9.6
     - 230 iterations, 2.5
   * - 1
     - 43
     - 324 iterations, 12.8
     - 324 iterations, 2.1

(140 × 140 cells, the layers stretched by 1.08 to about 750 m, line LU-SGS with two sweeps at a CFL cap of 16,
``init.alpha_v = init.alpha_h = 1``.) The tolerances differ: back end B stops on the :math:`J`-weighted residual
(the cell-volume, mass-imbalance norm), back end A on the plain one. On these Askervein runs, where the
divergence sits in the thin near-ground cells (:math:`J \ll 1`), back end B left 10–30 times less divergence after
each correction at the same ``simplec.pressure_reltol``; where the divergence sits in tall cells aloft
(:math:`J > 1`) the same value stops looser in the plain norm, by up to :math:`J_\max/J_\min`.

**Residuals that were exactly zero.** Back end A solves a horizontally uniform problem exactly, so some residuals
stayed exactly 0 (v across two periodic y cells, θ with no gravity). With back end B they sit at round-off
(about 1e-13 in double). ``steady.residual_tol`` times :math:`10^{-N}` (``steady.residual_drop_orders``) far below
that floor can then not be met; the defaults (5e-5, four orders) are far from it.

**The initializer's split coefficient.** With :math:`\alpha_v \ne \alpha_h` the initializer's operator is
:math:`\nabla\cdot(\mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2)\nabla\lambda)`, whose z-face flux needs a
separate coefficient on the cross part (``FaceCoeffs::z_cross``). Back end B has none, so that one projection runs
on back end A, and the log says so; the SIMPLEC stage keeps ``solver.poisson_backend``. With
:math:`\alpha_v = \alpha_h` the cross part takes the z coefficient and no split is used.

The solve stops at ``solver.poisson_reltol`` or ``solver.poisson_abstol``. With no Dirichlet face the mean of the
right-hand side is removed first. In single precision the rounding floor is about :math:`\varepsilon N_z^2/2`, so
use a tolerance near :math:`4\varepsilon N_z^2`.

Newton–Krylov stage
-------------------

Sub-phase 3.1 adds one Newton step solved by GMRES. The application does not run it yet (sub-phase 3.5).

.. math::

   (D + J)\,\delta q = -R(q), \qquad A = T^{-1}(D + J)\,S, \qquad b = -T^{-1}R(q), \qquad \delta q = S x .

- :math:`J = \partial R/\partial q`, applied by the Jacobian–vector product of :doc:`discretization`. With fixed
  diffusivities :math:`\nu_t` is frozen; with a closure the residual computes :math:`\nu_t(q)`, so :math:`J`
  differentiates it too.
- :math:`D` is :math:`1/\Delta\tau` on momentum, :math:`\theta` and :math:`k`, and 0 on continuity and Dirichlet
  rows. :math:`D = 0` is plain Newton.
- :math:`S` and :math:`T` scale the fields so GMRES weighs them alike.

GMRES is flexible, so the preconditioner (sub-phase 3.3) may change between iterations; it is the identity for now.

**Ghost fill.** :math:`R` here is :math:`R(\mathrm{fill}(q))`: the ghost cells come from the boundary conditions,
and :math:`Jv` refills :math:`q + \epsilon v`. On terrain with ``bc.zlo = surface_layer`` one fill is not a fixed
point: the surface layer reads the ground :math:`w` face, which the no-flux rule then resets from the
:math:`u, v` ghosts, so a second fill changes :math:`R` by about 1e-5. The step therefore evaluates :math:`R` on a
copy of :math:`q` filled once more, and every :math:`Jv` differences against that. On flat ground and on walls
the fill is idempotent and nothing changes.

**When an update is skipped.** If :math:`|b|` is not finite, or the Krylov residual comes back above :math:`|b|`
(or not finite), the update is not applied and the step prints a line saying so. A stagnated solve (a short
restart at its iteration cap) can end at :math:`|b|` plus round-off; its update is then about zero, so skipping
it does not change the answer. The closing ghost fill still runs,
so on terrain under the wall law the ground :math:`w` face can move slightly.

**Memory.** FGMRES(:math:`m`) keeps :math:`2m + 3` vectors of state size. The run stops, naming
``steady.newton_mem_budget_gb``, if they do not fit on the rank with the most points. The estimate leaves out the
step's own buffers (about 8 vectors more while GMRES runs, 2 more with a preconditioner) and the residual's
temporaries. On one 512 × 512 × 128 box restart 20 passes the check at 69.4 GB, but the Krylov vectors and the
step's buffers come to about 82 GB, before the state, the mesh and the residual's temporaries, so leave that
margin.

**Limits of GMRES here.** The forward difference has an error near 1e-8 in double (1.6e-3 in single). Below that
level GMRES still reports convergence (its Arnoldi estimate drops to round-off) while the true linear residual
stays at the floor: 1.4e-8 on the Ekman column, 4.6e-8 to 1.7e-6 on the 30 m hill case. A Krylov tolerance below
that buys nothing. Short restarts can stall: on the Ekman column restart 5 reaches the 2000-iteration cap on the
third step at the default 1e-2, and restarts 5 and 10 at 1e-6; on the hill case plain Newton with restart 30
stalls at the cap once :math:`|R|` is small (2.6e-10 at rtol 1e-2), while restart 200 converges to round-off.

**Defaults.** ``steady.newton_krylov_rtol = 1e-2`` and ``steady.newton_gmres_restart = 30`` are provisional
(sub-phase 3.6 tunes them). On the column they take 3, 5 and 8 iterations for the first three steps, but over ten
plain Newton steps restart 30 takes 1469 iterations against 346 for restart 200, and restart 10 reaches the cap on
step 4. A single step can also raise :math:`|R|` (by 25× at rtol 0.1 on the column); the line search and the
:math:`\Delta\tau` ramp of sub-phase 3.2 address that.

Reproducibility
---------------

Residual norms and monitors use ``deterministic_sum``: each column is summed in k order, columns in storage order,
and each box total goes into its own slot of a vector that is reduced over the ranks. The result has the same bits
for any number of ranks on a fixed BoxArray, so every rank reaches the same verdict at the same iteration. The
Newton stage's inner product and the finite-difference step of :math:`Jv` (sums over face and cell points with
shared faces counted once) are summed the same way, so the GMRES iterates are bit-identical on any rank count (with
the identity preconditioner, or a hook whose result does not depend on the distribution). A
different BoxArray adds in a different order: the GMRES solutions of two box layouts agree to the finite-difference
round-off carried through the solve's restarts, not to the last bit.

References
----------

- Van Doormaal, J. P., Raithby, G. D. (1984). *Numerical Heat Transfer* 7, 147–163.
- Mulder, W. A., van Leer, B. (1985). *J. Comput. Phys.* 59, 232–246.
