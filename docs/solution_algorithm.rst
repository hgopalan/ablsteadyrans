Solution algorithm
==================

The steady equations :math:`R(q) = 0` are nonlinear and strongly coupled, so they cannot be solved in one step.
ablsteady approaches the solution by a *pseudo-time* iteration: it marches an artificial time until nothing changes
any more. The march does not need to be time-accurate, so each cell may take its own, as-large-as-stable time step.
Only the end point matters.

Pseudo-time iteration
---------------------

Each iteration solves for a correction :math:`\delta q`:

.. math::
   :label: pseudo-time

   \left(\frac{1}{\Delta\tau} + A\right)\delta q = -R(q), \qquad q \leftarrow q + \delta q .

- :math:`R(q)`: the full second-order residual (:doc:`discretization`).
- :math:`\Delta\tau`: the local pseudo-time step.
- :math:`A`: an approximation of the Jacobian :math:`\partial R/\partial q`.

A *stage* chooses :math:`A` and the way :math:`\delta q` is computed. When :math:`\delta q = 0`, :math:`R(q) = 0`,
whatever :math:`A` is. A cheaper or rougher :math:`A` changes how fast the iteration converges, not the answer.

The production stage is SIMPLEC (``steady.stage = simplec``), optionally with grid sequencing (below). The
Newton–Krylov and multigrid stages are in development.

Local time step and CFL ramp
----------------------------

In every cell the pseudo-time step is set by a CFL number:

.. math::

   \Delta\tau = \frac{\mathrm{CFL}}{\dfrac{|u_c|}{\Delta x} + \dfrac{|v_c|}{\Delta y} + \dfrac{|\Omega_c|}{J\,\Delta\zeta}
   + 2\nu\left(\dfrac{1}{\Delta x^2} + \dfrac{1}{\Delta y^2}\right)} .

Vertical diffusion is treated implicitly, so it does not limit the step.

The CFL number starts small and grows as the residual falls (switched evolution relaxation, SER; Mulder and van Leer
1985):

.. math::

   \mathrm{CFL}_{n+1} = \min\bigl(\mathrm{CFL}_\mathrm{max},\ g\,\mathrm{CFL}_n,\ \mathrm{CFL}_0\,r_n^{-p}\bigr).

- :math:`\mathrm{CFL}_0`, :math:`\mathrm{CFL}_\mathrm{max}`: ``steady.cfl0`` and ``steady.cfl_max``.
- :math:`r_n`: the largest residual relative to its reference value (:ref:`convergence-criteria`).
- :math:`p`: ``steady.ser_exponent``; :math:`g`: ``steady.ser_max_growth``, the largest growth per iteration.

SIMPLEC stage
-------------

SIMPLEC (Van Doormaal and Raithby 1984) solves the momentum equations first, then corrects pressure and velocity so
that mass is conserved. One iteration:

#. Fill the ghost cells and evaluate the wall model at the current state :math:`q`.
#. Evaluate the residual :math:`R(q)` and the turbulence model.
#. **Momentum predictor.** Solve :math:`(1/\Delta\tau + A_1)\,\delta\mathbf u = -R_\mathbf u`, giving
   :math:`\mathbf u^*`.
#. **Pressure correction.** Solve a Poisson equation for :math:`\phi` and set
   :math:`\mathbf u = \mathbf u^* - \nabla\phi/d` and :math:`p \leftarrow p + \alpha_p\phi`.
#. **Scalars.** Take one step of the same form for :math:`\theta` (with its convection updated to the corrected
   velocity) and :math:`k`; bound :math:`k` below by its floor.

:math:`A_1` is the first-order operator of :doc:`discretization`.

Momentum predictor and line solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ABL meshes have cells much thinner than they are wide, so the vertical coupling is the strongest. The predictor
therefore solves each vertical column exactly (a tridiagonal solve with the Thomas algorithm). The coupling to the
horizontal neighbours is handled by sweeps; ``steady.sweep`` selects how:

- ``zjacobi`` (default): **z-line Jacobi.** All columns are solved at once, with the horizontal neighbours taken
  from the previous sweep. ``steady.zline_sweeps`` (at least 2) sets the sweeps per iteration. More sweeps allow a
  larger stable CFL; the terrain cases use 16 sweeps with a CFL cap of 4.
- ``lusgs``: **line LU-SGS** (lower–upper symmetric Gauss–Seidel). The columns of each box are visited in order: a
  forward pass uses the already-updated lower neighbours, and a backward pass the updated upper neighbours. One sweep
  is a forward and a backward pass. When the wind blows along the sweep direction, a pass resolves the upwind
  coupling exactly, which allows a much larger CFL. The hill cases use two sweeps with a CFL cap of 16. When
  :math:`u` and :math:`v` have opposite signs, one sweep is stable only to a CFL of about 2.5 (two sweeps to about
  33), so use at least two sweeps. Box faces are coupled with lagged values, so the iteration path (not the answer)
  depends on ``mesh.max_grid_size``, and boxes shorter than the CFL number slow the iteration.

Four options refine the LU-SGS sweeps. They change the path to the solution, not the solution.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Input
     - Effect
   * - ``steady.line_ghosts``
     - ``lagged`` (default): neighbours across a box face take the previous pass's values, exchanged before every
       pass. ``zero``: no communication inside the sweeps; each box sweeps on its own, with :math:`\delta q = 0`
       outside it (a restricted additive Schwarz method).
   * - ``steady.line_overlap``
     - With ``zero``: each box sweeps a region grown by this many columns, then keeps only its own points.
   * - ``steady.line_order``
     - ``sequential`` or ``diagonal``. ``diagonal`` visits the diagonal planes :math:`i + j = \mathrm{const}`, one
       column per GPU thread. It gives the same result as ``sequential``, bit for bit.
   * - ``steady.line_alternate``
     - ``1``: every second sweep reverses the y direction. When :math:`u` and :math:`v` have opposite signs, one of
       the two sweep pairs then follows the wind.

On CPUs, ``lagged`` needs the fewest iterations. ``zero`` avoids communication inside the sweeps, which is expected
to matter on GPUs and on many ranks. Measurements are in :doc:`results`.

Pressure correction
^^^^^^^^^^^^^^^^^^^

With the SIMPLEC coefficient :math:`d = a_P - \sum a_{nb}`, the pressure correction solves

.. math::

   \nabla\cdot\Bigl(\frac{1}{d}\,\nabla\phi\Bigr) = \nabla\cdot\mathbf u^*, \qquad
   \mathbf u = \mathbf u^* - \frac{1}{d}\,\nabla\phi, \qquad p \leftarrow p + \alpha_p\,\phi .

SIMPLE neglects the neighbours' velocity corrections and needs a pressure under-relaxation. SIMPLEC assumes they
equal the cell's own, so :math:`\alpha_p = 1` works (``simplec.alpha_p``).

The Poisson solve stops at the relative tolerance ``simplec.pressure_reltol`` (1e-3). A looser tolerance leaves
divergence that drives :math:`w` through buoyancy, and the iteration diverges in stratified flow.

Point-implicit buoyancy and rotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two couplings oscillate in pseudo-time where :math:`\nu_t` is small: :math:`w` with :math:`\theta` (buoyancy,
frequency :math:`N`) and :math:`u` with :math:`v` (Coriolis, frequency :math:`f`). An explicit treatment amplifies
each oscillation by :math:`(1 + \omega^2\Delta\tau^2)^{1/2}` per step. The predictor therefore treats both with a
linearised backward Euler step:

.. math::

   w:\ \ a_P \mathrel{+}= \frac{g}{\theta_0}\max(\partial_z\theta, 0)\,\tau_\theta,\quad
   \mathrm{rhs} \mathrel{-}= \frac{g}{\theta_0}\,\tau_\theta R_\theta ; \qquad
   u, v:\ \ a_P \mathrel{+}= f^2\max(\tau_c, \tau_\mathrm{own}),\quad
   \mathrm{rhs} \mathrel{\mp}= f\bigl(\tau_c\langle R\rangle + \tau_\mathrm{own}(\overline R - \langle R\rangle)\bigr).

- :math:`\tau_\theta = 1/(1/\Delta\tau + a_{P,\theta})`: the step allowed by the :math:`\theta` rows.
- :math:`\langle R\rangle`: the horizontal mean, on each level, of the other component's residual;
  :math:`\overline R`: its local value.
- :math:`\tau_c = 1/\max(|f|/r, a_P - \sum a_{nb})`, with :math:`r` = ``simplec.rotation_relax`` (0.3).
- :math:`\tau_\mathrm{own}`: the row's own step.

The inertial oscillation is horizontally uniform, so only the level mean uses the long inertial step. Every added
term vanishes when :math:`R = 0`, so the steady solution does not change.

Vertical-velocity cap
^^^^^^^^^^^^^^^^^^^^^

On meshes with thin first cells, :math:`w` converges slowly. The vertical neighbour weights
:math:`a_{ZM}, a_{ZP}` of the :math:`w` rows cancel in :math:`d`, so the pressure correction moves :math:`w` too far
by the vertical von Neumann number

.. math::

   \mathrm{VNN} = \Delta\tau\,(a_{ZM} + a_{ZP}),

and the next predictor undoes most of the move. ``simplec.w_vnn_cap`` = :math:`V \ge 16` adds a term on the
:math:`w` rows only (not on the ground and top rows):

.. math::

   \frac{1}{\Delta\tau} \;\rightarrow\; \frac{1}{\Delta\tau} + \frac{a_{ZM} + a_{ZP}}{V}
   \quad\text{in } a_P \text{ and in } d,
   \qquad \mathrm{VNN}_\mathrm{eff} = \frac{\mathrm{VNN}}{1 + \mathrm{VNN}/V} < V .

The cap changes the path, not the steady answer: the predictor is in delta form, and :math:`d` enters only the
pressure correction (there is no Rhie–Chow interpolation on the staggered grid), so the fixed point is still
:math:`R = 0`. It is off by default. Values below 16 are refused: on flat ground the iteration diverges for
:math:`V \le 2` at a CFL of 16 and for :math:`V = 4` at a CFL of 64.

.. note::

   Over terrain, the cap can make the pressure correction unstable where thin cells meet steep slopes (roughly when
   :math:`1 + \mathrm{VNN}/V > 1/h^2`, with :math:`h` the local slope). Over terrain, the run prints a warning when
   the cap is on and the mesh aspect ratio is above the tested 50.78. A larger :math:`V` is safer. See :doc:`results` for the
   tested cases.

Pressure Poisson solver
-----------------------

The pressure correction, the mass-consistent initial field and the Tucker wall distance all solve

.. math::

   \frac{1}{J}\Bigl[\delta_\xi\bigl(a_x\beta_x f_x\bigr) + \delta_\eta\bigl(a_y\beta_y f_y\bigr)
   + \delta_\zeta\bigl(a_z\beta_z f_z\bigr)\Bigr] = \mathrm{rhs},

where :math:`f` are ERF's face gradients, terrain cross terms included, and :math:`\beta` are face coefficients
(:math:`1/d` for SIMPLEC). ``solver.poisson_backend`` selects the solver:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Back end
     - Method
   * - ``gmres_mlmg`` (default)
     - GMRES preconditioned by AMReX's ``MLTerrainPoisson`` geometric multigrid. The multigrid coarsens in x and y
       only and solves each column exactly, which suits thin cells.
   * - ``mlmg``
     - The same multigrid alone.
   * - ``gmres_fft``
     - ERF's ``TerrainPoisson`` operator with GMRES and ERF's FFT preconditioner (the flat stretched Laplacian: FFT
       in x and y, tridiagonal in z).

The default needs a nearly constant number of iterations per solve as the first layer becomes thinner, where the FFT
preconditioner needs more (:doc:`results`).

The solve stops at ``solver.poisson_reltol`` or ``solver.poisson_abstol``; the pressure correction uses
``simplec.pressure_reltol`` instead. ``gmres_mlmg`` measures the residual weighted by cell volume (:math:`J`),
``gmres_fft`` the plain residual and ``mlmg`` the maximum norm. The same tolerance therefore does not mean the same
divergence where :math:`J` varies. With no Dirichlet face the mean of the right-hand side is removed first.

Grid sequencing
---------------

Most of the iterations of a run remove large-scale errors, which the fine mesh removes slowly. Grid sequencing
removes them first on coarse meshes, where an iteration is cheap, and then starts each finer mesh from the coarser
solution. It uses the grid hierarchy of :doc:`grid_hierarchy`.

``steady.sequencing = 1`` with ``mesh.hierarchy_levels`` = :math:`L \ge 2` runs, on a fresh start:

#. **Coarsest level** :math:`L - 1`: the initial field (``init.*``) on that mesh, then SIMPLEC to a *loose stop*.
#. **Levels** :math:`L - 2, \dots, 1`: the solution of the level above is prolonged (interpolated, then made
   divergence-free), then SIMPLEC runs to a loose stop.
#. **Level 0** (the run's mesh): prolonged from level 1, then SIMPLEC to the normal convergence criteria.

Each level has its own boundary conditions, turbulence model, wall distance and Poisson solver. An iteration on
level :math:`l` costs about :math:`4^{-l}` of a level-0 iteration. ``steady.max_iterations`` applies to each level
separately. The log ends with a table of every level's iterations and wall time.

**Loose stop.** A coarse level stops when, for every residual :math:`q`,

.. math::

   R_q^n \le \max\!\left(R_q^{\mathrm{ref}}, \tau_R\right) 10^{-d_\ell} \quad\text{or}\quad R_q^n < \tau_R ,

with :math:`d_\ell` = ``steady.sequencing_drop_orders`` (2), or when the full convergence criteria hold.
``steady.sequencing_monitor_tol`` adds an optional monitor test; it is off by default, because right after a
prolongation the monitors settle before the residuals have dropped.

**Residual reference.** A prolonged start is already close to the solution, so a level's own first residuals are
small. Measuring the drop from them would hold level 0 to a stricter standard than a run without sequencing. Every
finer level therefore starts with the reference residuals :math:`S_q` of the coarsest level, which, like a run
without sequencing, starts from the initial field.

**Reference hold.** A driver normally runs its first ``steady.residual_ref_iterations`` iterations at
``steady.cfl0`` while it takes its reference. A seeded level already has its reference, so by default
(``steady.sequencing_hold = 0``) it skips this hold. From the first iteration it uses

.. math::

   r_n = \max_q \frac{R_q^n}{\max\bigl(R_q^{\mathrm{ref},n}, \tau_R\bigr)}, \qquad
   R_q^{\mathrm{ref},n} = \max\bigl(S_q, R_q^0, \dots, R_q^{\min(n, K - 1)}\bigr),

in the CFL ramp, which then climbs at once (still limited by ``steady.ser_max_growth``). ``steady.sequencing_hold =
1`` keeps the hold on every level. ``steady.sequencing_cfl_start = previous`` starts each seeded level at the CFL the
previous level reached, instead of ``steady.cfl0``.

**What it does not change.** The converged state is the same with and without sequencing, to the convergence
tolerance. A restart continues level 0 only; the coarse levels run on a fresh start. While the coarse levels run,
the level-0 solver stays allocated, so the peak memory is somewhat above that of level 0 alone.

On the 3-D hill, three levels reduce the level-0 iterations from 2091 to 706 and the wall time from 145 s to 54 s
(:doc:`results`).

Newton–Krylov stage (in development)
------------------------------------

A Newton–Krylov stage will finish the iteration once SIMPLEC is close to the solution. Each Newton step solves

.. math::

   (D + J)\,\delta q = -R(q)

by GMRES, with the Jacobian–vector products of :doc:`discretization`.

- :math:`J = \partial R/\partial q`, including the dependence of :math:`\nu_t` on :math:`q`.
- :math:`D`: :math:`1/\Delta\tau` on the momentum, :math:`\theta` and :math:`k` rows; 0 on continuity and boundary
  rows. :math:`D = 0` is plain Newton.
- The fields are scaled so GMRES weighs them alike.

The step is implemented and tested; the application does not run it yet. Its inputs (``steady.newton_*``) are
listed in :doc:`inputs`.

.. _convergence-criteria:

Convergence criteria
--------------------

A run is declared steady when the residuals and the monitored quantities have both converged.

**Residuals.** With :math:`R_q^n` the scaled norm of equation :math:`q` at iteration :math:`n`, every equation must
satisfy

.. math::

   R_q^n < \tau_R \quad\text{and}\quad R_q^n \le 10^{-N} R_q^\mathrm{ref}, \qquad
   R_q^\mathrm{ref} = \max\Bigl(\max_{m < K} R_q^m,\ \tau_R\Bigr) .

- :math:`\tau_R`: ``steady.residual_tol`` (5e-5).
- :math:`N`: ``steady.residual_drop_orders`` (4), the orders of magnitude the residual must drop.
- :math:`K`: ``steady.residual_ref_iterations`` (5). The reference is the peak over the first :math:`K` iterations,
  so the drop does not depend on how low the first iteration happened to start.

**Monitors.** Quantities of interest, such as the hilltop speed-up or :math:`u_*`, must also settle. Each monitor
:math:`M_j` must stay within :math:`T_j = \tau_M\max(|M_j|, s_j)` over the last :math:`W` iterations, and its
extrapolated limit (Aitken's method) must lie within the same distance:

.. math::

   \max_{n-W \le m \le n}|M_j^m - M_j^n| \le T_j, \qquad
   \left|\frac{\rho\,d_2}{1 - \rho}\right| \le T_j, \qquad \rho = \frac{d_2}{d_1},\ -1 < \rho < 1,

where :math:`d_1, d_2` are the changes over the two halves of the window, :math:`\tau_M` is
``steady.monitor_rel_tol`` (1e-3) and :math:`W` is ``steady.monitor_window`` (50). The extrapolation matters: for
a slowly converging monitor, a small change over the window alone can stop 10 % away from the limit.

The outcome is ``Steady``, ``NotSteady`` (listing what has not converged) or ``Diverged`` (a NaN or infinity).
Oscillations with a period up to :math:`4W` are never declared steady.

Stage switching, restart and reproducibility
--------------------------------------------

- **Stage switching.** The driver moves to the next stage after ``steady.switch_window`` iterations below
  ``steady.switch_residual``, and back when the residual grows by ``steady.switch_back_growth`` over a window.
- **Restart.** The driver records every iteration (stage, CFL, residuals, monitors, wall time, allocation units) in
  the checkpoint. A restarted run reproduces the uninterrupted run bit for bit on the same ranks and box layout.
- **Reproducibility.** Norms and monitors are reduced over all ranks, so every rank reaches the same verdict at the
  same iteration. A run on a different number of ranks agrees with the original to round-off, not bit for bit.
  ``steady.line_order = diagonal`` reproduces ``sequential`` bit for bit.
