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

*Sloping ground.* On a :math:`z`-face the pressure correction scales the whole normal gradient, slope term included,
by :math:`\beta_z = 1/d_w`, while the :math:`x`- and :math:`y`-faces use :math:`\beta_h = 1/d_u, 1/d_v`. For one
Fourier mode (wavenumber :math:`k` along a slope :math:`h`, :math:`m` in :math:`\zeta`) the operator's quadratic form
is

.. math::

   \beta_h k^2 + \beta_z (1 + h^2)\, m^2 - h (\beta_h + \beta_z)\, k\, m ,

which is definite only while :math:`r = \beta_z/\beta_h` exceeds

.. math::

   r_c = \frac{h^2}{\left(\sqrt{1 + h^2} + 1\right)^2} .

Without the cap :math:`r \approx 1`. With it :math:`r = 1/(1 + \mathrm{VNN}/V)`, so thin cells on steep slopes cross
:math:`r_c`: on Askervein's steepest slope (23°, :math:`r_c = 0.041`) a cap of 16 fails once VNN exceeds 370, as with
a 0.4 m first cell. ``simplec.w_vnn_cap_slope_margin`` = :math:`s` (default 4) limits the term on each :math:`w` row
to

.. math::

   \min\left(\frac{a_{ZM} + a_{ZP}}{V},\; \frac{\max\bigl(0,\, 1/(s\,r_c) - 1\bigr)}{\Delta\tau}\right),

with :math:`h^2 = h_x^2 + h_y^2` from the face's node heights (the slope the operator's cross term uses). The limit
assumes that the horizontal coefficient and the vertical one without the cap are both :math:`1/\Delta\tau`. Then
:math:`r \ge s\,r_c` wherever :math:`s\,r_c < 1` (below 53° at :math:`s = 4`; steeper rows get no cap). On the
solver's own coefficients the margin is smaller where the vertical coefficient is stiffer (stable stratification,
the Rayleigh layer: at :math:`s = 4` and 32.5°, :math:`r \ge 3.0\,r_c` when it is twice :math:`1/\Delta\tau` and
:math:`1.7\,r_c` at five times), and larger in the first layer, where the surface drag stiffens the horizontal
coefficient. The definite band is also bounded above, :math:`r < 1/r_c`. On steep slopes the cap acts as a larger
local :math:`V`, and flat ground keeps :math:`(a_{ZM} + a_{ZP})/V` exactly. Giving the slope term :math:`\beta_h` instead (``FaceCoeffs::z_cross``) would
keep the operator definite for any :math:`r`, but the default Poisson back end has no split coefficient. With
:math:`s = 0` (the cap of #32) the run warns over terrain when the mesh aspect ratio is above the tested 50.78. See
:doc:`results` for the tested cases.

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

``steady.newton_gmres_max_iters`` bounds the GMRES iterations of one solve over all its restarts. A solve that
reaches it is unconverged; its update is used when the Krylov residual is finite and below :math:`|b|`, like any
other.

The Newton preconditioner
^^^^^^^^^^^^^^^^^^^^^^^^^

GMRES solves the Newton system with a right preconditioner :math:`M^{-1} \approx (D + J)^{-1}`, chosen by
``steady.newton_preconditioner`` (``SimplecPreconditioner``, ``Source/NewtonPreconditioner.H``). ``identity`` is
no preconditioner. ``simplec`` (the default) applies one SIMPLEC iteration to the linearized system
:math:`(D + J)\,\delta q = r`, starting from :math:`\delta q = 0`:

1. The momentum predictor of the SIMPLEC stage, :math:`(1/\Delta\tau + A_1)\,\delta u^* = r_u`, with the
   coefficients, rotation and buoyancy terms at :math:`q`.
2. The pressure correction with the continuity part :math:`r_p` as a source,

   .. math::

      \nabla\cdot\left(\frac{1}{a_P}\nabla p'\right) = \nabla\cdot\delta u^* - r_p, \qquad
      \delta u = \delta u^* - \frac{1}{a_P}\nabla p', \qquad \delta p = \alpha_p\,p',

   so that :math:`\nabla\cdot\delta u = r_p`. :math:`a_P` is the predictor's diagonal and :math:`\alpha_p` is
   ``simplec.alpha_p``.
3. The :math:`\theta` and :math:`k` line steps on their parts of :math:`r`. The change of the :math:`\theta`
   convection by :math:`\delta u` is a directional derivative,
   :math:`[C(F(u + \sigma\delta u)) - C(F(u))]/\sigma` with :math:`\sigma` as in the Jacobian-vector product, and
   the line operators use the transport of :math:`u`, so ``apply`` is linear in :math:`r` (to the finite-difference
   and Poisson tolerances). With the stage's denominator (below) the full :math:`\delta u` and the transport of
   :math:`u + \delta u` are used instead, as in the stage. The increment :math:`\delta k` is not clipped.
4. On the Dirichlet rows of :math:`u, v, w`, which are identity rows of the Newton system, :math:`\delta q = r`.

The sweeps are those of the SIMPLEC stage: ``steady.zline_sweeps`` of ``steady.sweep``. ``setup`` rebuilds the
coefficients once per Newton step, at the current state and :math:`\Delta\tau`. Inside the continuation loop
that is every attempt. At infinite :math:`\Delta\tau` the preconditioner uses the local time step at
``steady.cfl_max``.

**The denominator.** The SIMPLEC stage divides by :math:`d = a_P - \sum a_{nb}`. In the interior that is about
:math:`1/\Delta\tau`, so the pressure increment vanishes as :math:`\Delta\tau` grows. The preconditioner
therefore uses the SIMPLE denominator :math:`a_P`. On the 3-D hill of UT-3.3.3 the GMRES iterations per Newton step
then stay near 10 from CFL 10 to 1e4; with :math:`d` they grow from 18 to more than 350 (:doc:`results`). With :math:`d` the
preconditioner is exactly one iteration of the stage (UT-3.3.1).

**Shared faces.** A face shared by two boxes, or a periodic seam face and its image, is one unknown. :math:`r`
arrives with equal copies, and :math:`\delta q` is returned with every copy equal to the owner's.

**The pressure gauge.** Without a Dirichlet pressure face, a constant pressure is the null space of
:math:`D + J`. The preconditioner removes the net source of :math:`r_p` that the boundaries cannot carry, and
returns :math:`\delta p` with zero :math:`J`-weighted mean. The Newton update then keeps the mean of :math:`p`, and
a constant added to :math:`p` changes neither the preconditioner nor the solve (UT-3.3.2).

**Limits.**

- On the horizontally uniform Ekman column the pressure plays no part. The stiff part at large CFL is the balance
  of rotation and vertical diffusion, and ``simplec`` needs about as many GMRES iterations as ``identity``.
- The pressure solve's reductions depend on the rank order, so with ``simplec`` the Newton iterates are not bit
  for bit the same on different rank counts. With ``identity`` they are.
- Memory: about 25 state-size fields and the Poisson solver's hierarchy, on top of the Krylov vectors.

Pseudo-transient continuation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A continuation loop (``Continuation``, ``Source/Continuation.H``) wraps the Newton step. The application does not
run it yet. One continuation step from :math:`q_n`:

1. :math:`\Delta\tau` from the local time step at :math:`\mathrm{CFL}_n`, or no pseudo-time term (:math:`D = 0`)
   once the CFL has reached its cap and ``steady.newton_dtau_inf_at_cap = 1``.
2. The direction :math:`\delta q` from :math:`(D + J)\,\delta q = -R(q_n)` (one GMRES solve).
3. A backtracking line search on :math:`|R| = |T^{-1}R(\mathrm{fill}(q))|_2`, the norm of the Newton step:

   .. math::

      |R(q_n + \lambda\,\delta q)| \le (1 - c\,\lambda)\,|R(q_n)|, \qquad \lambda = 1, \beta, \beta^2, \ldots
      \ge \lambda_\mathrm{min},

   with :math:`c` = ``steady.newton_armijo_c``, :math:`\beta` = ``steady.newton_backtrack`` and
   :math:`\lambda_\mathrm{min}` = ``steady.newton_lambda_min``. Each trial point keeps :math:`k` at or above
   ``turb.tke_floor`` and goes through the ghost fill before :math:`R` is evaluated. At a finite
   :math:`\Delta\tau` the direction is a pseudo-time step and need not reduce :math:`|R|`. With :math:`\gamma` =
   ``steady.newton_ptc_growth`` > 1, when no :math:`\lambda` meets the condition there, the full step is accepted
   if :math:`|R(q_n + \delta q)| \le \gamma\,|R(q_n)|`; the CFL law then lowers the next CFL.
4. Accepted: :math:`q_{n+1} = q_n + \lambda\,\delta q`, and the CFL follows the SER law of the SIMPLEC stage,

   .. math::

      \mathrm{CFL}_{n+1} = \min\!\left(\mathrm{CFL}_\mathrm{max},\; g\,\mathrm{CFL}_n,\;
      \mathrm{CFL}_0 \left(\frac{|R_0|}{|R_{n+1}|}\right)^{p}\right),

   with ``steady.cfl0``, ``steady.cfl_max``, :math:`g` = ``steady.ser_max_growth`` and :math:`p` =
   ``steady.ser_exponent``. :math:`|R_0|` is the first finite :math:`|R|`, or one the caller supplies
   (``set_reference_residual``). While :math:`\Delta\tau` is infinite the law starts from
   :math:`\mathrm{CFL}_\mathrm{max}`.
5. A skipped update or a line search that accepts no :math:`\lambda` leaves :math:`q_n` unchanged, multiplies the
   CFL by ``steady.newton_cfl_cut`` (from infinite :math:`\Delta\tau`, :math:`\mathrm{CFL}_\mathrm{max}` times it)
   and retries, up to ``steady.newton_max_retries`` times. Then the step returns *failed*, and the Newton stage
   will fall back to SIMPLEC.

Every attempt prints one line: the step, :math:`|R|` before and after, the CFL (``inf`` for :math:`D = 0`),
:math:`\lambda`, the backtracks, the GMRES iterations and status, the Krylov relative residual and the outcome
(accepted, skipped or rejected). The continuation state (step count, CFL, the switch to infinite
:math:`\Delta\tau`, :math:`|R_0|` and the last :math:`|R|`) is written to and read from a checkpoint directory
(``Continuation``, hex floats); a restart gives the uninterrupted run's bits.

**Limits.**

- With ``steady.newton_ptc_growth = 1`` (the default) a pseudo-time direction that raises :math:`|R|` for every
  :math:`\lambda` is rejected, and a CFL cut need not help. On the Ekman column with ``steady.ser_exponent = 0.5``
  the retries run out at 3.5e-2 of :math:`|R_0|`; with 1.5 the column converges (:doc:`results`).
- Once :math:`\Delta\tau` is infinite it stays infinite through damped steps and unconverged solves; only a failed
  attempt returns it to a finite CFL. Returning to :math:`\mathrm{CFL}_\mathrm{max}` after such a step was tried
  and measured worse.
- With GMRES(30) and infinite :math:`\Delta\tau`, the 30 m hill stagnates near 6e-3 of :math:`|R_0|` for twelve
  steps at the 2000-iteration limit, then fails at 7e-8. With ``steady.newton_dtau_inf_at_cap = 0`` it reaches
  1e-8 in 23 steps, and with GMRES(200) in 8. The restart length and the cap will be tuned with the Newton stage.

Smoothness handling
^^^^^^^^^^^^^^^^^^^

The Jacobian-vector product is a forward difference of the residual. Across a kink of the residual it measures the
jump of the slope, not the derivative, and Newton loses its convergence rate. ``NewtonSmoothness``
(``Source/Smoothness.H``) handles the two kinks before the switch to Newton:

1. **The limiter.** The van Leer slope :math:`2ab/(a+b)` (0 when :math:`ab \le 0`) has a kink wherever :math:`a`
   or :math:`b` passes through 0: at every extremum, and in nearly uniform regions where the differences are
   round-off of either sign. The Newton stage's residual uses the smooth van Albada slope instead
   (``van_leer`` → ``van_albada``, with ``convection.smooth_eps_*``; ``upwind`` and ``van_albada`` are kept).
2. **The k floor.** With ``turb.frozen_switch = 1`` (the default) the floor decisions are frozen at the state where
   Newton starts. A cell at the floor there stays at :math:`k = k_\mathrm{floor}` for the whole stage, so
   :math:`k_f = \max(k, k_\mathrm{floor})` has no kink as its :math:`k` rises; its branch and :math:`R_t`
   smoothings are decided at that :math:`k`. A cell above the floor is floored again only if its :math:`k` falls to
   the floor or below, so no :math:`k \le 0` (from a line-search trial or a Jacobian perturbation) reaches the
   closure. The stable/unstable branch and the :math:`R_t` smoothing are not frozen: they are continuously
   differentiable (:doc:`turbulence_model`), and freezing them would replace the state's decisions by those of the
   switch point.

After the Newton stage the floor thaws and a final SIMPLEC pass with the run's own inputs (van Leer, live decisions)
converges the run's discrete equations from the Newton solution. ``leave`` reports how many floor decisions moved
during the stage. The module's state (active, the frozen decisions) goes to the checkpoint
(``<checkpoint>/Smoothness``, ``Level_0/closure_switches``) and is restored on any layout.

**Variation from the plan.** Plan rev. 6 also freezes the stability branch. In a neutral layer :math:`N^2` is
round-off and its sign flips: on the 2-D hill 286–507 of 512 cells changed branch during Newton. With the branch
frozen at the switch state Newton was 64 times further from convergence after 8 steps on the inversion and failed
at 5.7e-9 on the neutral van Albada run; without it both reach round-off (:doc:`results`).

**Cost.** The two limiters give different discrete solutions (0.4–0.6 % of the largest :math:`u` on the 2-D hill),
so the verification pass is not short: 3650–4316 SIMPLEC iterations, about as many as SIMPLEC needs from 1e-4 to
1e-9 on its own. On a van Leer run the swap therefore does not save work; with ``convection.limiter = van_albada``
for the whole run Newton solves the run's own equations, reaches round-off in about 8 steps, and needs no
verification pass (:doc:`results`). The application does not run the Newton stage yet (sub-phase 3.5).

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
