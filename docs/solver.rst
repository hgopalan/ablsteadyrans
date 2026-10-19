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
odd–even mode as the CFL grows. With opposite signs, each pass leaves one upwind neighbour lagged; one sweep is
stable only to a CFL of about 2.5 and two to about 33. Couplings across box faces use the previous pass, so the
iteration path depends on ``mesh.max_grid_size`` and boxes shorter than the CFL limit the solve. The hill results
are in :doc:`results`; ``tools/sweeps_cfl_scan.py --thresholds`` prints the scalar stability factors.

**Box-local sweeps, overlap and face ownership.** Sub-phase 2.2 adds two ways to treat the neighbours across a box
face (``Source/LineSweeps.H``):

- ``steady.line_ghosts = lagged`` (the default; sub-phase 2.1, bit for bit): ``delta.FillBoundary`` before every
  pass, so the neighbours across a box face or a periodic seam take the previous pass's values. Both copies of a
  shared face are solved; ``OverrideSync`` keeps the lowest box index's copy.
- ``steady.line_ghosts = zero``: no communication inside the sweeps. Each box sweeps the region it owns, grown by
  :math:`n_o` = ``steady.line_overlap`` columns in x and y, with :math:`\delta = 0` outside it, and keeps only the
  owned points. The system (:math:`A`, :math:`a_P`, the right-hand side) on the grown region comes from one
  ``FillBoundary`` of copies with :math:`n_o` ghost cells. Along a periodic axis one cell wide (y in a 2-D deck)
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
values as in the lexicographic order: the result is the sequential sweep's, bit for bit (UT-2.2.4). The GPU
timings and the per-machine defaults come in a later sub-phase.

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

Poisson solver
--------------

The initializer, the pressure correction and the Tucker wall distance solve

.. math::

   \frac{1}{J}\Bigl[\delta_\xi\bigl(a_x\beta_x f_x\bigr) + \delta_\eta\bigl(a_y\beta_y f_y\bigr)
   + \delta_\zeta\bigl(a_z\beta_z f_z\bigr)\Bigr] = \mathrm{rhs},

where :math:`f` are ERF's face gradients with the terrain cross terms and :math:`\beta` are face coefficients
(:math:`1/d` for SIMPLEC, :math:`\Sigma` for the initializer). ``solver.poisson_backend`` selects:

- ``gmres_fft`` (default): ERF's ``TerrainPoisson`` with ``amrex::GMRES`` and ERF's FFT preconditioner (the flat
  stretched Laplacian, FFT in x and y, tridiagonal in z);
- ``mlmg``, ``gmres_mlmg``: AMReX's ``MLTerrainPoisson`` with geometric multigrid, alone or under GMRES. It gives the
  same solution to 1e-9.

The solve stops at ``solver.poisson_reltol`` or ``solver.poisson_abstol``. With no Dirichlet face the mean of the
right-hand side is removed first. In single precision the rounding floor is about :math:`\varepsilon N_z^2/2`, so
use a tolerance near :math:`4\varepsilon N_z^2`.

Newton–Krylov stage
-------------------

Sub-phase 3.1 adds one Newton step solved by GMRES. The application does not run it yet (sub-phase 3.5).

.. math::

   (D + J)\,\delta q = -R(q), \qquad A = T^{-1}(D + J)\,S, \qquad b = -T^{-1}R(q), \qquad \delta q = S x .

- :math:`J = \partial R/\partial q` with :math:`\nu_t` frozen, applied by the Jacobian–vector product of
  :doc:`discretization`.
- :math:`D` is :math:`1/\Delta\tau` on momentum, :math:`\theta` and :math:`k`, and 0 on continuity and Dirichlet
  rows. :math:`D = 0` is plain Newton.
- :math:`S` and :math:`T` scale the fields so GMRES weighs them alike.

GMRES is flexible, so the preconditioner (sub-phase 3.3) may change between iterations; it is the identity for now.
FGMRES(:math:`m`) keeps :math:`2m + 3` vectors of state size. The run stops, naming
``steady.newton_mem_budget_gb``, if they do not fit.

Reproducibility
---------------

Residual norms and monitors use ``deterministic_sum``: each column is summed in k order, columns in storage order,
and each box total goes into its own slot of a vector that is reduced over the ranks. The result has the same bits
for any number of ranks on a fixed BoxArray, so every rank reaches the same verdict at the same iteration.

References
----------

- Van Doormaal, J. P., Raithby, G. D. (1984). *Numerical Heat Transfer* 7, 147–163.
- Mulder, W. A., van Leer, B. (1985). *J. Comput. Phys.* 59, 232–246.
