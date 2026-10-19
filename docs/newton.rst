Newton–Krylov core
==================

Sub-phase 3.1 (``Source/Newton.H``, ``Source/Newton.cpp``) adds one Newton step solved by a Jacobian-free Krylov
method. The application does not run it yet: the switch from SIMPLEC to Newton is sub-phase 3.5, pseudo-transient
continuation and the line search are sub-phase 3.2, and the preconditioner is sub-phase 3.3. The unit tests
UT-3.1.1 to UT-3.1.3 (:doc:`testing`) exercise it.

The Newton system
-----------------

One step solves

.. math::

   \left(D + J\right)\delta q = -R(q), \qquad q \leftarrow q + \delta q,

where :math:`q` is the state (the C-grid face velocities :math:`u, v, w` and the cell values :math:`p, \theta, k`),
:math:`R(q)` the full second-order residual of :doc:`residual` (per unit volume), :math:`J = \partial R/\partial q`
its Jacobian with :math:`\nu_t` frozen, and :math:`D` the pseudo-time term: :math:`1/\Delta\tau` on the momentum
faces (the face step of the predictor, :doc:`predictor`), on θ and on k, and 0 on the continuity rows and on the
Dirichlet rows of the normal velocities (identity rows of :math:`J`). :math:`D` is the plan's :math:`V/\Delta\tau`
divided by the cell volume :math:`V`, the convention of the residual. Without a :math:`\Delta\tau` field,
:math:`D = 0` (plain Newton).

GMRES works in scaled variables so that its inner products weigh the fields alike:

.. math::

   A x = b, \qquad A = T^{-1}(D + J)\,S, \qquad b = -T^{-1} R(q), \qquad \delta q = S x,

with :math:`S = \mathrm{diag}(s_u, s_u, s_u, s_p, s_\theta, s_k)` the scales of ``JvpScales`` (10 m/s,
100 m²/s², 10 K, 0.1 m²/s²) and :math:`T` the same with 1 on the continuity rows (already in 1/s). The solver is
``amrex::GMRES<NewtonVector, NewtonOperator>`` with ``setFlexible(true)``, so the preconditioner may change from one
iteration to the next.

The Jacobian-vector product
---------------------------

:math:`J v` is never formed. ``NewtonOperator::apply`` builds :math:`\delta q = S x` with ghost cells, fills them with
the homogeneous boundary fill, and calls ``jvp`` (:doc:`residual`):

.. math::

   J v \approx \frac{R(q + \epsilon v) - R(q)}{\epsilon}, \qquad
   \epsilon = h\;\frac{1 + \mathrm{mean}\,|\tilde q|}{\mathrm{rms}\,\tilde v},

where :math:`\tilde q = q/s` and :math:`\tilde v = v/s` are the scaled state and direction and
:math:`h` = ``steady.newton_fd_eps_scale``. For a residual with bounded second derivative :math:`B`, the error of
the forward difference is

.. math::

   \frac{\|J_\mathrm{FD} v - J v\|}{\|J v\|} \approx
   \epsilon\,\frac{\|B(v, v)\|}{\|J v\|} + \frac{\epsilon_\mathrm{mach}\,\|R\|}{\epsilon\,\|J v\|},

a truncation term growing with :math:`\epsilon` and a round-off term falling as :math:`1/\epsilon`. Their sum is
smallest near :math:`\epsilon \propto \sqrt{\epsilon_\mathrm{mach}}`, so the default
:math:`h = \sqrt{\epsilon_\mathrm{mach}}` of ``amrex::Real`` (1.49e-8 in double, 3.45e-4 in single) adapts to the
precision. UT-3.1.1 measures the curve; ``tools/newton_numerics.py`` re-derives it.

**Why.** The forward difference costs one residual per Krylov iteration and needs no Jacobian storage. **What.** The
step :math:`\epsilon` is the only free parameter; it is set relative to the scaled state so that every field moves by
about :math:`h` times its own scale, on any mesh and any box layout.

Krylov vectors, reductions and the preconditioner hook
------------------------------------------------------

The unknown spans fields of four index types (x-, y-, z-faces and cells), which one ``MultiFab`` cannot hold, so the
vector type is ``NewtonVector``: six ``MultiFab`` without ghost cells on the state's layout. A face shared by two
boxes, or by a periodic image, is one unknown; the inner products count it once (the owner mask). ``dotProduct`` and
``norm2`` sum each box in double, all-reduce a vector of box sums with one contributor per slot, and add the slots
in box order, so the result has the same bits on every rank and for any rank count with the same ``BoxArray``.

``NewtonPreconditioner`` is the hook for sub-phase 3.3: ``setup(q, dtau)`` once per step and ``apply(dq, r)`` once
per GMRES iteration, in physical units. Without a hook the preconditioner is the identity (a copy).

Memory
------

FGMRES(:math:`m`) keeps :math:`m + 1` Arnoldi vectors, :math:`m` preconditioned vectors and two work vectors:
:math:`2m + 3` vectors of state size. ``NewtonKrylov`` prints the estimate at construction (the bytes of one vector
times the count, for the rank with the most points) and aborts, naming ``steady.newton_mem_budget_gb`` and
``steady.newton_gmres_restart``, when it exceeds the budget (1 GB = 10\ :sup:`9` bytes). On one 80 GB H100 a
512 × 512 × 128 box needs 1.61 GB per vector in double: restart 30 (63 vectors, 102 GB) is refused under the 70 GB
default and restart 20 (43 vectors, 69.4 GB) fits (UT-3.1.3, beyond the plan). The residual
evaluations allocate their own temporaries outside the estimate.

Inputs
------

``read_newton_options()`` reads the ``steady.newton_*`` keys; each range check aborts with the key and the value.
The defaults are provisional: sub-phase 3.6 (PR-28, the Phase 3 checkpoint) tunes them. Until sub-phase 3.5 the
application does not run the Newton step, so it refuses these keys as unknown inputs (with a note saying so).

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Input
     - Default
     - Range and meaning
   * - ``steady.newton_gmres_restart``
     - ``30``
     - Integer, ≥ 5 and ≤ 200. The GMRES restart length :math:`m`.
   * - ``steady.newton_fd_eps_scale``
     - :math:`\sqrt{\epsilon_\mathrm{mach}}`
     - > 0 and ≤ 1e-2. The relative step :math:`h` of the Jacobian-vector product.
   * - ``steady.newton_krylov_rtol``
     - ``1e-2``
     - > 0 and < 1. GMRES stops when the scaled linear residual falls below this fraction of :math:`\|b\|`. Keep it
       well above :math:`\sqrt{\epsilon_\mathrm{mach}}`: below the accuracy of the forward difference (about 1e-8
       relative in double) GMRES stagnates.
   * - ``steady.newton_mem_budget_gb``
     - ``70``
     - > 0 and ≤ 80. Krylov memory per rank, in GB (10\ :sup:`9` bytes).

Diagnostics
-----------

``NewtonKrylov::step`` returns, and ``last()`` keeps, a ``NewtonStepInfo``: the GMRES iterations and status (0
converged, 1 iteration cap), :math:`\|b\|` before the step, the final Krylov residual (absolute and relative), the
smallest and largest :math:`\epsilon` used and the number of ``jvp`` calls. With ``verbose > 0`` the step prints one
line with these values.

Measured behaviour
------------------

On the Ekman column of UT-3.1.2 (2 × 2 × 32 cells, wall law, Coriolis), plain Newton from the geostrophic wind with
GMRES(200) to 1e-10 reduces the scaled residual by 0.403, 0.171, 0.0131, 2.09e-4 and 9.41e-8 in its first five
steps: the last two gain 3.7 and 7.0 orders, and :math:`r_{n+1}/r_n^2` stays between 0.25 and 0.65 (quadratic
convergence). GMRES iterations per step against the restart length on the same column, three steps:

.. list-table::
   :header-rows: 1
   :widths: 16 14 14 14 14 14 14

   * - rtol \\ restart
     - 5
     - 10
     - 20
     - **30**
     - 60
     - 100, 200
   * - 1e-2
     - 3, 5, 2000 (cap)
     - 3, 5, 8
     - 3, 5, 8
     - **3, 5, 8**
     - 3, 5, 8
     - 3, 5, 8
   * - 1e-6
     - cap
     - cap
     - 637, 743, 884
     - 614, 729, 556
     - 330, 355, 310
     - 94, 96, 96

At the default tolerance 1e-2 every restart from 10 up gives the same iterations; the default 30 leaves a margin.
The tight tolerance needs a restart close to the Krylov dimension of the problem.
