Definition of steady
====================

A run is steady when two things hold at the same iteration: every residual is small and has dropped far enough,
and every monitored quantity has stopped moving and is close to its limit. ``Source/SteadyCriterion.H`` implements
this definition. The steady driver (:doc:`driver`) uses it to stop, and grid sequencing (sub-phase 2.5) will use it
for the strict stop on the finest level.

Residuals
---------

Let :math:`R_q^n \ge 0` be the scaled residual norm of equation :math:`q` (momentum, θ, k, ...) at iteration
:math:`n`. Its reference is its largest value over the first :math:`K` iterations, floored at the tolerance:

.. math::
   :label: residual-ref

   R_q^{\mathrm{ref}} = \max\Bigl( \max_{0 \le m < \min(K,\, n+1)} R_q^m ,\; \tau_R \Bigr) .

The residual criterion holds at iteration :math:`n` when, for every :math:`q`,

.. math::
   :label: residual-criterion

   R_q^n < \tau_R
   \quad\text{and}\quad
   R_q^n \le 10^{-N}\, R_q^{\mathrm{ref}} .

- :math:`\tau_R` is ``steady.residual_tol`` (default :math:`5\times10^{-5}`).
- :math:`N` is ``steady.residual_drop_orders`` (default 4).
- :math:`K` is ``steady.residual_ref_iterations`` (default 5).

**Why a reference over several iterations.** While the pseudo-time step ramps up from CFL ≈ 1, a residual can rise
for a few iterations. Measured from the early peak, the required drop depends on the problem rather than on how low
the first iteration happened to start.

**Why the floor.** A residual that is exactly zero at the start (θ with a uniform initial θ and no heat flux) would
otherwise have a zero reference, and the drop test would fail for any later round-off value, forever. With the
floor, such a residual must fall to :math:`10^{-N}\tau_R`, which round-off does.

**What this does not cover.** After a large initial spike (a poor initial field), the drop test is easy to meet and
the absolute tolerance :math:`\tau_R` does the work. A residual whose round-off floor lies above
:math:`10^{-N}\max(R^{\mathrm{ref}},\tau_R)` never meets the drop test; the verdict then says so.

Monitors
--------

Let :math:`M_j^n` be monitor :math:`j` at iteration :math:`n`: the hilltop speed-up ΔS, the friction velocity
u\*, the domain-integrated k, and so on. Each monitor has a scale :math:`s_j > 0`, set by the code that creates it.
Write :math:`T_j^n = \tau_M \max(|M_j^n|, s_j)`, where :math:`\tau_M` is ``steady.monitor_rel_tol`` (default
:math:`10^{-3}`, i.e. 0.1 %). With :math:`W` = ``steady.monitor_window`` (default 50), the monitor criterion holds
at iteration :math:`n \ge W` when, for every :math:`j`, both

.. math::
   :label: monitor-window

   \max_{n-W \le m \le n} \left| M_j^m - M_j^n \right| \;\le\; T_j^n

and the distance to the limit, extrapolated from the window, is within the same amount. With :math:`h = \lfloor
W/2 \rfloor` and the changes over the two halves of the window,

.. math::
   :label: monitor-aitken

   d_1 = M_j^{n-h} - M_j^{n-2h}, \qquad
   d_2 = M_j^{n} - M_j^{n-h}, \qquad
   \rho = d_2 / d_1 ,

.. math::
   :label: monitor-remaining

   -1 < \rho < 1
   \quad\text{and}\quad
   \left| \frac{\rho\, d_2}{1 - \rho} \right| \;\le\; T_j^n .

:eq:`monitor-remaining` is Aitken's estimate of :math:`|M_j^\infty - M_j^n|` for a monitor that converges
geometrically. A ratio :math:`\rho` outside :math:`(-1, 1)` means the monitor is not converging over the window.
When both :math:`|d_1|` and :math:`|d_2|` are below :math:`10^{-3}\,T_j^n` (round-off), the ratio carries no
information and :eq:`monitor-remaining` is skipped.

**Why the extrapolation.** A small change over the window does not bound the distance to the limit. For
:math:`M^n = M^\infty + A r^n`, the window test alone stops with an error of about :math:`\tau_M / (W(1-r))`: with
:math:`A = 0.1` and :math:`r = 0.9999`, it stops at iteration 50 with a 10 % error. With
:eq:`monitor-remaining`, the error at the stop is about :math:`\tau_M` for every :math:`r` tested
(0.95 to 0.9999). The cost is that a slowly converging monitor needs more iterations (4602 at :math:`r = 0.999`).

**The scale.** Below :math:`|M_j| = s_j` the test is absolute, with size :math:`\tau_M s_j`. Choose :math:`s_j`
so that this size is the smallest change that matters. For the speed-up ratio ΔS, :math:`s = 1` means 0.1 % of the
reference wind speed; over flat terrain, where ΔS is zero, it keeps round-off of order :math:`10^{-17}` from
blocking the stop. For a quantity that never approaches zero, a scale far below its value makes the test purely
relative. There is no default scale.

**What this does not cover.** An oscillation whose period is much longer than the window looks like convergence
inside one window. Undamped oscillations with a period up to :math:`4W` (200 iterations at the default) are never
declared steady (tested); longer periods can be. Choose :math:`W` above the longest expected period. Detecting a
limit cycle from the whole history is sub-phase 4.2.

Verdict
-------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Verdict
     - Meaning
   * - ``Steady``
     - :eq:`residual-criterion`, :eq:`monitor-window` and :eq:`monitor-remaining` all hold.
   * - ``NotSteady``
     - Some criterion does not hold yet. The verdict lists each one, with its value and the input that sets it.
   * - ``Diverged``
     - A residual or monitor has been NaN or infinite at some iteration. This verdict is final.

Restart and MPI
---------------

The recorded history (iteration count, residual references, the monitor window, the ``Diverged`` state) is
available as ``SteadyHistory``, and ``restore()`` continues from it. The steady driver writes it to its checkpoint
(:doc:`driver`), so a restarted run reaches the same verdict at the same iteration as a run without the restart.

Every rank must pass bit-identical residuals and monitors (all-reduced, for example with ``deterministic_sum``).
Then every rank reaches the same verdict, and a steady loop exits on all ranks at the same iteration.

Cases that may have no steady state
-----------------------------------

Some flows have no steady RANS solution, or the iteration does not settle into one. These cases are labelled
*possibly non-steady* in :doc:`cases`:

- stable stratification;
- Coriolis forcing (inertial oscillations);
- lee-side separation.
