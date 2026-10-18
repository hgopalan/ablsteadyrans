Steady driver
=============

The steady driver (sub-phase 1.13, ``Source/Driver.H``) runs the pseudo-time iteration to the steady state defined in
:doc:`steady`: it hands each iteration to a *stage*, sets the local pseudo-time step, records the residual and
monitor histories with the wall time and allocation units of every stage, switches between stages by rules with
hysteresis, and stops on the verdict of ``SteadyCriterion``. Phase 1 has one stage, SIMPLEC (sub-phase 1.14); the
Newton finish (phase 3) and the FAS multigrid (phase 4) plug into the same interface. The application runs the driver
from sub-phase 1.14, when it has a stage; this version builds and tests the framework.

Stages
------

A stage implements ``init(q)`` (called when the driver enters it), ``iterate(q, dtau)`` (one pseudo-time iteration
with the cell-centred local step), ``residual_norms()`` (the scaled residual norms after that iteration, one per
equation, bit-identical on every rank) and ``monitors()`` (the monitored quantities, :doc:`steady`). Every stage of a
run reports the same residual names and monitor specifications; the driver refuses a set that does not.

One iteration
-------------

1. The local pseudo-time step from the current CFL (below).
2. The current stage's ``iterate``.
3. Its residual norms and monitors go to the steady criterion and to the history (:ref:`driver-history`).
4. The relative residual

   .. math::

      r_n = \max_q \frac{R_q^n}{\max(R_q^{\mathrm{ref}}, \tau_R)},

   with :math:`R_q^{\mathrm{ref}}` the criterion's reference (the largest value over the first
   ``steady.residual_ref_iterations`` iterations) and :math:`\tau_R` = ``steady.residual_tol``, drives the SER ramp
   and the switch rules. While the reference is still being taken, :math:`r_n = 1`.
5. The verdict. While it is ``NotSteady``, the CFL of the next iteration follows the SER ramp and the switch rules
   run.

Local pseudo-time step
----------------------

In every cell (``Source/PseudoTime.H``),

.. math::

   \Delta\tau = \frac{\mathrm{CFL}}{\dfrac{|u_c|}{\Delta x} + \dfrac{|v_c|}{\Delta y} + \dfrac{|\Omega_c|}{J\,\Delta\zeta}
   + 2\nu\left(\dfrac{1}{\Delta x^2} + \dfrac{1}{\Delta y^2}\right)},

with :math:`u_c, v_c` the means of the cell's x- and y-face velocities, :math:`\Omega_c` the mean of its two z-face
contravariant fluxes (:doc:`transport`: the velocity through the terrain-following layers, so :math:`|\Omega_c| / (J
\Delta\zeta)` is the rate at which it crosses the cell's physical thickness) and :math:`\nu` the momentum diffusivity
(cell-centred; none when not given). **Why no vertical diffusive limit.** Every stage treats vertical diffusion
implicitly (the z-line solves of sub-phase 1.14), so it does not bound the step; the horizontal diffusion is explicit
and does. The rate is floored at :math:`10^{-12}` 1/s, so a cell at rest without diffusion gets a large finite step.
A terrain-following field (:math:`\Omega = 0`) has no vertical term however steep the terrain. With vertical
advection implicit as well (sub-phase 1.14) the vertical term is not a stability bound either; it is kept as a
conservative choice. On the initializer's fields it costs at most a factor 1.9 on the flanks of the bump case
(first layer 2.5 m, 25 m cells) and 1.6 on Askervein, and nothing aloft; the horizontal term dominates everywhere.
The face-area factors :math:`a_x/J, a_y/J` of the horizontal terms are within 3 % of 1 on those meshes and are
left out.

SER ramp
--------

The CFL starts at ``steady.cfl0`` (≈ 1: a robust start) and stays there while the criterion takes its reference.
Then, with :math:`p` = ``steady.ser_exponent``,

.. math::

   \mathrm{CFL}_{n+1} = \min\Bigl( \mathrm{CFL}_{\max},\; g\,\mathrm{CFL}_n,\; \mathrm{CFL}_0\, r_n^{-p} \Bigr),

:math:`\mathrm{CFL}_{\max}` = ``steady.cfl_max`` and :math:`g` = ``steady.ser_max_growth`` the growth allowed per
iteration (switched evolution relaxation: the step grows as the residual falls). Mulder and van Leer's form,
:math:`\mathrm{CFL}_{n+1} = \mathrm{CFL}_n (r_{n-1}/r_n)^p`, telescopes to the same law without the caps; with them
the reference form keeps growing towards :math:`\mathrm{CFL}_0 r^{-p}` after the growth cap has held it back, and
drops to :math:`\mathrm{CFL}_0 r^{-p}` at once on a residual spike. Below :math:`r_{\min} =
(\mathrm{CFL}_0/\mathrm{CFL}_{\max})^{1/p}` the law already gives the cap and the power is not taken, so a residual
that vanishes gives the cap, not an overflow. ``steady.ser_exponent = 0`` holds the CFL at ``cfl0``.

Switch rules
------------

With :math:`W` = ``steady.switch_window``:

- *forward* (to the next stage) after :math:`W` consecutive iterations with :math:`r_n \le` ``steady.switch_residual``;
- *back* (to the previous stage) when :math:`r_n \ge` ``steady.switch_back_growth`` :math:`\times\, r_{n-W}`: the
  residual has grown by that factor over the last window, which must lie within the current stage;
- *hysteresis*: no switch within :math:`W` iterations of the last one, so any window of :math:`W` iterations holds at
  most one switch.

Every switch is recorded (iteration, stages, reason, the scaled residuals) and printed with its reason and the
residuals, and the stage entered gets ``init``. With one stage nothing ever switches. A stage that always grows the residual the previous one
lowered gives a forward/back cycle of period at least :math:`2W`; the history shows it, and the stall detector of
sub-phase 4.2 is where such a cycle will be caught.

.. _driver-history:

History and accounting
----------------------

``DriverHistory`` holds one record per iteration (stage, CFL, relative residual, wall seconds, residuals, monitors),
the switches, one account per stage (iterations, wall seconds, allocation units from ``steady.machine`` and
``steady.nodes``: Kestrel CPU 10 AU per node-hour, Kestrel GPU 100, local 0, as the reference-run tooling of
sub-phase 1.3) and the steady criterion's history. An iteration's wall time is the slowest rank's, so every rank holds
the same history. It is written to a checkpoint as the text file ``<checkpoint>/Driver`` with hexadecimal floats and
read back bit for bit; ``Driver::restore`` continues from it (the criterion, the SER state, the switch state and the
stage, whose ``init`` runs again on the first step), so a restarted run reproduces the records and the verdict of an
uninterrupted one (the tests show this with a two-stage run interrupted mid-way). A stage keeps its own state in the
solution state ``q``, which the checkpoint holds; nothing else of a stage is checkpointed. If ``steady.cfl0`` or
``steady.cfl_max`` changed across the restart, the CFL continues within the new bounds. A missing, truncated or
inconsistent file stops the run, naming it (text after a complete history is ignored).

Inputs
------

Every key is in ``steady`` (:doc:`inputs` has the table with the ranges). The application does not read them yet: the
driver needs a stage (sub-phase 1.14). The reader ``read_steady_options()`` and the checks are tested.
