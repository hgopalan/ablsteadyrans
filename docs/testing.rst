Testing
========

Conventions
-----------

These follow the development plan.

- **Unit tests (UT-p.s.n)** use GoogleTest. They are deterministic, finish in seconds on 1–4 ranks, and run on
  every commit: ``ctest -L unit``.
- **Regression tests (RT-p.s.n)** run the application and compare with stored gold plotfiles using
  ``amrex_fcompare``. Their check scripts print *measured / target / tolerance / pass* rows.
  ``ctest -L smoke`` runs the structural checks CI does on every pull request; ``ctest -L physics`` runs the
  converged checks nightly on Kestrel.
- **Performance tests (PT-p.s.n)** run nightly or weekly on Kestrel. They track iterations, wall time and
  allocation-unit cost; a slowdown of more than 10 % fails.
- **Gold files** change only in a reviewed commit.

Every test is first shown to fail on the bug it guards against, or on a one-line mutant of the code it covers.
Tests the plan does not list are marked *beyond the plan*.

The tests run with the GPU CI flags (see :doc:`running`). Each test runs in its own directory under
``<build>/Tests/run/``.

Current tests
-------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - CTest
     - What it checks
   * - ``unit_np1``, ``unit_np2``, ``unit_np4``
     - The GoogleTest suite on 1, 2 and 4 ranks (below).
   * - ``unit_notiling``
     - The same suite with tiling off, as on a GPU.
   * - ``smoke_clean``, ``smoke_run``, ``smoke_restart``, ``smoke_restart_matches``
     - A 2-rank run; a 3-rank restart on another box layout, whose log must show the restart; the two plotfiles
       identical under ``amrex_fcompare`` at zero tolerance.
   * - ``abort_unknown_key``, ``abort_range``, ``abort_repeated_key``, ``abort_restart_mesh``,
       ``abort_restart_overwrite``
     - Decks that must stop, with the input named in the message.
   * - ``drift_check_selftest``
     - The ``extern_erf`` drift check on a throw-away ERF repository (:doc:`extraction`).

GoogleTest suite
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - ID
     - Check
   * - UT-1.2.1
     - Steady-criterion verdicts on converging, stalling (including a residual stuck exactly at the tolerance) and
       oscillating histories. The converging history becomes
       steady at iteration 136, pinned from an independent Python implementation of the definition.
   * - UT-1.2.2
     - *Beyond the plan.* A NaN or infinite residual or monitor gives ``Diverged``, permanently.
   * - UT-1.2.3
     - *Beyond the plan.* A monitor at round-off level around zero is steady with scale 1 and never steady with a
       scale below the noise.
   * - UT-1.2.4
     - *Beyond the plan.* The monitor window needs :math:`W+1` values; a change of exactly the tolerance passes and
       one just above it does not (values exact in binary).
   * - UT-1.2.5
     - *Beyond the plan.* The residual drop is measured from the early peak (iteration 38), not from the first
       value (iteration 49).
   * - UT-1.2.6
     - *Beyond the plan.* Parameter range checks name their keys; bad or missing monitor scales and size
       mismatches stop.
   * - UT-1.2.7
     - *Beyond the plan.* A monitor converging at rate 0.999 stops at iteration 4602, within the tolerance of its
       limit; the window test alone would stop at 1615 with a 2 % error.
   * - UT-1.2.8
     - *Beyond the plan.* Undamped oscillations with periods of 37, 100 and 200 iterations are never steady.
   * - UT-1.2.9
     - *Beyond the plan.* A residual that is zero over the reference iterations does not block the stop.
   * - UT-1.2.10
     - *Beyond the plan.* A criterion restored from a saved history gives the same verdict at every later
       iteration as one never interrupted; a mismatched history is refused.
   * - UT-1.2.11
     - *Beyond the plan.* A slowly growing monitor is not steady while its change over the window is still small.
   * - UT-1.3.1
     - Input parsing; unknown, repeated, multi-valued and out-of-range keys reported with their names.
   * - UT-1.3.2
     - Plotfile round trip, bit-identical, including the nodal heights and their ghost nodes.
   * - UT-1.3.3
     - Checkpoint and restart onto another box layout, bit-identical; ghost nodes restored without overwriting
       valid nodes.
   * - UT-1.3.4
     - The sum gives the same bits on 1, 2 and 4 ranks and within 16 ε Σ\|x\| of the exact sum; the maximum
       handles NaN and −0.
