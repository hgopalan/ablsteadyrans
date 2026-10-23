Testing
========

Conventions
-----------

- **Unit tests** use GoogleTest; each is named ``Suite.Test``. They are deterministic, finish in seconds on
  1–4 ranks, and run on every commit: ``ctest -L unit``.
- **Regression tests** run the application and compare with stored gold plotfiles using
  ``amrex_fcompare``. Their check scripts print *measured / target / tolerance / pass* rows.
  ``ctest -L smoke`` runs the structural checks CI does on every pull request; ``ctest -L physics`` runs the
  converged checks nightly on Kestrel.
- **Performance tests** run nightly or weekly on Kestrel. They track iterations, wall time and
  allocation-unit cost; a slowdown of more than 10 % fails.
- **Gold files** change only in a reviewed commit.

Every test is first shown to fail on the bug it guards against, or on a one-line mutant of the code it covers.

The one-line mutants of the ``gmres_mlmg`` default (``Source/Poisson.H``, ``Source/Initializer.cpp``,
``Source/App.cpp``), each run against ``unit_np2``, ``rt_steady_askervein`` and (M1) ``rt_init_askervein``; every one
is killed:

.. list-table::
   :header-rows: 1
   :widths: 6 50 44

   * - Mutant
     - Change
     - Failing tests
   * - P1
     - No fallback: the projection with ``alpha_v != alpha_h`` keeps ``gmres_mlmg``
     - ``Initializer.ProjectionBackEndFollowsTheAlphaRatio`` and others in ``unit_np2``, ``rt_steady_askervein``,
       ``rt_init_askervein``
   * - P2
     - The old default ``gmres_fft``
     - ``unit_np2``, ``rt_steady_askervein`` (its log pattern ``pressure gmres_mlmg``)
   * - P3
     - The split coefficient set for every :math:`\alpha \ne 1`, also :math:`\alpha_v = \alpha_h`
     - ``Initializer.EqualAlphasAreThePlainProjection`` (``unit_np2``)
   * - P4
     - The fallback not reported
     - ``unit_np2``, ``rt_steady_askervein`` (its log line)
   * - P5
     - The stage's log line without the back end
     - ``rt_steady_askervein``
   * - P6
     - The fallback to ``mlmg`` instead of ``gmres_fft``
     - ``unit_np2``, ``rt_steady_askervein``
   * - P7
     - The ``solver.poisson_precond`` refusal without its hint about the default
     - ``PoissonB.UnsupportedCombinationsAreRefused`` (``unit_np2``)

The one-line mutants of the line-sweep variants (``Source/LineSweeps.cpp`` and ``Source/Driver.cpp``), each run against
``LineSweeps.*``, ``Driver.*`` and ``Predictor.Lusgs*`` on two ranks by ``tools/line_sweep_mutants.py``; every one
is killed:

.. list-table::
   :header-rows: 1
   :widths: 6 50 44

   * - Mutant
     - Change
     - Failing tests
   * - M1
     - Overlap ignored (swept region not grown)
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``, ``LineSweeps.OverlapsReachTheSameState``
   * - M2
     - Face updated twice (box keeps its high face)
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``, ``LineSweeps.OverlapsReachTheSameState``, ``LineSweeps.ConvergedStateDoesNotDependOnTheLayout``
   * - M3
     - Periodic seam updated twice (first box keeps the domain's low face)
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``, ``LineSweeps.OverlapsReachTheSameState``, ``LineSweeps.ConvergedStateDoesNotDependOnTheLayout``
   * - M4
     - Forward wavefront skips its last plane
     - ``LineSweeps.DiagonalOrderingEqualsTheSequentialSweep``
   * - M5
     - Backward wavefront skips plane 0
     - ``LineSweeps.DiagonalOrderingEqualsTheSequentialSweep``
   * - M6
     - Ghost δ not zeroed (neighbours' values exchanged between passes)
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``
   * - M7
     - Non-owned face copies keep their own value (owner mask all 1)
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``, ``LineSweeps.OverlapsReachTheSameState``, ``LineSweeps.ConvergedStateDoesNotDependOnTheLayout``
   * - M8
     - One-cell periodic axis not folded into the diagonal
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``
   * - M9
     - ``steady.line_overlap`` not read
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M10
     - ``steady.line_ghosts`` not read
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M11
     - ``steady.line_order`` not read
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M12
     - ``steady.line_overlap`` above 8 accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M13
     - A negative ``steady.line_overlap`` accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``, ``LineSweeps.RefusalsAndNames``
   * - M14
     - Overlap with lagged ghosts accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``, ``LineSweeps.RefusalsAndNames``
   * - M15
     - Overlap wider than the narrowest box accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``, ``LineSweeps.RefusalsAndNames``
   * - M16
     - Line key given with ``zjacobi`` accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M17
     - Odd sweeps not reversed (``steady.line_alternate`` ignored)
     - ``LineSweeps.AlternatingOrderingSweepsBothDiagonals``
   * - M18
     - Reversed forward wavefront skips its last plane
     - ``LineSweeps.AlternatingOrderingSweepsBothDiagonals``
   * - M19
     - Reversed sequential forward pass walks j ascending
     - ``LineSweeps.AlternatingOrderingSweepsBothDiagonals``
   * - M20
     - ``steady.line_alternate`` not read
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M21
     - ``steady.line_alternate = 2`` accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``, ``LineSweeps.AlternatingOrderingSweepsBothDiagonals``
   * - M22
     - ``steady.line_alternate`` given with ``zjacobi`` accepted
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
   * - M23
     - ``steady.line_alternate = 1`` with ``steady.zline_sweeps = 1`` accepted (the option read and ignored)
     - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``

The one-line mutants of the transfers between levels (``Source/Transfer.cpp``), each run against ``Transfer.*`` on two ranks with the
unit-test flags (SNaN fill, FPE trap) by ``tools/transfer_mutants.py``. 23 are killed. Three survive and the script
lists them as equivalent: T6 and T9 change no value, and T24 changes only the ``gmres_fft`` solver's path, so its result
agrees with the original's to the solver tolerance. A first run, without the SNaN and FPE flags and before ``TargetNodeHeightsComeFromTheHierarchy``,
``OneCellDirectionIsConstant`` and ``WrongLevelPairsAbort`` were added, left T10, T14, T18 and T19 alive as well. "Crash" marks a mutant that stops the run with a
signal in that test (SIGILL on macOS arm64), so the tests after it did not run.

.. list-table::
   :header-rows: 1
   :widths: 6 50 44

   * - Mutant
     - Change
     - Failing tests
   * - T1
     - Cell means not ``detJ``-weighted
     - ``RestrictionConservesFaceFlux``
   * - T2
     - x/y-face flux without the fine area factor
     - ``RestrictionConservesFaceFlux``, ``RoundTripIsTheClosedFormSmoother``
   * - T3
     - x/y-face flux not divided by r
     - ``RestrictionConservesFaceFlux``, ``RoundTripIsTheClosedFormSmoother``
   * - T4
     - z-face flux divided by r instead of :math:`r^2`
     - ``RestrictionConservesFaceFlux``, ``RoundTripIsTheClosedFormSmoother``
   * - T5
     - Coarse u not divided by the coarse area factor
     - ``RestrictionConservesFaceFlux``, ``RoundTripIsTheClosedFormSmoother``
   * - T6
     - Coarse Ω not divided by the coarse area factor
     - Equivalent: :math:`a_z = 1` on every mesh
   * - T7
     - w restricted in place of Ω
     - ``RestrictionConservesFaceFlux``
   * - T8
     - No zero-gradient rows for the slope term's u
     - ``RestrictionConservesFaceFlux`` (crash)
   * - T9
     - No periodic fill for the slope term's u
     - Equivalent: the slope term reads u on its own cell's x-faces, no horizontal ghost
   * - T10
     - Coarse ``z_phys_nd`` not set by ``restrict_state``
     - ``TargetNodeHeightsComeFromTheHierarchy``
   * - T11
     - Cell-centred position shifted by half a fine cell
     - ``ProlongationOfBilinearFieldsIsExact``, ``RoundTripIsTheClosedFormSmoother``, ``OneCellDirectionIsConstant``
   * - T12
     - No linear extrapolation at a non-periodic face (reads past the last centre)
     - ``ProlongationOfBilinearFieldsIsExact`` (crash)
   * - T13
     - Periodic direction extrapolates instead of reading the periodic image
     - ``RoundTripIsTheClosedFormSmoother``
   * - T14
     - One-cell direction reads a second coarse cell (zero weight, unfilled ghost)
     - ``OneCellDirectionIsConstant`` (crash)
   * - T15
     - x weights swapped on the upper y row
     - ``ProlongationOfBilinearFieldsIsExact``, ``ProlongedTkeIsAboveTheFloor``, ``RoundTripIsTheClosedFormSmoother``, ``OneCellDirectionIsConstant``
   * - T16
     - w interpolated with the nodal x axis
     - ``ProlongationOfBilinearFieldsIsExact`` (crash)
   * - T17
     - Coarse image filled without the periodic images
     - ``ProlongationOfBilinearFieldsIsExact`` (crash)
   * - T18
     - Fine ``z_phys_nd`` not set by ``prolong_fields``
     - ``TargetNodeHeightsComeFromTheHierarchy``
   * - T19
     - Fine-below-coarse level check dropped
     - ``WrongLevelPairsAbort``
   * - T20
     - k not clipped
     - ``ProlongedTkeIsAboveTheFloor``
   * - T21
     - Clipped-cell count not summed over ranks
     - ``ProlongedTkeIsAboveTheFloor``
   * - T22
     - Infinite ``turb.tke_floor`` accepted
     - ``ProlongedTkeIsAboveTheFloor``
   * - T23
     - No boundary fill before the projection
     - ``ProlongationAndProjectionIsDivergenceFree``
   * - T24
     - Projection corrects the Dirichlet rows too (d = 1 + mask)
     - Equivalent with ``gmres_mlmg`` and ``mlmg``: every Dirichlet row is a domain face with a Neumann condition on p′, so the flux there is 0. With ``gmres_fft`` only to the Poisson tolerance: its preconditioner's cell coefficients next to the face change
   * - T25
     - No projection off the Dirichlet rows (d = −mask)
     - ``ProlongationAndProjectionIsDivergenceFree``, ``ProlongedTkeIsAboveTheFloor``, ``GatheredLevelEqualsDistributed``
   * - T26
     - p′ not restored after the projection
     - ``ProlongationAndProjectionIsDivergenceFree``

The one-line mutants of the Newton stage (``Source/NewtonStage.cpp`` and the switch rules in ``Source/Driver.cpp``),
each run against ``Driver.*`` and ``NewtonStage.*`` on two ranks with the unit-test flags by
``tools/newton_stage_mutants.py`` (macOS arm64, AppleClang). 13 of 15 are killed. N14 and N15 survive: they matter only
where the SIMPLEC preconditioner loses GMRES at a large CFL (the full-size 2-D hill), which no unit test reaches in
seconds; a Newton study on Askervein (not yet run) would measure them.

.. list-table::
   :header-rows: 1
   :widths: 6 44 50

   * - Mutant
     - Change
     - Failing tests
   * - N1
     - The fallback waits for the hysteresis window
     - ``Driver.NewtonSwitchOnceAndFallback``
   * - N2
     - The forward switch ignores the monitors
     - ``Driver.NewtonSwitchOnceAndFallback``
   * - N3
     - The forward switch ignores ``accepts_switch``
     - ``Driver.NewtonSwitchOnceAndFallback``
   * - N4
     - ``leave`` not called on a switch
     - ``Driver.NewtonSwitchOnceAndFallback``, ``HandsBackForTheFinalPassAfterTheSwap``
   * - N5
     - The monitor test against the window's first record, off by one
     - ``Driver.NewtonSwitchOnceAndFallback``
   * - N6
     - The fallback after one failure whatever the key
     - ``FallsBackAfterFailedSteps``, ``RestartInsideTheStage``
   * - N7
     - A failed step does not count
     - ``FallsBackAfterFailedSteps``
   * - N8
     - Converged before any accepted step
     - ``FallsBackAfterFailedSteps``
   * - N9
     - No hand-back for the final pass
     - ``HandsBackForTheFinalPassAfterTheSwap``
   * - N10
     - The continuation's start CFL not reset on entry
     - every ``NewtonStage`` test
   * - N11
     - The continuation state not restored on restart
     - ``RestartInsideTheStage``
   * - N12
     - No limit on failed visits
     - ``FallsBackAfterFailedSteps``
   * - N13
     - No idle hand-back
     - ``HandsBackWhenConvergedWithoutAVerdict``
   * - N14
     - The Newton cap not used (the SIMPLEC cap, plain Newton after one step)
     - survives (above)
   * - N15
     - The preconditioner of an infinite :math:`\Delta\tau` at the SIMPLEC cap
     - survives (above)

The one-line mutants of the smoothness options (``Source/Smoothness.cpp`` and the floor freeze and stratification gate in
``Source/Closure.cpp``), each run against ``Smoothness.*`` and ``Closure.*`` on two ranks with the unit-test flags by
``tools/smoothness_mutants.py`` (macOS arm64, AppleClang). 12 of 13 are killed. S5 survives: it removes the gate
from the decisions that ``freeze_switches`` and ``changed_switches`` take, and the floor scope's held bit does not
depend on :math:`N^2` (``evaluate`` gates :math:`N^2` itself; S4 and S13 cover it there).

.. list-table::
   :header-rows: 1
   :widths: 6 44 50

   * - Mutant
     - Change
     - Failing tests
   * - S1
     - The floor scope holds every decision
     - ``FloorScopeHoldsOnlyTheFloor``, ``NewtonCycleReturnsTheSimplecSolution``
   * - S2
     - ``changed_switches`` compares every bit
     - ``FloorScopeHoldsOnlyTheFloor``, ``NewtonCycleReturnsTheSimplecSolution``
   * - S3
     - The floor scope decides from the live state (frozen bits ignored)
     - ``FloorScopeHoldsOnlyTheFloor``, ``ModuleInputsAbortsAndCheckpoint``
   * - S4
     - The gate with :math:`S^2` not squared
     - ``RiGateClosedFormAndLimits``
   * - S5
     - The gate not applied in ``Closure::stratification``
     - survives (above)
   * - S6
     - ``turb.ri_gate`` not read
     - ``Closure.InputsAreReadAndChecked``
   * - S7
     - No limiter swap
     - ``ModuleInputsAbortsAndCheckpoint``, ``NewtonCycleReturnsTheSimplecSolution``
   * - S8
     - ``leave`` does not thaw
     - ``ModuleInputsAbortsAndCheckpoint``, ``NewtonCycleReturnsTheSimplecSolution``
   * - S9
     - ``enter`` freezes every decision
     - ``ModuleInputsAbortsAndCheckpoint``, ``NewtonCycleReturnsTheSimplecSolution``
   * - S10
     - The checkpoint's scope ignored on read
     - ``ModuleInputsAbortsAndCheckpoint``
   * - S11
     - A held cell decided at its own :math:`k`, not the floor (negative :math:`\nu_t`)
     - ``FloorScopeHoldsOnlyTheFloor``, ``ModuleInputsAbortsAndCheckpoint``
   * - S12
     - ``evaluate`` holds every decision under the floor scope
     - ``FloorScopeHoldsOnlyTheFloor`` (crash: NaN at :math:`k \le 0`)
   * - S13
     - :math:`P_b` from the gated :math:`N^2`
     - ``RiGateClosedFormAndLimits``

The one-line mutants of the continuation (``Source/Continuation.cpp``, the SER law in ``Source/PseudoTime.cpp`` and the
GMRES limit in ``Source/Newton.cpp``), each run against ``Continuation.*``, ``Newton.*`` and ``Driver.*`` on two ranks
with the unit-test flags by ``tools/continuation_mutants.py`` (run on macOS arm64 with AppleClang). #36 reported 40
mutants killed but kept no list; these 40 are a new set. The first run left C20, C25 and C36 alive: C25 because
``LineSearchRejectsOvershoot`` no longer pins the trajectory, and C20 and C36 because no test covered them.
``LineSearchRejectsOvershoot`` now checks the
:math:`k` floor, and ``RestoreCapsTheCflAndTheFloorIsChecked`` covers the other two. All 40 are killed with GCC; with
Clang and AppleClang 39 are, and C2 cannot be told apart from the original.

.. list-table::
   :header-rows: 1
   :widths: 6 44 50

   * - Mutant
     - Change
     - Failing tests
   * - C1
     - Armijo without the sufficient decrease (any decrease accepted)
     - ``ArmijoOnKnownQuadratic``
   * - C2
     - Armijo compares a NaN trial (no isfinite guard)
     - ``ArmijoOnKnownQuadratic`` with GCC (verified with GCC 16 on arm64, one process), whose compare with NaN raises
       the invalid flag the test checks (``fcmpe`` on aarch64, ``comisd`` on x86-64). Survives with Clang and
       AppleClang, which emit a quiet compare; the script reports it as expected there
   * - C3
     - Backtracking factor fixed at 1/2
     - ``ArmijoOnKnownQuadratic``
   * - C4
     - No lower bound on lambda
     - ``ArmijoOnKnownQuadratic``, ``GrowthBoundAtFiniteDtau``, ``LineSearchRejectsOvershoot``
   * - C5
     - Backtracks not counted
     - ``ArmijoOnKnownQuadratic``, ``LineSearchRejectsOvershoot``
   * - C6
     - Full-step residual not kept (the growth bound never passes)
     - ``ArmijoOnKnownQuadratic``
   * - C7
     - Growth bound strict at the bound
     - ``ArmijoOnKnownQuadratic``
   * - C8
     - Growth-bound step not re-evaluated (last trial kept)
     - ``ArmijoOnKnownQuadratic``, ``GrowthBoundAtFiniteDtau``
   * - C9
     - Growth-bound flag not recorded
     - ``ArmijoOnKnownQuadratic``, ``GrowthBoundAtFiniteDtau``
   * - C10
     - SER restarts from the stored CFL after infinite dtau
     - ``DtauRampFollowsSer``
   * - C11
     - Never switches to infinite dtau
     - ``DtauRampFollowsSer``, ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``LineSearchRejectsOvershoot``, ``RestartIsBitIdentical``
   * - C12
     - Last residual not stored
     - ``DtauRampFollowsSer``, ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``ReferenceResidualFromCaller``, ``RestartIsBitIdentical``
   * - C13
     - Accepted steps not counted
     - ``DtauRampFollowsSer``, ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``ReferenceResidualFromCaller``, ``RestartIsBitIdentical``
   * - C14
     - No CFL cut after a failed attempt
     - ``DtauRampFollowsSer``, ``LineSearchRejectsOvershoot``, ``SkippedStepCutsCflAndKeepsState``
   * - C15
     - No floor on the cut CFL
     - ``DtauRampFollowsSer``
   * - C16
     - A failed attempt keeps infinite dtau
     - ``DtauRampFollowsSer``, ``LineSearchRejectsOvershoot``
   * - C17
     - Growth bound also at infinite dtau
     - ``DtauRampFollowsSer``
   * - C18
     - Checkpoint drops the infinite flag
     - ``RestartIsBitIdentical``
   * - C19
     - Reader ignores r_last
     - ``RestartIsBitIdentical``
   * - C20
     - Restore keeps a CFL above a lowered ``steady.cfl_max``
     - ``RestoreCapsTheCflAndTheFloorIsChecked``
   * - C21
     - ``set_reference_residual`` ignores its argument
     - ``ReferenceResidualFromCaller``
   * - C22
     - ``set_reference_residual`` accepts 0
     - ``ReferenceResidualFromCaller``
   * - C23
     - No layout check in step
     - ``RestartIsBitIdentical`` (crash)
   * - C24
     - :math:`|R_0|` never taken from the first step
     - ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``ReferenceResidualFromCaller``, ``SkippedStepCutsCflAndKeepsState``
   * - C25
     - Trial k not clipped at the floor
     - ``LineSearchRejectsOvershoot``
   * - C26
     - Accepted step not copied back to q
     - ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``LineSearchRejectsOvershoot``
   * - C27
     - Finite dtau even when infinite
     - ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``
   * - C28
     - One retry fewer
     - ``LineSearchRejectsOvershoot``, ``SkippedStepCutsCflAndKeepsState``
   * - C29
     - ``run()`` never stops at the tolerance inside the loop
     - ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``ReferenceResidualFromCaller``
   * - C30
     - ``steady.newton_armijo_c`` accepted up to 1
     - ``InputReadersAndRangeChecks``
   * - C31
     - ``steady.newton_max_retries`` = 0 accepted
     - ``InputReadersAndRangeChecks``
   * - C32
     - ``steady.newton_ptc_growth`` below 1 accepted
     - ``InputReadersAndRangeChecks``
   * - C33
     - ``steady.newton_ptc_growth`` not read
     - ``InputReadersAndRangeChecks``
   * - C34
     - ``steady.newton_dtau_inf_at_cap`` not read
     - ``InputReadersAndRangeChecks``
   * - C35
     - A repeated key accepted
     - ``InputReadersAndRangeChecks``
   * - C36
     - ``turb.tke_floor`` not checked
     - ``RestoreCapsTheCflAndTheFloorIsChecked``
   * - C37
     - SER without the growth limit
     - ``DtauRampFollowsSer``, ``HillSameOnAnyRankCount``, ``Driver.ReleasedHoldRampsFromTheFirstIteration``, ``Driver.SerRampFollowsTheLawAndItsCaps``
   * - C38
     - SER exponent sign flipped
     - ``DtauRampFollowsSer``, ``GrowthBoundAtFiniteDtau``, ``HillSameOnAnyRankCount``, ``ReferenceResidualFromCaller``, ``RestartIsBitIdentical``, ``Driver.ReleasedHoldRampsFromTheFirstIteration``, ``Driver.SerRampFollowsTheLawAndItsCaps``
   * - C39
     - GMRES iteration limit not applied
     - ``GmresIterationLimit``
   * - C40
     - ``steady.newton_gmres_max_iters`` = 0 accepted
     - ``InputReadersAndRangeChecks``

The tests run with the GPU CI flags (see :doc:`running`). Each test runs in its own directory under
``<build>/Tests/run/``.

Each case has its own folder ``Tests/<group>/<case>/`` with a ``README.md``, its input decks, its check script and
its gold files (``Gold/``); ``Tests/README.md`` lists the groups. Check scripts used by more than one case and the
test-only drivers are in ``Tests/Supporting_Files/``.

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
   * - ``abort_repeated_turb_key``, ``abort_repeated_surface_key``
     - The column deck with ``turb.max_geom_lscale`` and ``bc.surface_z0`` repeated on the command line: each stops,
       naming the key. Before, these two readers checked only the value count and the last value won.
   * - ``warn_w_vnn_cap_terrain``, ``warn_w_vnn_cap_terrain_margin``, ``warn_w_vnn_cap_terrain_below``,
       ``abort_w_vnn_cap_floor``, ``abort_w_vnn_cap_slope_margin``, ``abort_w_vnn_cap_slope_margin_alone``
     - The thin-cell strip on a 50 m ridge (``Tests/Infrastructure/Input_Validation/inputs_w_vnn_cap_terrain``, max
       aspect ratio about 134), one iteration. With ``simplec.w_vnn_cap = 16`` and the default slope margin the run
       prints no warning and its banner names the margin (4); a margin of 0.5, or a margin given without the cap, is
       refused, naming the key. With ``simplec.w_vnn_cap_slope_margin = 0`` and a 0.7 m first cell the run exits 0 and
       prints the warning that names the key, the tested aspect ratio 50.78 and the slow convergence without the cap;
       with a 4 m first cell (aspect ratio below 50.78) it runs without the warning; with 4 it stops on the floor of
       16. The flat strip of ``rt_thin_cell_w_vnn_cap`` (aspect ratio above 50.78, no slope) prints no warning.
   * - ``rt_askervein_mesh``
     - ERF's Askervein mesh coarsened to 60 × 60 × 18 (20 m × 1.08 layers, BTF, ERF's terrain
       file): statistics of ``z`` and ``detJ`` and the quality report against ``Tests/Discretization/Terrain_Mesh/Gold/askervein_mesh.json``.
       Double-precision builds only: the statistics are compared at double-precision tolerances.
   * - ``rt_steep_ridge_quality``
     - A 600 m ridge under STF: the report against its gold, including both warnings (STF capped, top
       not flat). Double-precision builds only.
   * - ``rt_hill2d_run``, ``rt_hill2d_matches``
     - A run without a solve writes the mapped-grid plotfile of the 2-D hill; ``amrex_fcompare`` against
       ``Tests/Discretization/Terrain_Mesh/Gold/hill2d_plt00000`` at relative tolerance 1e-12 (1e-5 in single precision). The plotfile holds
       cell fields only, so the nodal heights are not compared here; ``Mesh.MetricsEqualErf`` compares them with ERF.

       The golds of the three mesh tests above were written by GUST itself (``Tests/Supporting_Files/check_mesh.py --write``):
       they detect changes, they do not show correctness. Correctness of the same code paths comes from ``Mesh.MetricsEqualErf``
       (ERF's own output) and the closed-form unit tests.
   * - ``rt_askervein_hierarchy``, ``rt_askervein_hierarchy_level1_matches``, ``rt_askervein_hierarchy_level2_matches``
     - The Askervein mesh of ``rt_askervein_mesh`` with STF and a three-level hierarchy
       (``inputs_mesh_askervein_hierarchy``: 60, 30 and 15 cells a side, the same 18 layers, 1 and 3 ground passes).
       In double precision, level 0's ``z`` and ``detJ`` statistics and every level's quality report against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_mesh.json`` (the per-level headers and the summary table included); in both
       precisions, ``z`` and ``detJ`` of levels 1 and 2 against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_level{1,2}_plt00000`` (relative 1e-12; 1e-5 in single precision). The golds
       were written by gust: they detect changes. The coarse grounds were checked against an independent
       rebuild (``MeshHierarchy.GroundSmoothingPasses``).
   * - ``rt_steep_ridge_hierarchy_quality``
     - The steep ridge of ``rt_steep_ridge_quality`` on 32 × 4 × 16 cells with three levels: the report
       of every level against ``Tests/Discretization/Grid_Hierarchy/Gold/steep_ridge_hierarchy_mesh.json``, with the STF and top warnings (no
       coarse-level hint: levels 1 and 2 cap no more than level 0), and the per-level headers and summary table. Double precision only. The gold was written
       by gust: it detects changes.
   * - ``rt_askervein_hierarchy_hint``
     - The Askervein hierarchy with unsmoothed coarse grounds (``terrain.hierarchy_passes = 0``): levels 1
       and 2 cap STF at 9 and 10 z-levels against level 0's 8, so the report gives the coarse-level hint; against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_p0_mesh.json`` (written by GUST). Double precision only.
   * - ``rt_hierarchy_restart_first``, ``rt_hierarchy_restart``, ``rt_hierarchy_restart_level{1,2}_matches``
     - A checkpoint of the Askervein mesh written without a hierarchy (2 ranks), restarted with three levels on 3
       ranks: the rebuilt levels 1 and 2 against the golds of ``rt_askervein_hierarchy``.
   * - ``abort_hierarchy``, ``abort_hierarchy_passes_unused``, ``abort_hierarchy_tool``
     - ``MeshHierarchy.InvalidLayoutsAbort`` at the application: n_cell x not divisible by 4 for three levels, a box size not a multiple of 4
       and passes out of range; passes given without a hierarchy; a test-only tool (``wall_distance``) given a
       hierarchy. Each message names its input.
   * - ``abort_sequencing_range``, ``abort_sequencing_coarsest``, ``abort_sequencing_unused``, ``abort_sequencing_value``, ``abort_sequencing_no_stage``
     - Grid sequencing's inputs at the application (the truly neutral hill deck with overrides): a loose drop of 9
       orders against the strict 8, a monitor tolerance of 1 and sequencing without a hierarchy, all three named
       in one abort; a coarsest level of 4 × 1 cells (64 × 16 over five levels); ``steady.sequencing_drop_orders``
       without ``steady.sequencing = 1``; ``steady.sequencing = 2``; ``steady.sequencing`` without the stage.
       The range test also refuses ``steady.sequencing_cfl_start = previous`` with ``steady.sequencing_hold = 1``,
       and the unused test the two keys without sequencing.
   * - ``abort_sequencing_hold_value``, ``abort_sequencing_cfl_start_value``
     - ``steady.sequencing_hold = 2`` and ``steady.sequencing_cfl_start = last``: each refused, naming its input.
   * - ``rt_sequencing_hold_kept``, ``rt_sequencing_hold_released``, ``rt_sequencing_cfl_previous``
     - Old against new ``steady.sequencing_hold``, on the two-level restart deck of the truly
       neutral hill (only level 0 is seeded): with the hold (``1``) level 0's first three iterations run at
       ``steady.cfl0`` with r = 1; without it (``0``, the default) r < 1 from the first iteration and the CFL leaves
       ``steady.cfl0`` at the second; with ``steady.sequencing_cfl_start = previous`` level 0 starts at level 1's
       last CFL (4, after 150 iterations). Seconds. ``rt_sequencing_restart`` also checks that the restored level-0
       history keeps the released hold.
   * - ``abort_poisson_options``
     - ``poisson_projection`` with ``solver.poisson_reltol`` repeated and out of range: both named.
       (``abort_unknown_key`` checks that the application refuses ``solver.poisson_reltol``.)
   * - ``rt_poisson_askervein``, ``rt_poisson_askervein_matches``
     - ``poisson_projection`` projects a fixed wind on the Askervein mesh (30 × 30 × 18). Iterations at
       most 10 % above ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein.json``, divergence after below 10 × the printed tolerance of
       before, and
       :math:`\phi` against ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_plt00000`` at relative 1e-8. Double precision only (the
       golds are double-precision runs); the gold was written by GUST, so it detects changes.
   * - ``rt_poisson_askervein_mlmg``, ``rt_poisson_askervein_mlmg_matches``
     - The same projection with back end B (MLMG, 6 iterations): iterations at most 10 % above
       ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_mlmg.json`` (written by GUST's first run) and the divergence check, and
       :math:`\phi` against back end A's gold
       ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_plt00000`` at relative 1e-8, since both solve the same operator.
   * - ``rt_walldist_hill2d``, ``rt_walldist_hill2d_matches``
     - ``wall_distance`` on ERF's Neutral_Hill_2D mesh: ERF's wall-distance checks (``Tests/Discretization/Wall_Distance/check_walldist.py``:
       max relative ≤ 0.15, mean ≤ 0.03, ≤ 0.2 cells below 100 m, all finite), and the field against
       ``Tests/Discretization/Wall_Distance/Gold/walldist_hill2d_plt00000`` at relative 1e-8 (a gold written by GUST). Double precision only.
   * - ``rt_init_askervein``, ``rt_init_askervein_matches``
     - The mass-consistent initial field on the Askervein mesh (30 × 30 × 18; ``Tests/Initialization/Mass_Consistent_Askervein/inputs_init_askervein``:
       log law u* = 0.6 m/s, z₀ = 0.03 m along x, inflow at xlo, outflow at xhi, :math:`\alpha_h = 1`,
       :math:`\alpha_v = 0.5`; divergence 4.05e-1 → 6.1e-12 1/s). ``Tests/Supporting_Files/check_init.py askervein``: iterations at most 10 % above the gold's 13,
       divergence after at most 10 × ``solver.poisson_reltol`` times before, divergence before (the printed 10 digits),
       max :math:`|w|` and the largest horizontal speed at the gold's values (1e-6), :math:`\theta` = 300 K everywhere, :math:`p' = 0`,
       :math:`k > 0` (``Tests/Initialization/Mass_Consistent_Askervein/Gold/init_askervein.json``); the state against ``Tests/Initialization/Mass_Consistent_Askervein/Gold/init_askervein_plt00000``
       at relative 1e-8. Two ranks; double precision only. :math:`\lambda` is not in the plotfile: the state
       plotfile is the regression artefact, and the unit tests pin :math:`\lambda`.
   * - ``rt_init_bump_trend``
     - The hilltop speed-up over FastWindTerrain's bump (``Tests/Initialization/Bump_Speedup_Trend/inputs_init_bump``,
       ``Tests/Initialization/Bump_Speedup_Trend/terrain_bump_fwt.txt``: its ``regtests/phase6_solve_correction`` hill and power-law inflow) for
       :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (five runs on two ranks), against FastWindTerrain's own at commit
       648f8f1 (``Tests/Initialization/Bump_Speedup_Trend/Gold/init_bump_trend.json``, written from its plotfiles by ``Tests/Initialization/Bump_Speedup_Trend/fwt_bump_speedup.py``): the
       same ranking of the five values at 40, 80 and 150 m above the ground (at 20 m two values tie within 0.002 in
       both codes), every speed-up above 1 at 20, 40, 80 and 150 m, the two codes within 0.02, 0.01, 0.005 and 0.003
       of each other at those heights (measured 0.011, 0.0057, 0.0020, 0.0008), and GUST's own values and
       iterations pinned (1e-6, 10 %). Double precision only (:doc:`initial_conditions`).
   * - ``rt_bc_neutral_flat``, ``rt_bc_neutral_flat_matches``
     - ``boundary_fluxes`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, 1000 m) with a fixed log-law state
       (u* = 0.35 m/s, z₀ = 0.1 m, 15°), the surface layer on the ground and a slip top: u* recovered at every ground
       cell to 1e-10, the ground stress :math:`u_*^2(\cos, \sin)` to 1e-10, :math:`k_w = u_*^2/C_{\mu 0}^2`; the
       residual of every cell equals its closed form (``Tests/Discretization/Boundary_Fluxes/check_bc.py``: arithmetic face ν times the log law's
       differences, the wall flux below, no stress through the slip top) to 1e-6 of :math:`u_*^2/\Delta z` for u and
       v (1e-13 measured) and 1e-9 K/s for θ (the sounding's diffusion with the top gradient), with the first cell's
       truncation −0.0971 and the top cell's lost flux +1.000 of :math:`u_*^2\cos/\Delta z` pinned; R_k and R_w zero;
       every field finite. The plotfile against ``Tests/Discretization/Boundary_Fluxes/Gold/bc_neutral_flat_plt00000`` at relative 1e-10 (a gold
       written by GUST). Double precision only.
   * - ``rt_closure_column``, ``rt_closure_column_matches``
     - ``closure_column`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, neutral) with the wall law's
       equilibrium state (u* = 0.35 m/s, z₀ = 0.1 m, k = u*²/C_μ0²), the uncapped geometric length and the
       logarithmic face mean (``Tests/Discretization/Closure_Column/check_closure.py``): ν_t = κ u* (z + z₀) and L = κ (z + z₀) at every cell
       (1e-9), u* recovered (1e-10), R_u and R_v zero below the top cell (1e-10 of u*²/Δz; the top cell loses
       u*²/Δz through the slip top, 1e-8), P_s, ε and R_k = −(P_s − ε) equal to their closed forms from the discrete
       strain (1e-9 of ε; the second cell's P_s/ε = 1.45 is the centred difference of ln over 2Δz), R_θ and R_w zero.
       The plotfile against ``Tests/Discretization/Closure_Column/Gold/closure_column_plt00000`` at relative 1e-10 (a gold written by GUST).
       Double precision only.
   * - ``rt_steady_column``
     - The application on the 1-D neutral column (``Tests/Canonical_RANS/Neutral_Column/inputs_steady_column``: 4 × 4 × 64,
       1000 m, surface layer z₀ = 0.1 m, slip top, uniform θ, the force F = u*²/H with u* = 0.35 m/s, the log law with
       u* = 0.25 m/s as the first guess, CFL cap 50) to the steady verdict (6 orders, 1e-8).
       ``Tests/Canonical_RANS/Neutral_Column/check_steady_column.py``: the verdict; the residual of u down 6 orders from its reference (the
       ``Driver`` history of the checkpoint); mean u* = √(FH) to 1e-4 (the force balance); the first cell's speed the
       wall law's (1e-6); the second cell within 1 % of the log law and the lowest eight within 2 % (the stress of the
       force-driven column is u*²(1 − z/H) and the discrete strain overshoots, ``rt_closure_column``); the wall k = u*²/C_μ0²
       (1e-6); ν_t of the second cell within 6 % of κ u* (z + z₀); every column identical (1e-10); every field finite.
       Measured: 2231 iterations, u* = 0.3499999, second cell 0.9 % low, cells 2–7 1.3–1.9 % low. Double precision
       only (six orders and 1e-8 are below single's floor); the restart chain below runs in both.
   * - ``rt_steady_column_constants``
     - The column of ``rt_steady_column`` with ``turb.kappa = 0.40`` and ``turb.Cmu0 = 0.55``: the same checks, with
       the wall law, the wall k and ν_t taken with those constants (``--kappa``, ``--cmu0``). The wall model now
       takes the closure's constants; before, it kept 0.41 and 0.5562, and the first cell was 2.4 % off the wall law (0.40/0.41 − 1).
       Double precision only.
   * - ``rt_steady_column_rough``, ``rt_steady_column_rough_floor_one``
     - The column on rough ground (``inputs_steady_column_rough``: 4 × 4 × 64, 1 m layers, 64 m,
       z₀ = 0.5 m, F = u*²/H with u* = √0.125 = 0.354 m/s). The first centre, 0.5 m, is below (e − 1) z₀, where
       ln((d + z₀)/z₀) = ln 2 < 1. With the default ``bc.surface_denominator_floor = log_law``: the checks of
       ``rt_steady_column``, with the cells above the first measured against the profile of the shallow column's
       stress u*²(1 − z/H) (``--stress-profile``, a closed form). With ``one``: the first cell's speed is
       (u*/κ) max(ln((z + z₀)/z₀), 1) to 1e-6, 1/ln 2 = 1.44 times the log law's. Measured: 347 and 411 iterations;
       log_law second cell 0.70 % low, lowest eight within 0.92 %. The earlier wall law (the ``one`` floor) fails the first
       (first cell 0.862 m/s against 0.598). Double precision only.
   * - ``rt_steady_column_restart_first``, ``…_continued``, ``…_straight``, ``…_matches``
     - 40 iterations of the column (``inputs_steady_column_restart``, the first guess on the
       command line), continued to 80 from the checkpoint (the driver and closure histories restored, the inflow keys
       of the deck accepted on the restart) against 80 straight on the same two ranks: the plotfiles identical under
       ``amrex_fcompare`` at relative 1e-14.
   * - ``rt_steady_hill2d_restart_first``, ``…_continued``, ``…_straight``, ``…_matches``
     - The same chain over terrain: 5 iterations of the 2-D hill (``inputs_steady_hill2d_restart``), continued to 10
       from the checkpoint, against 10 straight on two ranks, at relative 1e-10. The restart fills the ghost cells
       before the first step, whose local Δτ reads u and v below the ground and above the top (Ω on those faces).
       Without that fill the runs differ by 8.5e-4 (v) and the test fails; with it, by 1.5e-12. The remaining
       round-off comes from the ground fill, which is not a fixed point over terrain. Labels ``regression;smoke``, about
       1.5 s on two ranks. Double precision only.
   * - ``rt_steady_neutral_flat``
     - The application on ERF's Neutral_ABL_Flat (``Tests/Canonical_RANS/Neutral_ABL_Flat/inputs_steady_neutral_flat``: 8 × 8 × 64,
       geostrophic wind 10 m/s at latitude 90°, f = 1e-4 1/s, the sounding's inversion from 700 m, slip top with
       0.03 K/m, z₀ = 0.1 m, ERF's closure defaults; four sweeps, CFL cap 2, the log law with u* = 0.35 m/s as the
       first guess) to the steady verdict. ``Tests/Canonical_RANS/Neutral_ABL_Flat/check_steady_neutral.py``: the verdict, then ERF's own
       ``check_neutral.py --physics`` (``extern_erf/``) on the turbulence plotfile: every field finite, k ≥ 0, the
       capped κ(z + z₀) length (1e-3), the AL01 stability functions of R_t (1e-6), K_mv = c_μ √k L (5 %), ε =
       C_μ0³ k^1.5/L (5 %), the wall cell's k balance (1 %), u* within 0.35 ± 0.15 m/s, the wall k within 5 % of
       u*²/C_μ0², the log law in the lowest four cells within 10 %, K_mv of the second cell within 30 % of κ u* z.
       Measured: the verdict after about 2850 iterations (13 s on two ranks), u* = 0.377 m/s, every ERF check
       passing, on one to three ranks with both compilers (the bursts of an earlier version of the predictor's
       rotation step are explained and removed in :doc:`solution_algorithm`). Label ``regression`` (the run time); double only.
       The 1 % comparison with an ERF reference run is not in this version.
   * - ``rt_thin_cell_w_vnn_cap``
     - The w rows' VNN cap on a flat inflow/outflow strip with a 0.7 m first cell
       (``Tests/Canonical_RANS/Thin_Cell_Flat/inputs_thin_cell_flat``: 16 × 2 × 40, dx/dz₀ = 134, line LU-SGS with
       two sweeps at CFL 16, the verdict at 6 orders and 5e-7). ``check_w_vnn_cap.py`` runs it with
       ``simplec.w_vnn_cap`` 0 and 16 to the verdict: both reach it, the banner names the cap, the capped run needs
       at most a fifth of the iterations, and u, v, w of the two steady states agree to 2e-4 m/s; then a capped run of
       40 iterations continued to 80 from its checkpoint equals 80 straight in every field, bit for bit. Measured:
       7038 and 322 iterations (22 times fewer; 7038 and 299 before the residual scales used the depth of the layers),
       u and w within 8.5e-9 and 5.8e-9 m/s (since the Rayleigh layer ends at the mesh's top level; before, 6487 and
       687 iterations and 5.6e-5 m/s). About 50 s on two ranks; label
       ``regression``; double only.
   * - ``rt_steep_ridge_w_vnn_cap``, ``rt_steep_ridge_w_vnn_cap_margin_off``
     - The slope margin of the cap on a steep ridge with a 0.4 m first cell
       (``Tests/Canonical_RANS/Steep_Ridge_Thin_Cell/inputs_steep_ridge_thin_cell``: 48 × 2 × 40, an Agnesi ridge
       150 m high with a 150 m half-width, max slope 32.5°, max aspect ratio 372, the verdict at 3 orders), 3500
       iterations at most. With ``simplec.w_vnn_cap = 16`` and the default margin the run is steady (1765
       iterations; 1652 before the residual scales used the depth of the layers, 1895 before #44); with the margin off (the cap of #32) the pressure correction stalls and it is not. The run
       without the cap takes 2881 iterations (3983 before the Rayleigh layer ended at the mesh's top level, when it
       was not steady at 3500 either). Measured to 6 orders before that change: cap off 14154, margin 2 4142, and the
       steady states agree to 3e-8 m/s.
   * - ``rt_steady_hill2d``
     - ERF's Neutral_Hill_2D (``Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d``: 128 × 1 × 64, the Witch of Agnesi of
       100 m and 500 m, the flat deck's sounding and forcing, Rayleigh damping of w and θ from 700 m, 16 sweeps, CFL
       cap 4) to the steady verdict with a two-order drop of the u residual; ``Tests/Supporting_Files/check_steady_hill.py`` then runs
       ERF's ``check_hill.py --physics`` on the turbulence plotfile: the wall distance against the exact distance to
       the ridge, finite fields, k ≥ 0, the length-scale bound, the wall cell's k balance, the crest speed-up in the
       three lowest cells within 0.2–0.8 (Jackson and Hunt's 0.4), the speed-up positive in the lowest ten, the
       upstream u* within 0.25–0.55 m/s, the upstream log law in the lowest three cells within 15 %. Measured: the
       verdict after 4576 iterations (49 s on two ranks), speed-up 0.59/0.50/0.43, u* = 0.354 m/s, every check
       passing (:doc:`solution_algorithm`). Label ``regression``; double only.
   * - ``rt_steady_hill2d_lusgs``
     - The 2-D hill with the line LU-SGS sweeps (``inputs_steady_hill2d_lusgs``: two symmetric sweeps,
       CFL cap 16) to the verdict with a two-order drop and ERF's ``check_hill.py --physics``. Measured: 3048
       iterations (23 s on two ranks), every check passing, the steady state within the convergence tolerance of
       the Jacobi deck's (:doc:`results`). Label ``regression``; double only; it runs in the column
       test's directory (``run/rt_steady_hill2d_column``) so the spin-up check can read both histories.
   * - ``rt_steady_hill2d_column``, ``rt_column_guess_spinup``
     - The 2-D hill LU-SGS deck started from the converged column (``inputs_steady_hill2d_column``) to
       the verdict with a two-order drop and ERF's ``check_hill.py --physics``; then ``Tests/Canonical_RANS/Neutral_Hill_2D/check_column_guess.py``
       on the Driver histories of the profile run (``rt_steady_hill2d_lusgs``, now in the same run directory) and
       the column run: at iteration 300 the column run's u and v residuals at most half the profile run's, and
       fewer iterations to the verdict. Measured: 1237 iterations (60 under ctest, machine shared s on two ranks) against 3048 (1231 under ctest), residual
       ratios at iteration 300 of 0.084 (u) and 0.060 (v), every ERF check passing. Label ``regression``;
       double only.
   * - ``rt_steady_hill3d``
     - ERF's Neutral_Hill_3D (``inputs_steady_hill3d``: 64 × 64 × 20, the hill of revolution, damping
       from 550 m, 16 sweeps, cap 4) to the verdict with ``check_hill3d.py --physics`` (the 3-D checks, speed-up
       against 1.6 h/L). Measured: 2092 iterations (148 s on four ranks), every check passing. Label ``nightly``
       (minutes); double only.
   * - ``rt_steady_hill3d_neutral_zjacobi``, ``rt_steady_hill3d_neutral_lusgs``, ``rt_steady_hill3d_neutral_matches``
     - A truly neutral Hill_3D (``inputs_steady_hill3d_neutral``: 32 × 32 × 16 over the same domain, θ =
       300 K, ``forcing.gravity = 0``, damping of w alone, four boxes) run to near round-off (every residual below 1e-11
       and eight orders down; ten orders put the θ and p residuals, which carry only the divergence, under the
       rounding floor once ``gmres_mlmg`` became the default) with 16 z-line Jacobi sweeps at a cap of 4 and with the box-local line LU-SGS sweeps
       (``inputs_steady_hill3d_neutral_lusgs``: zero ghosts, overlap 2, the diagonal ordering, two sweeps, cap 16);
       ``Tests/Canonical_RANS/Truly_Neutral_Hill_3D/check_line_match.py`` compares the two plotfiles at 1e-8 of each field's maximum
       (``rt_steady_hill3d_neutral_clean`` empties the run directories first, so a stale plotfile cannot be
       compared, and the checker requires exactly one plotfile per run). The conventionally neutral deck has no
       fixed point to this tolerance (it keeps entraining the inversion), hence the variant. Measured: 3495
       iterations (70 s on four ranks) against 1792 (32 s); u, v, w, p and k agree to 9e-10 absolute (at most
       3e-8 of the tolerance scale). θ has a free constant (no θ condition, no gravity) that drifts by 1e-4 K at
       round-off, so θ is compared after subtracting each run's mean: 3e-11 K. Label ``nightly``; double only.
   * - ``rt_poisson_default_gmres_fft``, ``rt_poisson_default_matches``
     - The default ``solver.poisson_backend = gmres_mlmg`` against the previous default: the line LU-SGS deck of
       ``rt_steady_hill3d_neutral_lusgs`` run again with ``gmres_fft`` (``run/rt_steady_hill3d_neutral/lusgs_fft``) to round-off;
       ``check_line_match.py`` (the case folder) compares it with ``rt_steady_hill3d_neutral_lusgs`` (now ``gmres_mlmg``, the log
       pattern ``pressure gmres_mlmg`` checked) at 1e-8 of each field's maximum. Measured: both stop at 1789 iterations; u, v, w, p and k agree to 2.2e-13 absolute or better, θ − mean to 2.4e-11 K. Label
       ``nightly``; double only.
   * - ``rt_steady_hill3d_neutral_sequencing``, ``rt_steady_hill3d_neutral_sequencing_matches``
     - The deck of ``rt_steady_hill3d_neutral_zjacobi`` with ``mesh.hierarchy_levels = 3`` and
       ``steady.sequencing = 1`` (8 × 8 and 16 × 16 to the loose stop, then 32 × 32 to round-off), its log checked
       for the loose stops of levels 2 and 1, the strict verdict and the summary table; ``check_line_match.py``
       compares its plotfile with ``rt_steady_hill3d_neutral_zjacobi`` at 1e-8 of each field's maximum (the
       per-level plotfiles ``…_plt_level<l>_<step>`` are not counted). Measured: see :doc:`solution_algorithm` (grid
       sequencing). Label ``nightly``; double only.
   * - ``rt_steady_hill3d_sequencing``
     - ``inputs_steady_hill3d`` with three levels (16 × 16, 32 × 32, 64 × 64) and
       sequencing, with the checks of ``rt_steady_hill3d``: the steady verdict, a two-order drop of the u residual
       from the criterion's reference (the coarsest level's; ``check_steady_hill.py`` reads it from the checkpoint,
       which for a run without sequencing is the largest of the first 100 residuals, as before) and ERF's
       ``check_hill3d.py --physics``. Label ``nightly``; double only.
   * - ``rt_steady_askervein``
     - Askervein on the ERF grid of ``rt_init_askervein`` (``inputs_steady_askervein``: 30 × 30
       × 18, the log-law inflow, outflow, slip sides, no Coriolis, four sweeps, cap 2) to the steady verdict.
       Measured: 175 iterations (156 with ``bc.inflow_tke_length = unbounded``). ``init.alpha_v = 0.5``, so the log
       must show the projection on ``gmres_fft`` and the stage on ``gmres_mlmg``. The 5 % comparison with ERF's URANS
       on lines A and AA is not in this version (Kestrel). Labels ``regression;smoke``; double only.
   * - ``rt_precursor_askervein_cnbl``, ``rt_precursor_askervein_cnbl_fixed``, ``rt_precursor_askervein_cnbl_table``
     - The precursor column of the conventionally neutral Askervein case (``inputs_precursor_askervein_cnbl``: 4 × 4 ×
       26 on Askervein's 20 m × 1.08 levels, z0 = 0.03 m, the flat case's sounding, Coriolis at 57.2 N, geostrophic
       wind (17.028798, −7.380863) m/s) to the steady verdict. ``rt_precursor_askervein_cnbl_fixed`` runs the same
       deck to the step in the table's header (11087), with an unreachable drop criterion (the deck is written at
       configure time; a regenerated table carries its step). ``tools/make_inflow_table.py --check`` on that run: the
       plotfile step equal to the table's, the 10 m wind within 0.01° of +x, u* within 1e-3 of 0.6 m/s, and
       ``Tests/Real_Terrain/Askervein_CNBL/askervein_cnbl_inflow.txt`` within 1e-3 of the planar means (u, v relative
       to the geostrophic speed, θ in K, k relative to its maximum). **Why a fixed step.** Under the capping inversion
       θ has no fixed point: the inversion base cools by 1.35e-4 K per iteration as it rises, and the mixed layer
       warms by 1.5e-5 K. The verdict comes when the k residual leaves a limit cycle of about 600 iterations (the
       Poisson iteration count switching between 5 and 6), at an iteration that depends on the compiler's rounding:
       11087 with Apple clang and 11073 with gcc-16, 1.7e-3 K apart in θ. At step 11087 the two builds agree to 2e-4
       K. A build that left the cycle one cycle later would fail the table check, and the verdict run would show its
       step. Label ``regression``; double only.
   * - ``rt_steady_askervein_cnbl``
     - Askervein in a conventionally neutral boundary layer (``inputs_steady_askervein_cnbl``: 30 ×
       30 × 26 to 1599 m, the precursor's table as inflow on xlo and yhi, outflow on xhi and ylo, Rayleigh damping
       of w and θ above 1200 m, four sweeps, cap 2) to the steady verdict. Measured: 220 iterations (175 for
       ``rt_steady_askervein``); u residual down four orders. Label ``regression``; double only.
   * - ``rt_askervein_210_inflow_table``
     - ``Tests/Real_Terrain/Askervein_210/askervein_210_inflow.txt`` is the table ``tools/erf_inflow_table.py`` makes
       from ERF's Askervein inflow (``extern_erf/inflow_powerlaw.txt``, wind from 210°; w dropped, θ = 300 K, a z = 0
       row added), to 1e-9. Label ``regression;smoke``.
   * - ``rt_steady_askervein_210``, ``rt_askervein_210_erf_lines``, ``rt_askervein_210_lines_gold``,
       ``rt_askervein_210_level_history``
     - Askervein with ERF's setup (``inputs_steady_askervein_210``: 300 × 300 ×
       18, the 210° inflow through xlo and ylo, outflow through xhi and yhi, Rayleigh damping of w in the top 250 m)
       with the line LU-SGS configuration (``steady.sweep = lusgs`` with its default CFL cap, the default reference window,
       grid sequencing over two levels, the last three on the command line) to the steady verdict; measured: 95 level-0 iterations and 34 on level 1 (four
       ranks). The 10 m speeds on lines A and AA (``tools/askervein_lines.py check``) must be within 5 % of ERF's
       URANS average (``Gold/erf_urans_lines_10m.csv``) up to 300 m past the hilltop, and within 10 % at the four
       lee masts beyond (measured: 4.00 % and 9.21 %; the lee margin is 0.8 points, so a change of the closure or
       the advection that crosses it needs a review, not a wider tolerance). ``rt_askervein_210_lines_gold`` checks
       the same lines against ``Gold/askervein_210_lines.csv`` to 1e-4 (the ERF tolerances alone would pass the run
       with the inflow k of ``bc.inflow_tke_length = unbounded``). The level history (cells, stop reasons, iterations within
       10 %) and the stage switches must match ``Gold/askervein_210_level_history.json``. Label ``nightly``; double
       only.
   * - ``rt_steady_askervein_210_unbounded``, ``rt_askervein_210_unbounded_matches``
     - Old against new for ``bc.inflow_tke_length``: the line LU-SGS run of ``rt_steady_askervein_210`` with the inflow
       length of before (``unbounded``) reaches that run's lines (``Gold/askervein_210_lines_unbounded.csv``) to 1e-4
       (measured: identical to the six digits written, on the same machine). Label ``nightly``; double only.
   * - ``rt_steady_askervein_210_zjacobi``, ``rt_askervein_210_zjacobi_matches``
     - Old against new: the same deck with the z-line Jacobi configuration of ``rt_steady_askervein`` (``zjacobi``, four sweeps, cap 2, a
       100-iteration reference, one level; measured: 533 iterations, 565 with ``bc.inflow_tke_length = unbounded``) reaches the line LU-SGS run's lines to 0.1 % (measured:
       1e-4 %). Label ``nightly``; double only.
   * - ``drift_check_selftest``
     - The ``extern_erf`` drift check on a throw-away ERF repository (:doc:`extraction`).
   * - ``erf_plotfile_selftest``
     - ``extern_erf/erf_plotfile.py`` (a local change, :doc:`extraction`) on hand-made one-box plotfiles in double
       and float, little- and big-endian: every value read back exactly. The upstream reader fails the two float
       cases. Label ``unit``.
   * - ``refavg_snapshots``, ``refavg_build``, ``refavg_verify``, ``refavg_matches``
     - ``gust_refavg`` end to end: three copies of the smoke plotfile at 0, 5 and 10 s are archived on 2
       ranks with two inputs files; the archive verifies on 3 ranks with the expected metadata, AU and a zero
       half-window change; the average equals the smoke plotfile at zero tolerance. This checks the round trip;
       identical snapshots cannot check the weights, which ``Reference.TrapezoidWeights`` and
       ``Reference.AverageIsTheTrapezoidTimeAverage`` do.
   * - ``refavg_damage``, ``abort_refavg_damaged``, ``abort_refavg_inputs``, ``abort_refavg_verify_keys``
     - A copy of the archive with one bit flipped is refused, naming the file; bad tool inputs stop the tool with
       each input named; a build key given with ``refavg.verify`` is reported.
   * - ``versions_check``, ``versions_check_selftest``
     - ``external/VERSIONS`` is well formed and its ``amrex`` entry matches the submodule; the self-test shows a
       mismatch with a throw-away submodule is reported.
   * - ``kestrel_parse_*``, ``kestrel_dryrun``
     - Not registered on Windows. The Kestrel scripts parse; a dry run with stub ``module``/``sbatch``/``srun`` checks the Askervein job
       script (run length, plot period, no 2-D plotfiles, ranks) and archives a fake finished run, averaging only
       the 3-D plotfiles in the final window and archiving the inputs but not the run's outputs; a second case
       with the default prefix ``plt`` keeps its ``plt2d_`` plotfiles out.

In a build without MPI (``GUST_MPI=OFF``, as in the Windows CI job) every test runs on one rank, so
``unit_np2`` and ``unit_np4`` repeat ``unit_np1`` and ``smoke_restart_matches`` does not change the rank count.
MSVC does not trap floating-point exceptions, so ``amrex.fpe_trap_invalid=1`` has no effect there; the Linux GCC
jobs keep that check.

GoogleTest suite
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - GoogleTest
     - Check
   * - ``Reference.ArchiveRoundTripVerifies``, ``Reference.DamagedArchivesAreRefused``
     - A reference archive written, verified and loaded: fields bit-identical to the average, metadata and field
       statistics unchanged; a changed byte, a missing file, an extra file, edited AU, and data changed together
       with its checksum are each found and named, and loading refuses.
   * - ``Reference.TrapezoidWeights``, ``Reference.AverageIsTheTrapezoidTimeAverage``, ``Reference.HalfWindowChangeMeasuresDrift``, ``Reference.NaNSnapshotIsRefused``, ``Reference.AverageRejectsBadSnapshotLists``
     - The average equals a host-side recomputation bit for bit and, for fields linear in
       time, the mid-window value to 8 ε; trapezoid weights for uneven spacing; the half-window change is exact
       for a linear drift and zero for a steady field; a snapshot with a NaN, fewer than 3 snapshots and times out of
       order are refused with their messages; a trailing ``/`` is handled.
   * - ``Reference.AllocationUnits``, ``Reference.MetadataRangeChecks``
     - AU = nodes × hours × rate (10 CPU, 100 GPU); metadata range checks name their inputs.
   * - ``Sha256.*``
     - SHA-256 against the NIST example vectors, at every padding boundary and chunking;
       ``hex()`` can be called repeatedly.
   * - ``Mesh.VerticalLevels``, ``Mesh.ZLevelsFile``
     - Vertical levels: uniform; solved ratio (top exact, ratio 1.0421553943273616 from an independent exact-rational bisection);
       fixed ratio (ERF's Askervein 20 m × 1.08, top 749.0048747962335 m); z-levels file and its errors.
   * - ``Mesh.FlatTerrainIsTheCartesianMesh``
     - Flat terrain gives the Cartesian levels: bit for bit for STF and Sullivan, within 2 ulp for BTF.
   * - ``Mesh.FreeStreamIsPreserved``
     - Free-stream preservation on round Witch-of-Agnesi and Gaussian hills, all three mappings: divergence below
       1e-12 of :math:`(|U|+|V|)/\Delta x`. The cancellation is algebraic for any node heights, so this test does
       not detect wrong heights; ``Mesh.MetricsEqualErf`` does.
   * - ``Mesh.VolumesSumToTheDomain``
     - Cell volumes sum to the column volumes above the ground to 1e-12, and to a 1 m quadrature of the domain
       above the hill to 1e-4 (the bilinear ground differs from the continuous hill by 2.4e-5 of the domain).
   * - ``Mesh.SmoothSurfaceKernel``, ``Mesh.SmoothSurfaceNonPeriodic``
     - Smoothing: a constant is unchanged; a Fourier mode is damped by the analytic factor per pass; β ≥ 1/4 is
       refused. On a non-periodic domain, edge nodes use themselves as the outside neighbour and
       ghost nodes take the nearest domain node.
   * - ``Mesh.TerrainInterpolation``
     - IDW returns input points and, with 4 neighbours, is exact on a plane at a symmetric query (with the
       default 6 the two extra neighbours break the symmetry and IDW is not exact on a plane); ERF's bilinear form
       is exact on planes anywhere, with ERF's edge fallback; ``terrain.idw_neighbors`` above the number of points is
       refused.
   * - ``Mesh.MetricsEqualErf``
     - ``z_phys``, ``detJ`` and nodal heights equal ERF's own output for the same 2-D hill deck (BTF, STF,
       Sullivan; STF and Sullivan bit for bit on macOS, BTF to rounding because GUST keeps the ground layer as
       given, all within 1e-13 elsewhere); STF capped at 2 of 16 levels, as ERF reports.
   * - ``Mesh.IdwBucketsEqualBruteForce``, ``Mesh.TerrainAtMapCoordinates``, ``Mesh.FileTerrainSeamIsPeriodic``, ``Mesh.GroundLayerIsTheGivenSurface``, ``Mesh.TerrainThatFillsTheDomainIsRefused``, ``Mesh.QualitySlopeOfDiagonalRidge``, ``Mesh.RestartWithOtherTerrainIsRefused``, ``Mesh.QualityReportFlagsSteepTerrain``
     - The IDW bucket search picks the brute-force neighbours; results agree to 1e-13 (random
       cloud with far queries, a regular grid with exact distance ties, a tie across a bucket edge, a transect
       queried from a 2-D domain, one point); UTM-sized coordinates are read and interpolated in double; file terrain on a periodic domain has one
       seam height; every mapping keeps the ground layer as given; ground at the top and inverted cells stop the
       run; the slope of a diagonal ridge is reported in full; a restart with changed terrain is refused; the report
       of a steep ridge carries both warnings.
   * - ``ConvectionMMS.UnlimitedMusclIsSecondOrder``, ``ConvectionMMS.LimitedMusclOrder``, ``ConvectionMMS.UpwindIsFirstOrder``
     - Convection on a manufactured flow over a 120 m cosine hill (BTF, geometric layers, 16 × 12 × 8 to 64 × 48 × 32
       cells), ghost nodes continuing the layers: the unlimited MUSCL reconstruction (van Albada slope with
       :math:`\epsilon` above every difference) is second order for q (conservative and bounded forms) and for the u, v
       and w control volumes, L2 order ≥ 1.9 on both pairs and max-norm order ≥ 1.9 on the finest. With ``van_leer``
       and ``van_albada``: L1 order ≥ 1.5 on the finest pair and an L1 error at most a quarter of first-order upwind's
       (a limiter is first order at smooth extrema). First-order upwind converges at order 1.0 ± 0.1.
   * - ``DiffusionMMS.ScalarIsSecondOrder``, ``DiffusionMMS.MomentumIsSecondOrder``
     - Scalar diffusion with a variable :math:`\gamma`, and momentum diffusion with the full stress tensor and a
       variable :math:`\nu`, on the same meshes: L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on the finest,
       for q, u, v and w.
   * - ``ConvectionProperties.BoundedUpwindKeepsKPositive``, ``ConvectionProperties.BoundedFormOfConstantIsZero``
     - Random :math:`k \ge 0` (a third of the points zero) and random transport fluxes of ±50 m/s: the step
       :math:`k - C(k)/a_P` with the bounded upwind form leaves every :math:`k \ge 0`; the conservative upwind form
       goes below −0.1 in the same step. The bounded form of a constant is zero for every limiter
       and location.
   * - ``ImplicitOperator.DiagonallyDominantForRandomWindAndNu``
     - Random transport (±20 m/s) and diffusivity (0 to 50 m²/s, a fifth zero) over the hill: every neighbour
       weight of :math:`A_1` is ≥ 0 and :math:`a_P` equals their sum (1e-14; single 1e-6), for cells and the three
       velocity components.
   * - ``SmoothLimiter.DerivativeIsContinuous``, ``SmoothLimiter.CloseToVanLeerOnSmoothData``, ``ConvectionOptions.InvalidEpsIsRejected``
     - The van Albada slope's derivative changes by less than 0.1 across every point where van Leer switches
       (:math:`a = 0`, :math:`b = 0`, :math:`a = -b`, the origin), in three directions; van Leer's jumps by more than 1.
       On the manufactured flow (64 × 48 × 32) the two limiters differ by less than 1 % (L2) for q, u, v and w.
       Invalid ``smooth_eps`` (0, negative, NaN, infinite) is refused.
   * - ``ImplicitOperator.ConvectionPartIsBoundedUpwind``, ``ImplicitOperator.DiffusionPartIsTheFlatOperator``, ``ImplicitOperator.WeightsAreOneDimensionalFluxes``, ``ImplicitOperator.IndependentOfBoxLayout``
     - :math:`A_1` is the explicit operators' first-order part: bounded upwind convection equals
       :math:`A_1 q` for random q and transport (1e-12; single 1e-5), and on a flat mesh :math:`-D(q) = A_1 q` for the
       scalar and :math:`-D_c(u) = A_1 u_c` for each velocity component with random :math:`\nu`. Over the hill, every weight
       equals the explicit face flux of a field that varies along one index direction only (this pins the
       :math:`1 + \delta_{cd}` and :math:`|\mathbf{n}|^2 + n_c^2` factors on sloping faces). MUSCL convection,
       both diffusions and :math:`A_1` agree to round-off (1e-13) on 16² and 8² column boxes, and boxes that share a
       face point agree there (not bit for bit: GCC on arm64 contracts multiply-adds inside vectorised loops only).
   * - ``MassFluxMMS.DivergenceAndGradientAreSecondOrder``, ``MassFlux.GradientIsMinusTheAdjointOfTheDivergence``, ``MassFlux.GradientIgnoresTheGhostBelowTheGround``
     - The divergence (through :math:`\Omega`) and the three components of :math:`G p` over the 120 m hill (BTF,
       geometric layers, 16 × 12 × 8 to 64 × 48 × 32 cells): L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on
       the finest (single precision 1.8). :math:`\langle Gp, U\rangle = -\langle p, DU\rangle` with control-volume
       weights to 1e-12 of either side (7e-12 against 3.3e4), for random p and a random velocity that vanishes near
       the ground and top. :math:`G` does not read the ghost cells below the ground.
   * - ``MassFlux.DivergenceIsThePoissonLayers``, ``MassFlux.OmegaAndWRoundTrip``
     - The continuity residual is bit for bit ``face_divergence`` of :math:`(u, v, \Omega)`; after a projection
       (Poisson fluxes added to :math:`u, v, \Omega`, w from ``w_from_omega``) the divergence falls from 1.19 to
       2.4e-11, the solve's residual. ``omega_from_w`` and ``w_from_omega`` are inverses.
   * - ``LineJacobian.MatchesFiniteDifferenceJacobian``
     - ``line_jacobian`` for u, v, w, θ and k in x, y and z against a finite-difference Jacobian of the first-order
       operator (transport frozen), per coefficient: relative difference < 1e-6 (1e-13 measured; single 1e-3), over
       the hill with convection, and on a flat mesh with convection and a random diffusivity.
   * - ``Residual.RestOverSteepTerrainIsBalanced``
     - Rest over a 300 m hill (slopes up to 0.94), θ = θ_ref(z) stratified, Coriolis on: the momentum, continuity and
       k residuals are below 1e-14. With a constant reference the hydrostatic :math:`p(z)` balances only to
       truncation (R_u 1.6e-3 m/s², 0.5 % of the buoyancy; R_w exact on interior faces).
   * - ``Jvp.*``
     - ``jvp`` on a laterally uniform column (affine residual) equals :math:`R(q+v) - R(q)` to 1e-5 (single 1e-2).
       On the hill with MUSCL convection, against a central difference, L1 error below 1e-5
       (2.7e-6) while a step 1e5 times larger misses (1.8e-4); the result is the same without ``R0``; the step is the
       same on 16² and 8² column boxes and scales with each field (k and its scale times 100 leave it unchanged).
   * - ``Residual.SourcesHaveTheirClosedForms``, ``Residual.IsTheSumOfItsParts``, ``Residual.InvalidInputsAreRefused``
     - Sources in closed form: Coriolis with the :math:`\cos\varphi` terms and geostrophic forcing
       on a uniform wind (:math:`f = 7.27221\times10^{-5}` 1/s at 30°, T = 86400 s), buoyancy of a uniform excess;
       the geostrophic wind is balanced. The residual equals its operators applied one by one (momentum, continuity,
       θ, k) on the hill with a general state. Invalid options (latitude, period, geostrophic forcing at latitude 0 or with
       the forcing off, θ₀, gravity, profiles that do not increase) and missing inputs stop with a message.
   * - ``Mesh.SharedNodesAgree``
     - Node heights shared by two boxes (or periodic images) are equal, for BTF, STF and Sullivan,
       periodic and not (under gcc each box's copy differed by up to 1e-13 m before ``Mesh`` synchronised them).
   * - ``Mesh.GhostNodesStandForTheirImages``
     - Every ghost node in x and y, at every k (the rows below the ground and above the top included), equals the node
       it stands for, bit for bit: its periodic image, or the edge node of its row at a non-periodic face. BTF, STF
       and Sullivan; periodic in x only, y only, both and neither, and a layout split in z; on 1, 2 and 4 ranks; every
       copy of a domain node agrees. Before ``fill_ghost_columns`` the corners outside a wall and a periodic side held
       the clamped corner column (up to 8 m off) and the row below the ground of periodic ghost columns was
       extrapolated from the clamped columns (STF, Sullivan).
   * - ``SimilarityFunctions.MatchTheReferenceTable``
     - ERF's similarity functions (the copied ``ERF_MOSTUtils.H``, byte-identical by the drift check): Businger–Dyer
       ψ and φ and Jimenez ψ against a Python table of the same formulas, z/L from −5 to 5, to 1e-14 (single 1e-6);
       :math:`d\psi/d\zeta = (1-\phi)/\zeta`.
   * - ``SurfaceLayer.ReferenceIntervalIsHalfOpen``
     - Reference-height lookup on the ERF#3982 column (dz₀ = 10 m, ratio 1.1, 40 layers): centres belong to the interval
       above, faces (10 m included) to the interval that holds them, the last centre to the last interval; heights below
       the first or above the last centre are refused.
   * - ``SurfaceLayer.LogLawFrictionVelocityOnAStretchedColumn``, ``SurfaceLayer.MostFrictionVelocityWithAHeatFlux``, ``SurfaceLayer.RoughWallFollowsTheLogLaw``, ``SurfaceLayer.DenominatorFloorKeyIsChecked``
     - u* from a log-law state (u* = 0.45, z₀ = 0.1 m) on that column: within 2e-5 at the first centre, at 10 m, at 55.3 m
       and on a centre. With a heat flux, u* matches an independent bisection to 1e-10 for unstable,
       stable and low-wind convective cases, and is positive. On rough ground (z₀ = 0.5 m, 1 m first layer,
       ln((d + z₀)/z₀) = ln 2): u* from the log law at the first centre and at 0.8 m with
       ``bc.surface_denominator_floor = log_law``, 0.45 ln 2 with ``one``; with a heat flux both floors match the
       bisection with their floor to 1e-10. An unknown floor stops the run, naming the key.
   * - ``SurfaceLayer.WallKAndWallViscosity``
     - :math:`k_w = u_*^2/C_{\mu0}^2` exactly (neutral); Axell's buoyant form with a heat flux; the wall eddy viscosity
       equals :math:`C_{\mu0}\sqrt{k_w}\,\kappa(d+z_0)`.
   * - ``SurfaceLayer.ColumnsOverTheHill``, ``SurfaceLayer.WallModelInTheResidual``, ``SurfaceLayer.WallRowsOfTheLineJacobian``, ``SurfaceLayer.WallKRowAtLowWind``, ``SurfaceLayer.InvalidOptionsAreRefused``
     - Over the hill: the columns' distances and areas from their nodes; the wall stress along the
       ground with magnitude :math:`u_*^2|A|`. In the residual, R_u, R_v and R_θ change only in the first cells, by the
       wall flux minus the interior flux over the cell volume; the k rows are the wall-k rows (both modes). The
       ``line_jacobian`` rows lose the ground neighbour and gain the wall terms. At low wind the wall-k row uses
       u* ≥ 0.01 m/s. Invalid options and a reference height outside the columns stop with a message.
   * - ``BoundaryTypes.PoissonTypesMatchErfGetFftBc``
     - The Poisson boundary types of every face type on every face it may take, periodic directions included, equal
       ERF's ``get_fft_bc`` (the copied ``ERF_SolverUtils.H``).
   * - ``BoundaryFill.InflowFollowsTheGroundHeight``, ``BoundaryFill.RulesPerFaceType``, ``InflowProfile.*``
     - On a ridge whose height varies along the inflow face: every inflow ghost point (both layers, u, v, θ, k) takes the
       profile at its own height above the ground, the inlet wind varies along the face by more than 0.05 m/s, and the
       outflow face points of u are not Dirichlet rows. The ghost values of every rule for every
       face type (even, odd, value, profile, gradient, the Dirichlet face points, the top's θ gradient); the analytic
       profiles and their equilibrium k and :math:`\nu_t`, and a file table with and without a k column.
   * - ``BoundaryFill.GroundGradientIsTheWallLaw``, ``BoundaryImplicit.FoldsMatchTheHomogeneousFill``, ``BoundaryResidual.DirichletRowsAndIdentityLines``, ``BoundaryDamping.*``, ``BoundaryInputs.KeysAreReadAndChecked``
     - Over a 120 m hill with a stable heat flux, the ground ghosts of u, v and θ give MOST's
       gradient at the first centre (u* φ_m/(κ(d₀ + z₀)) along the stress, −q φ_h/(κ u* (d₀ + z₀))) in both layers,
       on the two lateral ghost columns too (the periodic images of the neighbouring boxes; this catches a gradient
       slab with too few ghost columns under ``init_snan``); k is even; the homogeneous fill gives zero gradients
       and keeps even rules. The implicit folds: :math:`A_1` on the interior values with the folds equals the unfolded
       :math:`A_1` on the homogeneously filled state, for u, v, w, θ and k, every lateral type and a geostrophic top
       (1e-13; single 1e-5); Dirichlet face points carry identity rows. The residual's Dirichlet rows are
       :math:`q - V`; ``line_jacobian`` gives the identity there. The Rayleigh and sponge coefficients have their
       shapes (the :math:`\sin^2` ramp, zero below :math:`z_d`; the quadratic sponge, zero at its inner edge and on
       faces that are not inflow or outflow) in the residual and on the implicit diagonal. Every ``bc`` key is read
       into its field, and each out-of-range, ineffective (the other kinds' inflow keys, surface keys without a
       surface-layer ground, a layer without a coefficient) or inconsistent key is named.
   * - ``Closure.KernelsMatchErfBitForBit``
     - ``closure_point`` equals ERF's composition of the copied AL01 relations (``ComputeTurbulentViscosityRANS``,
       ρ = 1) bit for bit on 20 000 random (k, N², d): L, ε, R_t, c_μ, c_μ', ν_t, ν_t', across the neutral, stable,
       unstable, smoothed and floored cases (each exercised).
   * - ``Closure.MatchesKynemaKLAxellWhereTheyCoincide``
     - Kynema's ``KLAxell`` (transcribed) where the formulations coincide (neutral, z = d, λ = 30 m, no z₀, ρ = 1):
       L, ν_t, ε and P_s to 1e-12.
   * - ``Closure.SmoothingLimitsAndLargeRichardsonNumbers``
     - ``smooth_Rt``: identity above R_t^crit, continuous at it, within (R_t^min, R_t^crit] and monotone down to
       R_t = −1e30; the unstable length finite and below its bound for k down to 1e-16 under N² to −1e3; the same
       decision as ERF's.
   * - ``Closure.PointImplicitUpdateStaysPositive``
     - The point-implicit update k + S_k/(1/Δτ + ∂S_k/∂k) with the Patankar diagonal stays ≥ 0 (> 0 with production)
       for 1/Δτ = 1, 1e-3 and 0 on 5000 random states; without the diagonal a decaying state goes negative.
   * - ``Closure.SourceJacobianMatchesFiniteDifferences``
     - The dual-number dS_k/dk and dν_t/dk against central finite differences on 3000 random states (branches held):
       relative 1e-6 (single 1e-2); the dual value part equals the Real evaluation to a few ulp.
   * - ``Closure.InputsAreReadAndChecked``
     - Every ``turb`` key is read into its field and the constants follow (the floor, C_μ0³); each invalid value is
       refused with its key named, and the constructor refuses it before forming the constants.
   * - ``Closure.FaceMeanIsLogarithmicAcrossVerticalPairsOnly``, ``Closure.LogarithmicFaceMeanReproducesTheLogLawFlux``, ``Closure.ResidualWithClosureIsTheResidualWithItsFields``, ``Closure.RelaxationCheckpointAndFrozenSwitches``
     - The logarithmic mean against (a − b)/ln(a/b) on random pairs (1e-14), equal values, zero,
       continuity across its series switch; on a log-law column the logarithmic face mean gives the wall flux u*² at
       every interior u face (1e-12) where the arithmetic mean misses the first by 2–12 %; the implicit operator with
       the same mean satisfies A₁u = −D(u). Over a 120 m hill the residual with the closure equals the residual with
       the closure's own fields given minus the sources bit for bit, the k line Jacobian's diagonal differs by dSdk
       exactly, every field is finite with positive diffusivities, misuse (fields given with a closure, no wall
       distance) is refused, and ``height_above_ground`` is the node-mean height above the column's ground.
       Relaxation: ν_t = 0.6 new + 0.4 stored exactly, unrelaxed without history, and the history survives a
       checkpoint round trip onto other boxes (identical relaxed evaluation); frozen switches hold the floor decision
       after k drops below it.
   * - ``Poisson.ManufacturedSecondOrderAllBCs``
     - Manufactured :math:`\phi` on a stretched STF mesh over a 60 m hill: every lateral BC (periodic, Neumann,
       Dirichlet) with a Dirichlet or Neumann top, 2-D (y periodic) at three resolutions, and a 3-D case (Neumann in x,
       Dirichlet in y, two resolutions): order ≥ 1.9 (single precision: 1.8 coarse, 1.6 fine), residual below 10 ×
       tolerance. The solution has zero gradient at the ground and level side faces, where ERF's operator is exact
       for the mirror ghost cells (next row).
   * - ``Poisson.GroundNeumannIsFirstOrderOnSlopes``
     - With a tangential gradient at sloped ground: order 2 on flat ground (≥ 1.95), between 0.8 and
       1.2 over a 120 m hill (64 → 128 → 256 cells), and a fine-mesh error 3 to 6 times that over a 60 m hill (slope²).
       Double precision only. Pins ERF's ground treatment (:doc:`solution_algorithm`).
   * - ``Poisson.VariableCoefficientsSecondOrder``
     - The anisotropic operator with :math:`\beta_v/\beta_h` = 1, 0.5, 0.1 (``z_cross`` set) on flat ground and over
       the hill, and a spatially varying :math:`\beta` over the hill: order ≥ 1.9 (single: 1.8, 1.6). Without
       ``z_cross`` the hill error is more than 5 times larger.
   * - ``Poisson.UnitCoefficientsAreThePlainOperator``
     - :math:`\beta = 1` gives the no-coefficient path bit for bit (iterations, solution, fluxes); :math:`\beta = 2.5`
       on every face gives :math:`\phi/2.5` and the same fluxes (to 1e-9; single 1e-4). Both paths are this layer's;
       there is no separate plain projection to compare with.
   * - ``Poisson.NeumannFacesCarryNoFlux``
     - On a ramp (slope 0.15) crossing the Neumann faces, their fluxes and the ground's are exactly zero, and
       :math:`\sum J\,A\phi` equals the flux out of the Dirichlet top (1e-12; single 1e-4).
   * - ``Poisson.OperatorMatchesFluxDivergence``
     - Operator and fluxes agree cell by cell, :math:`A\phi = -\nabla_h\cdot F` (1e-12 of
       :math:`\max|A\phi|`; single 1e-4), with varying coefficients and ``z_cross`` on the ramp, Dirichlet and Neumann
       top.
   * - ``Poisson.ProjectedFieldIsDivergenceFree``
     - Projection: the divergence of :math:`u + F` is below 10 × tolerance of the divergence of :math:`u`, and equals
       the residual of the solve.
   * - ``Poisson.MirrorSymmetry``
     - Symmetric hill and right-hand side: :math:`\phi` mirror-symmetric to 1e-12 of :math:`\max|\phi|` (single 1e-5),
       Neumann and Dirichlet sides.
   * - ``Poisson.PreconditionerCutsIterations``, ``Poisson.SingularProblemRemovesTheMean``, ``Poisson.ChecksNameTheInput``
     - Preconditioner, with bounds 10 % above the counts measured: 19 iterations isotropic (more than
       10 times as many without it); 16 for :math:`\beta_v/\beta_h = 0.1` told the ratio (at most 3 more than
       isotropic; not told, more than 1.5 times as many); 29 for :math:`\beta` varying over 100 times. A singular
       problem reports and removes the mean. Option and BC checks, coefficients that are not positive, MultiFabs on
       another layout and a solve that does not converge each stop with the input named.
       ``tools/poisson_ms_check.py`` (sympy, not run by CI) checks the manufactured solutions' boundary conditions.
   * - ``PoissonB.ManufacturedSecondOrder``
     - Back end B (``mlmg`` and ``gmres_mlmg``): every lateral BC × Dirichlet/Neumann top × Neumann (zero-gradient
       solution) / Dirichlet ground, 2-D (y hidden) at three resolutions (two in single precision), and a 3-D case:
       order ≥ 1.9 (single: 1.8, 1.6).
   * - ``PoissonB.DirichletGroundSecondOrderBackendA``
     - Back end A with a Dirichlet ground: order ≥ 1.9 for every lateral BC and top.
   * - ``PoissonB.BackendsAgreeWithVariableCoefficients``
     - Spatially varying coefficients, Neumann and Dirichlet ground, 3-D with non-periodic sides (Neumann x,
       Dirichlet y and top): both back ends solved to 1e-12 give :math:`\phi` within 1e-9 of max :math:`|\phi|` (single:
       1e-5 and 1e-4), and the fluxes satisfy :math:`\nabla_h\cdot F = -\mathrm{rhs}` to the residual.
   * - ``PoissonB.UnsupportedCombinationsAreRefused``
     - ``z_cross`` with back end B, ``z_cross`` with a Dirichlet ground, ``apply`` with back end B,
       an unknown back end, ``solver.poisson_precond`` with back end B, ``solver.poisson_restart`` with ``mlmg`` and a
       non-converging MLMG solve each stop with the input named.
   * - ``WallDistance.FlatGroundIsTheHeight``, ``WallDistance.FlatGroundSecondOrder``, ``WallDistance.RidgeWithinErfTolerances``, ``WallDistance.Hill3dWithinErfTolerances``, ``PoissonBWallDistance.*``
     - Tucker distance against the exact distance with ERF's tolerances on ERF's decks: flat ground (max relative
       ≤ 3e-3, ≤ 0.1 m below 100 m; single 1e-2 and 0.5 m), 2-D ridge and 3-D hill (max relative ≤ 0.15, mean ≤ 0.03,
       ≤ 0.2 cells below 100 m); the 3-D hill with periodic and Neumann sides; both back ends (back end B on 64-cell
       boxes, under CI's trap). Single precision solves to 1e-4. On flat ground the error
       converges at second order (32, 64, 128 layers; double precision).
   * - ``Mesh.NoSentinelNodesLeft``, ``WallDistance.GhostFillFollowsTheBoundaryConditions``
     - No node of ``z_phys_nd`` (ghosts included) holds ERF's ``bogus_large_value``, for BTF, STF and
       Sullivan with periodic and non-periodic sides; ``fill_ghosts`` fills every ghost cell, corners and between-box
       ghosts included, with odd copies below a Dirichlet ground and even copies at Neumann faces.
   * - ``Initializer.EqualAlphasAreThePlainProjection``
     - :math:`\alpha_h = \alpha_v = 1`: the initial field equals the Poisson layer's plain projection of the same first
       guess (the solve without coefficients, applied the same way) bit for bit, :math:`\lambda` included;
       :math:`\alpha_h = \alpha_v = 2` gives the same velocity within 1e-9 of max :math:`|u|` and :math:`\lambda/4`
       (single: 1e-4); the scalar state is untouched by the correction. Scaling the whole z flux by
       :math:`\alpha_v^2` instead of the cross part by :math:`\alpha_h^2` gives a different field; with periodic sides
       the hill case is divergence-free as a singular problem and the seam ghost of :math:`u` is the corrected
       periodic image; a :math:`\lambda` handed in with garbage is zeroed first; the divergence after is measured,
       not assumed.
   * - ``Initializer.*``
     - The divergence after the correction is below 10 × ``solver.poisson_reltol`` times the divergence before, in the
       solver's own operator (``velocity_divergence``), in every initializer test.
   * - ``Initializer.FlatTerrainLeavesTheProfile``
     - Flat terrain, open and periodic sides: the profile is divergence-free, so :math:`\lambda`, :math:`w` and the
       divergence stay below 1e-12 and :math:`u` equals the log law at every x face's height (single: 1e-4, 1e-5).
   * - ``Initializer.D4SymmetryOnTheRoundHill``
     - D4 symmetry: on a square domain with a Gaussian hill, the wind along ±x and ±y (its own inflow and outflow
       faces each time) gives fields that are the eight transforms of each other: the cell-centred velocity vector
       maps with the transform's matrix, :math:`w` and :math:`\lambda` as scalars, within 1e-12 of the maxima (single
       1e-5), with the solve converged to round-off (``reltol`` 3e-15).
   * - ``Initializer.SummitSpeedUpDecreasesWithAlphaRatio``
     - The speed-up above the summit (four central cells a quarter of the depth up, over the profile at their
       height) falls by more than 0.01 at each step of :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (1.015, 0.998,
       0.975, 0.943, 0.913 on the 1 km test domain, whose open faces are 500 m from the summit), and
       :math:`\alpha_h = 2, \alpha_v = 1` gives the same value as :math:`\alpha_v = 0.5` (the
       ratio alone matters; 1e-9). The preconditioner is told the ratio: at most 21 GMRES iterations at every ratio
       (13 to 19 measured).
   * - ``Initializer.FirstGuessByHeightAboveGround``
     - The first guess over sloping terrain: :math:`u` at every x face, the two ghost layers included, and :math:`v`
       at every y face equal the profile at the face's height above its column's ground with the 3 m floor (the
       floor acts on the ground faces and the ghosts below); :math:`\theta` equals the sounding by absolute height;
       :math:`k` the profile's value with the decay above 200 m, and :math:`u_*^2/C_{\mu 0}^2` everywhere without it;
       :math:`p' = 0`; :math:`w_0 = 0` except on the ground faces, which carry no flux (:math:`\Omega_0 = 0` to
       round-off); every ghost cell of the state is written (no NaN from the allocation); the O'Brien pass refuses a
       layout whose boxes do not span the column.
   * - ``Initializer.OBrienColumnCorrection``
     - O'Brien: :math:`\Omega_0` vanishes at the ground and the top exactly, the divergence left in
       every cell is the weight difference :math:`f_{k+1}^2 a_{z,k+1} - f_k^2 a_{z,k}` (:math:`f` the height fraction,
       which the test shows differs from the index fraction on the stretched mesh) times the column's residual
       (1e-11), the column sum equals the residual flux, the residual the function returns equals the largest column's;
       through the initializer the correction starts from a different divergence than without O'Brien (not a lower
       one: the terrain-following first guess is already consistent at the ground) and gives a different :math:`w`.
   * - ``Initializer.InputsAndPoissonTypes``
     - Inputs: :math:`z_{agl,\min} < 0` or NaN, :math:`\alpha_h \le 0` or infinite,
       :math:`\alpha_v < 0` or NaN, :math:`z_i < 0`, and any key other than ``first_guess`` given without
       ``first_guess = profile`` are refused with the key named, by the checks and by both entry points;
       ``first_guess = none`` leaves the state untouched; the reader takes every key and refuses an unknown
       ``first_guess``, ``obrien = 2``, a repeated key and the keys without the first guess; the Poisson types of the
       correction are Dirichlet on inflow and outflow faces, Neumann elsewhere, periodic where the mesh is; a face set
       that does not fit the mesh is refused.
   * - ``Poisson.SplitOperatorKeepsTheMirrorSymmetries``
     - The split anisotropic operator (``z_cross`` set, :math:`\beta_v = 0.36`) on a Gaussian hill varying in x and y
       keeps the mesh's x and y mirror symmetries within 1e-12 of max :math:`|\phi|`; the cross kernel's first
       transcription broke the y symmetry at the 1e-4 level.
   * - ``PseudoTime.LocalStepMatchesTheAnalyticLimits``
     - The local pseudo-time step on known fields: uniform :math:`(u, v, w)` and :math:`\nu` on a flat uniform mesh
       give :math:`\mathrm{CFL} / (|u|/\Delta x + |v|/\Delta y + |w|/\Delta z + 2\nu(1/\Delta x^2 + 1/\Delta y^2))` in
       every cell to 1e-14 (single 1e-6; the hill cases 1e-12 and 1e-5), without :math:`\nu` the diffusive term is gone, a field at rest without
       diffusion gets the floored step, :math:`u` varying in :math:`x` gives the face mean per cell, a
       terrain-following field (:math:`\Omega = 0`) over a 120 m hill has no vertical term and a plain :math:`w` over
       it the term :math:`|w| / (J \Delta\zeta)` cell by cell; a non-positive or NaN CFL and a face layout are
       refused.
   * - ``Driver.SerRampFollowsTheLawAndItsCaps``
     - The SER ramp on a geometric residual (rate 1/4): CFL 1, 1, 1, 4, 16, 64, 64, 64 (``cfl0`` through the two
       reference iterations, then the law capped by the growth 4 and ``cfl_max`` 64), :math:`p = 2` without caps
       256, 4096, 65536, a growth cap of 1.5 binding every iteration, :math:`p = 0` holding the CFL, a residual of
       exactly 0 reaching the cap through the growth cap; the stage's :math:`\Delta\tau` follows the CFL.
   * - ``Driver.SwitchRulesHaveHysteresis``
     - Switch rules on scripted residuals (:math:`W` = 5, threshold 0.01, back growth 10): a residual alternating
       about the threshold never switches; five iterations below switch forward at iteration 5 with the reason
       naming ``steady.switch_residual``, the residuals and the stage's ``init`` recorded; a 20× growth of the
       second stage's residual switches back at iteration 11 (the window compared must lie within the second stage,
       entered at 6); every pair of switches is at
       least :math:`W` apart; the records carry each iteration's stage; one stage never switches.
   * - ``Driver.HistoryIsRestartableBitForBit``
     - The history written to ``<dir>/Driver`` and read back is equal field by field (records, switches with their
       reasons, accounts, the criterion's history); a two-stage run interrupted after 6 iterations and continued by a
       new driver from the file on the saved state gives the same records (stage, CFL, relative residual, residuals,
       monitors), switches and verdict iteration as the uninterrupted run; the AU of Kestrel CPU on 2 nodes equal
       20 × wall hours; the history is identical on every rank; a truncated file, a missing file and a history for
       another number of stages are refused.
   * - ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
     - A residual that becomes infinite or NaN stops the run as ``Diverged`` with no arithmetic on
       it (the traps); ``max_iterations`` bounds
       a run that never settles and a second ``run`` does nothing; no stages, a null stage and stages that disagree
       on residual names are refused (a null first stage too); ``restore`` refuses a NaN CFL, a switch state or
       iteration count that does not fit the records, and brings the CFL within changed bounds; every ``steady.*``
       range check names its key; the reader takes every key,
       its defaults, and refuses a repeated key, an unknown machine and a bad value.
   * - ``Predictor.APIsTheDiagonalPlusThePseudoTimeTerm``
     - :math:`a_P` on every x-, y- and z-face equals :math:`1/\Delta\tau` plus the diagonal of ``line_jacobian``
       exactly, with a random :math:`\Delta\tau` field over a 60 m hill with random diffusivities (both diagonals
       come from ``implicit_rows``, whose coefficients ``LineJacobian.MatchesFiniteDifferenceJacobian`` pins: the statement is the :math:`1/\Delta\tau`
       addition and the face step); the diagonal exceeds :math:`1/\Delta\tau`; the two copies of the periodic seam
       face have identical rows and identical corrections; the pressure and scalars are untouched and :math:`u`
       changes by exactly the z-line solve's correction.
   * - ``Predictor.ZLineJacobiConvergesToTheDenseSolution``
     - The z-line solve (the vertical upwind weights of :math:`A_1` included) against a dense direct solve of the same
       system,
       assembled independently from the coefficient fields: on a single column (1 × 1 × 12, periodic sides, the
       lagged neighbours are the column itself) to 1e-12 after 60 sweeps (single 1e-4); the ghosts outside the domain
       stay 0.
   * - ``Predictor.UnderRelaxationIsAPseudoTimeStep``
     - :math:`\alpha_u` relaxation equals pseudo-time relaxation: with :math:`\Delta\tau = \alpha_u/(1-\alpha_u)/
       \mathrm{diag}(A_1)` on every face, :math:`a_P = \mathrm{diag}(A_1)/\alpha_u` to 1e-12 (single 1e-5) and the
       converged correction of u, v and w equals the dense solve of the relaxed system to 1e-12 (single 1e-4).
   * - ``Predictor.ZLineJacobiConvergesToTheDenseSolution``
     - z-line Jacobi on a 3-D linear advection-diffusion problem (4 × 4 × 6, periodic, :math:`\Delta\tau` = 20 s)
       converges to the dense direct solution: the error falls at every doubling of the sweeps (1, 2, 4, 8, 16) and
       is below 1e-10 after 400 sweeps, with the last sweep's change below 1e-10.
   * - ``Predictor.AConvergedStateIsAFixedPoint``
     - Fixed point: a uniform wind over flat ground with uniform pressure has a residual of exactly zero and the
       predictor changes nothing (:math:`|\Delta q| = 0`; with a zero right-hand side the solve is zero for any
       coefficients, so this follows from the residual); a perturbation of :math:`u` at one face is pulled back: one
       predictor step with :math:`\Delta\tau` = 10⁴ s and 50 sweeps lowers max :math:`|R_u|` below 0.75 of its value
       (0.56 measured; the pressure and the transport are frozen in the step).
   * - ``Predictor.FaceStepRefusalsAndRankIndependence``
     - ``face_dtau`` is the mean of the two cells, the one cell at a wall, the ground and the top,
       and the periodic image across the seam; one sweep's change is the whole correction; ``sweeps < 1``, a
       non-positive step and wrong layouts are refused.
   * - ``Predictor.DirichletRowsAreSetExactly``
     - With boundary conditions (inflow at xlo, outflow at xhi, slip ground and top), the Dirichlet
       rows (the inflow face of u, the ground and top faces of w) have :math:`a_P = 1`, :math:`d = 0`, and one
       predictor step sets them to the boundary value exactly for any :math:`\Delta\tau` (the inflow faces then have
       a zero residual row); every other row has :math:`a_P = 1/\Delta\tau + \mathrm{diag}(A_1)` and :math:`d = a_P -
       \sum a_{nb} > 0`.
   * - ``Predictor.IndependentOfTheBoxLayout``
     - The same state on one box and on 2 × 2 boxes (the lagged neighbours cross the box seams
       through the ghost exchange) gives the same corrected u, v, w, :math:`a_P`, :math:`d` and face step to 1e-12
       (single 1e-5).
   * - ``Predictor.WVnnCapAddsToTheWDiagonalAndD``
     - The w rows' VNN cap (``simplec.w_vnn_cap``, :doc:`solution_algorithm`) on the boundary-condition state
       of ``Predictor.DirichletRowsAreSetExactly`` with :math:`\Delta\tau` of 400–1600 s (an uncapped VNN above the cap): with :math:`V = 16` every
       interior w row has :math:`a_P = 1/\Delta\tau + (a_{ZM} + a_{ZP})/V + \mathrm{diag}(A_1)` and :math:`d` larger by
       :math:`(a_{ZM} + a_{ZP})/V` (1e-13; single 1e-5), the Dirichlet ground and top rows keep :math:`a_P = 1` and
       :math:`d = 0`, the effective VNN stays below :math:`V`, and the w correction is the z-line solve with the
       capped diagonal; u and v, and every field with the cap 0, are bit-identical to the call without the argument;
       a negative, NaN, infinite or subnormal cap is refused. On 1, 2 and 4 ranks (2 × 1 boxes).
   * - ``Predictor.WVnnCapSlopeMarginKeepsThePressureOperatorDefinite``
     - The slope-aware cap (``simplec.w_vnn_cap_slope_margin``) on a steep manufactured hill
       (150 m over 1000 m, :math:`h^2` up to 0.12) with VNN up to about 70, :math:`V = 16`, margin 10: on every
       interior w row :math:`a_P` and :math:`d` gain :math:`\min(K/V, \max(0, 1/(s r_c) - 1)/\Delta\tau)` (the
       slope from the face's node heights, computed in the test; 1e-13, single 1e-5), the limit binds on some rows
       and not on others, :math:`1/(1 + \text{term}\,\Delta\tau) \ge s r_c` on every sloped row, the u and v
       rows are bit-identical to the call without the margin; with margin 100 (:math:`s r_c > 1` on the steepest
       rows) those rows get no cap and no row gets less than none; a margin of 0.5, NaN or infinity is refused; on
       1, 2 and 4 ranks. :math:`r_c` is pinned at :math:`a = 0.18` (0.04135501, the smaller root of
       :math:`a r^2 - (2a + 4) r + a`). ``Predictor.IndependentOfTheBoxLayout`` also runs the cap with margin 4 on one
       box and on 2 × 2 boxes (w, :math:`a_P` and :math:`d` of the w rows agree to 1e-12).
   * - ``Predictor.WVnnCapSlopeMarginIsBitIdenticalOnFlatGround``
     - On flat ground (several boxes) the cap :math:`V = 20` with margin 4 and with margin 0 gives
       the same w rows (:math:`a_P`, :math:`d`) and the same u, v and w, bit for bit; the cap itself changes
       :math:`a_P`. :math:`1/20` is inexact, so a regrouped term would differ.
   * - ``Simplec.OneCorrectionRemovesTheDivergence``
     - One pressure correction of a predicted velocity (random wind, slip ground and top, periodic sides) leaves a
       divergence below the Poisson tolerance of the predicted one (0.1 and 1e-10, 2-norms), the ground and top faces
       of w untouched (:math:`d = 0`), :math:`p` changed from a random start by :math:`\alpha_p` times the potential
       (:math:`\alpha_p` = 1 and 0.7); the stage's iterate reports the same divergences and six residual norms.
   * - ``Simplec.PoiseuilleToMachinePrecision``
     - Poiseuille flow between no-slip walls (uniform mesh, 16 layers, :math:`\nu` = 10 m²/s, F = 1e-4 m/s²): the stage
       at CFL 10 reaches the discrete closed form :math:`u_k = (F/2\nu)(z_k(H - z_k) + \Delta z^2/4)` (the odd mirror of
       the walls) to 1e-13 of :math:`u_{\max}` (single 1e-5) in under 400 iterations; w, v exactly 0, :math:`p'` at
       round-off, θ and k untouched; one more iteration changes u by less than 1e-14 (the fixed-point check, also run on a uniform state in
       ``Simplec.AUniformStateIsAFixedPoint``).
   * - ``Simplec.LidDrivenCavityRe100``
     - The lid-driven cavity at Re 100 (64 × 64 in x–z, two periodic cells in y, no-slip walls, the lid at U = 1,
       :math:`\nu` = 0.01) run by the driver (four sweeps, CFL 1 → 4) to its verdict: u on the vertical centreline and
       w on the horizontal one within 0.01 U of Ghia, Ghia and Shin (1982) at their tabulated points (measured 0.0045
       and 0.0083; 1760 iterations, 14 s); nothing varies in y.
   * - ``Simplec.CorrectionTendsToTheProjectionAsDtauGrows``
     - With :math:`d = 1/\Delta\tau` on every interior face and :math:`\Delta\tau` = 1e3, 1e6, 1e9, the corrected
       velocity equals the plain projection's to 1e-9 of U (single 1e-4), :math:`\Delta\tau\,p'` equals the
       projection's potential to 1e-8, and the divergence falls to the solver tolerance: no breakdown as
       :math:`\Delta\tau \to \infty`.
   * - ``Simplec.KStaysAboveTheFloor``
     - With the closure and a strong random wind, k stays at or above the floor after each of five iterations: a floor
       of 1.5 above most of the initial k (the clip acts and its count is reported) and a floor of 1e-6 below every
       value (the natural case: no cell is clipped), every field finite.
   * - ``Simplec.AUniformStateIsAFixedPoint``
     - A uniform wind with slip walls has a residual at round-off and a full iteration changes every field by less
       than 1e-14 of its scale (the Poisson solve of a zero right-hand side returns round-off); the monitors are the
       speed and the mean k.
   * - ``Simplec.PointImplicitCouplingsArePinned``
     - The point-implicit buoyancy and rotation of the predictor pinned against their formulas on a
       random stratified rotating state, with ``rotation_relax`` 1 and 0.3: :math:`a_P` and :math:`d` of u, v, w,
       the stiffness rows, and the corrections equal to the z-line solves of the modified rows (1e-13; single
       1e-4).
   * - ``Simplec.MovingLidAndLateralWalls``
     - The moving lid and lateral walls: the ghost of u above a moving top is 2U − u, the ghost of
       w beyond an x wall is −w, u is 0 on the x walls (2 n_y n_z Dirichlet rows), and the checks refuse a wall
       velocity on another face type, a normal component and a non-finite value, and a lateral surface layer.
   * - ``Simplec.InputReadersAndRefusals``
     - The readers of ``simplec.*``, ``forcing.*`` and ``convection.*``: values, defaults, range
       checks naming the key (``rotation_relax`` in (0, 1] and ``w_vnn_cap`` 0 or finite and ≥ 16 included), repeated keys, keys with no effect (``ug``
       without Coriolis and without geostrophic forcing, ``theta0`` without gravity, ``smooth_eps`` without the van
       Albada limiter), the value counts; the
       stage's refusals (sweeps < 1,
       no boundary conditions, a closure not the inputs', :math:`\alpha_p` out of range, a negative or non-finite
       SIMPLEC denominator, a non-positive step, a scalar sweep of a velocity, a negative floor).
   * - ``Simplec.WallKRowsAreFolded``
     - The wall-k rows folded into the implicit rows: with a surface layer the first cells' k
       rows are, for ``bc.surface_wall_k = dirichlet``, the diagonal :math:`\max(u_*, u_{*,\min})/d_0` alone with
       every neighbour zero, and for ``source`` the plain rows without the ground neighbour plus that diagonal.
   * - ``Terrain.GroundAndTopCarryNoFlux``
     - Over a 60 m hill (and over flat ground) the ground and top faces of w
       take the no-flux value after the fill: Ω vanishes on every ground and top face to 1e-12 (it is 0.1 or more
       before the fill on a random w), the two ghost layers mirror about the face value, the lateral ghost columns of
       the face rows are the periodic images of the valid faces, the Dirichlet rows of the residual vanish, and the
       homogeneous fill of a perturbation gives the same linear relation; over flat ground the value is 0.
   * - ``Terrain.RayleighFieldsAndHeight``
     - The Rayleigh layer's fields and height: with ``bc.rayleigh_fields = w theta`` the residual
       and the implicit diagonal of u, v and k stay 0 and those of w and θ do not; θ's residual is c(θ − θ_s) with
       the sounding at the absolute height (and at the height above the ground with ``bc.rayleigh_height =
       above_ground``) at every damped cell (1e-12; single 1e-4); the readers of the two keys, the refusal of an
       unknown field name and of a bad height word, and the no-effect check without a layer.
   * - ``Predictor.LusgsSweepMatchesTheDenseReference``, ``Predictor.LusgsOnTwoBoxesConvergesAndIsRankIndependent``
     - One and two symmetric line LU-SGS sweeps on a 6 × 3 × 4 box (cells and x faces; non-periodic, and periodic
       cells where the seam couplings are the previous pass's values) equal a dense block-triangular reference (the
       forward and backward solves with the lower or upper x/y couplings and every z coupling in the matrix, the rest
       on the right-hand side) to 1e-12 of the solution; 40 sweeps reach the dense solution of the full system to
       1e-11.
   * - ``Predictor.LusgsOnePassIsExactAlongTheWind``
     - Pure first-order upwind convection along x on a non-periodic box, wind in +x or −x: one symmetric sweep equals
       the dense solution to 1e-12, a second sweep changes nothing (1e-13), and two z-line Jacobi sweeps are still
       1e-3 off (the lag).
   * - ``Predictor.LusgsConvergesForEitherWindSign``
     - The real :math:`A_1` of u (implicit vertical advection and diffusion, the principal diffusion) on a periodic 4
       × 4 × 6 hill mesh with a uniform wind of either sign at Δτ = 20 s: the LU-SGS error falls with every sweep,
       is below the Jacobi error after four sweeps, and reaches 1e-10 of the dense solution after 40.
   * - ``Simplec.WallKRowsAreFolded`` (the scalar sweep)
     - The scalar sweep of k with ``lusgs`` is the ``zline_lusgs`` solve of the same rows (1e-13), differs from the
       Jacobi solve, and with a floor above the smallest result the clip holds k ≥ floor everywhere (cells clipped
       > 0).
   * - ``Predictor.LusgsFixedPointAndRefusals``
     - A zero right-hand side gives a zero correction and zero change for one and three sweeps; the dispatcher
       routes by kind (``line_sweeps``); sweeps < 1 and a wrong A are refused; the readers of ``steady.sweep``
       (values, the default, an unknown word refused, ``zline_sweeps`` ≥ 1 with ``lusgs`` and ≥ 2 with
       ``zjacobi``).
   * - ``Initializer.ColumnGuessAndProjectionSwitch``
     - The in-memory table profile evaluates its table (exact at the heights, linear between, constant beyond) and
       refuses a first height ≠ 0, non-increasing heights, θ ≤ 0, k < 0 and a single row; with the column override the
       first guess takes u, v and k from the table at the cell-centre heights of a flat mesh (1e-13) and θ from the
       sounding, and the override is required; ``init.projection = 0`` leaves the first guess untouched (no Poisson
       solve, divergence before = after) while the projection changes it over a hill; the readers of ``projection``,
       ``column_tol``, ``column_max_iterations`` (values, defaults, a bad value, no effect without the column, α
       without the projection) and ``check_column_guess`` (no forcing, no stage, no surface layer, a heat flux with the
       ``mol`` column mode, the source wall k, each naming its key).
   * - ``Initializer.ColumnGuessRunsTheColumn``
     - The 1-D column on 16 levels over 1000 m under the flat deck's sounding with the hill decks' Rayleigh layer on w
       and θ: converges within the limit, the ground row at z = 0, heights strictly increasing, u and v bounded, v > 0
       near the ground (the Ekman turning at latitude 90), the geostrophic wind at the top within 0.5 m/s, k at the
       floor in the inversion and largest near the ground, the first cell's k = u*²/C_μ0² (1e-12), u* in 0.2–0.6 m/s,
       the table identical on every rank.
   * - ``ColumnStability.*``
     - The stratified column (``init.column_stability``) on 64 × 1 m layers with a Rayleigh layer on θ above 48 m:
       the neutral default holds θ; ``mol`` (L = 500, −100 m), ``heat_flux`` and ``surface_temperature`` each meet
       MOST's closed forms at the first cell (1/L, q, u*, θ₁ − θ_s; Axell's wall k when unstable); ``mol`` reaches
       L = 418 m, past the flux peak (L = 424 m at θ_s = 298.9 K) and above the smallest L, 413.7 m (both from a
       sweep of ``surface_temperature`` columns); L = 50 m stops the run naming a smallest L of 414 m (2 %), and so
       does a tight ``init.column_rib_max``; the iteration limit names ``init.column_max_iterations``; the target
       wind is met between centres (linear in ln(z + z₀)) and below the first centre (MOST's profile), with a
       pressure gradient (:math:`|F| H = u_*^2`, 1e-6) and with geostrophic forcing; the wall gradients are MOST's at
       d/L; the surface solve against its Businger–Dyer closed form (Rib = 0.15, the 0.19 cap, a tighter
       ``init.column_rib_max``), QNSE at Rib = 1 and the −4 bound (``SurfaceSolveMatchesTheClosedForm``); QNSE's ψ on
       the stable side (``QnseStableColumn``); the input checks, including the Rayleigh layer on θ, the 3-D refusals
       of a target, the flux hand-off of ``mol`` and the Rib bounds.
   * - ``StabilityFunctions.*``, ``SurfaceLayer.QnseFrictionVelocity``, ``SurfaceLayer.StableFunctionsKeysAreRead``,
       ``BoundaryFill.GroundGradientIsTheWallLaw``
     - QNSE's φ and ψ against Python values at ζ = 0.1–100 (1e-14) and with distinct coefficients per function, ψ
       against a Simpson integral of (1 − φ)/x (1e-10); the unstable side and the default equal ERF's Businger–Dyer
       functions; the 3-D u* with QNSE matches a bisection (1e-10) and is 8 % above Businger–Dyer's at ζ = 0.087; the
       ground ghosts take QNSE's φ; the keys are read and checked by name. Twelve one-line mutants are caught.
   * - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``
     - On 1, 2 and 4 boxes (an uneven 7 × 6 split and a 2-D 8 × 1 mesh too), periodic and not, the owned regions cover
       every cell and every x, y and z face exactly once (periodic images identified, an image owning nothing); the
       box-local sweeps (zero ghosts, :math:`n_o` = 0 and 2, one and two sweeps) equal a host replay written from the
       ownership rule (1e-13 of the largest value); on one non-periodic box the zero and lagged ghosts give the same
       result, bit for bit.
   * - ``LineSweeps.OverlapsReachTheSameState``
     - A 16 × 8 mesh in 4 × 2 boxes with a strong x coupling (cells and x
       faces, periodic and not): the outer defect correction with two box-local sweeps per step reaches
       :math:`|r| \le 10^{-12}|b|` for :math:`n_o` = 0, 2 and 4, the three states agree to 1e-10, and the outer
       iterations fall strictly with :math:`n_o` (122, 54, 49 non-periodic). Two boxes of four columns with
       :math:`n_o = 4` reproduce the one-box sweep bit for bit, and differ from it with :math:`n_o = 0`.
   * - ``LineSweeps.ConvergedStateDoesNotDependOnTheLayout``
     - A 16 × 16 mesh in 1, 4 and 16 boxes (cells, x and y faces, periodic and not; zero ghosts with :math:`n_o` = 0
       and 2, and the lagged ghosts): the converged states agree to 1e-10. The box-local sweeps depend on the box
       shape, so the test is of the converged state, not of the iterates: one step is checked to differ between the
       layouts.
   * - ``LineSweeps.DiagonalOrderingEqualsTheSequentialSweep``
     - The diagonal wavefront ordering equals the sequential sweep bit for bit (the correction and the reported
       change): every location, 1, 2 and 4 boxes, periodic and not, one and three sweeps, the lagged ghosts and the
       zero ghosts with :math:`n_o` = 0 and 2. A per-rank launch counter (``line_sweep_launches``) pins which path
       ran: two launches per box and sweep for the sequential order, two per wavefront plane of the swept region for
       the diagonal one, so an ordering that is ignored fails. The stage tests
       (``Predictor.APIsTheDiagonalPlusThePseudoTimeTerm``, ``Simplec.WallKRowsAreFolded``) check the same counts
       through ``momentum_predictor`` and ``scalar_sweep``, so the options are shown to reach the sweeps.
   * - ``LineSweeps.RefusalsAndNames``, ``Driver.DivergenceStopsAndTheChecksNameTheirInput``
     - The refusals naming ``steady.line_overlap`` (< 0, > 8, with the lagged ghosts, wider than the narrowest box),
       no overlap along a one-cell axis, the readers of ``steady.line_ghosts``, ``line_overlap`` and ``line_order``
       (values, defaults, unknown words, a repeated key, a line key with ``zjacobi``), the dispatcher. The reader
       of ``steady.line_alternate`` (the value, the default 0, 2 and -1 refused, refused with ``zjacobi`` and with one sweep).
   * - ``LineSweeps.AlternatingOrderingSweepsBothDiagonals``
     - The alternating ordering (``steady.line_alternate = 1``). On an 8 × 8 × 4 system coupled only to the
       :math:`(i-1)` and :math:`(j+1)` neighbours (the upwind couplings of :math:`u > 0`, :math:`v < 0`), two
       alternating sweeps solve it to 1000 machine epsilons (the reversed forward pass is an exact solve), while two
       sweeps of the ascending pairs leave a residual above 1e-3 of the right-hand side. The diagonal ordering equals
       the sequential one bit for bit with the reversed pairs too, at the launch counts of
       ``LineSweeps.DiagonalOrderingEqualsTheSequentialSweep`` (three layouts, periodic and not, every location,
       lagged ghosts and zero ghosts with :math:`n_o = 2`, one to three sweeps). One sweep is unchanged by the option
       (only the odd sweeps are reversed); two and three sweeps differ. The range check names the key.
   * - ``MeshHierarchy.ZLevelsAreIdenticalOnEveryLevel``
     - Three levels of a stretched mesh (dz0 = 40 m, the ratio solved by Newton), flat and over a round STF hill:
       every level's z-levels equal level 0's bit for bit, with the same nz and z extent and x, y halved per level;
       on flat ground every node of every level sits at its z-level exactly.
   * - ``MeshHierarchy.GroundSmoothingPasses``
     - The pass counts :math:`p(2^l - 1)` pinned for p = 0 to 3 and l = 0 to 4; for p = 0, 1, 2 every coarse ground
       of the hill mesh (non-periodic x, periodic y, 8 boxes) equals an independent host rebuild (injection, then
       that many passes, β = 1/8, in double precision) to 1e-12 of the hill height and differs from one pass more
       or fewer by more than 1 mm; level 0 is the run's own mesh; a constant 37.5 m ground is unchanged on every
       level, bit for bit, ghost nodes included.
   * - ``MeshHierarchy.FreeStreamOnEveryLevel``
     - Free-stream preservation on every level for BTF, STF and Sullivan: the metric divergence of a uniform (10, 5)
       m/s wind below 1e-12 of :math:`(|U| + |V|)/\Delta x` of the level; the coarse grounds are not flat.
   * - ``MeshHierarchy.FlatTerrainColumnsAgree``
     - Flat terrain, three levels (16 × 8, 8 × 4, 4 × 2 × 24 cells, stretched): (a) each level's own Poisson solver
       gives the same φ column for a right-hand side in z alone; (b) the full residual of a conventionally neutral
       column state (the k-L closure, the surface layer, Coriolis with a geostrophic wind, a capping inversion with
       buoyancy) has the same u, v, w, p, θ and k columns on every level, so every level has the same column
       solutions; (c) ten SIMPLEC iterations keep every column of a level alike. Relative 1e-12. Part (b) compares the
       steady equations rather than iterates; the pull request that added the test gives the reason.
   * - ``MeshHierarchy.InvalidLayoutsAbort``
     - Every invalid hierarchy input gives one message naming the key (n_cell x or y not divisible, levels beyond
       the cells, a box size not a multiple of the coarsening, levels 0 and 11, passes −1 and 65, passes without a
       hierarchy), and ``validate`` reports it; ``hierarchy_box_array`` gives boxes of at most
       ``max_grid_size`` that coarsen exactly and span the column; a level 0 laid out by ``make_box_array`` with boxes
       of 6 cells is refused, naming ``mesh.max_grid_size``.
   * - ``MeshHierarchy.GatheredLevelEqualsDistributed``
     - A level gathered onto one box equals the distributed level (8 boxes) bit for bit: node heights with ghost
       nodes, ``detJ``, the three face areas, ``z_cc`` and the quality report, for BTF and STF, on 1, 2 and 4 ranks.
   * - ``Transfer.RestrictionConservesFaceFlux``
     - A 32 × 24 × 10 mesh in 6 boxes of unequal size (non-periodic x, periodic y), three
       levels, flat and over a round STF hill (150 m, half width 300 m), level pairs 0 → 1, 1 → 2 and 0 → 2: every
       coarse face flux equals the sum of the fine ones it covers (:math:`a_x u`, :math:`a_y v`, :math:`a_z \Omega`)
       to 1e-12 of the largest flux; θ, k and p′ equal the ``detJ``-weighted means; the coarse divergence of the
       restricted fluxes (u, v and the returned :math:`\Omega_c`) equals :math:`\sum J_f D_f / (r^2 J_c)` to 1e-12
       (2e-5 in single precision).
       Over the hill the test checks that the coarse areas are not sums of the fine ones and that J varies within
       a coarse cell, so the weights matter. 1, 2 and 4 ranks (6 boxes: uneven on 4).
   * - ``Transfer.ProlongationAndProjectionIsDivergenceFree``
     - Prolongation of a smooth level-1 state onto level 0 (inflow, outflow, surface layer, slip top, periodic y),
       flat and over the hill, with ``gmres_mlmg`` and ``gmres_fft``: the divergence after the projection is below
       ``reltol`` times the one before, in the J-weighted 2-norm for ``gmres_mlmg`` and the plain one for
       ``gmres_fft`` (the norms they stop on). The inflow u keeps its boundary value bit for bit, Ω is zero at the
       ground and the top, p′ and θ are the interpolated ones, and the projection moved the interior faces.
   * - ``Transfer.ProlongationOfBilinearFieldsIsExact``
     - Flat ground, periodic in x or in y, the three level pairs: a field bilinear in x and y (a tent of the period in
       the periodic direction, with its kinks on coarse points of the field's type) is reproduced at every fine point
       of each of the six fields, up to the non-periodic faces and across the seam, to 1e-13 (2e-6 in single
       precision). The coarse ghost cells are set to 1e30, so a read of them would show.
   * - ``Transfer.ProlongedTkeIsAboveTheFloor``
     - A coarse k between −6e-5 and 1.8e-4 against a floor of 1e-4 (the interpolant is negative somewhere and below
       the floor in more than 100 fine cells): after ``prolong_state`` every fine k is at least the floor, the
       reported count equals the number of cells clipped, and the other cells keep their interpolated values bit for
       bit. A floor of −1e-6, NaN or infinity aborts, naming ``turb.tke_floor``.
   * - ``Transfer.GatheredLevelEqualsDistributed``
     - With levels 1 and 2 gathered onto one box, restriction and interpolation between every pair, and
       ``prolong_state`` into level 0, give the distributed hierarchy's results bit for bit, over the hill, on 1, 2
       and 4 ranks.
   * - ``Transfer.RoundTripIsTheClosedFormSmoother``
     - restrict(prolong(c)) on flat ground, periodic in x or in y, against the closed form:
       :math:`(c_{I-1} + 6c_I + c_{I+1})/8` along a cell-centred direction (the identity in the first and last cells
       next to a non-periodic face), the identity along a nodal one, for a non-smooth coarse field; to 2e-14 (2e-6 in
       single precision).
   * - ``Transfer.TargetNodeHeightsComeFromTheHierarchy``
     - Over the hill, for the three level pairs: ``restrict_state`` and ``prolong_fields`` set the target State's
       ``z_phys_nd``, ghost nodes included, to its level's bit for bit, from a target that held −1.
   * - ``Transfer.OneCellDirectionIsConstant``
     - 32 × 4 cells, periodic x, walls in y, three levels: level 2 is one cell wide in y. A tent in x, constant in y,
       prolonged from level 2 to levels 1 and 0, is reproduced at every fine point of each field to 1e-13 (2e-6 in
       single precision). A read of the coarse image's unfilled ghost cells shows as a NaN under ``amrex.init_snan``,
       even with a zero weight.
   * - ``Transfer.WrongLevelPairsAbort``
     - Level pairs (1, 1), (2, 1), (−1, 1) and (0, 3) abort in ``restrict_state`` and ``prolong_fields``, with the
       message naming the rule 0 ≤ fine < coarse < levels.
   * - ``Sequencing.LooseResidualBranchStopsAtThePinnedIteration``, ``Sequencing.LooseResidualBranchWaitsForTheReference``, ``Sequencing.LooseMonitorBranchStopsAtThePinnedIteration``, ``Sequencing.NonFiniteHistoryNeverStops``
     - The loose stop on synthetic histories, pinned by an independent Python implementation: three residuals (one
       peaking at its third value, one below ``residual_tol`` from the start) stop at iteration 92 with two orders and
       137 with three, while the strict verdict is not met; residuals far below their first value do not stop before
       the reference is taken (iteration 4 of 5); two monitors (one with a scale above its value, an absolute test)
       stop at 88 with ``steady.sequencing_monitor_tol = 0.01`` and never with the default 0 (off); constant monitors
       stop when the window is full (50), and not before the reference is taken (a window of 10 with 30 reference
       iterations stops at 29); a history with an infinite residual never stops.
   * - ``Sequencing.SeededReferenceIsTheDropOrigin``, ``Sequencing.InputChecksNameTheirKeys``
     - ``Sequencing.SeededReferenceIsTheDropOrigin``: a residual of 1e-2 is two orders down from a seeded
       reference of 1 (stop at iteration 4) and never from its own reference of 1e-2; a larger value inside the
       reference window raises a seeded reference. ``Sequencing.InputChecksNameTheirKeys``: every range check of
       ``check_sequencing_options`` names its key; the equal-to-strict drop and a monitor tolerance of 0 are
       accepted; with sequencing off the defaults are not range-checked (a strict drop of 1 or 0 is accepted); a coarsest
       level of 4 × 2 is accepted and 2 × 1 refused. Mutants: 12 one-line mutants of
       ``loose_stop``, ``seed_reference`` and the checks; 11 killed, 1 equivalent (the floor of the reference at
       ``residual_tol``, which matters only where the below-tolerance branch already holds).
   * - ``Driver.ReleasedHoldRampsFromTheFirstIteration``, ``Sequencing.ReleasedHoldStopsFromTheFirstIteration``
     - ``Driver.ReleasedHoldRampsFromTheFirstIteration``: residuals 1/4, 1/16, … with a seed of 1 give the CFLs
       1, 1, 1, 4, 16, 64 with the hold and 1, 4, 16, 64, 64 without it (r = 0.25, 0.0625); a seed of 0.1 is
       raised to the first residual (r = 1, then 0.25: 1, 1, 4, 16, 64); ``start_cfl(8)`` gives 8, 4, 16, 64, and
       a start below ``steady.cfl0`` is raised to it; the released hold survives a write, read and restore (the
       continued run is the uninterrupted one), a format-1 history reads as the hold on, a hold other than −1 or 0
       is refused on read and on restore; the calls are refused after a step, without a seed and (``start_cfl``)
       with the hold on or out of range. ``Sequencing.ReleasedHoldStopsFromTheFirstIteration``: a residual of 1e-2
       under a seed of 1 stops at iteration 4 with the hold and at 2 without (the third-iteration minimum; at 0
       with one reference iteration); a falling and a rising residual stop at 3; the monitor branch waits for the
       reference iterations with either hold. ``Sequencing.InputChecksNameTheirKeys`` also checks
       ``steady.sequencing_cfl_start = previous`` against the hold and the has-no-effect messages. All pins were
       re-derived by an independent Python script. Mutants: 21 one-line mutants of the driver, the checkpoint
       format, ``loose_stop`` and the checks, all killed; 4 more of the application wiring are killed by the
       ``rt_sequencing_hold_*``, ``rt_sequencing_cfl_previous`` and ``rt_sequencing_restart`` tests.
   * - ``Simplec.PointImplicitCouplingsArePinned`` (the planar-mean split)
     - The rotation's pin holds the planar-mean split: :math:`a_P` with :math:`f^2\max(\tau_c, \tau_{\mathrm{own}})`
       and the right-hand side with the unique-face horizontal mean of the other residual per level at :math:`\tau_c`
       and the fluctuation at :math:`\tau_{\mathrm{own}}`.
   * - ``SteadyCriterion.ConvergingHistoryIsSteadyAtThePinnedIteration``, ``SteadyCriterion.StallingHistoryIsNeverSteady``, ``SteadyCriterion.ResidualAtTheToleranceIsNotBelowIt``, ``SteadyCriterion.OscillatingMonitorIsNeverSteady``, ``SteadyCriterion.OscillatingResidualIsNeverSteady``
     - Steady-criterion verdicts on converging, stalling (including a residual stuck exactly at the tolerance) and
       oscillating histories. The converging history becomes
       steady at iteration 136, pinned from an independent Python implementation of the definition.
   * - ``SteadyCriterion.NonFiniteValuesMeanDiverged``
     - A NaN or infinite residual or monitor gives ``Diverged``, permanently.
   * - ``SteadyCriterion.NearZeroMonitorUsesItsScale``
     - A monitor at round-off level around zero is steady with scale 1 and never steady with a
       scale below the noise.
   * - ``SteadyCriterion.WindowNeedsAllItsValues``
     - The monitor window needs :math:`W+1` values; a change of exactly the tolerance passes and
       one just above it does not (values exact in binary).
   * - ``SteadyCriterion.DropIsMeasuredFromTheEarlyPeak``
     - The residual drop is measured from the early peak (iteration 38), not from the first
       value (iteration 49).
   * - ``SteadyCriterion.ParameterChecks``
     - Parameter range checks name their keys; bad or missing monitor scales and size
       mismatches stop.
   * - ``SteadyCriterion.SlowMonitorStopsNearItsLimit``
     - A monitor converging at rate 0.999 stops at iteration 4602, within the tolerance of its
       limit; the window test alone would stop at 1615 with a 2 % error.
   * - ``SteadyCriterion.OscillationsUpToFourWindowsAreNeverSteady``
     - Undamped oscillations with periods of 37, 100 and 200 iterations are never steady.
   * - ``SteadyCriterion.ZeroResidualReferenceDoesNotBlock``
     - A residual that is zero over the reference iterations does not block the stop.
   * - ``SteadyCriterion.RestoredHistoryGivesTheSameVerdicts``
     - A criterion restored from a saved history gives the same verdict at every later
       iteration as one never interrupted; a mismatched history is refused.
   * - ``SteadyCriterion.SlowlyGrowingMonitorIsNotSteady``
     - A slowly growing monitor is not steady while its change over the window is still small.
   * - ``Inputs.*``
     - Input parsing; unknown, repeated, multi-valued and out-of-range keys reported with their names.
   * - ``IO.PlotfileRoundTripIsBitIdentical``
     - Plotfile round trip, bit-identical, including the nodal heights and their ghost nodes.
   * - ``IO.CheckpointRestartIsBitIdenticalOnANewLayout``, ``IO.CheckpointRestartRestoresNodalGhostsOnTheSameLayout``, ``IO.CheckpointWithOneGhostCellRestarts``
     - Checkpoint and restart onto another box layout, bit-identical; ghost nodes restored without overwriting
       valid nodes; a checkpoint written with one ghost cell on the solution fields (the earlier
       layout) restarts bit-identical.
   * - ``Reductions.*``
     - The sum gives the same bits on 1, 2 and 4 ranks and within 16 ε Σ\|x\| of the exact sum; the maximum
       handles NaN and −0.
   * - ``Newton.JvpEpsSweepAgainstQuadraticResidual``
     - On a residual that is exactly quadratic in q (upwind convection, Coriolis, slip walls), the analytic :math:`Jv`
       and :math:`B(v,v)` from central differences; the sweep of ``steady.newton_fd_eps_scale`` from 1e-2 to 3.2e-12
       (its lower bound is 1.8e-12) gives the V-shaped error, the truncation branch equal to
       :math:`\epsilon\|B\|/\|Jv\|` to 1e-3, and the default within a factor 3 of the minimum (1.5× measured); a
       uniform :math:`\Delta\tau` adds :math:`x/\Delta\tau` on the right rows; Dirichlet rows are identity rows.
   * - ``Newton.ColumnConvergesQuadratically``
     - Plain Newton on the Ekman column; the reductions 0.403, 0.171, 0.0131, 2.09e-4, 8.26e-8 pinned to a factor 2,
       the last two steps above the round-off floor gain more than three orders each, and :math:`r_{n+1} \le r_n^2`
       after the first step. In single precision the first three reductions are pinned and the last two steps sit on
       the forward difference's floor (they depend on FMA contraction); the single-precision CI job compiles only, so
       this branch runs on request.
   * - ``Newton.PoiseuilleInOneStep``, ``Newton.GmresSolvesKnownSystemOnAnyLayout``
     - ``Newton.PoiseuilleInOneStep``: one step gives the discrete Poiseuille profile to 1e-9.
       ``Newton.GmresSolvesKnownSystemOnAnyLayout``: ``amrex::GMRES<NewtonVector, NewtonOperator>`` recovers a known
       solution of a non-symmetric system to 1e-5 on one and four boxes (8 × 6 × 8) and on one and nine boxes
       (9 × 9 × 8); the layouts agree to 1e-6, a quarter of either solve's error (two box layouts add their sums
       in different orders and GMRES carries the rounding through its restarts; 6.5e-8 and 1.5e-8 measured, 1.4e-8
       to 3.1e-7 before jvp's step sums were made independent of the distribution), and
       the four boxes give the same iteration count and the same solution bits with every box on rank 0.
   * - ``Newton.*`` (vectors, options, Krylov)
     - ``Newton.VectorOperationsAndBitIdenticalInnerProduct``: the vector operations; the inner
       product counts shared and periodic faces once and has the same bits on every rank and for two distribution
       mappings. ``Newton.JvpBitIdenticalForAnyDistribution``: jvp's step and :math:`Jv` have the same bits with
       every box on rank 0 as on the distributed layout, for 16 rough perturbations on twelve boxes (on 2 and 4 ranks;
       the identity on one). ``Newton.OptionReadersAndRangeChecks``: the four ``steady.newton_*`` keys are read, and each
       range check names the key and value. ``Newton.KrylovMemoryEstimateAndBudget``: the :math:`2m + 3` vector
       count and bytes of a 512 × 512 × 128 box; the budget abort names both keys.
       ``Newton.PreconditionerHookAndOptionsReachTheSolver``: the hook is set up once per step and applied once
       per GMRES iteration, and its output is used; restart and tolerance reach GMRES.
       ``Newton.ColumnRestartAndToleranceSweep``: GMRES iterations against the restart length (:doc:`results`).
   * - ``Newton.*`` (``test_NewtonCoverage.cpp``)
     - Added after #19. ``Newton.TerrainWallLawOperatorMatchesCentralDifference``: on a 30 m hill
       under the wall law, where one ghost fill is not a fixed point, the operator matches a central difference of
       :math:`q \mapsto R(\mathrm{fill}(q))` to 1e-5 on one and six boxes (2.6e-8 to 4.8e-8 measured), and is off by
       more than 1 with the base point :math:`R(q)`. ``Newton.TerrainWallLawStepsReduceResidual``: four default steps
       there reduce :math:`|R|` at every step. ``Newton.UpdateSkippedAboveRhsNorm``: the rule that skips an update.
       ``Newton.KrylovBudgetIsPerRank``: one box per rank fits a budget the total exceeds.
       ``Newton.InnerProductSameBitsWithAllBoxesOnRankZero``: the inner product with every box on rank 0 has the bits
       of the distributed layout (on 2 and 4 ranks). ``Newton.StepDiagnosticsAndEpsRange``: the step's Krylov residual
       and eps range. ``Newton.NonUniformDtauOnFaces``: a non-uniform :math:`\Delta\tau` gives :math:`x/\Delta\tau_f`
       with the face mean of the two cells.
   * - ``Continuation.LineSearchRejectsOvershoot``
     - The 30 m hill with :math:`u` three times larger, plain Newton throughout, GMRES(30) stopped at 200 iterations.
       Plain ``NewtonKrylov::step`` raises :math:`|R|` more than twofold with a Krylov relative residual at most 1,
       and drives :math:`k` negative. Six continuation steps are accepted, each keeping :math:`k` at or above its
       floor; each meets the Armijo condition (written out) with
       :math:`\lambda = \beta^{\text{backtracks}} \ge \lambda_\mathrm{min}`, :math:`|R|` is monotone, and at least one
       step backtracks. From the state before the first backtracked step, with ``steady.newton_lambda_min = 0.9``, the
       full step fails the Armijo test and is rejected, and the retry runs at CFL 2.5 with a finite
       :math:`\Delta\tau`; if it is rejected too, the step fails with :math:`q` unchanged, bit for bit (the test
       prints which; ``SkippedStepCutsCflAndKeepsState`` covers :math:`q` unchanged after a failure everywhere). The
       solves end at their iteration limit, so the trajectory's digits depend on the platform's rounding: the test
       checks these properties rather than the trajectory (pinned to 1e-4 before, it held on Linux and failed on macOS
       arm64 and Windows MSVC).
   * - ``Continuation.DtauRampFollowsSer``
     - Given a sequence of :math:`|R|`, the CFL after each accepted step equals the closed form
       :math:`\min(\mathrm{CFL}_\mathrm{max}, g\,\mathrm{CFL}_n, \mathrm{CFL}_0(|R_0|/|R_n|)^p)` for four SER
       settings, with a flag per setting for whether the growth limit and the cap bind (each asserted); the switch to
       infinite :math:`\Delta\tau` at the cap (and not with ``steady.newton_dtau_inf_at_cap = 0``), the law from the
       cap while infinite, the growth bound at a finite :math:`\Delta\tau` only, and the CFL cut after a failure,
       floored at the smallest normal ``Real``. The defaults give 1.25, 2.5, 5, 10, then infinite.
   * - ``Continuation.*`` (failed steps, line search, inputs)
     - ``Continuation.SkippedStepCutsCflAndKeepsState``: a NaN preconditioner (the update is
       skipped) gives four skipped attempts at CFL 1, 1/4, 1/16, 1/64, a final CFL of 1/256, the status *failed*
       and :math:`q` unchanged; a NaN in :math:`q` is skipped the same way and does not set :math:`|R_0|`.
       ``Continuation.ArmijoOnKnownQuadratic``: on :math:`\phi(\lambda) = r(1 - \lambda + K\lambda^2)` the accepted
       :math:`\lambda` is the largest :math:`\beta^j \le (1 - c)/K` (six cases, one rejected); NaN and infinite
       trials fail without raising the invalid flag. ``Continuation.InputReadersAndRangeChecks``: the defaults, the
       seven new keys read, each range check refusing a value just outside its range (and NaN) and naming the
       key, the reader aborting on those and on repeated or two-valued keys, the SER checks of the driver, and
       the application's unknown-input hint. The line search's growth bound (``steady.newton_ptc_growth``) on the
       same quadratic with :math:`K` = 100: a bound of 100 accepts the full step, ``nextafter(100, 0)`` and 1 do
       not, a NaN full step never is. ``Continuation.GmresIterationLimit``: a limit of 3 stops GMRES at 3
       iterations with status 1.
   * - ``Continuation.*`` (rank counts, restart, reference)
     - ``Continuation.HillSameOnAnyRankCount``: the hill on six uneven boxes reaches
       :math:`10^{-8}|R_0|` in 8 steps at :math:`\lambda = 1`, each CFL from the SER law; the run with every box on
       rank 0 has the same :math:`|R|` history and :math:`q`, bit for bit (on 1, 2 and 4 ranks).
       ``Continuation.RestartIsBitIdentical``: 5 steps (already at infinite :math:`\Delta\tau`), the continuation
       state written and read back with every field compared, and ``run`` for 2 more give the bits of 7
       uninterrupted steps; a missing or corrupt file and an inconsistent state are refused, and ``step``
       refuses a :math:`q` on another BoxArray or DistributionMapping.
       ``Newton.TrueLinearResidualMatchesEstimate``: the true linear residual :math:`|b - Ax|/|b|` of a GMRES
       solve on the hill agrees with GMRES's estimate to 1e-6 (on one and six boxes).
       ``Continuation.ReferenceResidualFromCaller``: a caller's :math:`|R_0|` (ten times the first :math:`|R|`)
       stays, gives a CFL ten times that of the step's own :math:`|R_0|` and is the reference of ``run``'s
       tolerance; a zero, negative, infinite or NaN one is refused. ``Continuation.GrowthBoundAtFiniteDtau``: the
       Ekman column with ``steady.ser_exponent = 0.5`` fails after 7 steps without the growth bound and converges
       in 26 with 1.5, each step accepted by the bound at a finite :math:`\Delta\tau` with :math:`\lambda = 1`.
       ``Continuation.RestoreCapsTheCflAndTheFloorIsChecked``: a checkpointed CFL of 8 restored under
       ``steady.cfl_max = 5`` becomes 5, and 3 stays 3; a negative, infinite or NaN ``turb.tke_floor`` is refused,
       naming the key. ``Continuation.DISABLED_ParameterSweep`` (on request): the sweep of :doc:`results`.
   * - ``NewtonPreconditioner.ApplyIsOneSimplecIteration``
     - On a terrain hill split into six boxes (shared and periodic faces), with non-uniform fields and
       :math:`\Delta\tau`, ``apply`` on :math:`r = -R(q)` with the stage's denominator equals one SIMPLEC stage
       iteration from :math:`q` to 1e-9 (1e-12 measured). Both ``zjacobi`` and ``lusgs`` are checked, and ``lusgs``
       again with ``simplec.w_vnn_cap = 16`` and a slope margin of 1e6 (no cap on any sloped row; a preconditioner
       that drops the margin misses by 0.11); :math:`p` is compared modulo its mean. :math:`\delta q = r` on the
       Dirichlet rows, and every shared-face copy equals its owner's.
   * - ``NewtonPreconditioner.ConstantPressureShiftLeavesTheSolveUnchanged``
     - On the periodic hill (singular pressure problem), a constant added to :math:`r_p` or to :math:`q`'s :math:`p`
       leaves ``apply`` unchanged, and :math:`\delta p` has zero :math:`J`-weighted mean. A Newton step from :math:`q`
       and from :math:`q + c` takes the same GMRES iterations (53) and gives the same increment (to 1e-5; compared in
       double only, as :math:`p + c` keeps about four digits of :math:`p` in single precision).
   * - ``NewtonPreconditioner.IterationsBoundedOverCfl``
     - Four Newton steps on the open hill (9 × 6 × 8, six boxes) at CFL 10, 1e2, 1e3 and 1e4. With ``simplec`` the
       mean GMRES iterations per step grow by less than 50 % (10.25 to 10.5). The mutant check: the stage's
       denominator fails the bound (18 to 360–381, run to run). ``identity`` needs at least five times as many
       iterations at CFL 10 and 1e4.
   * - ``NewtonPreconditioner.*`` (continuation, linearity, inputs)
     - ``NewtonPreconditioner.ContinuationIdentityVersusSimplec``: the continuation loop on the
       open hill with each preconditioner. Both converge to :math:`10^{-8}|R_0|` in 8 steps to the same state;
       ``identity`` keeps the 8 steps of ``Continuation.HillSameOnAnyRankCount``, and ``simplec`` needs 54 GMRES iterations against 976 (pinned:
       8 steps, at most 70 iterations, at most a fifth of the identity's). ``NewtonPreconditioner.ApplyIsLinearInR``:
       on the stratified hill, :math:`\mathrm{apply}(s\,r)/s = \mathrm{apply}(r)` for
       :math:`s = 2, -1, 10^3` to 1e-4 with the default denominators (1e-5 measured; 0.5 in single precision, 0.21
       measured, where step 3's difference of the :math:`\theta` convections about 300 K keeps about one digit); the
       stage's denominators
       miss it (1.03), and so does the full-increment step 3 this test replaced (1.04, the mutant check). The
       preconditioner is set up once per attempt. ``NewtonPreconditioner.ReaderFactoryAndRangeChecks``: the
       default, ``identity``, an unknown value aborting with the key named, the range checks of the class, and
       ``apply`` before ``setup``. CTest ``abort_newton_preconditioner_value``, ``_no_newton`` and ``_no_stage``:
       the application's messages.
   * - ``Smoothness.NewtonCycleReturnsTheSimplecSolution``
     - The 2-D periodic hill with a capping inversion (8 × 2 × 8 cells, two boxes, the k–L closure, the wall law,
       Rayleigh damping of w and θ at the top). From the SIMPLEC state at :math:`r` = 1e-4: ``enter`` (the floor
       frozen), plain Newton on the smooth inputs to :math:`10^{-9}|R_0|` within 14 steps (10 measured), ``leave`` (no
       floor decision moved: no cell reaches the floor on this hill, so the floor freeze itself is tested in
       ``Smoothness.FloorScopeHoldsOnlyTheFloor``; the closure thawed), then the verification pass, SIMPLEC with the
       run's inputs to :math:`r` < 1e-8. Its solution equals SIMPLEC's own (run on from 1e-4 to 1e-8) to 1e-5 of the
       largest :math:`u` and :math:`k` and of 2 K in θ (6e-7 measured). With the standard limiter Newton is still
       above :math:`10^{-8}|R_0|` after the same 14 steps.
   * - ``Smoothness.NewtonCycleReturnsTheSimplecSolution`` (limiters)
     - The same test: the smooth-limiter Newton solution differs from the van Leer SIMPLEC solution by at most 2 % of
       the largest :math:`u` in :math:`u` and :math:`w`, 0.01 K in θ and 5 % of the largest :math:`k` in :math:`k`
       (1.1 %, 0.09 %, 7e-5 K and 3.2 % measured), and by more than 0.1 % in :math:`u` (the limiters differ).
   * - ``Smoothness.*`` (floor scope, gate, inputs)
     - ``Smoothness.FloorScopeHoldsOnlyTheFloor``: frozen with the floor scope, the closure
       evaluates as the live one after the θ gradient is reversed (every cell changes branch; bit for bit, no change
       counted), where the full freeze keeps the stable branch and counts every cell; with :math:`k` below the floor
       (1e-20, 0, −1e-10) the floor scope floors it as the live closure does (finite, bit for bit) and the full
       freeze keeps it unfloored; a cell held at the floor on an unstable state with :math:`k` back at 0.4 evaluates as
       the live closure at the floor (:math:`\nu_t > 0`). ``Smoothness.RiGateClosedFormAndLimits``: the gate's closed form
       at nine (:math:`N^2`, :math:`S^2`) points (sign kept, 0 at :math:`N^2 = 0`), the limits :math:`|Ri| = 10^3
       Ri_0` (kept to 1e-6) and :math:`10^{-3} Ri_0` (suppressed to 1e-6), ``turb.ri_gate = 0`` bit-identical to the
       default, a neutral state bit-identical with the gate on, a stratified one changed with :math:`P_b` of the
       ungated :math:`N^2`; the range check (0, or 1e-6 to 1) names the key. ``Smoothness.ModuleInputsAbortsAndCheckpoint``: the limiter swap (van Leer only; :math:`\epsilon` and
       the other inputs kept), ``turb.frozen_switch`` → the floor freeze, ``enter``/``leave`` (the moved floor
       decisions counted, the closure thawed; nothing frozen without a closure or with the switch off), the aborts
       (twice, out of order, another closure), and the checkpoint read on one box from four (the frozen evaluation
       bit for bit; a deck with ``turb.frozen_switch = 0`` and a corrupt file abort; no file reads as "not in the Newton
       stage"). ``Closure.InputsAreReadAndChecked`` reads ``turb.ri_gate`` and
       the new default of ``turb.frozen_switch``. CTest ``abort_frozen_switch_value``, ``abort_frozen_switch_no_newton``
       and ``abort_ri_gate_range``: the application's messages. ``Smoothness.DISABLED_Sweep`` (on request): the
       measurements of :doc:`results`.
   * - ``Driver.NewtonSwitchOnceAndFallback``
     - Synthetic histories with ``steady.switch_residual`` 1e-4 and a window of
       5: the relative residual crosses 1e-4 at iteration 3 and stays below; the switch happens once, at iteration 7,
       with the monitors flat. With a monitor still moving by 5 % per iteration (above ``steady.switch_monitor_tol``
       = 1 %) it waits until the monitor has been flat over a whole window (iteration 12); with the tolerance 0 it
       switches at 7 again. A stage that refuses (``accepts_switch`` false) is never entered.
   * - ``Driver.NewtonSwitchOnceAndFallback`` (fallback), ``NewtonStage.FallsBackAfterFailedSteps``
     - The same test: a second stage that asks to go back from its second iterate is left at once (iteration 9,
       inside the hysteresis window), its ``leave`` is called on every exit, the reason is logged, and the forward
       switch comes back only after the window. ``NewtonStage.FallsBackAfterFailedSteps``: on the hill, with a
       preconditioner that returns zero every continuation step fails; the stage asks to go back after
       ``steady.newton_fallback_failures`` = 2 failed steps, not after one, and after one with the key at 1; ``leave``
       thaws the closure.
   * - ``NewtonStage.*``
     - ``NewtonStage.SwitchesOnceAndConvergesInNewton``: the 2-D hill with a capping inversion
       (8 × 2 × 8, van Albada, 6 orders): one switch (at r ≤ 1e-4, iteration 2768), at most 10 accepted Newton steps
       (5 measured), steady in the Newton stage, at most 60 % of the SIMPLEC iterations of a SIMPLEC run (2769 of
       5417) and the same solution to 1e-4 of the largest :math:`u` and :math:`k` (7e-8 measured).
       ``NewtonStage.HandsBackForTheFinalPassAfterTheSwap``: with van Leer the stage hands the run back for the final
       SIMPLEC pass once its smooth equations have converged, thaws the closure, refuses any later switch, and the run
       stays in SIMPLEC over two more windows. ``NewtonStage.RestartInsideTheStage``: a checkpoint two steps into the
       stage restarts it; three more steps give the uninterrupted run's records and fields bit for bit with the
       identity preconditioner, and to 1e-3 of the relative residual with ``simplec`` (1e-4 measured; the cause of the
       difference is not found: the preconditioner's setup and apply give the same bits on the same input). That test
       found that the continuation's diffusivity must be the stage's own (the SIMPLEC stage's is not restored by a
       restart inside the Newton stage), and the fallback test that a failed first step must not count as converged.
   * - CTest ``rt_steady_hill2d_newton``, ``abort_newton_keys_no_newton_stage``
     - ``rt_steady_hill2d_newton`` (smoke, two ranks, about 10 s): the 2-D hill at 32 × 1 × 16
       without Coriolis (``inputs_steady_hill2d_newton``: truly neutral, a pressure gradient, van Albada, GMRES(200),
       the switch at 1e-3) switches to Newton once, finishes in the Newton stage without a fallback, and ends steady
       within 1800–2199 iterations (1868 measured on 1, 2 and 4 ranks; SIMPLEC alone takes 2906).
       ``abort_newton_keys_no_newton_stage``: a Newton key and a
       switch key with ``steady.stage = simplec`` abort, naming the key.
