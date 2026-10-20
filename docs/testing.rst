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

The one-line mutants of sub-phase 2.2 (``Source/LineSweeps.cpp`` and ``Source/Driver.cpp``), each run against
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
     - UT-2.2.1, UT-2.2.2
   * - M2
     - Face updated twice (box keeps its high face)
     - UT-2.2.1, UT-2.2.2, UT-2.2.3
   * - M3
     - Periodic seam updated twice (first box keeps the domain's low face)
     - UT-2.2.1, UT-2.2.2, UT-2.2.3
   * - M4
     - Forward wavefront skips its last plane
     - UT-2.2.4
   * - M5
     - Backward wavefront skips plane 0
     - UT-2.2.4
   * - M6
     - Ghost δ not zeroed (neighbours' values exchanged between passes)
     - UT-2.2.1
   * - M7
     - Non-owned face copies keep their own value (owner mask all 1)
     - UT-2.2.1, UT-2.2.2, UT-2.2.3
   * - M8
     - One-cell periodic axis not folded into the diagonal
     - UT-2.2.1
   * - M9
     - ``steady.line_overlap`` not read
     - UT-2.2.5 (``Driver``)
   * - M10
     - ``steady.line_ghosts`` not read
     - UT-2.2.5 (``Driver``)
   * - M11
     - ``steady.line_order`` not read
     - UT-2.2.5 (``Driver``)
   * - M12
     - ``steady.line_overlap`` above 8 accepted
     - UT-2.2.5 (``Driver``)
   * - M13
     - A negative ``steady.line_overlap`` accepted
     - UT-2.2.5 (``Driver`` and ``LineSweeps``)
   * - M14
     - Overlap with lagged ghosts accepted
     - UT-2.2.5 (``Driver`` and ``LineSweeps``)
   * - M15
     - Overlap wider than the narrowest box accepted
     - UT-2.2.5 (``Driver`` and ``LineSweeps``)
   * - M16
     - Line key given with ``zjacobi`` accepted
     - UT-2.2.5 (``Driver``)
   * - M17
     - Odd sweeps not reversed (``steady.line_alternate`` ignored)
     - UT-2.2.6
   * - M18
     - Reversed forward wavefront skips its last plane
     - UT-2.2.6
   * - M19
     - Reversed sequential forward pass walks j ascending
     - UT-2.2.6
   * - M20
     - ``steady.line_alternate`` not read
     - UT-2.2.5 (``Driver``)
   * - M21
     - ``steady.line_alternate = 2`` accepted
     - UT-2.2.5 (``Driver``), UT-2.2.6
   * - M22
     - ``steady.line_alternate`` given with ``zjacobi`` accepted
     - UT-2.2.5 (``Driver``)
   * - M23
     - ``steady.line_alternate = 1`` with ``steady.zline_sweeps = 1`` accepted (the option read and ignored)
     - UT-2.2.5 (``Driver``)

The one-line mutants of sub-phase 2.4 (``Source/Transfer.cpp``), each run against ``Transfer.*`` on two ranks with the
unit-test flags (SNaN fill, FPE trap) by ``tools/transfer_mutants.py``. 23 are killed. Three survive and the script
lists them as equivalent: T6 and T9 change no value, and T24 changes only the ``gmres_fft`` solver's path, so its result
agrees with the original's to the solver tolerance. A first run, without the SNaN and FPE flags and before UT-2.4.7 to
UT-2.4.9 were added, left T10, T14, T18 and T19 alive as well. "Crash" marks a mutant that stops the run with a
signal in that test (SIGILL on macOS arm64), so the tests after it did not run.

.. list-table::
   :header-rows: 1
   :widths: 6 50 44

   * - Mutant
     - Change
     - Failing tests
   * - T1
     - Cell means not ``detJ``-weighted
     - UT-2.4.1
   * - T2
     - x/y-face flux without the fine area factor
     - UT-2.4.1, UT-2.4.6
   * - T3
     - x/y-face flux not divided by r
     - UT-2.4.1, UT-2.4.6
   * - T4
     - z-face flux divided by r instead of :math:`r^2`
     - UT-2.4.1, UT-2.4.6
   * - T5
     - Coarse u not divided by the coarse area factor
     - UT-2.4.1, UT-2.4.6
   * - T6
     - Coarse Ω not divided by the coarse area factor
     - Equivalent: :math:`a_z = 1` on every mesh
   * - T7
     - w restricted in place of Ω
     - UT-2.4.1
   * - T8
     - No zero-gradient rows for the slope term's u
     - UT-2.4.1 (crash)
   * - T9
     - No periodic fill for the slope term's u
     - Equivalent: the slope term reads u on its own cell's x-faces, no horizontal ghost
   * - T10
     - Coarse ``z_phys_nd`` not set by ``restrict_state``
     - UT-2.4.7
   * - T11
     - Cell-centred position shifted by half a fine cell
     - UT-2.4.3, UT-2.4.6, UT-2.4.8
   * - T12
     - No linear extrapolation at a non-periodic face (reads past the last centre)
     - UT-2.4.3 (crash)
   * - T13
     - Periodic direction extrapolates instead of reading the periodic image
     - UT-2.4.6
   * - T14
     - One-cell direction reads a second coarse cell (zero weight, unfilled ghost)
     - UT-2.4.8 (crash)
   * - T15
     - x weights swapped on the upper y row
     - UT-2.4.3, UT-2.4.4, UT-2.4.6, UT-2.4.8
   * - T16
     - w interpolated with the nodal x axis
     - UT-2.4.3 (crash)
   * - T17
     - Coarse image filled without the periodic images
     - UT-2.4.3 (crash)
   * - T18
     - Fine ``z_phys_nd`` not set by ``prolong_fields``
     - UT-2.4.7
   * - T19
     - Fine-below-coarse level check dropped
     - UT-2.4.9
   * - T20
     - k not clipped
     - UT-2.4.4
   * - T21
     - Clipped-cell count not summed over ranks
     - UT-2.4.4
   * - T22
     - Infinite ``turb.tke_floor`` accepted
     - UT-2.4.4
   * - T23
     - No boundary fill before the projection
     - UT-2.4.2
   * - T24
     - Projection corrects the Dirichlet rows too (d = 1 + mask)
     - Equivalent with ``gmres_mlmg`` and ``mlmg``: every Dirichlet row is a domain face with a Neumann condition on p′, so the flux there is 0. With ``gmres_fft`` only to the Poisson tolerance: its preconditioner's cell coefficients next to the face change
   * - T25
     - No projection off the Dirichlet rows (d = −mask)
     - UT-2.4.2, UT-2.4.4, UT-2.4.5
   * - T26
     - p′ not restored after the projection
     - UT-2.4.2

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
   * - ``warn_w_vnn_cap_terrain``, ``warn_w_vnn_cap_terrain_below``, ``abort_w_vnn_cap_floor``
     - *Beyond the plan.* The thin-cell strip on a 50 m ridge
       (``Tests/Infrastructure/Input_Validation/inputs_w_vnn_cap_terrain``, max aspect ratio about 134), one
       iteration: with ``simplec.w_vnn_cap = 16`` and a 0.7 m first cell the run exits 0 and prints the warning that
       names the key, the tested aspect ratio 50.78 and the slow convergence without the cap; with a 4 m first cell
       (aspect ratio below 50.78) it runs without the warning; with 4 it stops on the floor of 16. The flat
       strip of ``rt_thin_cell_w_vnn_cap`` (aspect ratio above 50.78, no slope) prints no warning.
   * - ``rt_askervein_mesh``
     - RT-1.4.1 and RT-1.4.2: ERF's Askervein mesh coarsened to 60 × 60 × 18 (20 m × 1.08 layers, BTF, ERF's terrain
       file): statistics of ``z`` and ``detJ`` and the quality report against ``Tests/Discretization/Terrain_Mesh/Gold/askervein_mesh.json``.
       Double-precision builds only: the statistics use ERF's plotfile reader, which reads double-precision
       plotfiles.
   * - ``rt_steep_ridge_quality``
     - RT-1.4.2: a 600 m ridge under STF: the report against its gold, including both warnings (STF capped, top
       not flat). Double-precision builds only.
   * - ``rt_hill2d_run``, ``rt_hill2d_matches``
     - RT-1.4.3: a run without a solve writes the mapped-grid plotfile of the 2-D hill; ``amrex_fcompare`` against
       ``Tests/Discretization/Terrain_Mesh/Gold/hill2d_plt00000`` at relative tolerance 1e-12 (1e-5 in single precision). The plotfile holds
       cell fields only, so the nodal heights are not compared here; UT-1.4.7 compares them with ERF.

       The golds of the three mesh tests above were written by ablsteady itself (``Tests/Supporting_Files/check_mesh.py --write``):
       they detect changes, they do not show correctness. Correctness of the same code paths comes from UT-1.4.7
       (ERF's own output) and the closed-form unit tests.
   * - ``rt_askervein_hierarchy``, ``rt_askervein_hierarchy_level1_matches``, ``rt_askervein_hierarchy_level2_matches``
     - RT-2.3.1 and RT-2.3.2: the Askervein mesh of ``rt_askervein_mesh`` with STF and a three-level hierarchy
       (``inputs_mesh_askervein_hierarchy``: 60, 30 and 15 cells a side, the same 18 layers, 1 and 3 ground passes).
       In double precision, level 0's ``z`` and ``detJ`` statistics and every level's quality report against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_mesh.json`` (the per-level headers and the summary table included); in both
       precisions, ``z`` and ``detJ`` of levels 1 and 2 against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_level{1,2}_plt00000`` (relative 1e-12; 1e-5 in single precision). The golds
       were written by ablsteady: they detect changes. The coarse grounds were checked against an independent
       rebuild (UT-2.3.2).
   * - ``rt_steep_ridge_hierarchy_quality``
     - RT-2.3.2: the steep ridge of ``rt_steep_ridge_quality`` on 32 × 4 × 16 cells with three levels: the report
       of every level against ``Tests/Discretization/Grid_Hierarchy/Gold/steep_ridge_hierarchy_mesh.json``, with the STF and top warnings (no
       coarse-level hint: levels 1 and 2 cap no more than level 0), and the per-level headers and summary table. Double precision only. The gold was written
       by ablsteady: it detects changes.
   * - ``rt_askervein_hierarchy_hint``
     - RT-2.3.2: the Askervein hierarchy with unsmoothed coarse grounds (``terrain.hierarchy_passes = 0``): levels 1
       and 2 cap STF at 9 and 10 z-levels against level 0's 8, so the report gives the coarse-level hint; against
       ``Tests/Discretization/Grid_Hierarchy/Gold/askervein_hierarchy_p0_mesh.json`` (written by ablsteady). Double precision only.
   * - ``rt_hierarchy_restart_first``, ``rt_hierarchy_restart``, ``rt_hierarchy_restart_level{1,2}_matches``
     - A checkpoint of the Askervein mesh written without a hierarchy (2 ranks), restarted with three levels on 3
       ranks: the rebuilt levels 1 and 2 against the golds of ``rt_askervein_hierarchy``.
   * - ``abort_hierarchy``, ``abort_hierarchy_passes_unused``, ``abort_hierarchy_tool``
     - UT-2.3.5 at the application: n_cell x not divisible by 4 for three levels, a box size not a multiple of 4
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
     - Old against new ``steady.sequencing_hold`` (beyond the plan), on the two-level restart deck of the truly
       neutral hill (only level 0 is seeded): with the hold (``1``) level 0's first three iterations run at
       ``steady.cfl0`` with r = 1; without it (``0``, the default) r < 1 from the first iteration and the CFL leaves
       ``steady.cfl0`` at the second; with ``steady.sequencing_cfl_start = previous`` level 0 starts at level 1's
       last CFL (4, after 150 iterations). Seconds. ``rt_sequencing_restart`` also checks that the restored level-0
       history keeps the released hold.
   * - ``abort_poisson_options``
     - ``poisson_projection`` with ``solver.poisson_reltol`` repeated and out of range: both named.
       (``abort_unknown_key`` checks that the application refuses ``solver.poisson_reltol``.)
   * - ``rt_poisson_askervein``, ``rt_poisson_askervein_matches``
     - RT-1.10.1: ``poisson_projection`` projects a fixed wind on the Askervein mesh (30 × 30 × 18). Iterations at
       most 10 % above ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein.json``, divergence after below 10 × the printed tolerance of
       before, and
       :math:`\phi` against ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_plt00000`` at relative 1e-8. Double precision only (the
       golds are double-precision runs); the gold was written by ablsteady, so it detects changes.
   * - ``rt_poisson_askervein_mlmg``, ``rt_poisson_askervein_mlmg_matches``
     - RT-1.11.1: the same projection with back end B (MLMG, 6 iterations): iterations at most 10 % above
       ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_mlmg.json`` (written by ablsteady's first run) and the divergence check, and
       :math:`\phi` against back end A's gold
       ``Tests/Discretization/Poisson_Askervein/Gold/poisson_askervein_plt00000`` at relative 1e-8, since both solve the same operator.
   * - ``rt_walldist_hill2d``, ``rt_walldist_hill2d_matches``
     - RT-1.11.2: ``wall_distance`` on ERF's Neutral_Hill_2D mesh: ERF's wall-distance checks (``Tests/Discretization/Wall_Distance/check_walldist.py``:
       max relative ≤ 0.15, mean ≤ 0.03, ≤ 0.2 cells below 100 m, all finite), and the field against
       ``Tests/Discretization/Wall_Distance/Gold/walldist_hill2d_plt00000`` at relative 1e-8 (a gold written by ablsteady). Double precision only.
   * - ``rt_init_askervein``, ``rt_init_askervein_matches``
     - RT-1.12.1: the mass-consistent initial field on the Askervein mesh (30 × 30 × 18; ``Tests/Initialization/Mass_Consistent_Askervein/inputs_init_askervein``:
       log law u* = 0.6 m/s, z₀ = 0.03 m along x, inflow at xlo, outflow at xhi, :math:`\alpha_h = 1`,
       :math:`\alpha_v = 0.5`; divergence 4.05e-1 → 6.1e-12 1/s). ``Tests/Supporting_Files/check_init.py askervein``: iterations at most 10 % above the gold's 13,
       divergence after at most 10 × ``solver.poisson_reltol`` times before, divergence before (the printed 10 digits),
       max :math:`|w|` and the largest horizontal speed at the gold's values (1e-6), :math:`\theta` = 300 K everywhere, :math:`p' = 0`,
       :math:`k > 0` (``Tests/Initialization/Mass_Consistent_Askervein/Gold/init_askervein.json``); the state against ``Tests/Initialization/Mass_Consistent_Askervein/Gold/init_askervein_plt00000``
       at relative 1e-8. Two ranks; double precision only. :math:`\lambda` is not in the plotfile (a variation from
       the plan): the state plotfile is the regression artefact, and the unit tests pin :math:`\lambda`.
   * - ``rt_init_bump_trend``
     - RT-1.12.2: the hilltop speed-up over FastWindTerrain's bump (``Tests/Initialization/Bump_Speedup_Trend/inputs_init_bump``,
       ``Tests/Initialization/Bump_Speedup_Trend/terrain_bump_fwt.txt``: its ``regtests/phase6_solve_correction`` hill and power-law inflow) for
       :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (five runs on two ranks), against FastWindTerrain's own at commit
       648f8f1 (``Tests/Initialization/Bump_Speedup_Trend/Gold/init_bump_trend.json``, written from its plotfiles by ``Tests/Initialization/Bump_Speedup_Trend/fwt_bump_speedup.py``): the
       same ranking of the five values at 40, 80 and 150 m above the ground (at 20 m two values tie within 0.002 in
       both codes), every speed-up above 1 at 20, 40, 80 and 150 m, the two codes within 0.02, 0.01, 0.005 and 0.003
       of each other at those heights (measured 0.011, 0.0057, 0.0020, 0.0008), and ablsteady's own values and
       iterations pinned (1e-6, 10 %). Double precision only (:doc:`initial_conditions`).
   * - ``rt_bc_neutral_flat``, ``rt_bc_neutral_flat_matches``
     - RT-1.8.1: ``boundary_fluxes`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, 1000 m) with a fixed log-law state
       (u* = 0.35 m/s, z₀ = 0.1 m, 15°), the surface layer on the ground and a slip top: u* recovered at every ground
       cell to 1e-10, the ground stress :math:`u_*^2(\cos, \sin)` to 1e-10, :math:`k_w = u_*^2/C_{\mu 0}^2`; the
       residual of every cell equals its closed form (``Tests/Discretization/Boundary_Fluxes/check_bc.py``: arithmetic face ν times the log law's
       differences, the wall flux below, no stress through the slip top) to 1e-6 of :math:`u_*^2/\Delta z` for u and
       v (1e-13 measured) and 1e-9 K/s for θ (the sounding's diffusion with the top gradient), with the first cell's
       truncation −0.0971 and the top cell's lost flux +1.000 of :math:`u_*^2\cos/\Delta z` pinned; R_k and R_w zero;
       every field finite. The plotfile against ``Tests/Discretization/Boundary_Fluxes/Gold/bc_neutral_flat_plt00000`` at relative 1e-10 (a gold
       written by ablsteady). Double precision only.
   * - ``rt_closure_column``, ``rt_closure_column_matches``
     - RT-1.9.1: ``closure_column`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, neutral) with the wall law's
       equilibrium state (u* = 0.35 m/s, z₀ = 0.1 m, k = u*²/C_μ0²), the uncapped geometric length and the
       logarithmic face mean (``Tests/Discretization/Closure_Column/check_closure.py``): ν_t = κ u* (z + z₀) and L = κ (z + z₀) at every cell
       (1e-9), u* recovered (1e-10), R_u and R_v zero below the top cell (1e-10 of u*²/Δz; the top cell loses
       u*²/Δz through the slip top, 1e-8), P_s, ε and R_k = −(P_s − ε) equal to their closed forms from the discrete
       strain (1e-9 of ε; the second cell's P_s/ε = 1.45 is the centred difference of ln over 2Δz), R_θ and R_w zero.
       The plotfile against ``Tests/Discretization/Closure_Column/Gold/closure_column_plt00000`` at relative 1e-10 (a gold written by ablsteady).
       Double precision only.
   * - ``rt_steady_column``
     - RT-1.15.1: the application on the 1-D neutral column (``Tests/Canonical_RANS/Neutral_Column/inputs_steady_column``: 4 × 4 × 64,
       1000 m, surface layer z₀ = 0.1 m, slip top, uniform θ, the force F = u*²/H with u* = 0.35 m/s, the log law with
       u* = 0.25 m/s as the first guess, CFL cap 50) to the steady verdict (6 orders, 1e-8).
       ``Tests/Canonical_RANS/Neutral_Column/check_steady_column.py``: the verdict; the residual of u down 6 orders from its reference (the
       ``Driver`` history of the checkpoint); mean u* = √(FH) to 1e-4 (the force balance); the first cell's speed the
       wall law's (1e-6); the second cell within 1 % of the log law and the lowest eight within 2 % (the stress of the
       force-driven column is u*²(1 − z/H) and the discrete strain overshoots, RT-1.9.1); the wall k = u*²/C_μ0²
       (1e-6); ν_t of the second cell within 6 % of κ u* (z + z₀); every column identical (1e-10); every field finite.
       Measured: 2231 iterations, u* = 0.3499999, second cell 0.9 % low, cells 2–7 1.3–1.9 % low. Double precision
       only (six orders and 1e-8 are below single's floor); the restart chain below runs in both.
   * - ``rt_steady_column_rough``, ``rt_steady_column_rough_floor_one``
     - *Beyond the plan.* The column on rough ground (``inputs_steady_column_rough``: 4 × 4 × 64, 1 m layers, 64 m,
       z₀ = 0.5 m, F = u*²/H with u* = √0.125 = 0.354 m/s). The first centre, 0.5 m, is below (e − 1) z₀, where
       ln((d + z₀)/z₀) = ln 2 < 1. With the default ``bc.surface_denominator_floor = log_law``: the checks of
       ``rt_steady_column``, with the cells above the first measured against the profile of the shallow column's
       stress u*²(1 − z/H) (``--stress-profile``, a closed form). With ``one``: the first cell's speed is
       (u*/κ) max(ln((z + z₀)/z₀), 1) to 1e-6, 1/ln 2 = 1.44 times the log law's. Measured: 347 and 411 iterations;
       log_law second cell 0.70 % low, lowest eight within 0.92 %. Sub-phase 1.7's code fails the first
       (first cell 0.862 m/s against 0.598). Double precision only.
   * - ``rt_steady_column_restart_first``, ``…_continued``, ``…_straight``, ``…_matches``
     - *Beyond the plan.* 40 iterations of the column (``inputs_steady_column_restart``, the first guess on the
       command line), continued to 80 from the checkpoint (the driver and closure histories restored, the inflow keys
       of the deck accepted on the restart) against 80 straight on the same two ranks: the plotfiles identical under
       ``amrex_fcompare`` at relative 1e-14.
   * - ``rt_steady_neutral_flat``
     - RT-1.15.2: the application on ERF's Neutral_ABL_Flat (``Tests/Canonical_RANS/Neutral_ABL_Flat/inputs_steady_neutral_flat``: 8 × 8 × 64,
       geostrophic wind 10 m/s at latitude 90°, f = 1e-4 1/s, the sounding's inversion from 700 m, slip top with
       0.03 K/m, z₀ = 0.1 m, ERF's closure defaults; four sweeps, CFL cap 2, the log law with u* = 0.35 m/s as the
       first guess) to the steady verdict. ``Tests/Canonical_RANS/Neutral_ABL_Flat/check_steady_neutral.py``: the verdict, then ERF's own
       ``check_neutral.py --physics`` (``extern_erf/``) on the turbulence plotfile: every field finite, k ≥ 0, the
       capped κ(z + z₀) length (1e-3), the AL01 stability functions of R_t (1e-6), K_mv = c_μ √k L (5 %), ε =
       C_μ0³ k^1.5/L (5 %), the wall cell's k balance (1 %), u* within 0.35 ± 0.15 m/s, the wall k within 5 % of
       u*²/C_μ0², the log law in the lowest four cells within 10 %, K_mv of the second cell within 30 % of κ u* z.
       Measured: the verdict after about 2850 iterations (13 s on two ranks), u* = 0.377 m/s, every ERF check
       passing, on one to three ranks with both compilers (the bursts of an earlier version of the predictor's
       rotation step are explained and removed in :doc:`solver`). Label ``regression`` (the run time); double only.
       The 1 % comparison with an ERF reference run is not in this version.
   * - ``rt_thin_cell_w_vnn_cap``
     - *Beyond the plan.* The w rows' VNN cap on a flat inflow/outflow strip with a 0.7 m first cell
       (``Tests/Canonical_RANS/Thin_Cell_Flat/inputs_thin_cell_flat``: 16 × 2 × 40, dx/dz₀ = 134, line LU-SGS with
       two sweeps at CFL 16, the verdict at 6 orders and 5e-7). ``check_w_vnn_cap.py`` runs it with
       ``simplec.w_vnn_cap`` 0 and 16 to the verdict: both reach it, the banner names the cap, the capped run needs
       at most a fifth of the iterations, and u, v, w of the two steady states agree to 2e-4 m/s; then a capped run of
       40 iterations continued to 80 from its checkpoint equals 80 straight in every field, bit for bit. Measured:
       6487 and 687 iterations (9.4 times fewer), u and w within 5.6e-5 and 3.8e-5 m/s. That difference is the
       convergence error of the uncapped run at its verdict: the capped run at 687 iterations is within 5e-11 m/s of
       a 10000-iteration capped run. About 50 s on two ranks; label ``regression``; double only.
   * - ``rt_steady_hill2d``
     - RT-1.16.1: ERF's Neutral_Hill_2D (``Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d``: 128 × 1 × 64, the Witch of Agnesi of
       100 m and 500 m, the flat deck's sounding and forcing, Rayleigh damping of w and θ from 700 m, 16 sweeps, CFL
       cap 4) to the steady verdict with a two-order drop of the u residual; ``Tests/Supporting_Files/check_steady_hill.py`` then runs
       ERF's ``check_hill.py --physics`` on the turbulence plotfile: the wall distance against the exact distance to
       the ridge, finite fields, k ≥ 0, the length-scale bound, the wall cell's k balance, the crest speed-up in the
       three lowest cells within 0.2–0.8 (Jackson and Hunt's 0.4), the speed-up positive in the lowest ten, the
       upstream u* within 0.25–0.55 m/s, the upstream log law in the lowest three cells within 15 %. Measured: the
       verdict after 4590 iterations (49 s on two ranks), speed-up 0.59/0.50/0.43, u* = 0.354 m/s, every check
       passing (:doc:`solver`). Label ``regression``; double only.
   * - ``rt_steady_hill2d_lusgs``
     - RT-2.1.1: the 2-D hill with the line LU-SGS sweeps (``inputs_steady_hill2d_lusgs``: two symmetric sweeps,
       CFL cap 16) to the verdict with a two-order drop and ERF's ``check_hill.py --physics``. Measured: 3062
       iterations (23 s on two ranks), every check passing, the steady state within the convergence tolerance of
       the Jacobi deck's (:doc:`results`). Label ``regression``; double only; since 2.1b it runs in the column
       test's directory (``run/rt_steady_hill2d_column``) so the spin-up check can read both histories.
   * - ``rt_steady_hill2d_column``, ``rt_column_guess_spinup``
     - RT-2.1b.1: the 2-D hill LU-SGS deck started from the converged column (``inputs_steady_hill2d_column``) to
       the verdict with a two-order drop and ERF's ``check_hill.py --physics``; then ``Tests/Canonical_RANS/Neutral_Hill_2D/check_column_guess.py``
       on the Driver histories of the profile run (``rt_steady_hill2d_lusgs``, now in the same run directory) and
       the column run: at iteration 300 the column run's u and v residuals at most half the profile run's, and
       fewer iterations to the verdict. Measured: 1237 iterations (60 under ctest, machine shared s on two ranks) against 3062, residual
       ratios at iteration 300 of 0.084 (u) and 0.060 (v), every ERF check passing. Label ``regression``;
       double only.
   * - ``rt_steady_hill3d``
     - RT-1.16.2: ERF's Neutral_Hill_3D (``inputs_steady_hill3d``: 64 × 64 × 20, the hill of revolution, damping
       from 550 m, 16 sweeps, cap 4) to the verdict with ``check_hill3d.py --physics`` (the 3-D checks, speed-up
       against 1.6 h/L). Measured: 2092 iterations (148 s on four ranks), every check passing. Label ``nightly``
       (minutes); double only.
   * - ``rt_steady_hill3d_neutral_zjacobi``, ``rt_steady_hill3d_neutral_lusgs``, ``rt_steady_hill3d_neutral_matches``
     - RT-2.2.1: a truly neutral Hill_3D (``inputs_steady_hill3d_neutral``: 32 × 32 × 16 over the same domain, θ =
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
       RT-2.2.1 run again with ``gmres_fft`` (``run/rt_steady_hill3d_neutral/lusgs_fft``) to round-off;
       ``check_line_match.py`` (the case folder) compares it with ``rt_steady_hill3d_neutral_lusgs`` (now ``gmres_mlmg``, the log
       pattern ``pressure gmres_mlmg`` checked) at 1e-8 of each field's maximum. Measured: both stop at 1789 iterations; u, v, w, p and k agree to 2.2e-13 absolute or better, θ − mean to 2.4e-11 K. Label
       ``nightly``; double only.
   * - ``rt_steady_hill3d_neutral_sequencing``, ``rt_steady_hill3d_neutral_sequencing_matches``
     - RT-2.5.1 (the same-solution half): the RT-2.2.1 deck with ``mesh.hierarchy_levels = 3`` and
       ``steady.sequencing = 1`` (8 × 8 and 16 × 16 to the loose stop, then 32 × 32 to round-off), its log checked
       for the loose stops of levels 2 and 1, the strict verdict and the summary table; ``check_line_match.py``
       compares its plotfile with ``rt_steady_hill3d_neutral_zjacobi`` at 1e-8 of each field's maximum (the
       per-level plotfiles ``…_plt_level<l>_<step>`` are not counted). Measured: see :doc:`solver` (grid
       sequencing). Label ``nightly``; double only.
   * - ``rt_steady_hill3d_sequencing``
     - RT-2.5.1 (the ERF half): ``inputs_steady_hill3d`` with three levels (16 × 16, 32 × 32, 64 × 64) and
       sequencing, with the checks of ``rt_steady_hill3d``: the steady verdict, a two-order drop of the u residual
       from the criterion's reference (the coarsest level's; ``check_steady_hill.py`` reads it from the checkpoint,
       which for a run without sequencing is the largest of the first 100 residuals, as before) and ERF's
       ``check_hill3d.py --physics``. Label ``nightly``; double only.
   * - ``rt_steady_askervein``
     - RT-1.16.3, the convergence half: Askervein on the ERF grid of RT-1.12.1 (``inputs_steady_askervein``: 30 × 30
       × 18, the log-law inflow, outflow, slip sides, no Coriolis, four sweeps, cap 2) to the steady verdict.
       Measured: 156 iterations, unchanged by the ``gmres_mlmg`` default. ``init.alpha_v = 0.5``, so the log
       must show the projection on ``gmres_fft`` and the stage on ``gmres_mlmg``. The 5 % comparison with ERF's URANS
       on lines A and AA is not in this version (Kestrel). Labels ``regression;smoke``; double only.
   * - ``rt_precursor_askervein_cnbl``, ``rt_precursor_askervein_cnbl_fixed``, ``rt_precursor_askervein_cnbl_table``
     - Beyond the plan: the precursor column of the conventionally neutral Askervein case
       (``inputs_precursor_askervein_cnbl``: 4 × 4 × 26 on Askervein's 20 m × 1.08 levels, z0 = 0.03 m, the flat
       case's sounding, Coriolis at 57.2 N, geostrophic wind (17.028798, −7.380863) m/s) to the steady verdict.
       ``rt_precursor_askervein_cnbl_fixed`` runs the same deck to the step in the table's header (11087), with an
       unreachable drop criterion (the deck is written at configure time; a regenerated table carries its step). ``tools/make_inflow_table.py --check`` on
       that run: the plotfile step equal to the table's, the 10 m wind within 0.01° of +x, u* within 1e-3 of
       0.6 m/s, and ``Tests/Real_Terrain/Askervein_CNBL/askervein_cnbl_inflow.txt`` within 1e-3 of the planar means (u, v relative to the
       geostrophic speed, θ in K, k relative to its maximum). **Why a fixed step.** Under the capping inversion θ
       has no fixed point: the inversion base cools by 1.35e-4 K per iteration as it rises, and the mixed layer
       warms by 1.5e-5 K. The verdict comes when the k residual leaves a limit cycle of about 600 iterations (the
       Poisson iteration count switching between 5 and 6), at an iteration that depends on the compiler's
       rounding: 11087 with Apple clang and 11073 with gcc-16, 1.7e-3 K apart in θ. At step 11087 the two builds
       agree to 2e-4 K. A build that left the cycle one cycle later would fail the table check, and the verdict
       run would show its step. Label ``regression``; double only.
   * - ``rt_steady_askervein_cnbl``
     - Beyond the plan: Askervein in a conventionally neutral boundary layer (``inputs_steady_askervein_cnbl``: 30 ×
       30 × 26 to 1599 m, the precursor's table as inflow on xlo and yhi, outflow on xhi and ylo, Rayleigh damping
       of w and θ above 1200 m, four sweeps, cap 2) to the steady verdict. Measured: 220 iterations (156 for
       ``rt_steady_askervein``); u residual down four orders. Label ``regression``; double only.
   * - ``drift_check_selftest``
     - The ``extern_erf`` drift check on a throw-away ERF repository (:doc:`extraction`).
   * - ``refavg_snapshots``, ``refavg_build``, ``refavg_verify``, ``refavg_matches``
     - ``ablsteady_refavg`` end to end: three copies of the smoke plotfile at 0, 5 and 10 s are archived on 2
       ranks with two inputs files; the archive verifies on 3 ranks with the expected metadata, AU and a zero
       half-window change; the average equals the smoke plotfile at zero tolerance. This checks the round trip;
       identical snapshots cannot check the weights, which UT-1.1.2 does.
   * - ``refavg_damage``, ``abort_refavg_damaged``, ``abort_refavg_inputs``, ``abort_refavg_verify_keys``
     - A copy of the archive with one bit flipped is refused, naming the file; bad tool inputs stop the tool with
       each input named; a build key given with ``refavg.verify`` is reported.
   * - ``versions_check``, ``versions_check_selftest``
     - ``external/VERSIONS`` is well formed and its ``amrex`` entry matches the submodule; the self-test shows a
       mismatch with a throw-away submodule is reported.
   * - ``kestrel_parse_*``, ``kestrel_dryrun``
     - The Kestrel scripts parse; a dry run with stub ``module``/``sbatch``/``srun`` checks the Askervein job
       script (run length, plot period, no 2-D plotfiles, ranks) and archives a fake finished run, averaging only
       the 3-D plotfiles in the final window and archiving the inputs but not the run's outputs; a second case
       with the default prefix ``plt`` keeps its ``plt2d_`` plotfiles out.

GoogleTest suite
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - ID
     - Check
   * - UT-1.1.1
     - A reference archive written, verified and loaded: fields bit-identical to the average, metadata and field
       statistics unchanged; a changed byte, a missing file, an extra file, edited AU, and data changed together
       with its checksum are each found and named, and loading refuses.
   * - UT-1.1.2
     - *Beyond the plan.* The average equals a host-side recomputation bit for bit and, for fields linear in
       time, the mid-window value to 8 ε; trapezoid weights for uneven spacing; the half-window change is exact
       for a linear drift and zero for a steady field; a snapshot with a NaN, fewer than 3 snapshots and times out of
       order are refused with their messages; a trailing ``/`` is handled.
   * - UT-1.1.3
     - *Beyond the plan.* AU = nodes × hours × rate (10 CPU, 100 GPU); metadata range checks name their inputs.
   * - UT-1.1.4
     - *Beyond the plan.* SHA-256 against the NIST example vectors, at every padding boundary and chunking;
       ``hex()`` can be called repeatedly.
   * - UT-1.4.1
     - Vertical levels: uniform; solved ratio (top exact, ratio 1.0421553943273616 from an independent exact-rational bisection);
       fixed ratio (ERF's Askervein 20 m × 1.08, top 749.0048747962335 m); z-levels file and its errors.
   * - UT-1.4.2
     - Flat terrain gives the Cartesian levels: bit for bit for STF and Sullivan, within 2 ulp for BTF.
   * - UT-1.4.3
     - Free-stream preservation on round Witch-of-Agnesi and Gaussian hills, all three mappings: divergence below
       1e-12 of :math:`(|U|+|V|)/\Delta x`. The cancellation is algebraic for any node heights, so this test does
       not detect wrong heights; UT-1.4.7 does.
   * - UT-1.4.4
     - Cell volumes sum to the column volumes above the ground to 1e-12, and to a 1 m quadrature of the domain
       above the hill to 1e-4 (the bilinear ground differs from the continuous hill by 2.4e-5 of the domain).
   * - UT-1.4.5
     - Smoothing: a constant is unchanged; a Fourier mode is damped by the analytic factor per pass; β ≥ 1/4 is
       refused. *Beyond the plan:* on a non-periodic domain, edge nodes use themselves as the outside neighbour and
       ghost nodes take the nearest domain node.
   * - UT-1.4.6
     - IDW returns input points and, with 4 neighbours, is exact on a plane at a symmetric query (with the
       default 6 the two extra neighbours break the symmetry and IDW is not exact on a plane); ERF's bilinear form
       is exact on planes anywhere, with ERF's edge fallback; ``terrain.idw_neighbors`` above the number of points is
       refused.
   * - UT-1.4.7
     - ``z_phys``, ``detJ`` and nodal heights equal ERF's own output for the same 2-D hill deck (BTF, STF,
       Sullivan; STF and Sullivan bit for bit on macOS, BTF to rounding because ablsteady keeps the ground layer as
       given, all within 1e-13 elsewhere); STF capped at 2 of 16 levels, as ERF reports.
   * - UT-1.4.8
     - *Beyond the plan.* The IDW bucket search picks the brute-force neighbours; results agree to 1e-13 (random
       cloud with far queries, a regular grid with exact distance ties, a tie across a bucket edge, a transect
       queried from a 2-D domain, one point); UTM-sized coordinates are read and interpolated in double; file terrain on a periodic domain has one
       seam height; every mapping keeps the ground layer as given; ground at the top and inverted cells stop the
       run; the slope of a diagonal ridge is reported in full; a restart with changed terrain is refused; the report
       of a steep ridge carries both warnings.
   * - UT-1.5.1
     - Convection on a manufactured flow over a 120 m cosine hill (BTF, geometric layers, 16 × 12 × 8 to 64 × 48 × 32
       cells), ghost nodes continuing the layers: the unlimited MUSCL reconstruction (van Albada slope with
       :math:`\epsilon` above every difference) is second order for q (conservative and bounded forms) and for the u, v
       and w control volumes, L2 order ≥ 1.9 on both pairs and max-norm order ≥ 1.9 on the finest. With ``van_leer``
       and ``van_albada``: L1 order ≥ 1.5 on the finest pair and an L1 error at most a quarter of first-order upwind's
       (a limiter is first order at smooth extrema). First-order upwind converges at order 1.0 ± 0.1.
   * - UT-1.5.2
     - Scalar diffusion with a variable :math:`\gamma`, and momentum diffusion with the full stress tensor and a
       variable :math:`\nu`, on the same meshes: L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on the finest,
       for q, u, v and w.
   * - UT-1.5.3
     - Random :math:`k \ge 0` (a third of the points zero) and random transport fluxes of ±50 m/s: the step
       :math:`k - C(k)/a_P` with the bounded upwind form leaves every :math:`k \ge 0`; the conservative upwind form
       goes below −0.1 in the same step. *Beyond the plan:* the bounded form of a constant is zero for every limiter
       and location.
   * - UT-1.5.4
     - Random transport (±20 m/s) and diffusivity (0 to 50 m²/s, a fifth zero) over the hill: every neighbour
       weight of :math:`A_1` is ≥ 0 and :math:`a_P` equals their sum (1e-14; single 1e-6), for cells and the three
       velocity components.
   * - UT-1.5.5
     - The van Albada slope's derivative changes by less than 0.1 across every point where van Leer switches
       (:math:`a = 0`, :math:`b = 0`, :math:`a = -b`, the origin), in three directions; van Leer's jumps by more than 1.
       On the manufactured flow (64 × 48 × 32) the two limiters differ by less than 1 % (L2) for q, u, v and w.
       Invalid ``smooth_eps`` (0, negative, NaN, infinite) is refused.
   * - UT-1.5.6
     - *Beyond the plan.* :math:`A_1` is the explicit operators' first-order part: bounded upwind convection equals
       :math:`A_1 q` for random q and transport (1e-12; single 1e-5), and on a flat mesh :math:`-D(q) = A_1 q` for the
       scalar and :math:`-D_c(u) = A_1 u_c` for each velocity component with random :math:`\nu`. Over the hill, every weight
       equals the explicit face flux of a field that varies along one index direction only (this pins the
       :math:`1 + \delta_{cd}` and :math:`|\mathbf{n}|^2 + n_c^2` factors on sloping faces). MUSCL convection,
       both diffusions and :math:`A_1` agree to round-off (1e-13) on 16² and 8² column boxes, and boxes that share a
       face point agree there (not bit for bit: GCC on arm64 contracts multiply-adds inside vectorised loops only).
   * - UT-1.6.1
     - The divergence (through :math:`\Omega`) and the three components of :math:`G p` over the 120 m hill (BTF,
       geometric layers, 16 × 12 × 8 to 64 × 48 × 32 cells): L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on
       the finest (single precision 1.8). :math:`\langle Gp, U\rangle = -\langle p, DU\rangle` with control-volume
       weights to 1e-12 of either side (7e-12 against 3.3e4), for random p and a random velocity that vanishes near
       the ground and top. *Beyond the plan:* :math:`G` does not read the ghost cells below the ground.
   * - UT-1.6.2
     - The continuity residual is bit for bit ``face_divergence`` of :math:`(u, v, \Omega)`; after a projection
       (Poisson fluxes added to :math:`u, v, \Omega`, w from ``w_from_omega``) the divergence falls from 1.19 to
       2.4e-11, the solve's residual. *Beyond the plan:* ``omega_from_w`` and ``w_from_omega`` are inverses.
   * - UT-1.6.3
     - ``line_jacobian`` for u, v, w, θ and k in x, y and z against a finite-difference Jacobian of the first-order
       operator (transport frozen), per coefficient: relative difference < 1e-6 (1e-13 measured; single 1e-3), over
       the hill with convection, and on a flat mesh with convection and a random diffusivity.
   * - UT-1.6.4
     - Rest over a 300 m hill (slopes up to 0.94), θ = θ_ref(z) stratified, Coriolis on: the momentum, continuity and
       k residuals are below 1e-14. With a constant reference the hydrostatic :math:`p(z)` balances only to
       truncation (R_u 1.6e-3 m/s², 0.5 % of the buoyancy; R_w exact on interior faces).
   * - UT-1.6.5
     - ``jvp`` on a laterally uniform column (affine residual) equals :math:`R(q+v) - R(q)` to 1e-5 (single 1e-2).
       *Beyond the plan:* on the hill with MUSCL convection, against a central difference, L1 error below 1e-5
       (2.7e-6) while a step 1e5 times larger misses (1.8e-4); the result is the same without ``R0``; the step is the
       same on 16² and 8² column boxes and scales with each field (k and its scale times 100 leave it unchanged).
   * - UT-1.6.6
     - *Beyond the plan.* Sources in closed form: Coriolis with the :math:`\cos\varphi` terms and geostrophic forcing
       on a uniform wind (:math:`f = 7.27221\times10^{-5}` 1/s at 30°, T = 86400 s), buoyancy of a uniform excess;
       the geostrophic wind is balanced. The residual equals its operators applied one by one (momentum, continuity,
       θ, k) on the hill with a general state. Invalid options (latitude, period, geostrophic forcing at latitude 0 or with
       the forcing off, θ₀, gravity, profiles that do not increase) and missing inputs stop with a message.
   * - UT-1.4.9
     - *Beyond the plan.* Node heights shared by two boxes (or periodic images) are equal, for BTF, STF and Sullivan,
       periodic and not (under gcc each box's copy differed by up to 1e-13 m before ``Mesh`` synchronised them).
   * - UT-1.4.10
     - *Beyond the plan* (``Mesh.GhostNodesStandForTheirImages``). Every ghost node in x and y, at every k (the rows
       below the ground and above the top included), equals the node it stands for, bit for bit: its periodic image,
       or the edge node of its row at a non-periodic face. BTF, STF and Sullivan; periodic in x only, y only, both and
       neither, and a layout split in z; on 1, 2 and 4 ranks; every copy of a domain node agrees. Before ``fill_ghost_columns`` the corners outside a wall and a periodic side held
       the clamped corner column (up to 8 m off) and the row below the ground of periodic ghost columns was
       extrapolated from the clamped columns (STF, Sullivan).
   * - UT-1.7.1
     - ERF's similarity functions (the copied ``ERF_MOSTUtils.H``, byte-identical by the drift check): Businger–Dyer
       ψ and φ and Jimenez ψ against a Python table of the same formulas, z/L from −5 to 5, to 1e-14 (single 1e-6);
       :math:`d\psi/d\zeta = (1-\phi)/\zeta`.
   * - UT-1.7.2
     - Reference-height lookup on the ERF#3982 column (dz₀ = 10 m, ratio 1.1, 40 layers): centres belong to the interval
       above, faces (10 m included) to the interval that holds them, the last centre to the last interval; heights below
       the first or above the last centre are refused.
   * - UT-1.7.3
     - u* from a log-law state (u* = 0.45, z₀ = 0.1 m) on that column: within 2e-5 at the first centre, at 10 m, at 55.3 m
       and on a centre. *Beyond the plan:* with a heat flux, u* matches an independent bisection to 1e-10 for unstable,
       stable and low-wind convective cases, and is positive. On rough ground (z₀ = 0.5 m, 1 m first layer,
       ln((d + z₀)/z₀) = ln 2): u* from the log law at the first centre and at 0.8 m with
       ``bc.surface_denominator_floor = log_law``, 0.45 ln 2 with ``one``; with a heat flux both floors match the
       bisection with their floor to 1e-10. An unknown floor stops the run, naming the key.
   * - UT-1.7.4
     - :math:`k_w = u_*^2/C_{\mu0}^2` exactly (neutral); Axell's buoyant form with a heat flux; the wall eddy viscosity
       equals :math:`C_{\mu0}\sqrt{k_w}\,\kappa(d+z_0)`.
   * - UT-1.7.5
     - *Beyond the plan.* Over the hill: the columns' distances and areas from their nodes; the wall stress along the
       ground with magnitude :math:`u_*^2|A|`. In the residual, R_u, R_v and R_θ change only in the first cells, by the
       wall flux minus the interior flux over the cell volume; the k rows are the wall-k rows (both modes). The
       ``line_jacobian`` rows lose the ground neighbour and gain the wall terms. At low wind the wall-k row uses
       u* ≥ 0.01 m/s. Invalid options and a reference height outside the columns stop with a message.
   * - UT-1.8.1
     - The Poisson boundary types of every face type on every face it may take, periodic directions included, equal
       ERF's ``get_fft_bc`` (the copied ``ERF_SolverUtils.H``).
   * - UT-1.8.2
     - On a ridge whose height varies along the inflow face: every inflow ghost point (both layers, u, v, θ, k) takes the
       profile at its own height above the ground, the inlet wind varies along the face by more than 0.05 m/s, and the
       outflow face points of u are not Dirichlet rows. *Beyond the plan:* the ghost values of every rule for every
       face type (even, odd, value, profile, gradient, the Dirichlet face points, the top's θ gradient); the analytic
       profiles and their equilibrium k and :math:`\nu_t`, and a file table with and without a k column.
   * - UT-1.8.2b
     - *Beyond the plan.* Over a 120 m hill with a stable heat flux, the ground ghosts of u, v and θ give MOST's
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
   * - UT-1.9.1
     - ``closure_point`` equals ERF's composition of the copied AL01 relations (``ComputeTurbulentViscosityRANS``,
       ρ = 1) bit for bit on 20 000 random (k, N², d): L, ε, R_t, c_μ, c_μ', ν_t, ν_t', across the neutral, stable,
       unstable, smoothed and floored cases (each exercised).
   * - UT-1.9.2
     - Kynema's ``KLAxell`` (transcribed) where the formulations coincide (neutral, z = d, λ = 30 m, no z₀, ρ = 1):
       L, ν_t, ε and P_s to 1e-12.
   * - UT-1.9.3
     - ``smooth_Rt``: identity above R_t^crit, continuous at it, within (R_t^min, R_t^crit] and monotone down to
       R_t = −1e30; the unstable length finite and below its bound for k down to 1e-16 under N² to −1e3; the same
       decision as ERF's.
   * - UT-1.9.4
     - The point-implicit update k + S_k/(1/Δτ + ∂S_k/∂k) with the Patankar diagonal stays ≥ 0 (> 0 with production)
       for 1/Δτ = 1, 1e-3 and 0 on 5000 random states; without the diagonal a decaying state goes negative.
   * - UT-1.9.5
     - The dual-number dS_k/dk and dν_t/dk against central finite differences on 3000 random states (branches held):
       relative 1e-6 (single 1e-2); the dual value part equals the Real evaluation to a few ulp.
   * - UT-1.9.6
     - Every ``turb`` key is read into its field and the constants follow (the floor, C_μ0³); each invalid value is
       refused with its key named, and the constructor refuses it before forming the constants.
   * - UT-1.9.6b
     - *Beyond the plan.* The logarithmic mean against (a − b)/ln(a/b) on random pairs (1e-14), equal values, zero,
       continuity across its series switch; on a log-law column the logarithmic face mean gives the wall flux u*² at
       every interior u face (1e-12) where the arithmetic mean misses the first by 2–12 %; the implicit operator with
       the same mean satisfies A₁u = −D(u). Over a 120 m hill the residual with the closure equals the residual with
       the closure's own fields given minus the sources bit for bit, the k line Jacobian's diagonal differs by dSdk
       exactly, every field is finite with positive diffusivities, misuse (fields given with a closure, no wall
       distance) is refused, and ``height_above_ground`` is the node-mean height above the column's ground.
       Relaxation: ν_t = 0.6 new + 0.4 stored exactly, unrelaxed without history, and the history survives a
       checkpoint round trip onto other boxes (identical relaxed evaluation); frozen switches hold the floor decision
       after k drops below it.
   * - UT-1.10.1
     - Manufactured :math:`\phi` on a stretched STF mesh over a 60 m hill: every lateral BC (periodic, Neumann,
       Dirichlet) with a Dirichlet or Neumann top, 2-D (y periodic) at three resolutions, and a 3-D case (Neumann in x,
       Dirichlet in y, two resolutions): order ≥ 1.9 (single precision: 1.8 coarse, 1.6 fine), residual below 10 ×
       tolerance. The solution has zero gradient at the ground and level side faces, where ERF's operator is exact
       for the mirror ghost cells (next row).
   * - UT-1.10.1b
     - *Beyond the plan.* With a tangential gradient at sloped ground: order 2 on flat ground (≥ 1.95), between 0.8 and
       1.2 over a 120 m hill (64 → 128 → 256 cells), and a fine-mesh error 3 to 6 times that over a 60 m hill (slope²).
       Double precision only. Pins ERF's ground treatment (:doc:`solver`).
   * - UT-1.10.2
     - The anisotropic operator with :math:`\beta_v/\beta_h` = 1, 0.5, 0.1 (``z_cross`` set) on flat ground and over
       the hill, and a spatially varying :math:`\beta` over the hill: order ≥ 1.9 (single: 1.8, 1.6). Without
       ``z_cross`` the hill error is more than 5 times larger.
   * - UT-1.10.3
     - :math:`\beta = 1` gives the no-coefficient path bit for bit (iterations, solution, fluxes); :math:`\beta = 2.5`
       on every face gives :math:`\phi/2.5` and the same fluxes (to 1e-9; single 1e-4). Both paths are this layer's;
       there is no separate plain projection to compare with.
   * - UT-1.10.4
     - On a ramp (slope 0.15) crossing the Neumann faces, their fluxes and the ground's are exactly zero, and
       :math:`\sum J\,A\phi` equals the flux out of the Dirichlet top (1e-12; single 1e-4).
   * - UT-1.10.4b
     - *Beyond the plan.* Operator and fluxes agree cell by cell, :math:`A\phi = -\nabla_h\cdot F` (1e-12 of
       :math:`\max|A\phi|`; single 1e-4), with varying coefficients and ``z_cross`` on the ramp, Dirichlet and Neumann
       top.
   * - UT-1.10.5
     - Projection: the divergence of :math:`u + F` is below 10 × tolerance of the divergence of :math:`u`, and equals
       the residual of the solve.
   * - UT-1.10.6
     - Symmetric hill and right-hand side: :math:`\phi` mirror-symmetric to 1e-12 of :math:`\max|\phi|` (single 1e-5),
       Neumann and Dirichlet sides.
   * - UT-1.10.7
     - *Beyond the plan.* Preconditioner, with bounds 10 % above the counts measured: 19 iterations isotropic (more than
       10 times as many without it); 16 for :math:`\beta_v/\beta_h = 0.1` told the ratio (at most 3 more than
       isotropic; not told, more than 1.5 times as many); 29 for :math:`\beta` varying over 100 times. A singular
       problem reports and removes the mean. Option and BC checks, coefficients that are not positive, MultiFabs on
       another layout and a solve that does not converge each stop with the input named.
       ``tools/poisson_ms_check.py`` (sympy, not run by CI) checks the manufactured solutions' boundary conditions.
   * - UT-1.11.1
     - Back end B (``mlmg`` and ``gmres_mlmg``): every lateral BC × Dirichlet/Neumann top × Neumann (zero-gradient
       solution) / Dirichlet ground, 2-D (y hidden) at three resolutions (two in single precision), and a 3-D case:
       order ≥ 1.9 (single: 1.8, 1.6).
   * - UT-1.11.1b
     - *Beyond the plan.* Back end A with a Dirichlet ground: order ≥ 1.9 for every lateral BC and top.
   * - UT-1.11.2
     - Spatially varying coefficients, Neumann and Dirichlet ground, 3-D with non-periodic sides (Neumann x,
       Dirichlet y and top): both back ends solved to 1e-12 give :math:`\phi` within 1e-9 of max :math:`|\phi|` (single:
       1e-5 and 1e-4), and the fluxes satisfy :math:`\nabla_h\cdot F = -\mathrm{rhs}` to the residual.
   * - UT-1.11.2b
     - *Beyond the plan.* ``z_cross`` with back end B, ``z_cross`` with a Dirichlet ground, ``apply`` with back end B,
       an unknown back end, ``solver.poisson_precond`` with back end B, ``solver.poisson_restart`` with ``mlmg`` and a
       non-converging MLMG solve each stop with the input named.
   * - UT-1.11.3
     - Tucker distance against the exact distance with ERF's tolerances on ERF's decks: flat ground (max relative
       ≤ 3e-3, ≤ 0.1 m below 100 m; single 1e-2 and 0.5 m), 2-D ridge and 3-D hill (max relative ≤ 0.15, mean ≤ 0.03,
       ≤ 0.2 cells below 100 m); the 3-D hill with periodic and Neumann sides; both back ends (back end B on 64-cell
       boxes, under CI's trap). Single precision solves to 1e-4. *Beyond the plan:* on flat ground the error
       converges at second order (32, 64, 128 layers; double precision).
   * - UT-1.11.4
     - *Beyond the plan.* No node of ``z_phys_nd`` (ghosts included) holds ERF's ``bogus_large_value``, for BTF, STF and
       Sullivan with periodic and non-periodic sides; ``fill_ghosts`` fills every ghost cell, corners and between-box
       ghosts included, with odd copies below a Dirichlet ground and even copies at Neumann faces.
   * - UT-1.12.1
     - :math:`\alpha_h = \alpha_v = 1`: the initial field equals the Poisson layer's plain projection of the same first
       guess (the solve without coefficients, applied the same way) bit for bit, :math:`\lambda` included;
       :math:`\alpha_h = \alpha_v = 2` gives the same velocity within 1e-9 of max :math:`|u|` and :math:`\lambda/4`
       (single: 1e-4); the scalar state is untouched by the correction. *Beyond the plan:* scaling the whole z flux by
       :math:`\alpha_v^2` instead of the cross part by :math:`\alpha_h^2` gives a different field; with periodic sides
       the hill case is divergence-free as a singular problem and the seam ghost of :math:`u` is the corrected
       periodic image; a :math:`\lambda` handed in with garbage is zeroed first; the divergence after is measured,
       not assumed.
   * - UT-1.12.2
     - The divergence after the correction is below 10 × ``solver.poisson_reltol`` times the divergence before, in the
       solver's own operator (``velocity_divergence``), in every initializer test.
   * - UT-1.12.3
     - Flat terrain, open and periodic sides: the profile is divergence-free, so :math:`\lambda`, :math:`w` and the
       divergence stay below 1e-12 and :math:`u` equals the log law at every x face's height (single: 1e-4, 1e-5).
   * - UT-1.12.4
     - D4 symmetry: on a square domain with a Gaussian hill, the wind along ±x and ±y (its own inflow and outflow
       faces each time) gives fields that are the eight transforms of each other: the cell-centred velocity vector
       maps with the transform's matrix, :math:`w` and :math:`\lambda` as scalars, within 1e-12 of the maxima (single
       1e-5), with the solve converged to round-off (``reltol`` 3e-15).
   * - UT-1.12.5
     - The speed-up above the summit (four central cells a quarter of the depth up, over the profile at their
       height) falls by more than 0.01 at each step of :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (1.015, 0.998,
       0.975, 0.943, 0.913 on the 1 km test domain, whose open faces are 500 m from the summit), and
       :math:`\alpha_h = 2, \alpha_v = 1` gives the same value as :math:`\alpha_v = 0.5` (the
       ratio alone matters; 1e-9). The preconditioner is told the ratio: at most 21 GMRES iterations at every ratio
       (13 to 19 measured).
   * - UT-1.12.6
     - The first guess over sloping terrain: :math:`u` at every x face, the two ghost layers included, and :math:`v`
       at every y face equal the profile at the face's height above its column's ground with the 3 m floor (the
       floor acts on the ground faces and the ghosts below); :math:`\theta` equals the sounding by absolute height;
       :math:`k` the profile's value with the decay above 200 m, and :math:`u_*^2/C_{\mu 0}^2` everywhere without it;
       :math:`p' = 0`; :math:`w_0 = 0` except on the ground faces, which carry no flux (:math:`\Omega_0 = 0` to
       round-off); every ghost cell of the state is written (no NaN from the allocation); the O'Brien pass refuses a
       layout whose boxes do not span the column.
   * - UT-1.12.7
     - *Beyond the plan.* O'Brien: :math:`\Omega_0` vanishes at the ground and the top exactly, the divergence left in
       every cell is the weight difference :math:`f_{k+1}^2 a_{z,k+1} - f_k^2 a_{z,k}` (:math:`f` the height fraction,
       which the test shows differs from the index fraction on the stretched mesh) times the column's residual
       (1e-11), the column sum equals the residual flux, the residual the function returns equals the largest column's;
       through the initializer the correction starts from a different divergence than without O'Brien (not a lower
       one: the terrain-following first guess is already consistent at the ground) and gives a different :math:`w`.
   * - UT-1.12.8
     - *Beyond the plan.* Inputs: :math:`z_{agl,\min} < 0` or NaN, :math:`\alpha_h \le 0` or infinite,
       :math:`\alpha_v < 0` or NaN, :math:`z_i < 0`, and any key other than ``first_guess`` given without
       ``first_guess = profile`` are refused with the key named, by the checks and by both entry points;
       ``first_guess = none`` leaves the state untouched; the reader takes every key and refuses an unknown
       ``first_guess``, ``obrien = 2``, a repeated key and the keys without the first guess; the Poisson types of the
       correction are Dirichlet on inflow and outflow faces, Neumann elsewhere, periodic where the mesh is; a face set
       that does not fit the mesh is refused.
   * - UT-1.12.9
     - *Beyond the plan* (``Poisson.SplitOperatorKeepsTheMirrorSymmetries``). The split anisotropic operator
       (``z_cross`` set, :math:`\beta_v = 0.36`) on a Gaussian hill varying in x and y keeps the mesh's x and y mirror
       symmetries within 1e-12 of max :math:`|\phi|`; the cross kernel's first transcription broke the y symmetry at
       the 1e-4 level.
   * - UT-1.13.1
     - The local pseudo-time step on known fields: uniform :math:`(u, v, w)` and :math:`\nu` on a flat uniform mesh
       give :math:`\mathrm{CFL} / (|u|/\Delta x + |v|/\Delta y + |w|/\Delta z + 2\nu(1/\Delta x^2 + 1/\Delta y^2))` in
       every cell to 1e-14 (single 1e-6; the hill cases 1e-12 and 1e-5), without :math:`\nu` the diffusive term is gone, a field at rest without
       diffusion gets the floored step, :math:`u` varying in :math:`x` gives the face mean per cell, a
       terrain-following field (:math:`\Omega = 0`) over a 120 m hill has no vertical term and a plain :math:`w` over
       it the term :math:`|w| / (J \Delta\zeta)` cell by cell; a non-positive or NaN CFL and a face layout are
       refused.
   * - UT-1.13.2
     - The SER ramp on a geometric residual (rate 1/4): CFL 1, 1, 1, 4, 16, 64, 64, 64 (``cfl0`` through the two
       reference iterations, then the law capped by the growth 4 and ``cfl_max`` 64), :math:`p = 2` without caps
       256, 4096, 65536, a growth cap of 1.5 binding every iteration, :math:`p = 0` holding the CFL, a residual of
       exactly 0 reaching the cap through the growth cap; the stage's :math:`\Delta\tau` follows the CFL.
   * - UT-1.13.3
     - Switch rules on scripted residuals (:math:`W` = 5, threshold 0.01, back growth 10): a residual alternating
       about the threshold never switches; five iterations below switch forward at iteration 5 with the reason
       naming ``steady.switch_residual``, the residuals and the stage's ``init`` recorded; a 20× growth of the
       second stage's residual switches back at iteration 11 (the window compared must lie within the second stage,
       entered at 6); every pair of switches is at
       least :math:`W` apart; the records carry each iteration's stage; one stage never switches.
   * - UT-1.13.4
     - The history written to ``<dir>/Driver`` and read back is equal field by field (records, switches with their
       reasons, accounts, the criterion's history); a two-stage run interrupted after 6 iterations and continued by a
       new driver from the file on the saved state gives the same records (stage, CFL, relative residual, residuals,
       monitors), switches and verdict iteration as the uninterrupted run; the AU of Kestrel CPU on 2 nodes equal
       20 × wall hours; the history is identical on every rank; a truncated file, a missing file and a history for
       another number of stages are refused.
   * - UT-1.13.5
     - *Beyond the plan.* A residual that becomes infinite or NaN stops the run as ``Diverged`` with no arithmetic on
       it (the traps); ``max_iterations`` bounds
       a run that never settles and a second ``run`` does nothing; no stages, a null stage and stages that disagree
       on residual names are refused (a null first stage too); ``restore`` refuses a NaN CFL, a switch state or
       iteration count that does not fit the records, and brings the CFL within changed bounds; every ``steady.*``
       range check names its key; the reader takes every key,
       its defaults, and refuses a repeated key, an unknown machine and a bad value.
   * - UT-1.14.1
     - :math:`a_P` on every x-, y- and z-face equals :math:`1/\Delta\tau` plus the diagonal of ``line_jacobian``
       exactly, with a random :math:`\Delta\tau` field over a 60 m hill with random diffusivities (both diagonals
       come from ``implicit_rows``, whose coefficients UT-1.6.x pin: the statement is the :math:`1/\Delta\tau`
       addition and the face step); the diagonal exceeds :math:`1/\Delta\tau`; the two copies of the periodic seam
       face have identical rows and identical corrections; the pressure and scalars are untouched and :math:`u`
       changes by exactly the z-line solve's correction.
   * - UT-1.14.2
     - The z-line solve (the vertical upwind weights of :math:`A_1` included) against a dense direct solve of the same
       system,
       assembled independently from the coefficient fields: on a single column (1 × 1 × 12, periodic sides, the
       lagged neighbours are the column itself) to 1e-12 after 60 sweeps (single 1e-4); the ghosts outside the domain
       stay 0.
   * - UT-1.14.3
     - :math:`\alpha_u` relaxation equals pseudo-time relaxation: with :math:`\Delta\tau = \alpha_u/(1-\alpha_u)/
       \mathrm{diag}(A_1)` on every face, :math:`a_P = \mathrm{diag}(A_1)/\alpha_u` to 1e-12 (single 1e-5) and the
       converged correction of u, v and w equals the dense solve of the relaxed system to 1e-12 (single 1e-4).
   * - UT-1.14.4
     - z-line Jacobi on a 3-D linear advection-diffusion problem (4 × 4 × 6, periodic, :math:`\Delta\tau` = 20 s)
       converges to the dense direct solution: the error falls at every doubling of the sweeps (1, 2, 4, 8, 16) and
       is below 1e-10 after 400 sweeps, with the last sweep's change below 1e-10.
   * - UT-1.14.5
     - Fixed point: a uniform wind over flat ground with uniform pressure has a residual of exactly zero and the
       predictor changes nothing (:math:`|\Delta q| = 0`; with a zero right-hand side the solve is zero for any
       coefficients, so this follows from the residual); a perturbation of :math:`u` at one face is pulled back: one
       predictor step with :math:`\Delta\tau` = 10⁴ s and 50 sweeps lowers max :math:`|R_u|` below 0.75 of its value
       (0.56 measured; the pressure and the transport are frozen in the step).
   * - UT-1.14.6
     - *Beyond the plan.* ``face_dtau`` is the mean of the two cells, the one cell at a wall, the ground and the top,
       and the periodic image across the seam; one sweep's change is the whole correction; ``sweeps < 1``, a
       non-positive step and wrong layouts are refused.
   * - UT-1.14.7
     - *Beyond the plan.* With boundary conditions (inflow at xlo, outflow at xhi, slip ground and top), the Dirichlet
       rows (the inflow face of u, the ground and top faces of w) have :math:`a_P = 1`, :math:`d = 0`, and one
       predictor step sets them to the boundary value exactly for any :math:`\Delta\tau` (the inflow faces then have
       a zero residual row); every other row has :math:`a_P = 1/\Delta\tau + \mathrm{diag}(A_1)` and :math:`d = a_P -
       \sum a_{nb} > 0`.
   * - UT-1.14.8
     - *Beyond the plan.* The same state on one box and on 2 × 2 boxes (the lagged neighbours cross the box seams
       through the ghost exchange) gives the same corrected u, v, w, :math:`a_P`, :math:`d` and face step to 1e-12
       (single 1e-5).
   * - ``Predictor.WVnnCapAddsToTheWDiagonalAndD``
     - *Beyond the plan.* The w rows' VNN cap (``simplec.w_vnn_cap``, :doc:`solver`) on the boundary-condition state
       of UT-1.14.7 with :math:`\Delta\tau` of 400–1600 s (an uncapped VNN above the cap): with :math:`V = 16` every
       interior w row has :math:`a_P = 1/\Delta\tau + (a_{ZM} + a_{ZP})/V + \mathrm{diag}(A_1)` and :math:`d` larger by
       :math:`(a_{ZM} + a_{ZP})/V` (1e-13; single 1e-5), the Dirichlet ground and top rows keep :math:`a_P = 1` and
       :math:`d = 0`, the effective VNN stays below :math:`V`, and the w correction is the z-line solve with the
       capped diagonal; u and v, and every field with the cap 0, are bit-identical to the call without the argument;
       a negative, NaN, infinite or subnormal cap is refused. On 1, 2 and 4 ranks (2 × 1 boxes).
   * - UT-1.15.1
     - One pressure correction of a predicted velocity (random wind, slip ground and top, periodic sides) leaves a
       divergence below the Poisson tolerance of the predicted one (0.1 and 1e-10, 2-norms), the ground and top faces
       of w untouched (:math:`d = 0`), :math:`p` changed from a random start by :math:`\alpha_p` times the potential
       (:math:`\alpha_p` = 1 and 0.7); the stage's iterate reports the same divergences and six residual norms.
   * - UT-1.15.2
     - Poiseuille flow between no-slip walls (uniform mesh, 16 layers, :math:`\nu` = 10 m²/s, F = 1e-4 m/s²): the stage
       at CFL 10 reaches the discrete closed form :math:`u_k = (F/2\nu)(z_k(H - z_k) + \Delta z^2/4)` (the odd mirror of
       the walls) to 1e-13 of :math:`u_{\max}` (single 1e-5) in under 400 iterations; w, v exactly 0, :math:`p'` at
       round-off, θ and k untouched; one more iteration changes u by less than 1e-14 (UT-1.15.6 b).
   * - UT-1.15.3
     - The lid-driven cavity at Re 100 (64 × 64 in x–z, two periodic cells in y, no-slip walls, the lid at U = 1,
       :math:`\nu` = 0.01) run by the driver (four sweeps, CFL 1 → 4) to its verdict: u on the vertical centreline and
       w on the horizontal one within 0.01 U of Ghia, Ghia and Shin (1982) at their tabulated points (measured 0.0045
       and 0.0083; 1760 iterations, 14 s); nothing varies in y.
   * - UT-1.15.4
     - With :math:`d = 1/\Delta\tau` on every interior face and :math:`\Delta\tau` = 1e3, 1e6, 1e9, the corrected
       velocity equals the plain projection's to 1e-9 of U (single 1e-4), :math:`\Delta\tau\,p'` equals the
       projection's potential to 1e-8, and the divergence falls to the solver tolerance: no breakdown as
       :math:`\Delta\tau \to \infty`.
   * - UT-1.15.5
     - With the closure and a strong random wind, k stays at or above the floor after each of five iterations: a floor
       of 1.5 above most of the initial k (the clip acts and its count is reported) and a floor of 1e-6 below every
       value (the natural case: no cell is clipped), every field finite.
   * - UT-1.15.6
     - A uniform wind with slip walls has a residual at round-off and a full iteration changes every field by less
       than 1e-14 of its scale (the Poisson solve of a zero right-hand side returns round-off); the monitors are the
       speed and the mean k.
   * - UT-1.15.7
     - *Beyond the plan.* The point-implicit buoyancy and rotation of the predictor pinned against their formulas on a
       random stratified rotating state, with ``rotation_relax`` 1 and 0.3: :math:`a_P` and :math:`d` of u, v, w,
       the stiffness rows, and the corrections equal to the z-line solves of the modified rows (1e-13; single
       1e-4).
   * - UT-1.15.8
     - *Beyond the plan.* The moving lid and lateral walls: the ghost of u above a moving top is 2U − u, the ghost of
       w beyond an x wall is −w, u is 0 on the x walls (2 n_y n_z Dirichlet rows), and the checks refuse a wall
       velocity on another face type, a normal component and a non-finite value, and a lateral surface layer.
   * - UT-1.15.9
     - *Beyond the plan.* The readers of ``simplec.*``, ``forcing.*`` and ``convection.*``: values, defaults, range
       checks naming the key (``rotation_relax`` in (0, 1] and ``w_vnn_cap`` 0 or finite and ≥ 16 included), repeated keys, keys with no effect (``ug``
       without Coriolis and without geostrophic forcing, ``theta0`` without gravity, ``smooth_eps`` without the van
       Albada limiter), the value counts; the
       stage's refusals (sweeps < 1,
       no boundary conditions, a closure not the inputs', :math:`\alpha_p` out of range, a negative or non-finite
       SIMPLEC denominator, a non-positive step, a scalar sweep of a velocity, a negative floor).
   * - UT-1.15.10
     - *Beyond the plan.* The wall-k rows folded into the implicit rows: with a surface layer the first cells' k
       rows are, for ``bc.surface_wall_k = dirichlet``, the diagonal :math:`\max(u_*, u_{*,\min})/d_0` alone with
       every neighbour zero, and for ``source`` the plain rows without the ground neighbour plus that diagonal.
   * - UT-1.16.1
     - *Beyond the plan* (``test_Terrain.cpp``). Over a 60 m hill (and over flat ground) the ground and top faces of w
       take the no-flux value after the fill: Ω vanishes on every ground and top face to 1e-12 (it is 0.1 or more
       before the fill on a random w), the two ghost layers mirror about the face value, the lateral ghost columns of
       the face rows are the periodic images of the valid faces, the Dirichlet rows of the residual vanish, and the
       homogeneous fill of a perturbation gives the same linear relation; over flat ground the value is 0.
   * - UT-1.16.2
     - *Beyond the plan.* The Rayleigh layer's fields and height: with ``bc.rayleigh_fields = w theta`` the residual
       and the implicit diagonal of u, v and k stay 0 and those of w and θ do not; θ's residual is c(θ − θ_s) with
       the sounding at the absolute height (and at the height above the ground with ``bc.rayleigh_height =
       above_ground``) at every damped cell (1e-12; single 1e-4); the readers of the two keys, the refusal of an
       unknown field name and of a bad height word, and the no-effect check without a layer.
   * - UT-2.1.1
     - One and two symmetric line LU-SGS sweeps on a 6 × 3 × 4 box (cells and x faces; non-periodic, and periodic
       cells where the seam couplings are the previous pass's values) equal a dense block-triangular reference (the
       forward and backward solves with the lower or upper x/y couplings and every z coupling in the matrix, the rest
       on the right-hand side) to 1e-12 of the solution; 40 sweeps reach the dense solution of the full system to
       1e-11.
   * - UT-2.1.2
     - Pure first-order upwind convection along x on a non-periodic box, wind in +x or −x: one symmetric sweep equals
       the dense solution to 1e-12, a second sweep changes nothing (1e-13), and two z-line Jacobi sweeps are still
       1e-3 off (the lag).
   * - UT-2.1.3
     - The real :math:`A_1` of u (implicit vertical advection and diffusion, the principal diffusion) on a periodic 4
       × 4 × 6 hill mesh with a uniform wind of either sign at Δτ = 20 s: the LU-SGS error falls with every sweep,
       is below the Jacobi error after four sweeps, and reaches 1e-10 of the dense solution after 40.
   * - UT-2.1.4
     - The scalar sweep of k with ``lusgs`` is the ``zline_lusgs`` solve of the same rows (1e-13), differs from the
       Jacobi solve, and with a floor above the smallest result the clip holds k ≥ floor everywhere (cells clipped
       > 0).
   * - UT-2.1.5
     - A zero right-hand side gives a zero correction and zero change for one and three sweeps; the dispatcher
       routes by kind (``line_sweeps``); sweeps < 1 and a wrong A are refused; the readers of ``steady.sweep``
       (values, the default, an unknown word refused, ``zline_sweeps`` ≥ 1 with ``lusgs`` and ≥ 2 with
       ``zjacobi``).
   * - UT-2.1b.1
     - *Beyond the plan* (``Initializer.ColumnGuessAndProjectionSwitch``). The in-memory table profile evaluates
       its table (exact at the heights, linear between, constant beyond) and refuses a first height ≠ 0,
       non-increasing heights, θ ≤ 0, k < 0 and a single row; with the column override the first guess takes u, v
       and k from the table at the cell-centre heights of a flat mesh (1e-13) and θ from the sounding, and the
       override is required; ``init.projection = 0`` leaves the first guess untouched (no Poisson solve, divergence
       before = after) while the projection changes it over a hill; the readers of ``projection``,
       ``column_tol``, ``column_max_iterations`` (values, defaults, a bad value, no effect without the column, α
       without the projection) and ``check_column_guess`` (no forcing, no stage, no surface layer, a heat flux, the
       source wall k, each naming its key).
   * - UT-2.1b.2
     - *Beyond the plan* (``Initializer.ColumnGuessRunsTheColumn``). The 1-D column on 16 levels over 1000 m under
       the flat deck's sounding with the hill decks' Rayleigh layer on w and θ: converges within the limit, the
       ground row at z = 0, heights strictly increasing, u and v bounded, v > 0 near the ground (the Ekman turning
       at latitude 90), the geostrophic wind at the top within 0.5 m/s, k at the floor in the inversion and largest
       near the ground, the first cell's k = u*²/C_μ0² (1e-12), u* in 0.2–0.6 m/s, the table identical on every
       rank.
   * - UT-2.2.1
     - ``LineSweeps.EveryFaceIsOwnedAndSweptOnce``. On 1, 2 and 4 boxes (an uneven 7 × 6 split and a 2-D 8 × 1 mesh
       too), periodic and not, the owned regions cover every cell and every x, y and z face exactly once (periodic
       images identified, an image owning nothing); the box-local sweeps (zero ghosts, :math:`n_o` = 0 and 2, one
       and two sweeps) equal a host replay written from the ownership rule (1e-13 of the largest value); on one
       non-periodic box the zero and lagged ghosts give the same result, bit for bit.
   * - UT-2.2.2
     - ``LineSweeps.OverlapsReachTheSameState``. A 16 × 8 mesh in 4 × 2 boxes with a strong x coupling (cells and x
       faces, periodic and not): the outer defect correction with two box-local sweeps per step reaches
       :math:`|r| \le 10^{-12}|b|` for :math:`n_o` = 0, 2 and 4, the three states agree to 1e-10, and the outer
       iterations fall strictly with :math:`n_o` (122, 54, 49 non-periodic). Two boxes of four columns with
       :math:`n_o = 4` reproduce the one-box sweep bit for bit, and differ from it with :math:`n_o = 0`.
   * - UT-2.2.3
     - ``LineSweeps.ConvergedStateDoesNotDependOnTheLayout``. A 16 × 16 mesh in 1, 4 and 16 boxes (cells, x and y
       faces, periodic and not; zero ghosts with :math:`n_o` = 0 and 2, and the lagged ghosts): the converged
       states agree to 1e-10. The box-local sweeps depend on the box shape, so the test is of the converged state,
       not of the iterates: one step is checked to differ between the layouts.
   * - UT-2.2.4
     - ``LineSweeps.DiagonalOrderingEqualsTheSequentialSweep``. The diagonal wavefront ordering equals the
       sequential sweep bit for bit (the correction and the reported change): every location, 1, 2 and 4 boxes,
       periodic and not, one and three sweeps, the lagged ghosts and the zero ghosts with :math:`n_o` = 0 and 2.
       A per-rank launch counter (``line_sweep_launches``) pins which path ran: two launches per box and sweep for
       the sequential order, two per wavefront plane of the swept region for the diagonal one, so an ordering
       that is ignored fails. The stage tests (``Predictor.APIsTheDiagonalPlusThePseudoTimeTerm``,
       ``Simplec.WallKRowsAreFolded``) check the same counts through ``momentum_predictor`` and ``scalar_sweep``, so the options
       are shown to reach the sweeps.
   * - UT-2.2.5
     - *Beyond the plan* (``LineSweeps.RefusalsAndNames`` and ``Driver.DivergenceStopsAndTheChecksNameTheirInput``).
       The refusals naming ``steady.line_overlap`` (< 0, > 8, with the lagged ghosts, wider than the narrowest box),
       no overlap along a one-cell axis, the readers of ``steady.line_ghosts``, ``line_overlap`` and ``line_order``
       (values, defaults, unknown words, a repeated key, a line key with ``zjacobi``), the dispatcher. The reader
       of ``steady.line_alternate`` (the value, the default 0, 2 and -1 refused, refused with ``zjacobi`` and with one sweep).
   * - UT-2.2.6
     - *Beyond the plan* (``LineSweeps.AlternatingOrderingSweepsBothDiagonals``). The alternating ordering
       (``steady.line_alternate = 1``). On an 8 × 8 × 4 system coupled only to the :math:`(i-1)` and :math:`(j+1)`
       neighbours (the upwind couplings of :math:`u > 0`, :math:`v < 0`), two alternating sweeps solve it to 1000
       machine epsilons (the reversed forward pass is an exact solve), while two sweeps of the ascending pairs leave
       a residual above 1e-3 of the right-hand side. The diagonal ordering equals the sequential one bit for bit with
       the reversed pairs too, at the launch counts of UT-2.2.4 (three layouts, periodic and not, every location,
       lagged ghosts and zero ghosts with :math:`n_o = 2`, one to three sweeps). One sweep is unchanged by the option
       (only the odd sweeps are reversed); two and three sweeps differ. The range check names the key.
   * - UT-2.3.1
     - Three levels of a stretched mesh (dz0 = 40 m, the ratio solved by Newton), flat and over a round STF hill:
       every level's z-levels equal level 0's bit for bit, with the same nz and z extent and x, y halved per level;
       on flat ground every node of every level sits at its z-level exactly.
   * - UT-2.3.2
     - The pass counts :math:`p(2^l - 1)` pinned for p = 0 to 3 and l = 0 to 4; for p = 0, 1, 2 every coarse ground
       of the hill mesh (non-periodic x, periodic y, 8 boxes) equals an independent host rebuild (injection, then
       that many passes, β = 1/8, in double precision) to 1e-12 of the hill height and differs from one pass more
       or fewer by more than 1 mm; level 0 is the run's own mesh; a constant 37.5 m ground is unchanged on every
       level, bit for bit, ghost nodes included.
   * - UT-2.3.3
     - Free-stream preservation on every level for BTF, STF and Sullivan: the metric divergence of a uniform (10, 5)
       m/s wind below 1e-12 of :math:`(|U| + |V|)/\Delta x` of the level; the coarse grounds are not flat.
   * - UT-2.3.4
     - Flat terrain, three levels (16 × 8, 8 × 4, 4 × 2 × 24 cells, stretched): (a) each level's own Poisson solver
       gives the same φ column for a right-hand side in z alone; (b) the full residual of a conventionally neutral
       column state (the k-L closure, the surface layer, Coriolis with a geostrophic wind, a capping inversion with
       buoyancy) has the same u, v, w, p, θ and k columns on every level, so every level has the same column
       solutions; (c) ten SIMPLEC iterations keep every column of a level alike. Relative 1e-12. See *Variations
       from the plan* in the pull request for why (b) compares the steady equations rather than iterates.
   * - UT-2.3.5
     - Every invalid hierarchy input gives one message naming the key (n_cell x or y not divisible, levels beyond
       the cells, a box size not a multiple of the coarsening, levels 0 and 11, passes −1 and 65, passes without a
       hierarchy), and ``validate`` reports it; ``hierarchy_box_array`` gives boxes of at most
       ``max_grid_size`` that coarsen exactly and span the column; a level 0 laid out by ``make_box_array`` with boxes
       of 6 cells is refused, naming ``mesh.max_grid_size``.
   * - UT-2.3.6
     - A level gathered onto one box equals the distributed level (8 boxes) bit for bit: node heights with ghost
       nodes, ``detJ``, the three face areas, ``z_cc`` and the quality report, for BTF and STF, on 1, 2 and 4 ranks.
   * - UT-2.4.1
     - ``Transfer.RestrictionConservesFaceFlux``. A 32 × 24 × 10 mesh in 6 boxes of unequal size (non-periodic x, periodic y), three
       levels, flat and over a round STF hill (150 m, half width 300 m), level pairs 0 → 1, 1 → 2 and 0 → 2: every
       coarse face flux equals the sum of the fine ones it covers (:math:`a_x u`, :math:`a_y v`, :math:`a_z \Omega`)
       to 1e-12 of the largest flux; θ, k and p′ equal the ``detJ``-weighted means; the coarse divergence of the
       restricted fluxes (u, v and the returned :math:`\Omega_c`) equals :math:`\sum J_f D_f / (r^2 J_c)` to 1e-12
       (2e-5 in single precision).
       Over the hill the test checks that the coarse areas are not sums of the fine ones and that J varies within
       a coarse cell, so the weights matter. 1, 2 and 4 ranks (6 boxes: uneven on 4).
   * - UT-2.4.2
     - ``Transfer.ProlongationAndProjectionIsDivergenceFree``. Prolongation of a smooth level-1 state onto level 0
       (inflow, outflow, surface layer, slip top, periodic y), flat and over the hill, with ``gmres_mlmg`` and
       ``gmres_fft``: the divergence after the projection is below ``reltol`` times the one before, in the
       J-weighted 2-norm for ``gmres_mlmg`` and the plain one for ``gmres_fft`` (the norms they stop on). The inflow
       u keeps its boundary value bit for bit, Ω is zero at the ground and the top, p′ and θ are the interpolated
       ones, and the projection moved the interior faces.
   * - UT-2.4.3
     - ``Transfer.ProlongationOfBilinearFieldsIsExact``. Flat ground, periodic in x or in y, the three level pairs:
       a field bilinear in x and y (a tent of the period in the periodic direction, with its kinks on coarse points
       of the field's type) is reproduced at every fine point of each of the six fields, up to the non-periodic faces
       and across the seam, to 1e-13 (2e-6 in single precision). The coarse ghost cells are set to 1e30, so a read of them would show.
   * - UT-2.4.4
     - ``Transfer.ProlongedTkeIsAboveTheFloor``. A coarse k between −6e-5 and 1.8e-4 against a floor of 1e-4 (the
       interpolant is negative somewhere and below the floor in more than 100 fine cells): after ``prolong_state``
       every fine k is at least the floor, the reported count equals the number of cells clipped, and the other
       cells keep their interpolated values bit for bit. A floor of −1e-6, NaN or infinity aborts, naming
       ``turb.tke_floor``.
   * - UT-2.4.5
     - *Beyond the plan* (``Transfer.GatheredLevelEqualsDistributed``). With levels 1 and 2 gathered onto one box,
       restriction and interpolation between every pair, and ``prolong_state`` into level 0, give the distributed
       hierarchy's results bit for bit, over the hill, on 1, 2 and 4 ranks.
   * - UT-2.4.6
     - *Beyond the plan* (``Transfer.RoundTripIsTheClosedFormSmoother``). restrict(prolong(c)) on flat ground, periodic
       in x or in y, against the closed form: :math:`(c_{I-1} + 6c_I + c_{I+1})/8` along a cell-centred direction
       (the identity in the first and last cells next to a non-periodic face), the identity along a nodal one, for a
       non-smooth coarse field; to 2e-14 (2e-6 in single precision).
   * - UT-2.4.7
     - *Beyond the plan* (``Transfer.TargetNodeHeightsComeFromTheHierarchy``). Over the hill, for the three level
       pairs: ``restrict_state`` and ``prolong_fields`` set the target State's ``z_phys_nd``, ghost nodes included,
       to its level's bit for bit, from a target that held −1.
   * - UT-2.4.8
     - *Beyond the plan* (``Transfer.OneCellDirectionIsConstant``). 32 × 4 cells, periodic x, walls in y, three
       levels: level 2 is one cell wide in y. A tent in x, constant in y, prolonged from level 2 to levels 1 and 0, is
       reproduced at every fine point of each field to 1e-13 (2e-6 in single precision). A read of the coarse
       image's unfilled ghost cells shows as a NaN under ``amrex.init_snan``, even with a zero weight.
   * - UT-2.4.9
     - *Beyond the plan* (``Transfer.WrongLevelPairsAbort``). Level pairs (1, 1), (2, 1), (−1, 1) and (0, 3) abort
       in ``restrict_state`` and ``prolong_fields``, with the message naming the rule 0 ≤ fine < coarse < levels.
   * - UT-2.5.1
     - ``Sequencing.LooseResidualBranchStopsAtThePinnedIteration``, ``…WaitsForTheReference``,
       ``…LooseMonitorBranchStopsAtThePinnedIteration``, ``…NonFiniteHistoryNeverStops``. The loose stop on
       synthetic histories, pinned by an independent Python implementation: three residuals (one peaking at its
       third value, one below ``residual_tol`` from the start) stop at iteration 92 with two orders and 137 with
       three, while the strict verdict is not met; residuals far below their first value do not stop before the
       reference is taken (iteration 4 of 5); two monitors (one with a scale above its value, an absolute test)
       stop at 88 with ``steady.sequencing_monitor_tol = 0.01`` and never with the default 0 (off); constant
       monitors stop when the window is full (50), and not before the reference is taken (a window of 10 with 30
       reference iterations stops at 29); a history with an infinite residual never stops.
   * - UT-2.5.1 (beyond the plan)
     - ``Sequencing.SeededReferenceIsTheDropOrigin``: a residual of 1e-2 is two orders down from a seeded
       reference of 1 (stop at iteration 4) and never from its own reference of 1e-2; a larger value inside the
       reference window raises a seeded reference. ``Sequencing.InputChecksNameTheirKeys``: every range check of
       ``check_sequencing_options`` names its key; the equal-to-strict drop and a monitor tolerance of 0 are
       accepted; with sequencing off the defaults are not range-checked (a strict drop of 1 or 0 is accepted); a coarsest
       level of 4 × 2 is accepted and 2 × 1 refused. Mutants: 12 one-line mutants of
       ``loose_stop``, ``seed_reference`` and the checks; 11 killed, 1 equivalent (the floor of the reference at
       ``residual_tol``, which matters only where the below-tolerance branch already holds).
   * - UT-2.5.1 (beyond the plan, the reference hold)
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
   * - UT-1.15.7 (revised in 1.16)
     - The rotation's pin now holds the planar-mean split: :math:`a_P` with :math:`f^2\max(\tau_c, \tau_{\mathrm{own}})`
       and the right-hand side with the unique-face horizontal mean of the other residual per level at :math:`\tau_c`
       and the fluctuation at :math:`\tau_{\mathrm{own}}`.
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
       valid nodes; a checkpoint written with one ghost cell on the solution fields (before sub-phase 1.5)
       restarts bit-identical.
   * - UT-1.3.4
     - The sum gives the same bits on 1, 2 and 4 ranks and within 16 ε Σ\|x\| of the exact sum; the maximum
       handles NaN and −0.
   * - UT-3.1.1
     - ``Newton.JvpEpsSweepAgainstQuadraticResidual``: on a residual that is exactly quadratic in q (upwind
       convection, Coriolis, slip walls), the analytic :math:`Jv` and :math:`B(v,v)` from central differences; the
       sweep of ``steady.newton_fd_eps_scale`` from 1e-2 to 3.2e-12 (its lower bound is 1.8e-12) gives the V-shaped error, the truncation branch
       equal to :math:`\epsilon\|B\|/\|Jv\|` to 1e-3, and the default within a factor 3 of the minimum (1.5×
       measured); a uniform :math:`\Delta\tau` adds :math:`x/\Delta\tau` on the right rows; Dirichlet rows are
       identity rows.
   * - UT-3.1.2
     - ``Newton.ColumnConvergesQuadratically``: plain Newton on the Ekman column; the reductions 0.403, 0.171,
       0.0131, 2.09e-4, 8.26e-8 pinned to a factor 2, the last two steps above the round-off floor gain more than
       three orders each, and :math:`r_{n+1} \le r_n^2` after the first step. In single precision the first three
       reductions are pinned and the last two steps sit on the forward difference's floor (they depend on FMA
       contraction); the single-precision CI job compiles only, so this branch runs on request.
   * - UT-3.1.3
     - ``Newton.PoiseuilleInOneStep``: one step gives the discrete Poiseuille profile to 1e-9.
       ``Newton.GmresSolvesKnownSystemOnAnyLayout``: ``amrex::GMRES<NewtonVector, NewtonOperator>`` recovers a known
       solution of a non-symmetric system to 1e-5 on one and four boxes (8 × 6 × 8) and on one and nine boxes
       (9 × 9 × 8); the layouts agree to 1e-6, a quarter of either solve's error (two box layouts add their sums
       in different orders and GMRES carries the rounding through its restarts; 6.5e-8 and 1.5e-8 measured, 1.4e-8
       to 3.1e-7 before jvp's step sums were made independent of the distribution), and
       the four boxes give the same iteration count and the same solution bits with every box on rank 0.
   * - UT-3.1.3 (c)
     - *Beyond the plan.* ``Newton.VectorOperationsAndBitIdenticalInnerProduct``: the vector operations; the inner
       product counts shared and periodic faces once and has the same bits on every rank and for two distribution
       mappings. ``Newton.JvpBitIdenticalForAnyDistribution``: jvp's step and :math:`Jv` have the same bits with
       every box on rank 0 as on the distributed layout, for 16 rough perturbations on twelve boxes (on 2 and 4 ranks;
       the identity on one). ``Newton.OptionReadersAndRangeChecks``: the four ``steady.newton_*`` keys are read, and each
       range check names the key and value. ``Newton.KrylovMemoryEstimateAndBudget``: the :math:`2m + 3` vector
       count and bytes of a 512 × 512 × 128 box; the budget abort names both keys.
       ``Newton.PreconditionerHookAndOptionsReachTheSolver``: the hook is set up once per step and applied once
       per GMRES iteration, and its output is used; restart and tolerance reach GMRES.
       ``Newton.ColumnRestartAndToleranceSweep``: GMRES iterations against the restart length (:doc:`results`).
   * - UT-3.1.3 (d)
     - *Beyond the plan* (test_NewtonCoverage.cpp, from the post-merge review of #19).
       ``Newton.TerrainWallLawOperatorMatchesCentralDifference``: on a 30 m hill under the wall law, where one ghost
       fill is not a fixed point, the operator matches a central difference of :math:`q \mapsto R(\mathrm{fill}(q))`
       to 1e-5 on one and six boxes (2.6e-8 to 4.8e-8 measured), and is off by more than 1 with the base point
       :math:`R(q)`. ``Newton.TerrainWallLawStepsReduceResidual``: four default steps there reduce :math:`|R|` at
       every step. ``Newton.UpdateSkippedAboveRhsNorm``: the rule that skips an update. ``Newton.KrylovBudgetIsPerRank``:
       one box per rank fits a budget the total exceeds. ``Newton.InnerProductSameBitsWithAllBoxesOnRankZero``: the
       inner product with every box on rank 0 has the bits of the distributed layout (on 2 and 4 ranks).
       ``Newton.StepDiagnosticsAndEpsRange``: the step's Krylov residual and eps range.
       ``Newton.NonUniformDtauOnFaces``: a non-uniform :math:`\Delta\tau` gives :math:`x/\Delta\tau_f` with the
       face mean of the two cells.
