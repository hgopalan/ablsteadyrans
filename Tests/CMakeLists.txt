# ---------------------------------------------------------------------------
# Tests. Labels follow the development plan:
#   unit        GoogleTest and script self-tests, seconds, every commit
#   regression  whole-application runs
#   smoke       the regression runs CI does on every pull request
# Every test runs in its own directory under <build>/Tests/run/<name>, so
# tests can run in parallel and nothing is written into the source tree.
# ---------------------------------------------------------------------------

# --- GoogleTest (submodule) -------------------------------------------------
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
# MSVC: GoogleTest would link the static runtime (/MT) and AMReX the DLL one
# (/MD); the mismatch fails at link time.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest
                 ${CMAKE_BINARY_DIR}/external/googletest EXCLUDE_FROM_ALL)

add_executable(unit_tests
    Unit_Tests/main.cpp
    Unit_Tests/test_Inputs.cpp
    Unit_Tests/test_IO.cpp
    Unit_Tests/test_Reductions.cpp
    Unit_Tests/test_SteadyCriterion.cpp
    Unit_Tests/test_Sha256.cpp
    Unit_Tests/test_Reference.cpp
    Unit_Tests/test_Mesh.cpp
    Unit_Tests/test_MeshHierarchy.cpp
    Unit_Tests/test_Poisson.cpp
    Unit_Tests/test_PoissonB.cpp
    Unit_Tests/test_WallDistance.cpp
    Unit_Tests/test_Convection.cpp
    Unit_Tests/test_Diffusion.cpp
    Unit_Tests/test_ImplicitOperator.cpp
    Unit_Tests/test_MassFlux.cpp
    Unit_Tests/test_Residual.cpp
    Unit_Tests/test_SurfaceLayer.cpp
    Unit_Tests/test_BoundaryConditions.cpp
    Unit_Tests/test_Closure.cpp
    Unit_Tests/test_Initializer.cpp
    Unit_Tests/test_Driver.cpp
    Unit_Tests/test_Predictor.cpp
    Unit_Tests/test_LineSweeps.cpp
    Unit_Tests/test_Simplec.cpp
    Unit_Tests/test_Newton.cpp
    Unit_Tests/test_NewtonCoverage.cpp
    Unit_Tests/test_Continuation.cpp
    Unit_Tests/test_NewtonPreconditioner.cpp
    Unit_Tests/test_Terrain.cpp
    Unit_Tests/test_Transfer.cpp
    Unit_Tests/test_Sequencing.cpp
    Unit_Tests/test_Smoothness.cpp
    Unit_Tests/test_NewtonStage.cpp
)
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Unit_Tests)
target_link_libraries(unit_tests PRIVATE ablsteady_core GTest::gtest)
# Gold files (Tests/Unit_Tests/Gold) are read from the source tree.
target_compile_definitions(unit_tests PRIVATE ABLSTEADY_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# RT-1.10.1 projects a fixed wind with the Poisson layer: a test-only tool,
# since the application has no Poisson solve yet.
add_executable(poisson_projection Supporting_Files/Tools/poisson_projection.cpp)
target_link_libraries(poisson_projection PRIVATE ablsteady_core)
# RT-1.11.2 computes the Tucker wall distance: another test-only tool.
add_executable(wall_distance Supporting_Files/Tools/wall_distance.cpp)
target_link_libraries(wall_distance PRIVATE ablsteady_core)
add_executable(boundary_fluxes Supporting_Files/Tools/boundary_fluxes.cpp)
target_link_libraries(boundary_fluxes PRIVATE ablsteady_core)
# RT-1.9.1: the full residual of the equilibrium log-law column with the closure.
add_executable(closure_column Supporting_Files/Tools/closure_column.cpp)
target_link_libraries(closure_column PRIVATE ablsteady_core)
if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    setup_target_for_cuda_compilation(unit_tests)
    setup_target_for_cuda_compilation(poisson_projection)
    setup_target_for_cuda_compilation(wall_distance)
    setup_target_for_cuda_compilation(boundary_fluxes)
    setup_target_for_cuda_compilation(closure_column)
endif()

# --- How to launch on N ranks -------------------------------------------------
if(ABLSTEADY_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    separate_arguments(_extra_mpiexec UNIX_COMMAND "${ABLSTEADY_MPIEXEC_PREFLAGS}")
endif()

# Runtime flags of the GPU CI, used on the CPU: no arena preallocation,
# signalling-NaN fills, and a trap on invalid operations.
# A misspelt amrex.* flag would otherwise only print a warning at the end, and
# the check it was meant to switch on would be silently off.
set(ABLSTEADY_TEST_AMREX_FLAGS
    amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1
    amrex.abort_on_unused_inputs=1)

# ablsteady_add_test(NAME <name> NP <n> LABELS <labels...> COMMAND <exe> <args...>
#                    [FIXTURES_SETUP f] [FIXTURES_REQUIRED f] [RUN_DIR dir])
# In a serial build every test runs on one process.
function(ablsteady_add_test)
    cmake_parse_arguments(T "" "NAME;NP;FIXTURES_SETUP;RUN_DIR"
                          "LABELS;COMMAND;PREFIX;FIXTURES_REQUIRED" ${ARGN})
    if(NOT T_NP)
        set(T_NP 1)
    endif()
    if(NOT ABLSTEADY_MPI)
        set(T_NP 1)
    endif()
    if(NOT T_RUN_DIR)
        set(T_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/${T_NAME})
    endif()
    file(MAKE_DIRECTORY ${T_RUN_DIR})

    set(launch "")
    if(ABLSTEADY_MPI)
        set(launch ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${T_NP}
                   ${MPIEXEC_PREFLAGS} ${_extra_mpiexec})
    endif()
    add_test(NAME ${T_NAME} COMMAND ${T_PREFIX} ${launch} ${T_COMMAND}
             WORKING_DIRECTORY ${T_RUN_DIR})
    set_tests_properties(${T_NAME} PROPERTIES LABELS "${T_LABELS}" PROCESSORS ${T_NP})
    if(T_FIXTURES_SETUP)
        set_tests_properties(${T_NAME} PROPERTIES FIXTURES_SETUP ${T_FIXTURES_SETUP})
    endif()
    if(T_FIXTURES_REQUIRED)
        # Quoted: a list of fixtures ("a;b") must stay one argument.
        set_tests_properties(${T_NAME} PROPERTIES FIXTURES_REQUIRED "${T_FIXTURES_REQUIRED}")
    endif()
endfunction()

# One folder per case (see Tests/README.md): its decks, check script and golds.
set(SMOKE_CASE          ${CMAKE_CURRENT_SOURCE_DIR}/Infrastructure/Smoke_Restart)
set(VALIDATION_CASE     ${CMAKE_CURRENT_SOURCE_DIR}/Infrastructure/Input_Validation)
set(REFAVG_CASE         ${CMAKE_CURRENT_SOURCE_DIR}/Infrastructure/Reference_Archive)
set(KESTREL_CASE        ${CMAKE_CURRENT_SOURCE_DIR}/Infrastructure/Kestrel_Dryrun)
set(MESH_CASE           ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Terrain_Mesh)
set(HIERARCHY_CASE      ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Grid_Hierarchy)
set(BC_CASE             ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Boundary_Fluxes)
set(CLOSURE_CASE        ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Closure_Column)
set(POISSON_CASE        ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Poisson_Askervein)
set(WALLDIST_CASE       ${CMAKE_CURRENT_SOURCE_DIR}/Discretization/Wall_Distance)
set(INIT_CASE           ${CMAKE_CURRENT_SOURCE_DIR}/Initialization/Mass_Consistent_Askervein)
set(BUMP_CASE           ${CMAKE_CURRENT_SOURCE_DIR}/Initialization/Bump_Speedup_Trend)
set(COLUMN_CASE         ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Neutral_Column)
set(FLAT_CASE           ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Neutral_ABL_Flat)
set(THIN_CELL_CASE      ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Thin_Cell_Flat)
set(STEEP_RIDGE_CASE    ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Steep_Ridge_Thin_Cell)
set(HILL2D_CASE         ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Neutral_Hill_2D)
set(HILL3D_CASE         ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Neutral_Hill_3D)
set(HILL3D_NEUTRAL_CASE ${CMAKE_CURRENT_SOURCE_DIR}/Canonical_RANS/Truly_Neutral_Hill_3D)
set(ASKERVEIN_CASE      ${CMAKE_CURRENT_SOURCE_DIR}/Real_Terrain/Askervein)
set(CNBL_CASE           ${CMAKE_CURRENT_SOURCE_DIR}/Real_Terrain/Askervein_CNBL)
set(ASKERVEIN_210_CASE  ${CMAKE_CURRENT_SOURCE_DIR}/Real_Terrain/Askervein_210)
set(SUPPORT             ${CMAKE_CURRENT_SOURCE_DIR}/Supporting_Files)

# fcompare is built with the AMReX submodule; with an installed AMReX, look
# for its amrex_fcompare executable.
if(TARGET fcompare)
    set(FCOMPARE $<TARGET_FILE:fcompare>)
else()
    find_program(FCOMPARE amrex_fcompare REQUIRED
                 HINTS ${AMReX_ROOT}/bin ${AMReX_DIR}/../../../bin)
endif()
set(EXPECT_ABORT ${Python3_EXECUTABLE} ${SUPPORT}/expect_abort.py)

# --- Unit tests: UT-1.3.1 to UT-1.3.4 on 1, 2 and 4 ranks ---------------------
if(ABLSTEADY_MPI)
    set(_unit_nps 1 2 4)
else()
    set(_unit_nps 1)
endif()
foreach(np IN LISTS _unit_nps)
    ablsteady_add_test(NAME unit_np${np} NP ${np} LABELS unit
                       COMMAND $<TARGET_FILE:unit_tests> ${ABLSTEADY_TEST_AMREX_FLAGS})
endforeach()

# The same suite with tiling off (one tile per box, as on a GPU).
ablsteady_add_test(NAME unit_notiling NP 2 LABELS unit
                   COMMAND $<TARGET_FILE:unit_tests> ${ABLSTEADY_TEST_AMREX_FLAGS}
                           "fabarray.mfiter_tile_size=1024 1024 1024")

# --- Smoke run, then restart onto another layout and rank count ---------------
set(SMOKE_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/smoke)
# Start from an empty directory: a rerun would otherwise rename the previous
# outputs to *.old.NNNN and pile them up.
add_test(NAME smoke_clean
         COMMAND ${Python3_EXECUTABLE} -c
                 "import shutil, os, sys; shutil.rmtree(sys.argv[1], ignore_errors=True); os.makedirs(sys.argv[1])"
                 ${SMOKE_DIR})
set_tests_properties(smoke_clean PROPERTIES LABELS "regression;smoke" FIXTURES_SETUP smoke_clean)
ablsteady_add_test(NAME smoke_run NP 2 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke_clean FIXTURES_SETUP smoke
                   COMMAND $<TARGET_FILE:ablsteady> ${SMOKE_CASE}/inputs_smoke
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# The smoke state is mostly zeros, so the plotfile comparison below cannot
# tell a restart from a fresh start; the log line can.
ablsteady_add_test(NAME smoke_restart NP 3 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke FIXTURES_SETUP smoke_restart
                   PREFIX ${EXPECT_ABORT} --success
                          "restarted from smoke_chk00000 at step 0" --
                   COMMAND $<TARGET_FILE:ablsteady> ${SMOKE_CASE}/inputs_restart
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# Bitwise: fcompare's default tolerances are zero. The restart used another
# box layout, hence --allow_diff_grids.
ablsteady_add_test(NAME smoke_restart_matches NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke_restart
                   COMMAND ${FCOMPARE} --allow_diff_grids smoke_plt00000 restart_plt00000)

# --- Decks that must abort, with the input named ------------------------------
ablsteady_add_test(NAME abort_unknown_key NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.max_grid_sise.*not an input of namespace 'mesh'"
                          "turb\\.Cmu0 has no effect without steady\\.stage = simplec"
                          "mehs\\.n_cell.*unknown namespace 'mehs'"
                          "steady\\.residual_tol has no effect without steady\\.stage = simplec"
                          "solver\\.poisson_reltol.*not an input of namespace 'solver'" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_unknown_key
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_range NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.n_cell must be three integers >= 1"
                          "mesh\\.is_periodic: z cannot be periodic"
                          "mesh\\.prob_lo and mesh\\.prob_hi must be finite" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_range
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_no_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "forcing\\.gravity has no effect without steady\\.stage = simplec"
                          "convection\\.limiter has no effect without steady\\.stage = simplec"
                          "simplec\\.alpha_p has no effect without steady\\.stage = simplec"
                          "turb\\.Cmu0 has no effect without steady\\.stage = simplec"
                          "steady\\.cfl0 has no effect without steady\\.stage = simplec"
                          "io\\.turb_plot_file has no effect without steady\\.stage = simplec" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_no_stage
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_stage_no_guess NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "steady\\.stage = simplec needs a first guess" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_stage_no_guess
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# Grid sequencing (Sequencing.H): the range checks, a key without
# steady.sequencing = 1, a value other than 0 or 1 (or cfl0 or previous), and
# the keys without the stage.
ablsteady_add_test(NAME abort_sequencing_range NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "steady\\.sequencing_drop_orders = 9\\.000e\\+00: must be > 0 and <= steady\\.residual_drop_orders = 8\\.000e\\+00"
                          "steady\\.sequencing_monitor_tol = 1\\.000e\\+00: must be 0 \\(off\\) or > 0 and < 1"
                          "steady\\.sequencing = 1 needs mesh\\.hierarchy_levels >= 2 \\(it is 1\\)"
                          "steady\\.sequencing_cfl_start = previous needs steady\\.sequencing_hold = 0" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral steady.sequencing=1
                           steady.sequencing_drop_orders=9 steady.sequencing_monitor_tol=1 steady.sequencing_hold=1
                           steady.sequencing_cfl_start=previous ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_coarsest NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "steady\\.sequencing = 1 needs at least 2 cells in x and in y on the coarsest level: mesh\\.n_cell x and y / 2\\^\\(mesh\\.hierarchy_levels - 1\\) = 4 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral steady.sequencing=1
                           mesh.hierarchy_levels=5 "mesh.n_cell=64 16 16" ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_unused NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "steady\\.sequencing_drop_orders has no effect unless steady\\.sequencing = 1"
                          "steady\\.sequencing_hold has no effect unless steady\\.sequencing = 1"
                          "steady\\.sequencing_cfl_start has no effect unless steady\\.sequencing = 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.sequencing_drop_orders=1 steady.sequencing_hold=0 steady.sequencing_cfl_start=cfl0
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# bc.inflow_tke_length takes closure or unbounded only.
ablsteady_add_test(NAME abort_inflow_tke_length_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "bc\\.inflow_tke_length = capped: must be closure or unbounded" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           bc.inflow_tke_length=capped ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.sequencing = 2: must be 0 or 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.sequencing=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_hold_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.sequencing_hold = 2: must be 0 or 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.sequencing=1 mesh.hierarchy_levels=2 steady.sequencing_hold=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_cfl_start_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.sequencing_cfl_start = last: must be cfl0 or previous" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.sequencing=1 mesh.hierarchy_levels=2 steady.sequencing_cfl_start=last ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_sequencing_no_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.sequencing has no effect without steady\\.stage = simplec" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_no_stage steady.sequencing=1
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# steady.newton_preconditioner (sub-phase 3.3): a bad value aborts naming the
# key; a good one aborts without the Newton stage (steady.stage = simplec).
ablsteady_add_test(NAME abort_newton_preconditioner_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.newton_preconditioner = hypre: must be identity or simplec" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.newton_preconditioner=hypre ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_newton_preconditioner_no_newton NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.newton_preconditioner has no effect without steady\\.stage = simplec_newton" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.newton_preconditioner=identity ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_newton_preconditioner_no_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.newton_preconditioner has no effect without steady\\.stage = simplec" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_no_stage steady.newton_preconditioner=simplec
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# turb.frozen_switch (sub-phase 3.4): a bad value aborts naming the key; a good
# one aborts without the Newton stage. turb.ri_gate out of range aborts naming
# the key.
ablsteady_add_test(NAME abort_frozen_switch_value NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "turb\\.frozen_switch must be 0 or 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           turb.frozen_switch=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_frozen_switch_no_newton NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "turb\\.frozen_switch has no effect without steady\\.stage = simplec_newton" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           turb.frozen_switch=1 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_ri_gate_range NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "turb\\.ri_gate = 2\\.0+: must be 0 \\(off\\) or finite, >= 1e-6 and <= 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           turb.ri_gate=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
# Sub-phase 3.5: the Newton stage's keys and the switch rules have no effect
# with steady.stage = simplec (one stage), so they abort there.
ablsteady_add_test(NAME abort_newton_keys_no_newton_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.newton_gmres_restart has no effect without steady\\.stage = simplec_newton"
                          "steady\\.switch_residual has no effect without steady\\.stage = simplec_newton" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                           steady.newton_gmres_restart=50 steady.switch_residual=1e-3 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_newton_keys_no_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "steady\\.newton_gmres_restart has no effect without steady\\.stage = simplec_newton" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_no_stage steady.newton_gmres_restart=50
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# Sub-phase 3.5: the 2-D hill (32 x 1 x 16, no Coriolis, van Albada) with
# steady.stage = simplec_newton switches to Newton once (at r <= 1e-3, about
# iteration 3860) and ends steady in at most 4500 iterations (3919 measured;
# SIMPLEC alone takes 5381). During the Newton stage 8 k-floor decisions
# move, Newton stalls on that difference and falls back after two failed
# steps, and SIMPLEC finishes (the final pass). About 10 s on two ranks.
set(HILL2D_NEWTON_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_steady_hill2d_newton)
ablsteady_add_test(NAME rt_steady_hill2d_newton NP 2 LABELS "regression;smoke" RUN_DIR ${HILL2D_NEWTON_DIR}
                   PREFIX ${EXPECT_ABORT} --success
                          "switching from stage simplec to newton: relative residual"
                          "ablsteady: steady after (3[0-9][0-9][0-9]|4[0-4][0-9][0-9]) iterations"
                          "ablsteady: it [0-9]+ newton r" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL2D_CASE}/inputs_steady_hill2d_newton ${ABLSTEADY_TEST_AMREX_FLAGS})
# The LU-SGS overlap against the coarse levels' boxes: overlap 3 fits level 0's
# 16-wide boxes but not level 3's, 2 wide.
ablsteady_add_test(NAME abort_sequencing_overlap NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "\\(on grid hierarchy level 3, steady\\.sequencing = 1\\)" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_sequencing_overlap
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# A restart with steady.sequencing = 1 continues level 0 and says so: three
# sequenced iterations of the truly neutral hill, then three more from the
# checkpoint (seconds).
set(SEQ_RESTART_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_sequencing_restart)
ablsteady_add_test(NAME rt_sequencing_restart_first NP 2 LABELS "regression;smoke" RUN_DIR ${SEQ_RESTART_DIR}
                   FIXTURES_SETUP sequencing_restart
                   PREFIX ${EXPECT_ABORT} --success "ablsteady: not steady after 3 iterations" "grid sequencing summary"
                          "WARNING: level 1: stopped at steady\\.max_iterations = 3 before the loose stop" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_sequencing_restart steady.max_iterations=3
                           init.first_guess=profile io.chk_file=seq_chk ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_sequencing_restart NP 2 LABELS "regression;smoke" RUN_DIR ${SEQ_RESTART_DIR}
                   FIXTURES_REQUIRED sequencing_restart
                   PREFIX ${EXPECT_ABORT} --success "ablsteady: restarted from seq_chk00003"
                          "grid sequencing: the restart continues on level 0; the coarse levels are not run"
                          "ablsteady: not steady after 3 iterations \\(6 over the run\\)"
                          "driver history restored: 3 iterations, cfl .*, reference hold 0 iterations \\(released" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_sequencing_restart steady.max_iterations=6
                           io.restart=seq_chk00003 io.chk_file=seq_chk_b ${ABLSTEADY_TEST_AMREX_FLAGS})
# Old against new steady.sequencing_hold on the same two-level deck (only
# level 0 is seeded; seconds): with the hold (1, sub-phase 2.5) level 0 runs
# its first iterations at steady.cfl0 with r = 1; without it (0, the
# default) r is below 1 from the first iteration and the CFL leaves cfl0 at
# the second; with steady.sequencing_cfl_start = previous level 0 starts at
# the CFL level 1 ended with (steady.cfl_max = 4 after its 150 iterations).
ablsteady_add_test(NAME rt_sequencing_hold_kept NP 2 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success "every level holds steady\\.cfl0 for the reference iterations"
                          "ablsteady: it 1 cfl 1\\.000e\\+00 r 1\\.000e\\+00" "ablsteady: it 3 cfl 1\\.000e\\+00 r 1\\.000e\\+00" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_sequencing_restart steady.max_iterations=3
                           init.first_guess=profile steady.sequencing_hold=1 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_sequencing_hold_released NP 2 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success "the seeded levels without the reference hold"
                          "ablsteady: it 1 cfl 1\\.000e\\+00 r [1-9]\\.[0-9]+e-0[1-9]" "!ablsteady: it 2 cfl 1\\.000e\\+00"
                          "ablsteady: it 2 cfl [1-9]\\.[0-9]+e\\+00 r [1-9]\\.[0-9]+e-0[1-9]" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_sequencing_restart steady.max_iterations=3
                           init.first_guess=profile ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_sequencing_cfl_previous NP 2 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success "starting at the previous level's CFL"
                          "ablsteady: level 1: it 150 cfl 4\\.000e\\+00" "ablsteady: it 1 cfl 4\\.000e\\+00 r [1-9]\\.[0-9]+e-0[1-9]" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_sequencing_restart steady.max_iterations=150
                           init.first_guess=profile steady.sequencing_cfl_start=previous ${ABLSTEADY_TEST_AMREX_FLAGS})
# simplec.w_vnn_cap over terrain above the tested aspect ratio: without the
# slope margin (simplec.w_vnn_cap_slope_margin = 0, #32's cap) the run goes on
# with a warning; with the default margin there is none. Below its floor of 16
# the cap is refused, and so is a margin below 1 or given without the cap.
ablsteady_add_test(NAME warn_w_vnn_cap_terrain NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success
                          "WARNING: simplec\\.w_vnn_cap = 16 over terrain with max aspect ratio .* above the tested 50\\.78"
                          "simplec\\.w_vnn_cap_slope_margin = 0" "without the cap \\(simplec\\.w_vnn_cap = 0\\) w converges"
                          "ablsteady: not steady after 1 iterations" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain simplec.w_vnn_cap=16
                           simplec.w_vnn_cap_slope_margin=0 mesh.dz0=0.7 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME warn_w_vnn_cap_terrain_margin NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success
                          "max aspect ratio max\\(dx,dy\\) / thickness += [1-9]\\.[0-9]+e\\+02"
                          "!WARNING: simplec\\.w_vnn_cap" "w_vnn_cap 16 \\(slope margin 4\\)" "ablsteady: not steady after 1 iterations" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain simplec.w_vnn_cap=16
                           mesh.dz0=0.7 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_w_vnn_cap_slope_margin NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "simplec\\.w_vnn_cap_slope_margin = 0\\.50+: must be 0 \\(off\\), or finite and >= 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain simplec.w_vnn_cap=16
                           simplec.w_vnn_cap_slope_margin=0.5 mesh.dz0=0.7 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_w_vnn_cap_slope_margin_alone NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "simplec\\.w_vnn_cap_slope_margin has no effect unless simplec\\.w_vnn_cap > 0" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain
                           simplec.w_vnn_cap_slope_margin=2 mesh.dz0=0.7 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME warn_w_vnn_cap_terrain_below NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} --success
                          "max aspect ratio max\\(dx,dy\\) / thickness += [1-4]\\.[0-9]+e\\+01"
                          "!WARNING: simplec\\.w_vnn_cap"
                          "ablsteady: not steady after 1 iterations" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain simplec.w_vnn_cap=16
                           mesh.dz0=4. ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_w_vnn_cap_floor NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "simplec\\.w_vnn_cap = 4\\.0+: must be 0 \\(off\\), or finite and >= 16" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_w_vnn_cap_terrain simplec.w_vnn_cap=4 mesh.dz0=0.7
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_repeated_key NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.n_cell is given more than once"
                          "io\\.write_plotfile takes exactly one value"
                          "io\\.restart takes exactly one value" --
                   COMMAND $<TARGET_FILE:ablsteady> ${VALIDATION_CASE}/inputs_abort_repeated_key
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# The turb.* and bc.surface_* readers stop at their first error, so one
# test each: the deck sets the key and the command line repeats it.
ablsteady_add_test(NAME abort_repeated_turb_key NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "turb\\.max_geom_lscale is given more than once" --
                   COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column turb.max_geom_lscale=5
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_repeated_surface_key NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "bc\\.surface_z0 is given more than once" --
                   COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column bc.surface_z0=0.2
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_restart_mesh NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke
                   PREFIX ${EXPECT_ABORT}
                          "io\\.restart = smoke_chk00000: checkpoint mesh\\.n_cell" --
                   COMMAND $<TARGET_FILE:ablsteady> ${SMOKE_CASE}/inputs_abort_restart_mesh
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_restart_overwrite NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke
                   PREFIX ${EXPECT_ABORT}
                          "io\\.chk_file = smoke_chk would overwrite io\\.restart = smoke_chk00000" --
                   COMMAND $<TARGET_FILE:ablsteady> ${SMOKE_CASE}/inputs_abort_restart_overwrite
                           ${ABLSTEADY_TEST_AMREX_FLAGS})

# --- Mesh (sub-phase 1.4) -----------------------------------------------------
set(CHECK_MESH ${Python3_EXECUTABLE} ${SUPPORT}/check_mesh.py)
# The mesh golds are double-precision values: a single-precision build
# compares the plotfile at a single-precision tolerance and skips the
# statistics checks, whose golds and tolerances are double-precision ones.
if(ABLSTEADY_PRECISION STREQUAL "SINGLE")
    set(FCOMPARE_RTOL 1e-5)
else()
    set(FCOMPARE_RTOL 1e-12)
    # RT-1.4.1 and RT-1.4.2: ERF's Askervein mesh, coarsened; field
    # statistics of z and detJ and the quality report against
    # Gold/askervein_mesh.json.
    ablsteady_add_test(NAME rt_askervein_mesh NP 2 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${MESH_CASE}/Gold/askervein_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${MESH_CASE}/inputs_mesh_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.4.2: the quality report of a steep ridge, with the STF warning.
    ablsteady_add_test(NAME rt_steep_ridge_quality NP 1 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${MESH_CASE}/Gold/steep_ridge_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${MESH_CASE}/inputs_mesh_steep_ridge
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()
# RT-2.3.1 and RT-2.3.2 (sub-phase 2.3): the Askervein mesh with STF and a
# three-level x/y grid hierarchy (60, 30, 15 cells a side; the same 18 layers).
# In double precision, level 0's statistics and every level's quality report
# against Gold/askervein_hierarchy_mesh.json; in both precisions, the z and
# detJ of levels 1 and 2 against their gold plotfiles (fcompare).
set(ASKH_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_askervein_hierarchy)
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(_askh_prefix ${CHECK_MESH} ${HIERARCHY_CASE}/Gold/askervein_hierarchy_mesh.json --)
else()
    set(_askh_prefix "")
endif()
ablsteady_add_test(NAME rt_askervein_hierarchy NP 2 LABELS "regression;smoke" RUN_DIR ${ASKH_DIR}
                   FIXTURES_SETUP askervein_hierarchy PREFIX ${_askh_prefix}
                   COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_mesh_askervein_hierarchy
                           terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
foreach(lev 1 2)
    ablsteady_add_test(NAME rt_askervein_hierarchy_level${lev}_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${ASKH_DIR} FIXTURES_REQUIRED askervein_hierarchy
                       COMMAND ${FCOMPARE} --rel_tol ${FCOMPARE_RTOL} askervein_hier_plt_level${lev}_00000
                               ${HIERARCHY_CASE}/Gold/askervein_hierarchy_level${lev}_plt00000)
endforeach()
# A restart rebuilds the hierarchy: a checkpoint written without one,
# continued with three levels on another layout (4 boxes from 2 ranks to 3),
# gives the coarse levels of the fresh run's golds.
set(HRS_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_hierarchy_restart)
ablsteady_add_test(NAME rt_hierarchy_restart_first NP 2 LABELS "regression;smoke" RUN_DIR ${HRS_DIR}
                   FIXTURES_SETUP hierarchy_restart_first
                   COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_hierarchy_restart_first
                           terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_hierarchy_restart NP 3 LABELS "regression;smoke" RUN_DIR ${HRS_DIR}
                   FIXTURES_REQUIRED hierarchy_restart_first FIXTURES_SETUP hierarchy_restart
                   COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_hierarchy_restart
                           terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt ${ABLSTEADY_TEST_AMREX_FLAGS})
foreach(lev 1 2)
    ablsteady_add_test(NAME rt_hierarchy_restart_level${lev}_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${HRS_DIR} FIXTURES_REQUIRED hierarchy_restart
                       COMMAND ${FCOMPARE} --rel_tol ${FCOMPARE_RTOL} hier_restart_plt_level${lev}_00000
                               ${HIERARCHY_CASE}/Gold/askervein_hierarchy_level${lev}_plt00000)
endforeach()
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    # RT-2.3.2: the same mesh with unsmoothed coarse grounds: the coarse
    # levels cap STF more than level 0 and the report gives the hint.
    ablsteady_add_test(NAME rt_askervein_hierarchy_hint NP 2 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${HIERARCHY_CASE}/Gold/askervein_hierarchy_p0_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_mesh_askervein_hierarchy_p0
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.3.2: the steep ridge with a three-level hierarchy: the report of
    # every level, with the STF warnings and the coarse-level hint.
    ablsteady_add_test(NAME rt_steep_ridge_hierarchy_quality NP 1 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${HIERARCHY_CASE}/Gold/steep_ridge_hierarchy_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_mesh_steep_ridge_hierarchy
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()
# UT-2.3.5 at the application: a hierarchy that does not divide the mesh, a
# box size below the coarsening and passes out of range abort, naming the keys.
ablsteady_add_test(NAME abort_hierarchy NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.hierarchy_levels = 3 needs mesh\\.n_cell x and y divisible by 2\\^2 = 4 \\(got 30 32\\)"
                          "mesh\\.max_grid_size must be a multiple of 2\\^\\(mesh\\.hierarchy_levels - 1\\) = 4"
                          "terrain\\.hierarchy_passes must be an integer from 0 to 64 \\(got -1\\)" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_abort_hierarchy
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_hierarchy_passes_unused NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "terrain\\.hierarchy_passes has no effect unless mesh\\.hierarchy_levels > 1" --
                   COMMAND $<TARGET_FILE:ablsteady> ${HIERARCHY_CASE}/inputs_abort_hierarchy_unused
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# A test-only tool works on level 0 alone and refuses a hierarchy.
ablsteady_add_test(NAME abort_hierarchy_tool NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT} "mesh\\.hierarchy_levels has no effect in wall_distance" --
                   COMMAND $<TARGET_FILE:wall_distance> ${HIERARCHY_CASE}/inputs_abort_hierarchy_tool
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# RT-1.4.3: a run without a solve writes the mapped-grid plotfile of the 2-D
# hill; it matches Gold/hill2d_plt00000 (relative tolerance 1e-12: compilers
# that fuse multiply-adds differ from those that do not in the last bit).
set(HILL2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_hill2d)
ablsteady_add_test(NAME rt_hill2d_run NP 2 LABELS "regression;smoke" RUN_DIR ${HILL2D_DIR}
                   FIXTURES_SETUP hill2d
                   COMMAND $<TARGET_FILE:ablsteady> ${MESH_CASE}/inputs_mesh_hill2d ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_hill2d_matches NP 1 LABELS "regression;smoke" RUN_DIR ${HILL2D_DIR}
                   FIXTURES_REQUIRED hill2d
                   COMMAND ${FCOMPARE} --rel_tol ${FCOMPARE_RTOL} hill2d_plt00000 ${MESH_CASE}/Gold/hill2d_plt00000)

# --- Poisson layer (sub-phase 1.10) ---------------------------------------------
# RT-1.10.1: a fixed wind projected on the Askervein mesh; iterations and
# divergence against Gold/poisson_askervein.json, the potential against
# Gold/poisson_askervein_plt00000. The golds are double-precision runs, so a
# single-precision build skips both (its unit tests cover the layer).
# The reader of solver.poisson_* stops on a repeated key and on a value out
# of range (the deck sets poisson_reltol; the command line repeats it).
ablsteady_add_test(NAME abort_poisson_options NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "solver\\.poisson_reltol is given more than once"
                          "solver\\.poisson_reltol must be >=" --
                   COMMAND $<TARGET_FILE:poisson_projection> ${POISSON_CASE}/inputs_poisson_askervein
                           terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                           solver.poisson_reltol=0 ${ABLSTEADY_TEST_AMREX_FLAGS})
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(POISSON_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_poisson_askervein)
    ablsteady_add_test(NAME rt_poisson_askervein NP 2 LABELS "regression;smoke" RUN_DIR ${POISSON_DIR}
                       FIXTURES_SETUP poisson_askervein
                       PREFIX ${Python3_EXECUTABLE} ${POISSON_CASE}/check_poisson.py
                              ${POISSON_CASE}/Gold/poisson_askervein.json --
                       COMMAND $<TARGET_FILE:poisson_projection> ${POISSON_CASE}/inputs_poisson_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # The potential converges to solver.poisson_reltol = 1e-10; 1e-8 allows
    # for compilers that round differently (fused multiply-adds).
    ablsteady_add_test(NAME rt_poisson_askervein_matches NP 1 LABELS "regression;smoke" RUN_DIR ${POISSON_DIR}
                       FIXTURES_REQUIRED poisson_askervein
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 poisson_askervein_plt00000
                               ${POISSON_CASE}/Gold/poisson_askervein_plt00000)
    # RT-1.11.1: the same projection with back end B (MLMG). It solves the same
    # discrete operator, so its potential matches back end A's gold.
    set(POISSON_B_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_poisson_askervein_mlmg)
    ablsteady_add_test(NAME rt_poisson_askervein_mlmg NP 2 LABELS "regression;smoke" RUN_DIR ${POISSON_B_DIR}
                       FIXTURES_SETUP poisson_askervein_mlmg
                       PREFIX ${Python3_EXECUTABLE} ${POISSON_CASE}/check_poisson.py
                              ${POISSON_CASE}/Gold/poisson_askervein_mlmg.json --
                       COMMAND $<TARGET_FILE:poisson_projection> ${POISSON_CASE}/inputs_poisson_askervein_mlmg
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_poisson_askervein_mlmg_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${POISSON_B_DIR} FIXTURES_REQUIRED poisson_askervein_mlmg
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 poisson_askervein_mlmg_plt00000
                               ${POISSON_CASE}/Gold/poisson_askervein_plt00000)
    # RT-1.12.1: the mass-consistent initial field on the Askervein mesh (log-law
    # first guess, alpha_v = 0.5): iterations, divergence and field extrema
    # against Gold/init_askervein.json, the state against
    # Gold/init_askervein_plt00000 (1e-8 as the Poisson golds). The gold is a
    # double-precision run, so a single-precision build skips it.
    set(INIT_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_init_askervein)
    ablsteady_add_test(NAME rt_init_askervein NP 2 LABELS "regression;smoke" RUN_DIR ${INIT_DIR}
                       FIXTURES_SETUP init_askervein
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_init.py askervein
                              ${INIT_CASE}/Gold/init_askervein.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${INIT_CASE}/inputs_init_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_init_askervein_matches NP 1 LABELS "regression;smoke" RUN_DIR ${INIT_DIR}
                       FIXTURES_REQUIRED init_askervein
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 init_askervein_plt00000
                               ${INIT_CASE}/Gold/init_askervein_plt00000)
    # RT-1.12.2: the hilltop speed-up over FastWindTerrain's bump for alpha_v =
    # 1, 0.7, 0.5, 0.3, 0.1 (five runs) against FastWindTerrain's own on the
    # same hill and profile (Gold/init_bump_trend.json: ranking, sign, spread).
    ablsteady_add_test(NAME rt_init_bump_trend NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_init.py trend
                              ${BUMP_CASE}/Gold/init_bump_trend.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${BUMP_CASE}/inputs_init_bump
                               terrain.file=${BUMP_CASE}/terrain_bump_fwt.txt ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.11.2: the Tucker wall distance on ERF's Neutral_Hill_2D mesh: ERF's
    # wall-distance checks (Tests/Discretization/Wall_Distance/check_walldist.py), and the field against
    # Gold/walldist_hill2d_plt00000, at double-precision tolerances (double
    # builds only).
    set(WALLDIST_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_walldist_hill2d)
    ablsteady_add_test(NAME rt_walldist_hill2d NP 2 LABELS "regression;smoke" RUN_DIR ${WALLDIST_DIR}
                       FIXTURES_SETUP walldist_hill2d
                       PREFIX ${Python3_EXECUTABLE} ${WALLDIST_CASE}/check_walldist.py
                              --hmax 100 --half-width 500 --
                       COMMAND $<TARGET_FILE:wall_distance> ${WALLDIST_CASE}/inputs_walldist_hill2d
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_walldist_hill2d_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${WALLDIST_DIR} FIXTURES_REQUIRED walldist_hill2d
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 walldist_hill2d_plt00000
                               ${WALLDIST_CASE}/Gold/walldist_hill2d_plt00000)
endif()

# --- RT-1.8.1: boundary fluxes of a fixed log-law state (sub-phase 1.8) ------
# Double precision only (the gold is double).
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(BC_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_bc_neutral_flat)
    ablsteady_add_test(NAME rt_bc_neutral_flat NP 2 LABELS "regression;smoke" RUN_DIR ${BC_DIR}
                       FIXTURES_SETUP bc_neutral_flat
                       PREFIX ${Python3_EXECUTABLE} ${BC_CASE}/check_bc.py
                              --ustar 0.35 --z0 0.1 --direction 15 --dz 15.625 --
                       COMMAND $<TARGET_FILE:boundary_fluxes> ${BC_CASE}/inputs_bc_neutral_flat
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_bc_neutral_flat_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${BC_DIR} FIXTURES_REQUIRED bc_neutral_flat
                       COMMAND ${FCOMPARE} --rel_tol 1e-10 bc_neutral_flat_plt00000
                               ${BC_CASE}/Gold/bc_neutral_flat_plt00000)
endif()

# --- RT-1.9.1: the neutral log-law column with the k-L closure (sub-phase 1.9)
# Double precision only (the gold is double).
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(CL_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_closure_column)
    ablsteady_add_test(NAME rt_closure_column NP 2 LABELS "regression;smoke" RUN_DIR ${CL_DIR}
                       FIXTURES_SETUP closure_column
                       PREFIX ${Python3_EXECUTABLE} ${CLOSURE_CASE}/check_closure.py
                              --ustar 0.35 --z0 0.1 --direction 15 --dz 15.625 --
                       COMMAND $<TARGET_FILE:closure_column> ${CLOSURE_CASE}/inputs_closure_column
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_closure_column_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${CL_DIR} FIXTURES_REQUIRED closure_column
                       COMMAND ${FCOMPARE} --rel_tol 1e-10 closure_column_plt00000
                               ${CLOSURE_CASE}/Gold/closure_column_plt00000)
endif()

# --- RT-1.15.1 / RT-1.15.2: the SIMPLEC stage to steady (sub-phase 1.15) -----
# The 1-D neutral column driven by a pressure-gradient force: the steady
# verdict, u* = sqrt(F H) exactly, the wall law in the first cell, the log law
# within 1 % (second cell) and 2 % (lowest eight), the wall k, nu_t, a residual
# drop of six orders (the Driver history) and identical columns. Double
# precision only (six orders and 1e-8 are below single's floor); the restart
# chain below runs in both.
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    ablsteady_add_test(NAME rt_steady_column NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${COLUMN_CASE}/check_steady_column.py
                              --ustar 0.35 --z0 0.1 --force 1.225e-4 --height 1000 --dz 15.625 --orders 6 --
                       COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column ${ABLSTEADY_TEST_AMREX_FLAGS})
    # The column with turb.kappa and turb.Cmu0 off their defaults: the wall
    # model and the inflow profile take the closure's constants, so the first
    # cell is on the wall law with kappa = 0.40 and the wall k is u*^2/0.55^2.
    # Before, they kept 0.41 and 0.5562 (2.4 % and 2.2 % off). The pair is
    # chosen inside the swept range kappa 0.38-0.42, Cmu0 0.53-0.58, where the
    # discretization checks (log law over k < 8 within 2 %, nu_t within 6 %)
    # keep their margin; the wall law and wall k hold to round-off there.
    ablsteady_add_test(NAME rt_steady_column_constants NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${COLUMN_CASE}/check_steady_column.py
                              --ustar 0.35 --z0 0.1 --force 1.225e-4 --height 1000 --dz 15.625 --orders 6
                              --kappa 0.40 --cmu0 0.55 --
                       COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column turb.kappa=0.40
                               turb.Cmu0=0.55 ${ABLSTEADY_TEST_AMREX_FLAGS})
    # Beyond the plan: the rough column (z0 = 0.5 m, 1 m layers, the first
    # centre below (e - 1) z0). The default floor follows the log law in the
    # first cell, and the profile of the shallow column's stress u*^2 (1 - z/H)
    # above it; the floor one (sub-phase 1.7's form) gives 1/ln 2 = 1.44 times
    # the log law's speed in the first cell.
    ablsteady_add_test(NAME rt_steady_column_rough NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${COLUMN_CASE}/check_steady_column.py
                              --ustar 0.3535533905932738 --z0 0.5 --force 1.953125e-3 --height 64 --dz 1 --orders 6
                              --stress-profile --
                       COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column_rough ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_steady_column_rough_floor_one NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${COLUMN_CASE}/check_steady_column.py
                              --ustar 0.3535533905932738 --z0 0.5 --force 1.953125e-3 --height 64 --dz 1 --orders 6
                              --floor one --
                       COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column_rough
                               bc.surface_denominator_floor=one ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()
# A restart continues the stage bit for bit: 40 iterations, continued to 80
# from the checkpoint (the driver and closure histories), against 80 straight.
set(SCR_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_steady_column_restart)
ablsteady_add_test(NAME rt_steady_column_restart_first NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}
                   FIXTURES_SETUP steady_column_first
                   COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column_restart steady.max_iterations=40
                           init.first_guess=profile ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_steady_column_restart_continued NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}
                   FIXTURES_REQUIRED steady_column_first FIXTURES_SETUP steady_column_continued
                   COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column_restart steady.max_iterations=80
                           io.restart=steady_column_restart_chk00040 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_steady_column_restart_straight NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}/straight
                   FIXTURES_SETUP steady_column_straight
                   COMMAND $<TARGET_FILE:ablsteady> ${COLUMN_CASE}/inputs_steady_column_restart steady.max_iterations=80
                           init.first_guess=profile ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_steady_column_restart_matches NP 1 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}
                   FIXTURES_REQUIRED steady_column_continued steady_column_straight
                   COMMAND ${FCOMPARE} --rel_tol 1e-14 --abs_tol 1e-300 steady_column_restart_plt00080
                           straight/steady_column_restart_plt00080)
# ERF's Neutral_ABL_Flat to steady (about 13 s on two ranks); ERF's own
# check_neutral.py --physics on the turbulence plotfile (extern_erf/
# check_neutral.py). Regression label, not smoke (the run time).
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    ablsteady_add_test(NAME rt_steady_neutral_flat NP 2 LABELS "regression"
                       PREFIX ${Python3_EXECUTABLE} ${FLAT_CASE}/check_steady_neutral.py --
                       COMMAND $<TARGET_FILE:ablsteady> ${FLAT_CASE}/inputs_steady_neutral_flat ${ABLSTEADY_TEST_AMREX_FLAGS})
    # Beyond the plan: the w rows' VNN cap (simplec.w_vnn_cap = 16) on a flat
    # strip with a 0.7 m first cell (about 50 s on two ranks): the verdict
    # with and without the cap, at least 5 times fewer iterations with it,
    # u, v, w of the two steady states within 2e-4 m/s, and a restart of the
    # capped run bit for bit (check_w_vnn_cap.py runs all five).
    ablsteady_add_test(NAME rt_thin_cell_w_vnn_cap NP 2 LABELS "regression"
                       PREFIX ${Python3_EXECUTABLE} ${THIN_CELL_CASE}/check_w_vnn_cap.py
                              --cap 16 --min-ratio 5 --tol 2e-4 --
                       COMMAND $<TARGET_FILE:ablsteady> ${THIN_CELL_CASE}/inputs_thin_cell_flat ${ABLSTEADY_TEST_AMREX_FLAGS})
    set_tests_properties(rt_thin_cell_w_vnn_cap PROPERTIES ENVIRONMENT PYTHONDONTWRITEBYTECODE=1 TIMEOUT 300)
    # Beyond the plan: the slope margin of the cap on a steep ridge with a 0.4 m
    # first cell (Canonical_RANS/Steep_Ridge_Thin_Cell). With the default margin
    # the run is steady within 3500 iterations; with the margin off (#32's cap)
    # the pressure correction stalls and it is not.
    ablsteady_add_test(NAME rt_steep_ridge_w_vnn_cap NP 2 LABELS "regression"
                       PREFIX ${EXPECT_ABORT} --success "w_vnn_cap 16 \\(slope margin 4\\)"
                              "ablsteady: steady after [0-9]+ iterations" --
                       COMMAND $<TARGET_FILE:ablsteady> ${STEEP_RIDGE_CASE}/inputs_steep_ridge_thin_cell simplec.w_vnn_cap=16
                               steady.max_iterations=3500 ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_steep_ridge_w_vnn_cap_margin_off NP 2 LABELS "regression"
                       PREFIX ${EXPECT_ABORT} --success "w_vnn_cap 16 \\(slope margin off\\)"
                              "ablsteady: not steady after 3500 iterations" --
                       COMMAND $<TARGET_FILE:ablsteady> ${STEEP_RIDGE_CASE}/inputs_steep_ridge_thin_cell simplec.w_vnn_cap=16
                               simplec.w_vnn_cap_slope_margin=0 steady.max_iterations=3500 ${ABLSTEADY_TEST_AMREX_FLAGS})
    set_tests_properties(rt_steep_ridge_w_vnn_cap rt_steep_ridge_w_vnn_cap_margin_off PROPERTIES TIMEOUT 600)
    # RT-1.16.1: ERF's Neutral_Hill_2D to steady (about 80 s on two ranks);
    # the verdict, a two-order drop of the u residual and ERF's own
    # check_hill.py --physics (extern_erf/) on the turbulence plotfile.
    ablsteady_add_test(NAME rt_steady_hill2d NP 2 LABELS "regression"
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_steady_hill.py
                              --erf-check check_hill.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL2D_CASE}/inputs_steady_hill2d ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.1.1: the same case with the line LU-SGS sweeps (two symmetric
    # sweeps, CFL cap 16): the verdict, the residual drop and ERF's checks.
    set(HILL_COL_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_steady_hill2d_column)
    ablsteady_add_test(NAME rt_steady_hill2d_lusgs NP 2 LABELS "regression" RUN_DIR ${HILL_COL_DIR}
                       FIXTURES_SETUP hill2d_lusgs
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_steady_hill.py
                              --erf-check check_hill.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL2D_CASE}/inputs_steady_hill2d_lusgs ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.1b.1: the same case from the converged column (init.first_guess =
    # column): the verdict, the drop and ERF's checks; then the spin-up check
    # against the profile run's history (both checkpoints in HILL_COL_DIR).
    ablsteady_add_test(NAME rt_steady_hill2d_column NP 2 LABELS "regression" RUN_DIR ${HILL_COL_DIR}
                       FIXTURES_SETUP hill2d_column
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_steady_hill.py
                              --erf-check check_hill.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL2D_CASE}/inputs_steady_hill2d_column ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_column_guess_spinup
             COMMAND ${Python3_EXECUTABLE} ${HILL2D_CASE}/check_column_guess.py
                     steady_hill2d_lusgs_chk steady_hill2d_column_chk
             WORKING_DIRECTORY ${HILL_COL_DIR})
    set_tests_properties(rt_column_guess_spinup PROPERTIES LABELS "regression"
                         FIXTURES_REQUIRED "hill2d_lusgs;hill2d_column" ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # RT-1.16.2: ERF's Neutral_Hill_3D (64 x 64 x 20) to steady with
    # check_hill3d.py; minutes on four ranks, so label nightly (not in
    # ctest -L regression).
    ablsteady_add_test(NAME rt_steady_hill3d NP 4 LABELS "nightly"
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_steady_hill.py
                              --erf-check check_hill3d.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_CASE}/inputs_steady_hill3d ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.5.1 (the ERF half): the same deck with grid sequencing over three
    # levels (16 x 16, 32 x 32, then 64 x 64); the same checks and the same
    # two-order drop, measured from the coarsest level's reference.
    ablsteady_add_test(NAME rt_steady_hill3d_sequencing NP 4 LABELS "nightly"
                       PREFIX ${Python3_EXECUTABLE} ${SUPPORT}/check_steady_hill.py
                              --erf-check check_hill3d.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_CASE}/inputs_steady_hill3d mesh.hierarchy_levels=3
                               steady.sequencing=1 ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.2.1: a truly neutral Hill_3D (inputs_steady_hill3d_neutral) to
    # round-off with the z-line Jacobi sweeps and with the box-local line
    # LU-SGS sweeps (zero ghosts, overlap 2, diagonal ordering); then
    # check_line_match.py compares the two states to 1e-8: u, v, w, p and k
    # directly, theta after subtracting its mean (with no theta condition and
    # gravity 0 its constant is free and drifts by about 1e-4 K at round-off).
    # The setup step empties the run directories, so a plotfile left by an
    # earlier run cannot be compared.
    set(RT221_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_steady_hill3d_neutral)
    add_test(NAME rt_steady_hill3d_neutral_clean
             COMMAND ${Python3_EXECUTABLE} -c
                     "import shutil, os, sys; shutil.rmtree(sys.argv[1], ignore_errors=True); [os.makedirs(os.path.join(sys.argv[1], d)) for d in ('zjacobi', 'lusgs', 'lusgs_fft', 'sequencing')]"
                     ${RT221_DIR})
    set_tests_properties(rt_steady_hill3d_neutral_clean PROPERTIES LABELS "nightly" FIXTURES_SETUP hill3d_neutral_clean)
    ablsteady_add_test(NAME rt_steady_hill3d_neutral_zjacobi NP 4 LABELS "nightly" RUN_DIR ${RT221_DIR}/zjacobi
                       FIXTURES_REQUIRED hill3d_neutral_clean FIXTURES_SETUP hill3d_neutral_zjacobi
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # The default back end, gmres_mlmg, for the pressure correction and for the
    # alpha_v = alpha_h projection (both checked in the log).
    ablsteady_add_test(NAME rt_steady_hill3d_neutral_lusgs NP 4 LABELS "nightly" RUN_DIR ${RT221_DIR}/lusgs
                       FIXTURES_REQUIRED hill3d_neutral_clean FIXTURES_SETUP hill3d_neutral_lusgs
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" "pressure gmres_mlmg reltol"
                              "Poisson iterations \\(gmres_mlmg," --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral_lusgs
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_steady_hill3d_neutral_matches
             COMMAND ${Python3_EXECUTABLE} ${HILL3D_NEUTRAL_CASE}/check_line_match.py zjacobi lusgs
                     steady_hill3d_neutral_plt --rel_tol 1e-8
             WORKING_DIRECTORY ${RT221_DIR})
    set_tests_properties(rt_steady_hill3d_neutral_matches PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED "hill3d_neutral_zjacobi;hill3d_neutral_lusgs"
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # Old against new default of solver.poisson_backend: the same deck with
    # gmres_fft (the default before) reaches the gmres_mlmg run's steady state
    # to 1e-8 of each field's maximum (both run to round-off).
    ablsteady_add_test(NAME rt_poisson_default_gmres_fft NP 4 LABELS "nightly" RUN_DIR ${RT221_DIR}/lusgs_fft
                       FIXTURES_REQUIRED hill3d_neutral_clean FIXTURES_SETUP hill3d_neutral_lusgs_fft
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" "pressure gmres_fft reltol" --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral_lusgs
                               solver.poisson_backend=gmres_fft ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_poisson_default_matches
             COMMAND ${Python3_EXECUTABLE} ${HILL3D_NEUTRAL_CASE}/check_line_match.py lusgs lusgs_fft
                     steady_hill3d_neutral_plt --rel_tol 1e-8
             WORKING_DIRECTORY ${RT221_DIR})
    set_tests_properties(rt_poisson_default_matches PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED "hill3d_neutral_lusgs;hill3d_neutral_lusgs_fft"
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # RT-2.5.1: the same deck with grid sequencing over three levels (8 x 8,
    # 16 x 16, then 32 x 32 to round-off); the loose stops of levels 2 and 1
    # (by the residual branch) and the strict stop of level 0 are checked in
    # the log, level 0 in fewer than 3000 iterations (the z-line Jacobi run of
    # RT-2.2.1 takes 3495; sequencing 2688 (2773 with steady.sequencing_hold = 1), so a lost coarse start or
    # reference shows), and the steady state must match RT-2.2.1's to 1e-8;
    # p after subtracting each run's mean (no face fixes its level: the
    # prolonged start carries another constant, -1.8e-3 here).
    ablsteady_add_test(NAME rt_steady_hill3d_neutral_sequencing NP 4 LABELS "nightly" RUN_DIR ${RT221_DIR}/sequencing
                       FIXTURES_REQUIRED hill3d_neutral_clean FIXTURES_SETUP hill3d_neutral_sequencing
                       PREFIX ${EXPECT_ABORT} --success
                              "ablsteady: steady after ([0-9]{1,3}|[12][0-9]{3}) iterations"
                              "ablsteady: level 2: loose stop after [0-9]+ iterations \\(every residual 2\\.000e\\+00 orders below its reference\\)"
                              "ablsteady: level 1: loose stop after [0-9]+ iterations \\(every residual 2\\.000e\\+00 orders below its reference\\)"
                              "grid sequencing summary" --
                       COMMAND $<TARGET_FILE:ablsteady> ${HILL3D_NEUTRAL_CASE}/inputs_steady_hill3d_neutral
                               mesh.hierarchy_levels=3 steady.sequencing=1 ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_steady_hill3d_neutral_sequencing_matches
             COMMAND ${Python3_EXECUTABLE} ${HILL3D_NEUTRAL_CASE}/check_line_match.py zjacobi sequencing
                     steady_hill3d_neutral_plt --rel_tol 1e-8 --mean_free p
             WORKING_DIRECTORY ${RT221_DIR})
    set_tests_properties(rt_steady_hill3d_neutral_sequencing_matches PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED "hill3d_neutral_zjacobi;hill3d_neutral_sequencing"
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # RT-1.16.3 (the convergence half): Askervein on the ERF grid with an
    # inflow and an outflow face to the steady verdict (about 2 s on two
    # ranks), the Phase 1 sweep. The comparison with ERF's URANS on lines A
    # and AA is RT-2.6.1 (Askervein_210, below).
    # alpha_v != alpha_h: the projection falls back to gmres_fft (logged), the
    # stage keeps the default gmres_mlmg.
    ablsteady_add_test(NAME rt_steady_askervein NP 2 LABELS "regression;smoke"
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after"
                              "Poisson iterations \\(gmres_fft," "the projection used gmres_fft, the SIMPLEC stage keeps gmres_mlmg"
                              "pressure gmres_mlmg reltol" --
                       COMMAND $<TARGET_FILE:ablsteady> ${ASKERVEIN_CASE}/inputs_steady_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # Beyond the plan: Askervein in a conventionally neutral boundary layer.
    # The precursor column runs to steady; a companion run of the same deck
    # stops at the table's step, and tools/make_inflow_table.py checks its
    # 10 m wind is along +x (0.01 deg), u* = 0.6 m/s (1e-3) and the committed
    # inflow table matches it (1e-3, the same step). The 3-D case runs on that
    # table to the steady verdict.
    set(CNBL_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_precursor_askervein_cnbl)
    ablsteady_add_test(NAME rt_precursor_askervein_cnbl NP 2 LABELS "regression" RUN_DIR ${CNBL_DIR}
                       FIXTURES_SETUP precursor_askervein_cnbl
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${CNBL_CASE}/inputs_precursor_askervein_cnbl
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # The table check reads a companion run stopped at the table's own step
    # (read from the table's header, so a regenerated table carries its step
    # here), not the run above. Under the capping inversion theta has no
    # fixed point: the inversion base cools by 1.35e-4 K per iteration as it
    # rises and the mixed layer warms by 1.5e-5 K. The verdict comes when the
    # k residual leaves a limit cycle of about 600 iterations (the Poisson
    # iteration count switching between 5 and 6), at an iteration that depends
    # on the compiler's rounding: 11087 with Apple clang, 11073 with gcc-16,
    # theta 1.7e-3 K apart. At the fixed step the two builds agree to 2e-4 K,
    # so theta keeps 1e-3 K like u, v and k. A build that left the cycle one
    # cycle later would fail here, and the verdict run would show its step.
    # The deck is the precursor's with only the step limit and an unreachable
    # drop criterion changed, written at configure time.
    file(STRINGS ${CNBL_CASE}/askervein_cnbl_inflow.txt _cnbl_header LIMIT_COUNT 1 REGEX "planar means of")
    string(REGEX MATCH "_plt([0-9]+) " _cnbl_step_match "${_cnbl_header}")
    if(NOT _cnbl_step_match)
        message(FATAL_ERROR "askervein_cnbl_inflow.txt: no \"planar means of <plotfile><step> (\" header line")
    endif()
    set(CNBL_TABLE_STEP ${CMAKE_MATCH_1})
    file(READ ${CNBL_CASE}/inputs_precursor_askervein_cnbl _cnbl_deck)
    foreach(_kv "max_iterations          = ${CNBL_TABLE_STEP}" "residual_drop_orders    = 9")
        string(REGEX MATCH "^[a-z_]+" _k "${_kv}")
        string(REGEX MATCHALL "\nsteady\\.${_k} [^\n]*" _hits "${_cnbl_deck}")
        list(LENGTH _hits _n)
        if(NOT _n EQUAL 1)
            message(FATAL_ERROR "inputs_precursor_askervein_cnbl must set steady.${_k} exactly once (found ${_n})")
        endif()
        string(REGEX REPLACE "\nsteady\\.${_k} [^\n]*" "\nsteady.${_kv}" _cnbl_deck "${_cnbl_deck}")
    endforeach()
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/inputs_precursor_askervein_cnbl_fixed "${_cnbl_deck}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 ${CNBL_CASE}/inputs_precursor_askervein_cnbl ${CNBL_CASE}/askervein_cnbl_inflow.txt)
    set(CNBL_FIXED_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_precursor_askervein_cnbl_fixed)
    ablsteady_add_test(NAME rt_precursor_askervein_cnbl_fixed NP 2 LABELS "regression" RUN_DIR ${CNBL_FIXED_DIR}
                       FIXTURES_SETUP precursor_askervein_cnbl_fixed
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: not steady after ${CNBL_TABLE_STEP} iterations" --
                       COMMAND $<TARGET_FILE:ablsteady> ${CMAKE_CURRENT_BINARY_DIR}/inputs_precursor_askervein_cnbl_fixed
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_precursor_askervein_cnbl_table
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/make_inflow_table.py
                     precursor_askervein_cnbl_plt precursor_askervein_cnbl_turb
                     --ug 17.028798 --vg -7.380863 --check ${CNBL_CASE}/askervein_cnbl_inflow.txt
             WORKING_DIRECTORY ${CNBL_FIXED_DIR})
    set_tests_properties(rt_precursor_askervein_cnbl_table PROPERTIES LABELS "regression"
                         FIXTURES_REQUIRED precursor_askervein_cnbl_fixed)
    ablsteady_add_test(NAME rt_steady_askervein_cnbl NP 2 LABELS "regression"
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${CNBL_CASE}/inputs_steady_askervein_cnbl
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               bc.inflow_file=${CNBL_CASE}/askervein_cnbl_inflow.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.6.1 and RT-2.6.2 (sub-phase 2.6, the Phase 2 checkpoint): Askervein
    # with ERF's own setup (wind from 210 degrees, ERF's 300 x 300 x 18 grid,
    # Tests/Real_Terrain/Askervein_210). The inflow table must be the one
    # tools/erf_inflow_table.py makes from ERF's profile (extern_erf/).
    add_test(NAME rt_askervein_210_inflow_table
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/erf_inflow_table.py
                     ${PROJECT_SOURCE_DIR}/extern_erf/inflow_powerlaw.txt --theta 300
                     --check ${ASKERVEIN_210_CASE}/askervein_210_inflow.txt)
    set_tests_properties(rt_askervein_210_inflow_table PROPERTIES LABELS "regression;smoke"
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # The Phase 2 configuration: the CPU default sweep (lusgs, given here so
    # the gold holds on a GPU build too) with its default CFL cap, the default
    # reference window, grid sequencing over two levels (about 140
    # iterations, 2-3 min on four ranks). Its log is kept for RT-2.6.2. The
    # run directories are emptied first, so no check reads an earlier run's
    # plotfile or lines.
    set(ASK210_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_askervein_210)
    add_test(NAME rt_askervein_210_clean
             COMMAND ${Python3_EXECUTABLE} -c "import os, shutil, sys; [ (shutil.rmtree(d, ignore_errors=True), os.makedirs(d)) for d in sys.argv[1:] ]"
                     ${ASK210_DIR}/phase2 ${ASK210_DIR}/phase1 ${ASK210_DIR}/unbounded)
    set_tests_properties(rt_askervein_210_clean PROPERTIES LABELS "nightly" FIXTURES_SETUP askervein_210_clean)
    ablsteady_add_test(NAME rt_steady_askervein_210 NP 4 LABELS "nightly" RUN_DIR ${ASK210_DIR}/phase2
                       FIXTURES_SETUP askervein_210_phase2 FIXTURES_REQUIRED askervein_210_clean
                       PREFIX ${EXPECT_ABORT} --success --log ${ASK210_DIR}/phase2/run.log --timeout 3600
                              "ablsteady: steady after" "grid sequencing summary" --
                       COMMAND $<TARGET_FILE:ablsteady> ${ASKERVEIN_210_CASE}/inputs_steady_askervein_210
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               bc.inflow_file=${ASKERVEIN_210_CASE}/askervein_210_inflow.txt
                               steady.sweep=lusgs steady.sequencing=1 mesh.hierarchy_levels=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.6.1: 10 m wind speeds on lines A and AA within 5 % of ERF's URANS
    # average (2400-4800 s, Gold/erf_urans_lines_10m.csv), 10 % at the lee
    # masts more than 300 m past the hilltop (ablsteady's second-order
    # advection against ERF's fifth order sets the wake there).
    add_test(NAME rt_askervein_210_erf_lines
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/askervein_lines.py check
                     steady_askervein_210_plt ${ASKERVEIN_210_CASE}/Gold/erf_urans_lines_10m.csv
                     --terrain ${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                     --lines ${ASKERVEIN_210_CASE}/askervein_lines.csv --out lines.csv --lee_d 300 --lee_rel_tol 0.10
             WORKING_DIRECTORY ${ASK210_DIR}/phase2)
    set_tests_properties(rt_askervein_210_erf_lines PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED askervein_210_phase2 FIXTURES_SETUP askervein_210_lines
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # The same lines against their own gold to 1e-4 (Gold/askervein_210_lines.csv):
    # the ERF tolerances alone would pass the run with the inflow k of before.
    add_test(NAME rt_askervein_210_lines_gold
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/askervein_lines.py compare lines.csv
                     ${ASKERVEIN_210_CASE}/Gold/askervein_210_lines.csv --rel_tol 1e-4
             WORKING_DIRECTORY ${ASK210_DIR}/phase2)
    set_tests_properties(rt_askervein_210_lines_gold PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED askervein_210_lines ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # RT-2.6.2: the level history (cells, stop reasons, iterations within
    # 10 %) and the stage switches against the gold sequence.
    add_test(NAME rt_askervein_210_level_history
             COMMAND ${Python3_EXECUTABLE} ${SUPPORT}/check_sequencing_log.py ${ASK210_DIR}/phase2/run.log
                     ${ASKERVEIN_210_CASE}/Gold/askervein_210_level_history.json --iter_rel_tol 0.1)
    set_tests_properties(rt_askervein_210_level_history PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED askervein_210_phase2 ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # Old vs new: the Phase 1 configuration of RT-1.16.3 (zjacobi, four
    # sweeps, CFL cap 2, a 100-iteration reference, one level; about 540
    # iterations) reaches the same lines as the Phase 2 default to 0.1 %.
    ablsteady_add_test(NAME rt_steady_askervein_210_phase1 NP 4 LABELS "nightly" RUN_DIR ${ASK210_DIR}/phase1
                       FIXTURES_SETUP askervein_210_phase1 FIXTURES_REQUIRED askervein_210_clean
                       PREFIX ${EXPECT_ABORT} --success --timeout 3600 "ablsteady: steady after"
                              "4 zjacobi sweeps" --
                       COMMAND $<TARGET_FILE:ablsteady> ${ASKERVEIN_210_CASE}/inputs_steady_askervein_210
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               bc.inflow_file=${ASKERVEIN_210_CASE}/askervein_210_inflow.txt
                               steady.sweep=zjacobi steady.zline_sweeps=4 steady.cfl_max=2
                               steady.residual_ref_iterations=100 ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_askervein_210_phase1_matches
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/askervein_lines.py check
                     steady_askervein_210_plt ${ASK210_DIR}/phase2/lines.csv
                     --terrain ${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                     --lines ${ASKERVEIN_210_CASE}/askervein_lines.csv --out lines.csv --rel_tol 1e-3
             WORKING_DIRECTORY ${ASK210_DIR}/phase1)
    set_tests_properties(rt_askervein_210_phase1_matches PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED "askervein_210_phase1;askervein_210_lines"
                         ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    # Old against new for bc.inflow_tke_length: the Phase 2 run with the
    # inflow k of before (unbounded length) reaches the lines of that run
    # (Gold/askervein_210_lines_unbounded.csv, #43) to 1e-4.
    ablsteady_add_test(NAME rt_steady_askervein_210_unbounded NP 4 LABELS "nightly" RUN_DIR ${ASK210_DIR}/unbounded
                       FIXTURES_SETUP askervein_210_unbounded FIXTURES_REQUIRED askervein_210_clean
                       PREFIX ${EXPECT_ABORT} --success --timeout 3600 "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${ASKERVEIN_210_CASE}/inputs_steady_askervein_210
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               bc.inflow_file=${ASKERVEIN_210_CASE}/askervein_210_inflow.txt
                               bc.inflow_tke_length=unbounded
                               steady.sweep=lusgs steady.sequencing=1 mesh.hierarchy_levels=2 ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_askervein_210_unbounded_matches
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/askervein_lines.py check
                     steady_askervein_210_plt ${ASKERVEIN_210_CASE}/Gold/askervein_210_lines_unbounded.csv
                     --terrain ${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                     --lines ${ASKERVEIN_210_CASE}/askervein_lines.csv --out lines.csv --rel_tol 1e-4
             WORKING_DIRECTORY ${ASK210_DIR}/unbounded)
    set_tests_properties(rt_askervein_210_unbounded_matches PROPERTIES LABELS "nightly"
                         FIXTURES_REQUIRED askervein_210_unbounded ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
endif()

# --- Reference archives (sub-phase 1.1): ablsteady_refavg end to end --------
# Three copies of the smoke plotfile at t = 0, 5 and 10 s: their trapezoid
# average (weights 1/4, 1/2, 1/4, fma) reproduces the plotfile bit for bit.
# This checks the archive round trip; the weights themselves are tested in
# the unit tests (identical snapshots cannot tell weights apart).
set(REFAVG_FIX ${Python3_EXECUTABLE} ${REFAVG_CASE}/refavg_fixtures.py)
set(REFAVG_META refavg.case=Smoke refavg.erf_sha=92410d578f612ab9e844c450bf8e900194becfb5
                refavg.machine=kestrel-cpu refavg.nodes=2 refavg.ranks=8 refavg.wall_seconds=1800)
add_test(NAME refavg_snapshots
         COMMAND ${REFAVG_FIX} snapshots smoke_plt00000 refsnap_a:0 refsnap_b:5 refsnap_c:10
         WORKING_DIRECTORY ${SMOKE_DIR})
set_tests_properties(refavg_snapshots PROPERTIES LABELS "regression;smoke"
                     FIXTURES_REQUIRED smoke FIXTURES_SETUP refavg_snapshots
                     ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
ablsteady_add_test(NAME refavg_build NP 2 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED refavg_snapshots FIXTURES_SETUP refavg_archive
                   COMMAND $<TARGET_FILE:ablsteady_refavg> "refavg.plotfiles=refsnap_a refsnap_b refsnap_c/"
                           refavg.output=ref_archive
                           "refavg.inputs_files=${SMOKE_CASE}/inputs_smoke ${SMOKE_CASE}/inputs_restart"
                           ${REFAVG_META} ${ABLSTEADY_TEST_AMREX_FLAGS})
# 2 nodes x 0.5 h x 10 AU = 10 AU
ablsteady_add_test(NAME refavg_verify NP 3 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED refavg_archive
                   PREFIX ${EXPECT_ABORT} --success
                          "ref_archive verifies: case Smoke, ERF 92410d578, 3 snapshots, kestrel-cpu, 2 node\\(s\\), 8 rank\\(s\\), 10.000 AU"
                          "theta +half-window change 0.000e\\+00" --
                   COMMAND $<TARGET_FILE:ablsteady_refavg> refavg.verify=ref_archive
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME refavg_matches NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED refavg_archive
                   COMMAND ${FCOMPARE} --allow_diff_grids smoke_plt00000 ref_archive/average)
add_test(NAME refavg_damage
         COMMAND ${REFAVG_FIX} damage ref_archive ref_archive_damaged
         WORKING_DIRECTORY ${SMOKE_DIR})
set_tests_properties(refavg_damage PROPERTIES LABELS "regression;smoke"
                     FIXTURES_REQUIRED refavg_archive FIXTURES_SETUP refavg_damaged
                     ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
ablsteady_add_test(NAME abort_refavg_damaged NP 2 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED refavg_damaged
                   PREFIX ${EXPECT_ABORT}
                          "reference ref_archive_damaged does not verify"
                          "average/Level_0/Cell_D_00000: SHA-256 differs from MANIFEST" --
                   COMMAND $<TARGET_FILE:ablsteady_refavg> refavg.verify=ref_archive_damaged
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_refavg_inputs NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "refavg\\.machine must be kestrel-cpu, kestrel-gpu or local \\(got kestrel\\)"
                          "refavg\\.nodes must be >= 1"
                          "refavg\\.erf_sha must be a full 40-character"
                          "refavg\\.plotfiles needs at least 3 plotfiles"
                          "'refavg\\.typo': not an input of namespace 'refavg'" --
                   COMMAND $<TARGET_FILE:ablsteady_refavg> refavg.plotfiles=only_one refavg.output=x
                           refavg.case=C refavg.erf_sha=92410d5 refavg.machine=kestrel refavg.nodes=0
                           refavg.typo=1 ${ABLSTEADY_TEST_AMREX_FLAGS})

ablsteady_add_test(NAME abort_refavg_verify_keys NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED refavg_archive
                   PREFIX ${EXPECT_ABORT}
                          "refavg\\.output has no effect with refavg\\.verify; remove it" --
                   COMMAND $<TARGET_FILE:ablsteady_refavg> refavg.verify=ref_archive refavg.output=x
                           ${ABLSTEADY_TEST_AMREX_FLAGS})

# --- external/VERSIONS and the Kestrel scripts --------------------------------
add_test(NAME versions_check
         COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/check_versions.py)
set_tests_properties(versions_check PROPERTIES LABELS unit ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
add_test(NAME versions_check_selftest
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_check_versions
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tools)
set_tests_properties(versions_check_selftest PROPERTIES LABELS unit ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
# The Kestrel scripts cannot run here; at least they must parse. They are
# bash scripts for a Linux cluster, so they are not registered on Windows.
if(NOT WIN32)
    file(GLOB KESTREL_SCRIPTS ${PROJECT_SOURCE_DIR}/tools/kestrel/*.sh)
    foreach(script IN LISTS KESTREL_SCRIPTS)
        get_filename_component(sname ${script} NAME_WE)
        add_test(NAME kestrel_parse_${sname} COMMAND bash -n ${script})
        set_tests_properties(kestrel_parse_${sname} PROPERTIES LABELS unit)
    endforeach()
    # ... and a dry run with stub module/sbatch/srun exercises submit and archive.
    add_test(NAME kestrel_dryrun
             COMMAND bash ${KESTREL_CASE}/kestrel_dryrun.sh ${PROJECT_SOURCE_DIR}
                     $<TARGET_FILE:ablsteady_refavg> ${SMOKE_DIR}/smoke_plt00000
                     ${CMAKE_CURRENT_BINARY_DIR}/run/kestrel_dryrun)
    set_tests_properties(kestrel_dryrun PROPERTIES LABELS "regression;smoke" FIXTURES_REQUIRED smoke)
endif()

# The plotfile reader of the check scripts (extern_erf/erf_plotfile.py, a
# local change): hand-made one-box plotfiles in float and double, little- and
# big-endian, read back exactly.
add_test(NAME erf_plotfile_selftest
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_erf_plotfile
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tools)
set_tests_properties(erf_plotfile_selftest PROPERTIES LABELS unit ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)

# --- RT-1.3.1 helper: the extern_erf drift check's own tests -----------------
add_test(NAME drift_check_selftest
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_check_extern_erf
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tools)
set_tests_properties(drift_check_selftest PROPERTIES LABELS unit
                     ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
