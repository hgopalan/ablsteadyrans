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
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest
                 ${CMAKE_BINARY_DIR}/external/googletest EXCLUDE_FROM_ALL)

add_executable(unit_tests
    Unit/main.cpp
    Unit/test_Inputs.cpp
    Unit/test_IO.cpp
    Unit/test_Reductions.cpp
    Unit/test_SteadyCriterion.cpp
    Unit/test_Sha256.cpp
    Unit/test_Reference.cpp
    Unit/test_Mesh.cpp
    Unit/test_Poisson.cpp
    Unit/test_PoissonB.cpp
    Unit/test_WallDistance.cpp
    Unit/test_Convection.cpp
    Unit/test_Diffusion.cpp
    Unit/test_ImplicitOperator.cpp
    Unit/test_MassFlux.cpp
    Unit/test_Residual.cpp
    Unit/test_SurfaceLayer.cpp
    Unit/test_BoundaryConditions.cpp
    Unit/test_Closure.cpp
    Unit/test_Initializer.cpp
    Unit/test_Driver.cpp
    Unit/test_Predictor.cpp
    Unit/test_Simplec.cpp
    Unit/test_Newton.cpp
    Unit/test_Terrain.cpp
)
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Unit)
target_link_libraries(unit_tests PRIVATE ablsteady_core GTest::gtest)
# Gold files (Tests/Gold) are read from the source tree.
target_compile_definitions(unit_tests PRIVATE ABLSTEADY_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# RT-1.10.1 projects a fixed wind with the Poisson layer: a test-only tool,
# since the application has no Poisson solve yet.
add_executable(poisson_projection Tools/poisson_projection.cpp)
target_link_libraries(poisson_projection PRIVATE ablsteady_core)
# RT-1.11.2 computes the Tucker wall distance: another test-only tool.
add_executable(wall_distance Tools/wall_distance.cpp)
target_link_libraries(wall_distance PRIVATE ablsteady_core)
add_executable(boundary_fluxes Tools/boundary_fluxes.cpp)
target_link_libraries(boundary_fluxes PRIVATE ablsteady_core)
# RT-1.9.1: the full residual of the equilibrium log-law column with the closure.
add_executable(closure_column Tools/closure_column.cpp)
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

set(DECKS ${CMAKE_CURRENT_SOURCE_DIR}/Decks)

# fcompare is built with the AMReX submodule; with an installed AMReX, look
# for its amrex_fcompare executable.
if(TARGET fcompare)
    set(FCOMPARE $<TARGET_FILE:fcompare>)
else()
    find_program(FCOMPARE amrex_fcompare REQUIRED
                 HINTS ${AMReX_ROOT}/bin ${AMReX_DIR}/../../../bin)
endif()
set(EXPECT_ABORT ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/expect_abort.py)

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
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_smoke
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
# The smoke state is mostly zeros, so the plotfile comparison below cannot
# tell a restart from a fresh start; the log line can.
ablsteady_add_test(NAME smoke_restart NP 3 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke FIXTURES_SETUP smoke_restart
                   PREFIX ${EXPECT_ABORT} --success
                          "restarted from smoke_chk00000 at step 0" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_restart
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
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_unknown_key
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_range NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.n_cell must be three integers >= 1"
                          "mesh\\.is_periodic: z cannot be periodic"
                          "mesh\\.prob_lo and mesh\\.prob_hi must be finite" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_range
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_no_stage NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "forcing\\.gravity has no effect without steady\\.stage = simplec"
                          "convection\\.limiter has no effect without steady\\.stage = simplec"
                          "simplec\\.alpha_p has no effect without steady\\.stage = simplec"
                          "turb\\.Cmu0 has no effect without steady\\.stage = simplec"
                          "steady\\.cfl0 has no effect without steady\\.stage = simplec"
                          "io\\.turb_plot_file has no effect without steady\\.stage = simplec" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_no_stage
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_stage_no_guess NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "steady\\.stage = simplec needs a first guess" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_stage_no_guess
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_repeated_key NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.n_cell is given more than once"
                          "io\\.write_plotfile takes exactly one value"
                          "io\\.restart takes exactly one value" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_repeated_key
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_restart_mesh NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke
                   PREFIX ${EXPECT_ABORT}
                          "io\\.restart = smoke_chk00000: checkpoint mesh\\.n_cell" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_restart_mesh
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_restart_overwrite NP 1 LABELS "regression;smoke" RUN_DIR ${SMOKE_DIR}
                   FIXTURES_REQUIRED smoke
                   PREFIX ${EXPECT_ABORT}
                          "io\\.chk_file = smoke_chk would overwrite io\\.restart = smoke_chk00000" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_restart_overwrite
                           ${ABLSTEADY_TEST_AMREX_FLAGS})

# --- Mesh (sub-phase 1.4) -----------------------------------------------------
set(CHECK_MESH ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mesh.py)
set(GOLD ${CMAKE_CURRENT_SOURCE_DIR}/Gold)
# The mesh golds are double-precision values: a single-precision build
# compares the plotfile at a single-precision tolerance and skips the
# statistics checks, whose reader (ERF's extern_erf/erf_plotfile.py) handles
# double-precision plotfiles only.
if(ABLSTEADY_PRECISION STREQUAL "SINGLE")
    set(FCOMPARE_RTOL 1e-5)
else()
    set(FCOMPARE_RTOL 1e-12)
    # RT-1.4.1 and RT-1.4.2: ERF's Askervein mesh, coarsened; field
    # statistics of z and detJ and the quality report against
    # Gold/askervein_mesh.json.
    ablsteady_add_test(NAME rt_askervein_mesh NP 2 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${GOLD}/askervein_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_mesh_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.4.2: the quality report of a steep ridge, with the STF warning.
    ablsteady_add_test(NAME rt_steep_ridge_quality NP 1 LABELS "regression;smoke"
                       PREFIX ${CHECK_MESH} ${GOLD}/steep_ridge_mesh.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_mesh_steep_ridge
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()
# RT-1.4.3: a run without a solve writes the mapped-grid plotfile of the 2-D
# hill; it matches Gold/hill2d_plt00000 (relative tolerance 1e-12: compilers
# that fuse multiply-adds differ from those that do not in the last bit).
set(HILL2D_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_hill2d)
ablsteady_add_test(NAME rt_hill2d_run NP 2 LABELS "regression;smoke" RUN_DIR ${HILL2D_DIR}
                   FIXTURES_SETUP hill2d
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_mesh_hill2d ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_hill2d_matches NP 1 LABELS "regression;smoke" RUN_DIR ${HILL2D_DIR}
                   FIXTURES_REQUIRED hill2d
                   COMMAND ${FCOMPARE} --rel_tol ${FCOMPARE_RTOL} hill2d_plt00000 ${GOLD}/hill2d_plt00000)

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
                   COMMAND $<TARGET_FILE:poisson_projection> ${DECKS}/inputs_poisson_askervein
                           terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                           solver.poisson_reltol=0 ${ABLSTEADY_TEST_AMREX_FLAGS})
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(POISSON_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_poisson_askervein)
    ablsteady_add_test(NAME rt_poisson_askervein NP 2 LABELS "regression;smoke" RUN_DIR ${POISSON_DIR}
                       FIXTURES_SETUP poisson_askervein
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_poisson.py
                              ${GOLD}/poisson_askervein.json --
                       COMMAND $<TARGET_FILE:poisson_projection> ${DECKS}/inputs_poisson_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # The potential converges to solver.poisson_reltol = 1e-10; 1e-8 allows
    # for compilers that round differently (fused multiply-adds).
    ablsteady_add_test(NAME rt_poisson_askervein_matches NP 1 LABELS "regression;smoke" RUN_DIR ${POISSON_DIR}
                       FIXTURES_REQUIRED poisson_askervein
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 poisson_askervein_plt00000
                               ${GOLD}/poisson_askervein_plt00000)
    # RT-1.11.1: the same projection with back end B (MLMG). It solves the same
    # discrete operator, so its potential matches back end A's gold.
    set(POISSON_B_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_poisson_askervein_mlmg)
    ablsteady_add_test(NAME rt_poisson_askervein_mlmg NP 2 LABELS "regression;smoke" RUN_DIR ${POISSON_B_DIR}
                       FIXTURES_SETUP poisson_askervein_mlmg
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_poisson.py
                              ${GOLD}/poisson_askervein_mlmg.json --
                       COMMAND $<TARGET_FILE:poisson_projection> ${DECKS}/inputs_poisson_askervein_mlmg
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_poisson_askervein_mlmg_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${POISSON_B_DIR} FIXTURES_REQUIRED poisson_askervein_mlmg
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 poisson_askervein_mlmg_plt00000
                               ${GOLD}/poisson_askervein_plt00000)
    # RT-1.12.1: the mass-consistent initial field on the Askervein mesh (log-law
    # first guess, alpha_v = 0.5): iterations, divergence and field extrema
    # against Gold/init_askervein.json, the state against
    # Gold/init_askervein_plt00000 (1e-8 as the Poisson golds). The gold is a
    # double-precision run, so a single-precision build skips it.
    set(INIT_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_init_askervein)
    ablsteady_add_test(NAME rt_init_askervein NP 2 LABELS "regression;smoke" RUN_DIR ${INIT_DIR}
                       FIXTURES_SETUP init_askervein
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_init.py askervein
                              ${GOLD}/init_askervein.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_init_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_init_askervein_matches NP 1 LABELS "regression;smoke" RUN_DIR ${INIT_DIR}
                       FIXTURES_REQUIRED init_askervein
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 init_askervein_plt00000
                               ${GOLD}/init_askervein_plt00000)
    # RT-1.12.2: the hilltop speed-up over FastWindTerrain's bump for alpha_v =
    # 1, 0.7, 0.5, 0.3, 0.1 (five runs) against FastWindTerrain's own on the
    # same hill and profile (Gold/init_bump_trend.json: ranking, sign, spread).
    ablsteady_add_test(NAME rt_init_bump_trend NP 2 LABELS "regression;smoke"
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_init.py trend
                              ${GOLD}/init_bump_trend.json --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_init_bump
                               terrain.file=${DECKS}/terrain_bump_fwt.txt ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.11.2: the Tucker wall distance on ERF's Neutral_Hill_2D mesh: ERF's
    # wall-distance checks (Tests/check_walldist.py), and the field against
    # Gold/walldist_hill2d_plt00000. ERF's reader takes double-precision
    # plotfiles only.
    set(WALLDIST_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_walldist_hill2d)
    ablsteady_add_test(NAME rt_walldist_hill2d NP 2 LABELS "regression;smoke" RUN_DIR ${WALLDIST_DIR}
                       FIXTURES_SETUP walldist_hill2d
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_walldist.py
                              --hmax 100 --half-width 500 --
                       COMMAND $<TARGET_FILE:wall_distance> ${DECKS}/inputs_walldist_hill2d
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_walldist_hill2d_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${WALLDIST_DIR} FIXTURES_REQUIRED walldist_hill2d
                       COMMAND ${FCOMPARE} --rel_tol 1e-8 walldist_hill2d_plt00000
                               ${GOLD}/walldist_hill2d_plt00000)
endif()

# --- RT-1.8.1: boundary fluxes of a fixed log-law state (sub-phase 1.8) ------
# Double precision only (the gold is double).
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(BC_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_bc_neutral_flat)
    ablsteady_add_test(NAME rt_bc_neutral_flat NP 2 LABELS "regression;smoke" RUN_DIR ${BC_DIR}
                       FIXTURES_SETUP bc_neutral_flat
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_bc.py
                              --ustar 0.35 --z0 0.1 --direction 15 --dz 15.625 --
                       COMMAND $<TARGET_FILE:boundary_fluxes> ${DECKS}/inputs_bc_neutral_flat
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_bc_neutral_flat_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${BC_DIR} FIXTURES_REQUIRED bc_neutral_flat
                       COMMAND ${FCOMPARE} --rel_tol 1e-10 bc_neutral_flat_plt00000
                               ${GOLD}/bc_neutral_flat_plt00000)
endif()

# --- RT-1.9.1: the neutral log-law column with the k-L closure (sub-phase 1.9)
# Double precision only (the gold is double).
if(ABLSTEADY_PRECISION STREQUAL "DOUBLE")
    set(CL_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_closure_column)
    ablsteady_add_test(NAME rt_closure_column NP 2 LABELS "regression;smoke" RUN_DIR ${CL_DIR}
                       FIXTURES_SETUP closure_column
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_closure.py
                              --ustar 0.35 --z0 0.1 --direction 15 --dz 15.625 --
                       COMMAND $<TARGET_FILE:closure_column> ${DECKS}/inputs_closure_column
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    ablsteady_add_test(NAME rt_closure_column_matches NP 1 LABELS "regression;smoke"
                       RUN_DIR ${CL_DIR} FIXTURES_REQUIRED closure_column
                       COMMAND ${FCOMPARE} --rel_tol 1e-10 closure_column_plt00000
                               ${GOLD}/closure_column_plt00000)
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
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_column.py
                              --ustar 0.35 --z0 0.1 --force 1.225e-4 --height 1000 --dz 15.625 --orders 6 --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_column ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()
# A restart continues the stage bit for bit: 40 iterations, continued to 80
# from the checkpoint (the driver and closure histories), against 80 straight.
set(SCR_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_steady_column_restart)
ablsteady_add_test(NAME rt_steady_column_restart_first NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}
                   FIXTURES_SETUP steady_column_first
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_column_restart steady.max_iterations=40
                           init.first_guess=profile ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_steady_column_restart_continued NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}
                   FIXTURES_REQUIRED steady_column_first FIXTURES_SETUP steady_column_continued
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_column_restart steady.max_iterations=80
                           io.restart=steady_column_restart_chk00040 ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME rt_steady_column_restart_straight NP 2 LABELS "regression;smoke" RUN_DIR ${SCR_DIR}/straight
                   FIXTURES_SETUP steady_column_straight
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_column_restart steady.max_iterations=80
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
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_neutral.py --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_neutral_flat ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.16.1: ERF's Neutral_Hill_2D to steady (about 80 s on two ranks);
    # the verdict, a two-order drop of the u residual and ERF's own
    # check_hill.py --physics (extern_erf/) on the turbulence plotfile.
    ablsteady_add_test(NAME rt_steady_hill2d NP 2 LABELS "regression"
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_hill.py
                              --erf-check check_hill.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_hill2d ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-2.1.1: the same case with the line LU-SGS sweeps (two symmetric
    # sweeps, CFL cap 16): the verdict, the residual drop and ERF's checks.
    ablsteady_add_test(NAME rt_steady_hill2d_lusgs NP 2 LABELS "regression"
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_hill.py
                              --erf-check check_hill.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_hill2d_lusgs ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.16.2: ERF's Neutral_Hill_3D (64 x 64 x 20) to steady with
    # check_hill3d.py; minutes on four ranks, so label nightly (not in
    # ctest -L regression).
    ablsteady_add_test(NAME rt_steady_hill3d NP 4 LABELS "nightly"
                       PREFIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_hill.py
                              --erf-check check_hill3d.py --orders 2 --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_hill3d ${ABLSTEADY_TEST_AMREX_FLAGS})
    # RT-1.16.3 (the convergence half): Askervein on the ERF grid with an
    # inflow and an outflow face to the steady verdict (about 2 s on two
    # ranks). The comparison with ERF's URANS on lines A and AA needs the
    # ERF reference run and is not in this version.
    ablsteady_add_test(NAME rt_steady_askervein NP 2 LABELS "regression;smoke"
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_askervein
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    # Beyond the plan: Askervein in a conventionally neutral boundary layer.
    # The precursor column runs to steady, then tools/make_inflow_table.py
    # checks its 10 m wind is along +x (0.01 deg), u* = 0.6 m/s (1e-3) and the
    # committed inflow table matches it (1e-3). The 3-D case runs on that
    # table to the steady verdict.
    set(CNBL_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/rt_precursor_askervein_cnbl)
    ablsteady_add_test(NAME rt_precursor_askervein_cnbl NP 2 LABELS "regression" RUN_DIR ${CNBL_DIR}
                       FIXTURES_SETUP precursor_askervein_cnbl
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_precursor_askervein_cnbl
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
    add_test(NAME rt_precursor_askervein_cnbl_table
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/make_inflow_table.py
                     precursor_askervein_cnbl_plt precursor_askervein_cnbl_turb
                     --ug 17.028798 --vg -7.380863 --check ${DECKS}/askervein_cnbl_inflow.txt
             WORKING_DIRECTORY ${CNBL_DIR})
    set_tests_properties(rt_precursor_askervein_cnbl_table PROPERTIES LABELS "regression"
                         FIXTURES_REQUIRED precursor_askervein_cnbl)
    ablsteady_add_test(NAME rt_steady_askervein_cnbl NP 2 LABELS "regression"
                       PREFIX ${EXPECT_ABORT} --success "ablsteady: steady after" --
                       COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_steady_askervein_cnbl
                               terrain.file=${PROJECT_SOURCE_DIR}/extern_erf/askervein.txt
                               bc.inflow_file=${DECKS}/askervein_cnbl_inflow.txt
                               ${ABLSTEADY_TEST_AMREX_FLAGS})
endif()

# --- Reference archives (sub-phase 1.1): ablsteady_refavg end to end --------
# Three copies of the smoke plotfile at t = 0, 5 and 10 s: their trapezoid
# average (weights 1/4, 1/2, 1/4, fma) reproduces the plotfile bit for bit.
# This checks the archive round trip; the weights themselves are tested in
# the unit tests (identical snapshots cannot tell weights apart).
set(REFAVG_FIX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/refavg_fixtures.py)
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
                           "refavg.inputs_files=${DECKS}/inputs_smoke ${DECKS}/inputs_restart"
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
# The Kestrel scripts cannot run here; at least they must parse.
file(GLOB KESTREL_SCRIPTS ${PROJECT_SOURCE_DIR}/tools/kestrel/*.sh)
foreach(script IN LISTS KESTREL_SCRIPTS)
    get_filename_component(sname ${script} NAME_WE)
    add_test(NAME kestrel_parse_${sname} COMMAND bash -n ${script})
    set_tests_properties(kestrel_parse_${sname} PROPERTIES LABELS unit)
endforeach()
# ... and a dry run with stub module/sbatch/srun exercises submit and archive.
add_test(NAME kestrel_dryrun
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/kestrel_dryrun.sh ${PROJECT_SOURCE_DIR}
                 $<TARGET_FILE:ablsteady_refavg> ${SMOKE_DIR}/smoke_plt00000
                 ${CMAKE_CURRENT_BINARY_DIR}/run/kestrel_dryrun)
set_tests_properties(kestrel_dryrun PROPERTIES LABELS "regression;smoke" FIXTURES_REQUIRED smoke)

# --- RT-1.3.1 helper: the extern_erf drift check's own tests -----------------
add_test(NAME drift_check_selftest
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_check_extern_erf
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tools)
set_tests_properties(drift_check_selftest PROPERTIES LABELS unit
                     ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
