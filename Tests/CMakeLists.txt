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
)
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Unit)
target_link_libraries(unit_tests PRIVATE ablsteady_core GTest::gtest)
if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    setup_target_for_cuda_compilation(unit_tests)
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
    cmake_parse_arguments(T "" "NAME;NP;FIXTURES_SETUP;FIXTURES_REQUIRED;RUN_DIR"
                          "LABELS;COMMAND;PREFIX" ${ARGN})
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
        set_tests_properties(${T_NAME} PROPERTIES FIXTURES_REQUIRED ${T_FIXTURES_REQUIRED})
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
                          "turb\\.Cmu0.*not an input of namespace 'turb'"
                          "mehs\\.n_cell.*unknown namespace 'mehs'"
                          "steady\\.residual_tol.*not an input of namespace 'steady'" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_unknown_key
                           ${ABLSTEADY_TEST_AMREX_FLAGS})
ablsteady_add_test(NAME abort_range NP 1 LABELS "regression;smoke"
                   PREFIX ${EXPECT_ABORT}
                          "mesh\\.n_cell must be three integers >= 1"
                          "mesh\\.is_periodic: z cannot be periodic"
                          "mesh\\.prob_lo and mesh\\.prob_hi must be finite" --
                   COMMAND $<TARGET_FILE:ablsteady> ${DECKS}/inputs_abort_range
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

# --- RT-1.3.1 helper: the extern_erf drift check's own tests -----------------
add_test(NAME drift_check_selftest
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_check_extern_erf
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tools)
set_tests_properties(drift_check_selftest PROPERTIES LABELS unit
                     ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
