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
)
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Unit)
target_link_libraries(unit_tests PRIVATE ablsteady_core GTest::gtest)
# Gold files (Tests/Gold) are read from the source tree.
target_compile_definitions(unit_tests PRIVATE ABLSTEADY_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# RT-1.10.1 projects a fixed wind with the Poisson layer: a test-only tool,
# since the application has no Poisson solve yet.
add_executable(poisson_projection Tools/poisson_projection.cpp)
target_link_libraries(poisson_projection PRIVATE ablsteady_core)
if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    setup_target_for_cuda_compilation(unit_tests)
    setup_target_for_cuda_compilation(poisson_projection)
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
                          "steady\\.residual_tol.*not an input of namespace 'steady'"
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
