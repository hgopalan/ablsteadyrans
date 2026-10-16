# ---------------------------------------------------------------------------
# ablsteady -- top-level CMake build
#
#   git submodule update --init --recursive
#   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-rel -j 8
#   ctest --test-dir build-rel --output-on-failure
#
# AMReX comes from the submodule in external/amrex, pinned to the SHA the
# development plan was checked against. To use an installed AMReX instead:
#
#   cmake -S . -B build-rel -DABLSTEADY_USE_INTERNAL_AMREX=OFF -DAMReX_ROOT=/path
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.25)   # required by the pinned AMReX

project(ablsteady
        VERSION 0.1.0
        DESCRIPTION "Stand-alone AMReX steady RANS solver for the atmospheric boundary layer"
        LANGUAGES C CXX)

# AMReX at the pinned SHA uses C++20 in its headers.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(ABLSTEADY_MPI                "Build with MPI"                                   ON)
option(ABLSTEADY_OMP                "Build with OpenMP"                                OFF)
option(ABLSTEADY_FFT                "Build AMReX FFT (needed by the Poisson layer)"     ON)
option(ABLSTEADY_HYPRE              "Build AMReX with hypre (optional preconditioner)" OFF)
option(ABLSTEADY_TESTING            "AMReX_TESTING + assertions (SNaN fills, as GPU CI)" OFF)
option(ABLSTEADY_USE_INTERNAL_AMREX "Use the AMReX git submodule"                      ON)
option(ABLSTEADY_ENABLE_TESTS       "Build the unit tests and register all tests"      ON)

set(ABLSTEADY_PRECISION DOUBLE CACHE STRING "Floating-point precision: DOUBLE or SINGLE")
set_property(CACHE ABLSTEADY_PRECISION PROPERTY STRINGS DOUBLE SINGLE)

# GPU backend, passed to AMReX. Kestrel's GPU nodes are H100 (CUDA arch 9.0).
set(ABLSTEADY_GPU_BACKEND NONE CACHE STRING "GPU backend: NONE or CUDA")
set_property(CACHE ABLSTEADY_GPU_BACKEND PROPERTY STRINGS NONE CUDA)

# Extra flags for mpiexec in the tests (CI passes --oversubscribe to Open MPI
# so 4-rank tests run on 2-core runners).
set(ABLSTEADY_MPIEXEC_PREFLAGS "" CACHE STRING "Extra mpiexec flags for the tests")

foreach(opt PRECISION GPU_BACKEND)
    get_property(allowed CACHE ABLSTEADY_${opt} PROPERTY STRINGS)
    if(NOT ABLSTEADY_${opt} IN_LIST allowed)
        message(FATAL_ERROR "ABLSTEADY_${opt}=${ABLSTEADY_${opt}} is not one of: ${allowed}")
    endif()
endforeach()

if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
endif()

# --- AMReX ----------------------------------------------------------------
if(ABLSTEADY_USE_INTERNAL_AMREX)
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/amrex/CMakeLists.txt)
        message(FATAL_ERROR
            "The AMReX submodule is empty. Run:\n"
            "  git submodule update --init --recursive")
    endif()

    set(AMReX_SPACEDIM       3)
    set(AMReX_PRECISION      ${ABLSTEADY_PRECISION})
    set(AMReX_MPI            ${ABLSTEADY_MPI})
    set(AMReX_OMP            ${ABLSTEADY_OMP})
    set(AMReX_FFT            ${ABLSTEADY_FFT})
    set(AMReX_HYPRE          ${ABLSTEADY_HYPRE})
    set(AMReX_GPU_BACKEND    ${ABLSTEADY_GPU_BACKEND})
    set(AMReX_FORTRAN        OFF)
    set(AMReX_INSTALL        OFF)
    set(AMReX_LINEAR_SOLVERS ON)
    set(AMReX_PARTICLES      OFF)
    set(AMReX_AMRLEVEL       OFF)
    set(AMReX_EB             OFF)
    # fcompare is the comparison tool for the regression tests.
    set(AMReX_PLOTFILE_TOOLS ON)
    if(ABLSTEADY_TESTING)
        set(AMReX_TESTING    ON)
        set(AMReX_ASSERTIONS ON)
    endif()

    add_subdirectory(external/amrex)
    # For AMReXTargetHelpers (setup_target_for_cuda_compilation).
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/amrex/Tools/CMake)
else()
    find_package(AMReX REQUIRED 3D ${ABLSTEADY_PRECISION})
endif()

add_subdirectory(Source)

# --- Tests ----------------------------------------------------------------
if(ABLSTEADY_ENABLE_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_subdirectory(Tests)
endif()
