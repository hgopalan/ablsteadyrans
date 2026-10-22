# ---------------------------------------------------------------------------
# GUST -- top-level CMake build
#
#   git submodule update --init --recursive
#   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-rel -j 8
#   ctest --test-dir build-rel --output-on-failure
#
# AMReX comes from the submodule in external/amrex, pinned to a fixed SHA.
# To use an installed AMReX instead:
#
#   cmake -S . -B build-rel -DGUST_USE_INTERNAL_AMREX=OFF -DAMReX_ROOT=/path
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.25)   # required by the pinned AMReX

project(gust
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

option(GUST_MPI                "Build with MPI"                                   ON)
option(GUST_OMP                "Build with OpenMP"                                OFF)
option(GUST_FFT                "Build AMReX FFT (required by the Poisson layer)"   ON)
option(GUST_HYPRE              "Build AMReX with hypre (optional preconditioner)" OFF)
option(GUST_TESTING            "AMReX_TESTING + assertions (SNaN fills, as GPU CI)" OFF)
option(GUST_USE_INTERNAL_AMREX "Use the AMReX git submodule"                      ON)
option(GUST_ENABLE_TESTS       "Build the unit tests and register all tests"      ON)

set(GUST_PRECISION DOUBLE CACHE STRING "Floating-point precision: DOUBLE or SINGLE")
set_property(CACHE GUST_PRECISION PROPERTY STRINGS DOUBLE SINGLE)

# GPU backend, passed to AMReX. Kestrel's GPU nodes are H100 (CUDA arch 9.0).
# HIP (AMD) and SYCL (Intel) need the vendor compiler as CMAKE_CXX_COMPILER
# (ROCm clang++, icpx) and, for HIP, AMReX_AMD_ARCH. The FFT library follows
# the backend (AMReX FFT): FFTW on the CPU, cuFFT, rocFFT or oneMKL.
set(GUST_GPU_BACKEND NONE CACHE STRING "GPU backend: NONE, CUDA, HIP or SYCL")
set_property(CACHE GUST_GPU_BACKEND PROPERTY STRINGS NONE CUDA HIP SYCL)

# Extra flags for mpiexec in the tests (CI passes --oversubscribe to Open MPI
# so 4-rank tests run on 2-core runners).
set(GUST_MPIEXEC_PREFLAGS "" CACHE STRING "Extra mpiexec flags for the tests")

# The solver was called ablsteady until October 2026. The old cache options
# (-DABLSTEADY_MPI=OFF and so on) still work: each sets its GUST_* option, with
# a warning, and is then removed from the cache. If both names are given, the
# old one's value is used. (message(DEPRECATION) would be silent on a second
# configure: AMReX sets CMAKE_WARN_DEPRECATED=OFF in the cache.)
foreach(opt MPI OMP FFT HYPRE TESTING USE_INTERNAL_AMREX ENABLE_TESTS
            PRECISION GPU_BACKEND MPIEXEC_PREFLAGS)
    if(DEFINED ABLSTEADY_${opt})
        message(WARNING "ABLSTEADY_${opt} is deprecated; use GUST_${opt}")
        set_property(CACHE GUST_${opt} PROPERTY VALUE "${ABLSTEADY_${opt}}")
        unset(ABLSTEADY_${opt} CACHE)
    endif()
endforeach()

# The Poisson layer uses AMReX's FFT module for its
# preconditioner; stop before AMReX is configured without it.
if(NOT GUST_FFT)
    message(FATAL_ERROR "GUST needs GUST_FFT=ON: the Poisson layer uses AMReX's FFT module")
endif()

foreach(opt PRECISION GPU_BACKEND)
    get_property(allowed CACHE GUST_${opt} PROPERTY STRINGS)
    if(NOT GUST_${opt} IN_LIST allowed)
        message(FATAL_ERROR "GUST_${opt}=${GUST_${opt}} is not one of: ${allowed}")
    endif()
endforeach()

if(GUST_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
endif()

# --- AMReX ----------------------------------------------------------------
if(GUST_USE_INTERNAL_AMREX)
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/amrex/CMakeLists.txt)
        message(FATAL_ERROR
            "The AMReX submodule is empty. Run:\n"
            "  git submodule update --init --recursive")
    endif()

    set(AMReX_SPACEDIM       3)
    set(AMReX_PRECISION      ${GUST_PRECISION})
    set(AMReX_MPI            ${GUST_MPI})
    set(AMReX_OMP            ${GUST_OMP})
    set(AMReX_FFT            ${GUST_FFT})
    set(AMReX_HYPRE          ${GUST_HYPRE})
    set(AMReX_GPU_BACKEND    ${GUST_GPU_BACKEND})
    set(AMReX_FORTRAN        OFF)
    set(AMReX_INSTALL        OFF)
    set(AMReX_LINEAR_SOLVERS ON)
    set(AMReX_PARTICLES      OFF)
    set(AMReX_AMRLEVEL       OFF)
    set(AMReX_EB             OFF)
    # fcompare is the comparison tool for the regression tests.
    set(AMReX_PLOTFILE_TOOLS ON)
    if(GUST_TESTING)
        set(AMReX_TESTING    ON)
        set(AMReX_ASSERTIONS ON)
    endif()

    add_subdirectory(external/amrex)
    # For AMReXTargetHelpers (setup_target_for_cuda_compilation).
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/amrex/Tools/CMake)
else()
    find_package(AMReX REQUIRED 3D ${GUST_PRECISION} FFT)
endif()

# MSVC: M_PI and friends need _USE_MATH_DEFINES; NOMINMAX keeps <windows.h>
# (pulled in by MS-MPI and the CUDA runtime) from defining min and max macros.
if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
endif()

add_subdirectory(Source)

# --- Tests ----------------------------------------------------------------
if(GUST_ENABLE_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_subdirectory(Tests)
endif()
