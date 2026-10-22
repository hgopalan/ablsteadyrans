Installation
============

Requirements
------------

- CMake 3.25 or newer and a C++20 compiler (GCC, Clang, AppleClang, MSVC; nvcc for CUDA, ROCm ``clang++`` for HIP,
  ``icpx`` for SYCL).
- MPI (Open MPI, MPICH or Cray MPICH), unless ``GUST_MPI=OFF``.
- An FFT library for AMReX's FFT module: FFTW on the CPU, cuFFT (CUDA), rocFFT (HIP) or oneMKL (SYCL).
- Python 3, for the test scripts.

AMReX and GoogleTest are git submodules. AMReX is pinned to a fixed commit.

Quick start
-----------

.. code-block:: bash

   git clone --recursive https://github.com/hgopalan/gust-abl.git
   cd gust-abl
   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
   cmake --build build-rel -j 8
   ctest --test-dir build-rel --output-on-failure

The executable is ``build-rel/gust``. On a case-insensitive file system, name the build directory something
other than ``build``.

CMake options
-------------

An option outside its allowed values stops the configuration.

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Option
     - Default
     - Effect
   * - ``GUST_MPI``
     - ``ON``
     - Build with MPI.
   * - ``GUST_OMP``
     - ``OFF``
     - Build with OpenMP.
   * - ``GUST_PRECISION``
     - ``DOUBLE``
     - ``DOUBLE`` or ``SINGLE``. The steady stop needs double precision; single precision is kept compiling only.
   * - ``GUST_GPU_BACKEND``
     - ``NONE``
     - ``NONE``, ``CUDA``, ``HIP`` or ``SYCL``. For an H100 add ``-DCMAKE_CUDA_ARCHITECTURES=90``; HIP needs
       ``-DAMReX_AMD_ARCH`` (for example ``gfx90a``) and ROCm's ``clang++`` as the C++ compiler; SYCL needs ``icpx``.
   * - ``GUST_FFT``
     - ``ON``
     - AMReX FFT. Required: the Poisson solver uses its FFT preconditioner. Configuring with ``OFF`` stops
       with an error, and an external AMReX must have been built with FFT (``find_package`` asks for it).
   * - ``GUST_HYPRE``
     - ``OFF``
     - AMReX with hypre, for the optional Newton–Krylov preconditioner (in development).
   * - ``GUST_TESTING``
     - ``OFF``
     - ``AMReX_TESTING`` and assertions on, as in GPU CI.
   * - ``GUST_ENABLE_TESTS``
     - ``ON``
     - Build the unit tests and register all tests with CTest.
   * - ``GUST_MPIEXEC_PREFLAGS``
     - empty
     - Extra ``mpiexec`` flags for the tests. CI passes ``--oversubscribe`` to Open MPI.
   * - ``GUST_USE_INTERNAL_AMREX``
     - ``ON``
     - Use the AMReX submodule. ``OFF`` uses an installed AMReX found through ``AMReX_ROOT``.

The options were called ``ABLSTEADY_*`` before the rename to GUST (October 2026). The old names still work: each
sets its ``GUST_*`` option and prints a warning.

Windows
-------

On Windows the solver builds with MSVC and no MPI, and FFTW comes from vcpkg. It needs Visual Studio 2022 with
the "Desktop development with C++" workload (including "C++ CMake tools for Windows", which provides CMake and
Ninja), Git and Python 3. The vcpkg that ships with Visual Studio supports only manifest mode, so use a standalone
clone. From a "Developer PowerShell for VS":

.. code-block:: powershell

   git clone https://github.com/microsoft/vcpkg $HOME/vcpkg
   & $HOME/vcpkg/bootstrap-vcpkg.bat
   $env:VCPKG_ROOT = "$HOME/vcpkg"
   & $env:VCPKG_ROOT/vcpkg install fftw3:x64-windows

   git clone --recursive https://github.com/hgopalan/gust-abl.git
   cd gust-abl
   cmake -S . -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release `
     -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
     -DVCPKG_TARGET_TRIPLET=x64-windows -DGUST_MPI=OFF
   cmake --build build-rel
   $env:PATH = "$env:VCPKG_ROOT/installed/x64-windows/bin;$env:PATH"
   ctest --test-dir build-rel -L "unit|smoke" --output-on-failure

The FFTW DLLs are found through ``PATH``. Without MPI every test runs on one rank, and MSVC does not trap
floating-point exceptions (:doc:`testing`). The Kestrel shell-script tests are not registered on Windows.

A CUDA build on Windows uses ``cl`` with ``nvcc``. With CUDA 12.6 and C++20 it needs Windows SDK 10.0.22621.0:
install that SDK, open the prompt with ``-winsdk=10.0.22621.0`` (for example ``VsDevCmd.bat -arch=x64
-winsdk=10.0.22621.0``) and configure with ``-DCMAKE_SYSTEM_VERSION=10.0.22621.0``, as the CI job does.

Kestrel
-------

``tools/kestrel/build_kestrel.sh cpu`` (on a CPU login node) or ``tools/kestrel/build_kestrel.sh gpu`` (on a GPU
login node) builds into ``build-kestrel-<cpu|gpu>``. The modules come from ``tools/kestrel/kestrel_env.sh``, which
follows ERF's Kestrel profiles at the pinned commit (after Kestrel's RHEL9 upgrade, erf-model/ERF#4292). Kestrel's
module names change with system updates, so check them with ``module avail`` first. A CUDA build needs memory: on
a shared node, build inside an allocation with ``--mem=240G`` or ``--exclusive``.

The scripts have not yet been run on Kestrel. :doc:`references` describes the ERF reference runs.

Documentation
-------------

.. code-block:: bash

   python3 -m pip install -r docs/requirements.txt
   make -C docs html

Warnings are errors (``-W`` in ``docs/Makefile``), so the local build and the CI build are the same. The pages
are written to ``docs/_build/html``.
