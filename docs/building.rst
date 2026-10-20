Building
========

Requirements
------------

- CMake 3.25 or newer and a C++20 compiler (GCC, Clang, AppleClang, MSVC; nvcc for CUDA, ROCm ``clang++`` for HIP,
  ``icpx`` for SYCL).
- MPI (Open MPI, MPICH or Cray MPICH), unless ``ABLSTEADY_MPI=OFF``.
- An FFT library for AMReX's FFT module: FFTW on the CPU, cuFFT (CUDA), rocFFT (HIP) or oneMKL (SYCL).
- Python 3, for the test scripts.

AMReX and GoogleTest are git submodules. AMReX is pinned to the commit the development plan was checked against.

Quick start
-----------

.. code-block:: bash

   git clone --recursive https://github.com/hgopalan/ablsteadyrans.git
   cd ablsteadyrans
   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
   cmake --build build-rel -j 8
   ctest --test-dir build-rel --output-on-failure

The executable is ``build-rel/ablsteady``. On a case-insensitive file system, name the build directory something
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
   * - ``ABLSTEADY_MPI``
     - ``ON``
     - Build with MPI.
   * - ``ABLSTEADY_OMP``
     - ``OFF``
     - Build with OpenMP.
   * - ``ABLSTEADY_PRECISION``
     - ``DOUBLE``
     - ``DOUBLE`` or ``SINGLE``. The steady stop needs double precision; single precision is kept compiling only.
   * - ``ABLSTEADY_GPU_BACKEND``
     - ``NONE``
     - ``NONE``, ``CUDA``, ``HIP`` or ``SYCL``. For an H100 add ``-DCMAKE_CUDA_ARCHITECTURES=90``; HIP needs
       ``-DAMReX_AMD_ARCH`` (for example ``gfx90a``) and ROCm's ``clang++`` as the C++ compiler; SYCL needs ``icpx``.
   * - ``ABLSTEADY_FFT``
     - ``ON``
     - AMReX FFT. Required: the Poisson layer (sub-phase 1.10) uses its preconditioner. Configuring with ``OFF`` stops
       with an error, and an external AMReX must have been built with FFT (``find_package`` asks for it).
   * - ``ABLSTEADY_HYPRE``
     - ``OFF``
     - AMReX with hypre, for the optional Newton preconditioner (sub-phase 3.3).
   * - ``ABLSTEADY_TESTING``
     - ``OFF``
     - ``AMReX_TESTING`` and assertions on, as in GPU CI.
   * - ``ABLSTEADY_ENABLE_TESTS``
     - ``ON``
     - Build the unit tests and register all tests with CTest.
   * - ``ABLSTEADY_MPIEXEC_PREFLAGS``
     - empty
     - Extra ``mpiexec`` flags for the tests. CI passes ``--oversubscribe`` to Open MPI.
   * - ``ABLSTEADY_USE_INTERNAL_AMREX``
     - ``ON``
     - Use the AMReX submodule. ``OFF`` uses an installed AMReX found through ``AMReX_ROOT``.

Windows
-------

On Windows the solver builds with MSVC and no MPI. FFTW comes from vcpkg. From a "Developer PowerShell for VS"
(or after ``ilammy/msvc-dev-cmd`` in CI):

.. code-block:: powershell

   vcpkg install fftw3:x64-windows
   cmake -S . -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release `
     -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
     -DVCPKG_TARGET_TRIPLET=x64-windows -DABLSTEADY_MPI=OFF
   cmake --build build-rel
   $env:PATH = "C:/vcpkg/installed/x64-windows/bin;$env:PATH"
   ctest --test-dir build-rel -L "unit|smoke" --output-on-failure

The FFTW DLLs are found through ``PATH``. The Kestrel shell-script tests are not registered on Windows. A CUDA
build on Windows uses ``cl`` with ``nvcc``; with CUDA 12.6 and C++20 use Windows SDK 10.0.22621.0
(``-DCMAKE_SYSTEM_VERSION=10.0.22621.0``), as the CI job does.

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
