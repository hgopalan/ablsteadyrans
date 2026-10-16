Building
========

Requirements
------------

- CMake 3.25 or newer and a C++20 compiler (GCC, Clang, AppleClang; nvcc for CUDA).
- MPI (Open MPI, MPICH or Cray MPICH).
- FFTW, for AMReX's FFT module (cuFFT is used in CUDA builds).
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
     - ``NONE`` or ``CUDA``. For an H100 add ``-DAMReX_CUDA_ARCH=9.0``.
   * - ``ABLSTEADY_FFT``
     - ``ON``
     - AMReX FFT, used by the Poisson layer from sub-phase 1.10.
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
