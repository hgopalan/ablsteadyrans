Continuous integration
======================

Every pull request runs two workflows.

``.github/workflows/ci.yml``
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Job
     - What it does
   * - CPU (GCC Release, GCC Debug with ``AMReX_TESTING``, Clang Release)
     - Builds with Open MPI, runs ``ctest -L "unit|smoke"`` (4-rank tests on 2 cores with ``--oversubscribe``),
       and checks that the tests left the working tree clean. Linux GCC traps floating-point exceptions.
   * - Single precision
     - Compiles only. The steady stop needs double precision.
   * - CUDA sm_90
     - Compiles only, for the H100. Hosted runners have no GPU.
   * - Windows MSVC Release
     - Builds with MSVC, no MPI and FFTW from vcpkg, and runs ``ctest -L "unit|smoke"``. MSVC does not trap
       floating-point exceptions; that check stays with Linux GCC.
   * - Windows CUDA sm_90
     - Compiles only, with ``cl`` and ``nvcc`` (CUDA 12.6, Windows SDK 10.0.22621.0).
   * - HIP gfx90a
     - Compiles only, with ROCm 6.2.4 and rocFFT, for AMD GPUs (MI250X).
   * - SYCL
     - Compiles only, with Intel ``icpx`` and oneMKL, for Intel GPUs.
   * - ``extern_erf`` drift check
     - Runs the drift check's self-tests and the check itself against a clone of ERF, and lists upstream ERF
       changes to copied files (:doc:`extraction`).
   * - Spelling and whitespace
     - ``codespell`` (settings in ``.codespellrc``) and no trailing whitespace.

``.github/workflows/docs.yml``
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Job
     - What it does
   * - Build documentation
     - ``make -C docs html`` with warnings as errors, and uploads the HTML as a workflow artifact.
   * - Link check
     - ``make -C docs linkcheck``. It reaches third-party sites, so a failure can come from outside this
       repository.

``.github/workflows/static.yml``
--------------------------------

On every push to ``main`` (or by hand), the ``deploy`` job builds the documentation as above, with warnings as
errors, and publishes ``docs/_build/html`` to GitHub Pages. Only the built HTML is uploaded, not the repository.

Kestrel
-------

The regression tests labelled ``physics`` (``ctest -L physics``) and the performance tests will run nightly on
Kestrel CPU nodes and H100s once they exist (from sub-phase 1.1).
