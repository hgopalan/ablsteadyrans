ablsteady
=========

**ablsteady** computes steady-state wind fields in the atmospheric boundary layer (ABL) over complex terrain. It
solves the steady Reynolds-averaged Navier–Stokes (RANS) equations with a one-equation k–L turbulence model on a
terrain-following mesh. It is a stand-alone `AMReX <https://github.com/AMReX-Codes/amrex>`__ application and runs on
CPUs and GPUs.

A transient ABL code such as `ERF <https://github.com/erf-model/ERF>`__ reaches a steady wind field by integrating in
time for many hours of simulated flow. ablsteady solves for the steady state directly. It drives the residual of the
steady equations to zero with a pseudo-time iteration, and accelerates the iteration with line solvers, a good
initial field and coarse grids.

ablsteady reuses ERF components where they fit:

- the Axell–Liungman k–L turbulence model;
- the stretched terrain-following mesh with terrain smoothing;
- the terrain Poisson solver;
- the Monin–Obukhov surface layer.

The initial field comes from a mass-consistent model in the style of WindNinja, which follows
`FastWindTerrain <https://github.com/hgopalan/FastWindTerrain>`__.

Development status
------------------

The solver is developed in four phases. Each phase adds a faster solution method.

.. list-table::
   :header-rows: 1
   :widths: 10 70 20

   * - Phase
     - Solution method
     - Status
   * - 1
     - Pseudo-time SIMPLEC with implicit vertical lines (z-line Jacobi)
     - Complete
   * - 2
     - Line LU-SGS sweeps, a horizontal grid hierarchy and grid sequencing
     - In progress
   * - 3
     - Jacobian-free Newton–Krylov, preconditioned by SIMPLEC
     - In progress
   * - 4
     - Nonlinear FAS multigrid, GPU tuning, validation and release
     - Planned

The production solver today is the SIMPLEC stage (:doc:`solution_algorithm`). The target machine is NLR Kestrel:
CPU nodes with 104 Sapphire Rapids cores and GPU nodes with four NVIDIA H100.

Where to start
--------------

- **New users:** build the code (:doc:`building`), then work through the tutorials, starting with
  :doc:`tutorial_column`.
- **Setting up your own case:** :doc:`howto` and the :doc:`inputs`.
- **The model and the numerics:** the Mathematical model and Numerical method sections below.
- **Accuracy and speed:** :doc:`results`.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   building
   running
   output

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorial_column
   tutorial_flat_abl
   tutorial_hill
   tutorial_askervein

.. toctree::
   :maxdepth: 2
   :caption: How-to guides and reference

   howto
   inputs

.. toctree::
   :maxdepth: 2
   :caption: Mathematical model

   governing_equations
   turbulence_model
   terrain_coordinates
   boundary_conditions

.. toctree::
   :maxdepth: 2
   :caption: Numerical method

   discretization
   solution_algorithm
   initial_conditions
   grid_hierarchy

.. toctree::
   :maxdepth: 2
   :caption: Verification and validation

   test_cases
   results
   references
   bibliography

.. toctree::
   :maxdepth: 2
   :caption: Development

   testing
   ci
   extraction
   development
