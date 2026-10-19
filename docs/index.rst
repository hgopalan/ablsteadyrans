ablsteady
=========

**ablsteady** is a stand-alone `AMReX <https://github.com/AMReX-Codes/amrex>`__ application for fast steady-state
RANS of the atmospheric boundary layer (ABL) over complex terrain. It reuses `ERF <https://github.com/erf-model/ERF>`__
components where they fit: the Axell–Liungman k–L closure, the stretched terrain-following mesh with terrain
smoothing, the terrain Poisson solver, the MOST surface layer and the column tridiagonal solver. The starting field
comes from a mass-consistent (WindNinja-style) initializer that builds on
`FastWindTerrain <https://github.com/hgopalan/FastWindTerrain>`__.

The steady solver is built in four phases. Each phase ends with a faster solver than the one before.

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Phase
     - Solver at the end of the phase
   * - 1
     - Stand-alone application and pseudo-time SIMPLEC with implicit z (z-line Jacobi)
   * - 2
     - Line LU-SGS sweeps, an x/y grid hierarchy and grid sequencing
   * - 3
     - Jacobian-free Newton–Krylov finish, preconditioned by one SIMPLEC iteration
   * - 4
     - Nonlinear FAS multigrid, GPU tuning, scaling, benchmark, validation and release

The target machine is NLR Kestrel: CPU nodes with 104 Sapphire Rapids cores, and GPU nodes with four NVIDIA H100.

**Status.** Phase 1 is under way. The application reads and checks its inputs, builds the stretched
terrain-following mesh and the solution fields, and writes and restarts from plotfiles and checkpoints. The definition
of a steady solution, the Poisson layer, the convection and diffusion operators, the residual, the surface layer, the
domain boundary conditions, the k–L closure and the mass-consistent initial field are in place, and the SIMPLEC
stage (:doc:`solver`) runs under the steady driver.

.. toctree::
   :maxdepth: 2
   :caption: Using ablsteady

   building
   running
   inputs
   output

.. toctree::
   :maxdepth: 2
   :caption: Model formulation

   governing_cartesian
   governing_curvilinear
   discretization
   boundary_conditions
   initial_conditions
   solver

.. toctree::
   :maxdepth: 2
   :caption: Verification

   cases
   results
   references

.. toctree::
   :maxdepth: 2
   :caption: Development

   testing
   ci
   extraction
   development
