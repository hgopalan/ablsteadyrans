Poisson layer
=============

Three parts of ablsteady solve a Poisson problem on the terrain-following mesh: the mass-consistent initializer, the
SIMPLEC pressure correction and the Tucker wall distance. They share one layer, ``Source/Poisson.{H,cpp}``:

.. code-block:: c++

   PoissonSolver solver(mesh, geom, z_phys_nd, read_poisson_options());
   PoissonResult r = solver.solve(rhs, phi, face_coeffs, bc, &fluxes);

Two back ends solve the same discrete operator (``solver.poisson_backend``):

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Back end
     - What it is
   * - ``gmres_fft`` (A, default)
     - ERF's own operator, ``TerrainPoisson``, copied unchanged into ``extern_erf/`` (:doc:`extraction`), solved by
       ``amrex::GMRES`` with ERF's FFT preconditioner (sub-phase 1.10).
   * - ``mlmg``, ``gmres_mlmg`` (B)
     - AMReX's ``MLTerrainPoisson`` with geometric multigrid, alone or as GMRES's preconditioner (one V-cycle per
       iteration), sub-phase 1.11. Its stencil is ERF's term by term: with the same coefficients the two back ends
       give the same :math:`\phi` to 1e-9 relative when both are solved to 1e-12 (UT-1.11.2). It takes no
       ``z_cross`` (no split kernel) and ignores ``FaceCoeffs::precond_ratio``; ``solver.poisson_precond`` with
       back end B, and ``solver.poisson_restart`` with ``mlmg``, stop the run as keys with no effect. Its convergence
       test is MLMG's (maximum norm of the residual against that of the right-hand side; ``max_iter`` counts
       V-cycles), or GMRES's under ``gmres_mlmg``. A one-cell periodic direction is MLMG's hidden direction. It
       copies the mesh metrics when built, so a mesh changed in place needs a new ``PoissonSolver``.

On the Askervein projection test (30 × 30 × 18) back end A takes 10 GMRES iterations and back end B 6 MLMG
iterations.

**Single precision.** Rounding the solution to float leaves a residual floor of about
:math:`\varepsilon N_z^2/2` relative to the right-hand side (independent of the problem's scale): 2e-4 on 64 layers,
1e-3 on 128. The default ``solver.poisson_reltol`` in single precision (1e-5) is below that floor on fine meshes;
pick about :math:`4\varepsilon N_z^2`.

**AMReX finding (not fixed here).** ``MLTerrainPoisson``'s fused CPU red-black smoother (AMReX
``AMReX_MLTerrainPoisson.cpp``, the loop ``sa(i,j,k) += ta(i,j,k)`` and its neighbours) is vectorized over a stride-2
loop: it adds whole vector registers and stores only the active colour. The inactive lanes of its scratch array are
never written, so they can raise FE_INVALID (with ``amrex.fpe_trap_invalid = 1`` the run traps, for boxes wide
enough to vectorize: 32+ cells in double precision, fewer in single), while the discarded lanes never reach the
result. Back end B therefore switches the invalid trap off around the MLMG solve (``amrex::disableFPExcept``, then
``std::feclearexcept`` and ``amrex::setFPExcept``) and checks the solution for non-finite values afterwards. The
back end B tests run under CI's full flags on 64-cell boxes. Zeroing the scratch array (or not vectorizing the
stride-2 loops) in AMReX would remove the need.

Operator
--------

.. math::

   A\phi = \frac{1}{J}\Bigl[\partial_\xi\bigl(a_x\,\beta_x\,f_x\bigr) + \partial_\eta\bigl(a_y\,\beta_y\,f_y\bigr)
           + \partial_\zeta\bigl(a_z\,\beta_z\,f_z\bigr)\Bigr],

where :math:`J, a_x, a_y, a_z` are the mesh metrics (:doc:`mesh`), :math:`f = -\,\mathrm{flux}(\phi)` are ERF's face
gradients of :math:`\phi` in physical space, terrain cross terms included, and :math:`\beta` are optional face
coefficients (finite and > 0; 1 when not given). The coefficients enter as premultiplied areas :math:`a\beta`, so the
copied stencil is unchanged, and :math:`\beta = 1` gives ERF's operator bit for bit.

**Anisotropic coefficients.** On the terrain-following mesh the x- and y-face gradients are the physical
:math:`\partial_x\phi, \partial_y\phi`, but the z-face gradient is
:math:`n\cdot\nabla\phi = \partial_z\phi - (h_x\partial_x\phi + h_y\partial_y\phi)`. The initializer's operator
:math:`\nabla\cdot(\mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2)\nabla\phi)` needs
:math:`\alpha_v^2\partial_z\phi - \alpha_h^2(h_x\partial_x\phi + h_y\partial_y\phi)` on the z-faces, which a single
z-face multiplier cannot give. ``FaceCoeffs::z_cross`` is the coefficient of that cross part: the z-face then carries
:math:`\beta_z\,\partial_z\phi - \beta_{zc}(h_x\partial_x\phi + h_y\partial_y\phi)`, added to ERF's operator as
:math:`(\beta_z - \beta_{zc})` times the cross part by a separate kernel (``Source/PoissonCross_K.H``, ERF's formula;
sub-phase 1.12 corrected one node index of its low x-edge metric, which the first transcription read off a y-edge, so
the split operator now has the mirror symmetries of the mesh).
Without it, the anisotropic test over a 60 m hill does not converge (order 0); with it, it is second order. Scalar
coefficients (SIMPLEC's :math:`1/a_P`) leave it unset.

``solve`` also returns the face fluxes :math:`F = \beta\,\mathrm{flux}(\phi)` when asked. They are consistent with
``face_divergence``:

.. math::

   A\phi = -\,\nabla_h\!\cdot F, \qquad
   \nabla_h\!\cdot F = \frac{1}{J}\bigl[\delta_\xi(a_x F_x) + \delta_\eta(a_y F_y) + \delta_\zeta(a_z F_z)\bigr],

so a projection solves :math:`A\phi = \nabla_h\!\cdot u` and sets :math:`u \leftarrow u + F`; the divergence left is
the residual of the solve.

Boundary conditions
-------------------

All boundary conditions are homogeneous (``PoissonBC``, one entry per domain face):

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Face
     - Meaning
   * - ``periodic``
     - Exactly the periodic directions of the mesh.
   * - ``neumann``
     - Zero total flux through the face, the terrain cross term included (ERF#4272, merged as ERF#4279).
   * - ``dirichlet``
     - :math:`\phi = 0` on the face (an outflow or open face; the ground for the wall distance).

**Dirichlet ground.** ERF's kernel has no ground flux (its ground is a wall). With a Dirichlet ground, back end A
fills the ghost cells below the ground with :math:`-\phi` (so the bottom row's cross terms see :math:`\phi = 0` on the
ground), applies ERF's kernel, and adds the ground-face flux from ERF's own one-sided Dirichlet formula,
``terrpoisson_flux_zlo_dir`` (in the copied kernel file). The FFT preconditioner takes the ground as odd. The result
is second order for every lateral BC and top, and agrees with back end B's independent implementation to 1e-9
relative (UT-1.11.2).
``z_cross`` is not available with a Dirichlet ground.

With no Dirichlet face the problem is singular. The J-weighted mean of the right-hand side is then removed first
(and reported as ``removed_mean``), and :math:`\phi` comes back with zero J-weighted mean.

**Accuracy next to sloped ground and sloped side walls.** The ghost cell below the ground copies the bottom value of
:math:`\phi`, so the cross terms see :math:`\partial\phi/\partial\zeta = 0` at the ground. On sloped ground the
zero-normal-flux condition allows a nonzero value there, which leaves an O(1) truncation error in the bottom row and a
first-order solution. Its size grows with the square of the slope: over a 120 m hill (slopes up to 21°) the order is 1,
and the fine-mesh error is about 4 times that over a 60 m hill (:doc:`testing` pins both). The ghost cells beyond a
Neumann side face copy :math:`\phi` the same way, so a Neumann side face that the terrain crosses with a slope is first
order next to it too. Away from those faces the layer is second order for every boundary condition. The operator is
kept as ERF has it, so projections match the ERF reference solutions; a second-order treatment, setting those ghost
cells from the zero-flux condition, would be a separate option.

Solver
------

GMRES (restart length ``solver.poisson_restart``, zero initial guess) stops when the residual 2-norm is below
``solver.poisson_reltol`` times the initial one or below ``solver.poisson_abstol``. A solve that reaches
``solver.poisson_max_iter`` stops the run with the residuals and the inputs named. The reported final residual is
recomputed from the solution, not GMRES's estimate. The preconditioner is a fixed linear map, so standard (not
flexible) GMRES applies; flexible GMRES would keep ``solver.poisson_restart`` more solution-sized vectors (about twice
the memory).

**Preconditioner.** ``fft`` applies ERF's ``FFT::PoissonHybrid``: the Laplacian on the flat stretched grid, solved
by FFTs in x and y and a tridiagonal solve in z, with vertical spacing

.. math::

   \Delta z'_k = \Delta z_k\,\sqrt{\beta_h/\beta_v} = \Delta z_k\,\frac{\alpha_h}{\alpha_v},

where the ratio :math:`\beta_v/\beta_h` is ``FaceCoeffs::precond_ratio`` (1 by default; the typical value for
spatially varying coefficients; a new value rebuilds the FFT plan). The residual is divided by the cell's horizontal
coefficient (the mean of its four x- and y-face values) before the FFT solve. A constant factor does not change
right-preconditioned GMRES; a varying one is captured to first order: with :math:`\beta` varying over 10, 100 and
1000 times, the test case takes 24, 29 and 34 iterations, against 19 for :math:`\beta = 1`.

Inputs
------

The options are read by ``read_poisson_options()`` from ``solver.poisson_*`` (:doc:`inputs`). The application reads
them with the mass-consistent initializer (:doc:`initializer`), as do the regression tools (``poisson_projection``,
``wall_distance``); without the initializer the application has no Poisson solve yet, so setting them stops the run.

``PoissonSolver::fill_ghosts(phi, bc)`` fills :math:`\phi`'s ghost cells for the boundary conditions (periodic images,
even copies at Neumann faces, odd copies at Dirichlet faces), as a caller that differentiates :math:`\phi` needs.
