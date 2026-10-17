Poisson layer
=============

Three parts of ablsteady solve a Poisson problem on the terrain-following mesh: the mass-consistent initializer, the
SIMPLEC pressure correction and the Tucker wall distance. They are to share one layer, ``Source/Poisson.{H,cpp}``:

.. code-block:: c++

   PoissonSolver solver(mesh, geom, z_phys_nd, read_poisson_options());
   PoissonResult r = solver.solve(rhs, phi, face_coeffs, bc, &fluxes);

Back end A (sub-phase 1.10) is ERF's own operator, ``TerrainPoisson``, copied unchanged into ``extern_erf/``
(:doc:`extraction`), solved by ``amrex::GMRES`` with ERF's FFT preconditioner. Back end B (``MLTerrainPoisson``)
follows in sub-phase 1.11.

This version serves the initializer and SIMPLEC. The Tucker distance needs :math:`\phi = 0` on the ground, and here the
ground is always a zero-flux wall and every boundary condition is homogeneous; sub-phase 1.11 adds a Dirichlet
ground (ERF's ``terrpoisson_flux_zlo_dir`` is in the copied kernel file, unused so far).

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
:math:`(\beta_z - \beta_{zc})` times the cross part by a separate kernel (``Source/PoissonCross_K.H``, ERF's formula).
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
     - Zero total flux through the face, the terrain cross term included (ERF#4272, merged as ERF#4279). The ground
       is always Neumann.
   * - ``dirichlet``
     - :math:`\phi = 0` on the face (an outflow or open face).

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

The options are read by ``read_poisson_options()`` from ``solver.poisson_*`` (:doc:`inputs`). In this version only
the regression tool ``poisson_projection`` reads them; the application has no Poisson solve yet, and these keys stop
an ablsteady run like any other key it does not read.
