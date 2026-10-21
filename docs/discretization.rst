Spatial discretization
======================

The equations are discretized with a second-order finite-volume method. The domain is divided into cells; each
equation is integrated over a cell, which turns it into a balance of fluxes through the cell faces. The same
discretization serves flat ground and terrain: on flat ground :math:`z_x = z_y = 0` and every operator reduces to its
Cartesian form. The operators are in ``Source/Convection``, ``Source/Diffusion``, ``Source/MassFlux``,
``Source/Forcing`` and ``Source/Residual``.

Staggered grid
--------------

GUST uses ERF's staggered (Arakawa C) grid:

- :math:`u, v, w` sit on the x-, y- and z-faces of the cells;
- :math:`p, \theta, k` sit at the cell centres;
- the node heights ``z_phys_nd`` sit at the cell corners.

Each variable has its own control volume centred on its location. For a velocity component this is the two half
cells on either side of its face. The staggering couples pressure and velocity tightly and avoids the checkerboard
pressure modes of a collocated grid.

Boxes are split in x and y only, so every vertical column lies inside one box.

Every operator is a discrete form of :eq:`curvilinear-divergence`:

.. math::

   L(q) = \frac{1}{J}\left[\frac{\delta G_x}{\delta\xi} + \frac{\delta G_y}{\delta\eta}
          + \frac{\delta G_z}{\delta\zeta}\right],

where :math:`G_d` is the flux through the faces of the control volume and :math:`\delta/\delta\xi` is the difference
across the control volume divided by :math:`\Delta\xi`. Values between grid points are two-point averages; all
operators are second order on a smooth mesh.

Convection
----------

.. math::

   C(q) = \frac{1}{J}\left[\frac{\delta (F_x q_f)}{\delta\xi} + \frac{\delta (F_y q_f)}{\delta\eta}
          + \frac{\delta (F_z q_f)}{\delta\zeta}\right], \qquad
   F_x = a_x u,\quad F_y = a_y v,\quad F_z = \Omega .

The face value :math:`q_f` is a MUSCL reconstruction from the upwind side. For :math:`F \ge 0`,

.. math::

   q_f = q_0 + \tfrac12\,s(q_0 - q_{-1},\ q_1 - q_0),

where :math:`s(a, b)` is a slope limiter, set by ``convection.limiter``:

.. list-table::
   :header-rows: 1
   :widths: 18 44 38

   * - Limiter
     - Slope :math:`s(a, b)`
     - Typical use
   * - ``van_leer`` (default)
     - :math:`2ab/(a+b)` where :math:`ab > 0`, else 0
     - Momentum and :math:`\theta`; total-variation diminishing.
   * - ``van_albada``
     - :math:`\dfrac{(b^2+\epsilon^2)\,a + (a^2+\epsilon^2)\,b}{a^2 + b^2 + 2\epsilon^2}`
     - Smooth; suited to Newton methods.
   * - ``upwind``
     - 0
     - First order.

:math:`k` always uses first-order upwinding in the bounded form :math:`C(k) - k\,C(1)`, so :math:`k` stays
non-negative. The MUSCL stencil reaches two cells upwind, which is why the fields carry two ghost layers.

Diffusion
---------

.. math::

   D(q) = \nabla\cdot(\gamma\nabla q), \qquad
   D_c(\mathbf u) = \frac{\partial \tau_{cj}}{\partial x_j}, \qquad G_d = \boldsymbol\tau_c\cdot\mathbf A_d .

- The face flux is the physical gradient (chain rule, :doc:`terrain_coordinates`) dotted with the face area vector.
  All metric cross terms are included.
- The momentum stress keeps the transposed term :math:`\partial_c u_j`, which matters where :math:`\nu` varies.
- Cell diffusivities are averaged to the faces. ``turb.face_mean = logarithmic`` uses :math:`(a - b)/\ln(a/b)`
  across vertical pairs.

The logarithmic mean is exact for the wall layer: there the diffusivity :math:`\kappa u_*(z + z_0)` is linear in
height, and the logarithmic mean gives the exact face flux of the log profile. The arithmetic mean underestimates the
first face flux by about 10 %.

Divergence and pressure gradient
--------------------------------

.. math::

   D(\mathbf U) = \frac1J\left[\frac{\delta(a_x u)}{\delta\xi} + \frac{\delta(a_y v)}{\delta\eta}
   + \frac{\delta\Omega}{\delta\zeta}\right], \qquad
   \Omega = w - z_x\bar u - z_y\bar v ,

where :math:`\bar u, \bar v` are four-point averages around the z-face. The pressure gradient :math:`G` is the
negative adjoint of the divergence:

.. math::

   (Gp)_x = \frac{1}{J_u}\left[a_x\frac{\delta p}{\delta\xi} - \frac14\sum z_x \frac{\delta p}{\delta\zeta}\right],
   \qquad
   (Gp)_z = \frac{1}{J_w}\frac{\delta p}{\delta\zeta} .

With :math:`G = -D^T` (except in the ground row, below), the pressure operator :math:`DG` is symmetric and negative semidefinite, and only constants
lie in its null space. There are no spurious pressure modes. On the ground face the slope term uses the one-sided
difference :math:`-2p_0 + 3p_1 - p_2`.

Source terms
------------

- Buoyancy on a w-face uses the average of its two cells, :math:`g(\theta - \theta_\mathrm{ref})/\theta_0`.
- Coriolis uses four-point averages of the other velocity components around each face, as in ERF.

Approximate Jacobian
--------------------

The implicit solvers (:doc:`solution_algorithm`) need a cheap approximation of the Jacobian
:math:`\partial R/\partial q`. GUST uses the first-order operator

.. math::

   A_1 q = a_P\,q_P - \sum_{nb} a_{nb}\,q_{nb}, \qquad a_P = \sum_{nb} a_{nb},

over the six neighbours of each point:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Term
     - Weight of the neighbour across a face
   * - Convection
     - First-order upwind, :math:`\max(\pm F, 0)/(J\Delta\xi_d)`
   * - Diffusion, x- and y-faces
     - :math:`g\,z_\zeta\,(1 + \delta_{cd})/(J\,\Delta\xi_d^2)`
   * - Diffusion, z-faces
     - :math:`g\,(|\mathbf n|^2 + n_c^2)/(z_\zeta\,J\,\Delta\zeta^2)`

Here :math:`g` is the face diffusivity and :math:`\mathbf n = \mathbf A_\zeta`. The metric cross terms stay on the
right-hand side.

Every :math:`a_{nb} \ge 0`, so :math:`1/\Delta\tau + A_1` is an M-matrix and each grid line is a well-conditioned
tridiagonal system. :math:`A_1` only steers the iteration: a converged state satisfies the full second-order
:math:`R = 0`.

Jacobian–vector product
-----------------------

The Newton–Krylov solver (:doc:`solution_algorithm`) needs products :math:`Jv` with the exact Jacobian
:math:`J = \partial R/\partial q`, without forming :math:`J`. They are computed by a finite difference:

.. math::

   J v \approx \frac{R(q + \epsilon v) - R(q)}{\epsilon}, \qquad
   \epsilon = h\;\frac{1 + \mathrm{mean}\,|\tilde q|}{\mathrm{rms}\,\tilde v} .

- :math:`\tilde q = q/s` and :math:`\tilde v = v/s` are scaled by 10 m/s, 100 m²/s², 10 K and 0.1 m²/s², with
  :math:`\theta` taken as :math:`\theta - \theta_\mathrm{ref}`. Without the scaling, :math:`\theta \approx 300` K would
  dominate the norm.
- :math:`h` is ``steady.newton_fd_eps_scale``, by default :math:`\sqrt{\epsilon_\mathrm{mach}}`, which balances
  truncation and round-off errors.

Everything inside :math:`R` (boundary conditions, damping, wall model, turbulence model) is seen by :math:`Jv`.
The relative error of :math:`Jv` is 2e-8 to 5e-8 in double precision; the :math:`\theta` rows are 3 to 10 times
less accurate (1.6e-7), because :math:`q + \epsilon v` is formed with the absolute :math:`\theta` near 300 K.
