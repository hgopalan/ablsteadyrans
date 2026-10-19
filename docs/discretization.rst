Numerical discretization
========================

One finite-volume discretization serves the Cartesian and the curvilinear form. On flat ground the metric terms
are :math:`z_x = z_y = 0` and :math:`J = \Delta z_k/\Delta\zeta`, and every operator below reduces to its Cartesian
counterpart. The operators are in ``Source/Convection``, ``Source/Diffusion``, ``Source/MassFlux``,
``Source/Forcing`` and ``Source/Residual``.

Grid and control volumes
------------------------

ablsteady uses ERF's Arakawa C-grid:

- :math:`u, v, w` sit on the x-, y- and z-faces;
- :math:`p, \theta, k` sit at cell centres;
- the physical heights ``z_phys_nd`` sit at nodes.

Each field owns the control volume centred on its points. For a velocity component that is the two half cells
beside its face. Boxes are split in x and y only, so every z column lies inside one box.

Each operator is the discrete form of :eq:`curvilinear-divergence`:

.. math::

   L(q) = \frac{1}{J}\left[\frac{\delta G_x}{\delta\xi} + \frac{\delta G_y}{\delta\eta}
          + \frac{\delta G_z}{\delta\zeta}\right],

where :math:`G_d` is the flux through the control-volume faces and :math:`\delta/\delta\xi` is the difference across
the control volume divided by :math:`\Delta\xi`. Values between grid locations are two-point means; derivatives are
compact or centred differences. All are second order on a smooth mesh.

Convection
----------

.. math::

   C(q) = \frac{1}{J}\left[\frac{\delta (F_x q_f)}{\delta\xi} + \frac{\delta (F_y q_f)}{\delta\eta}
          + \frac{\delta (F_z q_f)}{\delta\zeta}\right], \qquad
   F_x = a_x u,\quad F_y = a_y v,\quad F_z = \Omega .

:math:`q_f` is the MUSCL face value on the upwind side, :math:`q_f = q_0 + \tfrac12 s(q_0 - q_{-1}, q_1 - q_0)` for
:math:`F \ge 0`. ``convection.limiter`` sets the slope :math:`s(a, b)`:

.. list-table::
   :header-rows: 1
   :widths: 18 44 38

   * - Limiter
     - Slope
     - Use
   * - ``van_leer``
     - :math:`2ab/(a+b)` where :math:`ab > 0`, else 0
     - Momentum and :math:`\theta` (TVD).
   * - ``van_albada``
     - :math:`\dfrac{(b^2+\epsilon^2)\,a + (a^2+\epsilon^2)\,b}{a^2 + b^2 + 2\epsilon^2}`
     - Smooth; for Newton (Phase 3).
   * - ``upwind``
     - 0
     - First order; the implicit side, and :math:`k`.

:math:`k` uses the bounded form :math:`C(k) - k\,C(1)` with upwind, so every neighbour has a non-negative weight
and :math:`k` stays non-negative. The MUSCL stencil reaches two cells upwind, so the fields carry two ghost layers.

Diffusion
---------

.. math::

   D(q) = \nabla\cdot(\gamma\nabla q), \qquad
   D_c(\mathbf u) = \frac{\partial \tau_{cj}}{\partial x_j}, \qquad G_d = \boldsymbol\tau_c\cdot\mathbf A_d .

- The face flux is the dot product of the physical gradient (chain rule, :doc:`governing_curvilinear`) with the face
  area vector. Every metric cross term is included.
- The momentum stress keeps the transposed term :math:`\partial_c u_j`, which matters where :math:`\nu` varies.
- The cell diffusivity is averaged to the faces. ``turb.face_mean = logarithmic`` uses
  :math:`(a - b)/\ln(a/b)` across vertical pairs.

**Why the logarithmic mean.** The wall-layer diffusivity :math:`\kappa u_*(z + z_0)` is linear in height. The
logarithmic mean then gives the exact face flux of the log profile; the arithmetic mean misses the first face by
about 10 %.

Divergence and pressure gradient
--------------------------------

.. math::

   D(\mathbf U) = \frac1J\left[\frac{\delta(a_x u)}{\delta\xi} + \frac{\delta(a_y v)}{\delta\eta}
   + \frac{\delta\Omega}{\delta\zeta}\right], \qquad
   \Omega = w - z_x\bar u - z_y\bar v ,

where :math:`\bar u, \bar v` are the four-point means around the z-face. The pressure gradient :math:`G` is the
negative adjoint of :math:`D`:

.. math::

   (Gp)_x = \frac{1}{J_u}\left[a_x\frac{\delta p}{\delta\xi} - \frac14\sum z_x \frac{\delta p}{\delta\zeta}\right],
   \qquad
   (Gp)_z = \frac{1}{J_w}\frac{\delta p}{\delta\zeta} .

**Why the adjoint.** :math:`DG` is symmetric negative semidefinite with only the constants in its null space, so
there are no spurious pressure modes. On the ground face the slope term uses the one-sided difference
:math:`-2p_0 + 3p_1 - p_2`.

Sources
-------

- Buoyancy on a w-face is the mean of its two cells, :math:`g(\theta - \theta_\mathrm{ref})/\theta_0`.
- Coriolis uses the four-point means of the other components around each face, as ERF does.

Implicit-side operator
----------------------

The solvers (:doc:`solver`) need a cheap approximate Jacobian. :math:`A_1` is the first-order operator

.. math::

   A_1 q = a_P\,q_P - \sum_{nb} a_{nb}\,q_{nb}, \qquad a_P = \sum_{nb} a_{nb},

over the six neighbours:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Part
     - Weight of the neighbour across a face
   * - Convection
     - Bounded first-order upwind, :math:`\max(\pm F, 0)/(J\Delta\xi_d)`.
   * - Diffusion, x- and y-faces
     - :math:`g\,z_\zeta\,(1 + \delta_{cd})/(J\,\Delta\xi_d^2)`
   * - Diffusion, z-faces
     - :math:`g\,(|\mathbf n|^2 + n_c^2)/(z_\zeta\,J\,\Delta\zeta^2)`

:math:`g` is the face diffusivity and :math:`\mathbf n = \mathbf A_\zeta`. Metric cross terms stay on the right-hand
side.

**Why it works.** Every :math:`a_{nb} \ge 0` and :math:`a_P = \sum a_{nb}`, so :math:`1/\Delta\tau + A_1` is an
M-matrix. Each grid line is a tridiagonal system. :math:`A_1` and the full residual differ, but the fixed point is
the residual's: a converged state has :math:`R = 0`.

Jacobian–vector product
-----------------------

Newton–Krylov (:doc:`solver`) needs :math:`Jv` without forming :math:`J`:

.. math::

   J v \approx \frac{R(q + \epsilon v) - R(q)}{\epsilon}, \qquad
   \epsilon = h\;\frac{1 + \mathrm{mean}\,|\tilde q|}{\mathrm{rms}\,\tilde v} .

- :math:`\tilde q = q/s` and :math:`\tilde v = v/s` are scaled by 10 m/s, 100 m²/s², 10 K and 0.1 m²/s², with
  :math:`\theta` taken as :math:`\theta - \theta_\mathrm{ref}`.
- :math:`h` is ``steady.newton_fd_eps_scale``, by default :math:`\sqrt{\epsilon_\mathrm{mach}}`.

**Why scaled.** Each field then moves by about :math:`h` times its own scale, on any mesh and box layout. A plain
norm would be dominated by :math:`\theta`'s 300 K.

Everything that enters :math:`R` (boundary rows, damping, wall model, closure) is evaluated inside the residual, so
:math:`Jv` sees it.
