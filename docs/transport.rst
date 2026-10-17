Convection and diffusion
========================

``Source/Convection.{H,cpp}``, ``Source/Diffusion.{H,cpp}`` and ``Source/ImplicitOperator.{H,cpp}`` hold the
convection and diffusion operators on the terrain-following C-grid and the first-order implicit-side operator built
from them (sub-phase 1.5). ``Source/Stencil_K.H`` holds the interpolations and differences between grid locations
that all three share. The operators are interior operators: they read ghost cells that the caller fills (periodic
images now, boundary conditions in sub-phases 1.7 and 1.8).

Control volumes
---------------

Every field owns the control volume centred on its points:

- θ, k and p′: the cells;
- velocity component :math:`u_c` (on the c-faces): the halves of the two cells that share its face.

Each operator is the divergence of fluxes :math:`G_d` through the control volume's faces,

.. math::

   L(q) = \frac{1}{J}\left[\frac{\delta G_x}{\delta\xi} + \frac{\delta G_y}{\delta\eta}
          + \frac{\delta G_z}{\delta\zeta}\right],

- :math:`\xi = x`, :math:`\eta = y` and :math:`\zeta` are the index-space coordinates, with spacings
  :math:`\Delta\xi = \Delta x`, :math:`\Delta\eta = \Delta y` and :math:`\Delta\zeta`;
- :math:`\delta/\delta\xi` is the difference across the control volume divided by :math:`\Delta\xi`;
- :math:`J` is the control volume's Jacobian: ERF's ``detJ`` for a cell, the mean of the two cells' values for a
  velocity component.

The faces of :math:`u`'s control volume sit at cell centres in x, on xy-edges in y and on xz-edges in z.

Geometry
--------

The mesh maps :math:`(\xi, \eta, \zeta)` to :math:`(x, y, z(\xi, \eta, \zeta))`. At any location the metric terms are
the slopes of the coordinate surface :math:`z_x = \partial z/\partial\xi`, :math:`z_y = \partial z/\partial\eta` and
the thickness factor :math:`z_\zeta = \partial z/\partial\zeta`. They are differences of the node heights
``z_phys_nd``. Physical derivatives follow from the chain rule,

.. math::

   \frac{\partial q}{\partial x} = q_\xi - \frac{z_x}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial y} = q_\eta - \frac{z_y}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial z} = \frac{q_\zeta}{z_\zeta}.

The flux through a face is the dot product with the face's area vector per unit index area:

.. math::

   \mathbf{A}_x = (z_\zeta, 0, 0), \qquad \mathbf{A}_y = (0, z_\zeta, 0), \qquad
   \mathbf{A}_z = \mathbf{n} = (-z_x, -z_y, 1).

x- and y-faces are vertical; z-faces follow the terrain.

**Index-space differences.** A field seen from another location takes the mean of the two points that bracket that
location in each direction where the two are staggered. A derivative along a direction is the compact difference of
the two bracketing points, or the centred difference :math:`(q_{i+1} - q_{i-1})/2` where field and location are not
staggered. Both are second order on a smooth mesh.

Convection
----------

.. math::

   C(q) = \frac{1}{J}\left[\frac{\delta (F_x q_f)}{\delta\xi} + \frac{\delta (F_y q_f)}{\delta\eta}
          + \frac{\delta (F_z q_f)}{\delta\zeta}\right]

- :math:`F_x = a_x u`, :math:`F_y = a_y v` and :math:`F_z = \Omega = w - z_x u - z_y v` are the transport (volume)
  fluxes through the cell faces per unit index area. They are the fluxes ``face_divergence`` sums (:doc:`poisson`).
  A velocity component's control volume takes the mean of the two cell fluxes it straddles, so a control volume made
  of divergence-free cells is divergence-free too.
- :math:`q_f` is the MUSCL value on the upwind side,
  :math:`q_f = q_0 + \tfrac12 s(q_0 - q_{-1},\, q_1 - q_0)` for :math:`F \ge 0`, and the mirror image otherwise.

The slope :math:`s(a, b)` is set by the limiter:

.. list-table::
   :header-rows: 1
   :widths: 18 44 38

   * - ``Limiter``
     - Slope
     - Use
   * - ``van_leer``
     - :math:`2ab/(a+b)` where :math:`ab > 0`, else 0
     - The standard limiter for momentum and θ (TVD).
   * - ``van_albada``
     - :math:`\dfrac{(b^2+\epsilon^2)\,a + (a^2+\epsilon^2)\,b}{a^2 + b^2 + 2\epsilon^2}`
     - The smooth limiter for Newton (Phase 3). Infinitely differentiable for :math:`\epsilon > 0`;
       :math:`\epsilon` is ``ConvectionOptions::smooth_eps`` (default :math:`10^{-5}`): an absolute value in the
       units of q, set per field and not scaled with the mesh. The default is far below the cell differences, so the
       limiter acts like van Albada with :math:`\epsilon = 0` and clips smooth extrema as van Leer does; Phase 3 sets
       :math:`\epsilon` per field (for example the Venkatakrishnan scaling :math:`\epsilon^2 = q_\mathrm{ref}^2 (K
       h/L_\mathrm{ref})^3`).
   * - ``upwind``
     - 0
     - First-order upwind: the implicit side, and k.

The **bounded form** :math:`C(q) - q\,C(1)` (``ConvectionOptions::bounded``) sums
:math:`F\,(q_f - q)` over the faces. With upwind every neighbour enters with a non-negative weight, so a
pseudo-time step keeps :math:`k \ge 0`. k uses it, with upwind. Where the transport is divergence-free the two forms
agree.

The MUSCL stencil reaches two cells upwind of a face, so the solution fields carry two ghost cells
(``ngrow_state = 2``). Every ghost layer must hold finite values, including where F = 0: a ground face with
:math:`F_z = 0` carries no flux, but its reconstruction still reads both layers below the ground, and
:math:`0 \times \mathrm{NaN}` is NaN. The boundary-condition code fills both layers.

Diffusion
---------

Scalars (θ, k):

.. math::

   D(q) = \nabla\cdot(\gamma\nabla q), \qquad G_d = \gamma\,(\nabla q\cdot\mathbf{A}_d).

- :math:`\gamma` is a cell-centred diffusivity, averaged to each face. The closure (sub-phase 1.9) supplies it:
  :math:`\nu_t/\mathrm{Pr}_t(R_t)` for θ and :math:`\nu_t/\sigma_k` for k.

Momentum, with the full stress tensor:

.. math::

   D_c(\mathbf{u}) = \partial_j\tau_{cj}, \qquad
   \tau_{cj} = \nu\,(\partial_j u_c + \partial_c u_j), \qquad G_d = \boldsymbol{\tau}_c\cdot\mathbf{A}_d.

- :math:`\nu` is the cell-centred viscosity (:math:`\nu + \nu_t`), averaged to each face.
- The transposed term :math:`\partial_c u_j` matters where :math:`\nu` varies; for constant :math:`\nu` it is
  :math:`\nu\,\partial_c(\nabla\cdot\mathbf{u})`.

Both operators include every metric cross term.

Implicit side
-------------

The defect-correction iteration of sub-phase 1.14 solves :math:`(V/\Delta\tau + A_1)\,\delta q = -R(q)`. The
first-order operator :math:`A_1` is

.. math::

   A_1 q = a_P\,q_P - \sum_{nb} a_{nb}\,q_{nb}, \qquad a_P = \sum_{nb} a_{nb},

over the six neighbours, with the weights split by direction:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Part
     - Weight of the neighbour across a face
   * - Convection
     - Bounded first-order upwind: :math:`\max(F, 0)/(J\Delta\xi_d)` on the low face and
       :math:`\max(-F, 0)/(J\Delta\xi_d)` on the high face.
   * - Diffusion, x- and y-faces
     - :math:`g\,z_\zeta\,(1 + \delta_{cd})/(J\,\Delta\xi_d^2)`
   * - Diffusion, z-faces
     - :math:`g\,(|\mathbf{n}|^2 + n_c^2)/(z_\zeta\,J\,\Delta\zeta^2)`

- :math:`g` is the face diffusivity (:math:`\gamma` or :math:`\nu`).
- The :math:`\delta_{cd}` and :math:`n_c^2` terms are the part of the transposed stress along :math:`u_c`'s own
  derivative. A scalar has neither.
- The metric cross terms and the transposed stress across components stay on the right-hand side.

Every :math:`a_{nb} \ge 0` and :math:`a_P = \sum a_{nb}` for any transport and diffusivity. So
:math:`V/\Delta\tau + A_1` is an M-matrix (diagonally dominant, off-diagonals not positive). Line d of a line solver
is tridiagonal in :math:`(-a_{d-},\, a_P,\, -a_{d+})`. The coefficients are stored as seven components
(``coef_P``, ``coef_XM`` … ``coef_ZP``). Weights that reach outside the domain are kept; the boundary-condition code
folds them in, and replaces the rows of boundary points (the ground face of w, whose control volume is half a cell
here) in both :math:`A_1` and the residual.

Notes for the code that uses :math:`A_1`:

- :math:`A_1` linearises the bounded convection form for every field; the residual uses the conservative form for
  momentum and θ. They differ by :math:`q\,C(1)`, which vanishes as the transport becomes divergence-free, so
  defect correction converges to the residual's solution.
- The vertical upwind weights are the implicit vertical advection of sub-phase 1.14 (its task 2).
- An independent reference for :math:`A_1` is first-order bounded convection plus the principal diffusion; on a flat
  mesh the principal part is the whole diffusion operator. The unit tests check both.

Accuracy
--------

On a manufactured flow over a 120 m cosine hill (slopes up to 0.38), BTF mapping and geometric layers, at three
resolutions (16 × 12 × 8 to 64 × 48 × 32 cells):

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Operator
     - L2 order (finest pair)
     - Max-norm order (finest pair)
   * - Convection, unlimited MUSCL: q, u, v, w
     - 1.99 – 2.10
     - 1.94 – 2.10
   * - Convection, ``van_leer`` / ``van_albada``
     - 1.5 – 2.2
     - about 1 (extrema)
   * - Convection, ``upwind``
     - 1.00
     - 0.96
   * - Scalar diffusion
     - 1.99
     - 1.99
   * - Momentum diffusion: u, v, w
     - 1.99 – 2.00
     - 1.96 – 1.99

A limiter cuts the slope at smooth extrema. The error there is first order, so the limited schemes have a max-norm
order of about 1; their L1 order tends to 2.

**Ground ghosts.** The test mesh continues the layers smoothly below the ground and above the top. ERF's terrain
code sets the ghost nodes by linear extrapolation, :math:`z_{-1} = 2z_0 - z_1`. That is a kink in the mapping on a
stretched mesh. With it, the centred differences that straddle the ground (for example :math:`\partial u/\partial z`
at the first cell centre, in the w control volume) carry a relative error of :math:`(r-1)/4` for a layer ratio
:math:`r`, which does not shrink when the mesh is refined at fixed ratio. On the test meshes with ERF's ghost nodes,
convection is unchanged, while momentum diffusion drops to a w max-norm order of 1.30 / 0.73 (L2 1.92 / 1.75, the
largest error at the first interior w point) and a u max-norm order of 1.86. Sub-phase 1.7 replaces the ground
fluxes with the wall model and sets the ghost values; it also has to supply :math:`\partial u/\partial z` where the
first-cell stencils use it, since a log profile differenced across the ground is not usable.
