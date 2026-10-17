Residual
========

``Source/Residual.{H,cpp}`` assembles the steady residual from the operators of :doc:`transport`, the divergence and
pressure gradient in ``Source/MassFlux.{H,cpp}`` and the sources in ``Source/Forcing.{H,cpp}`` (sub-phase 1.6). It
also gives the two linearisations the solvers need: the line coefficients of the first-order implicit side, and a
Jacobian–vector product.

Equations
---------

.. math::

   R_c &= C_c(\mathbf u) - D_c(\mathbf u) + (Gp)_c - S_c \qquad (c = x, y, z)\\
   R_p &= D(\mathbf U)\\
   R_\theta &= C(\theta) - \nabla\cdot(\gamma_\theta\nabla\theta)\\
   R_k &= C_b(k) - \nabla\cdot(\gamma_k\nabla k)

- :math:`C`, :math:`D_c` and the scalar diffusion are the operators of :doc:`transport`; :math:`C_b` is the bounded
  first-order upwind form, which k always uses.
- :math:`p` is the kinematic perturbation pressure :math:`(p - p_\mathrm{ref})/\rho_0`.
- :math:`S` is the source: buoyancy on w, Coriolis and geostrophic forcing.
- :math:`\gamma_\theta = \kappa + \nu_t/\mathrm{Pr}_t` and :math:`\gamma_k = \nu + \nu_t/\sigma_k` are cell-centred
  inputs, like :math:`\nu + \nu_t` for momentum; the closure supplies them (sub-phase 1.9).

**Conventions.**

- Every residual is per unit volume: m/s² for momentum, 1/s for continuity, K/s for θ, m²/s³ for k. A steady state
  has :math:`R = 0`; the solvers solve for :math:`\delta q` against :math:`-R`.
- The residual is an interior operator. The state's ghost cells (two layers) must be filled by the caller. The
  boundary rows (the ground and top faces of w, the faces on non-periodic sides) are left for the boundary-condition
  code (sub-phase 1.8). The closure sources of k come with sub-phase 1.9, the wall model with 1.7.
- The faces shared by two boxes, or periodic images, get one value (``OverrideSync``), so per-box updates cannot let
  the copies drift apart.
- **Everything that belongs to R enters inside** ``residual()`` (through ``ResidualInputs``), so that ``jvp``, which
  differentiates ``residual()``, sees it. The diffusivities are inputs held fixed: until the closure computes
  :math:`\nu_t` from the state inside ``residual()``, ``jvp`` is the Jacobian with :math:`\nu_t` frozen.

Divergence and pressure gradient
--------------------------------

The state holds the Cartesian velocity. A z-face follows the terrain, so the flux through it per unit horizontal
index area is the contravariant velocity

.. math::

   \Omega = w - z_x\,\bar u - z_y\,\bar v ,

- :math:`z_x`, :math:`z_y` are the face's slopes, from its four nodes (the slopes ERF's ``TerrainPoisson`` uses at
  z-faces);
- :math:`\bar u`, :math:`\bar v` are the means of the four u (v) points around the face.

ERF converts w to :math:`\Omega` the same way (``OmegaFromW``, with a fourth-order slope), takes the divergence with
:math:`\Omega`, adds the projection fluxes to :math:`(u, v, \Omega)` and converts back. ``w_from_omega`` is the
inverse here.

**Divergence.** :math:`D(\mathbf U) = \frac1J\left[\frac{\delta(a_x u)}{\delta\xi} + \frac{\delta(a_y v)}{\delta\eta}
+ \frac{\delta\Omega}{\delta\zeta}\right]`, computed by the Poisson layer's ``face_divergence``. The continuity
residual is therefore bit for bit the right-hand side the Poisson layer inverts. A projection that adds the Poisson
fluxes to :math:`(u, v, \Omega)` and recovers w with ``w_from_omega`` leaves exactly the solve's residual as
divergence (UT-1.6.2); this is how the SIMPLEC correction (sub-phase 1.15) keeps continuity.

**Pressure gradient.** :math:`G` is the negative adjoint of :math:`D` with the control-volume weights,

.. math::

   \sum_\text{faces} V_f\,(Gp)_f\,U_f = -\sum_\text{cells} V_c\,p_c\,D(\mathbf U)_c ,

which gives

.. math::

   (Gp)_x = \frac{1}{J_u}\left[a_x\frac{\delta p}{\delta\xi}
            - \frac14\sum z_x \frac{\delta p}{\delta\zeta}\right], \qquad
   (Gp)_z = \frac{1}{J_w}\frac{\delta p}{\delta\zeta}

(the sum runs over the four z-faces whose :math:`\bar u` contains this u). To second order it is the physical
gradient, :math:`\partial p/\partial x = p_\xi - (z_x/z_\zeta)\,p_\zeta`. :math:`D G` is symmetric negative
semidefinite with only the constants in its null space, so no spurious pressure modes (checked in review on periodic
hills with slopes up to 1).

- **Ground row.** On the ground face, :math:`\delta p/\delta\zeta` in the slope term is the second-order one-sided
  difference :math:`-2p_0 + 3p_1 - p_2` of the three cells above: the ghost below the ground is not read. An even
  ghost (the Poisson layer's Neumann fill) would halve the slope term in the bottom row, an error of 7–17 % of
  :math:`\partial p/\partial x` at slopes 0.3–0.9. Adjointness then holds everywhere but that row.
- **Not ERF's Poisson operator.** :math:`D G` and ERF's ``TerrainPoisson`` operator, which the SIMPLEC correction
  inverts, differ in the cross terms: :math:`DG` averages :math:`p_\zeta` over three faces and does not see the
  :math:`2\Delta\zeta` vertical mode. On a uniform slope s the eigenvalues of :math:`A_\mathrm{ERF}^{-1} D G` lie in
  :math:`[1/(1+s^2), 1]`, so a pressure correction with :math:`\alpha_p = 1` contracts the momentum–pressure
  mismatch by :math:`s^2/(1+s^2)` per outer iteration: 0.08 at s = 0.3, 0.26 at 0.6, 0.5 at 1. It slows the outer
  iteration on steep slopes; it does not destabilise it.

Sources
-------

**Buoyancy** about a hydrostatic reference: :math:`S_w = g\,(\theta - \theta_\mathrm{ref}(z))/\theta_0`, averaged
from the two cells of the w-face.

- :math:`\theta_\mathrm{ref}(z)` is a function of physical height alone, given as a 1-D profile
  (``ResidualInputs::theta_ref``: heights and temperatures, linear in between, constant beyond). The module builds the
  cell values itself: a 3-D reference would drop the horizontal gradient of :math:`p_\mathrm{ref}`. The profile
  should follow the far-field sounding; it is an input, so a restart rebuilds the same field.
- :math:`p_\mathrm{ref}(z)` balances gravity and the buoyancy of :math:`\theta_\mathrm{ref}`
  (:math:`dp_\mathrm{ref}/dz = -\rho_0 g\,[1 - (\theta_\mathrm{ref} - \theta_0)/\theta_0]`) and never enters the
  residual; the total pressure for output is :math:`p_\mathrm{ref} + \rho_0 p`.
- :math:`\theta_0` is the Boussinesq reference temperature (``BuoyancyOptions::theta0``, 300 K; ``gravity`` 9.81).
- A state at rest with :math:`\theta = \theta_\mathrm{ref}` and :math:`p = 0` is in exact discrete balance on any
  terrain. Where θ departs from :math:`\theta_\mathrm{ref}`, the hydrostatic part of :math:`p` is differenced along
  the sloping coordinate surfaces, with a truncation error that grows with the departure and the slope: 0.5 % of the
  departure's buoyancy over a 300 m hill (UT-1.6.4). On geometric layers with ratio r it does not shrink with
  :math:`\Delta x`: about :math:`s\,B\,(r-1)^2/(4r)` (review estimate), :math:`6.6\times10^{-4}` m/s² for a 5 K
  inversion at slope 0.3 and :math:`10^{-2}` m/s² at slope 0.9 (for scale, :math:`fU \approx 10^{-3}` m/s²).

**Coriolis**, as in ERF (``ERF_MakeMomSources.cpp``), with :math:`f_c = 4\pi/T`:

.. math::

   S_u = f_c\,(\bar v\sin\varphi - \bar w\cos\varphi), \qquad S_v = -f_c\,\bar u\sin\varphi, \qquad
   S_w = f_c\,\bar u\cos\varphi .

- :math:`\varphi` is the latitude, :math:`T` the rotation period (ERF's default, 86400 s); a bar is the mean of the
  four points around the face.
- The :math:`\cos\varphi` terms come with ``three_d`` (ERF's ``coriolis_3d``, default on).
- Geostrophic forcing, horizontal as in ERF, adds :math:`-f_c v_g \sin\varphi` to :math:`S_u` and
  :math:`f_c u_g\sin\varphi` to :math:`S_v`. It needs Coriolis on and :math:`\varphi \ne 0`; a geostrophic wind given
  with the forcing off is refused.
- Over terrain the u and v control volumes differ, so the discrete Coriolis force does a little work; it is second
  order (:math:`10^{-5}` of the Coriolis magnitude at 64 cells, review estimate) and harmless for a steady solve.

``CoriolisOptions``, ``BuoyancyOptions`` and the profile are range-checked; there are no input keys yet (the driver
adds them).

Linearisations
--------------

``line_jacobian(q, in, eq, dir, …)`` returns, for one equation (u, v, w, θ, k) and one direction, the tridiagonal
coefficients :math:`(-a_{d-},\; a_P,\; -a_{d+})` of the first-order implicit side :math:`A_1` (:doc:`transport`),
built with the state's transport fluxes (an overload takes precomputed fluxes). The neighbours across the other
directions are lagged. The caller adds :math:`1/\Delta\tau` (residuals and :math:`A_1` are per unit volume; the plan's
:math:`V/\Delta\tau` is the same term integrated over a cell) and source Jacobians to the diagonal.

``jvp(q, v, …)`` returns :math:`J v \approx [R(q + \epsilon v) - R(q)]/\epsilon` with

.. math::

   \epsilon = \sqrt{\epsilon_\mathrm{mach}}\;\frac{1 + \mathrm{mean}\,|\tilde q|}{\mathrm{rms}\,\tilde v},

over every owned point of every field (shared faces once), on scaled variables :math:`\tilde q = q/s`,
:math:`\tilde v = v/s` with θ taken as :math:`\theta - \theta_\mathrm{ref}`. The scales are ``JvpScales`` (10 m/s,
10 K, 0.1 m²/s², 100 m²/s²). A perturbation along a typical v moves each value by about
:math:`\sqrt{\epsilon_\mathrm{mach}}` times its own scale, whatever the field and the mesh size; a plain norm would be
dominated by θ's 300 K and grow with the number of cells. The step is the same on any box layout.
