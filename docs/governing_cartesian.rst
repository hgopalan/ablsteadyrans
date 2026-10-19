Governing equations in Cartesian form
=====================================

ablsteady solves the steady Reynolds-averaged equations of an incompressible, Boussinesq atmosphere. The unknowns
are the velocity :math:`\mathbf u = (u, v, w)`, the kinematic perturbation pressure :math:`p`, the potential
temperature :math:`\theta` and the turbulent kinetic energy :math:`k`.

Conservation equations
----------------------

.. math::
   :label: rans

   \nabla\cdot\mathbf u &= 0,\\
   \nabla\cdot(\mathbf u\,u_i) &= -\frac{\partial p}{\partial x_i} + \frac{\partial \tau_{ij}}{\partial x_j} + S_i,\\
   \nabla\cdot(\mathbf u\,\theta) &= \nabla\cdot(\gamma_\theta\nabla\theta),\\
   \nabla\cdot(\mathbf u\,k) &= \nabla\cdot(\gamma_k\nabla k) + P_s + P_b - \varepsilon,

with the stress

.. math::

   \tau_{ij} = \nu\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right),
   \qquad \nu = \nu_\mathrm{mol} + \nu_t .

- :math:`p = (p_\mathrm{tot} - p_\mathrm{ref}(z))/\rho_0` is the pressure about a hydrostatic reference.
- :math:`\gamma_\theta = \alpha_\mathrm{mol} + \nu_t'` and :math:`\gamma_k = \nu_\mathrm{mol} + \nu_t/\sigma_k` are
  the diffusivities of :math:`\theta` and :math:`k`.
- :math:`S_i` is the momentum source, and :math:`P_s, P_b, \varepsilon` are the sources of :math:`k` (below).

**What "steady" means.** Each equation is written as a residual, left side minus right side, per unit volume. A
steady solution has :math:`R(q) = 0` for the state :math:`q = (\mathbf u, p, \theta, k)`. :doc:`solver` gives the
stopping test.

Momentum sources
----------------

.. math::

   S_u &= f_c\,(v\sin\varphi - w\cos\varphi) - f_c\,v_g\sin\varphi + F_x,\\
   S_v &= -f_c\,u\sin\varphi + f_c\,u_g\sin\varphi + F_y,\\
   S_w &= f_c\,u\cos\varphi + g\,\frac{\theta - \theta_\mathrm{ref}(z)}{\theta_0} + F_z .

- :math:`f_c = 4\pi/T` with the rotation period :math:`T` (86 400 s), and :math:`\varphi` is the latitude. The
  :math:`\cos\varphi` terms are switched by ``forcing.coriolis_3d``.
- :math:`(u_g, v_g)` is the geostrophic wind. It needs Coriolis on and :math:`\varphi \ne 0`.
- :math:`\mathbf F` is a constant pressure-gradient force (ERF's ``abl_pressure_grad``).
- :math:`\theta_\mathrm{ref}(z)` is the sounding, a function of height alone, and :math:`\theta_0` (300 K) is the
  Boussinesq reference temperature.

**Why a reference profile.** :math:`p_\mathrm{ref}` balances gravity and the buoyancy of :math:`\theta_\mathrm{ref}`.
A state at rest with :math:`\theta = \theta_\mathrm{ref}` and :math:`p = 0` is then in exact balance on any terrain.

Turbulence closure
------------------

The closure is the one-equation k–L model of Axell and Liungman (2001, AL01), with ERF's point relations
(``extern_erf/ERF_RANSClosure.H``).

.. math::

   \nu_t = c_\mu\,k^{1/2}L, \qquad \nu_t' = c_\mu'\,k^{1/2}L, \qquad
   P_s = \nu_t S^2, \qquad P_b = -\nu_t' N^2, \qquad \varepsilon = C_{\mu 0}^3\,\frac{k^{3/2}}{L} .

- :math:`S^2 = 2S_{ij}S_{ij}` is the strain rate squared, and :math:`N^2 = (g/\theta_0)\,\partial\theta/\partial z` is
  the Brunt–Väisälä frequency squared.
- :math:`C_{\mu 0} = 0.5562`. The turbulent Prandtl number is :math:`c_\mu/c_\mu'`, 1 in neutral air.

The length scale :math:`L` starts from the geometric length :math:`l_g = \kappa(d + z_0)`, capped harmonically at
``turb.max_geom_lscale``, with :math:`d` the wall distance and :math:`z_0` the roughness. Stratification changes it:

.. math::

   \frac{1}{L^2} = \frac{1}{l_g^2} + \frac{N^2}{C_b^2 k}\ \ (N^2 > 0), \qquad
   L = l_g\left(1 - \frac{C_{\mu 0}^6}{C_b^2}R_t\right)^{1/2}\ \ (N^2 < 0), \qquad
   R_t = \frac{L^2 N^2}{C_{\mu 0}^6 k} .

The stability functions are (AL01 Eqs. 31–32)

.. math::

   c_\mu = \frac{C_{\mu 0} + 0.108\,R_t}{1 + 0.308\,R_t + 0.00837\,R_t^2}, \qquad
   c_\mu' = \frac{C_{\mu 0}}{1 + 0.277\,R_t} .

- Below ``turb.Rt_crit`` the Richardson number :math:`R_t` is smoothed towards ``turb.Rt_min`` (Burchard and
  Petersen's map). The poles of :math:`c_\mu` sit near :math:`R_t = -3.61`.
- The unstable length is evaluated once, from :math:`l_g`, as ERF does.
- :math:`k` is floored at ``turb.tke_floor`` (at least machine epsilon).

**Why this closure gives the log law.** With :math:`C_{\mu 0}^3` in :math:`\varepsilon`, the neutral wall layer has
:math:`u = (u_*/\kappa)\ln((z+z_0)/z_0)`, :math:`k = u_*^2/C_{\mu 0}^2` and :math:`P_s = \varepsilon` at every height.

Wall distance
-------------

``turb.wall_distance`` selects :math:`d`:

- ``terrain_height`` (default): the height of the cell centre above its column's ground.
- ``tucker``: Tucker's (2003) differential distance,

  .. math::

     \nabla^2\phi = -1,\quad \phi = 0 \text{ on the ground}, \qquad
     d = -|\nabla\phi| + \sqrt{|\nabla\phi|^2 + 2\phi} .

  On flat ground :math:`d` is the height. Over terrain it approximates the distance and is weakest at convex crests.

References
----------

- Axell, L. B., Liungman, O. (2001). *Environ. Fluid Mech.* 1, 71–106.
- Tucker, P. G. (2003). Differential equation-based wall distance computation for DES and RANS. *J. Comput. Phys.*
  190, 229–248.
