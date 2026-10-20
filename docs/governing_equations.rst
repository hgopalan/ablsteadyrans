Governing equations
===================

ablsteady models the atmosphere as an incompressible fluid under the Boussinesq approximation. Density changes
enter only through buoyancy. The unknowns are the mean velocity :math:`\mathbf u = (u, v, w)`, the kinematic
perturbation pressure :math:`p`, the potential temperature :math:`\theta` and the turbulent kinetic energy (TKE)
:math:`k`. This page gives the equations in Cartesian form; :doc:`terrain_coordinates` gives their form on the
terrain-following mesh.

Conservation equations
----------------------

The steady RANS equations are

.. math::
   :label: rans

   \nabla\cdot\mathbf u &= 0,\\
   \nabla\cdot(\mathbf u\,u_i) &= -\frac{\partial p}{\partial x_i} + \frac{\partial \tau_{ij}}{\partial x_j} + S_i,\\
   \nabla\cdot(\mathbf u\,\theta) &= \nabla\cdot(\gamma_\theta\nabla\theta),\\
   \nabla\cdot(\mathbf u\,k) &= \nabla\cdot(\gamma_k\nabla k) + P_s + P_b - \varepsilon,

with the viscous and turbulent stress

.. math::

   \tau_{ij} = \nu\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right),
   \qquad \nu = \nu_\mathrm{mol} + \nu_t .

- :math:`p = (p_\mathrm{tot} - p_\mathrm{ref}(z))/\rho_0`: the total pressure :math:`p_\mathrm{tot}` relative to a
  hydrostatic reference, divided by the constant reference density :math:`\rho_0`.
- :math:`\nu_\mathrm{mol}`, :math:`\nu_t`: molecular and turbulent viscosity.
- :math:`\gamma_\theta = \alpha_\mathrm{mol} + \nu_t'`: diffusivity of :math:`\theta`, with :math:`\alpha_\mathrm{mol}`
  the molecular and :math:`\nu_t'` the turbulent heat diffusivity.
- :math:`\gamma_k = \nu_\mathrm{mol} + \nu_t/\sigma_k`: diffusivity of :math:`k`, with the Schmidt number
  :math:`\sigma_k` (``turb.sigma_k``, 1).
- :math:`S_i`: momentum sources (below).
- :math:`P_s`, :math:`P_b`, :math:`\varepsilon`: shear production, buoyancy production and dissipation of :math:`k`
  (:doc:`turbulence_model`).

Momentum sources
----------------

.. math::

   S_u &= f_c\,(v\sin\varphi - w\cos\varphi) - f_c\,v_g\sin\varphi + F_x,\\
   S_v &= -f_c\,u\sin\varphi + f_c\,u_g\sin\varphi + F_y,\\
   S_w &= f_c\,u\cos\varphi + g\,\frac{\theta - \theta_\mathrm{ref}(z)}{\theta_0} + F_z .

- :math:`f_c = 4\pi/T`: twice the Earth's rotation rate, with rotation period :math:`T` (86 400 s).
- :math:`\varphi`: latitude. The :math:`\cos\varphi` terms can be switched off (``forcing.coriolis_3d``).
- :math:`(u_g, v_g)`: geostrophic wind. It requires Coriolis forcing and :math:`\varphi \ne 0`.
- :math:`\mathbf F`: a constant pressure-gradient force, as ERF's ``abl_pressure_grad``.
- :math:`g`: gravitational acceleration.
- :math:`\theta_\mathrm{ref}(z)`: the reference sounding, a function of height only.
- :math:`\theta_0`: the Boussinesq reference temperature (300 K).

The reference pressure :math:`p_\mathrm{ref}(z)` balances gravity and the buoyancy of :math:`\theta_\mathrm{ref}`.
A fluid at rest with :math:`\theta = \theta_\mathrm{ref}` and :math:`p = 0` is therefore in exact balance on any
terrain.

The forcing inputs are in the ``forcing`` namespace (:doc:`inputs`).

Steady-state problem
--------------------

Each equation of :eq:`rans` is written as a residual: the left side minus the right side, per unit volume. With the
state :math:`q = (\mathbf u, p, \theta, k)`, a steady solution satisfies

.. math::

   R(q) = 0 .

The solver drives :math:`R` to zero (:doc:`solution_algorithm`). The convergence criteria are given in
:ref:`convergence-criteria`.
