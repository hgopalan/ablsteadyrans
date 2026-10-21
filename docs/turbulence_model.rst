Turbulence model
================

RANS equations describe the mean flow. The effect of turbulent eddies on the mean flow appears as extra mixing,
which is modelled by a turbulent viscosity :math:`\nu_t`. A turbulence model computes :math:`\nu_t` from the flow.

ablsteady uses the one-equation k–L model of Axell and Liungman (2001), referred to below as AL01. The model needs
two quantities:

- the turbulent kinetic energy :math:`k`, how energetic the eddies are, from its own transport equation
  (:eq:`rans`);
- the length scale :math:`L`, how large the eddies are, from an algebraic formula. Near the ground the eddies are
  small; stable stratification makes them smaller still.

The model's point relations are ERF's (``extern_erf/ERF_RANSClosure.H``). The inputs are in the ``turb`` namespace
(:doc:`inputs`).

Eddy viscosity and TKE sources
------------------------------

.. math::

   \nu_t = c_\mu\,k^{1/2}L, \qquad \nu_t' = c_\mu'\,k^{1/2}L, \qquad
   P_s = \nu_t S^2, \qquad P_b = -\nu_t' N^2, \qquad \varepsilon = C_{\mu 0}^3\,\frac{k^{3/2}}{L} .

- :math:`S^2 = 2S_{ij}S_{ij}`: the squared strain rate.
- :math:`N^2 = (g/\theta_0)\,\partial\theta/\partial z`: the squared Brunt–Väisälä frequency.
- :math:`c_\mu`, :math:`c_\mu'`: stability functions (below). The turbulent Prandtl number is
  :math:`c_\mu/c_\mu'`, equal to 1 in neutral flow.
- :math:`C_{\mu 0} = 0.5562` (``turb.Cmu0``) and :math:`\kappa = 0.41` (``turb.kappa``).
- :math:`k` is bounded below by ``turb.tke_floor`` (at least machine epsilon).

Length scale
------------

The length scale starts from the geometric length

.. math::

   l_g = \kappa\,(d + z_0),

capped harmonically at ``turb.max_geom_lscale``. Here :math:`\kappa` is the von Kármán constant, :math:`d` the wall
distance and :math:`z_0` the roughness length. Stratification modifies it:

.. math::

   \frac{1}{L^2} = \frac{1}{l_g^2} + \frac{N^2}{C_b^2 k}\quad (N^2 > 0), \qquad
   L = l_g\left(1 - \frac{C_{\mu 0}^6}{C_b^2}R_t\right)^{1/2}\quad (N^2 < 0),

with the turbulent Richardson number and the constant :math:`C_b` (``turb.Cb``, 0.35)

.. math::

   R_t = \frac{L^2 N^2}{C_{\mu 0}^6 k} .

As in ERF, the unstable length is evaluated once, from :math:`l_g`.

Stability functions
-------------------

From AL01, Eqs. 31–32:

.. math::

   c_\mu = \frac{C_{\mu 0} + 0.108\,R_t}{1 + 0.308\,R_t + 0.00837\,R_t^2}, \qquad
   c_\mu' = \frac{C_{\mu 0}}{1 + 0.277\,R_t} .

The poles of :math:`c_\mu` and :math:`c_\mu'` lie near :math:`R_t = -3.6`. Below ``turb.Rt_crit``, :math:`R_t` is smoothly mapped
towards ``turb.Rt_min`` (the map of Burchard and Petersen) so that it never reaches them.

Neutral wall layer
------------------

Because :math:`\varepsilon` contains :math:`C_{\mu 0}^3`, the model reproduces the neutral log law exactly when
:math:`l_g = \kappa(z + z_0)`, that is, where the cap ``turb.max_geom_lscale`` has no effect (set it very large to
remove it, as tutorial 1 does). In the wall layer

.. math::

   u = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad k = \frac{u_*^2}{C_{\mu 0}^2}, \qquad P_s = \varepsilon

at every height, with :math:`u_*` the friction velocity.

Wall distance
-------------

``turb.wall_distance`` selects the wall distance :math:`d`:

- ``terrain_height`` (default): the height of the cell centre above the ground of its column.
- ``tucker``: the differential wall distance of Tucker (2003),

  .. math::

     \nabla^2\phi = -1,\quad \phi = 0 \text{ on the ground}, \qquad
     d = -|\nabla\phi| + \sqrt{|\nabla\phi|^2 + 2\phi} .

  On flat ground :math:`d` is the height. Over terrain it approximates the true distance; the error is largest at
  convex crests (:doc:`results`).

Model options
-------------

- ``turb.face_mean``: how the cell diffusivities are averaged to the faces (:doc:`discretization`).
- ``turb.nut_relax``: under-relaxation of :math:`\nu_t` between iterations.
- ``turb.frozen_switch``: holds at the floor, during the Newton stage, the cells whose :math:`k` was at the floor
  when Newton started (:doc:`solution_algorithm`). The stable/unstable branch and the :math:`R_t` smoothing are not
  frozen: the length and :math:`R_t` are continuously differentiable across them (below).
- ``turb.ri_gate``: the stratification gate (below).

Smoothness at neutral
---------------------

With :math:`y = l_g^2 N^2/(C_b^2 k)`, the stable length is :math:`l_g(1 + y)^{-1/2}` and the unstable one
:math:`l_g(1 - y)^{1/2}` (before the :math:`R_t` smoothing). Both are :math:`l_g(1 - y/2) + O(y^2)`, so :math:`L`,
:math:`R_t` and :math:`\nu_t` have the same value and the same first derivative in :math:`N^2` on both sides of
:math:`N^2 = 0`; only the second derivative jumps. The smoothing map of :math:`R_t` has slope 1 on both sides of
``turb.Rt_crit``. The kink of the closure is the floor :math:`k_f = \max(k, k_\mathrm{floor})`. (The neutral band
:math:`|N^2| \le` machine epsilon, where :math:`L = l_g`, is a jump in :math:`L` of relative size :math:`y/2` at
its edge: 8e-10 at :math:`k` = 1e-3, but close to 1 at :math:`k` = machine epsilon, where :math:`\nu_t` is
negligible.) Near the floor the closure is smooth but steep: at :math:`k` = 1e-6 the :math:`N^2` scale of
:math:`L`, :math:`C_b^2 k / l_g^2`, is about 4e-10 s⁻², close to the change of :math:`N^2` a finite-difference
Jacobian makes, which then misses :math:`\partial L/\partial N^2` by about 30 %.

Stratification gate
-------------------

*Beyond the plan.* In a nearly neutral layer, under an inversion for example, :math:`N^2` can be round-off or
numerical noise of either sign. With ``turb.ri_gate`` :math:`= Ri_0 > 0` the length :math:`L` and :math:`R_t`
(and through :math:`R_t` the stability functions and the heat diffusivity) use

.. math::

   N^2_\mathrm{eff} = N^2\,\frac{Ri^2}{Ri^2 + Ri_0^2} = N^2\,\frac{N^4}{N^4 + Ri_0^2 S^4}, \qquad Ri = \frac{N^2}{S^2},

with :math:`S^2 = 2 S_{ij}S_{ij}` (:math:`(\partial u/\partial z)^2` in simple shear, so :math:`Ri` is the gradient
Richardson number). Where the shear dominates (:math:`|Ri| \ll Ri_0`) the stratification is suppressed by
:math:`(Ri/Ri_0)^2`; in an inversion (:math:`|Ri| \gg Ri_0`) it is kept to within :math:`(Ri_0/Ri)^2`. The gate
keeps the sign of :math:`N^2`, is monotone in :math:`N^2` (its slope is at most 1.125, at :math:`|Ri| = \sqrt{3}\,
Ri_0`), is infinitely differentiable wherever :math:`S^2 > 0` (at :math:`N^2 = S^2 = 0` it is continuous but not
differentiable), and is exact at :math:`N^2 = 0`, so a neutral run is unchanged bit for bit. The buoyancy
production :math:`P_b = -\nu_t' N^2` keeps the ungated :math:`N^2`, so the buoyancy sink of :math:`k` stays that of
the heat flux :math:`\nu_t'\,\partial\theta/\partial z`; the buoyancy of the momentum equation is not gated
either.

It is a smooth, local form of the neutral switch of Kynema-SGF's ``KLAxell`` (``lengthscale_switch``: below
``meso_sponge_start``, when the deck's surface heat flux is below 1e-5, the neutral length and :math:`R_t = 0`).
That switch is a jump in value and depends on a global input: with KLAxell's defaults (``Cb_stable`` 0.25,
:math:`\lambda` = 30 m) at :math:`N^2` = 1e-4 s⁻², :math:`k` = 0.1 m²/s² and :math:`z` = 100 m it changes :math:`L`
from 7.2 m to 17.3 m. The gate changes stratified answers: in the stable
surface layer :math:`Ri` is about 0.05–0.2, so :math:`Ri_0` = 0.01 lowers :math:`N^2` there by 3.8 % (at
:math:`Ri` = 0.05) to 0.25 % (at 0.2).
It is off by default (:doc:`results`).
