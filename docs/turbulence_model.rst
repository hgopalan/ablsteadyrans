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
- ``turb.frozen_switch``: freezes the branch choices of the model (stable or unstable) before a Jacobian is formed.
