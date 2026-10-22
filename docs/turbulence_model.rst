Turbulence model
================

RANS equations describe the mean flow. The effect of turbulent eddies on the mean flow appears as extra mixing,
which is modelled by a turbulent viscosity :math:`\nu_t`. A turbulence model computes :math:`\nu_t` from the flow.

GUST uses the one-equation k–L model of Axell and Liungman (2001), referred to below as AL01. The model needs
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

In a nearly neutral layer, under an inversion for example, :math:`N^2` can be round-off or
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

.. _nut-cap:

Eddy-viscosity cap for separated flow
-------------------------------------

A one-equation closure with an algebraic length keeps :math:`\nu_t` large in the shear layer
that leaves a crest, so the flow reattaches too soon behind steep terrain (no separation behind the Bolund
escarpment or the Ishihara–Zhou ridge with the default closure). With ``turb.realizability`` :math:`= C_r > 0` the
new eddy viscosity is capped before the relaxation of ``turb.nut_relax``:

.. math::

   \nu_t = w\,\min\!\Bigl(\nu_t^\mathrm{new},\ C_r\,\frac{k}{S}\Bigr) + (1 - w)\,\nu_t^\mathrm{old},
   \qquad S = \sqrt{2 S_{ij} S_{ij}},

with :math:`w` = ``turb.nut_relax`` and :math:`\nu_t^\mathrm{old}` the stored value of the last iteration. The form
is that of the SST limiter (Bradshaw's :math:`-\overline{uw}/k = a_1 = 0.31`); Durbin's (1996) time-scale bound gives
:math:`\alpha/\sqrt{6} \approx 0.245` for :math:`\alpha = 0.6`. Unlike SST it acts everywhere, not only inside the
boundary layer. With the closure's relations, :math:`\nu_t S/k = C_{\mu 0}^2 \sqrt{P/\varepsilon}` in neutral flow, so the
cap acts where :math:`P/\varepsilon > (C_r/C_{\mu 0}^2)^2`. :math:`C_r = 0.309` is 0.1 % below :math:`C_{\mu 0}^2 = 0.30936`
for the default :math:`C_{\mu 0}`: the cap is active in an equilibrium surface layer too, where it lowers
:math:`\nu_t` by 0–3 % (most near the wall, where :math:`P/\varepsilon` is 1.04–1.05 in the first cells), and it limits
:math:`\nu_t` strongly where the strain far exceeds equilibrium, as in the shear layer leaving a crest. Applying the
cap before the relaxation lets :math:`w < 1` damp the cap's changes between iterations: an imposed two-iteration
oscillation of the cap passes with the factor :math:`w/(2 - w)` (1/3 at :math:`w = 0.5`) instead of 1. Without the
option the cap is the largest real number and the closure is unchanged bit for bit. The derivative the Newton stage
reads follows the cap (:math:`\mathrm{d}\,\mathrm{cap}/\mathrm{d}k = \mathrm{cap}/k`); the column first guess
(``init.first_guess = column``) does not apply it.

.. _hybrid-convection:

Near-ground blend towards upwind
--------------------------------

With the cap, the second-order limited convection can leave a steady run at a fixed,
unconverged state (the rough Ishihara–Zhou ridge stalls at a relative residual of 0.17 with
``van_leer`` and ``van_albada``; first-order upwind converges). With ``convection.hybrid_sigma_max``
:math:`= \sigma_\mathrm{max} > 0` the momentum and θ convection of the residual become

.. math::

   C = C_\mathrm{hi} + \sigma\,(C_\mathrm{up} - C_\mathrm{hi}), \qquad
   \sigma_t = \sigma_\mathrm{max}\,s\!\Bigl(\frac{a - 0.8}{0.2}\Bigr)\,
   s\!\Bigl(\frac{1.5 H_h - d}{0.5 H_h}\Bigr), \qquad a = \frac{\nu_t^\mathrm{new}}{C_r k / S},

- :math:`C_\mathrm{hi}`: the limited scheme of ``convection.limiter``; :math:`C_\mathrm{up}`: first-order upwind
  (the scheme of the implicit side);
- :math:`s(t) = t^2(3 - 2t)` on :math:`[0, 1]` (0 below, 1 above): a smooth step;
- :math:`a`: the ratio of the closure's eddy viscosity to the cap (:math:`a \ge 1` where the cap acts);
- :math:`d`: the closure's wall distance; :math:`H_h` = ``convection.hybrid_height``.

:math:`\sigma` is a cell field, interpolated to the faces of :math:`u, v, w`, updated once per iteration as
:math:`\sigma \leftarrow \sigma + \lambda(\sigma_t - \sigma)` with :math:`\lambda` = ``convection.hybrid_relax``, and
frozen once the relative residual falls below ``convection.hybrid_freeze`` (after at least :math:`2/\lambda` updates
on each grid level), so a converged state satisfies one fixed scheme. The blend is of the convection operators,
not of the face fluxes: it is not conservative where :math:`\sigma` varies (in a 1-D test with :math:`\sigma` rising
from 0 to 0.1 over a few cells, a net source of 0.2–0.4 % of the flux jump; mass is not affected). :math:`\sigma` is
not checkpointed, so ``io.restart`` with the blend aborts.

**Where it acts.** Because :math:`C_r` is at the equilibrium value, :math:`a` is 1.00–1.03 in an equilibrium surface
layer, so :math:`\sigma` reaches :math:`\sigma_\mathrm{max}` over most of the near-ground layer below :math:`H_h`, on the
windward side as well as in the lee; inside a recirculation, where :math:`k` is carried in (:math:`P < \varepsilon`),
:math:`a` can fall below 0.8 and the blend switches off. With :math:`\sigma_\mathrm{max} = 0.1` the added horizontal
numerical diffusion near the ground is of the order of the turbulent viscosity there (about 5 m²/s at 10 m on
Askervein) and changes the approach-flow speed-up by about 0.2 %. A gate that acts only above equilibrium (for
example :math:`a > 1.05`) is not in this version.

Recommended setting
~~~~~~~~~~~~~~~~~~~

For terrain with separated lees: ``turb.realizability = 0.309``, ``turb.nut_relax = 0.5``,
``convection.hybrid_sigma_max = 0.1``, ``convection.hybrid_height`` = 2.5 terrain heights and
``turb.max_geom_lscale = 100`` (on Big Southern Butte the length is uncapped with a moving lid, so that the inflow
does not drift over its 18.7 km fetch). :math:`C_r = 0.309` is :math:`C_{\mu 0}^2` for the default ``turb.Cmu0``
only. On the benchmarks (``Tests/Real_Terrain/WindNinja2019``, ``IshiharaZhou``,
``RUSHIL``) every case reaches the steady verdict with it:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Case
     - Default closure
     - Recommended setting
   * - Ishihara–Zhou ridge, smooth (tunnel reattachment 150–200 m)
     - no reversed flow; U rms 0.80 m/s
     - reversed 42–170 m; 0.30 m/s
   * - Ishihara–Zhou ridge, rough (250–300 m)
     - reversed 59–122 m; 0.70 m/s
     - reversed 47–173 m; 0.58 m/s
   * - RUSHIL Hill 3 (tunnel 0.27a–2.17a at 0.1 H)
     - 0.64a–1.27a; 0.50 m/s
     - 0.58a–1.45a; 0.48 m/s
   * - RUSHIL Hill 5 / Hill 8
     - not steady / attached, 0.23 m/s
     - attached, 0.21 m/s / attached, 0.20 m/s
   * - Askervein, all masts (RMSE)
     - 1.24 m/s
     - 0.97 m/s
   * - Bolund cases 1 / 3 / 4, all masts (RMSE)
     - 1.46 / 0.66 / 0.89 m/s
     - 1.25 / 0.64 / 1.23 m/s

Bolund case 4 (flow from the east, the lee over the escarpment) is slower than measured in the lee with the
cap. On Big Southern Butte the lee separates only on the unsmoothed terrain with a 4 m first layer
(``Tests/Real_Terrain/WindNinja2019/inputs_bsb_srtm``: lee RMSE 2.58 m/s against 3.53 m/s on the paper-like mesh,
with 68 m first cells and smoothed terrain).

The comparison uses data made openly available by N. Wagenbrenner (windninja-tools, the Big Southern Butte
database and Wagenbrenner et al. 2019), T. Ishihara and Q. Zhou (the ridge wind-tunnel data) and ERCOFTAC (the
RUSHIL data of Khurshudyan, Snyder and Nekrasov); we thank them. The case READMEs list the sources and licences.
