Surface layer
=============

``Source/SurfaceLayer.{H,cpp}`` is the wall model (sub-phase 1.7): Monin–Obukhov similarity (MOST) on every ground
column, the fluxes it sets through the ground, and the wall value of k. The residual uses it when
``ResidualInputs::surface`` is set (:doc:`residual`).

Friction velocity
-----------------

.. math::

   u_* = \frac{\kappa\,U}{\max\left(\ln\dfrac{d + z_0}{z_0} - \psi_m(d/L),\; 1\right)}, \qquad
   \frac1L = -\frac{\kappa\,g\,q}{u_*^3\,\theta_0} .

- :math:`U = \max(|\mathbf U_t|, 0.1\ \mathrm{m/s})`, with :math:`\mathbf U_t` the wind parallel to the ground at the
  reference height (ERF's 0.1 m/s floor).
- :math:`d = h/|A|` is the distance normal to the ground of a point at vertical height :math:`h` above it;
  :math:`|A| = \sqrt{1 + z_x^2 + z_y^2}` uses the column's ground slopes, and the ground is the mean of the column's
  four ground nodes.
- The log law in :math:`d + z_0` matches the closure's mixing length :math:`\kappa(d + z_0)` (sub-phase 1.9) and ERF's
  neutral check, :math:`u = (u_*/\kappa)\ln((z+z_0)/z_0)`. ERF's own MOST uses :math:`\ln(z/z_0)`; the two differ by
  about :math:`z_0/z` (1 % at 10 m for :math:`z_0` = 0.1 m).
- :math:`z_0` is the roughness, :math:`q` the kinematic surface heat flux (K m/s, upward > 0), :math:`g` and
  :math:`\theta_0` the buoyancy options' (:doc:`residual`), :math:`\psi_m` ERF's Businger–Dyer function
  (``extern_erf/ERF_MOSTUtils.H``, unchanged).
- :math:`z/L` is clamped to [-100, 100], and :math:`L` uses :math:`u_* \ge 0.001` m/s, as ERF. The denominator is kept
  at least 1: under strong convection at low wind the relation otherwise has only a negative root.
- ``neutral`` is the log law. ``heat_flux`` is solved by Newton's method to round-off, starting from the neutral value,
  so the residual is a smooth function of the state (ERF iterates a relaxed fixed point to an absolute 1e-5, which makes
  u* a step function of the wind and spoils a finite-difference Jacobian). A column with no positive, finite u* stops the
  run.
- Under stable stratification MOST can have two solutions below a critical wind (a coupled and a collapsed branch,
  e.g. 0.16 and 0.002 m/s near 2.4 m/s for q = −0.01 K m/s at 10 m). Newton from the neutral value finds the coupled one.

**Reference height.** The default is each column's first cell centre. A given ``zref`` is a vertical height above the
ground: the wind is taken on the column, in the interval :math:`[h_k, h_{k+1})` of the cell-centre heights (a height on a
centre belongs to the interval above it, as in ERF#3982, which is already in the pinned ERF), interpolated linearly in
:math:`\ln(d + z_0)`. That is exact for the neutral log law (ERF interpolates linearly in z, 3.5 % off at a 10 m face over
a 10 m first layer); for MOST profiles it is within 0.7 % unstable and 0.5–7 % stable (L = 100–20 m, layers 2–20 m).

Ground fluxes
-------------

The wall model replaces the diffusive fluxes of the ground faces:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Field
     - Ground flux (per unit horizontal index area)
   * - u, v
     - The kinematic wall stress :math:`u_*^2\,\mathbf U_t/U\;|A|`, along the ground wind; on the u (v) control volume
       the mean of the two columns that share the edge. Its component across the ground is zero; its vertical component
       belongs to the half control volume of w at the ground, a boundary row, and is dropped (as in ERF, whose MOST
       stress enters only the u and v rows; ERF leaves out :math:`|A|` unless ``most.terrain_rotate`` is set).
   * - θ
     - :math:`-q\,|A|` (zero for neutral)
   * - k
     - zero

MUSCL convection and the cross-metric and transposed-stress terms of the first cells still read the ghost cells below
the ground; the boundary conditions (sub-phase 1.8) fill them.

Wall k
------

.. math::

   k_w = \frac{\left(u_*^3 + \kappa\,B\,d_0\right)^{2/3}}{C_{\mu 0}^2}, \qquad
   B = \frac{g\,q}{\theta_0} \text{ for } q > 0, \text{ else } 0

(Axell and Liungman 2001, Eq. 16, as ERF; :math:`d_0` the first cell centre's normal distance; :math:`C_{\mu 0}` = 0.5562).
With :math:`a = \max(u_*, 0.01)/d_0`:

- ``dirichlet`` (default, ERF's ``dirichlet_k``): the first cells' k rows are :math:`R_k = a\,(k - k_w)`, and
  ``line_jacobian`` gives them :math:`(0, a, 0)`. ERF writes :math:`k_w` into the first cell instead; the two agree at a
  steady state.
- ``source``: :math:`a\,(k - k_w)` is added to the row and :math:`a` to the diagonal. The rate is close to the natural
  production rate of the first cell, so it only halves an error in first-cell production.

The wall eddy viscosity consistent with this log law is :math:`\kappa\,u_*\,(d + z_0)` (``wall_eddy_viscosity``).

Implicit side
-------------

With the wall model the first cells' rows of :math:`A_1` (:doc:`transport`) lose the ground neighbour (the ground flux no
longer depends on the ghost cell), and u and v gain the wall stress's Picard coefficient
:math:`u_*^2|A|/U\;/(J\,\Delta\zeta)` on the diagonal.

Options
-------

The keys are ``bc.surface_z0``, ``bc.surface_zref``, ``bc.surface_stability``, ``bc.surface_heat_flux`` and
``bc.surface_wall_k`` (``read_surface_layer_options()``, :doc:`boundary`; ERF's ``most.z0``, ``most.zref``,
``most.surf_temp_flux`` and ``dirichlet_k`` are the models). ``Cmu0`` and ``kappa`` are the closure's.

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - ``SurfaceLayerOptions``
     - Default
     - Range and effect
   * - ``z0``
     - 0.1 m
     - > 0, roughness length
   * - ``zref``
     - 0
     - ≤ 0: each column's first cell centre; otherwise a vertical height between the first and last cell-centre
       heights of every column
   * - ``stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``
   * - ``heat_flux``
     - 0 K m/s
     - finite; ``heat_flux`` only (refused with ``neutral``, where it would be ignored)
   * - ``wall_k``
     - ``dirichlet``
     - ``dirichlet`` or ``source``
   * - ``Cmu0``, ``kappa``
     - 0.5562, 0.41
     - > 0 (ERF's defaults)

The surface layer is recomputed from the state at every residual evaluation; there is nothing to checkpoint.
