Boundary conditions
===================

The computational domain is a box with six faces: west and east (``xlo``, ``xhi``), south and north (``ylo``,
``yhi``), the ground (``zlo``) and the top (``zhi``). Each face gets one boundary type. The typical ABL setup is
periodic sides, a rough ground described by a wall model, and a slip top. The boundary types are ERF's. The code is
in ``Source/BoundaryConditions``, ``Source/InflowProfile`` and ``Source/SurfaceLayer``; the inputs are in the ``bc``
namespace (:doc:`inputs`).

Boundary types
--------------

.. list-table::
   :header-rows: 1
   :widths: 18 14 68

   * - Type
     - Faces
     - Condition
   * - ``periodic``
     - x, y
     - The flow leaving one side enters the opposite side. Used where ``mesh.is_periodic = 1``.
   * - ``inflow``
     - x, y
     - :math:`u, v, \theta, k` from the inflow profile at their height above the ground.
   * - ``outflow``
     - x, y
     - Zero normal gradient for every field; :math:`p = 0`.
   * - ``slip_wall``
     - all
     - No flow through the face; no stress along it. At the top, :math:`\theta` may have a given gradient
       (``bc.zhi_theta_grad``).
   * - ``no_slip_wall``
     - all
     - The velocity equals the wall velocity (0, or a moving lid set by ``bc.<face>_wall_velocity``).
   * - ``surface_layer``
     - ground
     - Monin–Obukhov wall model (below). The default ground.
   * - ``geostrophic``
     - top
     - :math:`u, v` held at the geostrophic wind.

The defaults are periodic sides, a ``surface_layer`` ground and a ``slip_wall`` top, as in ERF's canonical ABL
cases.

**No flow through the ground and the top.** On a terrain-following mesh the velocity normal to the ground is the
contravariant velocity :math:`\Omega` (:doc:`terrain_coordinates`), not :math:`w`. The ground and top values of
:math:`w` are therefore

.. math::

   w_\mathrm{face} = z_x\,\bar u + z_y\,\bar v ,

which gives :math:`\Omega = 0`. Over flat ground :math:`w_\mathrm{face} = 0`.

Wall model
----------

The first cell above the ground is usually much thicker than the thin layer where the wind drops to zero. A wall
model bridges this layer: it gives the stress at the ground from the wind in the first cell, using Monin–Obukhov
similarity theory (MOST). In each column,

.. math::

   u_* = \frac{\kappa\,U}{\max\left(\ln\dfrac{d + z_0}{z_0} - \psi_m(d/L_O),\; m\right)}, \qquad
   \frac{1}{L_O} = -\frac{\kappa\,g\,q}{u_*^3\,\theta_0} .

- :math:`u_*`: friction velocity.
- :math:`U = \max(|\mathbf U_t|, 0.1~\mathrm{m/s})`: wind speed parallel to the ground at the reference height
  (``bc.surface_zref``; by default the first cell centre).
- :math:`d = h/|A|`: distance from the ground along its normal, with :math:`|A| = \sqrt{1 + z_x^2 + z_y^2}`.
- :math:`z_0`: roughness length (``bc.surface_z0``).
- :math:`q`: kinematic surface heat flux (``bc.surface_heat_flux``).
- :math:`L_O`: Obukhov length.
- :math:`\psi_m`: ERF's Businger–Dyer stability function, or QNSE's on the stable side (below).
- :math:`m`: a floor on the denominator (below).

With ``bc.surface_stability = heat_flux``, :math:`u_*` and :math:`L_O` are solved together by Newton's method to
round-off. (ERF uses a relaxed fixed point to 1e-5; a converged :math:`u_*` keeps the residual smooth, which a
finite-difference Jacobian needs.)

**Denominator floor.** Without a floor, strong convection at low wind gives :math:`u_* < 0`.
``bc.surface_denominator_floor`` selects the floor:

- ``log_law`` (default): :math:`m = \min(\ln((d + z_0)/z_0),\, 1)`. In neutral flow the floor never acts, so
  :math:`u_*` is the log-law value at every height.
- ``one``: :math:`m = 1`. In neutral flow this lowers :math:`u_*` when the reference point is below
  :math:`(e - 1)\,z_0 \approx 1.72\,z_0`, for example a 1 m first layer over :math:`z_0 = 0.5` m.

The two give identical results where :math:`\ln((d + z_0)/z_0) \ge 1`.

.. _stable-functions:

**Stable functions.** ``bc.surface_stable_functions`` selects the stable side
(:math:`\zeta = z/L > 0`) of the flux-gradient functions :math:`\phi_m, \phi_h` and their integrals
:math:`\psi = \int_0^\zeta (1 - \phi(x))/x\,dx`. Both use one form,

.. math::

   \phi(\zeta) = \frac{1 + a\,\zeta}{1 + b\,\zeta}, \qquad
   \psi(\zeta) = -\frac{a - b}{b}\,\ln(1 + b\,\zeta) .

- ``businger_dyer`` (default): :math:`a = 5`, :math:`b = 0`, so :math:`\phi = 1 + 5\zeta` and
  :math:`\psi = -5\zeta` (Dyer 1974, ERF's functions).
- ``qnse``: the quasi-normal scale elimination functions (Sukoriansky, Galperin and Perov 2005), as ERF's
  ``enable_qnse_stable_functions``: :math:`a_m, b_m` (``bc.surface_qnse_am``, ``_bm``) for momentum and
  :math:`a_h, b_h` (``_ah``, ``_bh``) for heat, ERF's defaults 2.5 and 0.2. :math:`\phi` tends to
  :math:`a/b = 12.5` in very stable air, so the surface stress does not collapse as fast as with
  Businger–Dyer's linear :math:`\phi`.

The unstable side is Businger–Dyer's with either choice. The selection applies to the wall model's
:math:`u_*` (``bc.surface_stability = heat_flux``), the ground gradients of the ghost cells, and the stratified
column (:ref:`stratified column <column-stability>`).

The wall model replaces the diffusive fluxes through the ground:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Field
     - Ground flux
   * - :math:`u, v`
     - :math:`u_*^2\,|A|\,\mathbf U_t/U`, along the ground wind
   * - :math:`\theta`
     - :math:`-q\,|A|`
   * - :math:`k`
     - 0; the first cell relaxes to :math:`k_w = (u_*^3 + \kappa B d_0)^{2/3}/C_{\mu 0}^2` (AL01 Eq. 16),
       with :math:`d_0` the first-cell-centre distance and :math:`B = gq/\theta_0` for :math:`q > 0`. With
       ``bc.surface_wall_k = dirichlet`` (default) this replaces the first cell's :math:`k` equation; with
       ``source`` it is added to it.

The ghost cells below the ground continue the MOST profile.

Inflow profile
--------------

The inflow profile gives the wind, :math:`\theta` and :math:`k` as functions of the height :math:`z` above the
ground. It is used on inflow faces, in the damping layers and for the initial field (:doc:`initial_conditions`). The
wind speed :math:`S(z)` points along ``bc.inflow_direction``:

.. math::

   \text{log\_law:}\quad S = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad
   \text{power\_law:}\quad S = S_\mathrm{ref}\left(\frac{z}{z_\mathrm{ref}}\right)^\alpha .

A ``file`` profile tabulates ``z u v theta [k]``. For the analytic profiles :math:`\theta` is the sounding
``bc.theta_profile``.

**Turbulence of the profile.** Without a :math:`k` column, :math:`k` and :math:`\nu_t` are in equilibrium with the
wind's shear: shear production equals dissipation in the neutral k–L closure (:doc:`turbulence_model`),

.. math::

   u_{*l} = l\,\frac{dS}{dz}, \qquad k = \frac{u_{*l}^2}{C_{\mu 0}^2}, \qquad \nu_t = l\,u_{*l}, \qquad
   \frac{1}{l} = \frac{1}{\kappa (z + z_0)} + \frac{1}{l_\mathrm{max}} .

- :math:`u_{*l}`: the local friction velocity; :math:`l`: the length; :math:`l_\mathrm{max}`: ``turb.max_geom_lscale``.
- ``bc.inflow_tke_length = closure`` (the default) uses this :math:`l`, the closure's own geometric length, so the
  inflow :math:`k` is the closure's equilibrium. ``unbounded`` drops the cap (:math:`l = \kappa(z + z_0)`); the log
  law then has :math:`u_{*l} = u_*` and :math:`k = u_*^2/C_{\mu 0}^2` at every height.
- **Why the cap.** Without it, :math:`\kappa z\,dS/dz` grows with height on a power law, and so does :math:`k`: on
  ERF's Askervein inflow it reaches about 9 m²/s² at 316 m, where ERF's flow has 0.31; with the cap, 0.31 at 316 m.
- A :math:`k` column is used as given, and ``bc.inflow_tke_length`` is refused (it would change nothing).
- The equilibrium is the neutral one, also above an inversion, and the height is the height above the ground (the
  closure's wall distance with ``turb.wall_distance = terrain_height``). The roughness is ``bc.inflow_z0``, which
  should equal ``bc.surface_z0``.
- A log-law wind with the capped length is not a momentum equilibrium (its stress falls with height as
  :math:`(l/\kappa(z + z_0))^2`); with the unbounded length the :math:`k` equation was out of balance instead. The
  domain adjusts either way; the capped form at least balances the :math:`k` equation.

Damping layers
--------------

Damping layers absorb waves before they reflect from the top or the sides. Both are off by default. They relax a
field :math:`q` towards the inflow profile :math:`q_\mathrm{ref}` (0 for :math:`w`) by adding a term to the
residual:

.. math::

   \text{Rayleigh (top):}\quad R \mathrel{+}= c_r \sin^2\!\left(\frac{\pi}{2}\,\frac{z - z_d}{z_\mathrm{top} - z_d}\right)
   (q - q_\mathrm{ref}),
   \qquad
   \text{sponge (sides):}\quad R \mathrel{+}= c_s\left(1 - \frac{s}{w_s}\right)^2 (q - q_\mathrm{ref}) .

- :math:`z_d`: start height of the Rayleigh layer (``bc.rayleigh_zdamp``); :math:`z_\mathrm{top}`: the top level of the
  mesh, which is not ``mesh.prob_hi`` z on a stretched grid (:doc:`terrain_coordinates`).
- :math:`c_r`: Rayleigh coefficient (``bc.rayleigh_dampcoef``); ``bc.rayleigh_fields`` selects the fields.
- :math:`s`: distance to an inflow or outflow face; :math:`w_s` is the sponge width (``bc.sponge_width``) and
  :math:`c_s` its strength (``bc.sponge_strength``).
- :math:`\theta` is relaxed to the sounding at the absolute height by default (``bc.rayleigh_height``).

In stratified flow over a hill, damp :math:`\theta` together with :math:`w`. Gravity waves exchange energy between
the two, and on ERF's 2-D hill damping :math:`w` alone let the pseudo-time iteration grow.

Pressure
--------

The pressure correction is periodic where the mesh is periodic, Dirichlet (:math:`p = 0`) on ``outflow`` faces and
Neumann on all other faces. The types follow ERF's ``get_fft_bc``.

Ghost cells
-----------

The solution has two layers of ghost cells outside each face. With interior values :math:`c_0, c_1` and ghost
values :math:`g_1, g_2` along the face normal:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Rule
     - Ghost values
   * - even
     - :math:`g_1 = c_0,\ g_2 = c_1`
   * - odd
     - :math:`g_1 = -c_0,\ g_2 = -c_1` (0 on the face)
   * - value :math:`V`
     - :math:`g_1 = 2V - c_0,\ g_2 = 2V - c_1`
   * - profile
     - :math:`g = P(h)` at the ghost's height above the ground
   * - gradient :math:`G`
     - the centred difference at :math:`c_0` equals :math:`G`

The directions are filled in the order x, y, z, so corner ghosts get values too.
