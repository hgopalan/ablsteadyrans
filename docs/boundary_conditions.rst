Boundary conditions
===================

Each of the six domain faces has one type (``bc.xlo`` … ``bc.zhi``). The types are ERF's, and the pressure types
follow ERF's ``get_fft_bc``. The input table is in :doc:`inputs`. The code is ``Source/BoundaryConditions``,
``Source/InflowProfile`` and ``Source/SurfaceLayer``.

Face types
----------

.. list-table::
   :header-rows: 1
   :widths: 18 14 68

   * - Type
     - Faces
     - Meaning
   * - ``periodic``
     - x, y
     - The directions with ``mesh.is_periodic = 1``.
   * - ``inflow``
     - x, y
     - :math:`u, v, \theta, k` from the inflow profile at their height above the ground; normal velocity given;
       :math:`p` even.
   * - ``outflow``
     - x, y
     - Zero gradient for every field; :math:`p = 0` on the face.
   * - ``slip_wall``
     - all
     - No normal flux; tangential velocity, :math:`\theta` and :math:`k` even. At the top :math:`\theta` may take
       the gradient ``bc.zhi_theta_grad``.
   * - ``no_slip_wall``
     - all
     - Velocity equal to the wall's (0, or a moving lid ``bc.<face>_wall_velocity``).
   * - ``surface_layer``
     - ground
     - The wall model (below). The default ground.
   * - ``geostrophic``
     - top
     - :math:`u, v` held at the geostrophic wind.

The defaults are periodic sides, a ``surface_layer`` ground and a ``slip_wall`` top, ERF's canonical ABL setup.

**No flux through the ground and the top.** On a terrain-following mesh the normal velocity there is
:math:`\Omega`, not :math:`w`. The Dirichlet value of :math:`w` is therefore

.. math::

   w_\mathrm{face} = z_x\,\bar u + z_y\,\bar v ,

so :math:`\Omega = 0` on the face. Over flat ground it is 0.

Ghost cells
-----------

The state has two ghost layers. With interior values :math:`c_0, c_1` and ghosts :math:`g_1, g_2` along the normal:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Rule
     - Ghosts
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

Directions are filled x, then y, then z, so every corner ghost gets a value. A perturbation (for :math:`Jv`) uses the
same rules with :math:`V = P = G = 0`. The first interior row of :math:`A_1` folds its ghost coupling the same way
(added to :math:`a_P` for value/odd, subtracted for even, dropped for profile).

Wall model
----------

The ground uses Monin–Obukhov similarity (MOST) on every column:

.. math::

   u_* = \frac{\kappa\,U}{\max\left(\ln\dfrac{d + z_0}{z_0} - \psi_m(d/L),\; m\right)}, \qquad
   \frac1L = -\frac{\kappa\,g\,q}{u_*^3\,\theta_0} .

- :math:`U = \max(|\mathbf U_t|, 0.1)` m/s is the wind parallel to the ground at the reference height
  (``bc.surface_zref``; by default the first cell centre).
- :math:`d = h/|A|` is the normal distance, with :math:`|A| = \sqrt{1 + z_x^2 + z_y^2}`.
- :math:`z_0` is ``bc.surface_z0``, :math:`q` the kinematic heat flux ``bc.surface_heat_flux`` and :math:`\psi_m`
  ERF's Businger–Dyer function.
- ``bc.surface_stability = heat_flux`` solves :math:`u_*` by Newton's method to round-off.
- :math:`m` is the floor of the denominator. Strong convection at low wind would otherwise give :math:`u_* < 0`.
  ``bc.surface_denominator_floor = log_law`` (the default) uses :math:`m = \min(\ln((d + z_0)/z_0),\, 1)`;
  ``one`` uses :math:`m = 1`.

**Why the floor is** :math:`\min(\ln((d+z_0)/z_0), 1)`. For neutral flow (:math:`\psi_m = 0`) the floor then never
acts, and :math:`u_*` is the log law's at every height, as in the inflow profile and the column guess. The floor
:math:`m = 1` acts on neutral flow when the reference point is below :math:`(e-1)\,z_0 \approx 1.72\,z_0`, i.e. a
first layer thinner than :math:`3.44\,z_0`. Then :math:`u_*` is :math:`\ln((d+z_0)/z_0)` times the log law's value.
For example, :math:`z_0 = 0.5` m over a 1 m first layer gives 0.69 of the log law's :math:`u_*` and 0.48 of its
stress. Rough ground (:math:`z_0` = 0.1–0.5 m) with a first layer of 1–4 m is in this range. Where
:math:`\ln((d+z_0)/z_0) \ge 1` the two floors give the same answer, bit for bit. ``one`` keeps the form of
sub-phase 1.7.

**Why Newton.** ERF iterates a relaxed fixed point to 1e-5. That makes :math:`u_*` a step function of the wind,
which spoils a finite-difference Jacobian.

The wall model replaces the diffusive ground fluxes:

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
     - 0; the first cells relax to :math:`k_w = (u_*^3 + \kappa B d_0)^{2/3}/C_{\mu 0}^2`
       (AL01 Eq. 16, :math:`B = gq/\theta_0` for :math:`q > 0`)

The ghost cells below the ground continue the MOST profile, with gradient :math:`u_*\phi_m/(\kappa(d_0 + z_0))` for
the wind and :math:`-q\,\phi_h/(\kappa u_*(d_0 + z_0))` for :math:`\theta`.

Inflow profile
--------------

The inflow profile gives the wind, :math:`\theta` and :math:`k` at height :math:`z` above the ground. The wind is
:math:`S(z)` along the direction ``bc.inflow_direction``:

.. math::

   \text{log\_law:}\quad S = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad
   \text{power\_law:}\quad S = S_\mathrm{ref}\left(\frac{z}{z_\mathrm{ref}}\right)^\alpha .

A ``file`` profile tabulates ``z u v theta [k]``. Without a :math:`k` column, :math:`k = u_{*l}^2/C_{\mu 0}^2` with
the local friction velocity :math:`u_{*l} = \kappa(z + z_0)\,dS/dz` (production equals dissipation). :math:`\theta`
is the sounding ``bc.theta_profile`` for the analytic kinds.

Damping layers
--------------

Both are off by default. They relax a field :math:`q` to the inflow profile :math:`q_\mathrm{ref}` (0 for :math:`w`):

.. math::

   \text{Rayleigh:}\quad R \mathrel{+}= c_r \sin^2\!\left(\frac{\pi}{2}\,\frac{z - z_d}{z_\mathrm{top} - z_d}\right)
   (q - q_\mathrm{ref}),
   \qquad
   \text{sponge:}\quad R \mathrel{+}= c_s\left(1 - \frac{s}{w}\right)^2 (q - q_\mathrm{ref}) .

- :math:`z_d` is ``bc.rayleigh_zdamp`` and :math:`c_r` is ``bc.rayleigh_dampcoef``; ``bc.rayleigh_fields`` selects
  the fields.
- :math:`s` is the distance to an inflow or outflow face and :math:`w` is ``bc.sponge_width``.
- :math:`\theta` is relaxed to the sounding at the absolute height by default (``bc.rayleigh_height``).

**Why damp θ as well as w.** A stratified layer over a hill carries internal gravity waves. Damping :math:`w` alone
leaves the :math:`\theta` half of the oscillator free, and the pseudo-time iteration over ERF's 2-D hill grew.

Pressure
--------

The pressure correction is periodic where the mesh is, Dirichlet (:math:`p = 0`) on ``outflow`` faces and Neumann
on every other face.
