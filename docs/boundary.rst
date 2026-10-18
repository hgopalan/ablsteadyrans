Boundary conditions
===================

``Source/BoundaryConditions.{H,cpp}`` and ``Source/InflowProfile.{H,cpp}`` are the domain boundary conditions
(sub-phase 1.8): one type per domain face, the ghost-cell fills that express it, the residual rows and the rows of
the implicit side it replaces, the pressure boundary types the Poisson layer needs, and the optional Rayleigh and
sponge layers. The residual (:doc:`residual`) applies them when ``ResidualInputs::bc`` is set.

Face types
----------

Each of the six faces (``bc.xlo``, ``bc.xhi``, ``bc.ylo``, ``bc.yhi``, ``bc.zlo``, ``bc.zhi``) has one type. ERF's name
for it is in parentheses; the Poisson types follow ERF's ``get_fft_bc`` for that name.

.. list-table::
   :header-rows: 1
   :widths: 16 12 72

   * - Type
     - Faces
     - Meaning
   * - ``periodic`` (Periodic)
     - x, y
     - Exactly the directions with ``mesh.is_periodic = 1``; the images come from ``FillBoundary``.
   * - ``inflow`` (Inflow)
     - x, y
     - u, v, θ and k take the inflow profile (below) at their own height above the ground; w = 0 on the face
       (mirrored inside), p' is even.
   * - ``outflow`` (Outflow)
     - x, y
     - Zero gradient for every field; p' = 0 on the face.
   * - ``slip_wall`` (SlipWall)
     - x, y, ground, top
     - Normal velocity 0; the tangential velocities, θ and k even. At the top θ's gradient is ``bc.zhi_theta_grad``
       instead (a slip-wall ground keeps θ even).
   * - ``no_slip_wall`` (NoSlipWall)
     - x, y, ground, top
     - Velocity 0 on the face; θ, k and p' even. *Not in the plan*: added so the wall-bounded verification cases of
       later sub-phases have it.
   * - ``surface_layer`` (surface_layer)
     - ground
     - The wall model (:doc:`surface`): w = 0; u, v and θ continued below the ground with the wall law's gradient at
       the first cell centre; k and p' even.
   * - ``geostrophic`` (SlipWall)
     - top
     - u and v held at the geostrophic wind of the Coriolis options, w = 0; θ as ``slip_wall``. Needs Coriolis with
       geostrophic forcing.

The defaults are periodic sides, a ``surface_layer`` ground and a ``slip_wall`` top: ERF's canonical ABL cases.

Ghost fills
-----------

The state has two ghost layers. Along the normal of a face, with interior values :math:`c_0, c_1` (first and second
inside) and ghosts :math:`g_1, g_2`:

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Rule
     - Ghosts
   * - even
     - :math:`g_1 = c_0,\ g_2 = c_1` (zero gradient)
   * - odd
     - :math:`g_1 = -c_0,\ g_2 = -c_1` (the value is 0 on the face)
   * - value :math:`V`
     - :math:`g_1 = 2V - c_0,\ g_2 = 2V - c_1`
   * - profile
     - :math:`g = P(h)` at each ghost point's own height :math:`h` above the ground
   * - gradient :math:`G`
     - :math:`g_1 = c_1 - (z_1 - z_{g1})\,G,\ g_2 = g_1 - (z_{g1} - z_{g2})\,G`: the centred difference at
       :math:`c_0` equals :math:`G`, with :math:`z` the points' heights

A face-normal velocity has a valid point on the face: a value is set there and the ghosts mirror the interior about
it (even for outflow). The directions are filled x, then y, then z, each over the ghosts the earlier ones filled, so
every corner ghost gets a value; this is the order of the Poisson layer's ``fill_bc_ghosts``. The homogeneous fill,
used for a perturbation vector (``jvp``), sets :math:`V = P = G = 0`.

The heights of the ghost points come from the nodes: the lateral indices clamped to the stored nodes (one ghost
layer) and the vertical ones extrapolated linearly, so the second ghost layer continues the first layer's thickness.
The second lateral ghost column therefore takes the height of the first one (they differ only over laterally
varying terrain at the face, and only the far MUSCL ghost sees it). The gradient fill makes the difference in height
across the first point equal to :math:`G`; the operator's centred derivative there (:doc:`transport`, index-centred
with the metric at the point) equals :math:`G` on uniform layers and :math:`G\,(0.75 + 0.25\,r)` with a layer ratio
:math:`r` (1.25 % high at :math:`r = 1.05`).

**The ground.** With ``surface_layer`` the surface state is computed from the state after the lateral fills (it reads
the lateral ghost columns) and the ground gradient of each column is MOST's flux-gradient relation at the first cell
centre, :math:`u_*\,\phi_m(d_0/L)/(\kappa (d_0 + z_0))` along the wall stress for u and v, and
:math:`-q\,\phi_h(d_0/L)/(\kappa u_* (d_0 + z_0))` for θ (0 when neutral), with :math:`d_0` the normal distance of the
first cell centre, :math:`L` the wall model's Obukhov length and :math:`\phi_m, \phi_h` ERF's Businger–Dyer gradient
functions (1 when neutral). The closure's turbulent Prandtl number does not enter here: the ghost continues the
profile the wall model solved, and the closure's :math:`\mathrm{Pr}_t(R_t) = c_\mu/c_\mu'` is 1 in neutral air
(:doc:`closure`); the ground face flux itself is the wall model's. The face points of u and v average the two columns beside them. The gradient
slab carries two ghost columns like the state: the first from the surface state's own ghost column, the rest from
the neighbouring boxes and the periodic images, and outside a non-periodic face beyond that, the nearest column
inside the domain, so every ghost is finite without depending on the wall model outside the domain.

Residual rows and the implicit side
-----------------------------------

The face points of a Dirichlet normal velocity (inflow, the walls, the ground, the top) get the row
:math:`R = q - V` (:math:`V = 0` when homogeneous) and the identity row on the implicit side. The damping sources
(below) are added before these rows, after the wall-k rows of :doc:`surface`.

The first interior row along a face folds its ghost coupling :math:`a_g` of :math:`A_1` (:doc:`transport`):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Fill of the ghost
     - Fold
   * - value or odd (:math:`g = 2V - c_0` or :math:`-c_0`)
     - :math:`a_P \mathrel{+}= a_g`
   * - even (:math:`g = c_0`)
     - :math:`a_P \mathrel{-}= a_g`
   * - profile
     - the coupling is dropped (the ghost does not depend on the state)
   * - gradient (:math:`g_1 = c_1 + \dots`), and a normal velocity's even mirror about its face point
     - the coupling moves to the opposite neighbour

The ground's own diffusive coupling is the wall model's and is already gone (:doc:`surface`). The damping
coefficient :math:`c` (below) is added to the diagonal first, so the line solves see the relaxation and stay stable
for any pseudo-time step. UT-1.8.2b checks the folds against the homogeneous fill for every equation.

Pressure
--------

``poisson_bc_from_faces`` gives the Poisson layer's types as ERF's ``get_fft_bc``: periodic where the geometry is,
odd (Dirichlet, p' = 0) for ``outflow``, even (Neumann) for every other type. UT-1.8.1 compares with ERF's function
on every type.

Inflow profile
--------------

``InflowProfile`` gives the wind, θ and the turbulence that enter through an inflow face and that the damping layers
relax towards, as functions of the height :math:`z` above the ground. The wind is :math:`S(z)` along a fixed direction
:math:`a` (degrees from the x axis): :math:`(u, v) = S(\cos a, \sin a)`.

.. math::

   \text{log\_law:}\quad S = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad
   \text{power\_law:}\quad S = S_{\mathrm{ref}}\left(\frac{z}{z_{\mathrm{ref}}}\right)^\alpha .

A ``file`` profile tabulates ``z u v theta [k]`` per line, :math:`z` strictly increasing from 0, linear in between
and constant beyond the ends. θ of the analytic kinds is the reference profile ``bc.theta_profile`` (also the
hydrostatic reference of :doc:`residual`).

**Turbulence.** In the wall law the eddy viscosity is :math:`\nu_t = \kappa u_{*l}(z + z_0)` with the local friction
velocity :math:`u_{*l} = \kappa (z + z_0)\,dS/dz` (shear production equals dissipation), so

.. math::

   k = \frac{u_{*l}^2}{C_{\mu 0}^2}, \qquad \nu_t = \kappa\,u_{*l}\,(z + z_0) .

For the log law :math:`u_{*l} = u_*` at every height, so :math:`k` is the wall value :math:`u_*^2/C_{\mu 0}^2`. A file
with a k column uses it as given.

Damping
-------

Both layers are off by default. With :math:`q_{\mathrm{ref}}` the inflow profile's wind and θ (and k), and 0 for w:

.. math::

   \text{Rayleigh:}\quad R \mathrel{+}= c(z)\,(q - q_{\mathrm{ref}}),\quad
   c = c_r \sin^2\!\left(\frac{\pi}{2}\,\frac{z - z_d}{z_{\mathrm{top}} - z_d}\right)\ \text{for}\ z > z_d ,

.. math::

   \text{sponge:}\quad R \mathrel{+}= c(s)\,(q - q_{\mathrm{ref}}),\quad
   c = c_s\left(1 - \frac{s}{w}\right)^2\ \text{for}\ s < w ,

with :math:`s` the distance to an inflow or outflow face and :math:`w` the width. Both relax u, v, θ and k; the
reference is the inflow profile at the point's height above the local ground, so over a hill the Rayleigh layer
relaxes towards a profile lifted by the terrain height (ERF relaxes to its level-indexed sounding; the difference
is of order :math:`h/z` of the profile's gradient). Without an inflow face the reference is still ``bc.inflow_*``
(the log law with its defaults unless set). ERF's Rayleigh layer
(``rayleigh_zdamp``, ``rayleigh_dampcoef``) has the same :math:`\sin^2` shape, and ERF's sponge zones
(``ERF_ApplySpongeZoneBCs.cpp``) the same quadratic ramp, 1 on the face and 0 at the inner edge. Both are sources
inside ``residual()``, so ``jvp`` sees them.

Inputs
------

Every key is in ``bc``. A key with no effect for the chosen types stops the run. The surface-layer keys are read
here too (``bc.surface_*``, :doc:`surface`).

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and effect
   * - ``bc.xlo``, ``bc.xhi``, ``bc.ylo``, ``bc.yhi``
     - ``periodic``
     - ``periodic`` exactly when ``mesh.is_periodic`` is 1 in that direction; otherwise ``inflow``, ``outflow``,
       ``slip_wall`` or ``no_slip_wall``.
   * - ``bc.zlo``
     - ``surface_layer``
     - ``surface_layer``, ``slip_wall`` or ``no_slip_wall``.
   * - ``bc.zhi``
     - ``slip_wall``
     - ``slip_wall``, ``no_slip_wall`` or ``geostrophic`` (needs Coriolis with geostrophic forcing).
   * - ``bc.inflow_profile``
     - ``log_law``
     - ``log_law``, ``power_law`` or ``file``. The inflow keys need an inflow face or a damping layer.
   * - ``bc.inflow_ustar``
     - 0.4 m/s
     - Finite, > 0 (``log_law`` only; refused with the other kinds).
   * - ``bc.inflow_z0``
     - 0.1 m
     - Finite, > 0: the log law's roughness, and the roughness of the equilibrium k and :math:`\nu_t` of every kind.
       Normally equal to ``bc.surface_z0``.
   * - ``bc.inflow_speed``, ``bc.inflow_zref``, ``bc.inflow_alpha``
     - 10 m/s, 10 m, 0.14
     - Speed finite and ≥ 0; height and exponent finite and > 0 (``power_law`` only; refused with the other kinds).
       Below :math:`z_0` the law's slope at :math:`z_0` gives the equilibrium k (the slope itself diverges).
   * - ``bc.inflow_direction``
     - 0
     - Finite, degrees from the x axis.
   * - ``bc.inflow_file``
     - empty
     - Required with ``file`` (refused with the other kinds): ``z u v theta [k]`` rows (4 or 5 columns, at least 2),
       z increasing from 0, θ > 0, k ≥ 0.
   * - ``bc.theta_profile_z``, ``bc.theta_profile``
     - ``0``, ``300``
     - The same number of values, heights increasing from 0, θ > 0: the reference profile.
   * - ``bc.zhi_theta_grad``
     - 0
     - Finite (K/m): θ's gradient through a ``slip_wall`` or ``geostrophic`` top; no effect with ``no_slip_wall``.
   * - ``bc.rayleigh_zdamp``
     - 0 (off)
     - In [0, ``mesh.prob_hi`` z): the Rayleigh layer's start height; > 0 needs ``bc.rayleigh_dampcoef`` > 0.
   * - ``bc.rayleigh_dampcoef``
     - 0
     - Finite, ≥ 0 (1/s); needs ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.sponge_width``
     - 0 (off)
     - Finite, ≥ 0 (m); > 0 needs an inflow or outflow face and ``bc.sponge_strength`` > 0.
   * - ``bc.sponge_strength``
     - 0
     - Finite, ≥ 0 (1/s); needs ``bc.sponge_width`` > 0.
   * - ``bc.surface_z0``, ``bc.surface_zref``, ``bc.surface_stability``, ``bc.surface_heat_flux``, ``bc.surface_wall_k``
     - see :doc:`surface`
     - The surface layer's options (``read_surface_layer_options``); refused unless ``bc.zlo = surface_layer``.

The application reads the ``bc`` keys with the mass-consistent initializer (:doc:`initializer`, which checks the face
types against the mesh and takes the inflow profile and the sounding from them), and so does the regression tool
``boundary_fluxes`` (its deck also has the tool's own ``fixed_state.ustar``, ``fixed_state.z0`` and
``fixed_state.direction``); the SIMPLEC stage of sub-phase 1.14 applies them in the solver. Without the initializer, setting
them in an ablsteady deck still stops the run.
