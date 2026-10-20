Input reference
===============

This page lists every input of ablsteady. The first table is a quick reference; the sections after it give the
allowed range and meaning of each input, grouped by namespace.

Input file format
-----------------

An input file has one input per line, ``namespace.key = value``. Lists are separated by spaces, and ``#`` starts a
comment:

.. code-block:: text

   mesh.n_cell  = 64 64 32     # cells in x, y, z
   bc.zlo       = surface_layer

Inputs can also be given on the command line after the input file (:doc:`running`).

The run stops before doing any work, naming the input, when:

- a key is misspelled or not used by this version;
- a key's namespace is unknown (the AMReX namespaces ``amrex``, ``fab``, ``fabarray``, ``vismf``, ``tiny_profiler``,
  ``DistributionMapping``, ``device`` and ``hypre`` are allowed);
- a key is given twice, or a single-value key is given several values;
- a value is outside its range;
- a key has no effect with the other inputs (for example ``terrain.hill_height`` with ``terrain.type = flat``).

All errors are reported together.

Quick reference
---------------

.. list-table::
   :header-rows: 1
   :widths: 34 20 46

   * - Input
     - Default
     - Description
   * - **Mesh** (`mesh`_)
     -
     -
   * - ``mesh.n_cell``
     - required
     - Number of cells in x, y, z
   * - ``mesh.prob_lo``
     - ``0 0 0``
     - Low corner of the domain (m)
   * - ``mesh.prob_hi``
     - required
     - High corner of the domain (m)
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Periodic directions
   * - ``mesh.max_grid_size``
     - ``64``
     - Largest box size in x and y
   * - ``mesh.dz0``
     - ``0``
     - First layer thickness of a stretched grid (m)
   * - ``mesh.stretch_ratio``
     - none
     - Fixed ratio between successive layers
   * - ``mesh.zlevels_file``
     - none
     - File with the vertical levels
   * - ``mesh.hierarchy_levels``
     - ``1``
     - Number of grid-hierarchy levels
   * - **Terrain** (`terrain`_)
     -
     -
   * - ``terrain.type``
     - ``flat``
     - Terrain source
   * - ``terrain.file``
     - none
     - Terrain file
   * - ``terrain.smoothing``
     - ``1``
     - Terrain-following mapping (0 BTF, 1 STF, 2 Sullivan)
   * - ``terrain.gamma_m``
     - ``0.5``
     - Smallest layer thickness under STF, relative to flat
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Hill height (m)
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Hill half-width (m)
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - Ridge or round hill
   * - ``terrain.idw_neighbors``
     - ``6``
     - Points used per height for ``xyz`` terrain
   * - ``terrain.idw_power``
     - ``2``
     - Inverse-distance weight exponent
   * - ``terrain.hierarchy_passes``
     - ``1``
     - Terrain smoothing on coarse hierarchy levels
   * - **Initial conditions** (`init`_)
     -
     -
   * - ``init.first_guess``
     - ``none``
     - Initial profile: ``none``, ``profile`` or ``column``
   * - ``init.projection``
     - ``1``
     - Mass-consistent adjustment on or off
   * - ``init.alpha_h``, ``init.alpha_v``
     - ``1``, ``1``
     - Horizontal and vertical transmissivity
   * - ``init.z_agl_min``
     - ``0``
     - Lowest height at which the profile is evaluated (m)
   * - ``init.k_decay_height``
     - ``0``
     - Height above which the initial TKE decays (m)
   * - ``init.obrien``
     - ``0``
     - O'Brien vertical-velocity correction
   * - ``init.column_tol``
     - ``1e-8``
     - Convergence tolerance of the single-column solve
   * - ``init.column_max_iterations``
     - ``5000``
     - Step limit of the single-column solve
   * - **Boundary conditions** (`bc`_)
     -
     -
   * - ``bc.xlo`` … ``bc.yhi``
     - ``periodic``
     - Side boundary types
   * - ``bc.zlo``
     - ``surface_layer``
     - Ground boundary type
   * - ``bc.zhi``
     - ``slip_wall``
     - Top boundary type
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Velocity of a no-slip wall (m/s)
   * - ``bc.inflow_profile``
     - ``log_law``
     - Inflow profile type
   * - ``bc.inflow_ustar``, ``bc.inflow_z0``
     - 0.4 m/s, 0.1 m
     - Log-law friction velocity and roughness
   * - ``bc.inflow_speed``, ``bc.inflow_zref``, ``bc.inflow_alpha``
     - 10 m/s, 10 m, 0.14
     - Power-law speed, reference height and exponent
   * - ``bc.inflow_direction``
     - ``0``
     - Wind direction (degrees from the x axis)
   * - ``bc.inflow_file``
     - none
     - Tabulated inflow profile
   * - ``bc.theta_profile_z``, ``bc.theta_profile``
     - ``0``, ``300``
     - Potential-temperature sounding (m, K)
   * - ``bc.zhi_theta_grad``
     - ``0``
     - θ gradient through the top (K/m)
   * - ``bc.rayleigh_zdamp``, ``bc.rayleigh_dampcoef``
     - ``0``, ``0``
     - Rayleigh layer start height (m) and strength (1/s)
   * - ``bc.rayleigh_fields``
     - ``u v w theta tke``
     - Fields the Rayleigh layer damps
   * - ``bc.rayleigh_height``
     - ``absolute``
     - Height used by the Rayleigh layer
   * - ``bc.sponge_width``, ``bc.sponge_strength``
     - ``0``, ``0``
     - Side sponge width (m) and strength (1/s)
   * - ``bc.surface_z0``
     - 0.1 m
     - Ground roughness length
   * - ``bc.surface_zref``
     - ``0``
     - Wall-model reference height (0: first cell centre)
   * - ``bc.surface_stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``
   * - ``bc.surface_heat_flux``
     - ``0``
     - Kinematic surface heat flux (K m/s)
   * - ``bc.surface_wall_k``
     - ``dirichlet``
     - Wall treatment of the TKE
   * - ``bc.surface_denominator_floor``
     - ``log_law``
     - Floor of the MOST denominator
   * - **Forcing** (`forcing`_)
     -
     -
   * - ``forcing.gravity``
     - ``9.81``
     - Gravitational acceleration (m/s²)
   * - ``forcing.theta0``
     - ``300``
     - Boussinesq reference temperature (K)
   * - ``forcing.coriolis``
     - ``0``
     - Coriolis force on or off
   * - ``forcing.latitude``
     - ``0``
     - Latitude (degrees)
   * - ``forcing.rotation_period``
     - ``86400``
     - Earth's rotation period (s)
   * - ``forcing.coriolis_3d``
     - ``1``
     - Include the cos(latitude) Coriolis terms
   * - ``forcing.geostrophic``
     - ``0``
     - Geostrophic forcing on or off
   * - ``forcing.ug``, ``forcing.vg``
     - ``0``, ``0``
     - Geostrophic wind (m/s)
   * - ``forcing.pressure_gradient``
     - ``0 0 0``
     - Constant pressure-gradient force (m/s²)
   * - **Turbulence model** (`turb`_)
     -
     -
   * - ``turb.Cmu0``
     - ``0.5562``
     - Model constant :math:`C_{\mu 0}`
   * - ``turb.Cb``
     - ``0.35``
     - Model constant :math:`C_b`
   * - ``turb.Rt_crit``, ``turb.Rt_min``
     - ``-1``, ``-3``
     - Richardson-number smoothing limits
   * - ``turb.max_geom_lscale``
     - 30 m
     - Cap on the geometric length scale
   * - ``turb.tke_floor``
     - ``0``
     - Lower bound on k (0: machine epsilon)
   * - ``turb.sigma_k``
     - ``1``
     - TKE diffusion Schmidt number
   * - ``turb.nu_molecular``, ``turb.alpha_molecular``
     - ``0``, ``0``
     - Molecular viscosity and heat diffusivity (m²/s)
   * - ``turb.kappa``
     - ``0.41``
     - von Kármán constant
   * - ``turb.nut_relax``
     - ``1``
     - Under-relaxation of :math:`\nu_t`
   * - ``turb.frozen_switch``
     - ``0``
     - Freeze the model's branch choices for a Jacobian
   * - ``turb.face_mean``
     - ``arithmetic``
     - Face average of the diffusivities
   * - ``turb.wall_distance``
     - ``terrain_height``
     - Wall-distance method
   * - **Convection** (`convection`_)
     -
     -
   * - ``convection.limiter``
     - ``van_leer``
     - Slope limiter for momentum and θ
   * - ``convection.smooth_eps_velocity``, ``convection.smooth_eps_theta``
     - ``1e-5``, ``1e-5``
     - Smoothing constants of ``van_albada``
   * - **Steady solver** (`steady`_)
     -
     -
   * - ``steady.stage``
     - ``none``
     - ``none`` (build and stop) or ``simplec``
   * - ``steady.max_iterations``
     - ``10000``
     - Iteration limit
   * - ``steady.residual_tol``
     - ``5e-5``
     - Absolute residual tolerance
   * - ``steady.residual_drop_orders``
     - ``4``
     - Required residual drop (orders of magnitude)
   * - ``steady.residual_ref_iterations``
     - ``5``
     - Iterations over which the reference residual is taken
   * - ``steady.monitor_rel_tol``, ``steady.monitor_window``
     - ``1e-3``, ``50``
     - Monitor tolerance and window
   * - ``steady.cfl0``, ``steady.cfl_max``
     - ``1``, ``10``
     - Initial and largest CFL number
   * - ``steady.ser_exponent``, ``steady.ser_max_growth``
     - ``1``, ``2``
     - CFL ramp exponent and largest growth per iteration
   * - ``steady.switch_window``, ``steady.switch_residual``, ``steady.switch_back_growth``
     - ``20``, ``1e-2``, ``10``
     - Stage-switching rules
   * - ``steady.sweep``
     - ``zjacobi``
     - Line solver: ``zjacobi`` or ``lusgs``
   * - ``steady.zline_sweeps``
     - ``2``
     - Sweeps per iteration
   * - ``steady.line_ghosts``, ``steady.line_overlap``
     - ``lagged``, ``0``
     - LU-SGS coupling across box faces
   * - ``steady.line_order``
     - ``sequential``
     - LU-SGS visiting order
   * - ``steady.line_alternate``
     - ``0``
     - Alternate the LU-SGS sweep direction
   * - ``steady.sequencing``
     - ``0``
     - Grid sequencing, coarse to fine
   * - ``steady.sequencing_drop_orders``
     - ``2``
     - Residual drop for the loose stop of a coarse level
   * - ``steady.sequencing_monitor_tol``
     - ``0``
     - Optional monitor test of the loose stop (0: off)
   * - ``steady.sequencing_hold``
     - ``0``
     - Keep the reference hold on the seeded levels
   * - ``steady.sequencing_cfl_start``
     - ``cfl0``
     - CFL of a seeded level's first iteration
   * - ``steady.machine``, ``steady.nodes``
     - ``local``, ``1``
     - Machine and node count for the allocation accounting
   * - ``steady.newton_*``
     - see `steady`_
     - Newton–Krylov inputs (not yet used by the application)
   * - ``steady.newton_preconditioner``
     - ``simplec``
     - Newton preconditioner: ``identity`` or ``simplec`` (not yet used by the application)
   * - **SIMPLEC** (`simplec`_)
     -
     -
   * - ``simplec.pressure_reltol``
     - ``1e-3``
     - Tolerance of the pressure-correction solve
   * - ``simplec.alpha_p``
     - ``1``
     - Pressure under-relaxation
   * - ``simplec.rotation_relax``
     - ``0.3``
     - Point-implicit rotation step
   * - ``simplec.w_vnn_cap``
     - ``0``
     - Vertical-velocity cap (0: off)
   * - ``simplec.w_vnn_cap_slope_margin``
     - ``4``
     - The cap's safety margin on sloping ground (0: off)
   * - **Poisson solver** (`solver`_)
     -
     -
   * - ``solver.poisson_backend``
     - ``gmres_mlmg``
     - Poisson solver
   * - ``solver.poisson_reltol``, ``solver.poisson_abstol``
     - ``1e-10``, ``0``
     - Relative and absolute tolerance
   * - ``solver.poisson_max_iter``
     - ``200``
     - Iteration limit
   * - ``solver.poisson_restart``
     - ``50``
     - GMRES restart length
   * - ``solver.poisson_precond``
     - ``fft``
     - Preconditioner of ``gmres_fft``
   * - ``solver.poisson_verbose``
     - ``0``
     - Output level
   * - **Input/output** (`io`_)
     -
     -
   * - ``io.plot_file``, ``io.chk_file``
     - ``plt``, ``chk``
     - Plotfile and checkpoint prefixes
   * - ``io.write_plotfile``, ``io.write_checkpoint``
     - ``1``, ``1``
     - Write the plotfile and checkpoint
   * - ``io.turb_plot_file``
     - none
     - Prefix of the turbulence plotfile
   * - ``io.restart``
     - none
     - Checkpoint to restart from

``mesh``
--------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``mesh.n_cell``
     - required
     - Three integers ≥ 1: cells in x, y and z.
   * - ``mesh.prob_lo``
     - ``0 0 0``
     - Three finite values: the low corner of the domain (m). The z value must be 0.
   * - ``mesh.prob_hi``
     - required
     - Three finite values, each greater than ``prob_lo``: the high corner of the domain (m).
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Three values, each 0 or 1. z must be 0.
   * - ``mesh.max_grid_size``
     - ``64``
     - One or two integers ≥ 1: the largest box size in x and y. Boxes are never split in z.
   * - ``mesh.dz0``
     - ``0``
     - Finite, ≥ 0 (m). First layer thickness of a stretched grid; 0 gives uniform layers. Without
       ``mesh.stretch_ratio`` it must not exceed ``prob_hi`` z / ``n_cell`` z (:doc:`terrain_coordinates`).
   * - ``mesh.stretch_ratio``
     - none
     - Finite, ≥ 1; needs ``mesh.dz0`` > 0. A fixed ratio between successive layers; the top is then where the
       layers end, and ``prob_hi`` z only sets the computational spacing.
   * - ``mesh.zlevels_file``
     - none
     - A file of ``n_cell`` z + 1 increasing heights starting at 0. Not together with ``mesh.dz0``.
   * - ``mesh.hierarchy_levels``
     - ``1``
     - Integer from 1 to 10: the levels :math:`L` of the grid hierarchy (:doc:`grid_hierarchy`). With
       :math:`L > 1`, ``n_cell`` and ``max_grid_size`` in x and y must be divisible by :math:`2^{L-1}`, and the
       coarsest level must be larger than 1 × 1 cells. The box layout of level 0 can then differ from :math:`L = 1`,
       which changes results in the last bits. The solver uses the coarse levels only with ``steady.sequencing =
       1``. May change across a restart.

``terrain``
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``terrain.type``
     - ``flat``
     - ``flat``, ``agnesi`` (analytic hill), ``xyz`` (point cloud) or ``erf_grid`` (ERF's terrain file)
       (:doc:`terrain_coordinates`).
   * - ``terrain.file``
     - required for ``xyz``, ``erf_grid``
     - The point file (``x y z`` per line) or ERF's structured terrain file.
   * - ``terrain.smoothing``
     - ``1``
     - 0 (BTF), 1 (STF) or 2 (Sullivan).
   * - ``terrain.gamma_m``
     - ``0.5``
     - > 0 and < 1; STF only. Smallest layer thickness as a fraction of the flat layer.
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Finite, > 0 (m).
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Finite, > 0 (m). Distance from the crest at which the hill is half as high.
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - ``ridge_x`` (a ridge along y, varying in x), ``ridge_y`` or ``round``; ``agnesi`` only.
   * - ``terrain.idw_neighbors``
     - ``6``
     - ≥ 1; ``xyz`` only. Points used for each height.
   * - ``terrain.idw_power``
     - ``2``
     - Finite, > 0; ``xyz`` only. Weight exponent :math:`p` in :math:`d^{-p}`.
   * - ``terrain.hierarchy_passes``
     - ``1``
     - Integer from 0 to 64; only with ``mesh.hierarchy_levels`` > 1. The ground of level :math:`l` is smoothed by
       this number times :math:`2^l - 1` passes (:doc:`grid_hierarchy`).

``init``
--------

See :doc:`initial_conditions`. The keys other than ``init.first_guess`` are refused with ``init.first_guess = none``.
A restart cannot be combined with ``init.*`` keys.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``init.first_guess``
     - ``none``
     - ``none`` (fluid at rest), ``profile`` (the inflow profile by height above the ground) or ``column`` (the
       steady single-column profile; needs a horizontal forcing and the SIMPLEC stage, see
       :doc:`initial_conditions` for the full list).
   * - ``init.projection``
     - ``1``
     - 0 or 1: the mass-consistent adjustment. With 0, ``init.alpha_h`` and ``init.alpha_v`` are refused.
   * - ``init.alpha_h``
     - ``1``
     - Finite, > 0. Horizontal transmissivity :math:`\alpha_h`.
   * - ``init.alpha_v``
     - ``1``
     - Finite, > 0. Vertical transmissivity :math:`\alpha_v`. Only :math:`\alpha_v/\alpha_h` matters. Ratios below
       0.1 need more than the default ``solver.poisson_max_iter``. With :math:`\alpha_v \ne \alpha_h` the
       adjustment uses the ``gmres_fft`` back end.
   * - ``init.z_agl_min``
     - ``0``
     - Finite, ≥ 0 (m). Lowest height above the ground at which the profile is evaluated.
   * - ``init.k_decay_height``
     - ``0``
     - Finite, ≥ 0 (m). Height :math:`z_i` above which the initial :math:`k` decays; 0: no decay. Not with
       ``column``.
   * - ``init.obrien``
     - ``0``
     - 0 or 1. The O'Brien correction of the vertical velocity before the adjustment.
   * - ``init.column_tol``
     - ``1e-8``
     - Finite, > 0; ``column`` only. Convergence tolerance of the single-column solve: the largest relative change
       of :math:`u, v, k` over a step. Raised to 100 machine epsilon in single precision (1.2e-5).
   * - ``init.column_max_iterations``
     - ``5000``
     - ≥ 1; ``column`` only. The run stops if the column has not converged after this many steps.

``bc``
------

See :doc:`boundary_conditions`. Read when ``init.first_guess`` is not ``none`` or with the SIMPLEC stage; refused
otherwise. A key that has no effect for the chosen boundary types is refused. The initial field uses the boundary
types, the inflow profile and the sounding; the other keys act in the SIMPLEC stage.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
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
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Three finite values (m/s): the velocity of a ``no_slip_wall`` face (``<face>`` is ``xlo`` … ``zhi``). The
       normal component must be 0.
   * - ``bc.inflow_profile``
     - ``log_law``
     - ``log_law``, ``power_law`` or ``file``. The inflow keys need an inflow face, a damping layer,
       ``init.first_guess`` or a restart.
   * - ``bc.inflow_ustar``
     - 0.4 m/s
     - Finite, > 0; ``log_law`` only.
   * - ``bc.inflow_z0``
     - 0.1 m
     - Finite, > 0. Roughness of the log law, and of the equilibrium :math:`k` of every profile type. Normally equal
       to ``bc.surface_z0``.
   * - ``bc.inflow_speed``, ``bc.inflow_zref``, ``bc.inflow_alpha``
     - 10 m/s, 10 m, 0.14
     - ``power_law`` only. Speed finite and ≥ 0; height and exponent finite and > 0.
   * - ``bc.inflow_direction``
     - ``0``
     - Finite (degrees). The direction the wind blows towards, counter-clockwise from the x axis. ``log_law`` and
       ``power_law`` only; a ``file`` profile gives u and v itself.
   * - ``bc.inflow_file``
     - none
     - Required with ``file``: rows of ``z u v theta [k]`` (4 or 5 columns, at least 2 rows), z increasing from
       0, θ > 0, k ≥ 0.
   * - ``bc.theta_profile_z``, ``bc.theta_profile``
     - ``0``, ``300``
     - Lists of equal length: heights increasing from 0 (m) and θ > 0 (K). The reference sounding.
   * - ``bc.zhi_theta_grad``
     - ``0``
     - Finite (K/m): θ gradient through a ``slip_wall`` or ``geostrophic`` top.
   * - ``bc.rayleigh_zdamp``
     - ``0`` (off)
     - In [0, ``mesh.prob_hi`` z): start height of the Rayleigh layer (m); > 0 needs ``bc.rayleigh_dampcoef`` > 0.
   * - ``bc.rayleigh_dampcoef``
     - ``0``
     - Finite, ≥ 0 (1/s); needs ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.rayleigh_fields``
     - ``u v w theta tke``
     - One or more of ``u v w theta tke``; needs ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.rayleigh_height``
     - ``absolute``
     - ``absolute`` or ``above_ground``: the height used for the layer's coefficient and for θ's reference
       (:math:`u, v, k` always use the height above the ground).
   * - ``bc.sponge_width``
     - ``0`` (off)
     - Finite, ≥ 0 (m); > 0 needs an inflow or outflow face and ``bc.sponge_strength`` > 0.
   * - ``bc.sponge_strength``
     - ``0``
     - Finite, ≥ 0 (1/s); needs ``bc.sponge_width`` > 0.
   * - ``bc.surface_z0``
     - 0.1 m
     - > 0. Roughness length of the wall model. The ``bc.surface_*`` keys need ``bc.zlo = surface_layer``.
   * - ``bc.surface_zref``
     - ``0``
     - ≤ 0: each column's first cell centre. Otherwise a height between the first and last cell-centre heights of
       every column (m).
   * - ``bc.surface_stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``.
   * - ``bc.surface_heat_flux``
     - ``0``
     - Finite (K m/s); ``heat_flux`` only.
   * - ``bc.surface_wall_k``
     - ``dirichlet``
     - ``dirichlet`` or ``source``.
   * - ``bc.surface_denominator_floor``
     - ``log_law``
     - ``log_law`` (:math:`m = \min(\ln((d + z_0)/z_0), 1)`) or ``one`` (:math:`m = 1`): the floor of the MOST
       denominator. They differ only where :math:`\ln((d + z_0)/z_0) < 1`.

The wall model and the inflow profile's equilibrium :math:`k` and :math:`\nu_t` use the turbulence model's
``turb.Cmu0`` and ``turb.kappa``. Without the SIMPLEC stage (``turb.*`` is then refused) they use the defaults, 0.5562
and 0.41.

``forcing``
-----------

See :doc:`governing_equations`. Read with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``forcing.gravity``
     - ``9.81``
     - Finite, ≥ 0 (m/s²). 0 switches buoyancy off.
   * - ``forcing.theta0``
     - ``300``
     - Finite, > 0 (K). Boussinesq reference temperature.
   * - ``forcing.coriolis``
     - ``0``
     - 0 or 1: the Coriolis force.
   * - ``forcing.latitude``
     - ``0``
     - In [−90, 90] (degrees). With ``coriolis = 1`` only.
   * - ``forcing.rotation_period``
     - ``86400``
     - Finite, > 0 (s); :math:`f = 4\pi\sin\varphi/T`. With ``coriolis = 1`` only.
   * - ``forcing.coriolis_3d``
     - ``1``
     - 0 or 1: the :math:`\cos\varphi` terms on :math:`u` and :math:`w`. With ``coriolis = 1`` only.
   * - ``forcing.geostrophic``
     - ``0``
     - 0 or 1: geostrophic forcing. Needs ``coriolis = 1`` and a non-zero latitude.
   * - ``forcing.ug``, ``forcing.vg``
     - ``0``
     - Finite (m/s): the geostrophic wind. With ``geostrophic = 1`` only.
   * - ``forcing.pressure_gradient``
     - ``0 0 0``
     - Three finite values (m/s²): a constant pressure-gradient force on :math:`u, v, w`.

``turb``
--------

See :doc:`turbulence_model`. Read with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``turb.Cmu0``
     - ``0.5562``
     - Finite, > 0: :math:`C_{\mu 0}`; also used by the wall k and the inflow profile's equilibrium :math:`k`.
   * - ``turb.Cb``
     - ``0.35``
     - Finite, > 0: :math:`C_b` of the stable and unstable length scales.
   * - ``turb.Rt_crit``
     - ``-1``
     - Finite, ≤ 0: :math:`R_t` is smoothed below this value.
   * - ``turb.Rt_min``
     - ``-3``
     - Finite, < ``turb.Rt_crit`` and > −3.6: the smoothed :math:`R_t` stays above it.
   * - ``turb.max_geom_lscale``
     - 30 m
     - Finite, > 0: the harmonic cap of :math:`l_g` (about :math:`\kappa \times 0.1` of the boundary-layer depth).
   * - ``turb.tke_floor``
     - ``0``
     - 0 (machine epsilon) or ≥ machine epsilon (m²/s²): lower bound on :math:`k`.
   * - ``turb.sigma_k``
     - ``1``
     - Finite, > 0: :math:`\gamma_k = \nu_\mathrm{mol} + \nu_t/\sigma_k`.
   * - ``turb.nu_molecular``, ``turb.alpha_molecular``
     - ``0``, ``0``
     - Finite, ≥ 0 (m²/s): molecular viscosity and heat diffusivity.
   * - ``turb.kappa``
     - ``0.41``
     - Finite, > 0: von Kármán constant; also used by the wall law and the inflow profile's :math:`\nu_t`.
   * - ``turb.nut_relax``
     - ``1``
     - In (0, 1]: under-relaxation of :math:`\nu_t` (1: none).
   * - ``turb.frozen_switch``
     - ``0``
     - 0 or 1: freeze the stable/unstable branch choices before a Jacobian is formed.
   * - ``turb.face_mean``
     - ``arithmetic``
     - ``arithmetic`` or ``logarithmic`` (:doc:`discretization`).
   * - ``turb.wall_distance``
     - ``terrain_height``
     - ``terrain_height`` or ``tucker``.

``convection``
--------------

See :doc:`discretization`. Read with the SIMPLEC stage only. :math:`k` always uses first-order upwinding.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``convection.limiter``
     - ``van_leer``
     - ``van_leer``, ``van_albada`` or ``upwind``, for momentum and θ.
   * - ``convection.smooth_eps_velocity``
     - ``1e-5``
     - > 0 (m/s), with a finite square: the :math:`\epsilon` of ``van_albada`` for momentum.
   * - ``convection.smooth_eps_theta``
     - ``1e-5``
     - The same for θ (K).

``steady``
----------

See :doc:`solution_algorithm`. With ``steady.stage = none`` the application builds the mesh and initial field,
writes them and stops; the other ``steady.*`` keys are then refused.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Input
     - Default
     - Range and meaning
   * - ``steady.stage``
     - ``none``
     - ``none`` or ``simplec``.
   * - ``steady.max_iterations``
     - ``10000``
     - ≥ 1. Iteration limit over the whole run, restarts included. With grid sequencing it applies to each level.
   * - ``steady.residual_tol``
     - ``5e-5``
     - Finite, > 0. Every scaled residual must be below it.
   * - ``steady.residual_drop_orders``
     - ``4``
     - ≥ 0 and ≤ the largest decimal exponent of the precision (308 double, 38 single). Every residual must also
       drop this many orders below its reference.
   * - ``steady.residual_ref_iterations``
     - ``5``
     - ≥ 1. The reference is the largest residual over this many first iterations. During them the CFL stays at
       ``steady.cfl0`` (for grid sequencing see ``steady.sequencing_hold``).
   * - ``steady.monitor_rel_tol``
     - ``1e-3``
     - Finite, > 0. Relative tolerance of the monitors.
   * - ``steady.monitor_window``
     - ``50``
     - ≥ 2. Monitor window (iterations).
   * - ``steady.cfl0``
     - ``1``
     - Finite, > 0. CFL number of the first iterations.
   * - ``steady.cfl_max``
     - ``10``
     - Finite, ≥ ``steady.cfl0``. Largest CFL number.
   * - ``steady.ser_exponent``
     - ``1``
     - Finite, ≥ 0. Exponent :math:`p` of the CFL ramp; 0 holds the CFL at ``cfl0``.
   * - ``steady.ser_max_growth``
     - ``2``
     - Finite, ≥ 1. Largest CFL growth factor per iteration.
   * - ``steady.switch_window``
     - ``20``
     - ≥ 1. Iterations below ``switch_residual`` before a switch to the next stage.
   * - ``steady.switch_residual``
     - ``1e-2``
     - > 0 and ≤ 1. Relative residual at which the next stage starts.
   * - ``steady.switch_back_growth``
     - ``10``
     - Finite, > 1. Residual growth over a window that returns the run to the previous stage.
   * - ``steady.sweep``
     - ``zjacobi``
     - ``zjacobi`` (z-line Jacobi) or ``lusgs`` (line LU-SGS).
   * - ``steady.zline_sweeps``
     - ``2``
     - Sweeps per iteration: 2 to 100 with ``zjacobi``, 1 to 100 with ``lusgs``. With ``zjacobi`` the stable CFL
       grows with the sweeps: 2.4 with 4, 3.7 with 8, 6.1 with 16.
   * - ``steady.line_ghosts``
     - ``lagged``
     - ``lagged`` or ``zero``; ``lusgs`` only. ``lagged``: neighbours across a box face take the previous pass's
       values. ``zero``: each box sweeps on its own, without communication.
   * - ``steady.line_overlap``
     - ``0``
     - Integer from 0 to 8, and at most the narrowest box width (on every hierarchy level with grid sequencing);
       > 0 only with ``line_ghosts = zero``. Columns by which each box's swept region is grown.
   * - ``steady.line_order``
     - ``sequential``
     - ``sequential`` or ``diagonal``; ``lusgs`` only. ``diagonal`` runs one column per GPU thread with the same
       result.
   * - ``steady.line_alternate``
     - ``0``
     - 0 or 1; ``lusgs`` only, with ``zline_sweeps`` ≥ 2. 1 reverses the y direction on every second sweep.
   * - ``steady.sequencing``
     - ``0``
     - 0 or 1. 1: grid sequencing over the ``mesh.hierarchy_levels`` levels (:doc:`solution_algorithm`). Needs
       ``mesh.hierarchy_levels`` ≥ 2 and at least 2 cells in x and y on the coarsest level. ``init.*`` then acts on
       the coarsest level only, and ``steady.max_iterations`` applies to each level. On a restart the coarse levels
       are skipped and the ``steady.sequencing*`` keys have no effect (the log says so).
   * - ``steady.sequencing_drop_orders``
     - ``2``
     - Finite, > 0 and ≤ ``steady.residual_drop_orders``; with ``sequencing = 1`` only. A coarse level stops once
       every residual is this many orders below its reference, or below ``steady.residual_tol``.
   * - ``steady.sequencing_monitor_tol``
     - ``0`` (off)
     - 0, or finite, > 0 and < 1; with ``sequencing = 1`` only. A coarse level also stops once every monitor changed
       by at most this fraction over ``steady.monitor_window`` iterations.
   * - ``steady.sequencing_hold``
     - ``0``
     - 0 or 1; with ``sequencing = 1`` only. 0: the seeded levels (all but the coarsest) skip the reference hold and
       ramp the CFL from their first iteration. 1: every level holds the CFL at ``steady.cfl0`` for its first
       ``steady.residual_ref_iterations`` iterations. A restart keeps the setting of the checkpoint.
   * - ``steady.sequencing_cfl_start``
     - ``cfl0``
     - ``cfl0`` or ``previous``; with ``sequencing = 1`` only. The CFL of a seeded level's first iteration:
       ``steady.cfl0``, or the CFL the previous level reached (at least ``steady.cfl0``). ``previous`` needs
       ``steady.sequencing_hold = 0``.
   * - ``steady.machine``
     - ``local``
     - ``local``, ``kestrel-cpu`` or ``kestrel-gpu``: allocation-unit rate (0, 10, 100 AU per node-hour).
   * - ``steady.nodes``
     - ``1``
     - ≥ 1. Nodes of the run, for the allocation units.
   * - ``steady.newton_gmres_restart``
     - ``30``
     - Integer from 5 to 200. GMRES restart length of the Newton step.
   * - ``steady.newton_fd_eps_scale``
     - :math:`\sqrt{\epsilon_\mathrm{mach}}`
     - ≥ :math:`\epsilon_\mathrm{mach}^{3/4}` and ≤ 1e-2. Relative step of the Jacobian–vector product.
   * - ``steady.newton_krylov_rtol``
     - ``1e-2``
     - > 0 and < 1. GMRES tolerance in each Newton step.
   * - ``steady.newton_mem_budget_gb``
     - ``70``
     - > 0 and ≤ 80. Krylov memory per rank (GB); the run stops if the estimate exceeds it.
   * - ``steady.newton_gmres_max_iters``
     - ``2000``
     - Integer, ≥ 1. GMRES iteration limit of one Newton step (all restarts). A solve that reaches it is unconverged;
       its update goes through the usual usability rule.
   * - ``steady.newton_preconditioner``
     - ``simplec``
     - ``identity`` or ``simplec``; any other value aborts and names the key. The preconditioner of each GMRES
       solve of the Newton step (:doc:`solution_algorithm`). ``simplec``: one SIMPLEC iteration with
       ``steady.sweep`` and ``steady.zline_sweeps``. ``identity``: none.
   * - ``steady.newton_armijo_c``
     - ``1e-4``
     - > 0 and < 0.5. Armijo constant :math:`c` of the continuation's line search.
   * - ``steady.newton_backtrack``
     - ``0.5``
     - > 0 and < 1. Backtracking factor :math:`\beta`: :math:`\lambda = 1, \beta, \beta^2, \ldots`
   * - ``steady.newton_lambda_min``
     - ``1e-2``
     - > 0 and < 1. The smallest :math:`\lambda` tried; below it the attempt is rejected.
   * - ``steady.newton_cfl_cut``
     - ``0.25``
     - > 0 and < 1. A rejected or skipped attempt multiplies the CFL by this factor (from infinite
       :math:`\Delta\tau`: ``steady.cfl_max`` times it) and retries.
   * - ``steady.newton_max_retries``
     - ``4``
     - Integer, ≥ 1. Retries of one continuation step; after the last the step returns *failed*.
   * - ``steady.newton_dtau_inf_at_cap``
     - ``1``
     - 0 or 1. 1: once the CFL reaches ``steady.cfl_max``, the next steps run without a pseudo-time term
       (:math:`D = 0`, plain Newton) until an attempt fails. 0: the CFL stays at the cap.
   * - ``steady.newton_ptc_growth``
     - ``1``
     - Finite, ≥ 1. At a finite :math:`\Delta\tau`, when no :math:`\lambda` meets the Armijo condition, accept the
       full step if :math:`|R|` grows by at most this factor. 1: off. Not used at infinite :math:`\Delta\tau`.

The ``steady.newton_*`` keys belong to the Newton–Krylov stage, which the application does not run yet; until then
they are refused as unknown keys. ``steady.newton_preconditioner`` is read and range-checked, then refused with
"has no effect until the run can switch to Newton". Their defaults are provisional. The last seven rows are the pseudo-transient
continuation of :doc:`solution_algorithm`; its CFL ramp reuses ``steady.cfl0``, ``steady.cfl_max``,
``steady.ser_exponent`` and ``steady.ser_max_growth``, with the same range checks.

``simplec``
-----------

See :doc:`solution_algorithm`. Read with ``steady.stage = simplec`` only.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``simplec.pressure_reltol``
     - ``1e-3``
     - ≥ 10 machine epsilon and < 1: relative tolerance of the pressure-correction solve. Looser values can be
       unstable in stratified flow.
   * - ``simplec.alpha_p``
     - ``1``
     - > 0 and ≤ 1: pressure under-relaxation.
   * - ``simplec.rotation_relax``
     - ``0.3``
     - > 0 and ≤ 1: step of the point-implicit rotation, in inertial times. 1 often converges faster (the hill
       cases use it).
   * - ``simplec.w_vnn_cap``
     - ``0`` (off)
     - 0, or finite and ≥ 16: the cap :math:`V` on the vertical VNN number of the :math:`w` rows. Changes the path,
       not the answer. May change across a restart.
   * - ``simplec.w_vnn_cap_slope_margin``
     - ``4``
     - 0 (off), or finite and ≥ 1; refused unless ``simplec.w_vnn_cap`` > 0. The safety factor :math:`s` that lowers the
       cap on sloping ground so the pressure operator stays definite (:doc:`solution_algorithm`). 0 is the cap of
       #32, which can stall where thin cells meet steep slopes; the run then warns over terrain when the mesh
       aspect ratio is above 50.78. Changes the path, not the answer. May change across a restart.

``solver``
----------

See :doc:`solution_algorithm`. Read when ``init.first_guess`` is not ``none`` or with the SIMPLEC stage; refused
otherwise.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``solver.poisson_backend``
     - ``gmres_mlmg``
     - ``gmres_mlmg`` (GMRES with multigrid), ``mlmg`` (multigrid alone) or ``gmres_fft`` (ERF's operator with an
       FFT preconditioner).
   * - ``solver.poisson_reltol``
     - ``1e-10`` (``1e-5`` single)
     - ≥ 10 machine epsilon and < 1. Converged when the residual is below this fraction of the initial one. In
       single precision use about :math:`4\varepsilon N_z^2`.
   * - ``solver.poisson_abstol``
     - ``0``
     - Finite, ≥ 0. Converged when the residual is below this value.
   * - ``solver.poisson_max_iter``
     - ``200``
     - ≥ 1. GMRES iterations, or multigrid cycles with ``mlmg``; the run stops when it is reached.
   * - ``solver.poisson_restart``
     - ``50``
     - ≥ 1. GMRES restart length; not with ``mlmg``.
   * - ``solver.poisson_precond``
     - ``fft``
     - ``fft`` or ``none``; only with ``gmres_fft``.
   * - ``solver.poisson_verbose``
     - ``0``
     - ≥ 0. 1 prints each solve; 2 also prints the residual at each GMRES restart.

``io``
------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``io.plot_file``
     - ``plt``
     - Plotfile prefix; not empty.
   * - ``io.chk_file``
     - ``chk``
     - Checkpoint prefix; not empty.
   * - ``io.write_plotfile``
     - ``1``
     - 0 or 1.
   * - ``io.write_checkpoint``
     - ``1``
     - 0 or 1.
   * - ``io.turb_plot_file``
     - none
     - Prefix of the turbulence plotfile with ERF's variable names (:doc:`output`). Needs the SIMPLEC stage.
   * - ``io.restart``
     - none
     - Checkpoint directory to restart from (:doc:`running`).

``refavg`` (reference-averaging tool)
-------------------------------------

The separate tool ``ablsteady_refavg`` (:doc:`references`) reads only the ``refavg`` namespace.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``refavg.plotfiles``
     - required
     - At least 3 single-level plotfiles with the same variables, domain and box layout, in increasing time.
   * - ``refavg.output``
     - required
     - Archive directory to create; must not exist.
   * - ``refavg.case``
     - required
     - Case name, one word.
   * - ``refavg.erf_sha``
     - required
     - Full 40-character ERF commit.
   * - ``refavg.machine``
     - ``local``
     - ``kestrel-cpu``, ``kestrel-gpu`` or ``local``.
   * - ``refavg.nodes``
     - ``1``
     - ≥ 1.
   * - ``refavg.ranks``
     - ``1``
     - ≥ 1. MPI ranks (CPU) or GPUs of the run.
   * - ``refavg.wall_seconds``
     - ``0``
     - Finite, ≥ 0.
   * - ``refavg.inputs_files``
     - none
     - Files to copy into the archive; each must exist, with distinct names.
   * - ``refavg.verify``
     - none
     - Instead of building, check this archive and print its metadata. The other keys are then refused.
