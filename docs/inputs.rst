Inputs
======

Inputs are ``namespace.key = value`` lines. ablsteady owns eleven namespaces: ``mesh``, ``terrain``, ``init``,
``turb``, ``bc``, ``steady``, ``solver``, ``io``, ``forcing``, ``convection`` and ``simplec``.

The run stops, naming the input, when:

- a key in an owned namespace is not read by this version (a misspelling, or a key from a later sub-phase);
- a key's namespace is neither owned nor AMReX's (``amrex``, ``fab``, ``fabarray``, ``vismf``, ``tiny_profiler``,
  ``DistributionMapping``, ``device``, ``hypre``);
- a key is given twice, or a single-value key is given more than one value;
- a value is outside its range.

All errors are reported together. Each row below matches the range check in ``Source/Inputs.cpp``.

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
     - Three finite values: the domain's low corner (m). The z value must be 0: node heights start at the flat
       ground.
   * - ``mesh.prob_hi``
     - required
     - Three finite values, each greater than ``prob_lo``: the domain's high corner (m).
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Three values, each 0 or 1. z must be 0 (ground and top are not periodic).
   * - ``mesh.max_grid_size``
     - ``64``
     - One or two integers ≥ 1: the largest box size in x and y. Boxes are never split in z, so the z-line solves
       stay inside one box.
   * - ``mesh.dz0``
     - ``0``
     - Finite, ≥ 0. First layer thickness of a stretched vertical grid; 0 means uniform layers. Without
       ``stretch_ratio`` it must not exceed prob_hi z / n_cell z (:doc:`governing_curvilinear`).
   * - ``mesh.stretch_ratio``
     - ``0`` (none)
     - Finite, ≥ 1 (an explicit 0 is refused), needs ``mesh.dz0`` > 0. A fixed ratio between successive layers;
       the physical top is then wherever the layers end, and prob_hi z only sets the computational spacing.
   * - ``mesh.zlevels_file``
     - empty
     - A file of n_cell z + 1 increasing heights starting at 0. Not together with ``mesh.dz0``.
   * - ``mesh.hierarchy_levels``
     - ``1``
     - An integer from 1 to 10: the levels L of the x/y grid hierarchy (:doc:`governing_curvilinear`). 1 is the mesh
       alone. With L > 1, n_cell x and y must be divisible by :math:`2^{L-1}` and ``mesh.max_grid_size`` x and y
       must be multiples of :math:`2^{L-1}`. Level 0 is then laid out so that every level coarsens it exactly, which
       can differ from the L = 1 layout (the last bits of a solve can differ). The levels are built, reported and
       written with the plotfile; the solver runs on level 0 alone unless ``steady.sequencing = 1`` (grid sequencing
       over all levels). The test-only tools refuse L > 1.
       May change across a restart.

``terrain``
-----------

A key that has no effect for the chosen ``terrain.type`` or ``terrain.smoothing`` stops the run.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``terrain.type``
     - ``flat``
     - ``flat``, ``agnesi``, ``xyz`` or ``erf_grid`` (:doc:`governing_curvilinear`).
   * - ``terrain.file``
     - required for ``xyz``, ``erf_grid``
     - The point file (``x y z`` per line) or ERF's structured terrain file.
   * - ``terrain.smoothing``
     - ``1``
     - 0 (BTF), 1 (STF, Klemp 2011) or 2 (Sullivan).
   * - ``terrain.gamma_m``
     - ``0.5``
     - > 0 and < 1; STF only. Smallest layer thickness as a fraction of the flat layer.
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Finite, > 0 (m); ``agnesi`` only.
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Finite, > 0 (m); ``agnesi`` only. The hill is at half height this far from its crest.
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - ``ridge_x``, ``ridge_y`` or ``round``; ``agnesi`` only.
   * - ``terrain.idw_neighbors``
     - ``6``
     - ≥ 1; ``xyz`` only. Points used for each height.
   * - ``terrain.idw_power``
     - ``2``
     - Finite, > 0; ``xyz`` only. Weight exponent, :math:`d^{-p}`.
   * - ``terrain.hierarchy_passes``
     - ``1``
     - An integer from 0 to 64; only with ``mesh.hierarchy_levels`` > 1. The ground of hierarchy level l is smoothed
       by this many times :math:`2^l - 1` passes (:doc:`governing_curvilinear`).

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
   * - ``io.restart``
     - empty
     - Checkpoint directory to restart from (see :doc:`running`).
   * - ``io.write_plotfile``
     - ``1``
     - 0 or 1.
   * - ``io.write_checkpoint``
     - ``1``
     - 0 or 1.
   * - ``io.turb_plot_file``
     - empty
     - Prefix of the turbulence diagnostics plotfile with ERF's variable names (:doc:`output`); empty: none. Needs
       ``steady.stage = simplec``.

``steady``
----------

The steady driver (:doc:`solver`) reads ``steady.*`` through ``read_steady_options()``: the five parameters of the
definition of a steady solution (:doc:`solver`) and its own. The application runs the driver with
``steady.stage = simplec`` (:doc:`solver`); with the default ``none`` it writes the state and stops, and every other
``steady.*`` key is refused (it would have no effect). Each row matches the check in ``Source/SteadyCriterion.cpp`` or
``Source/Driver.cpp``; a repeated key is refused.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Input
     - Default
     - Range and meaning
   * - ``steady.stage``
     - ``none``
     - ``none`` (write the state and stop) or ``simplec`` (run the SIMPLEC stage under the driver).
   * - ``steady.residual_tol``
     - ``5e-5``
     - Finite, > 0. Every scaled residual must be below it.
   * - ``steady.residual_drop_orders``
     - ``4``
     - ≥ 0 and ≤ the largest decimal exponent of the precision (308 in double, 38 in single). Every residual must
       also have dropped this many orders below its reference.
   * - ``steady.residual_ref_iterations``
     - ``5``
     - ≥ 1. A residual's reference is its largest value over this many first iterations.
   * - ``steady.monitor_rel_tol``
     - ``1e-3``
     - Finite, > 0. Every monitor must change by at most this fraction of max(\|value\|, scale) over the window, and be
       at most this far from its extrapolated limit.
   * - ``steady.monitor_window``
     - ``50``
     - ≥ 2. The window, in iterations.
   * - ``steady.max_iterations``
     - ``10000``
     - ≥ 1. Iterations over the whole run (restarts included) after which the driver stops. With
       ``steady.sequencing = 1`` it caps each level of the sequence separately.
   * - ``steady.cfl0``
     - ``1``
     - Finite, > 0. The CFL of the first iterations (while the residual reference is taken).
   * - ``steady.cfl_max``
     - ``10``
     - Finite, ≥ ``steady.cfl0``. The cap of the SER ramp.
   * - ``steady.ser_exponent``
     - ``1``
     - Finite, ≥ 0. :math:`p` in :math:`\mathrm{CFL} = \mathrm{CFL}_0\, r^{-p}`; 0 holds the CFL at ``cfl0``.
   * - ``steady.ser_max_growth``
     - ``2``
     - Finite, ≥ 1. The CFL may grow by at most this factor per iteration.
   * - ``steady.switch_window``
     - ``20``
     - ≥ 1. Iterations below ``switch_residual`` before a forward switch, the window of the back rule, and the
       hysteresis after any switch.
   * - ``steady.switch_residual``
     - ``1e-2``
     - > 0 and ≤ 1. The relative residual at or below which the next stage is entered.
   * - ``steady.switch_back_growth``
     - ``10``
     - Finite, > 1. Growth of the relative residual over a window that sends the run back a stage.
   * - ``steady.zline_sweeps``
     - ``2``
     - Line sweeps per momentum predictor step (:doc:`solver`): ≥ 2 and ≤ 100 with ``steady.sweep = zjacobi``
       (one lagged Jacobi sweep is unstable at large horizontal CFL), and ≥ 1 and ≤ 100 with ``lusgs`` (each sweep
       is a forward and a backward pass). With the second-order residual the lagged Jacobi sweeps are stable up to a
       horizontal CFL of 2.4 with 4 sweeps, 3.7 with 8 and 6.1 with 16 (:doc:`solver`, sub-phase 1.16); the terrain
       decks run 16 with ``steady.cfl_max = 4``. A
       horizontally uniform column has exact lagged neighbours and takes any count.
   * - ``steady.sweep``
     - ``zjacobi``
     - ``zjacobi`` or ``lusgs`` (sub-phase 2.1): the line solver of the sweeps, z-line Jacobi (independent columns,
       the x/y neighbours lagged) or line LU-SGS (a symmetric Gauss–Seidel over the columns of each box with the
       exact z-line solve, exact along a wind whose u and v share a sign and stable at every CFL there; with
       u and v of opposite sign the defect correction is stable with one sweep only to a CFL of about 2.5 and with
       two to about 33 (:doc:`solver`); a sequential
       thread per box on a GPU unless ``steady.line_order = diagonal``). The sweep's path depends on the box layout
       (``mesh.max_grid_size``); the converged state does not. The LU-SGS hill deck runs two sweeps at
       ``steady.cfl_max = 16`` (two boxes of 64 × 1 columns).
   * - ``steady.line_ghosts``
     - ``lagged``
     - ``lagged`` or ``zero`` (sub-phase 2.2); refused unless ``steady.sweep = lusgs``. ``lagged``: the x/y
       neighbours across a box face take the previous pass's values (sub-phase 2.1, bit-identical). ``zero``:
       box-local sweeps with δ = 0 outside the swept region, no communication inside the sweeps, each C-grid
       face swept by one box (its low face; the domain's high face, and on a periodic axis the seam, by the last
       box).
   * - ``steady.line_overlap``
     - ``0``
     - An integer ≥ 0 and ≤ 8, and ≤ the narrowest box width along each axis it applies to (not along a
       one-cell axis); > 0 only with ``steady.line_ghosts = zero``; refused unless ``steady.sweep = lusgs``. Each box sweeps its owned region grown by
       this many columns and keeps the owned part. With ``steady.sequencing = 1`` the bound holds on every
       hierarchy level, whose boxes are :math:`2^l` times narrower.
   * - ``steady.line_order``
     - ``sequential``
     - ``sequential`` or ``diagonal``; refused unless ``steady.sweep = lusgs``. ``diagonal`` sweeps the
       wavefronts i + j = const, one column per thread: the same result as ``sequential``, bit for bit.
   * - ``steady.line_alternate``
     - ``0``
     - ``0`` or ``1`` (beyond the plan of sub-phase 2.2); refused unless ``steady.sweep = lusgs``. ``1``: the odd
       sweeps (the second, fourth, ...) of each call take the j-reversed pair, (+i, −j) then (−i, +j), so a wind
       whose u and v have opposite signs is swept along its own direction once per two sweeps (:doc:`solver`).
       ``1`` is refused with ``steady.zline_sweeps = 1`` (there is no odd sweep to reverse); on a domain one cell
       wide in y only v (two y faces per column) changes.
   * - ``steady.sequencing``
     - ``0``
     - ``0`` or ``1`` (sub-phase 2.5). ``1``: grid sequencing, coarse to fine, over the ``mesh.hierarchy_levels``
       levels (:doc:`solver`): the initializer and the SIMPLEC stage to a loose stop on the coarsest level, then on
       each finer level the prolonged state to a loose stop, and on level 0 the strict stop. Needs
       ``mesh.hierarchy_levels`` ≥ 2 and at least 2 cells in x and in y on the coarsest level
       (``mesh.n_cell`` / 2^(L−1)). ``init.*`` then acts on the coarsest level only, and ``steady.max_iterations``
       caps each level separately. The coarse levels run on a fresh start only; a restart continues level 0, and the
       ``steady.sequencing*`` keys then have no effect (the log says so).
   * - ``steady.sequencing_drop_orders``
     - ``2``
     - Finite, > 0 and ≤ ``steady.residual_drop_orders``; refused unless ``steady.sequencing = 1``. A coarse level
       stops once every residual is this many orders below its reference (the coarsest level's, :doc:`solver`) or
       below ``steady.residual_tol``.
   * - ``steady.sequencing_monitor_tol``
     - ``0``
     - 0 (off) or finite, > 0 and < 1; refused unless ``steady.sequencing = 1``. When > 0 a coarse level also
       stops once every monitor changed by at most this fraction of max(\|value\|, scale) over
       ``steady.monitor_window`` iterations. Off by default: after a prolongation the monitors settle before the
       residuals drop (:doc:`solver`).
   * - ``steady.machine``
     - ``local``
     - ``local``, ``kestrel-cpu`` or ``kestrel-gpu``: the allocation-unit rate of the accounting (0, 10, 100 AU per
       node-hour).
   * - ``steady.nodes``
     - ``1``
     - ≥ 1. Nodes of the run, for the allocation units.
   * - ``steady.newton_gmres_restart``
     - ``30``
     - Integer, ≥ 5 and ≤ 200. GMRES restart length of the Newton step (:doc:`solver`).
   * - ``steady.newton_fd_eps_scale``
     - :math:`\sqrt{\epsilon_\mathrm{mach}}` of ``amrex::Real``
     - ≥ :math:`\epsilon_\mathrm{mach}^{3/4}` of ``amrex::Real`` (1.8189894035458565e-12 in double,
       6.41553051e-06 in single) and ≤ 1e-2. Relative step :math:`h` of the Jacobian-vector product
       (:doc:`solver`). Away from the default :math:`Jv` loses digits: relative error 7.2e-8 at the default, 9.4e-3
       at 1e-2 and about 1e-3 at the lower bound (double). Below the bound the error grows as :math:`1/h`; at
       :math:`\epsilon_\mathrm{mach}` it is about 1 and Newton stalls or diverges.
   * - ``steady.newton_krylov_rtol``
     - ``1e-2``
     - > 0 and < 1. Relative tolerance of GMRES in each Newton step. Below the forward difference's floor
       (1.4e-8 on the Ekman column, up to 1.7e-6 on the hill case, in double) GMRES reports convergence without
       the true residual following (:doc:`solver`).
   * - ``steady.newton_mem_budget_gb``
     - ``70``
     - > 0 and ≤ 80. Krylov memory per rank in GB (10\ :sup:`9` bytes); the run aborts when the estimate of
       ``2 steady.newton_gmres_restart + 3`` state-size vectors on the rank with the most points exceeds it. The
       step's own buffers (about 9 vectors more) are not counted: leave that margin.

The four ``steady.newton_*`` rows are read by ``read_newton_options()`` (sub-phase 3.1). Their defaults are
provisional until sub-phase 3.6 (PR-28) tunes them. The application runs no Newton step before sub-phase 3.5, so
until then it refuses these keys as unknown inputs.

``solver``
----------

The Poisson layer (:doc:`solver`) reads ``solver.poisson_*`` through ``read_poisson_options()``. The application
reads them with the initializer (``init.first_guess = profile``, :doc:`initial_conditions`), as do the regression tools
``poisson_projection`` and ``wall_distance``; without the initializer the application has no Poisson solve yet, so
setting these keys stops the run, like any other key it does not read. Each row matches the range check in
``Source/Poisson.cpp``.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``solver.poisson_reltol``
     - ``1e-10`` (``1e-5`` single)
     - ≥ 10 machine epsilon and < 1. Converged when the residual is below this fraction of the initial one: the
       2-norm under GMRES (``gmres_fft``, ``gmres_mlmg``), the maximum norm under ``mlmg``. Back end B measures
       the :math:`J`-weighted residual (the cell-volume norm), back end A the plain one, so the same value is not
       the same divergence where :math:`J` varies (:doc:`solver`). In single precision rounding sets a floor near
       :math:`\varepsilon N_z^2/2`; pure ``mlmg`` (maximum norm) may not reach the single-precision default on
       large meshes, ``gmres_mlmg`` does.
   * - ``solver.poisson_abstol``
     - ``0``
     - Finite, ≥ 0. Or when the residual (in the same norm) is below this value.
   * - ``solver.poisson_max_iter``
     - ``200``
     - ≥ 1. GMRES iterations, or MLMG V-cycles under ``mlmg``; reaching it stops the run.
   * - ``solver.poisson_restart``
     - ``50``
     - ≥ 1. GMRES restart length (ERF's value); ``gmres_fft`` and ``gmres_mlmg`` only (an error with ``mlmg``).
   * - ``solver.poisson_precond``
     - ``fft``
     - ``fft`` (ERF's FFT preconditioner on the flat stretched grid) or ``none``; given only with
       ``solver.poisson_backend = gmres_fft`` (an error with back end B, which includes the default). The
       initializer's ``gmres_fft`` fallback below uses the default ``fft``; to change it, set
       ``solver.poisson_backend = gmres_fft``, which moves the SIMPLEC stage to back end A too.
   * - ``solver.poisson_backend``
     - ``gmres_mlmg``
     - ``gmres_mlmg`` (back end B under GMRES), ``mlmg`` (back end B: AMReX's ``MLTerrainPoisson`` with multigrid)
       or ``gmres_fft`` (back end A: ERF's operator, GMRES, FFT preconditioner). Back end B has no split z-face
       coefficient, so with ``init.alpha_v`` ≠ ``init.alpha_h`` the initializer's projection uses ``gmres_fft``
       whatever this says, and logs it; the SIMPLEC stage keeps this value (:doc:`solver`).
   * - ``solver.poisson_verbose``
     - ``0``
     - ≥ 0. 1 prints each solve's iterations, residuals and (singular problems) removed mean, and AMReX GMRES's solve
       time; 2 also GMRES's residual at each restart cycle.

``bc``
------

The boundary conditions (:doc:`boundary_conditions`) read ``bc.*`` through ``read_boundary_inputs()`` and
``read_surface_layer_options()``. The application reads ``bc.*`` with the
initializer (``init.first_guess = profile``: the face types, the inflow profile and the sounding act on the initial
field; the damping, top-gradient and surface-layer keys are range-checked only, see :doc:`initial_conditions`) and with the
SIMPLEC stage (``steady.stage = simplec``, where every key acts), as does the regression tool ``boundary_fluxes``;
without either, setting them stops the run. The ``bc.inflow_*`` keys are accepted without an inflow face when the
initializer uses them as the first guess (and on a restart of such a run).
A key with no effect for the chosen face types stops the run.

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
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Three finite values (m/s), the velocity of a ``no_slip_wall`` face (``xlo``, ``xhi``, ``ylo``, ``yhi``,
       ``zlo``, ``zhi``); the normal component must be 0. Refused on another face type.
   * - ``bc.inflow_profile``
     - ``log_law``
     - ``log_law``, ``power_law`` or ``file``. The inflow keys need an inflow face, a damping layer or the
       initializer's first guess (``init.first_guess = profile``, which uses them on periodic domains too).
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
   * - ``bc.rayleigh_fields``
     - ``u v w theta tke``
     - One or more of ``u v w theta tke`` (any other name is refused): the fields the Rayleigh layer damps; needs
       ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.rayleigh_height``
     - ``absolute``
     - ``absolute`` or ``above_ground``: the height the layer's coefficient and θ's reference use (u, v and k
       always use the height above the ground); needs ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.sponge_width``
     - 0 (off)
     - Finite, ≥ 0 (m); > 0 needs an inflow or outflow face and ``bc.sponge_strength`` > 0.
   * - ``bc.sponge_strength``
     - 0
     - Finite, ≥ 0 (1/s); needs ``bc.sponge_width`` > 0.
   * - ``bc.surface_z0``, ``bc.surface_zref``, ``bc.surface_stability``, ``bc.surface_heat_flux``, ``bc.surface_wall_k``,
       ``bc.surface_denominator_floor``
     - below
     - The surface layer's options (``read_surface_layer_options``); refused unless ``bc.zlo = surface_layer``.

The surface-layer keys (refused unless ``bc.zlo = surface_layer``), read by ``read_surface_layer_options()``:

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - Input
     - Default
     - Range and effect
   * - ``bc.surface_z0``
     - 0.1 m
     - > 0, roughness length
   * - ``bc.surface_zref``
     - 0
     - ≤ 0: each column's first cell centre; otherwise a vertical height between the first and last cell-centre
       heights of every column
   * - ``bc.surface_stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``
   * - ``bc.surface_heat_flux``
     - 0 K m/s
     - finite; ``heat_flux`` only (refused with ``neutral``, where it would be ignored)
   * - ``bc.surface_wall_k``
     - ``dirichlet``
     - ``dirichlet`` or ``source``
   * - ``bc.surface_denominator_floor``
     - ``log_law``
     - ``log_law`` or ``one``: the floor :math:`m` of the MOST denominator (:doc:`boundary_conditions`).
       ``log_law``: :math:`m = \min(\ln((d + z_0)/z_0), 1)`, the neutral log law at every height. ``one``:
       :math:`m = 1`, sub-phase 1.7's form, which lowers the neutral :math:`u_*` below :math:`d = 1.72\,z_0`.
       The two are identical where :math:`\ln((d + z_0)/z_0) \ge 1`. ``init.first_guess = column`` uses the same floor.

The wall model uses the closure's ``turb.Cmu0`` and ``turb.kappa``.

``turb``
--------

The k–L closure (:doc:`governing_cartesian`) reads ``turb.*`` through ``read_closure_options()``. The application
reads them with the SIMPLEC stage (``steady.stage = simplec``, which always runs the closure), as does the
regression tool ``closure_column``; without the stage, setting them stops the run. Each row matches the check in
``check_closure_options()``.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and effect
   * - ``turb.Cmu0``
     - 0.5562
     - Finite, > 0: :math:`C_{\mu 0}` (AL01, ERF).
   * - ``turb.Cb``
     - 0.35
     - Finite, > 0: :math:`C_b` of the stable and unstable lengths.
   * - ``turb.Rt_crit``
     - −1
     - Finite, ≤ 0: the smoothing starts below it.
   * - ``turb.Rt_min``
     - −3
     - Finite, < ``turb.Rt_crit`` and > −3.6: the smoothed :math:`R_t` stays above it.
   * - ``turb.max_geom_lscale``
     - 30 m
     - Finite, > 0: the harmonic cap of :math:`l_g` (about :math:`\kappa \times 0.1 z_i`).
   * - ``turb.tke_floor``
     - 0
     - 0 (machine epsilon) or ≥ machine epsilon (m²/s²): the floor of :math:`k`; a value between is refused rather
       than raised.
   * - ``turb.sigma_k``
     - 1
     - Finite, > 0: :math:`\gamma_k = \nu_{mol} + \nu_t/\sigma_k`.
   * - ``turb.nu_molecular``, ``turb.alpha_molecular``
     - 0, 0
     - Finite, ≥ 0 (m²/s): molecular viscosity and heat diffusivity added to the turbulent ones.
   * - ``turb.kappa``
     - 0.41
     - Finite, > 0: von Kármán's constant in :math:`l_g`.
   * - ``turb.nut_relax``
     - 1
     - Finite, in (0, 1]: the relaxation weight of :math:`\nu_t` (1: none).
   * - ``turb.frozen_switch``
     - 0
     - 0 or 1: freeze the branch decisions at the state before a Jacobian (``Closure::freeze_switches``).
   * - ``turb.face_mean``
     - ``arithmetic``
     - ``arithmetic`` or ``logarithmic`` (above).
   * - ``turb.wall_distance``
     - ``terrain_height``
     - ``terrain_height`` or ``tucker`` (above).

``forcing``
-----------

The forcing of a run (:doc:`governing_cartesian`), read through ``read_forcing_options()`` (``Source/Forcing.cpp``); with the
SIMPLEC stage only. Each row matches ``check_forcing_options``; a repeated key and a key with no effect are refused.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``forcing.gravity``
     - ``9.81``
     - Finite, ≥ 0 m/s². 0 turns the buoyancy (and the point-implicit buoyancy of the predictor) off.
   * - ``forcing.theta0``
     - ``300``
     - Finite, > 0 K. The Boussinesq reference potential temperature.
   * - ``forcing.coriolis``
     - ``0``
     - 0 or 1. The Coriolis force (and the point-implicit rotation of the predictor).
   * - ``forcing.latitude``
     - ``0``
     - In [−90, 90] degrees. With ``coriolis = 1`` only.
   * - ``forcing.rotation_period``
     - ``86400``
     - Finite, > 0 s (:math:`f = 4\pi\sin\varphi/T`). With ``coriolis = 1`` only.
   * - ``forcing.coriolis_3d``
     - ``1``
     - 0 or 1: the :math:`\cos\varphi` terms on u and w (ERF's ``coriolis_3d``). With ``coriolis = 1`` only.
   * - ``forcing.geostrophic``
     - ``0``
     - 0 or 1: the geostrophic forcing (needs ``coriolis = 1`` and a non-zero latitude).
   * - ``forcing.ug``, ``forcing.vg``
     - ``0``
     - Finite m/s: the geostrophic wind. With ``geostrophic = 1`` only.
   * - ``forcing.pressure_gradient``
     - ``0 0 0``
     - Three finite values (m/s²): the constant kinematic pressure-gradient force on u, v, w (ERF's
       ``abl_pressure_grad``).

``convection``
--------------

The convection options of the momentum and :math:`\theta` equations (:doc:`discretization`; :math:`k` is always bounded
upwind), read through ``read_convection_options()`` (``Source/Convection.cpp``); with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``convection.limiter``
     - ``van_leer``
     - ``van_leer``, ``van_albada`` or ``upwind``, for momentum and :math:`\theta`.
   * - ``convection.smooth_eps_velocity``
     - ``1e-5``
     - > 0 m/s with its square a finite normal number: the van Albada eps of the momentum equations. With
       ``limiter = van_albada`` only.
   * - ``convection.smooth_eps_theta``
     - ``1e-5``
     - The same for :math:`\theta` (K).

``simplec``
-----------

The SIMPLEC stage's own keys (:doc:`solver`), read through ``read_simplec_options()``; with ``steady.stage =
simplec`` only. The sweep count is ``steady.zline_sweeps``.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``simplec.pressure_reltol``
     - ``1e-3``
     - ≥ 10 machine epsilon and < 1: the relative tolerance of the pressure correction's Poisson solve (the plan's 0.1
       is unstable with a stable stratification: :doc:`solver`).
   * - ``simplec.alpha_p``
     - ``1``
     - > 0 and ≤ 1: the pressure under-relaxation :math:`p' \leftarrow p' + \alpha_p\phi` (SIMPLEC needs none).
   * - ``simplec.rotation_relax``
     - ``0.3``
     - > 0 and ≤ 1: the step of the predictor's point-implicit rotation in inertial times, :math:`\tau_c =
       1/\max(|f|/r, a_P - \sum a_{nb})`, applied to the horizontal mean of the other component's residual
       (:doc:`solver`). Since sub-phase 1.16 the fluctuation about the mean takes the row's own step, and 1 converges
       as well as 0.3 (Neutral_ABL_Flat: 2507 against 2613 iterations); the default stays 0.3.
   * - ``simplec.w_vnn_cap``
     - ``0`` (off)
     - 0, or finite and ≥ 16 (dimensionless): the cap :math:`V` on the vertical VNN number of the :math:`w` rows. On those rows alone
       :math:`1/\Delta\tau` becomes :math:`1/\Delta\tau + (a_{ZM} + a_{ZP})/V` in :math:`a_P` and in the SIMPLEC
       denominator :math:`d` (not on the Dirichlet ground and top rows). It changes the pseudo-time path, not the
       steady answer (:doc:`solver`). 0 is bit-identical to the solver without the input. May change across a
       restart (nothing about it is in the checkpoint). Over terrain the cap can make the pressure correction
       unstable where thin cells meet slopes; the run warns when the mesh's max aspect ratio is above the tested
       50.78 (Askervein, 1 m first cell). The aspect ratio is a proxy for the condition :math:`1 + \mathrm{VNN}/V <
       1/h^2`, with :math:`h` the local slope (:doc:`solver`). Above it a larger cap is safer (64 converged on Askervein
       with a 0.4 m first cell, where 16 burst).

``init``
--------

The mass-consistent initializer (:doc:`initial_conditions`) reads ``init.*`` through ``read_init_options()``. Each row matches
the range check in ``Source/Initializer.cpp``; the keys other than ``first_guess`` are refused without
``first_guess = profile`` (they would have no effect). A restart with the initializer on is refused.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``init.first_guess``
     - ``none``
     - ``none`` (the state starts at rest), ``profile`` (the inflow profile by height above the ground) or
       ``column`` (sub-phase 2.1b: the deck's horizontally uniform column solved as a 1-D problem first, in a
       fraction of a second, its u, v and k broadcast by height above the ground; the requirements it refuses
       otherwise, naming this key, are in :doc:`initial_conditions`: the stage, a horizontal forcing, the neutral surface
       layer with the Dirichlet wall k at the first cell, the top type, the Rayleigh fields, no ``k_decay_height``),
       each corrected to a divergence-free field unless ``init.projection = 0`` (:doc:`initial_conditions`).
   * - ``init.projection``
     - ``1``
     - 0 or 1: the mass-consistent correction of the first guess; 0 keeps the first guess as it is (the stage's own
       pressure correction then removes its divergence). ``alpha_h`` and ``alpha_v`` have no effect with 0 (refused).
   * - ``init.column_tol``
     - ``1e-8``
     - Finite, > 0: the 1-D column's convergence tolerance, the largest change of u, v (scaled by the geostrophic
       or first-cell speed) and k (by its maximum) over one backward Euler step. Lifted to 100 epsilon of ``Real``
       (1.2e-5 in single precision, logged) since the closure's point kernels run in ``Real``. With
       ``first_guess = column`` only (refused otherwise).
   * - ``init.column_max_iterations``
     - ``5000``
     - ≥ 1: the column's step limit; the run aborts, naming the key and the last change, if the column has not
       converged by then. With ``first_guess = column`` only.
   * - ``init.z_agl_min``
     - ``0``
     - Finite, ≥ 0 m. Floor on the height above the ground the profile is evaluated at.
   * - ``init.alpha_h``
     - ``1``
     - Finite, > 0. Horizontal transmissivity :math:`\alpha_h` of the correction.
   * - ``init.alpha_v``
     - ``1``
     - Finite, > 0. Vertical transmissivity :math:`\alpha_v`; only :math:`\alpha_v/\alpha_h` matters for the velocity.
       Ratios below 0.1 need more than the default ``solver.poisson_max_iter`` (:doc:`initial_conditions`). With
       ``alpha_v`` ≠ ``alpha_h`` the projection runs on ``gmres_fft`` (``solver.poisson_backend``).
   * - ``init.k_decay_height``
     - ``0``
     - Finite, ≥ 0 m. :math:`z_i` of the decay of :math:`k` above :math:`z_i`; 0 for no decay.
   * - ``init.obrien``
     - ``0``
     - 0 or 1. The O'Brien column correction of :math:`w_0` before the projection.
