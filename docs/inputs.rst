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
       ``stretch_ratio`` it must not exceed prob_hi z / n_cell z (:doc:`mesh`).
   * - ``mesh.stretch_ratio``
     - ``0`` (none)
     - Finite, ≥ 1 (an explicit 0 is refused), needs ``mesh.dz0`` > 0. A fixed ratio between successive layers;
       the physical top is then wherever the layers end, and prob_hi z only sets the computational spacing.
   * - ``mesh.zlevels_file``
     - empty
     - A file of n_cell z + 1 increasing heights starting at 0. Not together with ``mesh.dz0``.

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
     - ``flat``, ``agnesi``, ``xyz`` or ``erf_grid`` (:doc:`mesh`).
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

The steady driver (:doc:`driver`) reads ``steady.*`` through ``read_steady_options()``: the five parameters of the
definition of a steady solution (:doc:`steady`) and its own. The application runs the driver with
``steady.stage = simplec`` (:doc:`simplec`); with the default ``none`` it writes the state and stops, and every other
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
     - ≥ 1. Iterations over the whole run (restarts included) after which the driver stops.
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
     - Line sweeps per momentum predictor step (:doc:`predictor`): ≥ 2 and ≤ 100 with ``steady.sweep = zjacobi``
       (one lagged Jacobi sweep is the configuration the sub-phase 1.5 analysis found unstable at large horizontal
       CFL), ≥ 1 and ≤ 100 with ``lusgs`` (each sweep is a forward and a backward pass). With the second-order
       residual the lagged Jacobi sweeps are stable up to a horizontal CFL of 2.4 with 4 sweeps, 3.7 with 8 and 6.1
       with 16 (:doc:`predictor`, sub-phase 1.16); the terrain decks run 16 with ``steady.cfl_max = 4``. A
       horizontally uniform column has exact lagged neighbours and takes any count.
   * - ``steady.sweep``
     - ``zjacobi``
     - ``zjacobi`` or ``lusgs`` (sub-phase 2.1): the line solver of the sweeps, z-line Jacobi (independent columns,
       the x/y neighbours lagged) or line LU-SGS (a symmetric Gauss–Seidel over the columns of each box with the
       exact z-line solve, exact along a wind whose u and v share a sign and stable at every CFL there; with
       u and v of opposite sign one sweep is stable only to a CFL of about 2.5 and two to about 33; a sequential
       thread per box on a GPU until sub-phase 2.2). The sweep's path depends on the box layout
       (``mesh.max_grid_size``); the converged state does not. The LU-SGS hill deck runs two sweeps at
       ``steady.cfl_max = 16`` (two boxes of 64 × 1 columns).
   * - ``steady.machine``
     - ``local``
     - ``local``, ``kestrel-cpu`` or ``kestrel-gpu``: the allocation-unit rate of the accounting (0, 10, 100 AU per
       node-hour).
   * - ``steady.nodes``
     - ``1``
     - ≥ 1. Nodes of the run, for the allocation units.

``solver``
----------

The Poisson layer (:doc:`poisson`) reads ``solver.poisson_*`` through ``read_poisson_options()``. The application
reads them with the initializer (``init.first_guess = profile``, :doc:`initializer`), as do the regression tools
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
       2-norm under GMRES (``gmres_fft``, ``gmres_mlmg``), the maximum norm under ``mlmg``. In single precision
       rounding sets a floor near :math:`\varepsilon N_z^2/2` (:doc:`poisson`).
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
     - ``fft`` (ERF's FFT preconditioner on the flat stretched grid) or ``none``; ``gmres_fft`` only (an error with
       back end B).
   * - ``solver.poisson_backend``
     - ``gmres_fft``
     - ``gmres_fft`` (back end A: ERF's operator, GMRES, FFT preconditioner), ``mlmg`` or ``gmres_mlmg`` (back end
       B: AMReX's ``MLTerrainPoisson`` with multigrid, alone or under GMRES).
   * - ``solver.poisson_verbose``
     - ``0``
     - ≥ 0. 1 prints each solve's iterations, residuals and (singular problems) removed mean, and AMReX GMRES's solve
       time; 2 also GMRES's residual at each restart cycle.

``bc``
------

The boundary conditions (:doc:`boundary`) read ``bc.*`` through ``read_boundary_inputs()`` and
``read_surface_layer_options()``; their table is in :doc:`boundary`. The application reads ``bc.*`` with the
initializer (``init.first_guess = profile``: the face types, the inflow profile and the sounding act on the initial
field; the damping, top-gradient and surface-layer keys are range-checked only, see :doc:`initializer`) and with the
SIMPLEC stage (``steady.stage = simplec``, where every key acts), as does the regression tool ``boundary_fluxes``;
without either, setting them stops the run. The ``bc.inflow_*`` keys are accepted without an inflow face when the
initializer uses them as the first guess (and on a restart of such a run).

``turb``
--------

The k–L closure (:doc:`closure`) reads ``turb.*`` through ``read_closure_options()``; its table is in :doc:`closure`.
The application reads them with the SIMPLEC stage (``steady.stage = simplec``, which always runs the closure), as does
the regression tool ``closure_column``; without the stage, setting them stops the run.

``forcing``
-----------

The forcing of a run (:doc:`residual`), read through ``read_forcing_options()`` (``Source/Forcing.cpp``); with the
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

The convection options of the momentum and :math:`	heta` equations (:doc:`transport`; :math:`k` is always bounded
upwind), read through ``read_convection_options()`` (``Source/Convection.cpp``); with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``convection.limiter``
     - ``van_leer``
     - ``van_leer``, ``van_albada`` or ``upwind``, for momentum and :math:`	heta`.
   * - ``convection.smooth_eps_velocity``
     - ``1e-5``
     - > 0 m/s with its square a finite normal number: the van Albada eps of the momentum equations. With
       ``limiter = van_albada`` only.
   * - ``convection.smooth_eps_theta``
     - ``1e-5``
     - The same for :math:`	heta` (K).

``simplec``
-----------

The SIMPLEC stage's own keys (:doc:`simplec`), read through ``read_simplec_options()``; with ``steady.stage =
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
       is unstable with a stable stratification: :doc:`simplec`).
   * - ``simplec.alpha_p``
     - ``1``
     - > 0 and ≤ 1: the pressure under-relaxation :math:`p' \leftarrow p' + \alpha_p\phi` (SIMPLEC needs none).
   * - ``simplec.rotation_relax``
     - ``0.3``
     - > 0 and ≤ 1: the step of the predictor's point-implicit rotation in inertial times, :math:`\tau_c =
       1/\max(|f|/r, a_P - \sum a_{nb})`, applied to the horizontal mean of the other component's residual
       (:doc:`predictor`). Since sub-phase 1.16 the fluctuation about the mean takes the row's own step, and 1 converges
       as well as 0.3 (Neutral_ABL_Flat: 2507 against 2613 iterations); the default stays 0.3.

``init``
--------

The mass-consistent initializer (:doc:`initializer`) reads ``init.*`` through ``read_init_options()``. Each row matches
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
     - ``none`` (the state starts at rest) or ``profile`` (the inflow profile by height above the ground, corrected
       to a divergence-free field).
   * - ``init.z_agl_min``
     - ``0``
     - Finite, ≥ 0 m. Floor on the height above the ground the profile is evaluated at.
   * - ``init.alpha_h``
     - ``1``
     - Finite, > 0. Horizontal transmissivity :math:`\alpha_h` of the correction.
   * - ``init.alpha_v``
     - ``1``
     - Finite, > 0. Vertical transmissivity :math:`\alpha_v`; only :math:`\alpha_v/\alpha_h` matters for the velocity.
       Ratios below 0.1 need more than the default ``solver.poisson_max_iter`` (:doc:`initializer`).
   * - ``init.k_decay_height``
     - ``0``
     - Finite, ≥ 0 m. :math:`z_i` of the decay of :math:`k` above :math:`z_i`; 0 for no decay.
   * - ``init.obrien``
     - ``0``
     - 0 or 1. The O'Brien column correction of :math:`w_0` before the projection.
