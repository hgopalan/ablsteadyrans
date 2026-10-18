Testing
========

Conventions
-----------

These follow the development plan.

- **Unit tests (UT-p.s.n)** use GoogleTest. They are deterministic, finish in seconds on 1–4 ranks, and run on
  every commit: ``ctest -L unit``.
- **Regression tests (RT-p.s.n)** run the application and compare with stored gold plotfiles using
  ``amrex_fcompare``. Their check scripts print *measured / target / tolerance / pass* rows.
  ``ctest -L smoke`` runs the structural checks CI does on every pull request; ``ctest -L physics`` runs the
  converged checks nightly on Kestrel.
- **Performance tests (PT-p.s.n)** run nightly or weekly on Kestrel. They track iterations, wall time and
  allocation-unit cost; a slowdown of more than 10 % fails.
- **Gold files** change only in a reviewed commit.

Every test is first shown to fail on the bug it guards against, or on a one-line mutant of the code it covers.
Tests the plan does not list are marked *beyond the plan*.

The tests run with the GPU CI flags (see :doc:`running`). Each test runs in its own directory under
``<build>/Tests/run/``.

Current tests
-------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - CTest
     - What it checks
   * - ``unit_np1``, ``unit_np2``, ``unit_np4``
     - The GoogleTest suite on 1, 2 and 4 ranks (below).
   * - ``unit_notiling``
     - The same suite with tiling off, as on a GPU.
   * - ``smoke_clean``, ``smoke_run``, ``smoke_restart``, ``smoke_restart_matches``
     - A 2-rank run; a 3-rank restart on another box layout, whose log must show the restart; the two plotfiles
       identical under ``amrex_fcompare`` at zero tolerance.
   * - ``abort_unknown_key``, ``abort_range``, ``abort_repeated_key``, ``abort_restart_mesh``,
       ``abort_restart_overwrite``
     - Decks that must stop, with the input named in the message.
   * - ``rt_askervein_mesh``
     - RT-1.4.1 and RT-1.4.2: ERF's Askervein mesh coarsened to 60 × 60 × 18 (20 m × 1.08 layers, BTF, ERF's terrain
       file): statistics of ``z`` and ``detJ`` and the quality report against ``Tests/Gold/askervein_mesh.json``.
       Double-precision builds only: the statistics use ERF's plotfile reader, which reads double-precision
       plotfiles.
   * - ``rt_steep_ridge_quality``
     - RT-1.4.2: a 600 m ridge under STF: the report against its gold, including both warnings (STF capped, top
       not flat). Double-precision builds only.
   * - ``rt_hill2d_run``, ``rt_hill2d_matches``
     - RT-1.4.3: a run without a solve writes the mapped-grid plotfile of the 2-D hill; ``amrex_fcompare`` against
       ``Tests/Gold/hill2d_plt00000`` at relative tolerance 1e-12 (1e-5 in single precision). The plotfile holds
       cell fields only, so the nodal heights are not compared here; UT-1.4.7 compares them with ERF.

       The golds of the three mesh tests above were written by ablsteady itself (``Tests/check_mesh.py --write``):
       they detect changes, they do not show correctness. Correctness of the same code paths comes from UT-1.4.7
       (ERF's own output) and the closed-form unit tests.
   * - ``abort_poisson_options``
     - ``poisson_projection`` with ``solver.poisson_reltol`` repeated and out of range: both named.
       (``abort_unknown_key`` checks that the application refuses ``solver.poisson_reltol``.)
   * - ``rt_poisson_askervein``, ``rt_poisson_askervein_matches``
     - RT-1.10.1: ``poisson_projection`` projects a fixed wind on the Askervein mesh (30 × 30 × 18). Iterations at
       most 10 % above ``Tests/Gold/poisson_askervein.json``, divergence after below 10 × the printed tolerance of
       before, and
       :math:`\phi` against ``Tests/Gold/poisson_askervein_plt00000`` at relative 1e-8. Double precision only (the
       golds are double-precision runs); the gold was written by ablsteady, so it detects changes.
   * - ``rt_poisson_askervein_mlmg``, ``rt_poisson_askervein_mlmg_matches``
     - RT-1.11.1: the same projection with back end B (MLMG, 6 iterations): iterations at most 10 % above
       ``Tests/Gold/poisson_askervein_mlmg.json`` (written by ablsteady's first run) and the divergence check, and
       :math:`\phi` against back end A's gold
       ``Tests/Gold/poisson_askervein_plt00000`` at relative 1e-8, since both solve the same operator.
   * - ``rt_walldist_hill2d``, ``rt_walldist_hill2d_matches``
     - RT-1.11.2: ``wall_distance`` on ERF's Neutral_Hill_2D mesh: ERF's wall-distance checks (``Tests/check_walldist.py``:
       max relative ≤ 0.15, mean ≤ 0.03, ≤ 0.2 cells below 100 m, all finite), and the field against
       ``Tests/Gold/walldist_hill2d_plt00000`` at relative 1e-8 (a gold written by ablsteady). Double precision only.
   * - ``rt_init_askervein``, ``rt_init_askervein_matches``
     - RT-1.12.1: the mass-consistent initial field on the Askervein mesh (30 × 30 × 18; ``Tests/Decks/inputs_init_askervein``:
       log law u* = 0.6 m/s, z₀ = 0.03 m along x, inflow at xlo, outflow at xhi, :math:`\alpha_h = 1`,
       :math:`\alpha_v = 0.5`; divergence 4.05e-1 → 6.1e-12 1/s). ``Tests/check_init.py askervein``: iterations at most 10 % above the gold's 13,
       divergence after at most 10 × ``solver.poisson_reltol`` times before, divergence before (the printed 10 digits),
       max :math:`|w|` and the largest horizontal speed at the gold's values (1e-6), :math:`\theta` = 300 K everywhere, :math:`p' = 0`,
       :math:`k > 0` (``Tests/Gold/init_askervein.json``); the state against ``Tests/Gold/init_askervein_plt00000``
       at relative 1e-8. Two ranks; double precision only. :math:`\lambda` is not in the plotfile (a variation from
       the plan): the state plotfile is the regression artefact, and the unit tests pin :math:`\lambda`.
   * - ``rt_init_bump_trend``
     - RT-1.12.2: the hilltop speed-up over FastWindTerrain's bump (``Tests/Decks/inputs_init_bump``,
       ``Tests/Decks/terrain_bump_fwt.txt``: its ``regtests/phase6_solve_correction`` hill and power-law inflow) for
       :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (five runs on two ranks), against FastWindTerrain's own at commit
       648f8f1 (``Tests/Gold/init_bump_trend.json``, written from its plotfiles by ``Tests/fwt_bump_speedup.py``): the
       same ranking of the five values at 40, 80 and 150 m above the ground (at 20 m two values tie within 0.002 in
       both codes), every speed-up above 1 at 20, 40, 80 and 150 m, the two codes within 0.02, 0.01, 0.005 and 0.003
       of each other at those heights (measured 0.011, 0.0057, 0.0020, 0.0008), and ablsteady's own values and
       iterations pinned (1e-6, 10 %). Double precision only (:doc:`initializer`).
   * - ``rt_bc_neutral_flat``, ``rt_bc_neutral_flat_matches``
     - RT-1.8.1: ``boundary_fluxes`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, 1000 m) with a fixed log-law state
       (u* = 0.35 m/s, z₀ = 0.1 m, 15°), the surface layer on the ground and a slip top: u* recovered at every ground
       cell to 1e-10, the ground stress :math:`u_*^2(\cos, \sin)` to 1e-10, :math:`k_w = u_*^2/C_{\mu 0}^2`; the
       residual of every cell equals its closed form (``Tests/check_bc.py``: arithmetic face ν times the log law's
       differences, the wall flux below, no stress through the slip top) to 1e-6 of :math:`u_*^2/\Delta z` for u and
       v (1e-13 measured) and 1e-9 K/s for θ (the sounding's diffusion with the top gradient), with the first cell's
       truncation −0.0971 and the top cell's lost flux +1.000 of :math:`u_*^2\cos/\Delta z` pinned; R_k and R_w zero;
       every field finite. The plotfile against ``Tests/Gold/bc_neutral_flat_plt00000`` at relative 1e-10 (a gold
       written by ablsteady). Double precision only.
   * - ``rt_closure_column``, ``rt_closure_column_matches``
     - RT-1.9.1: ``closure_column`` on ERF's Neutral_ABL_Flat mesh (8 × 8 × 64, neutral) with the wall law's
       equilibrium state (u* = 0.35 m/s, z₀ = 0.1 m, k = u*²/C_μ0²), the uncapped geometric length and the
       logarithmic face mean (``Tests/check_closure.py``): ν_t = κ u* (z + z₀) and L = κ (z + z₀) at every cell
       (1e-9), u* recovered (1e-10), R_u and R_v zero below the top cell (1e-10 of u*²/Δz; the top cell loses
       u*²/Δz through the slip top, 1e-8), P_s, ε and R_k = −(P_s − ε) equal to their closed forms from the discrete
       strain (1e-9 of ε; the second cell's P_s/ε = 1.45 is the centred difference of ln over 2Δz), R_θ and R_w zero.
       The plotfile against ``Tests/Gold/closure_column_plt00000`` at relative 1e-10 (a gold written by ablsteady).
       Double precision only.
   * - ``drift_check_selftest``
     - The ``extern_erf`` drift check on a throw-away ERF repository (:doc:`extraction`).
   * - ``refavg_snapshots``, ``refavg_build``, ``refavg_verify``, ``refavg_matches``
     - ``ablsteady_refavg`` end to end: three copies of the smoke plotfile at 0, 5 and 10 s are archived on 2
       ranks with two inputs files; the archive verifies on 3 ranks with the expected metadata, AU and a zero
       half-window change; the average equals the smoke plotfile at zero tolerance. This checks the round trip;
       identical snapshots cannot check the weights, which UT-1.1.2 does.
   * - ``refavg_damage``, ``abort_refavg_damaged``, ``abort_refavg_inputs``, ``abort_refavg_verify_keys``
     - A copy of the archive with one bit flipped is refused, naming the file; bad tool inputs stop the tool with
       each input named; a build key given with ``refavg.verify`` is reported.
   * - ``versions_check``, ``versions_check_selftest``
     - ``external/VERSIONS`` is well formed and its ``amrex`` entry matches the submodule; the self-test shows a
       mismatch with a throw-away submodule is reported.
   * - ``kestrel_parse_*``, ``kestrel_dryrun``
     - The Kestrel scripts parse; a dry run with stub ``module``/``sbatch``/``srun`` checks the Askervein job
       script (run length, plot period, no 2-D plotfiles, ranks) and archives a fake finished run, averaging only
       the 3-D plotfiles in the final window and archiving the inputs but not the run's outputs; a second case
       with the default prefix ``plt`` keeps its ``plt2d_`` plotfiles out.

GoogleTest suite
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - ID
     - Check
   * - UT-1.1.1
     - A reference archive written, verified and loaded: fields bit-identical to the average, metadata and field
       statistics unchanged; a changed byte, a missing file, an extra file, edited AU, and data changed together
       with its checksum are each found and named, and loading refuses.
   * - UT-1.1.2
     - *Beyond the plan.* The average equals a host-side recomputation bit for bit and, for fields linear in
       time, the mid-window value to 8 ε; trapezoid weights for uneven spacing; the half-window change is exact
       for a linear drift and zero for a steady field; a snapshot with a NaN, fewer than 3 snapshots and times out of
       order are refused with their messages; a trailing ``/`` is handled.
   * - UT-1.1.3
     - *Beyond the plan.* AU = nodes × hours × rate (10 CPU, 100 GPU); metadata range checks name their inputs.
   * - UT-1.1.4
     - *Beyond the plan.* SHA-256 against the NIST example vectors, at every padding boundary and chunking;
       ``hex()`` can be called repeatedly.
   * - UT-1.4.1
     - Vertical levels: uniform; solved ratio (top exact, ratio 1.0421553943273616 from an independent exact-rational bisection);
       fixed ratio (ERF's Askervein 20 m × 1.08, top 749.0048747962335 m); z-levels file and its errors.
   * - UT-1.4.2
     - Flat terrain gives the Cartesian levels: bit for bit for STF and Sullivan, within 2 ulp for BTF.
   * - UT-1.4.3
     - Free-stream preservation on round Witch-of-Agnesi and Gaussian hills, all three mappings: divergence below
       1e-12 of :math:`(|U|+|V|)/\Delta x`. The cancellation is algebraic for any node heights, so this test does
       not detect wrong heights; UT-1.4.7 does.
   * - UT-1.4.4
     - Cell volumes sum to the column volumes above the ground to 1e-12, and to a 1 m quadrature of the domain
       above the hill to 1e-4 (the bilinear ground differs from the continuous hill by 2.4e-5 of the domain).
   * - UT-1.4.5
     - Smoothing: a constant is unchanged; a Fourier mode is damped by the analytic factor per pass; β ≥ 1/4 is
       refused. *Beyond the plan:* on a non-periodic domain, edge nodes use themselves as the outside neighbour and
       ghost nodes take the nearest domain node.
   * - UT-1.4.6
     - IDW returns input points and, with 4 neighbours, is exact on a plane at a symmetric query (with the
       default 6 the two extra neighbours break the symmetry and IDW is not exact on a plane); ERF's bilinear form
       is exact on planes anywhere, with ERF's edge fallback; ``terrain.idw_neighbors`` above the number of points is
       refused.
   * - UT-1.4.7
     - ``z_phys``, ``detJ`` and nodal heights equal ERF's own output for the same 2-D hill deck (BTF, STF,
       Sullivan; STF and Sullivan bit for bit on macOS, BTF to rounding because ablsteady keeps the ground layer as
       given, all within 1e-13 elsewhere); STF capped at 2 of 16 levels, as ERF reports.
   * - UT-1.4.8
     - *Beyond the plan.* The IDW bucket search picks the brute-force neighbours; results agree to 1e-13 (random
       cloud with far queries, a regular grid with exact distance ties, a tie across a bucket edge, a transect
       queried from a 2-D domain, one point); UTM-sized coordinates are read and interpolated in double; file terrain on a periodic domain has one
       seam height; every mapping keeps the ground layer as given; ground at the top and inverted cells stop the
       run; the slope of a diagonal ridge is reported in full; a restart with changed terrain is refused; the report
       of a steep ridge carries both warnings.
   * - UT-1.5.1
     - Convection on a manufactured flow over a 120 m cosine hill (BTF, geometric layers, 16 × 12 × 8 to 64 × 48 × 32
       cells), ghost nodes continuing the layers: the unlimited MUSCL reconstruction (van Albada slope with
       :math:`\epsilon` above every difference) is second order for q (conservative and bounded forms) and for the u, v
       and w control volumes, L2 order ≥ 1.9 on both pairs and max-norm order ≥ 1.9 on the finest. With ``van_leer``
       and ``van_albada``: L1 order ≥ 1.5 on the finest pair and an L1 error at most a quarter of first-order upwind's
       (a limiter is first order at smooth extrema). First-order upwind converges at order 1.0 ± 0.1.
   * - UT-1.5.2
     - Scalar diffusion with a variable :math:`\gamma`, and momentum diffusion with the full stress tensor and a
       variable :math:`\nu`, on the same meshes: L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on the finest,
       for q, u, v and w.
   * - UT-1.5.3
     - Random :math:`k \ge 0` (a third of the points zero) and random transport fluxes of ±50 m/s: the step
       :math:`k - C(k)/a_P` with the bounded upwind form leaves every :math:`k \ge 0`; the conservative upwind form
       goes below −0.1 in the same step. *Beyond the plan:* the bounded form of a constant is zero for every limiter
       and location.
   * - UT-1.5.4
     - Random transport (±20 m/s) and diffusivity (0 to 50 m²/s, a fifth zero) over the hill: every neighbour
       weight of :math:`A_1` is ≥ 0 and :math:`a_P` equals their sum (1e-14; single 1e-6), for cells and the three
       velocity components.
   * - UT-1.5.5
     - The van Albada slope's derivative changes by less than 0.1 across every point where van Leer switches
       (:math:`a = 0`, :math:`b = 0`, :math:`a = -b`, the origin), in three directions; van Leer's jumps by more than 1.
       On the manufactured flow (64 × 48 × 32) the two limiters differ by less than 1 % (L2) for q, u, v and w.
       Invalid ``smooth_eps`` (0, negative, NaN, infinite) is refused.
   * - UT-1.5.6
     - *Beyond the plan.* :math:`A_1` is the explicit operators' first-order part: bounded upwind convection equals
       :math:`A_1 q` for random q and transport (1e-12; single 1e-5), and on a flat mesh :math:`-D(q) = A_1 q` for the
       scalar and :math:`-D_c(u) = A_1 u_c` for each velocity component with random :math:`\nu`. Over the hill, every weight
       equals the explicit face flux of a field that varies along one index direction only (this pins the
       :math:`1 + \delta_{cd}` and :math:`|\mathbf{n}|^2 + n_c^2` factors on sloping faces). MUSCL convection,
       both diffusions and :math:`A_1` agree to round-off (1e-13) on 16² and 8² column boxes, and boxes that share a
       face point agree there (not bit for bit: GCC on arm64 contracts multiply-adds inside vectorised loops only).
   * - UT-1.6.1
     - The divergence (through :math:`\Omega`) and the three components of :math:`G p` over the 120 m hill (BTF,
       geometric layers, 16 × 12 × 8 to 64 × 48 × 32 cells): L2 order ≥ 1.9 on both pairs, max-norm order ≥ 1.9 on
       the finest (single precision 1.8). :math:`\langle Gp, U\rangle = -\langle p, DU\rangle` with control-volume
       weights to 1e-12 of either side (7e-12 against 3.3e4), for random p and a random velocity that vanishes near
       the ground and top. *Beyond the plan:* :math:`G` does not read the ghost cells below the ground.
   * - UT-1.6.2
     - The continuity residual is bit for bit ``face_divergence`` of :math:`(u, v, \Omega)`; after a projection
       (Poisson fluxes added to :math:`u, v, \Omega`, w from ``w_from_omega``) the divergence falls from 1.19 to
       2.4e-11, the solve's residual. *Beyond the plan:* ``omega_from_w`` and ``w_from_omega`` are inverses.
   * - UT-1.6.3
     - ``line_jacobian`` for u, v, w, θ and k in x, y and z against a finite-difference Jacobian of the first-order
       operator (transport frozen), per coefficient: relative difference < 1e-6 (1e-13 measured; single 1e-3), over
       the hill with convection, and on a flat mesh with convection and a random diffusivity.
   * - UT-1.6.4
     - Rest over a 300 m hill (slopes up to 0.94), θ = θ_ref(z) stratified, Coriolis on: the momentum, continuity and
       k residuals are below 1e-14. With a constant reference the hydrostatic :math:`p(z)` balances only to
       truncation (R_u 1.6e-3 m/s², 0.5 % of the buoyancy; R_w exact on interior faces).
   * - UT-1.6.5
     - ``jvp`` on a laterally uniform column (affine residual) equals :math:`R(q+v) - R(q)` to 1e-5 (single 1e-2).
       *Beyond the plan:* on the hill with MUSCL convection, against a central difference, L1 error below 1e-5
       (2.7e-6) while a step 1e5 times larger misses (1.8e-4); the result is the same without ``R0``; the step is the
       same on 16² and 8² column boxes and scales with each field (k and its scale times 100 leave it unchanged).
   * - UT-1.6.6
     - *Beyond the plan.* Sources in closed form: Coriolis with the :math:`\cos\varphi` terms and geostrophic forcing
       on a uniform wind (:math:`f = 7.27221\times10^{-5}` 1/s at 30°, T = 86400 s), buoyancy of a uniform excess;
       the geostrophic wind is balanced. The residual equals its operators applied one by one (momentum, continuity,
       θ, k) on the hill with a general state. Invalid options (latitude, period, geostrophic forcing at latitude 0 or with
       the forcing off, θ₀, gravity, profiles that do not increase) and missing inputs stop with a message.
   * - UT-1.4.9
     - *Beyond the plan.* Node heights shared by two boxes (or periodic images) are equal, for BTF, STF and Sullivan,
       periodic and not (under gcc each box's copy differed by up to 1e-13 m before ``Mesh`` synchronised them).
   * - UT-1.7.1
     - ERF's similarity functions (the copied ``ERF_MOSTUtils.H``, byte-identical by the drift check): Businger–Dyer
       ψ and φ and Jimenez ψ against a Python table of the same formulas, z/L from −5 to 5, to 1e-14 (single 1e-6);
       :math:`d\psi/d\zeta = (1-\phi)/\zeta`.
   * - UT-1.7.2
     - Reference-height lookup on the ERF#3982 column (dz₀ = 10 m, ratio 1.1, 40 layers): centres belong to the interval
       above, faces (10 m included) to the interval that holds them, the last centre to the last interval; heights below
       the first or above the last centre are refused.
   * - UT-1.7.3
     - u* from a log-law state (u* = 0.45, z₀ = 0.1 m) on that column: within 2e-5 at the first centre, at 10 m, at 55.3 m
       and on a centre. *Beyond the plan:* with a heat flux, u* matches an independent bisection to 1e-10 for unstable,
       stable and low-wind convective cases, and is positive.
   * - UT-1.7.4
     - :math:`k_w = u_*^2/C_{\mu0}^2` exactly (neutral); Axell's buoyant form with a heat flux; the wall eddy viscosity
       equals :math:`C_{\mu0}\sqrt{k_w}\,\kappa(d+z_0)`.
   * - UT-1.7.5
     - *Beyond the plan.* Over the hill: the columns' distances and areas from their nodes; the wall stress along the
       ground with magnitude :math:`u_*^2|A|`. In the residual, R_u, R_v and R_θ change only in the first cells, by the
       wall flux minus the interior flux over the cell volume; the k rows are the wall-k rows (both modes). The
       ``line_jacobian`` rows lose the ground neighbour and gain the wall terms. At low wind the wall-k row uses
       u* ≥ 0.01 m/s. Invalid options and a reference height outside the columns stop with a message.
   * - UT-1.8.1
     - The Poisson boundary types of every face type on every face it may take, periodic directions included, equal
       ERF's ``get_fft_bc`` (the copied ``ERF_SolverUtils.H``).
   * - UT-1.8.2
     - On a ridge whose height varies along the inflow face: every inflow ghost point (both layers, u, v, θ, k) takes the
       profile at its own height above the ground, the inlet wind varies along the face by more than 0.05 m/s, and the
       outflow face points of u are not Dirichlet rows. *Beyond the plan:* the ghost values of every rule for every
       face type (even, odd, value, profile, gradient, the Dirichlet face points, the top's θ gradient); the analytic
       profiles and their equilibrium k and :math:`\nu_t`, and a file table with and without a k column.
   * - UT-1.8.2b
     - *Beyond the plan.* Over a 120 m hill with a stable heat flux, the ground ghosts of u, v and θ give MOST's
       gradient at the first centre (u* φ_m/(κ(d₀ + z₀)) along the stress, −q φ_h/(κ u* (d₀ + z₀))) in both layers,
       on the two lateral ghost columns too (the periodic images of the neighbouring boxes; this catches a gradient
       slab with too few ghost columns under ``init_snan``); k is even; the homogeneous fill gives zero gradients
       and keeps even rules. The implicit folds: :math:`A_1` on the interior values with the folds equals the unfolded
       :math:`A_1` on the homogeneously filled state, for u, v, w, θ and k, every lateral type and a geostrophic top
       (1e-13; single 1e-5); Dirichlet face points carry identity rows. The residual's Dirichlet rows are
       :math:`q - V`; ``line_jacobian`` gives the identity there. The Rayleigh and sponge coefficients have their
       shapes (the :math:`\sin^2` ramp, zero below :math:`z_d`; the quadratic sponge, zero at its inner edge and on
       faces that are not inflow or outflow) in the residual and on the implicit diagonal. Every ``bc`` key is read
       into its field, and each out-of-range, ineffective (the other kinds' inflow keys, surface keys without a
       surface-layer ground, a layer without a coefficient) or inconsistent key is named.
   * - UT-1.9.1
     - ``closure_point`` equals ERF's composition of the copied AL01 relations (``ComputeTurbulentViscosityRANS``,
       ρ = 1) bit for bit on 20 000 random (k, N², d): L, ε, R_t, c_μ, c_μ', ν_t, ν_t', across the neutral, stable,
       unstable, smoothed and floored cases (each exercised).
   * - UT-1.9.2
     - Kynema's ``KLAxell`` (transcribed) where the formulations coincide (neutral, z = d, λ = 30 m, no z₀, ρ = 1):
       L, ν_t, ε and P_s to 1e-12.
   * - UT-1.9.3
     - ``smooth_Rt``: identity above R_t^crit, continuous at it, within (R_t^min, R_t^crit] and monotone down to
       R_t = −1e30; the unstable length finite and below its bound for k down to 1e-16 under N² to −1e3; the same
       decision as ERF's.
   * - UT-1.9.4
     - The point-implicit update k + S_k/(1/Δτ + ∂S_k/∂k) with the Patankar diagonal stays ≥ 0 (> 0 with production)
       for 1/Δτ = 1, 1e-3 and 0 on 5000 random states; without the diagonal a decaying state goes negative.
   * - UT-1.9.5
     - The dual-number dS_k/dk and dν_t/dk against central finite differences on 3000 random states (branches held):
       relative 1e-6 (single 1e-2); the dual value part equals the Real evaluation to a few ulp.
   * - UT-1.9.6
     - Every ``turb`` key is read into its field and the constants follow (the floor, C_μ0³); each invalid value is
       refused with its key named, and the constructor refuses it before forming the constants.
   * - UT-1.9.6b
     - *Beyond the plan.* The logarithmic mean against (a − b)/ln(a/b) on random pairs (1e-14), equal values, zero,
       continuity across its series switch; on a log-law column the logarithmic face mean gives the wall flux u*² at
       every interior u face (1e-12) where the arithmetic mean misses the first by 2–12 %; the implicit operator with
       the same mean satisfies A₁u = −D(u). Over a 120 m hill the residual with the closure equals the residual with
       the closure's own fields given minus the sources bit for bit, the k line Jacobian's diagonal differs by dSdk
       exactly, every field is finite with positive diffusivities, misuse (fields given with a closure, no wall
       distance) is refused, and ``height_above_ground`` is the node-mean height above the column's ground.
       Relaxation: ν_t = 0.6 new + 0.4 stored exactly, unrelaxed without history, and the history survives a
       checkpoint round trip onto other boxes (identical relaxed evaluation); frozen switches hold the floor decision
       after k drops below it.
   * - UT-1.10.1
     - Manufactured :math:`\phi` on a stretched STF mesh over a 60 m hill: every lateral BC (periodic, Neumann,
       Dirichlet) with a Dirichlet or Neumann top, 2-D (y periodic) at three resolutions, and a 3-D case (Neumann in x,
       Dirichlet in y, two resolutions): order ≥ 1.9 (single precision: 1.8 coarse, 1.6 fine), residual below 10 ×
       tolerance. The solution has zero gradient at the ground and level side faces, where ERF's operator is exact
       for the mirror ghost cells (next row).
   * - UT-1.10.1b
     - *Beyond the plan.* With a tangential gradient at sloped ground: order 2 on flat ground (≥ 1.95), between 0.8 and
       1.2 over a 120 m hill (64 → 128 → 256 cells), and a fine-mesh error 3 to 6 times that over a 60 m hill (slope²).
       Double precision only. Pins ERF's ground treatment (:doc:`poisson`).
   * - UT-1.10.2
     - The anisotropic operator with :math:`\beta_v/\beta_h` = 1, 0.5, 0.1 (``z_cross`` set) on flat ground and over
       the hill, and a spatially varying :math:`\beta` over the hill: order ≥ 1.9 (single: 1.8, 1.6). Without
       ``z_cross`` the hill error is more than 5 times larger.
   * - UT-1.10.3
     - :math:`\beta = 1` gives the no-coefficient path bit for bit (iterations, solution, fluxes); :math:`\beta = 2.5`
       on every face gives :math:`\phi/2.5` and the same fluxes (to 1e-9; single 1e-4). Both paths are this layer's;
       there is no separate plain projection to compare with.
   * - UT-1.10.4
     - On a ramp (slope 0.15) crossing the Neumann faces, their fluxes and the ground's are exactly zero, and
       :math:`\sum J\,A\phi` equals the flux out of the Dirichlet top (1e-12; single 1e-4).
   * - UT-1.10.4b
     - *Beyond the plan.* Operator and fluxes agree cell by cell, :math:`A\phi = -\nabla_h\cdot F` (1e-12 of
       :math:`\max|A\phi|`; single 1e-4), with varying coefficients and ``z_cross`` on the ramp, Dirichlet and Neumann
       top.
   * - UT-1.10.5
     - Projection: the divergence of :math:`u + F` is below 10 × tolerance of the divergence of :math:`u`, and equals
       the residual of the solve.
   * - UT-1.10.6
     - Symmetric hill and right-hand side: :math:`\phi` mirror-symmetric to 1e-12 of :math:`\max|\phi|` (single 1e-5),
       Neumann and Dirichlet sides.
   * - UT-1.10.7
     - *Beyond the plan.* Preconditioner, with bounds 10 % above the counts measured: 19 iterations isotropic (more than
       10 times as many without it); 16 for :math:`\beta_v/\beta_h = 0.1` told the ratio (at most 3 more than
       isotropic; not told, more than 1.5 times as many); 29 for :math:`\beta` varying over 100 times. A singular
       problem reports and removes the mean. Option and BC checks, coefficients that are not positive, MultiFabs on
       another layout and a solve that does not converge each stop with the input named.
       ``tools/poisson_ms_check.py`` (sympy, not run by CI) checks the manufactured solutions' boundary conditions.
   * - UT-1.11.1
     - Back end B (``mlmg`` and ``gmres_mlmg``): every lateral BC × Dirichlet/Neumann top × Neumann (zero-gradient
       solution) / Dirichlet ground, 2-D (y hidden) at three resolutions (two in single precision), and a 3-D case:
       order ≥ 1.9 (single: 1.8, 1.6).
   * - UT-1.11.1b
     - *Beyond the plan.* Back end A with a Dirichlet ground: order ≥ 1.9 for every lateral BC and top.
   * - UT-1.11.2
     - Spatially varying coefficients, Neumann and Dirichlet ground, 3-D with non-periodic sides (Neumann x,
       Dirichlet y and top): both back ends solved to 1e-12 give :math:`\phi` within 1e-9 of max :math:`|\phi|` (single:
       1e-5 and 1e-4), and the fluxes satisfy :math:`\nabla_h\cdot F = -\mathrm{rhs}` to the residual.
   * - UT-1.11.2b
     - *Beyond the plan.* ``z_cross`` with back end B, ``z_cross`` with a Dirichlet ground, ``apply`` with back end B,
       an unknown back end, ``solver.poisson_precond`` with back end B, ``solver.poisson_restart`` with ``mlmg`` and a
       non-converging MLMG solve each stop with the input named.
   * - UT-1.11.3
     - Tucker distance against the exact distance with ERF's tolerances on ERF's decks: flat ground (max relative
       ≤ 3e-3, ≤ 0.1 m below 100 m; single 1e-2 and 0.5 m), 2-D ridge and 3-D hill (max relative ≤ 0.15, mean ≤ 0.03,
       ≤ 0.2 cells below 100 m); the 3-D hill with periodic and Neumann sides; both back ends (back end B on 64-cell
       boxes, under CI's trap). Single precision solves to 1e-4. *Beyond the plan:* on flat ground the error
       converges at second order (32, 64, 128 layers; double precision).
   * - UT-1.11.4
     - *Beyond the plan.* No node of ``z_phys_nd`` (ghosts included) holds ERF's ``bogus_large_value``, for BTF, STF and
       Sullivan with periodic and non-periodic sides; ``fill_ghosts`` fills every ghost cell, corners and between-box
       ghosts included, with odd copies below a Dirichlet ground and even copies at Neumann faces.
   * - UT-1.12.1
     - :math:`\alpha_h = \alpha_v = 1`: the initial field equals the Poisson layer's plain projection of the same first
       guess (the solve without coefficients, applied the same way) bit for bit, :math:`\lambda` included;
       :math:`\alpha_h = \alpha_v = 2` gives the same velocity within 1e-9 of max :math:`|u|` and :math:`\lambda/4`
       (single: 1e-4); the scalar state is untouched by the correction. *Beyond the plan:* scaling the whole z flux by
       :math:`\alpha_v^2` instead of the cross part by :math:`\alpha_h^2` gives a different field; with periodic sides
       the hill case is divergence-free as a singular problem and the seam ghost of :math:`u` is the corrected
       periodic image; a :math:`\lambda` handed in with garbage is zeroed first; the divergence after is measured,
       not assumed.
   * - UT-1.12.2
     - The divergence after the correction is below 10 × ``solver.poisson_reltol`` times the divergence before, in the
       solver's own operator (``velocity_divergence``), in every initializer test.
   * - UT-1.12.3
     - Flat terrain, open and periodic sides: the profile is divergence-free, so :math:`\lambda`, :math:`w` and the
       divergence stay below 1e-12 and :math:`u` equals the log law at every x face's height (single: 1e-4, 1e-5).
   * - UT-1.12.4
     - D4 symmetry: on a square domain with a Gaussian hill, the wind along ±x and ±y (its own inflow and outflow
       faces each time) gives fields that are the eight transforms of each other: the cell-centred velocity vector
       maps with the transform's matrix, :math:`w` and :math:`\lambda` as scalars, within 1e-12 of the maxima (single
       1e-5), with the solve converged to round-off (``reltol`` 3e-15).
   * - UT-1.12.5
     - The speed-up above the summit (four central cells a quarter of the depth up, over the profile at their
       height) falls by more than 0.01 at each step of :math:`\alpha_v` = 1, 0.7, 0.5, 0.3, 0.1 (1.015, 0.998,
       0.975, 0.943, 0.913 on the 1 km test domain, whose open faces are 500 m from the summit), and
       :math:`\alpha_h = 2, \alpha_v = 1` gives the same value as :math:`\alpha_v = 0.5` (the
       ratio alone matters; 1e-9). The preconditioner is told the ratio: at most 21 GMRES iterations at every ratio
       (13 to 19 measured).
   * - UT-1.12.6
     - The first guess over sloping terrain: :math:`u` at every x face, the two ghost layers included, and :math:`v`
       at every y face equal the profile at the face's height above its column's ground with the 3 m floor (the
       floor acts on the ground faces and the ghosts below); :math:`\theta` equals the sounding by absolute height;
       :math:`k` the profile's value with the decay above 200 m, and :math:`u_*^2/C_{\mu 0}^2` everywhere without it;
       :math:`p' = 0`; :math:`w_0 = 0` except on the ground faces, which carry no flux (:math:`\Omega_0 = 0` to
       round-off); every ghost cell of the state is written (no NaN from the allocation); the O'Brien pass refuses a
       layout whose boxes do not span the column.
   * - UT-1.12.7
     - *Beyond the plan.* O'Brien: :math:`\Omega_0` vanishes at the ground and the top exactly, the divergence left in
       every cell is the weight difference :math:`f_{k+1}^2 a_{z,k+1} - f_k^2 a_{z,k}` (:math:`f` the height fraction,
       which the test shows differs from the index fraction on the stretched mesh) times the column's residual
       (1e-11), the column sum equals the residual flux, the residual the function returns equals the largest column's;
       through the initializer the correction starts from a different divergence than without O'Brien (not a lower
       one: the terrain-following first guess is already consistent at the ground) and gives a different :math:`w`.
   * - UT-1.12.8
     - *Beyond the plan.* Inputs: :math:`z_{agl,\min} < 0` or NaN, :math:`\alpha_h \le 0` or infinite,
       :math:`\alpha_v < 0` or NaN, :math:`z_i < 0`, and any key other than ``first_guess`` given without
       ``first_guess = profile`` are refused with the key named, by the checks and by both entry points;
       ``first_guess = none`` leaves the state untouched; the reader takes every key and refuses an unknown
       ``first_guess``, ``obrien = 2``, a repeated key and the keys without the first guess; the Poisson types of the
       correction are Dirichlet on inflow and outflow faces, Neumann elsewhere, periodic where the mesh is; a face set
       that does not fit the mesh is refused.
   * - UT-1.12.9
     - *Beyond the plan* (``Poisson.SplitOperatorKeepsTheMirrorSymmetries``). The split anisotropic operator
       (``z_cross`` set, :math:`\beta_v = 0.36`) on a Gaussian hill varying in x and y keeps the mesh's x and y mirror
       symmetries within 1e-12 of max :math:`|\phi|`; the cross kernel's first transcription broke the y symmetry at
       the 1e-4 level.
   * - UT-1.13.1
     - The local pseudo-time step on known fields: uniform :math:`(u, v, w)` and :math:`\nu` on a flat uniform mesh
       give :math:`\mathrm{CFL} / (|u|/\Delta x + |v|/\Delta y + |w|/\Delta z + 2\nu(1/\Delta x^2 + 1/\Delta y^2))` in
       every cell to 1e-14 (single 1e-6; the hill cases 1e-12 and 1e-5), without :math:`\nu` the diffusive term is gone, a field at rest without
       diffusion gets the floored step, :math:`u` varying in :math:`x` gives the face mean per cell, a
       terrain-following field (:math:`\Omega = 0`) over a 120 m hill has no vertical term and a plain :math:`w` over
       it the term :math:`|w| / (J \Delta\zeta)` cell by cell; a non-positive or NaN CFL and a face layout are
       refused.
   * - UT-1.13.2
     - The SER ramp on a geometric residual (rate 1/4): CFL 1, 1, 1, 4, 16, 64, 64, 64 (``cfl0`` through the two
       reference iterations, then the law capped by the growth 4 and ``cfl_max`` 64), :math:`p = 2` without caps
       256, 4096, 65536, a growth cap of 1.5 binding every iteration, :math:`p = 0` holding the CFL, a residual of
       exactly 0 reaching the cap through the growth cap; the stage's :math:`\Delta\tau` follows the CFL.
   * - UT-1.13.3
     - Switch rules on scripted residuals (:math:`W` = 5, threshold 0.01, back growth 10): a residual alternating
       about the threshold never switches; five iterations below switch forward at iteration 5 with the reason
       naming ``steady.switch_residual``, the residuals and the stage's ``init`` recorded; a 20× growth of the
       second stage's residual switches back at iteration 11 (the window compared must lie within the second stage,
       entered at 6); every pair of switches is at
       least :math:`W` apart; the records carry each iteration's stage; one stage never switches.
   * - UT-1.13.4
     - The history written to ``<dir>/Driver`` and read back is equal field by field (records, switches with their
       reasons, accounts, the criterion's history); a two-stage run interrupted after 6 iterations and continued by a
       new driver from the file on the saved state gives the same records (stage, CFL, relative residual, residuals,
       monitors), switches and verdict iteration as the uninterrupted run; the AU of Kestrel CPU on 2 nodes equal
       20 × wall hours; the history is identical on every rank; a truncated file, a missing file and a history for
       another number of stages are refused.
   * - UT-1.13.5
     - *Beyond the plan.* A residual that becomes infinite or NaN stops the run as ``Diverged`` with no arithmetic on
       it (the traps); ``max_iterations`` bounds
       a run that never settles and a second ``run`` does nothing; no stages, a null stage and stages that disagree
       on residual names are refused (a null first stage too); ``restore`` refuses a NaN CFL, a switch state or
       iteration count that does not fit the records, and brings the CFL within changed bounds; every ``steady.*``
       range check names its key; the reader takes every key,
       its defaults, and refuses a repeated key, an unknown machine and a bad value.
   * - UT-1.14.1
     - :math:`a_P` on every x-, y- and z-face equals :math:`1/\Delta\tau` plus the diagonal of ``line_jacobian``
       exactly, with a random :math:`\Delta\tau` field over a 60 m hill with random diffusivities (both diagonals
       come from ``implicit_rows``, whose coefficients UT-1.6.x pin: the statement is the :math:`1/\Delta\tau`
       addition and the face step); the diagonal exceeds :math:`1/\Delta\tau`; the two copies of the periodic seam
       face have identical rows and identical corrections; the pressure and scalars are untouched and :math:`u`
       changes by exactly the z-line solve's correction.
   * - UT-1.14.2
     - The z-line solve (the vertical upwind weights of :math:`A_1` included) against a dense direct solve of the same
       system,
       assembled independently from the coefficient fields: on a single column (1 × 1 × 12, periodic sides, the
       lagged neighbours are the column itself) to 1e-12 after 60 sweeps (single 1e-4); the ghosts outside the domain
       stay 0.
   * - UT-1.14.3
     - :math:`\alpha_u` relaxation equals pseudo-time relaxation: with :math:`\Delta\tau = \alpha_u/(1-\alpha_u)/
       \mathrm{diag}(A_1)` on every face, :math:`a_P = \mathrm{diag}(A_1)/\alpha_u` to 1e-12 (single 1e-5) and the
       converged correction of u, v and w equals the dense solve of the relaxed system to 1e-12 (single 1e-4).
   * - UT-1.14.4
     - z-line Jacobi on a 3-D linear advection-diffusion problem (4 × 4 × 6, periodic, :math:`\Delta\tau` = 20 s)
       converges to the dense direct solution: the error falls at every doubling of the sweeps (1, 2, 4, 8, 16) and
       is below 1e-10 after 400 sweeps, with the last sweep's change below 1e-10.
   * - UT-1.14.5
     - Fixed point: a uniform wind over flat ground with uniform pressure has a residual of exactly zero and the
       predictor changes nothing (:math:`|\Delta q| = 0`; with a zero right-hand side the solve is zero for any
       coefficients, so this follows from the residual); a perturbation of :math:`u` at one face is pulled back: one
       predictor step with :math:`\Delta\tau` = 10⁴ s and 50 sweeps lowers max :math:`|R_u|` below 0.75 of its value
       (0.56 measured; the pressure and the transport are frozen in the step).
   * - UT-1.14.6
     - *Beyond the plan.* ``face_dtau`` is the mean of the two cells, the one cell at a wall, the ground and the top,
       and the periodic image across the seam; one sweep's change is the whole correction; ``sweeps < 1``, a
       non-positive step and wrong layouts are refused.
   * - UT-1.14.7
     - *Beyond the plan.* With boundary conditions (inflow at xlo, outflow at xhi, slip ground and top), the Dirichlet
       rows (the inflow face of u, the ground and top faces of w) have :math:`a_P = 1`, :math:`d = 0`, and one
       predictor step sets them to the boundary value exactly for any :math:`\Delta\tau` (the inflow faces then have
       a zero residual row); every other row has :math:`a_P = 1/\Delta\tau + \mathrm{diag}(A_1)` and :math:`d = a_P -
       \sum a_{nb} > 0`.
   * - UT-1.14.8
     - *Beyond the plan.* The same state on one box and on 2 × 2 boxes (the lagged neighbours cross the box seams
       through the ghost exchange) gives the same corrected u, v, w, :math:`a_P`, :math:`d` and face step to 1e-12
       (single 1e-5).
   * - UT-1.2.1
     - Steady-criterion verdicts on converging, stalling (including a residual stuck exactly at the tolerance) and
       oscillating histories. The converging history becomes
       steady at iteration 136, pinned from an independent Python implementation of the definition.
   * - UT-1.2.2
     - *Beyond the plan.* A NaN or infinite residual or monitor gives ``Diverged``, permanently.
   * - UT-1.2.3
     - *Beyond the plan.* A monitor at round-off level around zero is steady with scale 1 and never steady with a
       scale below the noise.
   * - UT-1.2.4
     - *Beyond the plan.* The monitor window needs :math:`W+1` values; a change of exactly the tolerance passes and
       one just above it does not (values exact in binary).
   * - UT-1.2.5
     - *Beyond the plan.* The residual drop is measured from the early peak (iteration 38), not from the first
       value (iteration 49).
   * - UT-1.2.6
     - *Beyond the plan.* Parameter range checks name their keys; bad or missing monitor scales and size
       mismatches stop.
   * - UT-1.2.7
     - *Beyond the plan.* A monitor converging at rate 0.999 stops at iteration 4602, within the tolerance of its
       limit; the window test alone would stop at 1615 with a 2 % error.
   * - UT-1.2.8
     - *Beyond the plan.* Undamped oscillations with periods of 37, 100 and 200 iterations are never steady.
   * - UT-1.2.9
     - *Beyond the plan.* A residual that is zero over the reference iterations does not block the stop.
   * - UT-1.2.10
     - *Beyond the plan.* A criterion restored from a saved history gives the same verdict at every later
       iteration as one never interrupted; a mismatched history is refused.
   * - UT-1.2.11
     - *Beyond the plan.* A slowly growing monitor is not steady while its change over the window is still small.
   * - UT-1.3.1
     - Input parsing; unknown, repeated, multi-valued and out-of-range keys reported with their names.
   * - UT-1.3.2
     - Plotfile round trip, bit-identical, including the nodal heights and their ghost nodes.
   * - UT-1.3.3
     - Checkpoint and restart onto another box layout, bit-identical; ghost nodes restored without overwriting
       valid nodes; a checkpoint written with one ghost cell on the solution fields (before sub-phase 1.5)
       restarts bit-identical.
   * - UT-1.3.4
     - The sum gives the same bits on 1, 2 and 4 ranks and within 16 ε Σ\|x\| of the exact sum; the maximum
       handles NaN and −0.
