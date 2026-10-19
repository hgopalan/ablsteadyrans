Mass-consistent initializer
===========================

The initial field of a run (sub-phase 1.12) is a mass-consistent, WindNinja-style wind: a first guess from the inflow
profile, made divergence-free by a variational correction with an anisotropy between the horizontal and the vertical
adjustment. The method follows `FastWindTerrain <https://github.com/hgopalan/FastWindTerrain>`__ (``external/VERSIONS``
pins the commit): no code is copied, since FastWindTerrain solves on a Cartesian immersed-boundary mesh and ablsteady on
the terrain-following one, but the meaning of every parameter is the same. ``Source/Initializer.H`` holds the code;
``init.first_guess = profile`` turns it on (:doc:`inputs`).

Method
------

**First guess.** :math:`u_0, v_0` at every x- and y-face point are the inflow profile ``bc.inflow_*``
(:doc:`boundary`) evaluated at the point's height above its column's ground,

.. math::

   h = \max\big(z_{\mathrm{face}} - z_{\mathrm{ground}},\; h_{\min}\big),

with :math:`z_{\mathrm{face}}` the node-mean height of the face and :math:`z_{\mathrm{ground}}` the node-mean ground under
its footprint (``Source/Heights_K.H``, the same helper as the boundary fills), and :math:`h_{\min}` =
``init.z_agl_min``. The ghost cells of :math:`u_0, v_0` are filled the same way (lateral ghosts from the stored ghost
nodes, vertical ghosts from the extrapolated slab). :math:`w_0 = 0`, except on the ground faces, which carry no
flux (:math:`\Omega_0 = 0` there, :math:`w_0 = u_0\,\partial_x z + v_0\,\partial_y z`): the correction's flux
through the ground is zero, so the no-penetration constraint of Sherman's problem (the minimisation of
:math:`\|u - u_0\|_\Sigma` over the interior faces with the ground flux fixed) enters through the first guess; a first
guess with a flux through the sloping ground would keep it. The potential temperature is the sounding
``bc.theta_profile`` by absolute height (the reference :math:`\theta_{\mathrm{ref}}(z)` itself); the turbulence energy
is the profile's equilibrium value :math:`k(h)` (for the log law :math:`u_*^2/C_{\mu 0}^2`) times a decay above
:math:`z_i` = ``init.k_decay_height``,

.. math::

   k = k_{\mathrm{prof}}(h)\, e^{-(h - z_i)/z_i} \quad (h > z_i),\qquad k = k_{\mathrm{prof}}(h) \quad (h \le z_i),

no decay when :math:`z_i = 0`. The pressure perturbation is :math:`p' = 0`: the hydrostatic :math:`p'(z) = -\int_z^{\mathrm{top}}
g\,(\theta - \theta_{\mathrm{ref}})/\theta_0\,dz` vanishes because :math:`\theta = \theta_{\mathrm{ref}}`; its
integration comes with a first guess whose :math:`\theta` differs from the sounding.

**O'Brien column correction** (``init.obrien = 1``, off by default as in FastWindTerrain). Before the projection,
:math:`\Omega_0` on the z-faces follows from the column continuity of :math:`(a_x u_0, a_y v_0)` from the ground up,

.. math::

   a_z\Omega_{k+1/2} = a_z\Omega_{k-1/2} - J_k D^h_k\,\Delta\zeta,\qquad \Omega_{1/2} = 0,

with :math:`D^h` the horizontal part of the divergence (``face_divergence`` of :math:`(u_0, v_0, 0)`). The top value
:math:`E = \Omega_{\mathrm{top}}` is the column's residual; it is removed with the weight :math:`f_k^2`, :math:`f_k =
(z_k - z_{\mathrm{ground}})/(z_{\mathrm{top}} - z_{\mathrm{ground}})` the face's height above the ground over the
column's depth (O'Brien 1970, quadratic in height; FastWindTerrain's ``Obrien.cpp`` weights by the face index, the same
on its uniform layers but not on a stretched mesh), so :math:`\Omega_0` vanishes at the top exactly, is the kinematic
value near the ground, and the divergence left in cell :math:`k` is :math:`-(f_{k+1}^2 a_{z,k+1} - f_k^2 a_{z,k})\,E/(J_k
\Delta\zeta)`, concentrated where the layers are thick. :math:`w_0` follows from :math:`\Omega_0` (``w_from_omega``).
The run prints :math:`\max|E|`. O'Brien is meant for a first guess from soundings or fields; with the profile by height
above the ground, :math:`u_0` is horizontally divergent through the whole column over terrain (the profile is shifted by
the ground height), so the kinematic :math:`\Omega_0` grows with height and the corrected field keeps a :math:`w` aloft
that the projection alone would not produce. It is off by default.

**Correction.** With :math:`\Sigma = \mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2)`,

.. math::

   \nabla\cdot(\Sigma\nabla\lambda) = \nabla\cdot u_0, \qquad u = u_0 - \Sigma\nabla\lambda,

solved by the Poisson layer (:doc:`poisson`) with ``FaceCoeffs`` :math:`x = y = z_{\mathrm{cross}} = \alpha_h^2`,
:math:`z = \alpha_v^2` and ``precond_ratio`` :math:`= \alpha_v^2/\alpha_h^2`; :math:`\lambda = 0` on the inflow and
outflow faces, no flux through the walls, the ground and the top, periodic where the mesh is. The fluxes :math:`F =
-\Sigma\nabla\lambda` the solver returns are added to :math:`u, v` and to :math:`\Omega` (the transport fluxes are
:math:`\Omega`-based), and :math:`w` follows from :math:`\Omega`. With :math:`\alpha_h = \alpha_v = 1` no coefficients
are passed and the correction is the plain projection of the Poisson layer, bit for bit. :math:`\alpha` is a
transmissivity, not a penalty: a smaller :math:`\alpha_v/\alpha_h` means less vertical adjustment, so the flow goes
around steep terrain rather than over it, and only the ratio matters for the velocity (:math:`\lambda` scales with
:math:`1/\alpha^2`). In Sherman's (1978) notation the :math:`\alpha_i^2` are penalties, :math:`u = u_0 +
\partial_x\lambda/(2\alpha_1^2)`: ablsteady's :math:`\alpha_h^2, \alpha_v^2` are :math:`1/(2\alpha_1^2),
1/(2\alpha_2^2)`, and :math:`\alpha_v/\alpha_h` here is Sherman's :math:`\alpha_1/\alpha_2` (WindNinja's stability
parameter: small for stable flow, around the hill). The run prints the Poisson iterations and the 2-norm of the
divergence before and after. GMRES iterations grow as :math:`\alpha_v/\alpha_h` falls (13 at 1, 105 at 0.1 on the
bump case; 0.05 does not converge in the default ``solver.poisson_max_iter`` = 200): the FFT preconditioner is told the
ratio but not the terrain cross term, whose vertical coefficient :math:`\alpha_h^2(h_x^2 + h_y^2)` dominates on slopes
when :math:`\alpha_v` is small. Ratios of 0.1 and above are practical.

What the field is, and is not
-----------------------------

The result is divergence-free in the solver's own operator (``velocity_divergence``) to the solver's tolerance, has
the profile on the inflow faces up to the correction's normal component there (:math:`\lambda = 0` on an open face
fixes the potential, not the normal velocity), and :math:`\theta`, :math:`k`, :math:`p'` consistent with the sounding.
It is a starting point for the steady solver, not a RANS solution: there is no momentum balance in it. The boundary
fill (``BoundaryConditions::fill_ghosts``) is not applied here; the solver does that at its first iteration. Field
inputs for the first guess (WindNinja or NWP soundings, FastWindTerrain's own output) are not read in this version.

Comparison with FastWindTerrain
-------------------------------

RT-1.12.2 (``Tests/check_init.py trend``) runs both codes on FastWindTerrain's bump regression case
(``regtests/phase6_solve_correction/inputs_bump``: a 100 m Gaussian hill on a plane sloping 0.1 in x over 1000 m ×
1000 m, a power-law inflow of (8, 6) m/s at 10 m with exponent 0.14, :math:`\alpha_h = 1`) for :math:`\alpha_v` = 1,
0.7, 0.5, 0.3, 0.1. FastWindTerrain's numbers come from its plotfiles (``Tests/fwt_bump_speedup.py``) and are stored
in ``Tests/Gold/init_bump_trend.json`` with the commit. The hilltop speed-up :math:`|u_h|/S_0(h)` (:math:`S_0` the
inflow profile at the height :math:`h` above the ground):

.. list-table::
   :header-rows: 1
   :widths: 12 22 22 22 22

   * - :math:`\alpha_v`
     - 20 m (ablsteady / FWT)
     - 40 m
     - 80 m
     - 150 m
   * - 1.0
     - 1.346 / 1.335
     - 1.253 / 1.247
     - 1.157 / 1.155
     - 1.084 / 1.083
   * - 0.7
     - 1.390 / 1.382
     - 1.260 / 1.255
     - 1.143 / 1.141
     - 1.068 / 1.067
   * - 0.5
     - 1.398 / 1.394
     - 1.240 / 1.236
     - 1.118 / 1.117
     - 1.052 / 1.052
   * - 0.3
     - 1.335 / 1.336
     - 1.176 / 1.174
     - 1.079 / 1.079
     - 1.037 / 1.037
   * - 0.1
     - 1.154 / 1.158
     - 1.088 / 1.092
     - 1.052 / 1.052
     - 1.032 / 1.032

The two codes agree to 0.011 at 20 m, 0.006 at 40 m, 0.002 at 80 m and 0.0008 at 150 m above the ground (0.3–0.8 %
of the speed-up), with the same ranking of the five :math:`\alpha_v` values at 40 m and above (monotone at 80 m and
above, the maximum at :math:`\alpha_v = 0.5`–0.7 near the ground, where a smaller :math:`\alpha_v` first steers the
flow around the hill and then, as it removes the vertical adjustment altogether, flattens the near-ground speed-up; at
20 m two of the values are within 0.002 of each other in both codes, so the ranking there is not compared). Every
speed-up is above 1. The agreement is as close as it is because both codes impose no flux through the ground on the
first guess and solve the same variational problem; the meshes differ (FastWindTerrain: Cartesian, first layer 2 m,
immersed boundary, four passes of an approximate nodal projection; ablsteady: terrain-following, first layer 2.5 m, one
exact projection). The test pins the ranking at 40 m and above, the sign, the per-height spread (0.02, 0.01, 0.005,
0.003) and ablsteady's own numbers to 1e-6.

The converged column as the first guess
---------------------------------------

``init.first_guess = column`` (sub-phase 2.1b, ``Source/ColumnGuess.H``) solves the deck's horizontally uniform
column as a 1-D problem first: on the cells of the deck's vertical grid (its nodal heights from the mesh builder
over flat ground), with :math:`u + iv`, :math:`k` at the cell centres and :math:`\theta` held at the sounding,

.. math::

   \frac{\partial (u + iv)}{\partial t} = \frac{\partial}{\partial z}\Bigl(\nu\,\frac{\partial (u + iv)}{\partial z}\Bigr)
   - i f\,\bigl((u + iv) - (u_g + i v_g)\bigr) + F_x + i F_y, \qquad
   \frac{\partial k}{\partial t} = \frac{\partial}{\partial z}\Bigl(\gamma_k\,\frac{\partial k}{\partial z}\Bigr)
   + P_s + P_b - \varepsilon ,

by backward Euler steps of growing size (10 s to :math:`10^5` s): the diffusion, the Coriolis term (implicit through
the complex unknown), the wall drag and the Rayleigh layer implicit (a complex and a real tridiagonal solve per
step), the :math:`k` sinks point-implicit (Patankar), the closure's point kernels (:doc:`closure`, the same algebra as
the 3-D closure with :math:`N^2` from the sounding and :math:`S^2` from the cell's two faces) lagged one step and
:math:`\nu_t` relaxed by a half against the previous step. The ground is the neutral wall law (:math:`u_* = \kappa
|U_1| / \ln((d_1 + z_0)/z_0)`, the stress :math:`u_*^2` along :math:`U_1` as an implicit drag, the first cell's
:math:`k = u_*^2/C_{\mu 0}^2`), the top the deck's (slip: zero flux; geostrophic: :math:`u_g, v_g`; no slip: 0),
the Rayleigh layer the deck's on the fields ``bc.rayleigh_fields`` selects (:math:`\nu_t` is relaxed by a half
between steps and the face diffusivities are arithmetic means whatever ``turb.nut_relax`` and ``turb.face_mean``
say: the column is a first guess). The asymptotic contraction is about 0.88 per step, set by the lagged
:math:`k`–:math:`\nu_t`–:math:`S^2` loop (the review's Jacobian of the step map), not by the step size. The column
converges when the largest
change of :math:`u, v` (scaled by the geostrophic or first-cell speed) and :math:`k` (by its maximum) over a step
is below ``init.column_tol`` (1e-8; lifted to 100 epsilon of ``Real``, 1.2e-5 in single precision, since the closure's
point kernels run in ``Real``), typically in 100–200 steps and a few hundredths of a second; it aborts, naming
``init.column_max_iterations`` and the last change, otherwise. The converged :math:`u, v, k` per level with a
ground row (:math:`u = v = 0`, the first cell's :math:`k`) become the first guess's profile in place of the inflow
profile, applied by height above the local ground exactly as the inflow profile is (``init.z_agl_min`` and
``init.obrien`` act as for the profile; ``init.k_decay_height`` is refused, it would decay the converged k again;
below the first cell centre the table is linear from the ground row, where the log law is logarithmic, a
difference confined to a crest's thinnest cells); :math:`\theta` stays the sounding at the absolute height. The
projection then runs as for ``profile`` unless ``init.projection = 0``. An inflow face keeps the deck's inflow
profile on the face while the interior starts from the column; the stage reconciles the two as it does for any
first guess. The application writes the table to ``column_guess.dat`` in the run directory.

**Why.** The inflow profile has u along the geostrophic wind and v = 0 at every height, and the mass-consistent
correction knows nothing of the Coriolis balance, so a Coriolis-driven case starts out of Ekman balance and the
stage spends the inertial adjustment in pseudo-time: the point-implicit rotation damps that mode by about
:math:`f\Delta\tau` per iteration, and on ERF's 2-D hill the cross-wind residual sits on a plateau for the first
1500 of 3062 iterations (:doc:`simplec`). The column has the balance (and the closure's equilibrium k) built in,
to round-off on the flat case. The two switches are independent: ``profile`` with the projection (the earlier
default), ``column`` alone (the balanced profile on the terrain-following mesh, the stage's first pressure correction
removing the divergence) and ``column`` with the projection (the balanced profile bent over the terrain by the
α-weighted correction). Over a hill the lifted spiral carries an advective imbalance on the flanks that neither the
projection nor the profile guess has, so the first iteration's residual is larger than the profile guess's; the
balance shows after the CFL ramp (:doc:`testing`, the spin-up check).

**Where it does not help.** On the flat Neutral_ABL_Flat deck the column seed is the 1-D steady state of u and v to
round-off, but its k, equilibrated with the sounding's sharp inversion base, differs from the 3-D closure's at the
two cells there (the first iteration's k residual is 8e-2), and the 3-D run then falls into a slow oscillation of
k at the inversion (the k residual swinging between 1e-5 and 5e-4 with a period of about 2500 iterations, the mean
speed between 11.03 and 11.16 m/s) that neither the rotation relaxation, the LU-SGS sweeps nor skipping the
projection removes, and the verdict is not reached in 20 000 iterations where the profile start takes 2839. The
flat deck therefore keeps ``profile``; the terrain cases, whose residual is set by the hill rather than by the
inversion base, converge in a third of the iterations (:doc:`simplec`). Seeding k at the inversion base is left
for a later sub-phase.

**Why θ stays the sounding.** A column with an inversion warms without end (the entrainment of :doc:`simplec`), so
its θ has no fixed point; the Rayleigh layer and the top condition reference the sounding; and the column's θ
applied by height above the ground would lift the inversion over a crest against them. The price is a
:math:`k`–θ mismatch in the two cells at the inversion base, where the column's k equilibrated with the sounding's
sharper base, gone within about a hundred iterations.

The column needs a horizontal forcing (``forcing.geostrophic = 1`` or a non-zero ``forcing.pressure_gradient``),
the stage (``steady.stage = simplec``, whose closure and forcing it uses), the neutral surface layer with the
Dirichlet wall k at the first cell centre (``bc.surface_zref`` is not supported), a slip, geostrophic or no-slip
top, and ``bc.rayleigh_fields`` holding u and v together or neither; a pressure gradient without Coriolis is allowed
only with a neutral sounding or a Rayleigh layer holding u and v, because above an inversion the stress
:math:`F(H - z)` falls on a layer whose :math:`k` sits at the floor and the column has no steady state (the review's
finding); ``check_column_guess`` refuses the option otherwise, naming ``init.first_guess``, and a column whose
change stops decreasing aborts after 200 stalled steps. An inflow face in the deck is no obstacle: the column is a
profile by height. A surface heat flux and the source wall k are not in this version.

Inputs
------

Every key is in ``init`` (:doc:`inputs` has the table). The keys other than ``first_guess`` have no effect without
``first_guess = profile`` and are refused then. With the initializer on, the application also reads ``bc.*`` and
``solver.poisson_*``. Of the ``bc`` keys the initial field uses the face types (which must fit the mesh's periodicity),
the inflow profile and the sounding; the top gradient, the Rayleigh and sponge damping and the surface layer's options
(``bc.zhi_theta_grad``, ``bc.rayleigh_*``, ``bc.sponge_*``, ``bc.surface_*``) are range-checked but have no effect on
the initial field: the SIMPLEC stage of sub-phase 1.14 applies them. The ``init``, ``bc`` and ``solver`` readers each stop at
their first error, after the mesh, terrain and io inputs have been reported together. A restart (``io.restart``) with
the initializer on is refused: the checkpoint holds the state.

References
----------

O'Brien, J. J. (1970). Alternative solutions to the classical vertical velocity problem. *J. Appl. Meteorol.* 9,
197–203.

Sherman, C. A. (1978). A mass-consistent model for wind fields over complex terrain. *J. Appl. Meteorol.* 17, 312–319.

Forthofer, J. M., Butler, B. W., Wagenbrenner, N. S. (2014). A comparison of three approaches for simulating fine-scale
surface winds in support of wildland fire management. Part I. Model formulation and comparison against measurements
(WindNinja). *Int. J. Wildland Fire* 23, 969–981.
