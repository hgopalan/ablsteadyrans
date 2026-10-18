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
