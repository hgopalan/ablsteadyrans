Initial conditions
==================

``init.first_guess = profile`` builds a mass-consistent, WindNinja-style starting field: a first guess from the
inflow profile, made divergence-free by a variational correction (``init.projection = 0`` keeps the first guess as it
is); ``init.first_guess = column`` takes the first guess from a converged 1-D column instead (below). The method follows
`FastWindTerrain <https://github.com/hgopalan/FastWindTerrain>`__. No code is copied (FastWindTerrain uses a
Cartesian immersed-boundary mesh), but every parameter has the same meaning. The code is ``Source/Initializer``.

First guess
-----------

- :math:`u_0, v_0` are the inflow profile (:doc:`boundary_conditions`) at each face's height above the ground,
  :math:`h = \max(z_\mathrm{face} - z_\mathrm{ground}, h_\mathrm{min})`, with :math:`h_\mathrm{min}` =
  ``init.z_agl_min``.
- :math:`w_0` is 0, except on the ground faces, where it takes the no-flux value :math:`u_0 z_x + v_0 z_y`.
- :math:`\theta` is the sounding ``bc.theta_profile`` at the absolute height, and :math:`p = 0`.
- :math:`k` is the profile's equilibrium value, decaying above :math:`z_i` = ``init.k_decay_height``:

  .. math::

     k = k_\mathrm{prof}(h)\, e^{-(h - z_i)/z_i} \quad (h > z_i), \qquad k = k_\mathrm{prof}(h) \quad (h \le z_i).

Mass-consistent correction
--------------------------

.. math::

   \nabla\cdot(\Sigma\nabla\lambda) = \nabla\cdot\mathbf u_0, \qquad
   \mathbf u = \mathbf u_0 - \Sigma\nabla\lambda, \qquad
   \Sigma = \mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2) .

- :math:`\alpha_h, \alpha_v` are ``init.alpha_h`` and ``init.alpha_v``.
- :math:`\lambda = 0` on inflow and outflow faces; no flux through walls, the ground and the top.
- The equation is solved by the Poisson solver of :doc:`solver`, terrain cross terms included.

**What α does.** :math:`\alpha` is a transmissivity. A smaller :math:`\alpha_v/\alpha_h` means less vertical
adjustment, so the flow goes around steep terrain rather than over it. Only the ratio matters. In Sherman's (1978)
notation :math:`\alpha_v/\alpha_h` is :math:`\alpha_1/\alpha_2`, WindNinja's stability parameter. Ratios of 0.1 and
above are practical; smaller ratios need more Poisson iterations.

O'Brien column correction
-------------------------

``init.obrien = 1`` (off by default, as in FastWindTerrain) first sets :math:`\Omega_0` from the column continuity
of :math:`(u_0, v_0)` from the ground up, and removes the column's top residual :math:`E` with the weight
:math:`f_k^2`:

.. math::

   a_z\Omega_{k+1/2} = a_z\Omega_{k-1/2} - J_k D^h_k\,\Delta\zeta, \qquad
   f_k = \frac{z_k - z_\mathrm{ground}}{z_\mathrm{top} - z_\mathrm{ground}} ,

where :math:`D^h` is the horizontal part of the divergence (O'Brien 1970). It is meant for first guesses from
soundings or fields.

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
step), the :math:`k` sinks point-implicit (Patankar), the closure's point kernels (:doc:`governing_cartesian`, the same algebra as
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
1500 of 3062 iterations (:doc:`results`). The column has the balance (and the closure's equilibrium k) built in,
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
inversion base, converge in a third of the iterations (:doc:`results`). Seeding k at the inversion base is left
for a later sub-phase.

**Why θ stays the sounding.** A column with an inversion warms without end (the entrainment of :doc:`results`), so
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

What the field is
-----------------

The field is divergence-free to the Poisson tolerance and has :math:`\theta, k, p` consistent with the sounding. It
is a starting point for the steady solver, not a RANS solution: it has no momentum balance. A restart with the
initializer on is refused, since the checkpoint holds the state.

References
----------

- O'Brien, J. J. (1970). *J. Appl. Meteorol.* 9, 197–203.
- Sherman, C. A. (1978). *J. Appl. Meteorol.* 17, 312–319.
- Forthofer, J. M., Butler, B. W., Wagenbrenner, N. S. (2014). *Int. J. Wildland Fire* 23, 969–981.
