Initial conditions
==================

The steady solver starts from an initial field. A good initial field is already close to the solution, so the solver
needs fewer iterations. GUST builds the initial field in two steps:

#. **Initial profile.** A wind profile is placed on the mesh by height above the local ground.
#. **Mass-consistent adjustment.** The field is corrected so that it conserves mass, which bends the flow over and
   around the terrain.

The method is that of mass-consistent wind models such as WindNinja (Forthofer et al. 2014) and follows
FastWindTerrain. No FastWindTerrain code is copied (it uses a Cartesian immersed-boundary mesh), but every parameter
has the same meaning. The code is in ``Source/Initializer`` and ``Source/ColumnGuess``; the inputs are in the
``init`` namespace (:doc:`inputs`).

Options
-------

``init.first_guess`` selects the initial profile:

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Value
     - Initial profile
     - Recommended for
   * - ``none`` (default)
     - Fluid at rest.
     - Tests only.
   * - ``profile``
     - The inflow profile of :doc:`boundary_conditions` (log law, power law or file).
     - Cases without Coriolis forcing; flat cases with an inversion.
   * - ``column``
     - The steady solution of a single vertical column, computed first (below).
     - Terrain cases with Coriolis forcing.

``init.projection`` (default 1) switches the mass-consistent adjustment on or off. With 0, the solver's own pressure
correction removes the divergence in its first iteration.

Initial profile from the inflow profile
---------------------------------------

With ``init.first_guess = profile``:

- :math:`u_0, v_0` are the inflow profile at each face's height above the ground,
  :math:`h = \max(z_\mathrm{face} - z_\mathrm{ground}, h_\mathrm{min})`, where :math:`h_\mathrm{min}` is
  ``init.z_agl_min``.
- :math:`w_0 = 0`, except on the ground faces, where :math:`w_0 = u_0 z_x + v_0 z_y` (no flow through the ground).
- :math:`\theta` is the sounding ``bc.theta_profile`` at the absolute height, and :math:`p = 0`.
- :math:`k` is the profile's equilibrium value. Above :math:`z_i` = ``init.k_decay_height`` (0: no decay) it decays:

  .. math::

     k = k_\mathrm{prof}(h)\, e^{-(h - z_i)/z_i} \quad (h > z_i), \qquad k = k_\mathrm{prof}(h) \quad (h \le z_i).

Initial profile from a single column
------------------------------------

With Coriolis forcing, the steady wind turns with height (the Ekman spiral). The inflow profile has no such turning,
so the solver must first spend many iterations on the slow inertial adjustment. ``init.first_guess = column`` avoids
this: it first solves the horizontally uniform problem of the case, a single vertical column, and uses its steady
profile as the initial profile.

The column has the vertical grid of the case over flat ground. It solves for the horizontal wind, written as the
complex number :math:`u + iv`, and for :math:`k`:

.. math::

   \frac{\partial (u + iv)}{\partial t} = \frac{\partial}{\partial z}\Bigl(\nu\,\frac{\partial (u + iv)}{\partial z}\Bigr)
   - i f\,\bigl((u + iv) - (u_g + i v_g)\bigr) + F_x + i F_y, \qquad
   \frac{\partial k}{\partial t} = \frac{\partial}{\partial z}\Bigl(\gamma_k\,\frac{\partial k}{\partial z}\Bigr)
   + P_s + P_b - \varepsilon .

- :math:`\theta` is held at the sounding.
- :math:`\nu`, :math:`\gamma_k`, :math:`P_s`, :math:`P_b`, :math:`\varepsilon` come from the turbulence model of
  :doc:`turbulence_model`.
- The ground uses the neutral log law; the top and the Rayleigh layer are those of the case.

The column is marched to steady state with implicit (backward Euler) steps that grow from 10 s to :math:`10^5` s. It
has converged when the largest change of :math:`u, v` and :math:`k` over a step, relative to their scale, is below
``init.column_tol``. This typically takes 100–200 steps and a few hundredths of a second. If it has not converged
after ``init.column_max_iterations`` steps, the run stops.

The converged :math:`u, v, k` are then placed on the mesh by height above the local ground, exactly as the inflow
profile is; :math:`\theta` stays the sounding. The table is written to ``column_guess.dat`` in the run directory.
An inflow face keeps the case's inflow profile; the solver reconciles the two.

**Requirements.** The column option needs:

- a horizontal forcing: ``forcing.geostrophic = 1`` or a non-zero ``forcing.pressure_gradient``;
- the SIMPLEC stage (``steady.stage = simplec``);
- a ``surface_layer`` ground with the default wall :math:`k` (``bc.surface_zref`` is not supported), neutral
  unless ``init.column_stability = heat_flux`` (then ``bc.surface_stability = heat_flux``) or ``mol`` with a
  non-zero ``init.column_mol`` (``bc.surface_stability = heat_flux`` allowed: the run takes the column's flux);
- with a non-neutral column, a Rayleigh layer on :math:`\theta` (below);
- a ``slip_wall``, ``geostrophic`` or ``no_slip_wall`` top;
- ``bc.rayleigh_fields`` either containing both ``u`` and ``v`` or neither;
- with a pressure gradient but no Coriolis force: a neutral sounding or a Rayleigh layer on ``u`` and ``v``.

``init.k_decay_height`` cannot be used with it. Otherwise the run stops with a message naming
``init.first_guess``.

.. _column-stability:

**Stratified column.** With a non-neutral surface (``init.column_stability``, the forms of
amrTerrain's 1-D column solver) the column also solves for :math:`\theta`,

.. math::

   \frac{\partial \theta}{\partial t} = \frac{\partial}{\partial z}\Bigl(\gamma_\theta\,\frac{\partial \theta}{\partial z}\Bigr)
   - c_r(z)\,(\theta - \theta_p),

with the deck's gradient ``bc.zhi_theta_grad`` at the top and, at the ground, Monin–Obukhov similarity for a
surface temperature :math:`\theta_s` (ERF's Businger–Dyer :math:`\psi_m`, :math:`\psi_h`, or QNSE's on the stable
side with ``bc.surface_stable_functions = qnse``, :ref:`stable functions <stable-functions>`):

.. math::

   u_* = \frac{\kappa U_1}{\ln\frac{z_1 + z_0}{z_0} - \psi_m(\zeta)}, \qquad
   q = -\frac{\kappa u_*\,(\theta_1 - \theta_s)}{\ln\frac{z_1 + z_0}{z_0} - \psi_h(\zeta)}, \qquad
   \zeta = \frac{z_1}{L}, \quad \frac{1}{L} = -\frac{\kappa g q}{u_*^3 \theta_0},

- :math:`U_1`, :math:`\theta_1`: the first cell's speed and potential temperature at height :math:`z_1`;
- :math:`\zeta` from the bulk Richardson number, as ERF's surface-temperature solve:
  :math:`\zeta = \mathrm{Ri}_b\,(\ln - \psi_m)^2/(\ln - \psi_h)`,
  :math:`\mathrm{Ri}_b = z_1 g (\theta_1 - \theta_s)/(U_1^2 \theta_0)` held in
  [``init.column_rib_min``, ``init.column_rib_max``] (default ±4, ERF's), solved by ERF's fixed point
  under-relaxed by a half. Businger–Dyer's linear stable :math:`\psi` has a root only for
  :math:`\mathrm{Ri}_b < 0.2` (:math:`\zeta = \mathrm{Ri}_b \ln/(1 - 5\,\mathrm{Ri}_b)`), so with them the bound
  is at most 0.19 (:math:`\zeta = 3.8\ln`); QNSE's :math:`\phi` levels off and needs no further cap. The flux takes
  the full :math:`\theta_1 - \theta_s`, as in ERF and WRF;
- :math:`q`: the kinematic surface heat flux (K m/s, upward > 0), implicit in :math:`\theta_1`;
- the wall :math:`k` the surface layer's (Axell and Liungman 2001, Eq. 16):
  :math:`k_1 = (u_*^3 + \kappa B z_1)^{2/3}/C_{\mu 0}^2`, :math:`B = g q/\theta_0` for an upward flux only.

``surface_temperature`` fixes :math:`\theta_s` (``init.column_surface_theta``). ``mol`` and ``heat_flux`` fix
:math:`1/L` (``init.column_mol``) or :math:`q` (``bc.surface_heat_flux``), and the column finds the :math:`\theta_s`
that gives it by the secant method. The search is on :math:`1/L` for ``mol``, so it reaches every :math:`L` the
column allows, also those past its largest downward flux. :math:`\zeta` is taken at the first cell centre,
:math:`z_1/L`, as in the 3-D wall model.

**Heat budget.** At steady state the surface flux has to leave (or enter) the column somewhere. The top passes
almost nothing: its flux is the deck's gradient times a diffusivity at the floor. The Rayleigh layer on
:math:`\theta` is therefore the column's heat source or sink, and a non-neutral column requires it
(``bc.rayleigh_zdamp`` > 0 and ``theta`` in ``bc.rayleigh_fields``). The most stable column the deck allows (the
smallest :math:`L`, the largest downward flux) depends on that layer's height and strength and on the wind. It is a
property of the deck, not of the turbulence model. A request beyond it stops the run, naming the smallest
:math:`L` (or the largest flux) reached. On the unit test's 64 m column (Rayleigh layer above 48 m) the smallest
:math:`L` is 414 m; on a 120 m deck with the Bolund vertical grid the largest downward flux is −0.029 K m/s with the
Rayleigh layer above 80 m and −0.022 K m/s above 100 m.

The log prints the flux, :math:`L` and :math:`\theta_s` reached, and the column's total steps and columns. With
``mol`` and ``bc.surface_stability = heat_flux`` the 3-D run takes the column's flux as ``bc.surface_heat_flux``.
The wall :math:`k` is the surface layer's, :math:`(u_*^3 + \kappa B z_1)^{2/3}/C_{\mu 0}^2`; amrTerrain's column
uses the additive form :math:`(u_*^2 + (\kappa B z_1)^{2/3})/C_{\mu}^2`, up to :math:`2^{1/3}` larger for an
upward flux.

**Target wind.** ``init.column_target_u``, ``init.column_target_v`` at ``init.column_target_height`` (amrTerrain's
met-mast forcing): the geostrophic wind needed to keep a given wind at a given height under each stratification.
After each converged column the geostrophic wind and the column's Rayleigh reference wind are multiplied by
:math:`R = U_\mathrm{target}/U(z_t)` (complex: scaled and turned), the pressure gradient by :math:`R\,|R|` (the
stress goes with the square of the speed), and the column is solved again, until :math:`|R - 1| < 10^{-4}`.
:math:`U(z_t)` is MOST's profile below the first cell centre and linear in :math:`\ln(z + z_0)` above. The 3-D run
then takes that forcing. It reaches the target only if nothing holds its wind to the deck's inflow profile, so a
target is refused with a Rayleigh layer on :math:`u, v`, a sponge or an inflow face (for an inflow case, run the
column on a periodic deck and give its ``column_guess.dat`` as ``bc.inflow_file``), and with a stratified column
over a neutral 3-D surface.

**Restart.** The checkpoint holds the state, not the adjusted forcing or the column's flux, and ``init.*`` keys are
refused on a restart. The log lists the inputs a restart deck needs in their place (``forcing.ug``,
``forcing.vg``, ``forcing.pressure_gradient``, ``bc.surface_heat_flux``), with all their digits.

On the Bolund case 4 deck with :math:`L` = 250 m and 7.58 m/s at 5 m, the column gives :math:`u_*` = 0.509 m/s
(0.5085 m/s measured). That agreement mainly reflects MOST at 5 m with the chosen :math:`z_0` and target. The flux,
:math:`q` = −0.039 K m/s, is sustained by the Rayleigh layer and has no measured counterpart.

The column's :math:`\theta` goes into ``column_guess.dat`` (a ``bc.inflow_file`` for a run with the same
stratification). The 3-D first guess keeps :math:`\theta` from the sounding.

**When to use it.** On the 2-D and 3-D hills the column start halves the iterations (:doc:`results`). On the flat
``Neutral_ABL_Flat`` case it does not converge: the column's :math:`k` at the sharp inversion base differs from the
3-D solution, and a slow oscillation of :math:`k` develops there. Use ``profile`` for that case.

Mass-consistent adjustment
--------------------------

The initial profile :math:`\mathbf u_0` does not conserve mass over terrain. The adjustment finds the field
:math:`\mathbf u` closest to :math:`\mathbf u_0` that does (Sherman 1978):

.. math::

   \nabla\cdot(\Sigma\nabla\lambda) = \nabla\cdot\mathbf u_0, \qquad
   \mathbf u = \mathbf u_0 - \Sigma\nabla\lambda, \qquad
   \Sigma = \mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2) .

- :math:`\lambda`: a Lagrange multiplier (a potential), with :math:`\lambda = 0` on inflow and outflow faces and no
  flux through walls, the ground and the top.
- :math:`\alpha_h, \alpha_v`: horizontal and vertical transmissivity (``init.alpha_h``, ``init.alpha_v``).

The equation is solved by the Poisson solver of :doc:`solution_algorithm`, terrain cross terms included.

**Choosing α.** Only the ratio :math:`\alpha_v/\alpha_h` matters. A smaller ratio allows less vertical motion, so
the flow goes around steep terrain rather than over it, as in stable air. In Sherman's notation
:math:`\alpha_v/\alpha_h` is :math:`\alpha_1/\alpha_2`, WindNinja's stability parameter. Ratios of 0.1 and above are
practical; smaller ratios need more Poisson iterations.

With :math:`\alpha_v \ne \alpha_h` this one solve always uses the ``gmres_fft`` back end, whatever
``solver.poisson_backend`` says, and the log says so. The SIMPLEC stage keeps ``solver.poisson_backend``.

O'Brien correction
------------------

``init.obrien = 1`` (off by default) first sets the vertical flux :math:`\Omega_0` from column continuity, integrating
from the ground up, and removes the remaining flux at the top with a weight that grows with height (O'Brien 1970):

.. math::

   a_z\Omega_{k+1/2} = a_z\Omega_{k-1/2} - J_k D^h_k\,\Delta\zeta, \qquad
   f_k = \frac{z_k - z_\mathrm{ground}}{z_\mathrm{top} - z_\mathrm{ground}} ,

where :math:`D^h` is the horizontal part of the divergence and the top residual is removed with the weight
:math:`f_k^2`. It is meant for initial profiles taken from soundings or other models.

The resulting field
-------------------

The initial field is divergence-free to the Poisson tolerance, and :math:`\theta`, :math:`k`, :math:`p` are
consistent with the sounding. It is a starting point, not a RANS solution: it does not satisfy the momentum
equations.

A restart (``io.restart``) cannot be combined with ``init.first_guess``: the checkpoint already holds the state.
