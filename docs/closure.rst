Turbulence closure
==================

``Source/Closure.{H,cpp}`` and ``Source/Closure_K.H`` are the one-equation k–L closure of Axell and Liungman (2001,
*Environ. Fluid Mech.* 1, 71–106; AL01), sub-phase 1.9: the eddy viscosity and diffusivities from the state, and the
sources of :math:`k`. The point relations are ERF's (``extern_erf/ERF_RANSClosure.H``, copied unchanged, the
reference of UT-1.9.1), composed as ERF's ``ComputeTurbulentViscosityRANS`` composes them with :math:`\rho = 1`.
The residual (:doc:`residual`) uses the closure when ``ResidualInputs::closure`` is set, and then computes the
diffusivities from the state inside ``residual()``, so ``jvp`` differentiates through them.

Eddy viscosity and length scale
-------------------------------

.. math::

   \nu_t = c_\mu\,k^{1/2}\,L, \qquad \nu_t' = c_\mu'\,k^{1/2}\,L = \frac{c_\mu'}{c_\mu}\,\nu_t .

- :math:`L` is the turbulent length: the geometric length :math:`l_g = \kappa(d + z_0)`, capped harmonically at
  ``turb.max_geom_lscale`` (:math:`l_{max} l_g/(l_{max} + l_g)`; AL01 Eq. 22), when :math:`|N^2| \le \epsilon`;

  .. math::

     \frac{1}{L^2} = \frac{1}{l_g^2} + \frac{N^2}{C_b^2 k} \quad (N^2 > 0,\ \text{Eq. 26}), \qquad
     L = l_g\left(1 - \frac{C_{\mu 0}^6}{C_b^2}\,R_t^s(l_g)\right)^{1/2} \quad (N^2 < 0,\ \text{Eq. 28}),

  with :math:`N^2 = (g/\theta_0)\,\partial\theta/\partial z` the Brunt–Väisälä frequency squared (the physical
  gradient at the cell centre, :doc:`transport`), :math:`d` the wall distance and :math:`z_0` the surface layer's
  roughness (0 without one, as ERF). The unstable length uses the Richardson number of the geometric length,
  smoothed (next item), and is evaluated once, as ERF and Kynema do: iterating Eq. 28 has no fixed point in
  strongly convective, weakly turbulent air.
- :math:`R_t = L^2 N^2/(C_{\mu 0}^6 k)` is the turbulent Richardson number, smoothed below ``turb.Rt_crit`` by
  Burchard and Petersen's map onto :math:`(R_t^{min}, R_t^{crit})`, written as ERF writes it,
  :math:`R_t^{crit} + a x/(x + a)` with :math:`x = R_t - R_t^{crit}`, :math:`a = R_t^{min} - R_t^{crit}`, which stays
  finite for :math:`|R_t|` beyond :math:`10^{15}` (UT-1.9.3).
- The stability functions (AL01 Eq. 31, 32):

  .. math::

     c_\mu = \frac{C_{\mu 0} + 0.108\,R_t}{1 + 0.308\,R_t + 0.00837\,R_t^2}, \qquad c_\mu' = \frac{C_{\mu 0}}{1 + 0.277\,R_t}.

  Their poles sit near :math:`R_t = -3.61`, so ``turb.Rt_min`` must be above −3.6.
- :math:`k` is floored at :math:`\max(\texttt{turb.tke\_floor}, \epsilon_{machine})` before everything (ERF's rule).

**Diffusivities.** :math:`\nu = \nu_{mol} + \nu_t`, :math:`\gamma_\theta = \alpha_{mol} + \nu_t'`, so the
turbulent Prandtl number is :math:`\mathrm{Pr}_t(R_t) = c_\mu/c_\mu'`, 1 when neutral, in every direction. ERF's
k-eqn uses :math:`\nu_t'` for the vertical heat diffusivity (``Khv``) but :math:`\nu_t/\mathrm{Pr}_t` with
``erf.Pr_t`` = 1/3 horizontally unless ``erf.rans_consistent_diffusivities`` is set; ablsteady is that consistent
choice, so the horizontal heat flux differs from an ERF default run. And :math:`\gamma_k = \nu_{mol} + \nu_t/\sigma_k` with :math:`\sigma_k = 1` (AL01 Table I, ERF's value for
the k-eqn). The ground ghost gradient of θ (:doc:`boundary`) therefore needs no Prandtl factor.

Sources of :math:`k`
--------------------

.. math::

   R_k = C(k) - D_{\gamma_k}(k) - S_k, \qquad S_k = P_s + P_b - \varepsilon, \qquad
   P_s = \nu_t S^2, \quad P_b = -\nu_t' N^2, \quad \varepsilon = C_{\mu 0}^3\,\frac{k^{3/2}}{L},

with :math:`S^2 = 2 S_{ij} S_{ij}` from the cell-centred physical velocity gradients and :math:`N^2` with the
constant :math:`\theta_0` of the buoyancy options (ERF's ``theta_ref`` branch; the local-θ variant is not provided).
ERF forms :math:`P_b` from the mean of the two face heat fluxes, which equals this centred form on uniform layers and
differs by :math:`(r - 1)^2` on stretched ones (0.8 % at a ratio 1.05 in the first cell).

- **What makes the log law hold:** :math:`C_{\mu 0}^3` in :math:`\varepsilon` (AL01 Sect. 4.1): in the equilibrium
  state :math:`u = (u_*/\kappa)\ln((z + z_0)/z_0)`, :math:`k = u_*^2/C_{\mu 0}^2` with the uncapped :math:`l_g`,
  :math:`\nu_t = \kappa u_*(z + z_0)` and :math:`P_s = \varepsilon = u_*^3/(\kappa(z + z_0))` pointwise. RT-1.9.1 pins
  this on ERF's Neutral_ABL_Flat mesh: the only residual left is the discrete strain's truncation (the centred
  difference of :math:`\ln` over :math:`2\Delta z`): :math:`P_s/\varepsilon = (r/2 \ln((r+1)/(r-1)))^2` with
  :math:`r = (z + z_0)/\Delta z`, 1.45, 1.12, 1.06, 1.03 in cells 2–5 and :math:`1 + \tfrac23 (\Delta z/z)^2`
  aloft; the first cell is exact through the ground ghost's gradient. ERF (edge stress averaged to the centre) and
  Kynema (centred strain) share this stencil. What the steady solver inherits on the 15.6 m mesh, from a 1-D column
  solve of the coupled equations: :math:`k` 14 % above :math:`u_*^2/C_{\mu 0}^2` in the second cell (8 % in the
  third), :math:`\nu_t` 7 % high there and :math:`u` 1 % low aloft. A strain from the product of the two face
  gradients would overshoot by :math:`(\Delta z/z)^2/12` only; not in this PR.
- **Point-implicit treatment.** The line solves (``line_jacobian``) add the Patankar diagonal
  :math:`(\varepsilon + \max(-P_b, 0))/k` to :math:`k`'s rows: the sinks linearised as (sink/k) k, so the point
  update :math:`k + S_k/(1/\Delta\tau + \text{diag})` stays non-negative for any pseudo-time step (UT-1.9.4). ERF's
  ``implicit_tke_dissipation`` is the same idea for the dissipation alone. Against the exact derivative
  :math:`-\partial S_k/\partial k = (1.5\,\varepsilon - 0.5\,P_s)/k` (neutral, :math:`L` fixed) the Patankar value is
  exact at equilibrium, an under-relaxation where :math:`P_s < \varepsilon` or the air is strongly stable (0.3 of the
  exact value at :math:`R_t = 3`), and the only positive choice where :math:`P_s > 3\varepsilon` or near the floor
  (the exact derivative is dominated by :math:`0.5 P_s/k`, which diverges there). ``ClosureFields::dSdk_exact``
  holds the exact derivative for a Newton diagonal, :math:`\max(\text{Patankar}, -\partial S_k/\partial k)`.
- **Analytic Jacobian.** ``Closure_K.H`` is templated on the scalar type; evaluated on dual numbers (``Dual.H``) it
  gives :math:`dS_k/dk` and :math:`d\nu_t/dk` of exactly the arithmetic it evaluates, with the branch decisions
  held (UT-1.9.5 against finite differences, relative 1e-6). The *frozen-switch* mode of ``Closure`` records every
  point's decisions (neutral/stable/unstable, the two smoothings, the floor) at a state and keeps them while a
  Jacobian is taken; a frozen evaluation is meant for perturbations of that state (past a switch the smoothing is
  clamped at :math:`R_t^{crit}`), not for a line search.
- **Relaxation.** With ``turb.nut_relax`` :math:`w < 1`, :math:`\nu_t = w\,\nu_t^{new} + (1 - w)\,\nu_t^{stored}`,
  the stored field being the last one ``Closure::accept`` took; the residual stays a function of the state and of
  that field, which goes to the checkpoint (``Level_0/closure_nut``, ghost ring included) and is restored exactly
  on any box layout. The frozen switches are transient (one Jacobian) and are not checkpointed.

Trap safety
-----------

The kernels take their branch decisions once (``switches_for``) and then evaluate every arm with an argument that
is valid for every input: the stable arm with :math:`\max(N^2, 0)`, the unstable arm with :math:`\min(N^2, 0)`,
the smoothing with :math:`\min(x, 0)`, all written as :math:`\tfrac12(x \pm |x|)` rather than as a select. The
compiler evaluates both arms of a select and drops a clamp inside an arm where it knows the sign, so the square
root of :math:`1 - 16\zeta` or :math:`1 - C_{\mu 0}^6 R_t/C_b^2` for the other branch would otherwise raise the
invalid-operation trap on a value that is never used (``amrex.fpe_trap_invalid``). Where the arm applies these
forms equal the raw value exactly (doubling and halving are exact), so the results stay bit-identical to ERF's
(UT-1.9.1). The same holds for the logarithmic face mean below.

Face average of the diffusivities
---------------------------------

The diffusion operators average the cell diffusivity to the faces (:doc:`transport`). ``turb.face_mean`` chooses
``arithmetic`` (the default, bit-identical to before) or ``logarithmic``: the logarithmic mean
:math:`(a - b)/\ln(a/b)` across each vertical cell pair, arithmetic across lateral pairs. For a diffusivity linear
in the vertical coordinate, the wall law's :math:`\kappa u_*(z + z_0)`, the logarithmic mean gives the face flux of
the logarithmic profile exactly: in RT-1.9.1 every interior momentum flux is :math:`u_*^2` to 1e-12, where the
arithmetic mean misses the first face by 10 % (RT-1.8.1's 9.7 % truncation). The implicit operator uses the same
mean, so :math:`A_1 u = -D(u)` holds for either choice. The mean is written without a select feeding a division or
a logarithm (series for :math:`|a - b| < 10^{-4}(a + b)`, the direct form otherwise, each arm's divisor non-zero).

Wall distance
-------------

``height_above_ground()`` gives the vertical height of each cell centre above its column's ground (the node-mean
height minus the node-mean ground; ERF's ``terrain_height``), with one ghost cell. ``turb.wall_distance = tucker``
names the Tucker distance of :doc:`walldistance` instead; ``closure_column`` honours it (the Tucker solve with the
Poisson defaults), as the SIMPLEC stage (sub-phase 1.14) will. Likewise ``turb.frozen_switch = 1`` makes the tool record
the branch decisions at the state before the residual. The
surface layer's normal distance :math:`d_0` and this vertical height agree on flat ground and differ by the slope
factor :math:`|A|` over terrain.

Inputs
------

Every key is in ``turb``. In this version only the regression tool ``closure_column`` reads them; the driver of
sub-phase 1.14 adds them to the application, where setting them still stops the run. Each row matches the check in
``check_closure_options()``, which names the key.

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

The closure's :math:`C_{\mu 0}` and :math:`\kappa` are also the surface layer's (``SurfaceLayerOptions``) and the
inflow profile's; the driver will set them together.
