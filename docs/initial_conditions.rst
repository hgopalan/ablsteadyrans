Initial conditions
==================

``init.first_guess = profile`` builds a mass-consistent, WindNinja-style starting field: a first guess from the
inflow profile, made divergence-free by a variational correction. The method follows
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
