Governing equations in curvilinear form
=======================================

Over terrain the equations of :doc:`governing_cartesian` are solved on ERF's terrain-following mesh. A uniform
computational grid :math:`(\xi, \eta, \zeta)` maps to physical space :math:`(x, y, z(\xi, \eta, \zeta))`. The
horizontal coordinates are not distorted; only the heights follow the ground. The code that builds the mapping is
ERF's, copied into ``extern_erf/`` (:doc:`extraction`).

Vertical levels
---------------

Over flat ground the node heights are the reference levels :math:`0 = z_0 < z_1 < \dots < z_{n_z}`:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Inputs
     - Levels
   * - none (default)
     - Uniform, :math:`z_k = k\,\Delta\zeta`, :math:`\Delta\zeta = z_\mathrm{top}/n_z`.
   * - ``mesh.dz0``
     - Geometric, :math:`z_k - z_{k-1} = \Delta z_0\,r^{k-1}`, with :math:`r` solved so that the top is
       ``mesh.prob_hi`` z.
   * - ``mesh.dz0``, ``mesh.stretch_ratio``
     - Geometric with the given :math:`r`; the top is where the layers end.
   * - ``mesh.zlevels_file``
     - The :math:`n_z + 1` heights in the file.

Terrain-following mapping
-------------------------

With the ground height :math:`h(x, y)`, ``terrain.smoothing`` selects one of ERF's three mappings:

.. math::

   \text{0, BTF:}\quad & z = \frac{h\,z_\mathrm{top} + (z_\mathrm{top} - h)\,z_k}{z_\mathrm{top}},\\
   \text{1, STF (default):}\quad & z = z_k + A(z_k)\,h_s(x, y, k), \qquad
     A(z) = \cos^6\!\Bigl(\frac{\pi}{2}\frac{z}{z_H}\Bigr),\quad z_H = \frac{2.44\,h_\mathrm{max}}{1 - \gamma_m},\\
   \text{2, Sullivan:}\quad & z = z_k + \Bigl(1 - \frac{z_k}{z_\mathrm{top}}\Bigr)^3 h .

- :math:`z_\mathrm{top}` is the physical top, the last level.
- STF (Klemp 2011) smooths the terrain imprint :math:`h_s` level by level and is 0 above :math:`z_H`. A smoothing pass
  that would make a layer thinner than :math:`\gamma_m` (``terrain.gamma_m``, 0.5) of its flat thickness is
  rejected, so STF never inverts a cell.
- The run stops if the ground reaches the top or any layer has zero or negative thickness.

``terrain.type`` sets :math:`h`: ``flat``, ``agnesi`` (Witch of Agnesi ridge or round hill), ``xyz`` (a point cloud,
inverse-distance weighted as in FastWindTerrain) or ``erf_grid`` (ERF's terrain file, interpolated as ERF does).

Metrics
-------

The metric terms are differences of the node heights:

.. math::

   z_x = \frac{\partial z}{\partial \xi}, \qquad z_y = \frac{\partial z}{\partial \eta}, \qquad
   J = z_\zeta = \frac{\partial z}{\partial \zeta} .

- :math:`z_x, z_y` are the slopes of a coordinate surface.
- :math:`J` is the cell volume divided by :math:`\Delta x\,\Delta y\,\Delta\zeta`. The x- and y-face area factors
  :math:`a_x, a_y` are :math:`z_\zeta` on those faces.

Physical derivatives follow from the chain rule:

.. math::

   \frac{\partial q}{\partial x} = q_\xi - \frac{z_x}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial y} = q_\eta - \frac{z_y}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial z} = \frac{q_\zeta}{z_\zeta} .

Conservation form
-----------------

The divergence of any flux :math:`\mathbf G` becomes

.. math::
   :label: curvilinear-divergence

   \nabla\cdot\mathbf G = \frac{1}{J}\left[\frac{\partial(\mathbf G\cdot\mathbf A_\xi)}{\partial\xi}
   + \frac{\partial(\mathbf G\cdot\mathbf A_\eta)}{\partial\eta}
   + \frac{\partial(\mathbf G\cdot\mathbf A_\zeta)}{\partial\zeta}\right],
   \qquad
   \mathbf A_\xi = (z_\zeta, 0, 0),\quad \mathbf A_\eta = (0, z_\zeta, 0),\quad \mathbf A_\zeta = (-z_x, -z_y, 1).

:math:`\mathbf A` are the face area vectors per unit index area. x- and y-faces are vertical; z-faces follow the
terrain. Applied to the velocity, the flux through a z-face is the contravariant velocity

.. math::

   \Omega = w - z_x\,u - z_y\,v ,

and continuity reads

.. math::

   \frac{1}{J}\left[\frac{\partial(a_x u)}{\partial\xi} + \frac{\partial(a_y v)}{\partial\eta}
   + \frac{\partial\Omega}{\partial\zeta}\right] = 0 .

The momentum, :math:`\theta` and :math:`k` equations take the same form, with the convective flux
:math:`\mathbf u\,q` and the diffusive flux :math:`\gamma\nabla q` (or :math:`\boldsymbol\tau`) in :math:`\mathbf G`.
The state keeps the Cartesian velocity components; only the fluxes are contravariant.

**Why this form.** The metrics satisfy the geometric conservation law: a uniform wind has zero discrete divergence on
any terrain (tested to 1e-12 for all three mappings). Over flat ground :math:`z_x = z_y = 0`, :math:`\Omega = w` and
the equations are the Cartesian ones.

Mesh quality
------------

Every run prints the smallest cell and layer thickness relative to flat, the largest ground and coordinate-surface
slope, the largest aspect ratio, the lowest first-cell centre, the largest top-node displacement and the number of
STF levels capped by :math:`\gamma_m`. A restart checks that the node heights rebuilt from the inputs match the
checkpoint.

References
----------

- Klemp, J. B. (2011). A terrain-following coordinate with smoothed coordinate surfaces. *Mon. Weather Rev.* 139,
  2163–2169.
