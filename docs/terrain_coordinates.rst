Terrain-following coordinates
=============================

Over hills the ground is not flat, so a regular box-shaped grid does not fit it. GUST uses a
terrain-following mesh instead: the grid lines run parallel to the ground near the surface and become flat
higher up. The grid is still regular in computational space, which keeps the solver simple and fast.

A uniform computational grid :math:`(\xi, \eta, \zeta)` is mapped to physical space :math:`(x, y, z(\xi, \eta,
\zeta))`. Only the heights follow the ground; the horizontal coordinates are not distorted. The mapping code is
ERF's, copied into ``extern_erf/`` (:doc:`extraction`).

Vertical levels
---------------

Over flat ground the node heights are the reference levels :math:`0 = z_{(0)} < z_{(1)} < \dots < z_{(n_z)}`. In
the formulas below :math:`z_k` stands for :math:`z_{(k)}`; it is not the roughness length :math:`z_0`. They are set by
the ``mesh`` inputs:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Inputs
     - Levels
   * - none (default)
     - Uniform: :math:`z_k = k\,\Delta\zeta`, with :math:`\Delta\zeta = z_\mathrm{top}/n_z`.
   * - ``mesh.dz0``
     - Geometric: :math:`z_k - z_{k-1} = \Delta z_0\,r^{k-1}`, with the ratio :math:`r` chosen so that the top is at
       ``mesh.prob_hi`` z.
   * - ``mesh.dz0``, ``mesh.stretch_ratio``
     - Geometric with the given :math:`r`; the top is where the layers end.
   * - ``mesh.zlevels_file``
     - The :math:`n_z + 1` heights listed in the file.

Stretched levels put thin cells near the ground, where the wind changes fastest, and thick cells aloft.

Terrain-following mapping
-------------------------

Let :math:`h(x, y)` be the ground height and :math:`z_\mathrm{top}` the top of the domain. ``terrain.smoothing``
selects one of ERF's three mappings:

.. math::

   \text{0, BTF:}\quad & z = \frac{h\,z_\mathrm{top} + (z_\mathrm{top} - h)\,z_k}{z_\mathrm{top}},\\
   \text{1, STF (default):}\quad & z = z_k + A(z_k)\,h_s(x, y, k), \qquad
     A(z) = \cos^6\!\Bigl(\frac{\pi}{2}\frac{z}{z_H}\Bigr),\quad z_H = \frac{2.44\,h_\mathrm{max}}{1 - \gamma_m},\\
   \text{2, Sullivan:}\quad & z = z_k + \Bigl(1 - \frac{z_k}{z_\mathrm{top}}\Bigr)^3 h .

- **BTF** (basic terrain following) shrinks every column linearly; the terrain is felt up to the top.
- **STF** (smoothed terrain following, Klemp 2011) smooths the terrain imprint :math:`h_s` level by level, so the
  coordinate surfaces become flat above the height :math:`z_H`. A smoothing pass that would make a layer thinner than
  the fraction :math:`\gamma_m` (``terrain.gamma_m``, 0.5) of its flat thickness is rejected, so STF never folds a
  cell.
- **Sullivan** decays the terrain imprint with the cube of the distance to the top.

The run stops if the ground reaches the top or a layer has zero or negative thickness.

``terrain.type`` sets :math:`h`:

- ``flat``;
- ``agnesi``: a Witch of Agnesi ridge or round hill;
- ``xyz``: a point cloud, interpolated by inverse-distance weighting as in FastWindTerrain;
- ``erf_grid``: ERF's terrain file, interpolated as ERF does.

Metric terms
------------

The mapping is described by the slopes and the stretching of the coordinate surfaces, computed from the node
heights:

.. math::

   z_x = \frac{\partial z}{\partial \xi}, \qquad z_y = \frac{\partial z}{\partial \eta}, \qquad
   J = z_\zeta = \frac{\partial z}{\partial \zeta} .

- :math:`z_x, z_y`: slopes of a coordinate surface.
- :math:`J`: the Jacobian, the cell volume divided by :math:`\Delta x\,\Delta y\,\Delta\zeta`. The area factors
  :math:`a_x, a_y` of the x- and y-faces are :math:`z_\zeta` on those faces.

Physical derivatives follow from the chain rule:

.. math::

   \frac{\partial q}{\partial x} = q_\xi - \frac{z_x}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial y} = q_\eta - \frac{z_y}{z_\zeta}\,q_\zeta, \qquad
   \frac{\partial q}{\partial z} = \frac{q_\zeta}{z_\zeta} .

Ghost nodes outside the domain copy a domain node: the periodic image across a periodic side, or the edge node of
the same row across any other side. Below the ground and above the top the heights are extrapolated
(:math:`2z_{(0)} - z_{(1)}` below the ground).

Equations in conservation form
------------------------------

In the new coordinates the divergence of any flux :math:`\mathbf G` is

.. math::
   :label: curvilinear-divergence

   \nabla\cdot\mathbf G = \frac{1}{J}\left[\frac{\partial(\mathbf G\cdot\mathbf A_\xi)}{\partial\xi}
   + \frac{\partial(\mathbf G\cdot\mathbf A_\eta)}{\partial\eta}
   + \frac{\partial(\mathbf G\cdot\mathbf A_\zeta)}{\partial\zeta}\right],
   \qquad
   \mathbf A_\xi = (z_\zeta, 0, 0),\quad \mathbf A_\eta = (0, z_\zeta, 0),\quad \mathbf A_\zeta = (-z_x, -z_y, 1).

The vectors :math:`\mathbf A` are the face areas per unit index area. The x- and y-faces are vertical; the z-faces
follow the terrain. The flux of the velocity through a z-face is the contravariant velocity

.. math::

   \Omega = w - z_x\,u - z_y\,v ,

and continuity becomes

.. math::

   \frac{1}{J}\left[\frac{\partial(a_x u)}{\partial\xi} + \frac{\partial(a_y v)}{\partial\eta}
   + \frac{\partial\Omega}{\partial\zeta}\right] = 0 .

The momentum, :math:`\theta` and :math:`k` equations take the same form, with the convective flux
:math:`\mathbf u\,q` and the diffusive flux :math:`\gamma\nabla q` (or the stress :math:`\boldsymbol\tau`) as
:math:`\mathbf G`. The solution keeps the Cartesian velocity components :math:`(u, v, w)`; only the fluxes are
contravariant.

Two properties follow:

- **Free-stream preservation.** The metrics satisfy the geometric conservation law: a uniform wind has zero discrete
  divergence on any terrain.
- **Flat-ground limit.** With :math:`z_x = z_y = 0`, :math:`\Omega = w` and the equations reduce to the Cartesian
  ones.

Mesh quality report
-------------------

Every run prints a mesh-quality summary:

- the smallest cell and layer thickness relative to flat ground;
- the largest ground slope and coordinate-surface slope;
- the largest cell aspect ratio;
- the lowest first-cell centre and the largest displacement of the top nodes;
- the number of STF levels limited by :math:`\gamma_m`.

On restart, the node heights rebuilt from the inputs are checked against the checkpoint.
