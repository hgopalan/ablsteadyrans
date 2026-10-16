Mesh
====

ablsteady uses ERF's terrain-following mesh: a uniform computational grid :math:`(\xi, \eta, \zeta)` mapped to physical
heights :math:`z(\xi, \eta, \zeta)` stored at nodes (``z_phys_nd``). The horizontal coordinates are not distorted.
The code that builds the mapping and its metrics is ERF's, copied into ``extern_erf/`` (:doc:`extraction`).

Vertical levels
---------------

The reference levels :math:`z_0 = 0 < z_1 < \dots < z_{n_z}` are the node heights of the mesh over flat ground
(``Source/ZLevels.cpp``, after ERF's ``init_zlevels``):

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Inputs
     - Levels
   * - none (default)
     - Uniform, :math:`z_k = k\,\Delta\zeta` with :math:`\Delta\zeta = z_\mathrm{top}/n_z` and
       :math:`z_\mathrm{top}` = ``mesh.prob_hi`` z.
   * - ``mesh.dz0``
     - Geometric, :math:`z_k - z_{k-1} = \Delta z_0\,r^{k-1}`, with :math:`r` solved by Newton's method so that
       :math:`z_{n_z} = z_\mathrm{top}` exactly (ERF's ``build_geometric_zlevels``).
   * - ``mesh.dz0`` and ``mesh.stretch_ratio``
     - Geometric with the given :math:`r`; the top node is wherever the layers end (ERF's
       ``grid_stretching_ratio`` with ``initial_dz``, as in its Askervein deck).
   * - ``mesh.zlevels_file``
     - The :math:`n_z + 1` heights in the file.

Two tops are involved. ``mesh.prob_hi`` z is the *computational* top :math:`z_\mathrm{top}`: the computational grid
stays uniform, and the metric terms below use :math:`\Delta\zeta = z_\mathrm{top}/n_z` whatever the levels. The
*physical* top is the last level :math:`z_{n_z}`. The two are equal for uniform and solved-ratio levels; with
``mesh.stretch_ratio`` or ``mesh.zlevels_file`` the physical top is wherever the levels put it. The metrics stay
consistent either way, since :math:`J` and the areas are differences of node heights divided by the same
:math:`\Delta\zeta`. In the BTF and Sullivan formulas below, :math:`z_\mathrm{top}` stands for the physical top
:math:`z_{n_z}`, as in ERF.

Terrain-following coordinates
-----------------------------

With the ground height :math:`h(x, y)` in the bottom layer of nodes, the node heights follow one of ERF's three
mappings (``terrain.smoothing``):

**0, basic terrain following (BTF).**

.. math::

   z = \frac{h\,z_\mathrm{top} + (z_\mathrm{top} - h)\,z_k}{z_\mathrm{top}} ,

which follows the terrain all the way up and flattens only at the top.

**1, smoothed terrain following (STF, Klemp 2011), the default.** The terrain imprint at level :math:`k` is a
smoothed copy :math:`h_s` of the one below, weighted by a decay function:

.. math::

   z = z_k + A(z_k)\, h_s(x, y, k),
   \qquad
   A(z) = \cos^6\!\Bigl(\frac{\pi}{2}\frac{z}{z_H}\Bigr)\ (z < z_H),\quad 0 \text{ above},
   \qquad
   z_H = \frac{2.44\, h_\mathrm{max}}{1 - \gamma_m} .

:math:`h_s` is obtained by up to 50 passes of the five-point smoothing
:math:`h \leftarrow h + \beta (h_{i\pm1} + h_{j\pm1} - 4h)` with :math:`\beta = 0.2 \min(z_k / 2h_\mathrm{max}, 1)`.
A pass is rejected, and the smoothing at that level stops, if it would make any layer thinner than
:math:`\gamma_m` (``terrain.gamma_m``, default 0.5) times its flat thickness. The mesh-quality report counts the
levels where this happened and warns: the terrain is then steep for the vertical grid.

**2, Sullivan.**

.. math::

   z = z_k + \Bigl(1 - \frac{z_k}{z_\mathrm{top}}\Bigr)^3 h .

On flat ground STF and Sullivan reproduce the levels :math:`z_k` exactly; BTF reproduces them to within the last
bit of each value (its formula reduces to :math:`z_\mathrm{top} z_k / z_\mathrm{top}`).

``Source/SmoothSurface.cpp`` provides the same five-point smoothing as a separate function. The mesh does not call
it yet (ERF's STF loop smooths inline); it is meant for the coarse levels of the grid hierarchy (sub-phase 2.3).

When :math:`z_H` exceeds the physical top, the decay function is not zero there and the top nodes follow the terrain.
The quality report gives the largest top-node displacement and warns; a later rigid-lid solver assumes a flat top.

The mesh is refused (the run stops, naming the inputs to change) if the ground reaches the physical top or if any
layer has zero or negative thickness. BTF and Sullivan can invert cells under terrain that is high relative to the
domain; STF cannot, since it keeps every layer at least :math:`\gamma_m` of its flat thickness.

Metric terms
------------

From the node heights (ERF's ``make_J``, ``make_areas``, ``make_zcc``):

.. math::

   J = \frac{1}{4\,\Delta\zeta} \sum_{\text{4 corners}} \bigl[z_{k+1} - z_k\bigr],
   \qquad
   a_x = \frac{1}{2\,\Delta\zeta}\bigl[(z_{k+1} - z_k)_{j} + (z_{k+1} - z_k)_{j+1}\bigr],

with :math:`a_y` analogous and :math:`a_z = 1`. :math:`J` is the cell volume divided by
:math:`\Delta x\,\Delta y\,\Delta\zeta`; :math:`a_x, a_y` are the x- and y-face area factors. The slopes of the
coordinate surfaces at z-faces, :math:`\partial z/\partial\xi` and :math:`\partial z/\partial\eta`, come from ERF's
inline helpers (``Compute_h_xi_AtKface``, ``Compute_h_eta_AtKface``).

These metrics satisfy the geometric conservation law: a uniform wind has zero discrete divergence on any terrain
(tested to 1e-12 of :math:`U/\Delta x` on round Witch-of-Agnesi and Gaussian hills, for all three mappings).

Ground
------

``terrain.type`` sets :math:`h`:

- ``flat``: :math:`h = 0`;
- ``agnesi``: ERF's Witch-of-Agnesi hill, :math:`h = H / (1 + (r/L)^2)`, centred in the domain, as a ridge in x, a
  ridge in y, or round (``terrain.hill_shape``);
- ``xyz``: a point cloud, interpolated by inverse distance weighting from the ``terrain.idw_neighbors`` nearest
  points with weight :math:`d^{-p}` (FastWindTerrain's method and defaults, 6 points, :math:`p = 2`). IDW returns
  the input points exactly but is not exact on a plane between them: on a regular point grid the error with 6
  neighbours is a fraction of the height change across one point spacing. ``terrain.idw_neighbors`` may not exceed
  the number of points;
- ``erf_grid``: ERF's structured terrain file, interpolated bilinearly exactly as ERF does it (including ERF's
  1 mm offset of the query point), so ablsteady and the ERF reference solutions see the same ground.

Ghost nodes outside a non-periodic domain take the height of the nearest domain node, as in ERF.

Mesh-quality report
-------------------

Every run prints:

.. code-block:: text

   ablsteady: mesh quality
     min cell thickness / flat layer thickness  = ...
     min nodal layer thickness / flat thickness = ...
     max ground height (m)                      = ...
     max ground slope (deg)                     = ...
     max coordinate-surface slope (deg)         = ...
     max aspect ratio max(dx,dy) / thickness    = ...
     min first-cell centre above ground (m)     = ...
     max top-node displacement (m)              = ...
     STF levels capped by terrain.gamma_m       = n of nz

The cell thickness comes from :math:`J` (a mean over the cell's four node columns); the nodal thickness is the
quantity :math:`\gamma_m` bounds. Slopes are :math:`\arctan|\nabla z|` from one-sided differences at the nodes, so
a ridge diagonal to the grid reports its full slope. A warning line follows when STF capped any level, and another
when the top nodes are not flat. The first-cell centre should lie at least 10 roughness lengths
above the ground; the roughness length arrives with the surface layer (sub-phase 1.7).

Restart
-------

A restart reads the node heights from the checkpoint, rebuilds them from the inputs and requires the two to agree
on every node; a restart with changed mesh or terrain inputs is refused. The metrics are then recomputed.
