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

Every ghost node of the node heights stands for a domain node: its periodic image across a periodic side, the
edge node of its row across any other side, at every height including the extrapolated rows below the ground
(:math:`2z_0 - z_1`) and above the top. A restart takes the ghost nodes from the rebuilt mesh.

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

Grid hierarchy
--------------

``mesh.hierarchy_levels`` = L > 1 adds L − 1 coarser meshes of the same domain (``Source/MeshHierarchy.H``). They
are for grid sequencing and nonlinear multigrid (later sub-phases); in this version they are built, reported and
written, and the solver runs on level 0.

- Level l has :math:`n_x / 2^l \times n_y / 2^l` cells: x and y are coarsened by 2 per level, **z never**. Every
  level has the z-levels of level 0, bit for bit, so the first cell and the wall treatment are the same on all of
  them.
- The ground of level l is level 0's ground at its nodes (node :math:`(I, J)` is level 0's node
  :math:`(2^l I, 2^l J)`), smoothed on level l's own grid by :math:`n_l = p\,(2^l - 1)` passes of the STF stencil

  .. math::

     h \leftarrow h + \beta\,(h_{i+1,j} + h_{i-1,j} + h_{i,j+1} + h_{i,j-1} - 4h), \qquad \beta = \tfrac18,

  with :math:`p` = ``terrain.hierarchy_passes``. Level 0 is not smoothed. A Fourier mode
  :math:`\cos(k_x x)\cos(k_y y)` is multiplied per pass by :math:`1 - \beta(4 - 2\cos k_x\Delta x - 2\cos k_y\Delta
  y) = \tfrac12[\cos^2(k_x\Delta x/2) + \cos^2(k_y\Delta y/2)] \in [0, 1]`: no mode changes sign, the 2-D
  checkerboard of the level's grid is removed in one pass and a 1-D odd-even mode is halved. After :math:`n_l`
  passes the filter is close to a Gaussian of width :math:`\sigma_l = \Delta x_0\,2^l\sqrt{p(2^l-1)/4}` (1.0, 3.5
  and 10.6 :math:`\Delta x_0` on levels 1–3 for p = 1), so the coarse grounds keep less relief per level: on the
  Askervein test mesh 124, 108 and 63 m on levels 0–2. The ground is sampled before it is smoothed, so fine scales
  above the coarse Nyquist wavenumber fold onto coarse ones; on Askervein this differs from smoothing on the fine
  grid first by 0.3 m RMS. A constant ground is unchanged.
- Aloft, each level follows its own ground with the run's ``terrain.smoothing`` and ``terrain.gamma_m``, and its
  metrics come from its own node heights, never from a finer level's. Each level is a full-domain mesh with its own
  Poisson solvers.
- Level l's boxes are level 0's coarsened by :math:`(2^l, 2^l, 1)`, on the same ranks. With L > 1, level 0's boxes
  are laid out so that this is exact: the coarsest domain cut into boxes of at most ``max_grid_size`` /
  :math:`2^{L-1}` cells, refined (``max_grid_size`` must be a multiple of :math:`2^{L-1}`). This layout can differ
  from the one without a hierarchy, so a level-0 solve can differ from an L = 1 run in the last bits. A level can also be gathered onto one box (``MeshHierarchy``'s ``gather_from``;
  no input yet), its node heights copied there with ``ParallelCopy``; the gathered level equals the distributed one
  bit for bit.
- A restart rebuilds the coarse levels from level 0's checkpointed node heights; nothing more is checkpointed, and
  ``mesh.hierarchy_levels`` and ``terrain.hierarchy_passes`` may change across a restart.

**Coarse cells over steep terrain.** The horizontal size doubles per level while the layers do not, so the aspect
ratio doubles too (on the Askervein test mesh: 8.5, 15.6, 31 on levels 0, 1, 2). The report prints every level's
quality and a summary table. When STF stops early on a coarse level at more z-levels than on level 0, it adds a
hint: raise ``terrain.hierarchy_passes``, lower ``terrain.gamma_m`` or use fewer levels. (Capping that level 0
shares is level 0's, and its own warning names it.) More passes reduce capping and slopes but raise the aspect
ratio, since STF then leaves thinner layers near the ground (Askervein level 2: aspect 25, 31, 34 for p = 0, 1,
2). A level of 1 × 1 cells is refused (the FFT preconditioner needs two directions longer than one cell). ERF's STF note, printed while each
level is built, names its own level 0 for every hierarchy level.

References
----------

- Klemp, J. B. (2011). A terrain-following coordinate with smoothed coordinate surfaces. *Mon. Weather Rev.* 139,
  2163–2169.
