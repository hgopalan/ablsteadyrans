Grid hierarchy
==============

Iterative solvers remove short-wavelength errors quickly and long-wavelength errors slowly. Coarser grids help: a
long wave on the fine grid is a short wave on a coarse grid, and coarse iterations are cheap. GUST builds a
hierarchy of coarser meshes for this purpose. Grid sequencing (``steady.sequencing``, :doc:`solution_algorithm`)
uses it: it solves on the coarse meshes first and interpolates each result to the next finer mesh. Nonlinear
multigrid, in development, will also use it. The code is in
``Source/MeshHierarchy`` and ``Source/Transfer``.

Coarse levels
-------------

``mesh.hierarchy_levels`` = :math:`L > 1` adds :math:`L - 1` coarser meshes of the same domain.

- **Horizontal coarsening only.** Level :math:`l` has :math:`n_x/2^l \times n_y/2^l` cells. The vertical levels are
  never coarsened: every level has the vertical levels of level 0, so the first cell and the wall model are the same
  on all levels.
- **Smoothed ground.** The ground of level :math:`l` is level 0's ground sampled at its nodes, then smoothed on
  level :math:`l`'s own grid by :math:`n_l = p\,(2^l - 1)` passes of the filter

  .. math::

     h \leftarrow h + \beta\,(h_{i+1,j} + h_{i-1,j} + h_{i,j+1} + h_{i,j-1} - 4h), \qquad \beta = \tfrac18,

  where :math:`p` is ``terrain.hierarchy_passes``. Level 0 is not smoothed. The filter damps every Fourier mode
  without changing its sign and removes the grid-scale checkerboard in one pass. Coarser levels therefore keep less
  terrain relief; on the Askervein test mesh the relief is 124, 108 and 63 m on levels 0, 1 and 2.
- **Own mesh.** Each level follows its own ground with the run's ``terrain.smoothing`` and ``terrain.gamma_m``, and
  computes its own metrics and Poisson solvers.
- **Same ranks.** Level :math:`l`'s boxes are level 0's boxes coarsened by :math:`(2^l, 2^l, 1)`, on the same MPI
  ranks. For this, ``mesh.n_cell`` and ``mesh.max_grid_size`` in x and y must be divisible by :math:`2^{L-1}`, and
  level 0 is laid out so that every level coarsens it exactly. That layout can differ from the one without a
  hierarchy, so a level-0 result can differ from an :math:`L = 1` run in the last bits.
- **Restart.** The coarse levels are rebuilt from level 0. ``mesh.hierarchy_levels`` and
  ``terrain.hierarchy_passes`` may change across a restart.

**Coarse cells over steep terrain.** The horizontal cell size doubles per level while the vertical size does not, so
the aspect ratio doubles too. The run prints a mesh-quality report per level. If STF smoothing is limited on a coarse
level, the report suggests raising ``terrain.hierarchy_passes``, lowering ``terrain.gamma_m`` or using fewer levels.

Transfer operators
------------------

Two operators move a solution between a fine level :math:`l_f` and a coarse level :math:`l_c`, with ratio
:math:`r = 2^{l_c - l_f}` in x and y. Since z is never coarsened, both work level by level in the vertical. Grid sequencing uses the prolongation;
the restriction is for nonlinear multigrid.

**Restriction** (fine to coarse) conserves mass:

- :math:`\theta`, :math:`k` and :math:`p'` at cell centres: the volume-weighted mean of the :math:`r \times r` fine
  cells, :math:`q_c = \sum J_f q_f / \sum J_f`.
- :math:`u` and :math:`v`: the face flux. A coarse x-face covers :math:`r` fine x-faces:

  .. math::

     a_{x,c}\, u_c\, r = \sum_{b=0}^{r-1} a_{x,f}\, u_f .

- :math:`w`: the vertical flux :math:`\Omega` of the :math:`r^2` fine z-faces,
  :math:`a_{z,c}\,\Omega_c\, r^2 = \sum a_{z,f}\,\Omega_f`, then :math:`w_c` from :math:`\Omega_c` with the coarse
  slopes.

Each coarse face carries the total flux of the fine faces it covers. The coarse divergence is therefore the
volume-weighted mean of the fine one,

.. math::

   r^2 J_c D_c = \sum J_f D_f ,

and a divergence-free fine field restricts to a divergence-free coarse field, over flat ground and terrain alike.

**Prolongation** (coarse to fine):

#. Each field is interpolated bilinearly at its own staggered location, level by level. Next to a non-periodic
   domain face the interpolation is extended linearly. A field :math:`a + bx + cy + dxy` is reproduced exactly
   when each non-periodic direction has at least two coarse cells.
#. :math:`k` is bounded below by ``turb.tke_floor``.
#. The fine level's boundary conditions are applied.
#. A SIMPLEC pressure correction on the fine level makes the field divergence-free. The interpolated :math:`p'` is
   kept.

Bilinear interpolation is not conservative, so restriction after prolongation does not return the original coarse
field exactly; it applies a mild smoothing.

Over terrain the coarse ground is smoothed, so a prolonged field is an initial guess on the fine mesh, not an
interpolant of a fine solution.

The transfer operators have no inputs. Their tests are listed in :doc:`testing`.
