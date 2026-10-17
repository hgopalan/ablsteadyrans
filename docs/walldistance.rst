Wall distance
=============

The turbulence closures need the distance to the ground. ``Source/WallDistance.{H,cpp}`` provides Tucker's (2003)
differential wall distance, the algorithm of ERF's ``ERF::poisson_wall_dist`` lifted into a free function
(sub-phase 1.11):

.. code-block:: c++

   PoissonResult r = tucker_wall_distance(mesh, geom, z_phys_nd, dist, read_poisson_options());

Method
------

Solve, on the terrain-following mesh with the Poisson layer (:doc:`poisson`, back end A by default, full terrain cross
terms),

.. math::

   \nabla^2\phi = -1, \qquad \phi = 0 \text{ on the ground}, \qquad
   \partial\phi/\partial n = 0 \text{ through the top and non-periodic sides},

with periodic directions kept periodic, then set in every cell (Tucker 2003, Eq. 2)

.. math::

   d = -|\nabla\phi| + \sqrt{|\nabla\phi|^2 + 2\phi}.

:math:`\nabla\phi` is ERF's cell-centred physical gradient: centred differences in computational space and the
cell-centre metrics, :math:`\partial_z\phi = \partial_\zeta\phi/h_\zeta`,
:math:`\partial_x\phi = \partial_\xi\phi - h_\xi\partial_z\phi`, except in the first cell above the ground. There
:math:`\partial_z\phi` is the derivative of the quadratic along the column through :math:`\phi = 0` on the ground and
the two lowest cell centres, at heights :math:`a` and :math:`b` above the ground:

.. math::

   \partial_z\phi = \frac{\phi_0\, b\,(b - 2a) + \phi_1\, a^2}{a\,b\,(b - a)}
   \qquad\Bigl(= \frac{\phi_0 + \phi_1/3}{\Delta z} \text{ on uniform levels}\Bigr).

ERF uses the centred difference with the odd ghost, :math:`(\phi_1 + \phi_0)/(2\Delta z)` on uniform levels. Why the
change helps: on flat ground the discrete solution is the exact quadratic :math:`Hz - z^2/2` plus a constant
:math:`\Delta z^2/8` (from the odd ghost at the ground). The centred difference is then exact and the offset alone
gives a first-cell error of :math:`\Delta z/(4H)` of the height; the quadratic fit overshoots the gradient by
:math:`\Delta z/6`, which cancels part of the offset's effect and leaves :math:`\Delta z/(12H)`. On ERF's flat-fitted
deck that is 0.39 % → 0.13 %. Fitting along the physical column keeps the gain on stretched levels, where a fit in
:math:`\zeta` does not (first-cell error, 20 m first layer, 18 layers; a one-dimensional check made in review):

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - Stretch ratio
     - ERF centred
     - Fit in :math:`\zeta`
     - Fit along the column (used)
   * - 1.08
     - 1.26 %
     - 1.06 %
     - 0.20 %
   * - 1.2
     - 4.54 %
     - 3.14 %
     - 0.05 %

Above the first cell the error is the offset's, :math:`\Delta z^2/(8H)` in absolute terms, which converges at second
order (tested). The centred interior gradient divided by :math:`h_\zeta` is first order on stretched levels (0.7 % at
the second cell for ratio 1.2), unchanged from ERF.

For a flat ground :math:`d` is the height above it. Over terrain it approximates the distance, weakest at convex
crests.

ERF's code comment (``ERF_PoissonWallDist.cpp`` at 92410d57) reports that its ``MLABecLaplacian`` solve, with a
diagonal metric tensor and deferred correction of the cross terms, diverges on 3-D fitted meshes with
:math:`\Delta x \ne \Delta z`, and ERF uses ``terrain_height`` there. Here the operator carries the cross terms
directly; this work did not reproduce ERF's divergence.

Accuracy
--------

Against the exact distance, with ERF's own check tolerances (``Neutral_Hill_2D/check_hill.py``,
``Neutral_Hill_3D/check_hill3d.py``) on ERF's decks:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - Case (ERF deck)
     - Max relative error
     - Mean relative error
     - Max error, d < 100 m
   * - Flat ground (128 × 1 × 64)
     - 0.13 % (≤ 0.3 %)
     - —
     - 0.03 m (≤ 0.1 m)
   * - 2-D ridge, 100 m, L = 500 m (128 × 1 × 64)
     - 5.3 % (≤ 15 %)
     - 1.0 % (≤ 3 %)
     - 0.16 cells (≤ 0.2)
   * - 3-D hill, 100 m, L = 500 m (64 × 64 × 20)
     - 9.2 % (≤ 15 %)
     - 0.29 % (≤ 3 %)
     - 0.13 cells (≤ 0.2)

The 3-D values hold with periodic and with Neumann sides, and for both back ends. ERF's 3-D deck comment says
"dx = dy = dz = 40 m"; its cells are 80 m × 80 m × 40 m.

What the tolerances test: the ridge's 5.3 % maximum is at the domain top over the crest (exact distance 893 m, an
artefact of the Neumann top), not near the ground. The near-ground error is Tucker's own approximation, about 2.5 m
above the crest, and does not shrink with the mesh, so in cells it grows: 0.10, 0.16 and 0.33 cells on 64 × 32,
128 × 64 and 256 × 128. ERF's 0.2-cell tolerance therefore holds at ERF's deck resolution, not on finer meshes. The
first-cell error over the crest converges at second order.

In single precision the Tucker solve needs ``solver.poisson_reltol`` of about 1e-4 (the tests use it); the default
1e-5 is below the rounding floor (:doc:`poisson`).

Use
---

No part of the application uses the distance yet; the regression tool ``wall_distance`` computes it (RT-1.11.2).
The plan's choice between this and the height above the ground (``terrain_height``) arrives with the turbulence model.
