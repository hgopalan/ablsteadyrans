Numerics
========

Grid and fields
---------------

ablsteady uses ERF's Arakawa C-grid. The normal velocities ``u``, ``v``, ``w`` sit on x-, y- and z-faces.
The perturbation pressure p′, the potential temperature θ and the turbulent kinetic energy k sit at cell centres.
The physical height of each node, ``z_phys_nd``, is stored at nodes. Until terrain arrives (sub-phase 1.4) the
heights are those of a flat, uniform grid,

.. math::

   z_k = z_{\mathrm{lo}} + k\,\Delta z ,

computed as a fused multiply-add so every compiler and the GPU produce the same bits.

Boxes are split in x and y only: every box spans the full z column, so the implicit z-line solves of later
sub-phases stay inside one box.

Deterministic sum
-----------------

Residual norms and steady-state monitors are sums over the domain. A floating-point sum depends on the order of
its additions, so a sum in an order that changes with the number of ranks would give residual histories that
differ between runs. ``deterministic_sum`` (``Source/Reductions.H``) fixes the order:

.. math::
   :label: deterministic-sum

   S = \sum_{b=0}^{N_b-1} \Biggl( \sum_{j} \sum_{i} \Bigl( \sum_{k=k_{\mathrm{lo}}}^{k_{\mathrm{hi}}} a_{ijk} \Bigr) \Biggr)_b .

- :math:`a_{ijk}` is the cell value.
- :math:`b` is the box index in the BoxArray, and :math:`N_b` is the number of boxes.
- The innermost sum is one column, added in k order by one GPU thread.
- The column sums are added on the host in storage order (i fastest, then j).

Each box total is placed in its own slot of a vector of length :math:`N_b`, and the vector is summed over the
ranks. Each slot has a single non-zero contributor, and adding zeros is exact, so the order in which MPI combines
the ranks does not matter.

**Guarantee:** the same bits for any number of ranks and any DistributionMapping of a *fixed* BoxArray.

**Not guaranteed:** a different BoxArray (for example after a restart with another ``mesh.max_grid_size``) adds
in a different order, and the last bits can differ.

The sum accepts cell-centred data only: face and nodal boxes share their boundary points, which would be counted
twice.

``global_max`` returns NaN if any value is NaN, and returns zero as +0. Without these two rules the result could
depend on the decomposition, since ``max(x, NaN)`` and ``max(+0, -0)`` depend on argument order.
