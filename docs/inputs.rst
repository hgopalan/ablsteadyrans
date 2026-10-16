Inputs
======

Inputs are ``namespace.key = value`` lines. ablsteady owns eight namespaces: ``mesh``, ``terrain``, ``init``,
``turb``, ``bc``, ``steady``, ``solver`` and ``io``.

The run stops, naming the input, when:

- a key in an owned namespace is not read by this version (a misspelling, or a key from a later sub-phase);
- a key's namespace is neither owned nor AMReX's (``amrex``, ``fab``, ``fabarray``, ``vismf``, ``tiny_profiler``,
  ``DistributionMapping``, ``device``, ``hypre``);
- a key is given twice, or a single-value key is given more than one value;
- a value is outside its range.

All errors are reported together. Each row below matches the range check in ``Source/Inputs.cpp``.

``mesh``
--------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``mesh.n_cell``
     - required
     - Three integers ≥ 1: cells in x, y and z.
   * - ``mesh.prob_lo``
     - ``0 0 0``
     - Three finite values: the domain's low corner (m).
   * - ``mesh.prob_hi``
     - required
     - Three finite values, each greater than ``prob_lo``: the domain's high corner (m).
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Three values, each 0 or 1. z must be 0 (ground and top are not periodic).
   * - ``mesh.max_grid_size``
     - ``64``
     - One or two integers ≥ 1: the largest box size in x and y. Boxes are never split in z, so the z-line solves
       stay inside one box.

``io``
------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``io.plot_file``
     - ``plt``
     - Plotfile prefix; not empty.
   * - ``io.chk_file``
     - ``chk``
     - Checkpoint prefix; not empty.
   * - ``io.restart``
     - empty
     - Checkpoint directory to restart from (see :doc:`running`).
   * - ``io.write_plotfile``
     - ``1``
     - 0 or 1.
   * - ``io.write_checkpoint``
     - ``1``
     - 0 or 1.

``steady``
----------

The definition of a steady solution (:doc:`steady`) has five parameters. They will be read as ``steady.*`` inputs
from sub-phase 1.13, when the steady driver exists to use them. Until then the defaults below are the definition,
and setting these keys stops the run, like any other key this version does not read.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Input (from sub-phase 1.13)
     - Default
     - Range and meaning
   * - ``steady.residual_tol``
     - ``5e-5``
     - Finite, > 0. Every scaled residual must be below it.
   * - ``steady.residual_drop_orders``
     - ``4``
     - ≥ 0 and ≤ the largest decimal exponent of the precision (308 in double, 38 in single). Every residual must
       also have dropped this many orders below its reference.
   * - ``steady.residual_ref_iterations``
     - ``5``
     - ≥ 1. A residual's reference is its largest value over this many first iterations.
   * - ``steady.monitor_rel_tol``
     - ``1e-3``
     - Finite, > 0. Every monitor must change by at most this fraction of max(\|value\|, scale) over the window, and be
       at most this far from its extrapolated limit.
   * - ``steady.monitor_window``
     - ``50``
     - ≥ 2. The window, in iterations.

The other namespaces (``terrain``, ``init``, ``turb``, ``bc``, ``solver``) have no inputs yet.
