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
     - Three finite values: the domain's low corner (m). The z value must be 0: node heights start at the flat
       ground.
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
   * - ``mesh.dz0``
     - ``0``
     - Finite, ≥ 0. First layer thickness of a stretched vertical grid; 0 means uniform layers. Without
       ``stretch_ratio`` it must not exceed prob_hi z / n_cell z (:doc:`mesh`).
   * - ``mesh.stretch_ratio``
     - ``0`` (none)
     - Finite, ≥ 1 (an explicit 0 is refused), needs ``mesh.dz0`` > 0. A fixed ratio between successive layers;
       the physical top is then wherever the layers end, and prob_hi z only sets the computational spacing.
   * - ``mesh.zlevels_file``
     - empty
     - A file of n_cell z + 1 increasing heights starting at 0. Not together with ``mesh.dz0``.

``terrain``
-----------

A key that has no effect for the chosen ``terrain.type`` or ``terrain.smoothing`` stops the run.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``terrain.type``
     - ``flat``
     - ``flat``, ``agnesi``, ``xyz`` or ``erf_grid`` (:doc:`mesh`).
   * - ``terrain.file``
     - required for ``xyz``, ``erf_grid``
     - The point file (``x y z`` per line) or ERF's structured terrain file.
   * - ``terrain.smoothing``
     - ``1``
     - 0 (BTF), 1 (STF, Klemp 2011) or 2 (Sullivan).
   * - ``terrain.gamma_m``
     - ``0.5``
     - > 0 and < 1; STF only. Smallest layer thickness as a fraction of the flat layer.
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Finite, > 0 (m); ``agnesi`` only.
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Finite, > 0 (m); ``agnesi`` only. The hill is at half height this far from its crest.
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - ``ridge_x``, ``ridge_y`` or ``round``; ``agnesi`` only.
   * - ``terrain.idw_neighbors``
     - ``6``
     - ≥ 1; ``xyz`` only. Points used for each height.
   * - ``terrain.idw_power``
     - ``2``
     - Finite, > 0; ``xyz`` only. Weight exponent, :math:`d^{-p}`.

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

``solver``
----------

The Poisson layer (:doc:`poisson`) reads ``solver.poisson_*`` through ``read_poisson_options()``. In this version only
the regression tool ``poisson_projection`` calls it: the application has no Poisson solve yet, so setting these keys
in an ablsteady deck stops the run, like any other key this version does not read. Each row matches the range check
in ``Source/Poisson.cpp``.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``solver.poisson_reltol``
     - ``1e-10`` (``1e-5`` single)
     - ≥ 10 machine epsilon and < 1. Converged when the residual 2-norm is below this fraction of the initial one.
   * - ``solver.poisson_abstol``
     - ``0``
     - Finite, ≥ 0. Or when the residual 2-norm is below this value.
   * - ``solver.poisson_max_iter``
     - ``200``
     - ≥ 1. GMRES iterations; reaching it stops the run.
   * - ``solver.poisson_restart``
     - ``50``
     - ≥ 1. GMRES restart length (ERF's value).
   * - ``solver.poisson_precond``
     - ``fft``
     - ``fft`` (ERF's FFT preconditioner on the flat stretched grid) or ``none``.
   * - ``solver.poisson_verbose``
     - ``0``
     - ≥ 0. 1 prints each solve's iterations, residuals and (singular problems) removed mean, and AMReX GMRES's solve
       time; 2 also GMRES's residual at each restart cycle.

The other namespaces (``init``, ``turb``, ``bc``) have no inputs yet.
