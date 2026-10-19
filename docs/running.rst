Running
========

.. code-block:: bash

   mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_smoke

The first argument is the inputs file; further ``key=value`` arguments override it. A run today:

#. reads and checks the inputs (:doc:`inputs`), and stops with every error listed if any check fails;
#. builds the mesh and the solution fields (the first guess of :doc:`initial_conditions` with ``init.first_guess =
   profile``), or reads them from ``io.restart``;
#. with ``steady.stage = simplec``, runs the SIMPLEC stage under the steady driver (:doc:`solver`) to
   its verdict or ``steady.max_iterations``, printing one line per iteration (the first ten, then every tenth);
#. writes a plotfile, the turbulence diagnostics plotfile if ``io.turb_plot_file`` is set, and a checkpoint
   (:doc:`output`).

For example, the 1-D neutral column of RT-1.15.1 and ERF's Neutral_ABL_Flat:

.. code-block:: bash

   mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_steady_column
   mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_steady_neutral_flat

A key the deck already sets cannot be overridden on the command line: a repeated key is refused.

Restarting
----------

Set ``io.restart`` to a checkpoint directory. The checkpoint's mesh and floating-point precision must match the
inputs exactly. The box layout (``mesh.max_grid_size``) and the number of ranks may differ. The run stops if it
would write a checkpoint with the same name as the one it read, since AMReX would rename that checkpoint to
``<name>.old.NNNN``. With the SIMPLEC stage the checkpoint also holds the driver history and the closure's
:math:`\nu_t` history, and the continued run is bit for bit the uninterrupted one on the same rank count
(``rt_steady_column_restart_*``); a restart with the stage needs a checkpoint the stage wrote (its ``Driver`` file;
a checkpoint from a stage-less run is refused). ``init.first_guess = profile`` and ``io.restart`` are exclusive (a restart keeps its
state): the restart test gives the first guess on the command line of the first run only.

AMReX runtime flags
-------------------

AMReX reads its own ``amrex.*`` flags. The tests use these, which match the GPU CI settings:

.. code-block:: text

   amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1 amrex.abort_on_unused_inputs=1

``amrex.abort_on_unused_inputs=1`` makes a misspelt ``amrex.*`` flag stop the run instead of printing a warning
at the end. AppleClang builds do not trap floating-point exceptions; use a Linux GCC build to check for them.
