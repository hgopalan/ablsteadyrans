Running
========

.. code-block:: bash

   mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_smoke

The first argument is the inputs file; further ``key=value`` arguments override it. A run today:

#. reads and checks the inputs (:doc:`inputs`), and stops with every error listed if any check fails;
#. builds the mesh and the solution fields, or reads them from ``io.restart``;
#. writes a plotfile and a checkpoint (:doc:`output`).

There is no solver yet, so the run does not iterate.

Restarting
----------

Set ``io.restart`` to a checkpoint directory. The checkpoint's mesh and floating-point precision must match the
inputs exactly. The box layout (``mesh.max_grid_size``) and the number of ranks may differ. The run stops if it
would write a checkpoint with the same name as the one it read, since AMReX would rename that checkpoint to
``<name>.old.NNNN``.

AMReX runtime flags
-------------------

AMReX reads its own ``amrex.*`` flags. The tests use these, which match the GPU CI settings:

.. code-block:: text

   amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1 amrex.abort_on_unused_inputs=1

``amrex.abort_on_unused_inputs=1`` makes a misspelt ``amrex.*`` flag stop the run instead of printing a warning
at the end. AppleClang builds do not trap floating-point exceptions; use a Linux GCC build to check for them.
