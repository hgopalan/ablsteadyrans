Running a simulation
====================

GUST takes one input file and, optionally, extra inputs on the command line:

.. code-block:: bash

   mpirun -n 2 build-rel/gust inputs_file [key=value ...]

For a first run, follow :doc:`tutorial_column`.

What a run does
---------------

#. **Check the inputs** (:doc:`inputs`). If any check fails, the run lists every error and stops.
#. **Build the mesh** and print its quality report (:doc:`terrain_coordinates`).
#. **Build the initial field** (:doc:`initial_conditions`), or read it from a checkpoint (``io.restart``).
#. **Solve** with the SIMPLEC stage when ``steady.stage = simplec`` (:doc:`solution_algorithm`), until the
   convergence criteria are met or ``steady.max_iterations`` is reached. One log line is printed per iteration for
   the first ten, then every tenth.
#. **Write the output**: a plotfile, the turbulence plotfile if ``io.turb_plot_file`` is set, and a checkpoint
   (:doc:`output`).

With the default ``steady.stage = none`` the run skips step 4: it writes the mesh and initial field and stops.
This is useful for checking a new mesh or initial field.

Command-line inputs
-------------------

Inputs after the file name are added to the file's inputs. A key that the file already sets cannot be given again
on the command line: a repeated key is refused.

Restarting
----------

Set ``io.restart`` to a checkpoint directory (:doc:`howto` has a step-by-step recipe).

- The mesh and terrain inputs and the floating-point precision must match the checkpoint exactly. The box layout
  (``mesh.max_grid_size``) and the number of ranks may differ.
- ``init.*`` keys cannot be combined with ``io.restart``: the checkpoint holds the state.
- With grid sequencing, a restart continues level 0 only.
- A restart with the SIMPLEC stage needs a checkpoint written by the stage. The checkpoint holds the solver history
  and the turbulence model's state, so the continued run reproduces an uninterrupted one on the same ranks and box
  layout: bit for bit over flat ground, to round-off over terrain (1.5e-12 measured on the 2-D hill).
- The run stops if it would write a checkpoint with the name of the one it read.

AMReX runtime flags
-------------------

AMReX reads its own ``amrex.*`` flags. The tests use these, which match the GPU CI settings:

.. code-block:: text

   amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1 amrex.abort_on_unused_inputs=1

``amrex.abort_on_unused_inputs=1`` makes a misspelled ``amrex.*`` flag fail the run when it finishes (after the
output is written). AppleClang builds do not trap
floating-point exceptions; use a Linux GCC build for that check.
