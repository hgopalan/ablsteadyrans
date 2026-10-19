Development
===========

The development plan (rev. 6, October 2026) splits the work into four phases and 36 sub-phases, one pull request
each. Each phase ends at a gate, and the next phase starts only after the gate passes.

Pull requests
-------------

``CLAUDE.md`` at the top of the repository holds the rules. In short:

- one worktree and one build directory per pull request;
- every new input has a range check that stops the run and names the input, and a row in :doc:`inputs`;
- every test is shown to fail on a mutant before it is committed;
- three reviews before a pull request is ready: user-facing and maintainer, implementation, and numerics and
  physics;
- the description follows a fixed layout (Summary, Theory, Changes, Results, Tests, Not in this PR) in a neutral,
  factual tone;
- the documentation builds without warnings.

Notes on floating-point traps
-----------------------------

The CI runs with ``amrex.fpe_trap_invalid = 1``. Two places need care:

- **Closure kernels.** A compiler evaluates both arms of a select, so a square root meant for the other branch can
  raise the invalid trap on a value that is never used. The kernels take their branch decisions once and feed every
  arm a valid argument, written as :math:`\tfrac12(x \pm |x|)`. The results stay bit-identical to ERF's.
- **AMReX** ``MLTerrainPoisson`` **(not fixed here).** Its fused CPU red-black smoother is vectorized over a
  stride-2 loop and reads scratch lanes that are never written, which can raise the invalid trap on wide boxes. The
  multigrid back end therefore switches the trap off around the solve and checks the solution for non-finite values.
  Zeroing the scratch array in AMReX would remove the need.

Layout
------

.. code-block:: text

   Source/          application and ablsteady_refavg: inputs, state, I/O, reductions,
                    steady criterion, reference archives
   Tests/           one folder per test case (Tests/README.md), the GoogleTest suite
                    (Tests/Unit_Tests/) and shared check scripts (Tests/Supporting_Files/)
   docs/            this documentation (Sphinx)
   extern_erf/      files copied from ERF, listed in EXTRACTED.md
   tools/           drift check, VERSIONS check, Kestrel scripts (tools/kestrel/)
   external/        AMReX and GoogleTest submodules
