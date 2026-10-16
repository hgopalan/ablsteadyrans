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

Layout
------

.. code-block:: text

   Source/          application and ablsteady_refavg: inputs, state, I/O, reductions,
                    steady criterion, reference archives
   Tests/Unit/      GoogleTest unit tests
   Tests/Decks/     input decks for the smoke and abort tests
   docs/            this documentation (Sphinx)
   extern_erf/      files copied from ERF, listed in EXTRACTED.md
   tools/           drift check, VERSIONS check, Kestrel scripts (tools/kestrel/)
   external/        AMReX and GoogleTest submodules
