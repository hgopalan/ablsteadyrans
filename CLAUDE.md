# CLAUDE.md — ablsteady

Rules for working on **ablsteady** (`hgopalan/ablsteadyrans`), a stand-alone AMReX steady RANS solver for the
atmospheric boundary layer. This is the subset of the multi-code rules (written for ERF, AMR-Wind and related codes)
that applies here, plus what is specific to this repository.

- Development plan: `ablsteady_plan_rev6.md` (4 phases, 36 sub-phases, one PR each). Tests beyond the plan are
  welcome (mark them "beyond the plan"); variations from the plan go in a "Variations from the plan" list in the
  PR body, with the reason.
- PR tracker (order, status, inputs, tests): https://claude.ai/artifact/2DtP1MQHaQ91HFHAbPuWE6. Update the PR's row
  when its state changes.

## Tone and style

Every PR body, comment and review reply uses a neutral, factual tone.
- State what the code did and what it does now. No blame, superlatives or alarm words ("critical", "broken").
- Credit others by name and issue/PR number. Present related work constructively.
- Disclose costs and limits: performance loss, changed golds, what is "Not in this PR".
- Short sentences. The `## Summary` must be readable by a non-specialist.

### PR description layout
- **Line 1:** `Fixes #N` or `Related to #N` (or the plan sub-phase, e.g. `Plan rev. 6, sub-phase 1.3`).
- `## Summary`: problem, solution, relation to other PRs, bit-identity status.
- `## Theory` (physics or numerics only): LaTeX blocks defining every symbol.
- `## Changes`: numbered bold headings (`file:line`); an `Input | Default | Effect` table for every new input.
- `## Results`: compact tables; recommended configuration in bold.
- `## Test(s)`: exact CTest names and status; explain any expected change.
- `## Variations from the plan` (when any): each change from the plan and why.
- `## Not in this PR`: deferred items.

## Three review lenses (parallel agents)

No PR is marked ready until all three have checked the code paths (callers, callees, configurations not exercised)
and every finding is resolved.
1. **Ann (user-facing / maintainer):** inputs vs docs (names, defaults, units, ranges); keys read but ignored;
   defaults that change existing answers; restart completeness; vacuous tests.
2. **Aaron (implementation):** index arithmetic, periodic seams, domain faces, box splits; MPI consistency and
   collectives; GPU portability; memory lifetimes; single precision.
3. **Pressel (numerics and physics):** re-derive formulas in a short script; parameter sweeps for solvers and
   root-finders; stability of the pseudo-time iteration and of solver coupling.

## Workflow

### Before starting
- Check for other active worktrees or uncommitted edits on the same sub-phase.
- Pull `main`; check open PRs and issues.
- Work in a worktree (`.worktrees/<name>`), never in the main checkout. One build directory per worktree, named
  `build-<something>` (never plain `build` on a case-insensitive filesystem).

### Design and coding
- New capabilities live in their own class or module, which owns its diagnostics and checkpoint I/O.
- Reuse existing selectors; never add a parallel input for something an input already controls.
- Keep old behavior reachable behind a named option when a default changes; ship an old-vs-new regression test.
- Every input has a range check that **aborts** (never clamps) and names the input in the message.
- Unknown keys in the `mesh`, `terrain`, `init`, `turb`, `bc`, `steady`, `solver` and `io` namespaces abort.
  Add each new key to the reader in the same PR as the code that uses it.
- Document ghost-cell fills and NaN guards in the commit message.
- New options get their own `ParallelFor` kernel; no option branches inside device lambdas.
- Code copied from ERF lives in `extern_erf/` and gets a row in `extern_erf/EXTRACTED.md` (source path, ERF SHA,
  local changes). The drift check in CI fails on a copied file that differs from its record.

### Testing
- Unit tests (GoogleTest, `ctest -L unit`) are deterministic and run in seconds on 1–4 ranks.
- Regression tests compare against gold files; golds change only in a reviewed commit.
- **Mutant proofing:** each new test is shown to fail against the bug or a one-line mutant before it is committed.
- Pin closed-form or published numbers; plotfile existence or self-consistency is not a test.
- Test uneven decompositions (several boxes, 1/2/4 ranks) and non-uniform fields.
- A pre-existing bug found on the way: reproduce it on `main`, find the introducing commit, fix it in the PR.

### Portability (before pushing)
- No `this` or host-pointer captures in device lambdas; no device lambdas in private/protected members or `TEST_F`
  bodies; `Gpu::streamSynchronize()` before host reads of device data.
- Clean single-precision compile (`Real` literals, `amrex::Real` vs `double`).
- CPU runs that mimic GPU CI: `amrex.the_arena_init_size=0`, one box per rank
  (`fabarray.mfiter_tile_size = 1024 1024 1024`) and again with tiling, and
  `amrex.init_snan=1 amrex.fpe_trap_invalid=1`. AppleClang does not trap FPEs; use Linux/GCC for that check.

### Documentation and clean-up
- Detail goes in the Sphinx documentation (`docs/`); `README.md` stays short and points there.
- An input-table row in `docs/inputs.rst` for every new key, matching its range check.
- `make -C docs html` builds without warnings (`-W`) before pushing.
- `codespell`, no trailing whitespace, `git diff --check` clean.
- Never commit `docs/notes/`, plotfiles or checkpoints.

## Build

```bash
git submodule update --init --recursive
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
```
