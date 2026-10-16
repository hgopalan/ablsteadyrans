# ablsteady

**ablsteady** is a stand-alone [AMReX](https://github.com/AMReX-Codes/amrex) application for fast steady-state RANS
of the atmospheric boundary layer (ABL) over complex terrain. It reuses [ERF](https://github.com/erf-model/ERF)
components where they fit: the Axell–Liungman k–L closure, the stretched terrain-following mesh with terrain
smoothing, the terrain Poisson solver, the MOST surface layer and the column tridiagonal solver. The starting field
comes from a mass-consistent (WindNinja-style) initializer that builds on
[FastWindTerrain](https://github.com/hgopalan/FastWindTerrain).

The steady solver is built in four phases, each ending with a faster solver than the one before:

| Phase | Solver at the end of the phase |
|---|---|
| 1 | Stand-alone app + pseudo-time SIMPLEC with implicit z (z-line Jacobi) |
| 2 | Line LU-SGS sweeps + x/y grid hierarchy and grid sequencing |
| 3 | Jacobian-free Newton–Krylov finish, preconditioned by one SIMPLEC iteration |
| 4 | Nonlinear FAS multigrid, GPU tuning, scaling, benchmark, validation and release |

Target machine: NLR Kestrel (CPU nodes with 104 Sapphire Rapids cores; GPU nodes with 4 NVIDIA H100).

## Status

Sub-phase **1.3** (application skeleton, build, CI and I/O) is in progress. The application reads and checks its
inputs, builds the mesh and the C-grid state, and writes a plotfile and a checkpoint, or restarts from a checkpoint.
There is no solver yet. The mesh is flat and uniform until sub-phase 1.4.

## Build

Requirements: CMake ≥ 3.25, a C++20 compiler, MPI, and FFTW (for AMReX's FFT module). AMReX and GoogleTest are git
submodules.

```bash
git clone --recursive https://github.com/hgopalan/ablsteadyrans.git
cd ablsteadyrans
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
```

| CMake option | Default | Effect |
|---|---|---|
| `ABLSTEADY_MPI` | `ON` | Build with MPI |
| `ABLSTEADY_OMP` | `OFF` | Build with OpenMP |
| `ABLSTEADY_PRECISION` | `DOUBLE` | `DOUBLE` or `SINGLE`. The steady stop needs `DOUBLE`; `SINGLE` is kept compiling only |
| `ABLSTEADY_GPU_BACKEND` | `NONE` | `NONE` or `CUDA` (add `-DAMReX_CUDA_ARCH=9.0` for H100) |
| `ABLSTEADY_FFT` | `ON` | AMReX FFT, used by the Poisson layer from sub-phase 1.10 |
| `ABLSTEADY_HYPRE` | `OFF` | AMReX with hypre, for the optional Newton preconditioner (sub-phase 3.3) |
| `ABLSTEADY_TESTING` | `OFF` | `AMReX_TESTING` and assertions on, as in GPU CI |
| `ABLSTEADY_ENABLE_TESTS` | `ON` | Build the unit tests and register all tests with CTest |
| `ABLSTEADY_MPIEXEC_PREFLAGS` | empty | Extra `mpiexec` flags for the tests (CI uses `--oversubscribe`) |
| `ABLSTEADY_USE_INTERNAL_AMREX` | `ON` | Use the AMReX submodule; `OFF` uses an installed AMReX via `AMReX_ROOT` |

On Kestrel, `tools/kestrel/build_kestrel.sh cpu` or `tools/kestrel/build_kestrel.sh gpu` loads the modules from ERF's
Kestrel guide and builds. Check module names with `module avail` first.

## Run

```bash
mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_smoke
```

### Inputs

Inputs are `namespace.key = value`. The namespaces are `mesh`, `terrain`, `init`, `turb`, `bc`, `steady`, `solver`
and `io`. A key in one of them that this version does not read aborts the run, and so does a key in any other
namespace except AMReX's own (`amrex`, `fab`, `fabarray`, `vismf`, `tiny_profiler`, `DistributionMapping`, `device`,
`hypre`). A key given twice, or a single-value key given more than one value, also aborts. Every range check aborts
and names the input. All errors are reported together.

| Input | Default | Range / meaning |
|---|---|---|
| `mesh.n_cell` | required | Three integers ≥ 1: cells in x, y, z |
| `mesh.prob_lo` | `0 0 0` | Three finite values: domain low corner (m) |
| `mesh.prob_hi` | required | Three finite values, each greater than `prob_lo`: domain high corner (m) |
| `mesh.is_periodic` | `0 0 0` | Three values, each 0 or 1; z must be 0 |
| `mesh.max_grid_size` | `64` | One or two integers ≥ 1: largest box size in x and y. Boxes are never split in z |
| `io.plot_file` | `plt` | Plotfile prefix, not empty |
| `io.chk_file` | `chk` | Checkpoint prefix, not empty |
| `io.restart` | empty | Checkpoint directory to restart from. Its mesh and precision must match exactly; the box layout and rank count may differ. The run aborts if it would write a checkpoint of the same name |
| `io.write_plotfile` | `1` | 0 or 1 |
| `io.write_checkpoint` | `1` | 0 or 1 |

### Output

- **Plotfile** `<io.plot_file>NNNNN`: cell-centred `u v w p theta tke z`. The velocities are face averages; `z` is the
  mean of the cell's eight node heights. The nodal heights themselves, ghost nodes included, are in
  `Level_0/z_phys_nd` (a VisMF MultiFab).
- **Checkpoint** `<io.chk_file>NNNNN`: a text `Header` (format version, bytes per real, step, mesh with hex-float
  coordinates, BoxArray) and each field (`u v w p theta tke z_phys_nd`) with its ghost cells. A restart restores
  valid values exactly and the ghost nodes of `z_phys_nd`; other ghost cells are refilled by the solver.

## Tests

`ctest -L unit` runs in seconds; `ctest -L smoke` runs the application.

| Test | What it checks |
|---|---|
| `unit_np1`, `unit_np2`, `unit_np4` | GoogleTest suite on 1, 2 and 4 ranks: input parsing and unknown-key detection (UT-1.3.1), plotfile round trip bit-identical (UT-1.3.2), checkpoint/restart round trip bit-identical onto another box layout (UT-1.3.3), sum and max identical to the last bit on every rank count for a fixed box layout (UT-1.3.4) |
| `unit_notiling` | The same suite with tiling off, as on a GPU |
| `smoke_clean`, `smoke_run`, `smoke_restart`, `smoke_restart_matches` | Run on 2 ranks; restart on 3 ranks with another box layout (the log must show the restart); the two plotfiles are identical (`amrex_fcompare`, zero tolerance) |
| `abort_unknown_key`, `abort_range`, `abort_repeated_key`, `abort_restart_mesh`, `abort_restart_overwrite` | Decks that must abort, with the input named in the message |
| `drift_check_selftest` | The `extern_erf` drift check (RT-1.3.1) on a throw-away ERF repository |

Every run of `ablsteady` and `unit_tests` uses `amrex.the_arena_init_size=0 amrex.init_snan=1
amrex.fpe_trap_invalid=1`, as GPU CI does, plus `amrex.abort_on_unused_inputs=1` so a misspelt flag cannot switch a
check off unnoticed. Each test was shown to fail on a one-line mutant of the code it covers. In a build without MPI, every test runs on
one process.

## Continuous integration

`.github/workflows/ci.yml` runs on every pull request:

- GCC (Release, and Debug with `ABLSTEADY_TESTING`) and Clang, with Open MPI: build, unit and smoke tests, clean
  working tree;
- single precision, compile only;
- CUDA for sm_90 (H100), compile only;
- the `extern_erf` drift check, with a report of upstream ERF changes;
- codespell and trailing whitespace.

Kestrel runs the `--physics` and performance tests nightly once those exist (sub-phase 1.1 onward).

## Repository layout

```
Source/          application: inputs, state, I/O, reductions
Tests/Unit/      GoogleTest unit tests
Tests/Decks/     input decks for the smoke and abort tests
extern_erf/      files copied from ERF, listed in EXTRACTED.md
tools/           drift check (check_extern_erf.py) and Kestrel build script
external/        AMReX and GoogleTest submodules
```

## Development

The development plan (rev. 6) splits the work into 36 sub-phases, one pull request each. `CLAUDE.md` holds the rules
for pull requests, reviews and tests.

## Licence

BSD 2-Clause; see `LICENSE`. Files in `extern_erf/` keep ERF's licence.
