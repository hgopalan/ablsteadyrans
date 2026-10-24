# GUST

[![CI](https://github.com/hgopalan/gust-abl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/hgopalan/gust-abl/actions/workflows/ci.yml)
[![Docs](https://github.com/hgopalan/gust-abl/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/hgopalan/gust-abl/actions/workflows/docs.yml)

**GUST** (Grid-sequenced Unified Steady Terrain solver) computes steady wind fields over complex terrain. It is a
steady Reynolds-averaged Navier–Stokes (RANS) solver for the atmospheric boundary layer, built on
[AMReX](https://github.com/AMReX-Codes/amrex), and one code base builds for CPUs and GPUs (CUDA, HIP, SYCL).

## What it does

- **Two answers from one input file.** A mass-consistent wind in seconds (stop after the initializer), or the
  converged steady RANS solution in minutes.
- **Steady solver.** SIMPLEC pressure–velocity coupling with implicit vertical line relaxation (line LU-SGS),
  local pseudo-time stepping with switched-evolution relaxation, and grid sequencing over coarser horizontal
  levels. An optional Jacobian-free Newton–Krylov stage (`steady.stage = simplec_newton`) finishes to tight
  tolerances.
- **Physics.** The Axell–Liungman k–L turbulence closure, a Monin–Obukhov surface layer, buoyancy, Coriolis, and
  inflow from a log law or a measured profile.
- **Terrain.** A terrain-following mesh from an elevation file (`x y z` points or ERF's grid format), with
  stretched vertical layers for thin near-ground cells. The pressure equation uses a terrain-aware multigrid
  solver that handles high cell aspect ratios.
- **A clear definition of steady.** Residuals must fall by a set number of orders and the monitored quantities
  must settle; the run says why it stopped. Restarts continue bit for bit on the same ranks and box layout.

GUST reuses components of [ERF](https://github.com/erf-model/ERF) (the k–L closure, the terrain-following mesh and
Poisson operator, the surface layer), and its initializer builds on
[FastWindTerrain](https://github.com/hgopalan/FastWindTerrain). It was called ablsteady before October 2026.

## Quick start

```bash
git clone --recursive https://github.com/hgopalan/gust-abl.git
cd gust-abl
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel -L unit --output-on-failure
mpirun -n 2 build-rel/gust Tests/Infrastructure/Smoke_Restart/inputs_smoke
```

The [tutorials](https://hgopalan.github.io/gust-abl/) walk through a single column, a flat boundary layer, a hill
and Askervein Hill.

## Platforms

| Platform | Status in CI |
|---|---|
| Linux CPU (GCC Release and Debug, Clang; MPI) | built, unit and smoke tests run |
| Windows CPU (MSVC, serial) | built, unit and smoke tests run |
| Single precision (Linux) | compiled |
| NVIDIA CUDA sm_90 (Linux, Windows) | compiled |
| AMD HIP gfx90a (Linux) | compiled |
| Intel SYCL (Linux) | compiled |

Hosted CI runners have no GPUs, so the GPU builds are compiled but not yet run (see the
[CI page](https://hgopalan.github.io/gust-abl/ci.html)).

## Documentation

The full documentation is at <https://hgopalan.github.io/gust-abl/>: building and running, every input with its
range, the governing equations, discretization and solution algorithm, test cases and results. The sources are
in [`docs/`](docs/); to build them locally:

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html     # output in docs/_build/html
```

## Contributing

Issues and pull requests are welcome. The development rules (tests, input checks, review) are in
[`docs/development.rst`](docs/development.rst).

## Licence

BSD 2-Clause ([`LICENSE`](LICENSE)). Code copied from ERF keeps ERF's licence; see
[`docs/extraction.rst`](docs/extraction.rst).
