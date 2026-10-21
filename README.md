# GUST

[![CI](https://github.com/hgopalan/gust-abl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/hgopalan/gust-abl/actions/workflows/ci.yml)
[![Docs](https://github.com/hgopalan/gust-abl/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/hgopalan/gust-abl/actions/workflows/docs.yml)

GUST, the Grid-sequenced Unified Steady Terrain solver, is a steady RANS solver for atmospheric boundary layer
simulations using [AMReX](https://github.com/AMReX-Codes/amrex). It computes steady wind fields over complex
terrain. It reuses [ERF](https://github.com/erf-model/ERF) components (Axell–Liungman k–L closure,
terrain-following mesh, terrain Poisson solver, MOST surface layer) and starts from a mass-consistent initializer
that builds on [FastWindTerrain](https://github.com/hgopalan/FastWindTerrain). It was called ablsteady before
October 2026, and the repository was called ablsteadyrans.

**Status:** phases 1 to 3 are merged: SIMPLEC with line LU-SGS and grid sequencing, and the Newton stage
(`steady.stage = simplec_newton`, #48). Phase 4 (FAS multigrid) was cancelled on 2026-10-21: the target tolerance
is r = 1e-4, which LU-SGS with grid sequencing already reaches quickly. Phase 5 (buildings) is planned. Details
are in the [documentation](https://hgopalan.github.io/gust-abl/).

## Quick start

```bash
git clone --recursive https://github.com/hgopalan/gust-abl.git
cd gust-abl
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
mpirun -n 2 build-rel/gust Tests/Infrastructure/Smoke_Restart/inputs_smoke
```

## Platforms

The CI badge covers every job below. Hosted runners have no GPU, so the GPU builds are compiled, not run
(details in the [CI page](https://hgopalan.github.io/gust-abl/ci.html)).

| Platform | CI |
|---|---|
| Linux CPU (GCC Release and Debug, Clang; MPI) | built; unit and smoke tests run |
| Windows CPU (MSVC, serial, no MPI) | built; unit and smoke tests run |
| Single precision (Linux) | compiled |
| CUDA sm_90 (Linux and Windows) | compiled |
| HIP gfx90a (Linux) | compiled |
| SYCL (Linux) | compiled |

GPU runs are not yet tested.

## Documentation

Building, inputs, output, the definition of steady, the case ladder, numerics, testing and CI are documented in
[`docs/`](docs/) (reStructuredText). The documentation of `main` is published at
<https://hgopalan.github.io/gust-abl/>, and each CI run keeps the built HTML as the `docs-html` artifact. To
build the HTML:

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html     # output in docs/_build/html
```

## Licence

BSD 2-Clause (`LICENSE`). Code copied from ERF keeps ERF's licence; see
[`docs/extraction.rst`](docs/extraction.rst).
