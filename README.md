# ablsteady

[![CI](https://github.com/hgopalan/ablsteadyrans/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/hgopalan/ablsteadyrans/actions/workflows/ci.yml)
[![Docs](https://github.com/hgopalan/ablsteadyrans/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/hgopalan/ablsteadyrans/actions/workflows/docs.yml)

A stand-alone [AMReX](https://github.com/AMReX-Codes/amrex) solver for fast steady-state RANS of the atmospheric
boundary layer over complex terrain. It reuses [ERF](https://github.com/erf-model/ERF) components (Axell–Liungman
k–L closure, terrain-following mesh, terrain Poisson solver, MOST surface layer) and starts from a mass-consistent
initializer that builds on [FastWindTerrain](https://github.com/hgopalan/FastWindTerrain).

**Status:** phase 1 of 4 under way. Inputs, the terrain-following mesh (ERF's BTF, STF and Sullivan mappings),
I/O, restart, ERF reference archives, the definition of a steady solution, the Poisson layer (ERF's terrain
Poisson operator with GMRES and an FFT preconditioner, or AMReX's MLTerrainPoisson with multigrid) and the Tucker
wall distance are in place, and the SIMPLEC stage runs under the steady driver (`steady.stage = simplec`).

## Quick start

```bash
git clone --recursive https://github.com/hgopalan/ablsteadyrans.git
cd ablsteadyrans
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
mpirun -n 2 build-rel/ablsteady Tests/Infrastructure/Smoke_Restart/inputs_smoke
```

## Platforms

The CI badge covers every job below. Hosted runners have no GPU, so the GPU builds are compiled, not run
(details in the [CI page](https://hgopalan.github.io/ablsteadyrans/ci.html)).

| Platform | CI |
|---|---|
| Linux CPU (GCC Release and Debug, Clang; MPI) | built and tested |
| Windows CPU (MSVC, no MPI) | built and tested |
| Single precision (Linux) | compiled |
| CUDA, H100 (Linux and Windows) | compiled |
| HIP, AMD MI250X (Linux) | compiled |
| SYCL, Intel GPUs (Linux) | compiled |

## Documentation

Building, inputs, output, the definition of steady, the case ladder, numerics, testing and CI are documented in
[`docs/`](docs/) (reStructuredText). The documentation of `main` is published at
<https://hgopalan.github.io/ablsteadyrans/>, and each CI run keeps the built HTML as the `docs-html` artifact. To
build the HTML:

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html     # output in docs/_build/html
```

## Licence

BSD 2-Clause (`LICENSE`). Code copied from ERF keeps ERF's licence; see
[`docs/extraction.rst`](docs/extraction.rst).
