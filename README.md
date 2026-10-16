# ablsteady

A stand-alone [AMReX](https://github.com/AMReX-Codes/amrex) solver for fast steady-state RANS of the atmospheric
boundary layer over complex terrain. It reuses [ERF](https://github.com/erf-model/ERF) components (Axell–Liungman
k–L closure, terrain-following mesh, terrain Poisson solver, MOST surface layer) and starts from a mass-consistent
initializer that builds on [FastWindTerrain](https://github.com/hgopalan/FastWindTerrain).

**Status:** phase 1 of 4 under way. Inputs, the terrain-following mesh (ERF's BTF, STF and Sullivan mappings),
I/O, restart, ERF reference archives and the definition of a steady solution are in place; there is no solver yet.

## Quick start

```bash
git clone --recursive https://github.com/hgopalan/ablsteadyrans.git
cd ablsteadyrans
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel -j 8
ctest --test-dir build-rel --output-on-failure
mpirun -n 2 build-rel/ablsteady Tests/Decks/inputs_smoke
```

## Documentation

Building, inputs, output, the definition of steady, the case ladder, numerics, testing and CI are documented in
[`docs/`](docs/) (reStructuredText; each CI run also keeps the built HTML as the `docs-html` artifact). To build
the HTML:

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html     # output in docs/_build/html
```

## Licence

BSD 2-Clause (`LICENSE`). Code copied from ERF keeps ERF's licence; see
[`docs/extraction.rst`](docs/extraction.rst).
