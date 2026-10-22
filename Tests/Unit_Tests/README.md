# Unit_Tests

## Purpose
The GoogleTest suite (`unit_tests`), deterministic and seconds long, run on every commit. Each `test_*.cpp` covers
one layer of `Source/`; the tests (GoogleTest names `Suite.Test`) and what each checks are listed in
`docs/testing.rst`.

## Files
| File | Role |
|---|---|
| `main.cpp` | AMReX initialization for GoogleTest. |
| `*TestHelpers.H` | Shared fixtures for the mesh, Poisson and transport tests. |
| `test_*.cpp` | The tests. |
| `Gold/` | ERF's mesh output for `Mesh.MetricsEqualErf` (see `Gold/README.md`). |

## CTest
`unit_np1`, `unit_np2`, `unit_np4` (1, 2 and 4 ranks in an MPI build), `unit_notiling` (tiling off, as on a GPU).

## How to run
```bash
ctest --test-dir <build> -L unit
```
