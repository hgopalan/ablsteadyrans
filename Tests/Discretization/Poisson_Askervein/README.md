# Poisson_Askervein

## Purpose
A fixed wind projected onto a divergence-free field on the Askervein mesh (30 x 30 x 18) by the Poisson layer.
Back end A (RT-1.10.1) and back end B, MLMG (RT-1.11.1), solve the same discrete operator, so both potentials
match the same gold. The test-only driver is `poisson_projection` (`Supporting_Files/Tools/poisson_projection.cpp`).

## Files
| File | Role |
|---|---|
| `inputs_poisson_askervein` | Back end A. Terrain file on the command line: `terrain.file=extern_erf/askervein.txt`. |
| `inputs_poisson_askervein_mlmg` | Back end B (MLMG). |
| `check_poisson.py` | Runs the driver, then checks iterations and divergence against a JSON gold (`--write` rewrites it). |
| `Gold/` | `poisson_askervein.json`, `poisson_askervein_mlmg.json`, `poisson_askervein_plt00000/`. |

## CTest
`abort_poisson_options` (both precisions); `rt_poisson_askervein`, `rt_poisson_askervein_matches`,
`rt_poisson_askervein_mlmg`, `rt_poisson_askervein_mlmg_matches` (double precision only).

## Expected Results
- Iterations at most 10 % above the gold; divergence after below 10 x the printed tolerance times before.
- The potential matches `Gold/poisson_askervein_plt00000` at relative 1e-8 for both back ends.
- `abort_poisson_options`: `solver.poisson_reltol` repeated and out of range are both named.

The golds were written by ablsteady, so they detect changes.

## How to run
```bash
ctest --test-dir <build> -R "poisson"
```
