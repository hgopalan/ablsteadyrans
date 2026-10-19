# Neutral_Hill_3D

## Purpose
ERF's Neutral_Hill_3D run to steady: 64 x 64 x 20, the hill of revolution, damping from 550 m, 16 sweeps, CFL
cap 4. It takes minutes on four ranks, so it carries the `nightly` label and is not in `ctest -L regression`.

## Files
| File | Role |
|---|---|
| `inputs_steady_hill3d` | The deck. |

Check script: `Supporting_Files/check_steady_hill.py --erf-check check_hill3d.py --orders 2` (shared with
`Neutral_Hill_2D`), which runs ERF's `check_hill3d.py --physics` (`extern_erf/`).

## CTest
`rt_steady_hill3d` (label `nightly`, four ranks, double precision only).

## Expected Results
- The steady verdict with a two-order drop of the u residual.
- ERF's 3-D checks, including the speed-up against 1.6 h/L.
- Measured: 2092 iterations (148 s on four ranks), every check passing.

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_hill3d$"
```
