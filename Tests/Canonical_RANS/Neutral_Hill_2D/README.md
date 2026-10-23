# Neutral_Hill_2D

## Purpose
ERF's Neutral_Hill_2D run to steady: 128 x 1 x 64, the Witch of Agnesi of height 100 m and half-width 500 m,
the flat deck's sounding and forcing, Rayleigh damping of w and theta from 700 m. Three decks run the same case
with different solver settings and first guesses.

## Files
| File | Role |
|---|---|
| `inputs_steady_hill2d` | Jacobi sweeps (16 sweeps, CFL cap 4), log-law first guess. |
| `inputs_steady_hill2d_lusgs` | Line LU-SGS sweeps (two symmetric sweeps, CFL cap 16). |
| `inputs_steady_hill2d_column` | The LU-SGS deck started from the converged column (`init.first_guess = column`). |
| `inputs_steady_hill2d_restart` | The Jacobi deck for the restart chain (iterations and first guess on the command line). |
| `check_column_guess.py` | Compares the Driver histories of the LU-SGS profile run and the column run. |

The three runs use `Supporting_Files/check_steady_hill.py --erf-check check_hill.py --orders 2` (shared with
`Neutral_Hill_3D`): the verdict, a two-order drop of the u residual, then ERF's `check_hill.py --physics`.

## CTest
`rt_steady_hill2d`, `rt_steady_hill2d_lusgs`, `rt_steady_hill2d_column`, `rt_column_guess_spinup` (label
`regression`, double precision only). The LU-SGS and column runs share `<build>/Tests/run/rt_steady_hill2d_column`
so the spin-up check can read both histories. `rt_steady_hill2d_restart_first`, `rt_steady_hill2d_restart_continued`,
`rt_steady_hill2d_restart_straight`, `rt_steady_hill2d_restart_matches` (label `regression;smoke`, double
precision only): 5 iterations continued to 10 from the checkpoint against 10 straight, at relative 1e-10
(measured 1.5e-12; 8.5e-4 without the restart's ghost fill).

## Expected Results
- ERF's checks: the wall distance against the exact distance to the ridge, finite fields, k >= 0, the length-scale
  bound, the wall cell's k balance, the crest speed-up in the three lowest cells within 0.2-0.8 (Jackson and Hunt's
  0.4), the speed-up positive in the lowest ten, the upstream u* within 0.25-0.55 m/s, the upstream log law in the
  lowest three cells within 15 %.
- Jacobi: 4576 iterations (49 s on two ranks), speed-up 0.59 / 0.50 / 0.43, u* = 0.354 m/s.
- LU-SGS: 3048 iterations (23 s on two ranks).
- Column first guess: 1237 iterations; at iteration 300 the u and v residuals at most half the profile run's
  (measured 0.084 and 0.060).

## How to run
```bash
ctest --test-dir <build> -R "rt_steady_hill2d|rt_column_guess_spinup"
```
