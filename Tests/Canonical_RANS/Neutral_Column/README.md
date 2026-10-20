# Neutral_Column

## Purpose
The application on the 1-D neutral column: 4 x 4 x 64, 1000 m, surface layer z0 = 0.1 m, slip top, uniform
theta, driven by the force F = u*^2/H with u* = 0.35 m/s. The first guess is the log law with u* = 0.25 m/s; CFL
cap 50. The run goes to the steady verdict (6 orders, 1e-8).

A second deck runs the same column for 40 iterations, continues it to 80 from the checkpoint, and compares with
80 straight iterations.

## Files
| File | Role |
|---|---|
| `inputs_steady_column` | The column to steady. |
| `inputs_steady_column_restart` | The restart chain (iterations and first guess on the command line). |
| `inputs_steady_column_rough` | *Beyond the plan.* The column on rough ground: z0 = 0.5 m, 1 m layers, 64 m. |
| `check_steady_column.py` | Runs the application and checks the verdict, residual drop, force balance and wall law. |

## CTest
`rt_steady_column`, `rt_steady_column_constants`, `rt_steady_column_rough`, `rt_steady_column_rough_floor_one` (double precision only); `rt_steady_column_restart_first`, `rt_steady_column_restart_continued`,
`rt_steady_column_restart_straight`, `rt_steady_column_restart_matches` (both precisions). `abort_repeated_turb_key` and
`abort_repeated_surface_key` run the column deck with a repeated `turb.*` or `bc.surface_*` key, which must stop the run.

## Expected Results
- The steady verdict; the u residual down 6 orders from its reference (the Driver history of the checkpoint).
- Mean u* = sqrt(F H) to 1e-4 (the force balance); the first cell's speed the wall law's (1e-6).
- The second cell within 1 % of the log law, the lowest eight within 2 %; the wall k = u*^2/Cmu0^2 (1e-6);
  nu_t of the second cell within 6 % of kappa u* (z + z0); every column identical (1e-10); every field finite.
- Measured: 2231 iterations, u* = 0.3499999, second cell 0.9 % low, cells 2-7 1.3-1.9 % low.
- Restart: the continued and straight plotfiles identical under `amrex_fcompare` at relative 1e-14.
- Rough column, `log_law` floor (default): the same checks, with the cells above the first measured against the
  profile of the shallow column's stress u*^2 (1 - z/H) (`--stress-profile`). Measured: 347 iterations, first cell
  0.5977 m/s on the wall law, second cell 0.70 % low, lowest eight within 0.92 %.
- Rough column, `bc.surface_denominator_floor = one`: the first cell 1/ln 2 = 1.44 times the log law's speed
  ((u*/kappa) max(ln((z + z0)/z0), 1), 1e-6). Measured: 411 iterations, 0.8623 m/s.

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_column"
```
