# Closure_Column

## Purpose
The full residual of the equilibrium log-law column with the k-L closure, on ERF's Neutral_ABL_Flat mesh
(8 x 8 x 64, neutral): u* = 0.35 m/s, z0 = 0.1 m, k = u*^2/Cmu0^2, the uncapped geometric length and the
logarithmic face mean. The test-only driver `closure_column` (`Supporting_Files/Tools/closure_column.cpp`)
evaluates it once.

## Files
| File | Role |
|---|---|
| `inputs_closure_column` | The equilibrium column. |
| `check_closure.py` | Runs the driver and checks the closure fields and residuals against their closed forms. |
| `Gold/closure_column_plt00000/` | The driver's plotfile. |

## CTest
`rt_closure_column`, `rt_closure_column_matches` (double precision only).

## Expected Results
- nu_t = kappa u* (z + z0) and L = kappa (z + z0) at every cell (1e-9); u* recovered (1e-10).
- R_u and R_v zero below the top cell (1e-10 of u*^2/dz); the top cell loses u*^2/dz through the slip top (1e-8).
- P_s, epsilon and R_k = -(P_s - epsilon) equal their closed forms from the discrete strain (1e-9 of epsilon).
  R_theta and R_w are zero.
- The plotfile matches `Gold/closure_column_plt00000` at relative 1e-10.

## Gold files
| Gold | How it was made |
|---|---|
| `closure_column_plt00000/` | `closure_column inputs_closure_column` |

## How to run
```bash
ctest --test-dir <build> -R "^rt_closure_column"
```
