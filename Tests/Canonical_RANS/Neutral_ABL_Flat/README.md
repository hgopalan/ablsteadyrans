# Neutral_ABL_Flat

## Purpose
ERF's Neutral_ABL_Flat run to steady: 8 x 8 x 64, geostrophic wind 10 m/s at latitude 90 degrees
(f = 1e-4 1/s), the sounding's inversion from 700 m, slip top with 0.03 K/m, z0 = 0.1 m, ERF's closure defaults;
four sweeps, CFL cap 2, the log law with u* = 0.35 m/s as the first guess.

## Files
| File | Role |
|---|---|
| `inputs_steady_neutral_flat` | The deck. |
| `check_steady_neutral.py` | Requires the steady verdict, then runs ERF's `check_neutral.py --physics` (`extern_erf/`) on the turbulence plotfile. |

## CTest
`rt_steady_neutral_flat` (label `regression`, double precision only).

## Expected Results
- Every field finite, k >= 0, the capped kappa (z + z0) length (1e-3), the AL01 stability functions (1e-6),
  K_mv = c_mu sqrt(k) L (5 %), epsilon = Cmu0^3 k^1.5 / L (5 %), the wall cell's k balance (1 %).
- u* within 0.35 +- 0.15 m/s, the wall k within 5 % of u*^2/Cmu0^2, the log law in the lowest four cells within
  10 %, K_mv of the second cell within 30 % of kappa u* z.
- Measured: the verdict after about 2850 iterations (13 s on two ranks), u* = 0.377 m/s, every ERF check passing.

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_neutral_flat$"
```
