# Askervein

## Purpose
The convergence half of RT-1.16.3: Askervein on the ERF grid of RT-1.12.1 (30 x 30 x 18), the log-law inflow,
outflow, slip sides, no Coriolis, four sweeps, CFL cap 2, to the steady verdict.

## Files
| File | Role |
|---|---|
| `inputs_steady_askervein` | The deck. Terrain file on the command line: `terrain.file=extern_erf/askervein.txt`. |

Check script: `Supporting_Files/expect_abort.py --success "gust: steady after"` (shared).

## CTest
`rt_steady_askervein` (labels `regression;smoke`, double precision only).

## Expected Results
- The steady verdict. Measured: 175 iterations (2 s on two ranks; 156 with `bc.inflow_tke_length = unbounded`).

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_askervein$"
```
