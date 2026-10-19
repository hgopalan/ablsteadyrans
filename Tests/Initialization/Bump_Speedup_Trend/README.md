# Bump_Speedup_Trend

## Purpose
The hilltop speed-up over FastWindTerrain's bump (its `regtests/phase6_solve_correction` hill and power-law
inflow) for alpha_v = 1, 0.7, 0.5, 0.3 and 0.1, compared with FastWindTerrain's own values at commit 648f8f1 on
the same hill and profile. One CTest runs all five cases on two ranks.

## Files
| File | Role |
|---|---|
| `inputs_init_bump` | The 1000 m x 1000 m domain with FastWindTerrain's hill. |
| `terrain_bump_fwt.txt` | The hill, given on the command line as `terrain.file=`. |
| `fwt_bump_speedup.py` | Writes the FastWindTerrain section of the gold from its plotfiles. |
| `Gold/init_bump_trend.json` | Both codes' speed-ups and ablsteady's iterations. |

Check script: `Supporting_Files/check_init.py trend` (shared with `Mass_Consistent_Askervein`).

## CTest
`rt_init_bump_trend` (double precision only).

## Expected Results
- The same ranking of the five alpha_v values at 40, 80 and 150 m above the ground in both codes (at 20 m two
  values tie within 0.002 in both codes).
- Every speed-up above 1 at 20, 40, 80 and 150 m.
- The two codes within 0.02, 0.01, 0.005 and 0.003 of each other at those heights (measured 0.011, 0.0057,
  0.0020, 0.0008).
- ablsteady's own values (1e-6) and iterations (10 %) pinned.

## Gold files
| Gold | How it was made |
|---|---|
| `init_bump_trend.json` (FastWindTerrain section) | `fwt_bump_speedup.py Gold/init_bump_trend.json SHA ALPHA=PLOTFILE ...` on FastWindTerrain's plotfiles |
| `init_bump_trend.json` (ablsteady section) | `check_init.py trend Gold/init_bump_trend.json --write -- ablsteady inputs_init_bump terrain.file=terrain_bump_fwt.txt` |

## How to run
```bash
ctest --test-dir <build> -R "^rt_init_bump_trend$"
```
