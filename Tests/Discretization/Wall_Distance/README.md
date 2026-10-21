# Wall_Distance

## Purpose
The Tucker wall distance on ERF's Neutral_Hill_2D mesh: a 2-D Witch of Agnesi ridge (height 100 m, half-width
500 m), periodic in x and y. The test-only driver is `wall_distance` (`Supporting_Files/Tools/wall_distance.cpp`).

## Files
| File | Role |
|---|---|
| `inputs_walldist_hill2d` | ERF's Neutral_Hill_2D mesh. |
| `check_walldist.py` | Runs the driver, then ERF's wall-distance checks from `Neutral_Hill_2D/check_hill.py` against the exact distance to the ridge. |
| `Gold/walldist_hill2d_plt00000/` | The driver's plotfile. |

## CTest
`rt_walldist_hill2d`, `rt_walldist_hill2d_matches` (double precision only).

## Expected Results
- Maximum relative error at most 0.15, mean at most 0.03, at most 0.2 cells below 100 m, all values finite.
- The field matches `Gold/walldist_hill2d_plt00000` at relative 1e-8 (a gold written by GUST).

## How to run
```bash
ctest --test-dir <build> -R "^rt_walldist_hill2d"
```
