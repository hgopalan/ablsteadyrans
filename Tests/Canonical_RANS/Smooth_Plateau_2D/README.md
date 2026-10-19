# Smooth_Plateau_2D

## Purpose
A comparison case, not a CTest. A 2-D smooth plateau: a tanh ramp up from flat ground, a flat top and a ramp down
(100 m high, ramps of 250 m at +-750 m), with the settings of `Neutral_Hill_2D/inputs_steady_hill2d`. ERF runs the
same ground from the same file with its Neutral_Hill_2D deck minus the `prob.*` hill keys. The comparison is in
`docs/results.rst`.

## Files
| File | Role |
|---|---|
| `inputs_steady_ramp2d` | The deck. The terrain file is given on the command line as `terrain.file=ramp2d.txt`. |

The terrain file is written by `tools/make_ramp_terrain.py`.

## How to run
```bash
python3 tools/make_ramp_terrain.py ramp2d.txt
```
```bash
mpirun -n 2 <build>/ablsteady Tests/Canonical_RANS/Smooth_Plateau_2D/inputs_steady_ramp2d terrain.file=ramp2d.txt
```
