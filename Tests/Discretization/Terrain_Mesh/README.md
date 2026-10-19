# Terrain_Mesh

## Purpose
The terrain-following mesh on three grounds, with no solve: ERF's Askervein terrain, a 600 m ridge steep enough
to trigger the STF warnings, and the 2-D Witch of Agnesi hill of ERF's Neutral_Hill_2D.

## Files
| File | Role |
|---|---|
| `inputs_mesh_askervein` | ERF's Askervein mesh coarsened to 60 x 60 x 18 (20 m x 1.08 layers, BTF). Terrain file on the command line: `terrain.file=extern_erf/askervein.txt`. |
| `inputs_mesh_steep_ridge` | A 600 m ridge under STF. |
| `inputs_mesh_hill2d` | The 2-D hill; writes the mapped-grid plotfile. |
| `Gold/` | `askervein_mesh.json`, `steep_ridge_mesh.json`, `hill2d_plt00000/`. |

Check script: `Supporting_Files/check_mesh.py` (shared with `Grid_Hierarchy`). It runs the application, then compares
field statistics of `z` and `detJ` and the quality report with a JSON gold (`--write` rewrites it).

## CTest
`rt_askervein_mesh`, `rt_steep_ridge_quality` (double precision only), `rt_hill2d_run`, `rt_hill2d_matches`.

## Expected Results
- Askervein: statistics and quality report match `Gold/askervein_mesh.json`.
- Steep ridge: the report matches its gold, including both warnings (STF capped, top not flat).
- 2-D hill: the plotfile matches `Gold/hill2d_plt00000` at relative 1e-12 (1e-5 in single precision). The
  plotfile holds cell fields only; UT-1.4.7 compares the nodal heights with ERF.

The three golds were written by ablsteady: they detect changes, they do not show correctness. Correctness of
these code paths comes from UT-1.4.7 (ERF's own output) and the closed-form unit tests.

## Gold files
| Gold | How it was made |
|---|---|
| `askervein_mesh.json` | `Supporting_Files/check_mesh.py Gold/askervein_mesh.json --write -- ablsteady inputs_mesh_askervein terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_mesh.json` | `Supporting_Files/check_mesh.py Gold/steep_ridge_mesh.json --write -- ablsteady inputs_mesh_steep_ridge` |
| `hill2d_plt00000/` | `ablsteady inputs_mesh_hill2d` |

`hill2d_plt00000/Level_0/z_phys_nd_D_00000` (the node heights with ghost nodes; no test reads it) was rewritten when
the ghost columns were made to match their periodic images: 8 values of the row below the ground changed.

## How to run
```bash
ctest --test-dir <build> -R "rt_askervein_mesh|rt_steep_ridge|rt_hill2d_"
```
