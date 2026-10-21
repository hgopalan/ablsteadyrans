# Grid_Hierarchy

## Purpose
The horizontal-only grid hierarchy (sub-phase 2.3): coarser x/y levels over the same vertical layers, each with its
own smoothed ground. The case builds the hierarchy on the Askervein mesh and on a steep ridge. It checks each
level's `z` and `detJ` and the quality report of every level. It also checks that a restart rebuilds the
hierarchy, and that invalid hierarchy inputs stop the run.

## Files
| File | Role |
|---|---|
| `inputs_mesh_askervein_hierarchy` | The Askervein mesh of `Terrain_Mesh/inputs_mesh_askervein` with STF and three levels: 60 x 60, 30 x 30 and 15 x 15 cells, the same 18 layers; coarse grounds smoothed by 1 and 3 passes. Terrain file on the command line: `terrain.file=extern_erf/askervein.txt`. |
| `inputs_mesh_askervein_hierarchy_p0` | The same with unsmoothed coarse grounds (`terrain.hierarchy_passes = 0`). |
| `inputs_mesh_steep_ridge_hierarchy` | The steep ridge of `Terrain_Mesh/inputs_mesh_steep_ridge` on 32 x 4 x 16 cells with three levels (32, 16 and 8 cells in x). |
| `inputs_hierarchy_restart_first` | The Askervein mesh without a hierarchy; writes the checkpoint. |
| `inputs_hierarchy_restart` | Restarts that checkpoint with three levels on another box layout. |
| `inputs_abort_hierarchy` | n_cell not divisible by 4 for three levels, a box size not a multiple of 4, passes out of range: must stop. |
| `inputs_abort_hierarchy_unused` | `terrain.hierarchy_passes` without a hierarchy: must stop. |
| `inputs_abort_hierarchy_tool` | A hierarchy given to the test-only `wall_distance` tool, which works on level 0 alone: must stop. |
| `Gold/` | `askervein_hierarchy_mesh.json`, `askervein_hierarchy_p0_mesh.json`, `steep_ridge_hierarchy_mesh.json`, `askervein_hierarchy_level1_plt00000/`, `askervein_hierarchy_level2_plt00000/`. |

Check scripts: `Supporting_Files/check_mesh.py` (shared with `Terrain_Mesh`) and `Supporting_Files/expect_abort.py`.

## CTest
- `rt_askervein_hierarchy`, `rt_askervein_hierarchy_level1_matches`, `rt_askervein_hierarchy_level2_matches`
  (both precisions; the JSON check in double precision only).
- `rt_hierarchy_restart_first` (2 ranks), `rt_hierarchy_restart` (3 ranks), `rt_hierarchy_restart_level1_matches`,
  `rt_hierarchy_restart_level2_matches`.
- `rt_askervein_hierarchy_hint`, `rt_steep_ridge_hierarchy_quality` (double precision only).
- `abort_hierarchy`, `abort_hierarchy_passes_unused`, `abort_hierarchy_tool`.

All carry the labels `regression;smoke`.

## Expected Results
- RT-2.3.1 / RT-2.3.2, Askervein: level 0's `z` and `detJ` statistics and every level's quality report (per-level
  headers and the summary table included) match `Gold/askervein_hierarchy_mesh.json`. Levels 1 and 2 match their
  gold plotfiles at relative 1e-12 (1e-5 in single precision).
- Unsmoothed coarse grounds: levels 1 and 2 cap STF at 9 and 10 z-levels against level 0's 8, so the report gives
  the coarse-level hint.
- Steep ridge: the report of every level with the STF and top warnings. No coarse-level hint: levels 1 and 2 cap
  no more than level 0.
- Restart: the rebuilt levels 1 and 2 match the golds of `rt_askervein_hierarchy`.
- Each abort message names its input.

The golds were written by gust: they detect changes. The coarse grounds were checked against an independent
rebuild in the unit tests (UT-2.3.2).

## Gold files
| Gold | How it was made |
|---|---|
| `askervein_hierarchy_mesh.json` | `mpiexec -n 2 Supporting_Files/check_mesh.py Gold/askervein_hierarchy_mesh.json --write -- gust inputs_mesh_askervein_hierarchy terrain.file=extern_erf/askervein.txt` |
| `askervein_hierarchy_level1_plt00000/`, `askervein_hierarchy_level2_plt00000/` | The same run's `askervein_hier_plt_level{1,2}_00000`, without `Level_0/z_phys_nd*` |
| `askervein_hierarchy_p0_mesh.json` | `mpiexec -n 2 Supporting_Files/check_mesh.py Gold/askervein_hierarchy_p0_mesh.json --write -- gust inputs_mesh_askervein_hierarchy_p0 terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_hierarchy_mesh.json` | `Supporting_Files/check_mesh.py Gold/steep_ridge_hierarchy_mesh.json --write -- gust inputs_mesh_steep_ridge_hierarchy` |

## How to run
```bash
ctest --test-dir <build> -R "hierarchy"
```
