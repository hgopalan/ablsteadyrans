# Supporting_Files

## Purpose
Files shared by more than one case folder, and the test-only drivers.

## Files
| File | Used by |
|---|---|
| `expect_abort.py` | `Infrastructure/Smoke_Restart`, `Infrastructure/Input_Validation`, `Infrastructure/Reference_Archive`, `Discretization/Grid_Hierarchy`, `Discretization/Poisson_Askervein`, `Canonical_RANS/Truly_Neutral_Hill_3D`, `Real_Terrain/*`: runs a command and requires each regex in its output; with `--success`, the command must also exit 0. |
| `check_mesh.py` | `Discretization/Terrain_Mesh`, `Discretization/Grid_Hierarchy`: runs the application, then compares field statistics of `z` and `detJ` and the quality report of every level with a JSON gold. |
| `check_init.py` | `Initialization/Mass_Consistent_Askervein` (mode `askervein`), `Initialization/Bump_Speedup_Trend` (mode `trend`). |
| `check_steady_hill.py` | `Canonical_RANS/Neutral_Hill_2D`, `Canonical_RANS/Neutral_Hill_3D`: the verdict, the residual drop, then ERF's `check_hill.py` or `check_hill3d.py` from `extern_erf/`. |
| `Tools/poisson_projection.cpp` | `Discretization/Poisson_Askervein` |
| `Tools/boundary_fluxes.cpp` | `Discretization/Boundary_Fluxes` |
| `Tools/closure_column.cpp` | `Discretization/Closure_Column` |
| `Tools/wall_distance.cpp` | `Discretization/Wall_Distance`, `Discretization/Grid_Hierarchy` (`abort_hierarchy_tool`) |
