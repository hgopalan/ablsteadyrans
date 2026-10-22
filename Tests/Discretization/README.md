# Discretization

## Purpose
Cases that check one layer of the solver on a fixed state, before any iteration: the terrain mesh and its grid
hierarchy, the boundary
fluxes, the k-L closure, the Poisson layer and the wall distance. Each compares against a closed form, ERF's own
checks, or a gold written by GUST.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Terrain_Mesh` | `rt_askervein_mesh`, `rt_steep_ridge_quality`, `rt_hill2d_matches` | The Askervein mesh, a steep ridge under STF and the 2-D hill: statistics, quality report, mapped-grid plotfile. |
| `Grid_Hierarchy` | `rt_askervein_hierarchy`, `rt_hierarchy_restart`, `rt_steep_ridge_hierarchy_quality` | The horizontal-only grid hierarchy on the Askervein mesh and a steep ridge: per-level fields, quality reports, restart, invalid inputs. |
| `Boundary_Fluxes` | `rt_bc_neutral_flat` | The surface layer and slip top on a fixed log-law state; the residual against its closed form. |
| `Closure_Column` | `rt_closure_column` | The full residual of the equilibrium log-law column with the k-L closure. |
| `Poisson_Askervein` | `rt_poisson_askervein`, `rt_poisson_askervein_mlmg` | A fixed wind projected on the Askervein mesh with back ends A and B (MLMG). |
| `Wall_Distance` | `rt_walldist_hill2d` | The Tucker wall distance on ERF's Neutral_Hill_2D mesh. |

The drivers `poisson_projection`, `boundary_fluxes`, `closure_column` and `wall_distance` are test-only
executables in `Supporting_Files/Tools/`.
