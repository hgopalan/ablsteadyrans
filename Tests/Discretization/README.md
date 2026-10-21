# Discretization

## Purpose
Cases that check one layer of the solver on a fixed state, before any iteration: the terrain mesh and its grid
hierarchy, the boundary
fluxes, the k-L closure, the Poisson layer and the wall distance. Each compares against a closed form, ERF's own
checks, or a gold written by GUST.

## Cases
| Folder | Plan ID | Description |
|---|---|---|
| `Terrain_Mesh` | RT-1.4.1 to RT-1.4.3 | The Askervein mesh, a steep ridge under STF and the 2-D hill: statistics, quality report, mapped-grid plotfile. |
| `Grid_Hierarchy` | RT-2.3.1, RT-2.3.2 | The horizontal-only grid hierarchy on the Askervein mesh and a steep ridge: per-level fields, quality reports, restart, invalid inputs. |
| `Boundary_Fluxes` | RT-1.8.1 | The surface layer and slip top on a fixed log-law state; the residual against its closed form. |
| `Closure_Column` | RT-1.9.1 | The full residual of the equilibrium log-law column with the k-L closure. |
| `Poisson_Askervein` | RT-1.10.1, RT-1.11.1 | A fixed wind projected on the Askervein mesh with back ends A and B (MLMG). |
| `Wall_Distance` | RT-1.11.2 | The Tucker wall distance on ERF's Neutral_Hill_2D mesh. |

The drivers `poisson_projection`, `boundary_fluxes`, `closure_column` and `wall_distance` are test-only
executables in `Supporting_Files/Tools/`.
