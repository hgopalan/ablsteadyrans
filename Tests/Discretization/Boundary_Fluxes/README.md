# Boundary_Fluxes

## Purpose
The boundary conditions on ERF's Neutral_ABL_Flat mesh (8 x 8 x 64, 1000 m) with a fixed log-law state
(u* = 0.35 m/s, z0 = 0.1 m, 15 degrees): the surface layer on the ground and a slip top. The test-only driver
`boundary_fluxes` (`Supporting_Files/Tools/boundary_fluxes.cpp`) evaluates the fluxes and the residual once.

## Files
| File | Role |
|---|---|
| `inputs_bc_neutral_flat` | ERF's Neutral_ABL_Flat mesh and sounding with the fixed log-law state. |
| `check_bc.py` | Runs the driver and checks the fluxes and the residual against their closed forms. |
| `Gold/bc_neutral_flat_plt00000/` | The driver's plotfile. |

## CTest
`rt_bc_neutral_flat`, `rt_bc_neutral_flat_matches` (double precision only).

## Expected Results
- u* recovered at every ground cell to 1e-10; the ground stress u*^2 (cos, sin) to 1e-10; k_w = u*^2 / Cmu0^2.
- The residual of every cell equals its closed form to 1e-6 of u*^2/dz for u and v (1e-13 measured), and to
  1e-9 K/s for theta. The first cell's truncation (-0.0971) and the top cell's lost flux (+1.000 of
  u*^2 cos/dz) are pinned. R_k and R_w are zero. Every field is finite.
- The plotfile matches `Gold/bc_neutral_flat_plt00000` at relative 1e-10.

## Gold files
| Gold | How it was made |
|---|---|
| `bc_neutral_flat_plt00000/` | `boundary_fluxes inputs_bc_neutral_flat` |

## How to run
```bash
ctest --test-dir <build> -R "^rt_bc_neutral_flat"
```
