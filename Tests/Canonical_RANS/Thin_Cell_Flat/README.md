# Thin_Cell_Flat

## Purpose
The w rows' VNN cap (`simplec.w_vnn_cap`, `docs/solution_algorithm.rst`) on a flat inflow/outflow strip with a thin first cell.
16 x 2 x 40 cells over 1500 m x 187.5 m x 742 m: dx = 94 m, a 0.7 m first cell (dx/dz0 = 134) stretched by 1.2;
y periodic, the Askervein power-law inflow on xlo, outflow on xhi, the surface layer (z0 = 0.03 m), a slip top,
Rayleigh damping of w above 250 m, no gravity. Line LU-SGS with two sweeps at a CFL cap of 16; the verdict at 6
orders below the reference and 5e-7.

Without the cap the w residual decays slowly on thin first cells: the SIMPLEC denominator of the w rows is about
1/dtau, so the pressure correction moves w much further than the implicit predictor does where the vertical
diffusion is stiff. The cap bounds that ratio (the VNN number) and leaves the steady state unchanged.

## Files
| File | Role |
|---|---|
| `inputs_thin_cell_flat` | The strip (the cap, the iterations, the first guess and the restart on the command line). |
| `check_w_vnn_cap.py` | Runs the deck with the cap off and on, and a capped restart chain; checks the results. |

## CTest
`rt_thin_cell_w_vnn_cap` (double precision only; label `regression`).

## Expected Results
- Both runs reach the verdict; the stage banner says `w_vnn_cap off` and `w_vnn_cap 16 (slope margin 4)`.
- The capped run needs at most a fifth of the iterations. Measured: 7038 without, 322 with (22 times fewer;
  7038 and 299 before the residual scales used the depth of the layers, 5.1 km here, instead of `mesh.prob_hi` z;
  6487 and 687 before the Rayleigh layer ended at the mesh's top level).
- u, v and w of the two steady states within 2e-4 m/s. Measured: 8.5e-9 (u), 5.8e-9 (w), 3e-15 (v) (before the
  Rayleigh change: 5.6e-5 and 3.8e-5).
- A capped run of 40 iterations continued to 80 from its checkpoint equals 80 straight in every field, bit for bit.
- About 50 s on two ranks.

## How to run
```bash
ctest --test-dir <build> -R rt_thin_cell_w_vnn_cap -V
```
