# Steep_Ridge_Thin_Cell

## Purpose
The slope margin of the w rows' VNN cap (`simplec.w_vnn_cap_slope_margin`, `docs/solution_algorithm.rst`) on a steep
2-D ridge with a thin first cell. 48 x 2 x 40 cells over 1500 m x 187.5 m (dx = 31 m, dy = 94 m), a 0.4 m first layer
stretched by 1.2 (columns about 2.9 km high; `mesh.prob_hi` z = 742 m sets only the computational spacing), an
Agnesi ridge 150 m high with a 150 m half-width (max ground slope 32.5 deg; max aspect ratio max(dx, dy)/thickness
372, from dy and the 0.25 m first layer at the crest); y periodic, the Askervein power-law inflow on xlo, outflow on xhi, the surface layer (z0 = 0.03 m), a slip
top, Rayleigh damping of w above 250 m, no gravity. Line LU-SGS with two sweeps at a CFL cap of 16; the verdict at
3 orders below the reference and the default 5e-5.

With the cap on, the pressure correction's z-face coefficient is much smaller than the horizontal one near the
ground; on steep slopes the pressure operator then turns indefinite and the iteration stalls. The margin lowers the
cap on sloping ground just enough to keep the operator definite.

## Files
| File | Role |
|---|---|
| `inputs_steep_ridge_thin_cell` | The ridge from the profile first guess (the cap, the margin and the iterations on the command line). |

## CTest
- `rt_steep_ridge_w_vnn_cap`: V = 16 with the default margin (4) is steady within 3500 iterations.
- `rt_steep_ridge_w_vnn_cap_margin_off`: V = 16 with the margin off (#32's cap) is not steady after 3500. The test
  checks only that verdict; it does not tell a stall from a slow run (the run without the cap is not steady at 3500
  either).

Both have the label `regression`.

## Expected Results
- With the margin: steady in 1765 iterations (the tests stop at 3500); 1652 before the residual scales used the
  depth of the layers, 1895 before the inflow k took the closure's length (#44). Measured before #44 (not
  re-measured): margin 2: 2035; 8: 1940; cap off: 2881.
  Before the Rayleigh layer ended at the mesh's top level (the layers end at 2.9 km, above `mesh.prob_hi` z): 2568
  (margin 2: 2714; 8: 3111; cap off: 3983).
- Without it: the residual ratio stays between 0.1 and 0.7 after the transient and the CFL keeps being cut.
- Measured before the Rayleigh change (not re-measured): the first 2300 iterations are the same transient in every configuration (the lee separation forms while the
  residual rises to about 11 times its reference and the SER law lowers the CFL to about 0.08).
- To the 6-order verdict (not run here; before the Rayleigh change): 14154 iterations with the cap off, 4142 with V = 16 and margin 2, 3169
  with V = 64 without a margin; the three steady states agree to 3e-8 m/s. Margin 1 (the bound itself) stalls.
- About 55 s per run on two ranks.

## How to run
```bash
ctest --test-dir <build> -R rt_steep_ridge -V
```
