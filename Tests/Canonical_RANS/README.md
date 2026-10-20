# Canonical_RANS

## Purpose
The SIMPLEC stage run to the steady verdict on the canonical RANS cases. The folder names follow ERF's
`Exec/CanonicalTests/Canonical_RANS/`, and each case runs ERF's own check script (`extern_erf/`) on the
turbulence plotfile where ERF has one.

## Cases
| Folder | Plan ID | Description |
|---|---|---|
| `Neutral_Column` | RT-1.15.1 | The 1-D neutral column driven by a pressure-gradient force; a restart chain that must continue bit for bit. |
| `Neutral_ABL_Flat` | RT-1.15.2 | ERF's Neutral_ABL_Flat to steady, with ERF's `check_neutral.py --physics`. |
| `Neutral_Hill_2D` | RT-1.16.1, RT-2.1.1, RT-2.1b.1 | ERF's Neutral_Hill_2D with Jacobi sweeps, line LU-SGS sweeps and the column first guess. |
| `Neutral_Hill_3D` | RT-1.16.2 | ERF's Neutral_Hill_3D (64 x 64 x 20); label `nightly`. |
| `Truly_Neutral_Hill_3D` | RT-2.2.1 | A truly neutral Hill_3D to round-off with z-line Jacobi and box-local line LU-SGS sweeps; the two states agree to 1e-8. Label `nightly`. |
| `Thin_Cell_Flat` | (beyond the plan) | A flat inflow/outflow strip with a 0.7 m first cell, with and without `simplec.w_vnn_cap`: fewer iterations, the same steady state, a restart bit for bit. |
| `Steep_Ridge_Thin_Cell` | (beyond the plan) | A steep 2-D ridge (32.5°) with a 0.4 m first cell: `simplec.w_vnn_cap = 16` with the default slope margin is steady within 3500 iterations, without the margin it is not. |
| `Smooth_Plateau_2D` | (comparison case) | A 2-D plateau run in both codes; not a CTest. |

The hill cases use `Supporting_Files/check_steady_hill.py` (shared by the 2-D and 3-D folders).
