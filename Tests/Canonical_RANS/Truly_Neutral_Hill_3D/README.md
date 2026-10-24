# Truly_Neutral_Hill_3D

## Purpose
A truly neutral variant of ERF's Neutral_Hill_3D, run to round-off with two different sweeps, whose
steady states must agree to 1e-8. The conventionally neutral deck (`Neutral_Hill_3D`) keeps entraining the
inversion and has no fixed point to that tolerance, hence the variant: theta = 300 K, `forcing.gravity = 0`,
Rayleigh damping of w alone, 32 x 32 x 16 over the same domain, four boxes.

## Files
| File | Role |
|---|---|
| `inputs_steady_hill3d_neutral` | 16 z-line Jacobi sweeps, CFL cap 4. |
| `inputs_steady_hill3d_neutral_lusgs` | The box-local line LU-SGS sweeps: zero ghosts, overlap 2, the diagonal ordering, two sweeps, CFL cap 16. |
| `inputs_sequencing_z0map`, `z0_two_points.xyz` | The roughness map (`bc.surface_z0_file`) on two sequencing levels: two points (0.0003 m at x = −1000 m, 0.015 m at x = 1000 m), two IDW neighbours. |
| `check_line_match.py` | Compares the two final plotfiles: u, v, w, p and k to `--rel_tol` of each field's maximum; theta after subtracting each run's mean. |

The two runs require the steady verdict through `Supporting_Files/expect_abort.py --success` (shared).

It also carries the old-against-new test of the default `solver.poisson_backend`: the LU-SGS deck runs a second time
with `solver.poisson_backend = gmres_fft` (the default before `gmres_mlmg`), and `check_line_match.py` compares it
with `rt_steady_hill3d_neutral_lusgs` (the default, `gmres_mlmg`, checked in its log) at 1e-8.

`rt_steady_hill3d_neutral_sequencing` runs the z-line Jacobi deck a third time with grid sequencing (`mesh.hierarchy_levels=3
steady.sequencing=1` on the command line: 8 x 8 and 16 x 16 to the loose stop, then 32 x 32 to the strict stop),
and `check_line_match.py` compares it with `rt_steady_hill3d_neutral_zjacobi` at 1e-8. The checker counts only
`PREFIX<step>`, not the per-level plotfiles `PREFIX_level<l>_<step>` a hierarchy run also writes.

## CTest
`rt_steady_hill3d_neutral_clean`, `rt_steady_hill3d_neutral_zjacobi`, `rt_steady_hill3d_neutral_lusgs`,
`rt_poisson_default_gmres_fft`, `rt_steady_hill3d_neutral_sequencing` (four ranks each), `rt_steady_hill3d_neutral_matches`,
`rt_poisson_default_matches`, `rt_steady_hill3d_neutral_sequencing_matches`. Label `nightly`, double precision only. The clean step empties
the run directories first, so a stale plotfile cannot be compared, and the checker requires exactly one plotfile
per run.

`rt_sequencing_z0map` (two ranks, three iterations per level, label `regression;smoke`) pins the logged summary
of each level's map, computed independently in Python from the two points (ln z0 inverse-distance weights at the
32 x 32 cell centres, then the geometric mean of each 2 x 2 block): level 0 min 0.0003023, max 0.01489; level 1
min 0.0003134, max 0.01436; both geometric means 0.002121 m. Level 1 sampled by IDW at its own centres instead
would give min 0.0003105, max 0.01449; the default terrain.idw_neighbors = 6 would abort (two points).

## Expected Results
- Both runs reach near round-off: every residual below 1e-11 and eight orders down (ten orders put theta and p under the rounding floor; see the deck).
- Measured: 3495 iterations (70 s on four ranks) for z-line Jacobi against 1792 (32 s) for LU-SGS.
- u, v, w, p and k agree to 9e-10 absolute (at most 3e-8 of the tolerance scale).
- `gmres_mlmg` against `gmres_fft` on the LU-SGS deck: both stop at 1789 iterations; u, v, w, p and k agree to 2.2e-13 absolute or better, theta minus its mean to 2.4e-11 K.
- theta has a free constant (no theta condition, no gravity) that drifts by 1e-4 K at round-off. After subtracting
  each run's mean the two agree to 3e-11 K.

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_hill3d_neutral"
```
