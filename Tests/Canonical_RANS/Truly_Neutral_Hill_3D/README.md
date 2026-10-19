# Truly_Neutral_Hill_3D

## Purpose
RT-2.2.1. A truly neutral variant of ERF's Neutral_Hill_3D, run to round-off with two different sweeps, whose
steady states must agree to 1e-8. The conventionally neutral deck (`Neutral_Hill_3D`) keeps entraining the
inversion and has no fixed point to that tolerance, hence the variant: theta = 300 K, `forcing.gravity = 0`,
Rayleigh damping of w alone, 32 x 32 x 16 over the same domain, four boxes.

## Files
| File | Role |
|---|---|
| `inputs_steady_hill3d_neutral` | 16 z-line Jacobi sweeps, CFL cap 4. |
| `inputs_steady_hill3d_neutral_lusgs` | The box-local line LU-SGS sweeps: zero ghosts, overlap 2, the diagonal ordering, two sweeps, CFL cap 16. |
| `check_line_match.py` | Compares the two final plotfiles: u, v, w, p and k to `--rel_tol` of each field's maximum; theta after subtracting each run's mean. |

The two runs require the steady verdict through `Supporting_Files/expect_abort.py --success` (shared).

## CTest
`rt_steady_hill3d_neutral_clean`, `rt_steady_hill3d_neutral_zjacobi`, `rt_steady_hill3d_neutral_lusgs` (four
ranks each), `rt_steady_hill3d_neutral_matches`. Label `nightly`, double precision only. The clean step empties
the run directories first, so a stale plotfile cannot be compared, and the checker requires exactly one plotfile
per run.

## Expected Results
- Both runs reach round-off: every residual below 1e-11 and ten orders down.
- Measured: 3495 iterations (70 s on four ranks) for z-line Jacobi against 1792 (32 s) for LU-SGS.
- u, v, w, p and k agree to 9e-10 absolute (at most 3e-8 of the tolerance scale).
- theta has a free constant (no theta condition, no gravity) that drifts by 1e-4 K at round-off. After subtracting
  each run's mean the two agree to 3e-11 K.

## How to run
```bash
ctest --test-dir <build> -R "^rt_steady_hill3d_neutral"
```
