# Askervein_CNBL

## Purpose
Askervein in a conventionally neutral boundary layer. A precursor column runs to steady and
supplies the inflow table; the 3-D case then runs on that table to the steady verdict.

## Files
| File | Role |
|---|---|
| `inputs_precursor_askervein_cnbl` | The precursor column: 4 x 4 x 26 on Askervein's 20 m x 1.08 levels, z0 = 0.03 m, the flat case's sounding, Coriolis at 57.2 N, geostrophic wind (17.028798, -7.380863) m/s. |
| `askervein_cnbl_inflow.txt` | The committed inflow table (z u v theta k), written from the precursor by `tools/make_inflow_table.py`. |
| `inputs_steady_askervein_cnbl` | The 3-D case: 30 x 30 x 26 to 1599 m, the table as inflow on xlo and yhi, outflow on xhi and ylo, Rayleigh damping of w and theta above 1200 m, four sweeps, CFL cap 2. The table and terrain file are given on the command line. |

The table check is `tools/make_inflow_table.py --check`; the run verdicts use `Supporting_Files/expect_abort.py`
(shared). At configure time, `Tests/CMakeLists.txt` writes `<build>/Tests/inputs_precursor_askervein_cnbl_fixed`:
the precursor deck with `steady.max_iterations` set to the step in the table's header and an unreachable
`steady.residual_drop_orders`.

## CTest
`rt_precursor_askervein_cnbl`, `rt_precursor_askervein_cnbl_fixed`, `rt_precursor_askervein_cnbl_table`,
`rt_steady_askervein_cnbl` (label
`regression`, double precision only).

## Expected Results
- Precursor: the steady verdict (11087 iterations with Apple clang, 11073 with gcc-16).
- Fixed-step run: stops at the table's step (11087) without the verdict.
- Table check, on the fixed-step run: the plotfile step equal to the table's, the 10 m wind within 0.01 degrees
  of +x, u* within 1e-3 of 0.6 m/s, and `askervein_cnbl_inflow.txt` within 1e-3 of the planar means.
- 3-D case: the steady verdict. Measured: 220 iterations; u residual down four orders.

Why a fixed step: under the capping inversion theta has no fixed point. The verdict comes when the k residual leaves
a limit cycle of about 600 iterations, at an iteration that depends on the compiler's rounding (11087 with Apple
clang, 11073 with gcc-16, 1.7e-3 K apart in theta). At step 11087 the two builds agree to 2e-4 K.

## How to run
```bash
ctest --test-dir <build> -R "askervein_cnbl"
```
