# Input_Validation

## Purpose
Decks with mistakes in them. Each must stop the application, and the message must name the input. Every
message listed in `Tests/CMakeLists.txt` must appear, so one deck checks several mistakes at once.

## Files
| File | Mistakes |
|---|---|
| `inputs_abort_unknown_key` | A misspelt key, an unknown namespace, keys that only apply with `steady.stage = simplec`, a Poisson key the application does not read. |
| `inputs_abort_range` | `mesh.n_cell` below 1, a periodic z, non-finite `mesh.prob_lo` / `mesh.prob_hi`. |
| `inputs_abort_no_stage` | Solver keys without `steady.stage = simplec`. |
| `inputs_abort_stage_no_guess` | `steady.stage = simplec` without a first guess. |
| `inputs_abort_repeated_key` | A key given twice; a key given more than one value. |

Check script: `Supporting_Files/expect_abort.py` (shared).

## CTest
`abort_unknown_key`, `abort_range`, `abort_no_stage`, `abort_stage_no_guess`, `abort_repeated_key`.

## How to run
```bash
ctest --test-dir <build> -R "^abort_(unknown_key|range|no_stage|stage_no_guess|repeated_key)$"
```
