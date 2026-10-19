# Smoke_Restart

## Purpose
A short run of the application on two ranks, then a restart from its checkpoint on three ranks and another box
layout. The restart must continue the run bit for bit. Two further decks restart from the same checkpoint with a
mistake and must stop.

## Files
| File | Role |
|---|---|
| `inputs_smoke` | The 2-rank run; writes `smoke_plt00000` and `smoke_chk00000`. |
| `inputs_restart` | Restarts from `smoke_chk00000`; writes `restart_plt00000`. |
| `inputs_abort_restart_mesh` | Restart with a different `mesh.n_cell`: must stop. |
| `inputs_abort_restart_overwrite` | `io.chk_file` that would overwrite the restart checkpoint: must stop. |

Check script: `Supporting_Files/expect_abort.py` (shared), which runs a command and requires the given messages.

## CTest
`smoke_clean`, `smoke_run`, `smoke_restart`, `smoke_restart_matches`, `abort_restart_mesh`,
`abort_restart_overwrite`. All share the run directory `<build>/Tests/run/smoke`, which the reference-archive
and Kestrel cases also read.

## Expected Results
- The restart log shows `restarted from smoke_chk00000 at step 0` (the smoke state is mostly zeros, so the
  plotfile comparison alone cannot tell a restart from a fresh start).
- `smoke_plt00000` and `restart_plt00000` are identical under `amrex_fcompare` at zero tolerance
  (`--allow_diff_grids`, since the layouts differ).
- The two abort decks stop with the input named.

## How to run
```bash
ctest --test-dir <build> -R "^smoke_|^abort_restart_"
```
