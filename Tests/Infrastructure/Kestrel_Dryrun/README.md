# Kestrel_Dryrun

## Purpose
The Kestrel scripts (`tools/kestrel/`) cannot run here. This case checks that they parse, and runs them with stub
`module`, `sbatch` and `srun` commands to exercise the submit and archive paths.

## Files
| File | Role |
|---|---|
| `kestrel_dryrun.sh` | `kestrel_dryrun.sh SRC REFAVG PLOTFILE WORK`: builds a fake ERF tree and finished run in `WORK`, then submits and archives it. |

## CTest
`kestrel_parse_*` (one per script, `bash -n`), `kestrel_dryrun` (needs `smoke_run` for its plotfile).

## Expected Results
- The Askervein job script has the expected run length, plot period, rank count and no 2-D plotfiles.
- The archive averages only the 3-D plotfiles in the final window and keeps the inputs but not the run's outputs.
- A second case with the default prefix `plt` keeps its `plt2d_` plotfiles out.

## How to run
```bash
ctest --test-dir <build> -R "^kestrel_"
```
