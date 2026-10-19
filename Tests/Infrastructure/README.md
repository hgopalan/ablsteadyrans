# Infrastructure

## Purpose
Cases that check the application around the solver: a run and its restart, input decks that must stop, reference
archives and the Kestrel scripts.

## Cases
| Folder | Description |
|---|---|
| `Smoke_Restart` | A 2-rank run, then a 3-rank restart on another box layout; the two plotfiles identical. Two decks that must stop on a bad restart. |
| `Input_Validation` | Decks that must stop, with the input named in the message. |
| `Reference_Archive` | `ablsteady_refavg` end to end on copies of the smoke plotfile. |
| `Kestrel_Dryrun` | A dry run of the Kestrel submit and archive scripts with stub `module`/`sbatch`/`srun`. |
