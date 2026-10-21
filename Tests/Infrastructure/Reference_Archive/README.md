# Reference_Archive

## Purpose
`gust_refavg` end to end. Three copies of the smoke plotfile at t = 0, 5 and 10 s are archived on two ranks
with two inputs files. The archive must verify on three ranks, and its trapezoid average (weights 1/4, 1/2, 1/4)
must reproduce the smoke plotfile bit for bit. A damaged archive and bad tool inputs must be refused.

This checks the archive round trip. Identical snapshots cannot tell weights apart; the weights are tested in the
unit tests (UT-1.1.2).

## Files
| File | Role |
|---|---|
| `refavg_fixtures.py` | `snapshots`: writes the three time-stamped copies. `damage`: copies an archive and flips one bit. |

The archived inputs files are `Infrastructure/Smoke_Restart/inputs_smoke` and `inputs_restart`.

## CTest
`refavg_snapshots`, `refavg_build`, `refavg_verify`, `refavg_matches`, `refavg_damage`, `abort_refavg_damaged`,
`abort_refavg_inputs`, `abort_refavg_verify_keys`. They run in `<build>/Tests/run/smoke` and need `smoke_run`.

## Expected Results
- The archive verifies with case `Smoke`, ERF 92410d578, 3 snapshots, `kestrel-cpu`, 2 nodes, 8 ranks, 10.000 AU
  and a zero half-window change.
- The average equals `smoke_plt00000` at zero tolerance.
- The damaged copy is refused, naming the file whose SHA-256 differs.
- Bad `refavg.*` inputs stop the tool with each input named.

## How to run
```bash
ctest --test-dir <build> -R "refavg"
```
