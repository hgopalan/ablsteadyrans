# Askervein_210

## Purpose
RT-2.6.1 and RT-2.6.2 (sub-phase 2.6, the Phase 2 checkpoint): Askervein with ERF's own setup, so that ERF's URANS
is a reference for the same case. ERF's deck is `Exec/CanonicalTests/Real_Terrain/Askervein/inputs_anel` at the
pinned commit: the power-law inflow from 210° through the west and south faces, outflow through the east and north,
a Rayleigh layer on w in the top 250 m, neutral, no Coriolis, on a 300 x 300 x 18 mesh (20 m, a 20 m first layer
stretched by 1.08).

## Files
| File | Role |
|---|---|
| `inputs_steady_askervein_210` | The deck. Terrain and inflow on the command line: `terrain.file=extern_erf/askervein.txt`, `bc.inflow_file=Tests/Real_Terrain/Askervein_210/askervein_210_inflow.txt`. The deck sets no sweep, CFL cap or sequencing; RT-2.6.1 adds `steady.sweep=lusgs steady.sequencing=1 mesh.hierarchy_levels=2`, the old-vs-new run the Phase 1 settings. |
| `askervein_210_inflow.txt` | z u v θ rows from `extern_erf/inflow_powerlaw.txt` (`tools/erf_inflow_table.py --theta 300`). |
| `askervein_lines.csv` | The masts of lines A and AA (Taylor and Teunissen 1987): positions in the frame of the terrain file, as tabulated in N. Wagenbrenner's windninja-tools. |
| `Gold/erf_urans_lines_10m.csv` | ERF's 10 m speeds at the masts, averaged over 2400–4800 s. |
| `Gold/askervein_210_level_history.json` | The grid-sequencing history of the Phase 2 run (RT-2.6.2). |

## The ERF reference
ERF at the pinned commit ran `inputs_anel` with the overrides of `tools/kestrel/erf_references.sh`
(`stop_time=4800 erf.plot_per_1=100 erf.plot_int_1=-1 erf.plot2d_int_1=-1`) on 8 ranks of one workstation
(10 378 steps, 8.6 h), was archived with `ablsteady_refavg` (25 snapshots, 2400–4800 s; the half-window averages
differ by at most 7.6e-4 m/s) and sampled with `tools/askervein_lines.py sample`. The archive itself is not in the
repository (100 MB).

## CTest
`rt_askervein_210_inflow_table` (`regression;smoke`); `rt_steady_askervein_210`, `rt_askervein_210_erf_lines`,
`rt_askervein_210_level_history`, `rt_steady_askervein_210_phase1`, `rt_askervein_210_phase1_matches` (`nightly`,
double precision only).

## Expected Results
- Phase 2 configuration (`lusgs` and its default cap 8, the default reference window, two levels): steady after 96
  level-0 iterations (34 on level 1; 96 and 42 with `bc.inflow_tke_length = unbounded`). Phase 1 configuration (`zjacobi`, four sweeps, cap 2, a 100-iteration
  reference): 512 (538 with `unbounded`). The two agree on the lines to 1e-4 %.
- Against ERF: within 4.00 % up to 300 m past the hilltop; the lee masts beyond are 6.6–9.2 % slower (second-order
  against ERF's fifth-order advection, and less inflow turbulence mixing into the wake than with the unbounded
  inflow length, which gave 4.48 % and 4.9–7.5 %). See `docs/results.rst`, "Phase 2 checkpoint".

## How to run
```bash
ctest --test-dir <build> -R "askervein_210"
```
