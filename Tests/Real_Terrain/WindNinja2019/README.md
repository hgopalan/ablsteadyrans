# WindNinja2019

## Purpose
Beyond the plan: GUST on the three sites of Wagenbrenner, Forthofer, Page and Butler (2019), "Development and
evaluation of a Reynolds-averaged Navier–Stokes solver in WindNinja for operational wildland fire applications",
*Atmosphere* 10, 672. The transect plots of the paper (Figures 8, 12 and 15) are redrawn with GUST - Fast (mass
consistent) and GUST - CFD (steady RANS) next to the observations and the paper's WindNinja curves.

## Files
| File | Role |
|---|---|
| `inputs_askervein` | Askervein, 8.9 m/s from 210° at 10 m; 256 × 256 × 64 cells, 23.4 m, 1 m first layer; masts on lines A, AA, B at 10 m. |
| `inputs_bolund_270`, `inputs_bolund_239`, `inputs_bolund_90` | Bolund cases 1, 3, 4 (10.9 m/s from 270°, 8.7 m/s from 239°, 7.6 m/s from 90°, at 5 m); 256 × 128 × 48 cells, 1.6 m, 0.5 m first layer; line B (cases 1, 4) or A (case 3) at 5 m. |
| `inputs_bsb_wn_consistent` | Big Southern Butte, 8.3 m/s from 222° at 3 m, on the paper's 68 m cells and smoothed terrain (`bsb_wn.xyz`), 272 × 296 × 20 cells; an uncapped length and a no-slip lid moving at the inflow's top speed keep the inflow from drifting over the 18.7 km fetch; transects TSW, TWSW at 3 m. |
| `inputs_askervein_stable`, `inputs_bsb`, `inputs_bsb_wn`, `inputs_bsb_wn_flat`, `inputs_bsb_wn_flat_consistent`, `inputs_bsb_wn_precursor` | Diagnostics, not part of the comparison: Askervein with a lapse rate; the Butte on the SRTM terrain with 1 m first cells (does not reach the steady verdict); the paper-like Butte with a slip top (its inflow drifts by −18.5 %); the horizontal-homogeneity tests of the two Butte inflows on flat ground; a precursor column. |
| `windninja2019.py` | `prepare` (masts and terrain from windninja-tools), `smooth` (the paper-like Butte terrain), `sample` (speed-up at the masts and along the transects), `plot` (the paper's figures), `stats` (windward / ridge / lee RMSE, bias and MAPE). |
| `run_all.sh` | Runs the cases one at a time (GUST - Fast and GUST - CFD), samples and plots them. |
| `sweep.sh`, `sweep_variants.txt` | One-change variants of a case (a file of `name key=value ...` lines). |
| `wagenbrenner2019_digitized.csv` | The paper's Figures 8, 12 and 15, digitized (figure, panel, series, x, y). |
| `bsb_observed_fig15.csv` | The Butte observations of Figure 15 (speed-up and error-bar ends along TSW and TWSW). |
| `bsb_obs_1700.csv` | The 40 Butte sites at 17:00 local on 18 July 2010 (10-min mean speed and direction). |
| `bsb_profiles_fig14.csv`, `bsb_transects_fig14.csv` | The paper's Figure 14 (terrain along the transects), used to fit the smoothed Butte terrain. |

The comparison decks carry the recommended GUST - CFD setting for the lee: `turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, `turb.max_geom_lscale = 100` (the Butte: uncapped), `convection.hybrid_sigma_max = 0.1`
and `convection.hybrid_height` = 2.5 terrain heights (`docs/turbulence_model.rst`). The mass-consistent runs
drop the solver keys. All runs are truly neutral (the paper's CFD solver is neutral), use the log law with the
paper's roughness (0.01 m; Askervein 0.03 m), line LU-SGS with grid sequencing over 3 levels, and stop at the
steady verdict. Speed-up is (S − U_ref)/U_ref with U_ref the input wind, as in the paper.

## Data and acknowledgements
Downloaded by the user (not in the repository):
- **windninja-tools** (N. Wagenbrenner, <https://github.com/nwagenbrenner/windninja-tools>):
  `data/observed_data/points.asc` (the Askervein masts), `data/bolund/Bolund_Measurements.zip` and
  `data/bolund/bolund_hill_extended.tif` (the Bolund masts and 0.25 m terrain), `test_cases/big_butte.asc`
  (the Butte, SRTM 30 m). The repository states no licence, so none of its files is copied here;
  `windninja2019.py prepare` reads them from a checkout. `stats` for the Butte reads `big_butte.asc` from the
  environment variable `BSB_ASC`.

In the repository (CC BY 4.0 sources):
- Wagenbrenner, N.S., Forthofer, J.M., Page, W.G., Butler, B.W. (2019), *Atmosphere* 10, 672,
  doi:10.3390/atmos10110672: `wagenbrenner2019_digitized.csv`, `bsb_observed_fig15.csv`, `bsb_profiles_fig14.csv`
  and `bsb_transects_fig14.csv`, digitized from the figures. The observed points match the measurements to 0.002
  RMS in speed-up; where a curve is hidden under another in a figure it is traced along the covering curve.
  `plot` draws the paper's KELU (its recommended CFD setting) and COM curves; the LES curves are in the file but
  not plotted (Bolund case 3 has two LES curves, which the series mixes).
- Wagenbrenner, N.S., Big Southern Butte surface wind data (`bsb.sqlite`), Zenodo, doi:10.5281/zenodo.8226638:
  `bsb_obs_1700.csv`.

Field campaigns behind the data: Askervein (Taylor and Teunissen 1987, *Boundary-Layer Meteorol.* 39, 15),
Bolund (Berg et al. 2011, *Boundary-Layer Meteorol.* 141, 219; Bechmann et al. 2011, the blind comparison, DTU
Wind Energy), Big Southern Butte (Butler et al. 2015, *Atmos. Chem. Phys.* 15, 3785; Forthofer et al. 2014,
*Int. J. Wildland Fire* 23, 969; USDA Forest Service, Missoula Fire Sciences Laboratory).

We thank Natalie Wagenbrenner for windninja-tools, the Big Southern Butte database and the open-access paper; the
Missoula Fire Sciences Laboratory for the Big Southern Butte campaign; DTU Wind Energy and the Bolund team for the
Bolund data; and the Askervein Hill project for the Askervein data.

## How to run
```bash
Tests/Real_Terrain/WindNinja2019/run_all.sh build-rel ~/windninja-tools wn2019-runs 8
CASES="bolund_270" MODES="rans" Tests/Real_Terrain/WindNinja2019/run_all.sh build-rel ~/windninja-tools wn2019-runs
BSB_ASC=~/windninja-tools/test_cases/big_butte.asc python3 Tests/Real_Terrain/WindNinja2019/windninja2019.py \
    stats wn2019-runs/*_mc.csv wn2019-runs/*_rans_van_leer.csv --masts wn2019-runs/masts.csv
```
`run_all.sh` writes one directory per run, `<case>_<mode>.csv` samples and `<case>.png` plots in the run
directory. Runs are sequential (one case on the machine at a time).

## Expected Results
With the recommended setting, on 8 ranks of an Apple M2 Pro laptop (not idle; indicative):

| Case | GUST - CFD iterations (finest level) | Time | RMSE all masts, Fast / CFD (m/s) |
|---|---|---|---|
| Askervein | 150 | 3.5 min | 1.71 / 0.97 |
| Bolund case 1 | 365 | 2.9 min | 2.27 / 1.25 |
| Bolund case 3 | 361 | 2.9 min | 1.47 / 0.64 |
| Bolund case 4 | 815 | 6.0 min | 1.01 / 1.23 |
| Big Southern Butte | 158 | 1.3 min | 3.34 / 3.22 |

The mass-consistent runs take 2–5 s. Known differences: the Butte lee is faster than measured with both modes,
and the Bolund case 4 lee mast is slower than measured with the cap (bias −2.4 m/s).

## CTest
None (the runs take minutes; the unit tests cover the new options).
