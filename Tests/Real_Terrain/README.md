# Real_Terrain

## Purpose
Askervein hill run to steady on ERF's grid (`extern_erf/askervein.txt`). The folder name follows ERF's
`Exec/CanonicalTests/Real_Terrain/Askervein`.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Askervein` | `rt_steady_askervein` | Log-law inflow, outflow and slip sides, no Coriolis, to the steady verdict. |
| `Askervein_210` | `rt_askervein_210_erf_lines`, `rt_askervein_210_level_history` | ERF's own setup (210° power-law inflow, ERF's full mesh) against ERF's URANS on lines A and AA. |
| `Askervein_CNBL` | `rt_precursor_askervein_cnbl`, `rt_steady_askervein_cnbl` | A conventionally neutral boundary layer: a precursor column writes the inflow table, then the 3-D case runs on it. |
| `WindNinja2019` | (comparison case) | Askervein, Bolund cases 1, 3, 4 and Big Southern Butte against the observations and WindNinja curves of Wagenbrenner et al. (2019); GUST - Fast and GUST - CFD. |
| `IshiharaZhou` | (comparison case) | The steep 2-D ridge of Ishihara and Zhou (2025), smooth and rough, against the wind tunnel. |
| `RUSHIL` | (comparison case) | The RUSHIL 2-D hills (a/H = 3, 5, 8) against the wind tunnel (ERCOFTAC AC5-05). |

The first three use `Supporting_Files/expect_abort.py --success "gust: steady after"` to require the verdict.
