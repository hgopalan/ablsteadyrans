# Real_Terrain

## Purpose
Askervein hill run to steady on ERF's grid (`extern_erf/askervein.txt`). The folder name follows ERF's
`Exec/CanonicalTests/Real_Terrain/Askervein`.

## Cases
| Folder | Plan ID | Description |
|---|---|---|
| `Askervein` | RT-1.16.3 (convergence half) | Log-law inflow, outflow and slip sides, no Coriolis, to the steady verdict. |
| `Askervein_210` | RT-2.6.1, RT-2.6.2 | ERF's own setup (210° power-law inflow, ERF's full mesh) against ERF's URANS on lines A and AA. |
| `Askervein_CNBL` | beyond the plan | A conventionally neutral boundary layer: a precursor column writes the inflow table, then the 3-D case runs on it. |

Both use `Supporting_Files/expect_abort.py --success "gust: steady after"` to require the verdict.
