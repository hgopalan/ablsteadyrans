# Initialization

## Purpose
Cases that check the mass-consistent initial field: a first guess projected onto a divergence-free field with the
horizontal and vertical weights alpha_h and alpha_v.

## Cases
| Folder | Plan ID | Description |
|---|---|---|
| `Mass_Consistent_Askervein` | RT-1.12.1 | The initial field on the Askervein mesh: iterations, divergence, field extrema and the state. |
| `Bump_Speedup_Trend` | RT-1.12.2 | The hilltop speed-up over FastWindTerrain's bump for five values of alpha_v, against FastWindTerrain's own. |

Both cases use `Supporting_Files/check_init.py` (modes `askervein` and `trend`).
