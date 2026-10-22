# Column_Target

## Purpose
The target wind of the column first guess (`init.column_target_u`, `init.column_target_height`) carried into the
3-D run. The 1-D column of `Neutral_Column` (4 x 4 x 64, 1000 m, z0 = 0.1 m, slip top, pressure-gradient forcing)
starts from the converged column with a target of 8 m/s at 10 m. The column adjusts `forcing.pressure_gradient`
until its wind at 10 m is the target; the 3-D run takes that forcing. The Rayleigh layer above 700 m acts on w and
theta only: a layer on u and v would hold the 3-D wind to the deck's inflow profile, and is refused with a target.

The second run adds a stable column, `init.column_mol = 1000`, with `bc.surface_stability = heat_flux`: the 3-D run
takes the column's surface flux and starts from its theta.

## Files
| File | Role |
|---|---|
| `inputs_column_target` | The neutral run; the stable run adds its keys on the command line. |
| `check_column_target.py` | Runs the application and checks the 3-D wind at the target height and the surface's L. |

## How to run
```bash
ctest --test-dir <build> -R "^rt_column_target"
```

## Expected Results
- The column reaches the target, and the log lists the inputs a restart deck needs (`forcing.pressure_gradient`;
  stable: also `bc.surface_heat_flux`).
- The 3-D wind at 10 m (planar mean, linear in ln(z + z0) between centres) is 8 m/s to 0.5 %, along x.
  Measured after 1500 iterations: 8.004 m/s neutral, 8.008 m/s stable.
- Stable: the 3-D surface's Obukhov length, from its mean u* and the flux it took, is 1000 m to 2 %.
  Measured: 1001 m.
- The runs do not reach the steady verdict in 1500 iterations: they start from the converged column, so the
  residuals' reference is already small and six orders of drop are out of reach. The wind is steady to the
  tolerance checked.

## CTest
`rt_column_target`, `rt_column_target_mol`.
