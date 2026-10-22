# RUSHIL

## Purpose
The RUSHIL two-dimensional hills (Khurshudyan, Snyder and Nekrasov 1981) at the full scale of
Manning, Hancock and Whiting (2010): H = 117 m, a/H = 3, 5 and 8. Hill 3 separates in the tunnel, Hill 8 stays
attached.

## Files
| File | Role |
|---|---|
| `inputs_rushil` | One deck for the three hills (`terrain.file = hill3.xyz`, `hill5.xyz` or `hill8.xyz`, summit at x = 4680 m): 936 × 8 × 64 cells over ±40 H, 10 m horizontally, a 1 m first layer, slip ceiling at 1603 m, z0 = 0.157 m, the tunnel's flat-floor profile as inflow. |
| `rushil.py` | `prepare` (the three terrain files, the inflow table and the station files from the ERCOFTAC data), `compare` (reversed-flow extent at 0.1 H, U rms against the tunnel), `plot` (speed-up planes and station profiles; with `--manning`, Meteodyn WT curves from a user-supplied digitized file). |

The deck carries the recommended GUST - CFD setting (`../WindNinja2019/README.md`) with
`convection.hybrid_height = 292` m. `turb.max_geom_lscale = 100` m: with the default 30 m cap the near-ground
speed-up of the tunnel profile falls 12 % over the 40 H fetch (flat-ground test); 100 m holds it within 4 % and
keeps the Hill 3 bubble that an uncapped length shortens.

## Data and acknowledgements
Downloaded by the user (not in the repository): the files of ERCOFTAC Classic Collection case 69 / Application
Challenge AC5-05, "Flow over 2D hills (RUSHIL)" (`AC5-05_Flatrefcase.dat`, `AC5-05_h3??s.dat`, `AC5-05_h5??s.dat`,
`AC5-05_h8??s.dat`), from the ERCOFTAC Knowledge Base. The measurements are those of Khurshudyan, L.H., Snyder,
W.H., Nekrasov, I.V. (1981), "Flow and dispersion of pollutants over two-dimensional hills", EPA-600/4-81-067, made
in the RUSHIL tunnel of the US EPA Fluid Modeling Facility.

Not distributed: the Meteodyn WT and WAsP curves of Manning, J., Hancock, P., Whiting, P. (2010), "Wind flow over
two-dimensional hills: a comparison of Meteodyn WT and WAsP with wind tunnel data", *Wind Engineering* 34, 477.
`rushil.py plot --manning CSV` draws them from a file the user digitizes (figure, panel, series, x, y; Figs.
4.3–4.6).

We thank ERCOFTAC and the QNET-CFD Knowledge Base for keeping the RUSHIL data available, and the US EPA Fluid
Modeling Facility for the original measurements.

## How to run
```bash
python3 Tests/Real_Terrain/RUSHIL/rushil.py prepare ~/Downloads/rushil rushil
cd rushil && mpirun -np 8 ../build-rel/gust ../Tests/Real_Terrain/RUSHIL/inputs_rushil \
    terrain.file=hill3.xyz bc.inflow_file=inflow.txt
python3 ../Tests/Real_Terrain/RUSHIL/rushil.py compare plt --hill 3 --outdir .
python3 ../Tests/Real_Terrain/RUSHIL/rushil.py plot plt --hill 3 --outdir . --out hill_3.png
```

## Expected Results
With the recommended setting (8 ranks, Apple M2 Pro):

| Hill | Reversed flow at 0.1 H | Tunnel | U rms vs tunnel | Iterations, time |
|---|---|---|---|---|
| 3 | 0.58a to 1.45a | 0.27a to 2.17a (Meteodyn WT 0.8a–1.8a) | 0.48 m/s | 298, 40 s |
| 5 | none | — | 0.21 m/s | 241, 32 s |
| 8 | none | attached | 0.20 m/s | 229, 32 s |

With the default closure Hill 3 separates from 0.64a to 1.27a and Hill 5 does not reach the steady verdict. The
reversed flow on Hill 3 reaches about 0.2 H above the ground, against about 0.4 H in the tunnel.

## CTest
None.
