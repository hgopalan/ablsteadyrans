# IshiharaZhou

## Purpose
Beyond the plan: the steep two-dimensional ridge of Ishihara and Zhou (2025) at full scale (tunnel × 1000), a
test of lee separation and reattachment, on a smooth (z0 = 0.01 m) and a rough (z0 = 0.3 m) surface. The setup
follows H. Gopalan's Kynema-SGF cases (kynema/kynema-sgf#2031, #2036).

## Files
| File | Role |
|---|---|
| `inputs_iz_smooth` | Smooth surface: z = 40 cos²(π (x − 800)/200) m for \|x − 800\| ≤ 100 m (maximum slope about 32°) in 2400 m × 896 m under a slip ceiling, periodic in y; 600 × 8 × 64 cells, 4 m, 0.5 m first layer; the tunnel's approach flow as inflow (u* 0.21 m/s). Terrain and inflow on the command line. |
| `inputs_iz_rough` | The same with z0 = 0.3 m (u* 0.32 m/s). |
| `ishihara_zhou.py` | `prepare` (terrain, inflow tables and station files from the article's data), `compare` (reversed-flow extent on the near-surface line, U rms against the tunnel stations), `plot` (near-surface line and station profiles). |

Both decks carry the recommended GUST - CFD setting (`../WindNinja2019/README.md`) with
`convection.hybrid_height = 100` m (2.5 ridge heights).

## Data and acknowledgements
Downloaded by the user (not in the repository): the supplementary spreadsheets mmc1.xls–mmc4.xls (smooth and
rough approach flow, smooth and rough ridge stations) of Ishihara, T., Zhou, Q. (2025), "Wind tunnel data of
turbulent flows over a two-dimensional steep ridge", *Data in Brief* 63, 112260 (CC BY 4.0), from the article page.
`ishihara_zhou.py prepare DATA OUTDIR` writes `ridge.xyz`, `inflow_smooth.txt`, `inflow_rough.txt` and
`stations_*.csv`.

The LES reattachment lengths quoted by `compare` are those of Liu, Z., Diao, Z., Ishihara, T. (2019), *Renewable
Energy* 136, 968, as digitized in kynema/kynema-sgf#2036.

We thank Takeshi Ishihara and Qiang Zhou for publishing the wind-tunnel data openly.

## How to run
```bash
python3 Tests/Real_Terrain/IshiharaZhou/ishihara_zhou.py prepare ~/Downloads/ishihara_zhou iz
cd iz && mpirun -np 8 ../build-rel/gust ../Tests/Real_Terrain/IshiharaZhou/inputs_iz_smooth \
    terrain.file=ridge.xyz bc.inflow_file=inflow_smooth.txt
python3 ../Tests/Real_Terrain/IshiharaZhou/ishihara_zhou.py compare plt --surface smooth --outdir .
python3 ../Tests/Real_Terrain/IshiharaZhou/ishihara_zhou.py plot plt --surface smooth --outdir . --out ridge_smooth.png
```

## Expected Results
With the recommended setting (8 ranks, Apple M2 Pro):

| Surface | Reversed flow near the ground | Tunnel reattachment | LES | U rms vs tunnel | Iterations, time |
|---|---|---|---|---|---|
| smooth (line at 5 m) | 42 to 170 m behind the crest | 150–200 m | 45–180 m | 0.30 m/s | 600, 40 s |
| rough (line at 7 m) | 47 to 173 m | 250–300 m | 37–239 m | 0.58 m/s | 417, 31 s |

With the default closure the smooth ridge does not separate (U rms 0.80 m/s) and the rough one separates from 59
to 122 m (0.70 m/s). On the rough surface the tunnel's lowest points at the foot of the ridge (7 m) are 60 % below
the tunnel's own approach profile at that height (21 % on the smooth surface, which GUST reproduces); they may lie
in the roughness layer of the tunnel's roughness elements.

## CTest
None.
