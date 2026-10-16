# Gold files

| Gold | What | How it was made |
|---|---|---|
| `erf_hill2d_btf/`, `erf_hill2d_stf/`, `erf_hill2d_sullivan/` | ERF's `z_phys`, `detJ` and nodal heights (`Level_0/Nu_nd`, component 2 = z - k dz) for UT-1.4.7 | `tools/make_erf_mesh_gold.sh` with ERF at the commit in `external/VERSIONS` |
| `askervein_mesh.json` | Statistics and quality report for RT-1.4.1 / RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/askervein_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_askervein terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_mesh.json` | Quality report for RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/steep_ridge_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_steep_ridge` |
| `hill2d_plt00000/` | Plotfile for RT-1.4.3 | `ablsteady Tests/Decks/inputs_mesh_hill2d` |

The ERF golds are a cross-code check: they come from ERF, not from ablsteady. The others are change detectors made
by ablsteady itself. A gold changes only in a reviewed commit that says why.
