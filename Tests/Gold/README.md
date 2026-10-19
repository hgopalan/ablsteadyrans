# Gold files

| Gold | What | How it was made |
|---|---|---|
| `erf_hill2d_btf/`, `erf_hill2d_stf/`, `erf_hill2d_sullivan/` | ERF's `z_phys`, `detJ` and nodal heights (`Level_0/Nu_nd`, component 2 = z - k dz) for UT-1.4.7 | `tools/make_erf_mesh_gold.sh` with ERF at the commit in `external/VERSIONS` |
| `askervein_mesh.json` | Statistics and quality report for RT-1.4.1 / RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/askervein_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_askervein terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_mesh.json` | Quality report for RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/steep_ridge_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_steep_ridge` |
| `hill2d_plt00000/` | Plotfile for RT-1.4.3 | `ablsteady Tests/Decks/inputs_mesh_hill2d` |
| `bc_neutral_flat_plt00000/` | Boundary fluxes and residual of the fixed log-law state, RT-1.8.1 | `boundary_fluxes Tests/Decks/inputs_bc_neutral_flat` |
| `closure_column_plt00000/` | Full residual of the equilibrium log-law column with the k-L closure, RT-1.9.1 | `closure_column Tests/Decks/inputs_closure_column` |

`hill2d_plt00000/Level_0/z_phys_nd_D_00000` (the node heights with ghost nodes; no test reads it) was rewritten when
the ghost columns were made to match their periodic images: 8 values of the row below the ground changed.

The ERF golds are a cross-code check: they come from ERF, not from ablsteady. The others are change detectors made
by ablsteady itself. A gold changes only in a reviewed commit that says why.
