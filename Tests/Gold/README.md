# Gold files

| Gold | What | How it was made |
|---|---|---|
| `erf_hill2d_btf/`, `erf_hill2d_stf/`, `erf_hill2d_sullivan/` | ERF's `z_phys`, `detJ` and nodal heights (`Level_0/Nu_nd`, component 2 = z - k dz) for UT-1.4.7 | `tools/make_erf_mesh_gold.sh` with ERF at the commit in `external/VERSIONS` |
| `askervein_mesh.json` | Statistics and quality report for RT-1.4.1 / RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/askervein_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_askervein terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_mesh.json` | Quality report for RT-1.4.2 | `Tests/check_mesh.py Tests/Gold/steep_ridge_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_steep_ridge` |
| `askervein_hierarchy_mesh.json` | Level 0's statistics and every level's quality report for RT-2.3.1 / RT-2.3.2 | `mpiexec -n 2 Tests/check_mesh.py Tests/Gold/askervein_hierarchy_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_askervein_hierarchy terrain.file=extern_erf/askervein.txt` |
| `askervein_hierarchy_level1_plt00000/`, `askervein_hierarchy_level2_plt00000/` | z and detJ of hierarchy levels 1 and 2, RT-2.3.1 | the same run's `askervein_hier_plt_level{1,2}_00000`, without `Level_0/z_phys_nd*` |
| `askervein_hierarchy_p0_mesh.json` | The same with unsmoothed coarse grounds (the coarse-level hint), RT-2.3.2 | `mpiexec -n 2 Tests/check_mesh.py Tests/Gold/askervein_hierarchy_p0_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_askervein_hierarchy_p0 terrain.file=extern_erf/askervein.txt` |
| `steep_ridge_hierarchy_mesh.json` | Every level's quality report for RT-2.3.2 | `Tests/check_mesh.py Tests/Gold/steep_ridge_hierarchy_mesh.json --write -- ablsteady Tests/Decks/inputs_mesh_steep_ridge_hierarchy` |
| `hill2d_plt00000/` | Plotfile for RT-1.4.3 | `ablsteady Tests/Decks/inputs_mesh_hill2d` |
| `bc_neutral_flat_plt00000/` | Boundary fluxes and residual of the fixed log-law state, RT-1.8.1 | `boundary_fluxes Tests/Decks/inputs_bc_neutral_flat` |
| `closure_column_plt00000/` | Full residual of the equilibrium log-law column with the k-L closure, RT-1.9.1 | `closure_column Tests/Decks/inputs_closure_column` |

The ERF golds are a cross-code check: they come from ERF, not from ablsteady. The others are change detectors made
by ablsteady itself. A gold changes only in a reviewed commit that says why.
