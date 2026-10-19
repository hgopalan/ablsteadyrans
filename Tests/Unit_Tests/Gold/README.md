# Gold files

| Gold | What | How it was made |
|---|---|---|
| `erf_hill2d_btf/`, `erf_hill2d_stf/`, `erf_hill2d_sullivan/` | ERF's `z_phys`, `detJ` and nodal heights (`Level_0/Nu_nd`, component 2 = z - k dz) for UT-1.4.7 | `tools/make_erf_mesh_gold.sh` with ERF at the commit in `external/VERSIONS` |

These golds are a cross-code check: they come from ERF, not from ablsteady. The regression golds written by
ablsteady live in each case folder (`Tests/<group>/<case>/Gold/`), with a "Gold files" table in that case's
README. A gold changes only in a reviewed commit that says why.
