# Tests

## Overview
The GUST tests. Each case has its own folder with a `README.md`, its input decks, its check script and its
gold files (`Gold/`). The layout follows ERF's `Exec/CanonicalTests/<group>/<case>/`. `CMakeLists.txt` in this
folder registers every case with CTest; each test runs in its own directory under `<build>/Tests/run/<name>`, so
nothing is written into the source tree.

## Groups
| Folder | Cases | Character |
|---|---|---|
| `Infrastructure` | Smoke run and restart, input validation, reference archives, the Kestrel scripts. | unit / regression / smoke |
| `Discretization` | The terrain mesh and its grid hierarchy, boundary fluxes, the k-L closure, the Poisson layer and the wall distance, each on a fixed state. | regression |
| `Initialization` | The mass-consistent initial field: Askervein, and the hilltop speed-up trend against FastWindTerrain. | regression |
| `Canonical_RANS` | The SIMPLEC stage to steady on ERF's canonical RANS cases: the neutral column, Neutral_ABL_Flat, Neutral_Hill_2D, Neutral_Hill_3D and its truly neutral variant, and the smooth plateau comparison case. | regression / nightly |
| `Real_Terrain` | Askervein to steady, with a log-law inflow and in a conventionally neutral boundary layer. | regression |
| `Unit_Tests` | The GoogleTest suite (`unit_tests`) and its ERF golds. | unit |
| `Supporting_Files` | Check scripts shared by more than one case, and the test-only drivers (`Tools/`). | assets |

## Running
```bash
ctest --test-dir <build> -L "unit|smoke"
```
`-L regression` adds the longer runs; `-L nightly` runs the 3-D hills. `-R <name>` runs one test; each case
README lists its CTest names.

## Notes
- Gold files change only in a reviewed commit that says why. The ERF golds (`Unit_Tests/Gold/erf_hill2d_*`) come
  from ERF; the others are change detectors written by GUST itself.
- `Gold/` folders are exempt from the whitespace and spelling checks (`.gitattributes`, `.codespellrc`).
- The full list of tests and what each checks is in `docs/testing.rst`.
- `inputs_reference` is the master input reference: every input with its default, range and meaning, as comments
  to copy into a deck. It is generated from `docs/inputs.rst` (`tools/input_reference.py --write`) and checked by
  the CTest `input_reference_current`.
