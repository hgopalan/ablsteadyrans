# Code copied from ERF

GUST copies some ERF files rather than linking ERF (development plan rev. 6, extraction audit).
Every copied file lives under `extern_erf/` and has one row below. `tools/check_extern_erf.py` reads this table:

- **Status `unchanged`:** the local file must be byte-identical to `ERF source` at `ERF SHA`.
- **Status `modified`:** the local changes are recorded in `extern_erf/patches/<local file name>.patch`, a unified
  diff from the ERF file to the local file. The check regenerates the diff and fails if the stored patch is out of date.
- Every file under `extern_erf/` (except this file, `LICENSE.md` and `patches/`) must have a row, and every row must name a file
  that exists.

The CI drift job runs the check against an ERF clone and also lists upstream ERF changes to each source file since
its pinned SHA. That list is a report; it never fails the job.

ERF is distributed under a BSD 3-clause licence with the LBNL "Enhancements" clause (ERF `LICENSE.md`). Condition (1)
of that licence requires redistributed source to keep ERF's copyright notice, conditions and disclaimer, so a copy of
ERF's `LICENSE.md` goes into `extern_erf/` with the first copied file. Copied files keep their original headers.

<!-- manifest:start -->
| Local file | ERF source | ERF SHA | Status | Notes |
|---|---|---|---|---|
| `extern_erf/ERF_GridUtils.H` | `Source/Utils/ERF_GridUtils.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Geometric vertical stretch (build_geometric_zlevels) |
| `extern_erf/ERF_Constants.H` | `Source/ERF_Constants.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | bogus_large_value and constants used by ERF_TerrainMetrics.cpp |
| `extern_erf/ERF_NumericalConstants.H` | `Source/ERF_NumericalConstants.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | zero, one, myhalf, PI |
| `extern_erf/ERF_Interpolation_1D.H` | `Source/Utils/ERF_Interpolation_1D.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Included by ERF_TerrainMetrics.cpp |
| `extern_erf/ERF_TerrainMetrics.H` | `Source/Utils/ERF_TerrainMetrics.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | modified | Terrain-fitted coordinate builders take terrain_smoothing, terrain_gamma_m and return the STF capped count |
| `extern_erf/ERF_TerrainMetrics.cpp` | `Source/Utils/ERF_TerrainMetrics.cpp` | 92410d578f612ab9e844c450bf8e900194becfb5 | modified | Same signatures; no ParmParse reads; lev > 0 (AMR) paths removed; STF note names terrain.gamma_m; BTF leaves the ground layer inside the domain as given (ERF's in-place rewrite is a GPU race) and sets its ghost columns outside non-periodic sides from the clamped domain node |
| `extern_erf/ERF_ChopGrids.cpp` | `Source/Utils/ERF_ChopGrids.cpp` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | ChopGrids2D, referenced by make_terrain_fitted_coords |
| `extern_erf/erf_plotfile.py` | `Exec/CanonicalTests/Canonical_RANS/erf_plotfile.py` | 92410d578f612ab9e844c450bf8e900194becfb5 | modified | Plotfile reader for the regression check scripts; the real size and byte order from the FAB's byte-order descriptor (the upstream reader takes the format's field count, 8 for float too, as the size, and reads float plotfiles as byte-swapped doubles) |
| `extern_erf/check_neutral.py` | `Exec/CanonicalTests/Canonical_RANS/Neutral_ABL_Flat/check_neutral.py` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | ERF's Neutral_ABL_Flat checks, run on the stage's turbulence plotfile (RT-1.15.2) |
| `extern_erf/check_hill.py` | `Exec/CanonicalTests/Canonical_RANS/Neutral_Hill_2D/check_hill.py` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | ERF's Neutral_Hill_2D checks, run on the stage's turbulence plotfile (RT-1.16.1) |
| `extern_erf/check_hill3d.py` | `Exec/CanonicalTests/Canonical_RANS/Neutral_Hill_3D/check_hill3d.py` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | ERF's Neutral_Hill_3D checks (RT-1.16.2) |
| `extern_erf/rans_checks.py` | `Exec/CanonicalTests/Canonical_RANS/rans_checks.py` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Shared report and structural checks of ERF's RANS check scripts |
| `extern_erf/askervein.txt` | `Exec/CanonicalTests/Real_Terrain/Askervein/askervein.txt` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Askervein terrain (ERF's structured terrain format) |
| `extern_erf/inflow_powerlaw.txt` | `Exec/CanonicalTests/Real_Terrain/Askervein/inflow_powerlaw.txt` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Askervein power-law inflow, wind from 210° (z u v w); source of `Tests/Real_Terrain/Askervein_210/askervein_210_inflow.txt` (RT-2.6.1) |
| `extern_erf/ERF_TerrainPoisson.H` | `Source/LinearSolvers/ERF_TerrainPoisson.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Terrain Poisson operator for amrex::GMRES with the FFT preconditioner (sub-phase 1.10) |
| `extern_erf/ERF_TerrainPoisson.cpp` | `Source/LinearSolvers/ERF_TerrainPoisson.cpp` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Includes the zero-flux Neumann faces with the terrain cross term (ERF#4272, merged as ERF#4279) |
| `extern_erf/ERF_TerrainPoisson_3D_K.H` | `Source/LinearSolvers/ERF_TerrainPoisson_3D_K.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | Operator and face-flux kernels |
| `extern_erf/ERF_SolverUtils.H` | `Source/LinearSolvers/ERF_SolverUtils.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | get_fft_bc: maps the domain boundary names to FFT::Boundary (Outflow/Open: odd, periodic: periodic, others: even) |
| `extern_erf/ERF_MOSTUtils.H` | `Source/BoundaryConditions/ERF_MOSTUtils.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | similarity_funs (Businger-Dyer and Jimenez stability functions) for the surface layer (sub-phase 1.7); ERF#3982 (reference-height lookup) is already in this ERF |
| `extern_erf/ERF_RANSClosure.H` | `Source/Diffusion/ERF_RANSClosure.H` | 92410d578f612ab9e844c450bf8e900194becfb5 | unchanged | AL01 k-L closure relations (geom_length, smooth_Rt, turb_length, cmu, cmu_prime): the reference of UT-1.9.1 |
<!-- manifest:end -->

The shims in `Source/erf_shims/` stand in for two ERF headers the copied files include (`ERF_IndexDefines.H`: only
the `ERF_BC` enum; `ERF_Utils.H`: only declarations of functions defined here).
