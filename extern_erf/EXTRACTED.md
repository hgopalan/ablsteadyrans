# Code copied from ERF

ablsteady copies some ERF files rather than linking ERF (development plan rev. 6, extraction audit).
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
<!-- manifest:end -->

No files are copied yet. The first copies arrive with sub-phase 1.4 (mesh and terrain metrics).
