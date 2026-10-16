Code from other projects
========================

ablsteady copies some ERF files rather than linking ERF, uses AMReX and GoogleTest as libraries, and
re-implements methods from FastWindTerrain. This page lists what comes from where and how copies are kept honest.

The extraction manifest
-----------------------

Every file copied from ERF lives under ``extern_erf/`` and has one row in ``extern_erf/EXTRACTED.md``: the local
file, its ERF source path, the ERF commit it was copied from, its status, and notes.

- Status ``unchanged``: the local file is byte-identical to the ERF source at that commit.
- Status ``modified``: the local changes are recorded in ``extern_erf/patches/<file>.patch``, a unified diff from
  the ERF file to the local file.

Copied so far (sub-phase 1.4): the mesh and terrain-metric code (``ERF_GridUtils.H``, ``ERF_TerrainMetrics.{H,cpp}``,
``ERF_ChopGrids.cpp`` and the small headers they include), ERF's Python plotfile reader for the check scripts, and
the Askervein terrain file. ``ERF_TerrainMetrics.{H,cpp}`` are modified: the coordinate builders take the smoothing
choice and :math:`\gamma_m` as arguments instead of reading ERF's inputs, return the number of STF levels capped,
and lose the refined-level (AMR) paths, which ablsteady does not use. BTF no longer rewrites the ground layer: ERF
writes :math:`h\,z_\mathrm{top}/z_\mathrm{top}` there (within 1 ulp of :math:`h`) while other GPU threads read the
same nodes, so the result depended on thread order. Two shim headers in ``Source/erf_shims/`` stand
in for ERF headers the copies include.

The drift check
---------------

``tools/check_extern_erf.py --erf <ERF clone>`` fails when:

- an ``unchanged`` file differs from its ERF source;
- a ``modified`` file has no patch, or its patch no longer matches the local changes;
- a file under ``extern_erf/`` (other than ``EXTRACTED.md``, ``LICENSE.md`` and ``patches/``) has no row, or a row
  names a missing file;
- the manifest table lacks its header row.

With ``--report-upstream <ref>`` it also lists the ERF commits that touched each copied file after its pinned
commit. That list is a report: upstream fixes (for example ERF#4272 and ERF#4279) are reviewed and copied over
deliberately, not automatically. CI runs the check on every pull request.

Licences and attribution
------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Project
     - How ablsteady uses it
     - Licence and what it requires
   * - ERF
     - Copied files in ``extern_erf/``; ``Tests/Gold/erf_hill2d_*`` are ERF outputs at the pinned commit.
     - BSD-style licence with the LBNL "Enhancements" clause (ERF ``LICENSE.md``). Redistributed source must keep
       ERF's copyright notice, conditions and disclaimer: ERF's ``LICENSE.md`` is copied into ``extern_erf/`` with
       the first copied file, and copied files keep their headers.
   * - AMReX
     - Git submodule, built as a library.
     - BSD-style licence (``external/amrex/LICENSE``), kept with the submodule.
   * - GoogleTest
     - Git submodule, used by the tests only.
     - BSD 3-Clause (``external/googletest/LICENSE``), kept with the submodule.
   * - FastWindTerrain
     - Method knowledge only (mass-consistent solve, anisotropy, inflow profiles, O'Brien correction),
       re-implemented. No code is copied.
     - —
   * - Kynema-SGF (formerly AMR-Wind)
     - Cross-check of the Axell–Liungman closure in unit tests (sub-phase 1.9). No code is copied.
     - —

ablsteady itself is distributed under the BSD 2-Clause licence (``LICENSE``).
