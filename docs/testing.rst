Testing
========

Conventions
-----------

These follow the development plan.

- **Unit tests (UT-p.s.n)** use GoogleTest. They are deterministic, finish in seconds on 1–4 ranks, and run on
  every commit: ``ctest -L unit``.
- **Regression tests (RT-p.s.n)** run the application and compare with stored gold plotfiles using
  ``amrex_fcompare``. Their check scripts print *measured / target / tolerance / pass* rows.
  ``ctest -L smoke`` runs the structural checks CI does on every pull request; ``ctest -L physics`` runs the
  converged checks nightly on Kestrel.
- **Performance tests (PT-p.s.n)** run nightly or weekly on Kestrel. They track iterations, wall time and
  allocation-unit cost; a slowdown of more than 10 % fails.
- **Gold files** change only in a reviewed commit.

Every test is first shown to fail on the bug it guards against, or on a one-line mutant of the code it covers.
Tests the plan does not list are marked *beyond the plan*.

The tests run with the GPU CI flags (see :doc:`running`). Each test runs in its own directory under
``<build>/Tests/run/``.

Current tests
-------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - CTest
     - What it checks
   * - ``unit_np1``, ``unit_np2``, ``unit_np4``
     - The GoogleTest suite on 1, 2 and 4 ranks (below).
   * - ``unit_notiling``
     - The same suite with tiling off, as on a GPU.
   * - ``smoke_clean``, ``smoke_run``, ``smoke_restart``, ``smoke_restart_matches``
     - A 2-rank run; a 3-rank restart on another box layout, whose log must show the restart; the two plotfiles
       identical under ``amrex_fcompare`` at zero tolerance.
   * - ``abort_unknown_key``, ``abort_range``, ``abort_repeated_key``, ``abort_restart_mesh``,
       ``abort_restart_overwrite``
     - Decks that must stop, with the input named in the message.
   * - ``rt_askervein_mesh``
     - RT-1.4.1 and RT-1.4.2: ERF's Askervein mesh coarsened to 60 × 60 × 18 (20 m × 1.08 layers, BTF, ERF's terrain
       file): statistics of ``z`` and ``detJ`` and the quality report against ``Tests/Gold/askervein_mesh.json``.
       Double-precision builds only: the statistics use ERF's plotfile reader, which reads double-precision
       plotfiles.
   * - ``rt_steep_ridge_quality``
     - RT-1.4.2: a 600 m ridge under STF: the report against its gold, including both warnings (STF capped, top
       not flat). Double-precision builds only.
   * - ``rt_hill2d_run``, ``rt_hill2d_matches``
     - RT-1.4.3: a run without a solve writes the mapped-grid plotfile of the 2-D hill; ``amrex_fcompare`` against
       ``Tests/Gold/hill2d_plt00000`` at relative tolerance 1e-12 (1e-5 in single precision). The plotfile holds
       cell fields only, so the nodal heights are not compared here; UT-1.4.7 compares them with ERF.

       The golds of the three mesh tests above were written by ablsteady itself (``Tests/check_mesh.py --write``):
       they detect changes, they do not show correctness. Correctness of the same code paths comes from UT-1.4.7
       (ERF's own output) and the closed-form unit tests.
   * - ``drift_check_selftest``
     - The ``extern_erf`` drift check on a throw-away ERF repository (:doc:`extraction`).
   * - ``refavg_snapshots``, ``refavg_build``, ``refavg_verify``, ``refavg_matches``
     - ``ablsteady_refavg`` end to end: three copies of the smoke plotfile at 0, 5 and 10 s are archived on 2
       ranks with two inputs files; the archive verifies on 3 ranks with the expected metadata, AU and a zero
       half-window change; the average equals the smoke plotfile at zero tolerance. This checks the round trip;
       identical snapshots cannot check the weights, which UT-1.1.2 does.
   * - ``refavg_damage``, ``abort_refavg_damaged``, ``abort_refavg_inputs``, ``abort_refavg_verify_keys``
     - A copy of the archive with one bit flipped is refused, naming the file; bad tool inputs stop the tool with
       each input named; a build key given with ``refavg.verify`` is reported.
   * - ``versions_check``, ``versions_check_selftest``
     - ``external/VERSIONS`` is well formed and its ``amrex`` entry matches the submodule; the self-test shows a
       mismatch with a throw-away submodule is reported.
   * - ``kestrel_parse_*``, ``kestrel_dryrun``
     - The Kestrel scripts parse; a dry run with stub ``module``/``sbatch``/``srun`` checks the Askervein job
       script (run length, plot period, no 2-D plotfiles, ranks) and archives a fake finished run, averaging only
       the 3-D plotfiles in the final window and archiving the inputs but not the run's outputs; a second case
       with the default prefix ``plt`` keeps its ``plt2d_`` plotfiles out.

GoogleTest suite
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - ID
     - Check
   * - UT-1.1.1
     - A reference archive written, verified and loaded: fields bit-identical to the average, metadata and field
       statistics unchanged; a changed byte, a missing file, an extra file, edited AU, and data changed together
       with its checksum are each found and named, and loading refuses.
   * - UT-1.1.2
     - *Beyond the plan.* The average equals a host-side recomputation bit for bit and, for fields linear in
       time, the mid-window value to 8 ε; trapezoid weights for uneven spacing; the half-window change is exact
       for a linear drift and zero for a steady field; a snapshot with a NaN, fewer than 3 snapshots and times out of
       order are refused with their messages; a trailing ``/`` is handled.
   * - UT-1.1.3
     - *Beyond the plan.* AU = nodes × hours × rate (10 CPU, 100 GPU); metadata range checks name their inputs.
   * - UT-1.1.4
     - *Beyond the plan.* SHA-256 against the NIST example vectors, at every padding boundary and chunking;
       ``hex()`` can be called repeatedly.
   * - UT-1.4.1
     - Vertical levels: uniform; solved ratio (top exact, ratio 1.0421553943273616 from an independent exact-rational bisection);
       fixed ratio (ERF's Askervein 20 m × 1.08, top 749.0048747962335 m); z-levels file and its errors.
   * - UT-1.4.2
     - Flat terrain gives the Cartesian levels: bit for bit for STF and Sullivan, within 2 ulp for BTF.
   * - UT-1.4.3
     - Free-stream preservation on round Witch-of-Agnesi and Gaussian hills, all three mappings: divergence below
       1e-12 of :math:`(|U|+|V|)/\Delta x`. The cancellation is algebraic for any node heights, so this test does
       not detect wrong heights; UT-1.4.7 does.
   * - UT-1.4.4
     - Cell volumes sum to the column volumes above the ground to 1e-12, and to a 1 m quadrature of the domain
       above the hill to 1e-4 (the bilinear ground differs from the continuous hill by 2.4e-5 of the domain).
   * - UT-1.4.5
     - Smoothing: a constant is unchanged; a Fourier mode is damped by the analytic factor per pass; β ≥ 1/4 is
       refused. *Beyond the plan:* on a non-periodic domain, edge nodes use themselves as the outside neighbour and
       ghost nodes take the nearest domain node.
   * - UT-1.4.6
     - IDW returns input points and, with 4 neighbours, is exact on a plane at a symmetric query (with the
       default 6 the two extra neighbours break the symmetry and IDW is not exact on a plane); ERF's bilinear form
       is exact on planes anywhere, with ERF's edge fallback; ``terrain.idw_neighbors`` above the number of points is
       refused.
   * - UT-1.4.7
     - ``z_phys``, ``detJ`` and nodal heights equal ERF's own output for the same 2-D hill deck (BTF, STF,
       Sullivan; STF and Sullivan bit for bit on macOS, BTF to rounding because ablsteady keeps the ground layer as
       given, all within 1e-13 elsewhere); STF capped at 2 of 16 levels, as ERF reports.
   * - UT-1.4.8
     - *Beyond the plan.* The IDW bucket search picks the brute-force neighbours; results agree to 1e-13 (random
       cloud with far queries, a regular grid with exact distance ties, a tie across a bucket edge, a transect
       queried from a 2-D domain, one point); UTM-sized coordinates are read and interpolated in double; file terrain on a periodic domain has one
       seam height; every mapping keeps the ground layer as given; ground at the top and inverted cells stop the
       run; the slope of a diagonal ridge is reported in full; a restart with changed terrain is refused; the report
       of a steep ridge carries both warnings.
   * - UT-1.2.1
     - Steady-criterion verdicts on converging, stalling (including a residual stuck exactly at the tolerance) and
       oscillating histories. The converging history becomes
       steady at iteration 136, pinned from an independent Python implementation of the definition.
   * - UT-1.2.2
     - *Beyond the plan.* A NaN or infinite residual or monitor gives ``Diverged``, permanently.
   * - UT-1.2.3
     - *Beyond the plan.* A monitor at round-off level around zero is steady with scale 1 and never steady with a
       scale below the noise.
   * - UT-1.2.4
     - *Beyond the plan.* The monitor window needs :math:`W+1` values; a change of exactly the tolerance passes and
       one just above it does not (values exact in binary).
   * - UT-1.2.5
     - *Beyond the plan.* The residual drop is measured from the early peak (iteration 38), not from the first
       value (iteration 49).
   * - UT-1.2.6
     - *Beyond the plan.* Parameter range checks name their keys; bad or missing monitor scales and size
       mismatches stop.
   * - UT-1.2.7
     - *Beyond the plan.* A monitor converging at rate 0.999 stops at iteration 4602, within the tolerance of its
       limit; the window test alone would stop at 1615 with a 2 % error.
   * - UT-1.2.8
     - *Beyond the plan.* Undamped oscillations with periods of 37, 100 and 200 iterations are never steady.
   * - UT-1.2.9
     - *Beyond the plan.* A residual that is zero over the reference iterations does not block the stop.
   * - UT-1.2.10
     - *Beyond the plan.* A criterion restored from a saved history gives the same verdict at every later
       iteration as one never interrupted; a mismatched history is refused.
   * - UT-1.2.11
     - *Beyond the plan.* A slowly growing monitor is not steady while its change over the window is still small.
   * - UT-1.3.1
     - Input parsing; unknown, repeated, multi-valued and out-of-range keys reported with their names.
   * - UT-1.3.2
     - Plotfile round trip, bit-identical, including the nodal heights and their ghost nodes.
   * - UT-1.3.3
     - Checkpoint and restart onto another box layout, bit-identical; ghost nodes restored without overwriting
       valid nodes.
   * - UT-1.3.4
     - The sum gives the same bits on 1, 2 and 4 ranks and within 16 ε Σ\|x\| of the exact sum; the maximum
       handles NaN and −0.
