ERF reference solutions
=======================

ablsteady is validated against ERF. ERF solves the same equations with the same turbulence model, but in time
(unsteady RANS). Running ERF long enough and averaging over the final hour gives a reference steady solution for each
test case (:doc:`test_cases`). The wall time and cost of each ERF run is also the baseline that ablsteady aims to
beat.

.. note::

   The reference runs need NLR's Kestrel and have not been run yet. Everything that can be tested without Kestrel
   (the averaging, the archive format and its checks) is tested locally and in CI (:doc:`testing`).

Pinned versions
---------------

``external/VERSIONS`` pins the commits that ablsteady builds on or compares with (the file holds the full
40-character hashes):

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Project
     - Commit
     - Role
   * - ERF
     - ``92410d57``
     - Reference solutions; source of the files in ``extern_erf/``.
   * - AMReX
     - ``347a3a2d``
     - The ``external/amrex`` submodule. The CTest ``versions_check`` fails if the two differ.
   * - FastWindTerrain
     - ``648f8f1e``
     - Method of the mass-consistent initial field; no code copied.
   * - Kynema-SGF
     - ``8500c568``
     - Cross-check of the turbulence model; no code copied.

Producing the references
------------------------

``tools/kestrel/erf_references.sh`` builds ERF at the pinned commit, runs the cases and archives the results:

.. code-block:: bash

   export ABLSTEADY_ACCOUNT=<allocation>
   tools/kestrel/build_kestrel.sh cpu                 # ablsteady and ablsteady_refavg
   tools/kestrel/erf_references.sh build cpu          # ERF at the pinned commit
   salloc -A $ABLSTEADY_ACCOUNT -N 1 -t 01:00:00
   tools/kestrel/erf_references.sh smoke cpu          # ERF's own RANS tests, inside the allocation
   tools/kestrel/erf_references.sh submit cpu all
   tools/kestrel/erf_references.sh archive cpu Neutral_ABL_Flat   # after each job ends

Use ``gpu`` instead of ``cpu`` for the H100 nodes. ``archive`` always uses the CPU build, so run it from a CPU login
node. A run that reaches its time limit is rerun from the start.

Run lengths
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 24 14 14 14 34

   * - Case
     - Simulated time
     - Averaging window
     - CPU ranks / GPUs
     - Note
   * - ``Neutral_ABL_Flat``
     - 48 h
     - 1 h
     - 8 / 1
     - ERF's input file runs 12 h (see below).
   * - ``Neutral_Hill_2D``
     - 48 h
     - 1 h
     - 32 / 1
     - ERF's input file runs 6 h.
   * - ``Neutral_Hill_3D``
     - 48 h
     - 1 h
     - 64 / 1
     - ERF's input file runs 4 h.
   * - Askervein (``inputs_anel``)
     - 4800 s
     - 2400 s
     - 104 / 4
     - ERF's input file stops at 600 s; the lee flow takes about 2000 s to cross the 7 km domain.
   * - ``Stable_ABL_Flat``
     - 9 h
     - 1 h
     - 8 / 1
     - Transient by design; not a steady reference.
   * - ``Convective_ABL_Flat``
     - 4 h
     - 1 h
     - 8 / 1
     - Transient by design; not a steady reference.

**Why 48 hours?** The flat case and both hills start from a uniform geostrophic wind with Coriolis forcing. This
starts an inertial oscillation with a period of :math:`2\pi/f` = 17.5 h (:math:`f = 10^{-4}` s\ :sup:`-1`), which
a one-hour average does not remove. A 1-D estimate puts the wind 1.3–1.4 m/s from steady after 4–6 h, 0.8 m/s after
12 h and 0.03 m/s after 48 h.

Time averaging
--------------

For snapshots at times :math:`t_0 < t_1 < \dots < t_n` in the window, the reference is the trapezoidal time average

.. math::
   :label: trapezoid-average

   \bar f = \sum_{i=0}^{n} w_i f_i ,
   \qquad
   w_i = \frac{t_{i+1} - t_{i-1}}{2\,(t_n - t_0)} ,

with :math:`t_{-1} = t_0` and :math:`t_{n+1} = t_n`. The weights allow unevenly spaced snapshots. The sum is formed
in a fixed order, so the result is the same on every machine and rank count. At least 3 snapshots are needed, and a
snapshot containing NaN is refused. ERF writes a snapshot every window/24 of simulated time.

**Is the run stationary?** The same average is also taken over each half of the window. For each variable the
archive records the largest difference between the two half-window averages, and the largest :math:`|\bar f|` for
scale. A run that is still drifting, or that oscillates with a period close to the window, shows a large difference.

Archive
-------

The tool ``ablsteady_refavg`` writes one archive directory per case:

.. code-block:: text

   <case>_<cpu|gpu>/
     average/     AMReX plotfile of the time average
     inputs/      the ERF input file, its data files and the command-line overrides
     MANIFEST     case, ERF commit, machine, nodes, ranks, wall time, allocation units, the snapshots,
                  per-variable sum, minimum, maximum and half-window change, and a SHA-256 of every file

The allocation units (AU) are nodes × wall hours × the machine's rate: 10 AU per node-hour on Kestrel CPU nodes and
100 on GPU nodes. Partial nodes count as whole nodes, for ERF and ablsteady alike.

Before an archive is used, it is checked. The check fails, naming the file or variable, if a file's checksum
differs, a file is missing or not listed, the AU do not match, the manifest does not parse, or the fields do not
have exactly the recorded sum, minimum and maximum.

Only cell-centred variables are averaged. On terrain cases the cell-centre height ``z_phys`` is kept when the ERF
input file writes it (both hills and Askervein do).

The inputs of ``ablsteady_refavg`` are listed in :doc:`inputs`.
