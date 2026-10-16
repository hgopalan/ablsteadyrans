Reference solutions
===================

ablsteady is checked against ERF. ERF runs each case of the ladder (:doc:`cases`) as an unsteady RANS (URANS)
simulation on Kestrel; the time average over the final window is the reference a steady ablsteady solution is
compared with. The wall time and allocation-unit (AU) cost of each ERF run is the baseline ablsteady has to beat.

Pinned versions
---------------

``external/VERSIONS`` pins the commits ablsteady builds on or compares with:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Name
     - Commit
     - Role
   * - ``erf``
     - ``92410d57``
     - Reference solutions; source of files copied to ``extern_erf/``.
   * - ``amrex``
     - ``347a3a2d``
     - The ``external/amrex`` submodule. ``tools/check_versions.py`` (CTest ``versions_check``) fails if the two
       differ.
   * - ``fastwindterrain``
     - ``648f8f1e``
     - Method of the mass-consistent initializer; no code copied.
   * - ``kynema-sgf``
     - ``8500c568``
     - Cross-check of the Axell–Liungman closure; no code copied.

The file holds the full 40-character commits.

Reference runs
--------------

``tools/kestrel/erf_references.sh`` builds ERF at the pinned commit, runs the cases and archives them:

.. code-block:: bash

   export ABLSTEADY_ACCOUNT=<allocation>
   tools/kestrel/build_kestrel.sh cpu                 # ablsteady and ablsteady_refavg
   tools/kestrel/erf_references.sh build cpu          # ERF at the pinned commit
   salloc -A $ABLSTEADY_ACCOUNT -N 1 -t 01:00:00
   tools/kestrel/erf_references.sh smoke cpu          # RT-1.1.1, inside the allocation
   tools/kestrel/erf_references.sh submit cpu all
   tools/kestrel/erf_references.sh archive cpu Neutral_ABL_Flat   # after each job ends

The same commands with ``gpu`` build and run on H100 nodes; ``archive`` always runs the CPU build of
``ablsteady_refavg`` under the CPU modules, so run it from a CPU login node. Each run writes a plotfile every
window/24 seconds of physical time (``erf.plot_per_1``) and no 2-D plotfiles. A run that reaches its time limit is
rerun from scratch; there is no restart support.

.. list-table::
   :header-rows: 1
   :widths: 22 12 11 13 42

   * - Case
     - Run (s)
     - Window (s)
     - CPU ranks / GPUs
     - Note
   * - ``Neutral_ABL_Flat``
     - 172 800
     - 3 600
     - 8 / 1
     - 48 h instead of the deck's 12 h (inertial oscillation, below).
   * - ``Stable_ABL_Flat``
     - 32 400
     - 3 600
     - 8 / 1
     - **Transient by design** (GABLS1 surface cooling): not a steady reference. Deck length.
   * - ``Convective_ABL_Flat``
     - 14 400
     - 3 600
     - 8 / 1
     - **Transient by design** (the mixed layer grows and warms): not a steady reference. Deck length.
   * - ``Neutral_Hill_2D``
     - 172 800
     - 3 600
     - 32 / 1
     - 48 h instead of 6 h.
   * - ``Neutral_Hill_3D``
     - 172 800
     - 3 600
     - 64 / 1
     - 48 h instead of 4 h.
   * - Askervein (``inputs_anel``)
     - 4 800
     - 2 400
     - 104 / 4
     - The deck stops at 600 s; the lee flow takes about 2 000 s to pass through the 7 km path. 4 800 s with the
       last 2 400 s averaged.

**Why 48 h.** The neutral flat case and both hills start from a uniform geostrophic wind with Coriolis forcing on.
That sets off an inertial oscillation with period :math:`2\pi/f` = 17.5 h (:math:`f = 10^{-4}` s\ :sup:`-1`), which
a one-hour average does not remove: averaged over one hour it keeps 99.5 % of its amplitude. A 1-D Ekman spin-up
estimate puts the wind in the last hour 1.3–1.4 m/s from steady after 4–6 h, 0.8 m/s after 12 h and 0.03 m/s after
48 h (about 2.75 periods).

RT-1.1.1 (``erf_references.sh smoke``) runs ERF's own ``ctest -L rans`` at the pinned commit on Kestrel before any
reference is produced.

Time average
------------

For snapshots at times :math:`t_0 < t_1 < \dots < t_n` in the window, the reference is the trapezoid-rule time
average

.. math::
   :label: trapezoid-average

   \bar f = \sum_{i=0}^{n} w_i f_i ,
   \qquad
   w_i = \frac{t_{i+1} - t_{i-1}}{2\,(t_n - t_0)} ,

with :math:`t_{-1} = t_0` and :math:`t_{n+1} = t_n`. The weights handle unevenly spaced snapshots, which a run with a
CFL-controlled time step produces. The sum is accumulated in snapshot order with an explicit fused multiply-add, so
its values are the same on every machine and rank count. A snapshot containing NaN is refused. At least 3 snapshots
are needed.

**Stationarity.** The same average is also taken over each half of the window, :math:`[t_0, t_m]` and
:math:`[t_m, t_n]` with :math:`m = \lfloor n/2 \rfloor`. The archive records, per variable, the largest change
:math:`\max |\bar f_2 - \bar f_1|` over the domain and the largest :math:`|\bar f|` for scale. A run still drifting
(or oscillating with a period comparable to the window) shows a large change; a steady one shows none.
``ablsteady_refavg refavg.verify=...`` prints both for every variable.

**What the average cannot remove.** An oscillation with a period of exactly window/24 (or a divisor of it) is
sampled at the same phase every time and stays in the average; other periods below half the window leave at most a
few per cent; periods longer than half the window leave 20–99 %. The half-window change makes the last case
visible.

Archive
-------

``ablsteady_refavg`` writes an archive directory:

.. code-block:: text

   <case>_<cpu|gpu>/
     average/     AMReX plotfile of the time average (level 0, all variables of the snapshots)
     inputs/      the deck, its input files (sounding, terrain, inflow) and the command-line overrides
     MANIFEST     case, ERF commit, machine, nodes, ranks, wall time, AU, the snapshots, per-variable
                  sum, minimum, maximum and half-window change, and the SHA-256 of every other file

The AU are nodes × wall hours × the machine's rate: 10 AU per node-hour on Kestrel CPU nodes and 100 on GPU nodes,
as in the development plan (``local`` is 0). The wall time is that of ``srun``, so it excludes queueing and job
start-up. A run on part of a node is counted as a whole node, so the AU are an upper bound on the charge; ablsteady
runs are counted the same way. The AU cover the whole run, including the spin-up length chosen above.

Before any use, ``load_reference`` checks the archive. It fails, naming the file or variable, if:

- a file differs from its SHA-256 in ``MANIFEST``;
- a listed file is missing, or a file is present but not listed;
- the AU do not equal nodes × wall hours × rate;
- ``MANIFEST`` does not parse;
- the fields as read do not have the recorded sum, minimum and maximum (exactly).

The last check ties the data to the manifest independently of the files: the averaged plotfile is written in one
file per writing rank, so an archive regenerated on another rank count has the same values but different files and
checksums.

The average keeps only cell-centred variables. On terrain cases the nodal heights ERF writes beside its plotfiles
are not kept; the cell-centred ``z_phys`` is, when the deck plots it (both hills and Askervein do). In a
single-precision build of ablsteady the snapshot times are rounded to single precision.

``ablsteady_refavg``
--------------------

The tool owns the ``refavg`` namespace; any other key stops it, as in :doc:`inputs`.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``refavg.plotfiles``
     - required
     - At least 3 single-level plotfiles with the same variables, domain and box layout, times strictly
       increasing. A trailing ``/`` is allowed.
   * - ``refavg.output``
     - required
     - Archive directory to create; must not exist.
   * - ``refavg.case``
     - required
     - One word.
   * - ``refavg.erf_sha``
     - required
     - Full 40-character lowercase hex commit.
   * - ``refavg.machine``
     - ``local``
     - ``kestrel-cpu``, ``kestrel-gpu`` or ``local``.
   * - ``refavg.nodes``
     - ``1``
     - ≥ 1.
   * - ``refavg.ranks``
     - ``1``
     - ≥ 1. MPI ranks (CPU) or GPUs of the run.
   * - ``refavg.wall_seconds``
     - ``0``
     - Finite, ≥ 0.
   * - ``refavg.inputs_files``
     - empty
     - Files to copy into ``inputs/``; each must exist, and no two may share a name.
   * - ``refavg.verify``
     - empty
     - Instead of building: check this archive (files and field values), print its metadata and half-window
       changes, and stop with an error if it does not verify. The other ``refavg.*`` keys have no effect here, so
       giving one stops the tool.

Status
------

The scripts have not yet run on Kestrel, so no archive exists yet. Everything that does not need Kestrel (the
averaging, the archive format and its checks, the tool) is tested locally and in CI (:doc:`testing`).
