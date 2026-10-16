Case ladder
===========

Each step up the ladder adds one difficulty to the one below it, so a failure points to the feature just added.
A new solver or option is run up the ladder in order. Earlier sub-phases use some of these cases for checks on a
fixed state (for example RT-1.8.1 on the flat case); the last column is the first run to a steady solution.

.. list-table::
   :header-rows: 1
   :widths: 5 22 37 16 20

   * - Step
     - Case
     - What it adds
     - Steady expected?
     - First run to steady
   * - 1
     - 1-D neutral column
     - Surface layer, k–L closure and vertical diffusion only; the log law and equilibrium k are known exactly.
     - Yes
     - 1.15 (RT-1.15.1)
   * - 2
     - Flat 3-D periodic (ERF ``Neutral_ABL_Flat``)
     - Horizontal directions, pressure, boundary fluxes.
     - Yes
     - 1.15 (RT-1.15.2)
   * - 3
     - 2-D hill (ERF ``Neutral_Hill_2D``)
     - Terrain-following mesh and metric terms; wall distance.
     - Yes
     - 1.16 (RT-1.16.1)
   * - 4
     - 3-D hill (ERF ``Neutral_Hill_3D``)
     - Flow around as well as over the hill.
     - Yes
     - 1.16 (RT-1.16.2)
   * - 5
     - Askervein
     - Real terrain and field measurements (lines A and AA).
     - Yes
     - 1.16 (RT-1.16.3)
   * - 6
     - Bolund, Perdigão
     - Steep escarpment (Bolund); double ridge (Perdigão).
     - **Possibly not** (lee-side separation)
     - 4.3 (RT-4.3.3)

The ERF Canonical_RANS set also has stratified flat cases:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Case
     - Steady expected?
     - Reason
   * - ``Stable_ABL_Flat``
     - **Possibly not**
     - Stable stratification.
   * - ``Convective_ABL_Flat``
     - Yes, in RANS
     - The closure carries the convective mixing; recheck if a run oscillates.

A case run with Coriolis forcing is labelled *possibly non-steady* whatever its row says, because of inertial
oscillations. A *possibly non-steady* case is still run. If it does not meet the definition of steady
(:doc:`steady`), it is reported as non-steady, not tuned until it passes.
