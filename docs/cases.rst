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
     - Real terrain and field measurements (lines A and AA). The ERF reference runs 3600 s, longer than the
       deck's 600 s (:doc:`references`).
     - Yes
     - 1.16 (RT-1.16.3)
   * - 5b
     - Askervein, conventionally neutral (beyond the plan)
     - Step 5 with the potential temperature, gravity, Coriolis (57.2 N) and the flat case's sounding kept: an
       inversion from 700 m, 0.03 K/m above, the top raised to 1599 m with Rayleigh damping above 1200 m. The
       inflow is the steady planar mean of a precursor column (``inputs_precursor_askervein_cnbl``), whose
       geostrophic wind turns the 10 m wind onto +x and sets u* = 0.6 m/s, the inflow of step 5. The wind veers to
       −y with height, so the flow enters through xlo and yhi and leaves through xhi and ylo. No ERF reference
       yet: step 5 stays the comparison case.
     - Yes
     - Beyond the plan (``rt_steady_askervein_cnbl``)
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
     - **No**
     - Transient by design: GABLS1 cools the surface throughout, and its low-level jet is an inertial
       oscillation. Compared through quasi-steady measures (u*, jet height, boundary-layer depth), not as a
       steady field.
   * - ``Convective_ABL_Flat``
     - **No**
     - Transient by design: the mixed layer deepens and warms throughout. Compared through quasi-steady
       measures, not as a steady field.

A case run with Coriolis forcing from a uniform geostrophic start oscillates with the inertial period
:math:`2\pi/f` (17.5 h at :math:`f = 10^{-4}` s\ :sup:`-1`) before it settles; its ERF reference runs 48 h
(:doc:`references`). A *possibly non-steady* case is still run. If it does not meet the definition of steady
(:doc:`steady`), it is reported as non-steady, not tuned until it passes.
