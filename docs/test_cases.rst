Test cases
==========

ablsteady is verified on a sequence of cases of increasing complexity. Each case adds one new feature to the case
before it, so a failure points to the feature just added. A new solver or option is tested on the cases in this
order. Most cases are ERF's Canonical_RANS cases, so ERF provides a reference solution (:doc:`references`). The input
files are in ``Tests/Canonical_RANS/``.

Steady cases
------------

.. list-table::
   :header-rows: 1
   :widths: 5 25 50 20

   * - No.
     - Case
     - What it tests
     - Test
   * - 1
     - 1-D neutral column
     - Wall model, turbulence model and vertical diffusion. The log law and the equilibrium :math:`k` are known
       exactly.
     - RT-1.15.1
   * - 2
     - Flat 3-D periodic (ERF ``Neutral_ABL_Flat``)
     - Horizontal directions, pressure, Coriolis forcing.
     - RT-1.15.2
   * - 3
     - 2-D hill (ERF ``Neutral_Hill_2D``)
     - Terrain-following mesh, metric terms, wall distance.
     - RT-1.16.1
   * - 4
     - 3-D hill (ERF ``Neutral_Hill_3D``)
     - Flow around as well as over a hill.
     - RT-1.16.2
   * - 5
     - Askervein hill
     - Real terrain and field measurements (masts on lines A and AA).
     - RT-1.16.3
   * - 5b
     - Askervein, conventionally neutral
     - Case 5 with stratification aloft, Coriolis forcing and a Rayleigh layer (below).
     - ``rt_steady_askervein_cnbl``
   * - 6
     - Bolund, Perdigão
     - A steep escarpment (Bolund) and a double ridge (Perdigão). Flow separation in the lee may prevent a steady
       solution.
     - RT-4.3.3 (planned)

**Case 5b.** Case 5 keeps the potential temperature, gravity, Coriolis forcing at 57.2° N and the flat case's
sounding: an inversion from 700 m, 0.03 K/m above, the top raised to 1599 m and Rayleigh damping above 1200 m. The
inflow is the steady mean of a precursor column (``inputs_precursor_askervein_cnbl``), whose geostrophic wind turns
the 10 m wind onto +x with :math:`u_* = 0.6` m/s, the inflow of case 5. The wind veers with height, so the flow enters
through the west and north faces and leaves through the east and south faces. There is no ERF reference yet.

Transient cases
---------------

ERF's Canonical_RANS set also has two stratified flat cases. They have no steady solution, so ablsteady does not
use them as steady references:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Case
     - Why it is not steady
   * - ``Stable_ABL_Flat``
     - GABLS1 cools the surface throughout, and its low-level jet is an inertial oscillation.
   * - ``Convective_ABL_Flat``
     - The mixed layer deepens and warms throughout.

They can be compared through quasi-steady quantities (:math:`u_*`, jet height, boundary-layer depth).

Notes
-----

- A case with Coriolis forcing started from a uniform geostrophic wind oscillates with the inertial period
  :math:`2\pi/f` (17.5 h at :math:`f = 10^{-4}` s\ :sup:`-1`) before it settles. Its ERF reference therefore runs
  48 h (:doc:`references`).
- A case that may not be steady (case 6) is still run. If it does not meet the convergence criteria
  (:ref:`convergence-criteria`), it is reported as not steady.

Results for these cases are in :doc:`results`.
