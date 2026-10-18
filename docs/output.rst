Output
======

Plotfile
--------

``<io.plot_file>NNNNN``, where ``NNNNN`` is the step, is an AMReX plotfile with the cell-centred variables below.
It can be read with yt, ParaView, VisIt or AMReX's ``amrex_fcompare``.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Variable
     - Meaning
   * - ``u``, ``v``, ``w``
     - Velocity components: the mean of the two face values around the cell.
   * - ``p``
     - Perturbation pressure p′.
   * - ``theta``
     - Potential temperature.
   * - ``tke``
     - Turbulent kinetic energy k.
   * - ``z``
     - Height of the cell centre: the mean of the cell's eight node heights.

The node heights themselves, ghost nodes included, are written to ``Level_0/z_phys_nd`` inside the plotfile as a
VisMF MultiFab.

Turbulence diagnostics plotfile
-------------------------------

With ``io.turb_plot_file`` set and the SIMPLEC stage on, ``<io.turb_plot_file>NNNNN`` holds the closure's fields at
the final state under ERF's variable names, so ERF's Canonical_RANS check scripts read it unchanged
(``Source/TurbulenceOutput.H``): ``x_velocity``, ``y_velocity``, ``z_velocity`` (cell means of the face velocities),
``theta``, ``KE`` (k), ``Kmv`` (:math:`\nu_t`; the density is 1), ``Khv`` (:math:`\gamma_\theta` less the molecular
diffusivity), ``Lturb``, ``walldist``, ``diss`` (:math:`\varepsilon`), ``density`` (1), ``Rt``, ``cmu`` and
``cmu_prime`` (the stability functions of :math:`R_t`, AL01 Eqs. 31–32, as ERF writes them), ``ustar`` (the ground
cells; 0 above), ``z_phys`` (the cell-centre height), and the residual at the final state, ``R_u``, ``R_v``, ``R_w`` (cell means of the face
residuals), ``R_p``, ``R_theta``, ``R_tke``: where the iteration has not converged. Beside it ``surf_hist.dat`` holds ``t u_star`` and one row with the iteration count and the
mean :math:`u_*` over the ground, the file ERF's ``--physics`` checks read.

Checkpoint
----------

``<io.chk_file>NNNNN`` holds a text ``Header`` and every field with its ghost cells: two for the solution fields,
one for the node heights. (Up to sub-phase 1.11 the solution fields had one ghost cell; those checkpoints restart
unchanged, since only valid values are read back. The second ghost layer makes a checkpoint of 8³-cell boxes
about 1.5 times larger; larger boxes grow less.)

The ``Header`` lists the format version, the bytes per real (8 or 4), the step, the mesh (``n_cell``, and
``prob_lo`` and ``prob_hi`` in hexadecimal floating point so they read back exactly), ``is_periodic`` and the
BoxArray. The fields are ``u v w`` (faces), ``p theta tke`` (cell centres) and ``z_phys_nd`` (nodes).

A restart restores every valid value exactly, and the ghost nodes of ``z_phys_nd``, which come from the terrain.
Other ghost cells are refilled by the solver. With the SIMPLEC stage the checkpoint also holds the driver history
(``Driver``, :doc:`driver`) and the closure's :math:`\nu_t` history (``Level_0/closure_nut``, :doc:`closure`).
