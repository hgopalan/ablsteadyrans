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
Other ghost cells are refilled by the solver.
