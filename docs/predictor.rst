SIMPLEC momentum predictor
==========================

The momentum predictor (sub-phase 1.14, ``Source/Predictor.H``) advances the three velocity components by one
pseudo-time step with defect correction. In delta form,

.. math::
   :label: predictor

   \Bigl(\frac{1}{\Delta\tau} + A_1\Bigr)\,\delta u = -R(q), \qquad a_P = \frac{1}{\Delta\tau} + \mathrm{diag}(A_1),

with :math:`A_1` the first-order implicit-side operator of the component (:doc:`transport`: bounded upwind
convection and the principal diffusion, folded rows from the wall model, the damping layers, the Dirichlet identity
rows of the open faces and the ghost couplings), :math:`R` the full second-order residual at the state (:doc:`residual`)
and :math:`\Delta\tau` the local pseudo-time step on the component's faces. Everything is per unit volume (the plan's
:math:`V/\Delta\tau` and :math:`a_P` are these times the cell volume). The pressure gradient stays in :math:`R`: the
pressure correction that makes the predicted velocity divergence-free is the SIMPLEC stage's (:doc:`simplec`). For it
the predictor stores
on the faces both :math:`a_P` and the SIMPLEC denominator :math:`d = a_P - \sum a_{nb}`, which is :math:`1/\Delta\tau`
plus whatever the folds added to the diagonal (the wall drag, the damping, a value fold) and cannot be recovered from
:math:`a_P` alone; :math:`d = 0` on the Dirichlet rows of the open faces, where the velocity is not corrected.

**Why defect correction.** :math:`A_1` is first order and diagonally dominant (an M-matrix with the pseudo-time
term), so it can be solved cheaply by line sweeps, while the fixed point is the second-order residual's: a converged
state has :math:`R = 0` and the predictor leaves it unchanged, whatever :math:`A_1` is.

Pseudo-time step on the faces
-----------------------------

The driver's step (:doc:`driver`) is cell-centred; a face takes the arithmetic mean of its two cells (``face_dtau``),
the one cell on a non-periodic domain face, the ground and the top. (The volume-weighted mean of :math:`1/\Delta\tau`
would be the control-volume-consistent choice; for the layer ratios in use the two differ by 0.1–0.2 %, and the step
is not a stability bound for the M-matrix solve, so the simpler mean is kept.) **Relation to under-relaxation.** The classic SIMPLE
relaxation solves :math:`(a_P/\alpha_u)\,\delta u - \sum a_{nb}\,\delta u_{nb} = -R`. With :math:`1/\Delta\tau =
a_P (1 - \alpha_u)/\alpha_u`, i.e. :math:`\Delta\tau = \alpha_u/(1 - \alpha_u)/\mathrm{diag}(A_1)`, :eq:`predictor`
is the same system, so the pseudo-time form contains the relaxed form and sets the step from a CFL instead of from a
relaxation factor (which near the ground, where :math:`\mathrm{diag}(A_1)` is large, would hold the effective
horizontal CFL to a few units). The tests show the two forms agree to 1e-12.

z-line Jacobi sweeps
--------------------

Every column :math:`(i, j)` is solved exactly for its tridiagonal part :math:`(-a_{z-}, a_P, -a_{z+})` with the
Thomas algorithm, the :math:`x` and :math:`y` neighbours taken from the previous sweep (zero at the first); the
columns are independent, one thread per column on a GPU (``zline_jacobi``). The vertical upwind weights of
:math:`A_1` make the vertical advection implicit along with the vertical diffusion (:doc:`transport`), so the step is
not bounded by the stiff vertical terms. The sweeps converge to the solution of :eq:`predictor` for any
:math:`\Delta\tau` (the system is an M-matrix), and the error falls with every sweep; the number per predictor step
is ``steady.zline_sweeps`` (default 2).

Line LU-SGS sweeps
------------------

``steady.sweep = lusgs`` (sub-phase 2.1; ``zjacobi`` is the default) keeps the exact z-line solve but visits the
columns of a box in order: a forward pass over :math:`(i, j)` in which the :math:`x` and :math:`y` neighbours below
the column in the ordering are the pass's own updated values, then a backward pass over the columns in reverse
with the updated upper neighbours; a symmetric line Gauss–Seidel in delta form, in place (``zline_lusgs``). One
"sweep" of ``steady.zline_sweeps`` is the pair of passes, and one is allowed (the lag of a single Jacobi sweep is
the reason that solver needs two). For a wind whose :math:`u` and :math:`v` have the same sign in the ordering a
pass resolves the first-order upwind coupling of :math:`A_1` exactly: in 1-D the system is triangular and one pass
solves it, so the defect correction becomes :math:`1 - A_1^{-1}A_2`, whose scalar pure-advection factor the scan
``tools/sweeps_cfl_scan.py --thresholds`` finds below 1 at every CFL (the odd–even mode's factor is
:math:`(2c - 1)/(2c + 1)`, tending to 1 as the CFL :math:`c` grows, so the gain per iteration flattens instead of
turning into growth). For :math:`u` and :math:`v` of opposite sign each pass leaves one upwind neighbour lagged
and the per-sweep factor is :math:`c_u c_v / ((1 + c_u)(1 + c_v))`: one sweep is then unstable beyond a CFL of about
2.5 and two sweeps beyond about 33 (the review's derivation), so the decks keep two sweeps and a cap of 16;
alternating the ordering between sweeps is for sub-phase 2.2. The iteration path depends on the box layout
(``mesh.max_grid_size``: the couplings across box faces are lagged, and the seam error decays downstream over
about :math:`c` cells, so boxes shorter than the CFL limit the solve), unlike Jacobi's; the converged state does
not. What the sweeps leave otherwise is the horizontal diffusion (at a cell Péclet number of 2 the per-sweep
factor is 0.57 at CFL 4 and 0.85 at 16; the outer iteration stays stable). Sub-phase 2.2 adds the box-local
handling, the overlap and the GPU diagonal ordering; here a box is one sequential thread on a GPU, so ``lusgs``
is a CPU solver until then. The stage reports the change of the last sweep's backward pass as ``last_sweep``,
which is not comparable with Jacobi's last-sweep change. The tests (UT-2.1.1–5): one and two symmetric sweeps equal a dense
block-triangular reference to 1e-12, and on two boxes the sweeps reach the dense solution with a rank-independent
path; one symmetric sweep is exact for pure upwind convection along the wind (either sign);
the sweeps converge to the dense solution of the real :math:`A_1` for a wind of either sign, faster than Jacobi;
:math:`k` stays at or above its floor; a zero right-hand side gives a zero correction. On ERF's 2-D hill two
symmetric sweeps at a CFL cap of 16 reach the same steady state as sixteen Jacobi sweeps at a cap of 4 in a third
fewer iterations and half the time (:doc:`simplec`).

**The flag from sub-phase 1.5.** The Fourier analysis of that review found that *one* lagged z-line sweep per outer
iteration is unstable with the second-order right-hand side at large horizontal CFL (an amplification factor of 1.09
at CFL 1, growing with the CFL), while symmetric Gauss–Seidel in :math:`x, y` is stable. That is why the default is
two sweeps (and one is not an input), why ``steady.cfl_max`` should stay moderate until the line LU-SGS of phase 2
(which adds the :math:`x, y` Gauss–Seidel ordering), and why the usable CFL of this stage is measured in 1.16 rather
than assumed. A Fourier scan of the outer iteration made for this review (first-order upwind :math:`A_1`, an unlimited
second-order right-hand side, pure advection, the vertical direction solved exactly) gives the worst-mode amplification
per outer iteration: with 2 sweeps 1.000 at CFL 0.8, 1.003 at 1, 1.05 at 2, 1.22 at 5; with 4 sweeps 1.000 up to CFL
4 and 1.001 at 5. Horizontal diffusion helps (at a cell Péclet number of 20, 2 sweeps hold to CFL 1.2 and 4 sweeps to
5). **Sub-phase 1.16 rescan** (``tools/sweeps_cfl_scan.py``: the same first-order upwind :math:`A_1`, the lagged
x neighbours zero at the first sweep, a 2-D staggered u–w–p model with the exact SIMPLEC projection, which reduces to
the scalar factor because the lagged solve commutes with the projection): the first CFL with an amplification above 1
is 0.9 with 2 sweeps, **2.4 with 4, 3.7 with 8 and 6.1 with 16** (second-order upwind residual; 2.3, 3.4 and 5.4
with central; ``--thresholds`` prints the table), the worst mode a short x wave (:math:`k_x\Delta x \approx
2.3`–2.9). The 1.15 figure of 4 for four sweeps
is not what this operator gives; the 1.16 thresholds are the ones the terrain runs confirm. Over ERF's 2-D hill at a
nominal cap of 2 the local CFL passes 2.4 near the crest (the speed-up, the thinner terrain-following cells) and the
run grows slowly after 2000–4000 iterations with four sweeps (an e-folding of about 400 iterations, in the w
residual first), converges with first-order upwind (no defect correction) and converges with 8 or 16 sweeps; a 1 cm
hill grows the same way, so any horizontal variation seeds it (:doc:`simplec`). The hill decks run 16 sweeps with
a cap of 4. The flat ABL deck keeps four sweeps and a cap of 2, below the threshold (at 4 the stratified
k–:math:`\theta`–:math:`w` system and the SER ramp fall into bursts); the driver's default cap of 10 is for phase 2.
A column with nothing varying horizontally has exact lagged neighbours and takes any cap (the 1-D column runs at 50).

Point-implicit buoyancy and rotation
------------------------------------

Two couplings of the residual are conservative in pseudo-time and undamped where :math:`\nu_t` vanishes: the
buoyancy–:math:`w` oscillator of a stable stratification (frequency :math:`N`) and the Coriolis rotation of
:math:`u, v` (frequency :math:`f`). An explicit pseudo-time step amplifies each by :math:`(1 + \omega^2\Delta\tau^2)^{1/2}`
per iteration whatever the CFL. With gravity on, the :math:`w` rows gain on the diagonal the stiffness
:math:`(g/\theta_0)\max(\partial_z\theta, 0)\,\tau_\theta` and on the right-hand side the buoyancy of the change the
:math:`\theta` equation predicts, :math:`-(g/\theta_0)\,\tau_\theta R_\theta` (the face mean of the two cells), with
:math:`\tau_\theta = 1/(1/\Delta\tau + a_{P,\theta})` the step :math:`\theta`'s own implicit rows allow (the face mean of
its diagonal). With Coriolis on, the :math:`u` and :math:`v` rows gain :math:`f^2\max(\tau_c, \tau_{\mathrm{own}})`
on the diagonal and on the right-hand side

.. math::

   -f\,\bigl(\tau_c\,\langle R_v\rangle + \tau_{\mathrm{own}}\,(\overline{R_v} - \langle R_v\rangle)\bigr),
   \qquad
   +f\,\bigl(\tau_c\,\langle R_u\rangle + \tau_{\mathrm{own}}\,(\overline{R_u} - \langle R_u\rangle)\bigr),

with :math:`\overline{R}` the four-face mean, :math:`\langle R\rangle` the horizontal mean of the other component's
residual over the level (one reduction per level and iteration), :math:`\tau_c = 1/\max(|f|/r, a_P - \sum a_{nb})`,
:math:`r` = ``simplec.rotation_relax`` (0.3), and :math:`\tau_{\mathrm{own}} = 1/(1/\Delta\tau + \max(a_P - \sum
a_{nb}, 0))` the row's own step. The inertial oscillation is horizontally uniform (:math:`u, v` of the whole level
rotate together), so the inertial time belongs to the mean alone; the fluctuation about it is the hill's wake, a
gravity wave or a roll, which the plain pseudo-time step treats at its own rate. :math:`\tau_c` is a fraction
:math:`r` of an inertial time where the drag folds are weak (a backward Euler step of the rotation alone damps it by
only :math:`(f\Delta\tau)^2/2` per iteration; with :math:`r/|f|` by about :math:`r f\Delta\tau`), the inverse of the
net diagonal where a drag dominates. Sub-phase 1.15 applied :math:`\tau_c` to the whole four-face mean; that
injected the other component's local residual at up to :math:`1/f` = 10 000 s and fed an internal gravity wave of an
inviscid inversion layer (growth 1.4 per iteration at :math:`r = 1`), and over ERF's 2-D hill grew at :math:`r =
0.3` as well. With the split both :math:`r = 0.3` and :math:`r = 1` converge on the flat case and the hill
(:doc:`simplec`). The diagonal carries the larger of the two steps so that the mean and the fluctuation are each a
backward Euler step of the rotation, a contraction for every :math:`f\Delta\tau` (the 2 × 2 analysis of the 1.16
review: with :math:`f^2\tau_c` alone the fluctuation's step is a contraction only for :math:`f\Delta\tau < r +
\sqrt{1 + r^2}`, 1.3 at :math:`r = 0.3`; the runs sit at :math:`f\Delta\tau` = 0.002–0.02). Both additions enter the
SIMPLEC denominator :math:`d` too,
vanish at a fixed point (the delta form and :math:`R = 0`), and are pinned by UT-1.15.7; why each is needed is in
:doc:`simplec`. The :math:`\cos\varphi` terms of the 3-D Coriolis force stay explicit.

What the predictor does with the boundaries
-------------------------------------------

The Dirichlet rows of the normal velocities on inflow and outflow faces (:doc:`boundary`: :math:`R = q - V`, identity
row in :math:`A_1`) get no :math:`1/\Delta\tau`: the step sets the boundary value exactly. **Why.** The identity row
is dimensionless while :math:`1/\Delta\tau` is a rate; adding them would remove a fraction :math:`1/(1 +
1/\Delta\tau)` of the boundary defect per step that depends on the time unit, not on the physics. The boundary
conditions supply the rows (``dirichlet_rows``); :math:`a_P = 1` there, which the pressure correction of 1.15 must
read as "no correction on this face". The ghost couplings of the first interior rows and the wall model's rows are
already folded into the coefficients. Outside the domain :math:`\delta u = 0`, and the two copies of a periodic seam
face are kept identical (``OverrideSync``). The caller refills the state's ghost cells
(``BoundaryConditions::fill_ghosts``) before the next residual.

Tests
-----

UT-1.14.1 to UT-1.14.5 in ``Tests/Unit/test_Predictor.cpp`` (:doc:`testing`): :math:`a_P` equals
:math:`1/\Delta\tau` plus the diagonal of ``line_jacobian`` exactly; the z-line solve and the Jacobi iteration
converge to a dense direct solution of the same system; the relaxed and the pseudo-time forms agree; a converged
state is a fixed point. The dense reference is assembled independently in the test from the coefficient fields.

Inputs
------

``steady.zline_sweeps`` (:doc:`inputs`). The SIMPLEC stage (:doc:`simplec`) combines it with the pressure correction
and the scalar sweeps.
