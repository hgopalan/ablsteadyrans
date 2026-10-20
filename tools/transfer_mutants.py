#!/usr/bin/env python3
"""
Sub-phase 2.4: the one-line mutants of the transfer operators
(Source/Transfer.cpp).

Run against Transfer.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/transfer_mutants.py BUILD_DIR [T1 T4 ...]
"""
import mutants

F = "Source/Transfer.cpp"
M = [
 # restriction
 ("T1", "cell means not detJ-weighted", F,
  "const amrex::Real w = jj(r * I + a, r * Jc + b, k);", "const amrex::Real w = 1;"),
 ("T2", "x/y-face flux without the fine area factor", F,
  "s += aa(i, j, k) * qq(i, j, k);", "s += qq(i, j, k);"),
 ("T3", "x/y-face flux not divided by r", F,
  "o(I, J, k) = s * rinv;", "o(I, J, k) = s;"),
 ("T4", "z-face flux divided by r instead of r^2", F,
  "const amrex::Real r2inv = amrex::Real(1) / amrex::Real(r * r);",
  "const amrex::Real r2inv = amrex::Real(1) / amrex::Real(r);"),
 ("T5", "coarse u not divided by the coarse area factor", F,
  "amrex::MultiFab::Divide(coarse.u, mc.ax, 0, 0, 1, 0);", ";"),
 ("T6", "coarse Omega not divided by the coarse area factor", F,
  "amrex::MultiFab::Divide(om_c, mc.az, 0, 0, 1, 0);", ";"),
 ("T7", "w restricted in place of Omega", F,
  "restrict_zfaces(mf.az, om_f, r, img);", "restrict_zfaces(mf.az, fine.w, r, img);"),
 ("T8", "no zero-gradient rows for the slope term's u", F,
  "zero_gradient_rows(u1, gc.Domain());", ";"),
 ("T9", "no periodic fill for the slope term's u", F,
  "u1.FillBoundary(gc.periodicity());", ";"),
 ("T10", "coarse z_phys_nd not set by restrict_state", F,
  "amrex::MultiFab::Copy(coarse.z_phys_nd, h.z_phys_nd(lc), 0, 0, 1, ngrow_nodal);", ";"),
 # interpolation
 ("T11", "cell-centred position shifted by half a fine cell", F,
  "const int m = Nodal ? (i - a.flo) : (2 * (i - a.flo) + 1 - a.r);",
  "const int m = Nodal ? (i - a.flo) : (2 * (i - a.flo) - a.r);"),
 ("T12", "no linear extrapolation at a non-periodic face (reads past the last centre)", F,
  "a.imax = N - 2;", "a.imax = N - 1;"),
 ("T13", "periodic direction extrapolates instead of reading the periodic image", F,
  "a.imin = -1;", "a.imin = 0;"),
 ("T14", "one-cell direction reads a second coarse cell", F,
  "a.step = 0;", "a.step = 1;"),
 ("T15", "x weights swapped on the upper y row", F,
  "+ ty * ((one - tx) * q(I, J1, k) + tx * q(I1, J1, k));",
  "+ ty * (tx * q(I, J1, k) + (one - tx) * q(I1, J1, k));"),
 ("T16", "w interpolated with the nodal x axis", F,
  "interpolate<stencil::ZF>(img(coarse.w), fine.w, cx, cy);",
  "interpolate<stencil::ZF>(img(coarse.w), fine.w, nx, cy);"),
 ("T17", "coarse image filled without the periodic images", F,
  "c.ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), gc.periodicity());",
  "c.ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), amrex::Periodicity::NonPeriodic());"),
 ("T18", "fine z_phys_nd not set by prolong_fields", F,
  "amrex::MultiFab::Copy(fine.z_phys_nd, h.z_phys_nd(lf), 0, 0, 1, ngrow_nodal);", ";"),
 ("T19", "fine-below-coarse level check dropped", F,
  "if (!(lf >= 0 && lf < lc && lc < h.num_levels())) {", "if (!(lf >= 0 && lc < h.num_levels())) {"),
 # prolong_state
 ("T20", "k not clipped", F,
  "if (low) { f(i, j, k) = tke_floor; }", "if (false) { f(i, j, k) = tke_floor; }"),
 ("T21", "clipped-cell count not summed over ranks", F,
  "amrex::ParallelDescriptor::ReduceLongSum(res.tke_clipped);", ";"),
 ("T22", "infinite tke_floor accepted", F,
  "if (!std::isfinite(tke_floor) || !(tke_floor >= 0)) {", "if (!(tke_floor >= 0)) {"),
 ("T23", "no boundary fill before the projection", F,
  "bc.fill_ghosts(fine, surface, buoyancy);", ";"),
 ("T24", "projection corrects the Dirichlet rows too (d = 1 + mask)", F,
  "dd[c]->mult(amrex::Real(-1), 0, 1, 0);", ";"),
 ("T25", "no projection off the Dirichlet rows (d = -mask)", F,
  "dd[c]->plus(amrex::Real(1), 0, 1, 0);", ";"),
 ("T26", "p' not restored after the projection", F,
  "amrex::MultiFab::Copy(fine.p, p_keep, 0, 0, 1, fine.p.nGrow());", ";"),
]

# Equivalent mutants.
#  T6: a_z = 1 on every mesh (make_areas), so the division by it is 1.
#  T9: the slope term interpolates u from the x-faces of its own cell
#      (stencil::interp<XF, ZF>), so it reads no horizontal ghost cell.
#  T24: every Dirichlet row is a domain face where p' has a Neumann
#      condition, so the projection's flux there is 0 whatever d is: exact
#      with gmres_mlmg and mlmg. gmres_fft's preconditioner averages the face
#      coefficients into the cells next to it (fill_beta_cell), so there the
#      result changes within the Poisson tolerance.
EQUIVALENT = {"T6", "T9", "T24"}



if __name__ == "__main__":
    mutants.run(M, "Transfer.*", EQUIVALENT, __doc__)
