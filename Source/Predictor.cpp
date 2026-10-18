#include "Predictor.H"

#include "BoundaryConditions.H"
#include "CVFluxes.H"
#include "ImplicitOperator.H"

#include <cmath>
#include <string>

namespace ablsteady {

namespace {

Loc loc_of (Equation eq)
{
    return eq == Equation::u ? Loc::xface : eq == Equation::v ? Loc::yface : eq == Equation::w ? Loc::zface : Loc::cell;
}

amrex::IntVect index_type (Loc loc)
{
    switch (loc) {
    case Loc::xface: return amrex::IntVect(1, 0, 0);
    case Loc::yface: return amrex::IntVect(0, 1, 0);
    case Loc::zface: return amrex::IntVect(0, 0, 1);
    default: return amrex::IntVect(0);
    }
}

} // namespace

void face_dtau (const amrex::MultiFab& dtau_cc, const Mesh& mesh, const amrex::Geometry& geom,
                amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces)
{
    const char* op = "face_dtau";
    check_operator_field(op, "dtau_cc", dtau_cc, Loc::cell, mesh, 0);
    {
        const amrex::Real lo = dtau_cc.min(0);   // collective: every rank aborts together
        if (std::isnan(lo) || !(lo > 0) || !std::isfinite(dtau_cc.max(0))) {
            amrex::Abort("ablsteady: face_dtau: dtau_cc must be finite and > 0 everywhere");
        }
    }
    const amrex::MultiFab* areas[3] = {&mesh.ax, &mesh.ay, &mesh.az};
    for (int d = 0; d < 3; ++d) {
        if (dtau_faces[d].empty()) { dtau_faces[d].define(areas[d]->boxArray(), mesh.detJ.DistributionMap(), 1, 0); }
    }
    check_operator_field(op, "dtau_x", dtau_faces[0], Loc::xface, mesh, 0);
    check_operator_field(op, "dtau_y", dtau_faces[1], Loc::yface, mesh, 0);
    check_operator_field(op, "dtau_z", dtau_faces[2], Loc::zface, mesh, 0);
    // One ghost cell with the periodic images; outside a non-periodic face
    // the ghost stays unread (the face takes its one cell).
    amrex::MultiFab dc(dtau_cc.boxArray(), dtau_cc.DistributionMap(), 1, 1);
    dc.setVal(0.0);
    amrex::MultiFab::Copy(dc, dtau_cc, 0, 0, 1, 0);
    dc.FillBoundary(geom.periodicity());
    const amrex::Dim3 dlo = amrex::lbound(geom.Domain()), dhi = amrex::ubound(geom.Domain());
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    for (amrex::MFIter mfi(dc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto d = dc.const_array(mfi);
        const auto fx = dtau_faces[0].array(mfi);
        const auto fy = dtau_faces[1].array(mfi);
        const auto fz = dtau_faces[2].array(mfi);
        amrex::ParallelFor(mfi.tilebox(amrex::IntVect(1, 0, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool lo = (!perx && i == dlo.x), hi = (!perx && i == dhi.x + 1);
            fx(i, j, k) = lo ? d(i, j, k) : hi ? d(i - 1, j, k) : amrex::Real(0.5) * (d(i - 1, j, k) + d(i, j, k));
        });
        amrex::ParallelFor(mfi.tilebox(amrex::IntVect(0, 1, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool lo = (!pery && j == dlo.y), hi = (!pery && j == dhi.y + 1);
            fy(i, j, k) = lo ? d(i, j, k) : hi ? d(i, j - 1, k) : amrex::Real(0.5) * (d(i, j - 1, k) + d(i, j, k));
        });
        amrex::ParallelFor(mfi.tilebox(amrex::IntVect(0, 0, 1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool lo = (k == dlo.z), hi = (k == dhi.z + 1);
            fz(i, j, k) = lo ? d(i, j, k) : hi ? d(i, j, k - 1) : amrex::Real(0.5) * (d(i, j, k - 1) + d(i, j, k));
        });
    }
}

amrex::Real zline_jacobi (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                          int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
    const char* op = "zline_jacobi";
    if (sweeps < 1) { amrex::Abort("ablsteady: zline_jacobi: sweeps must be >= 1 (got " + std::to_string(sweeps) + ")"); }
    check_operator_field(op, "A", A, loc, mesh, 0);
    {
        // Every box must span the column, or the line solve would be a block
        // Jacobi in z (make_box_array never splits z; Mesh::define would take it).
        const amrex::Box dom = amrex::convert(geom.Domain(), index_type(loc));
        for (int b = 0; b < A.boxArray().size(); ++b) {
            if (A.boxArray()[b].smallEnd(2) != dom.smallEnd(2) || A.boxArray()[b].bigEnd(2) != dom.bigEnd(2)) {
                amrex::Abort("ablsteady: zline_jacobi: every box must span the whole column (mesh.max_grid_size is horizontal)");
            }
        }
    }
    if (A.nComp() != n_implicit_coefs) { amrex::Abort("ablsteady: zline_jacobi: A needs the implicit operator's components"); }
    check_operator_field(op, "ap", ap, loc, mesh, 0);
    check_operator_field(op, "rhs", rhs, loc, mesh, 0);
    if (delta.empty()) { delta.define(A.boxArray(), A.DistributionMap(), 1, 1); }
    check_operator_field(op, "delta", delta, loc, mesh, 1);
    delta.setVal(0.0);   // the ghosts outside the domain stay 0: the couplings are folded
    // The Thomas recurrences and the new iterate, one column per thread.
    amrex::MultiFab cprime(A.boxArray(), A.DistributionMap(), 1, 0), dprime(A.boxArray(), A.DistributionMap(), 1, 0),
                    next(A.boxArray(), A.DistributionMap(), 1, 0);
    amrex::Real change = 0;
    for (int s = 0; s < sweeps; ++s) {
        delta.FillBoundary(geom.periodicity());
        for (amrex::MFIter mfi(next, false); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
            amrex::Box cols = bx;
            cols.setRange(2, klo, 1);
            const auto a = A.const_array(mfi);
            const auto dg = ap.const_array(mfi);
            const auto r = rhs.const_array(mfi);
            const auto d = delta.const_array(mfi);
            const auto c = cprime.array(mfi);
            const auto e = dprime.array(mfi);
            const auto n = next.array(mfi);
            amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                // Forward: with the lagged x/y neighbours on the right-hand side.
                for (int k = klo; k <= khi; ++k) {
                    const amrex::Real rk = r(i, j, k) + a(i, j, k, coef_XM) * d(i - 1, j, k) + a(i, j, k, coef_XP) * d(i + 1, j, k)
                                           + a(i, j, k, coef_YM) * d(i, j - 1, k) + a(i, j, k, coef_YP) * d(i, j + 1, k);
                    const amrex::Real lower = -a(i, j, k, coef_ZM), upper = -a(i, j, k, coef_ZP);
                    if (k == klo) {
                        const amrex::Real m = dg(i, j, k);
                        c(i, j, k) = upper / m;
                        e(i, j, k) = rk / m;
                    } else {
                        const amrex::Real m = dg(i, j, k) - lower * c(i, j, k - 1);
                        c(i, j, k) = upper / m;
                        e(i, j, k) = (rk - lower * e(i, j, k - 1)) / m;
                    }
                }
                // Backward.
                n(i, j, khi) = e(i, j, khi);
                for (int k = khi - 1; k >= klo; --k) { n(i, j, k) = e(i, j, k) - c(i, j, k) * n(i, j, k + 1); }
            });
        }
        // The last sweep's change (one reduction), then the new iterate into delta.
        if (s == sweeps - 1) {
            amrex::MultiFab::Copy(cprime, next, 0, 0, 1, 0);   // cprime is scratch now
            amrex::MultiFab::Subtract(cprime, delta, 0, 0, 1, 0);
            change = cprime.norm0();
        }
        amrex::MultiFab::Copy(delta, next, 0, 0, 1, 0);
    }
    // The two copies of a periodic seam face are the same unknown: keep them
    // bit-identical whatever the per-box arithmetic did.
    delta.OverrideSync(geom.periodicity());
    return change;
}

MomentumPredictorResult momentum_predictor (State& q, const ResidualFields& R,
                                            const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                                            const ResidualInputs& in,
                                            const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces,
                                            const Mesh& mesh, const amrex::Geometry& geom, int sweeps,
                                            PredictorFields& ap, const SurfaceState* surface, const ClosureFields* fields)
{
    const char* op = "momentum_predictor";
    if (sweeps < 1) { amrex::Abort("ablsteady: momentum_predictor: sweeps must be >= 1 (got " + std::to_string(sweeps) + ")"); }
    check_operator_field(op, "dtau_x", dtau_faces[0], Loc::xface, mesh, 0);
    check_operator_field(op, "dtau_y", dtau_faces[1], Loc::yface, mesh, 0);
    check_operator_field(op, "dtau_z", dtau_faces[2], Loc::zface, mesh, 0);
    for (int d = 0; d < 3; ++d) {
        const amrex::Real lo = dtau_faces[d].min(0);
        if (std::isnan(lo) || !(lo > 0) || !std::isfinite(dtau_faces[d].max(0))) {
            amrex::Abort("ablsteady: momentum_predictor: dtau_faces must be finite and > 0 everywhere");
        }
    }
    MomentumPredictorResult out;
    const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
    amrex::MultiFab* comps[3] = {&q.u, &q.v, &q.w};
    const amrex::MultiFab* res[3] = {&R.u, &R.v, &R.w};
    amrex::MultiFab* aps[3] = {&ap.ap_u, &ap.ap_v, &ap.ap_w};
    amrex::MultiFab* ds[3] = {&ap.d_u, &ap.d_v, &ap.d_w};
    amrex::Real* dmax[3] = {&out.max_du, &out.max_dv, &out.max_dw};
    amrex::Real* last[3] = {&out.last_sweep_du, &out.last_sweep_dv, &out.last_sweep_dw};
    for (int c = 0; c < 3; ++c) {
        const Loc loc = loc_of(eqs[c]);
        check_operator_field(op, "q", *comps[c], loc, mesh, 0);
        check_operator_field(op, "R", *res[c], loc, mesh, 0);
        amrex::MultiFab A;
        implicit_rows(F, q.z_phys_nd, in, eqs[c], mesh, geom, A, surface, fields);
        if (aps[c]->empty()) { aps[c]->define(A.boxArray(), A.DistributionMap(), 1, 0); }
        if (ds[c]->empty()) { ds[c]->define(A.boxArray(), A.DistributionMap(), 1, 0); }
        check_operator_field(op, "ap", *aps[c], loc, mesh, 0);
        check_operator_field(op, "d", *ds[c], loc, mesh, 0);
        // a_P = 1/dtau + diag(A1), except on the Dirichlet rows of the open
        // faces (identity rows, R = q - V): there a_P = 1 and the step sets the
        // boundary value exactly (1/dtau has the units of a rate, the identity
        // row none; adding them would relax the value by a unit-dependent
        // fraction).
        amrex::MultiFab mask;
        if (in.bc != nullptr) { in.bc->dirichlet_rows(eqs[c], loc, mesh, mask); }
        const bool has_mask = !mask.empty();
        for (amrex::MFIter mfi(*aps[c], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto p = aps[c]->array(mfi);
            const auto dd = ds[c]->array(mfi);
            const auto a = A.const_array(mfi);
            const auto dt = dtau_faces[c].const_array(mfi);
            const auto m = has_mask ? mask.const_array(mfi) : A.const_array(mfi);   // read only with has_mask
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real keep = has_mask ? amrex::Real(1) - m(i, j, k) : amrex::Real(1);
                p(i, j, k) = keep / dt(i, j, k) + a(i, j, k, coef_P);
                const amrex::Real nb = a(i, j, k, coef_XM) + a(i, j, k, coef_XP) + a(i, j, k, coef_YM) + a(i, j, k, coef_YP)
                                       + a(i, j, k, coef_ZM) + a(i, j, k, coef_ZP);
                dd(i, j, k) = keep * (p(i, j, k) - nb);   // 0 on a Dirichlet row
            });
        }
        amrex::MultiFab rhs(A.boxArray(), A.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(rhs, *res[c], 0, 0, 1, 0);
        rhs.mult(amrex::Real(-1), 0, 1, 0);
        amrex::MultiFab delta;
        *last[c] = zline_jacobi(loc, A, *aps[c], rhs, sweeps, mesh, geom, delta);
        *dmax[c] = delta.norm0();
        amrex::MultiFab::Add(*comps[c], delta, 0, 0, 1, 0);
    }
    return out;
}

} // namespace ablsteady
