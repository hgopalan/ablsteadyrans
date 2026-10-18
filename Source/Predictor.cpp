#include "Predictor.H"

#include "Forcing.H"
#include "Stencil_K.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

#include <vector>

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

namespace {

//! The checks shared by the line solvers; defines delta if empty and zeroes it.
void check_line_system (const char* op, Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                        int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
    if (sweeps < 1) { amrex::Abort(std::string("ablsteady: ") + op + ": sweeps must be >= 1 (got " + std::to_string(sweeps) + ")"); }
    check_operator_field(op, "A", A, loc, mesh, 0);
    {
        // Every box must span the column, or the line solve would be a block
        // Jacobi in z (make_box_array never splits z; Mesh::define would take it).
        const amrex::Box dom = amrex::convert(geom.Domain(), index_type(loc));
        for (int b = 0; b < A.boxArray().size(); ++b) {
            if (A.boxArray()[b].smallEnd(2) != dom.smallEnd(2) || A.boxArray()[b].bigEnd(2) != dom.bigEnd(2)) {
                amrex::Abort(std::string("ablsteady: ") + op + ": every box must span the whole column (mesh.max_grid_size is horizontal)");
            }
        }
    }
    if (A.nComp() != n_implicit_coefs) { amrex::Abort(std::string("ablsteady: ") + op + ": A needs the implicit operator's components"); }
    check_operator_field(op, "ap", ap, loc, mesh, 0);
    check_operator_field(op, "rhs", rhs, loc, mesh, 0);
    if (delta.empty()) { delta.define(A.boxArray(), A.DistributionMap(), 1, 1); }
    check_operator_field(op, "delta", delta, loc, mesh, 1);
    delta.setVal(0.0);   // the ghosts outside the domain stay 0: the couplings are folded
}

//! The exact z-line solve of column (i, j): the tridiagonal (-a_zm, dg,
//! -a_zp) with the x/y neighbour terms of d on the right-hand side, the
//! Thomas recurrences in c, e, the result into out. The in-place sweep
//! passes the same mutable view as d and out (Array4's pointer is restrict-
//! qualified, so two views of one array must not be mixed in a kernel); the
//! column's own values are not read before they are written.
template <class SourceView>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve_column (amrex::Array4<amrex::Real const> const& a, amrex::Array4<amrex::Real const> const& dg,
                   amrex::Array4<amrex::Real const> const& r, SourceView const& d,
                   amrex::Array4<amrex::Real> const& c, amrex::Array4<amrex::Real> const& e,
                   amrex::Array4<amrex::Real> const& out, int i, int j, int klo, int khi) noexcept
{
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
    out(i, j, khi) = e(i, j, khi);
    for (int k = khi - 1; k >= klo; --k) { out(i, j, k) = e(i, j, k) - c(i, j, k) * out(i, j, k + 1); }
}

} // namespace

const char* sweep_name (SweepKind kind) noexcept { return kind == SweepKind::lusgs ? "lusgs" : "zjacobi"; }

amrex::Real line_sweeps (SweepKind kind, Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap,
                         const amrex::MultiFab& rhs, int sweeps, const Mesh& mesh, const amrex::Geometry& geom,
                         amrex::MultiFab& delta)
{
    return (kind == SweepKind::lusgs) ? zline_lusgs(loc, A, ap, rhs, sweeps, mesh, geom, delta)
                                      : zline_jacobi(loc, A, ap, rhs, sweeps, mesh, geom, delta);
}

amrex::Real zline_lusgs (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                         int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
    check_line_system("zline_lusgs", loc, A, ap, rhs, sweeps, mesh, geom, delta);
    amrex::MultiFab cprime(A.boxArray(), A.DistributionMap(), 1, 0), dprime(A.boxArray(), A.DistributionMap(), 1, 0),
                    before(A.boxArray(), A.DistributionMap(), 1, 0);
    amrex::Real change = 0;
    for (int s = 0; s < sweeps; ++s) {
        for (const bool forward : {true, false}) {
            // The reported change is the last sweep's backward pass alone
            // (with one sweep the whole sweep's change is |delta| itself).
            if (s == sweeps - 1 && !forward) { amrex::MultiFab::Copy(before, delta, 0, 0, 1, 0); }
            // The neighbours across box faces and periodic seams: the
            // previous pass's values (a block Jacobi across boxes; the two
            // copies of a periodic seam face differ within a pass and
            // OverrideSync keeps the owner's at the end).
            delta.FillBoundary(geom.periodicity());
            for (amrex::MFIter mfi(delta, false); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.validbox();
                const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0), jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
                const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
                const auto a = A.const_array(mfi);
                const auto dg = ap.const_array(mfi);
                const auto r = rhs.const_array(mfi);
                const auto d = delta.array(mfi);   // one view for the reads and the writes
                const auto c = cprime.array(mfi);
                const auto e = dprime.array(mfi);
                // The columns in order: one sequential thread per box (the
                // GPU diagonal ordering is sub-phase 2.2).
                amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
                {
                    if (forward) {
                        for (int j = jlo; j <= jhi; ++j) {
                            for (int i = ilo; i <= ihi; ++i) { solve_column(a, dg, r, d, c, e, d, i, j, klo, khi); }
                        }
                    } else {
                        for (int j = jhi; j >= jlo; --j) {
                            for (int i = ihi; i >= ilo; --i) { solve_column(a, dg, r, d, c, e, d, i, j, klo, khi); }
                        }
                    }
                });
            }
        }
        if (s == sweeps - 1) {
            amrex::MultiFab::Subtract(before, delta, 0, 0, 1, 0);
            change = before.norm0();
        }
    }
    delta.OverrideSync(geom.periodicity());
    return change;
}

amrex::Real zline_jacobi (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                          int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
    check_line_system("zline_jacobi", loc, A, ap, rhs, sweeps, mesh, geom, delta);
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
                // The lagged x/y neighbours on the right-hand side.
                solve_column(a, dg, r, d, c, e, n, i, j, klo, khi);
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
                                            PredictorFields& ap, const SurfaceState* surface, const ClosureFields* fields,
                                            amrex::Real rotation_relax, SweepKind sweep)
{
    const char* op = "momentum_predictor";
    if (std::isnan(rotation_relax) || !(rotation_relax > 0 && rotation_relax <= 1)) {
        amrex::Abort("ablsteady: momentum_predictor: rotation_relax must be > 0 and <= 1 (simplec.rotation_relax)");
    }
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
    if (in.bc == nullptr && (in.buoyancy.gravity > 0 || (in.coriolis.enabled && coriolis_parameter(in.coriolis) != 0))) {
        amrex::Abort("ablsteady: momentum_predictor: the point-implicit buoyancy and rotation need the boundary conditions "
                     "(the Dirichlet rows they skip); give in.bc, or gravity 0 and Coriolis off");
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
        if (c < 2 && in.coriolis.enabled && coriolis_parameter(in.coriolis) != 0) {
            // Point-implicit rotation (the header): the diagonal gains
            // f^2 max(tau_c, tau_own) and the right-hand side the Coriolis
            // force of the change the other component's residual predicts,
            // its horizontal mean per level at the inertial step tau_c and
            // the fluctuation at the row's own step tau_own: -f (...) on u,
            // +f (...) on v.
            const amrex::Real fs = coriolis_parameter(in.coriolis);
            const amrex::MultiFab& other = (c == 0) ? R.v : R.u;
            amrex::MultiFab og(other.boxArray(), other.DistributionMap(), 1, 1);
            amrex::MultiFab::Copy(og, other, 0, 0, 1, 0);
            og.FillBoundary(geom.periodicity());
            // The horizontal mean of the other residual per level: the
            // inertial oscillation is horizontally uniform, so the inertial
            // time step applies to the mean alone and the fluctuation gets
            // the row's own step (the header).
            const amrex::Box& dom = geom.Domain();
            const int klo = dom.smallEnd(2), nz = dom.length(2);   // other is cell-centred in z
            amrex::Gpu::DeviceVector<amrex::Real> dsum(static_cast<std::size_t>(nz), amrex::Real(0));
            amrex::Real* sum_ptr = dsum.data();
            // Every face once: a box's high face plane in its nodal direction
            // is the next box's low plane (or, on a periodic domain, the
            // image of the domain's low plane), kept only on a non-periodic
            // domain edge; the count per level follows the same rule.
            const int d_o = 1 - c;
            const bool periodic_o = geom.isPeriodic(d_o);
            for (amrex::MFIter mfi(other, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box bx = mfi.tilebox();
                const amrex::Box& vb = mfi.validbox();
                if (bx.bigEnd(d_o) == vb.bigEnd(d_o) && (periodic_o || vb.bigEnd(d_o) != dom.bigEnd(d_o) + 1)) {
                    bx.growHi(d_o, -1);
                }
                const auto ro = other.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::Gpu::Atomic::AddNoRet(sum_ptr + (k - klo), ro(i, j, k));
                });
            }
            amrex::Gpu::streamSynchronize();
            std::vector<amrex::Real> hsum(static_cast<std::size_t>(nz));
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, dsum.begin(), dsum.end(), hsum.begin());
            amrex::ParallelDescriptor::ReduceRealSum(hsum.data(), nz);
            const amrex::Real npts_level = static_cast<amrex::Real>(dom.length(d_o) + (periodic_o ? 0 : 1))
                                           * static_cast<amrex::Real>(dom.length(c));
            for (amrex::Real& v : hsum) { v /= npts_level; }
            amrex::Gpu::DeviceVector<amrex::Real> dmean(static_cast<std::size_t>(nz));
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, hsum.begin(), hsum.end(), dmean.begin());
            const amrex::Real* mean_ptr = dmean.data();
            // The four faces around a row are inside the domain except for a
            // row on a non-periodic domain face that is not Dirichlet (an
            // outflow's even mirror): those rows get no coupling.
            const int d = c, lo_edge = dom.smallEnd(d), hi_edge = dom.bigEnd(d) + 1;
            const bool periodic = geom.isPeriodic(d);
            for (amrex::MFIter mfi(*aps[c], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto p = aps[c]->array(mfi);
                const auto dd = ds[c]->array(mfi);
                const auto r = rhs.array(mfi);
                const auto a = A.const_array(mfi);
                const auto dt = dtau_faces[c].const_array(mfi);
                const auto ro = og.const_array(mfi);
                const auto m = has_mask ? mask.const_array(mfi) : A.const_array(mfi);
                const bool is_u = (c == 0);
                const amrex::Real relax = rotation_relax;
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const amrex::Real keep = has_mask ? amrex::Real(1) - m(i, j, k) : amrex::Real(1);
                    if (keep == 0) { return; }
                    const int n = is_u ? i : j;
                    if (!periodic && (n == lo_edge || n == hi_edge)) { return; }
                    // The rotation's step: rotation_relax inertial times, or
                    // the inverse of the net diagonal where a drag dominates
                    // (the header).
                    const amrex::Real nb = a(i, j, k, coef_XM) + a(i, j, k, coef_XP) + a(i, j, k, coef_YM) + a(i, j, k, coef_YP)
                                           + a(i, j, k, coef_ZM) + a(i, j, k, coef_ZP);
                    const amrex::Real net = amrex::max(a(i, j, k, coef_P) - nb, amrex::Real(0));
                    const amrex::Real tau_c = amrex::Real(1) / amrex::max(std::abs(fs) / relax, net);
                    const amrex::Real tau_own = amrex::Real(1) / (amrex::Real(1) / dt(i, j, k) + net);
                    // The diagonal carries the larger of the two steps so that
                    // both the mean (tau_c) and the fluctuation (tau_own) are
                    // backward Euler steps of the rotation: contractions for
                    // every f dtau (the review of 1.16). In practice tau_c =
                    // r/|f| >> dtau >= tau_own and the max is tau_c.
                    const amrex::Real stiff = fs * fs * amrex::max(tau_c, tau_own);
                    p(i, j, k) += stiff;
                    dd(i, j, k) += stiff;
                    const amrex::Real avg = is_u ? stencil::interp<stencil::YF, stencil::XF>(ro, i, j, k)
                                                 : stencil::interp<stencil::XF, stencil::YF>(ro, i, j, k);
                    const amrex::Real mean = mean_ptr[k - klo];
                    r(i, j, k) += (is_u ? -fs : fs) * (tau_c * mean + tau_own * (avg - mean));
                });
            }
            amrex::Gpu::streamSynchronize();   // dmean is read by the kernels above
        }
        if (c == 2 && in.buoyancy.gravity > 0) {
            // Point-implicit buoyancy (the header): the w rows gain the
            // stiffness (g/theta_0) max(d theta/dz, 0) dtau on the diagonal
            // (into d as well) and the buoyancy of the theta change the theta
            // equation predicts, -(g/theta_0) dtau R_theta, on the right-hand
            // side; nothing on the Dirichlet rows.
            check_operator_field(op, "R_theta", R.theta, Loc::cell, mesh, 0);
            check_operator_field(op, "theta", q.theta, Loc::cell, mesh, 1);
            const amrex::Real gt = in.buoyancy.gravity / in.buoyancy.theta0;
            // The theta change its own step predicts, -R_theta / (1/dtau +
            // a_P,theta): the diagonal of theta's implicit rows bounds the
            // change where theta's operator is stiff (the diffusion of a
            // kink), where dtau R_theta alone would be far too large.
            amrex::MultiFab At;
            implicit_rows(F, q.z_phys_nd, in, Equation::theta, mesh, geom, At, surface, fields);
            for (amrex::MFIter mfi(*aps[c], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto p = aps[c]->array(mfi);
                const auto dd = ds[c]->array(mfi);
                const auto r = rhs.array(mfi);
                const auto dt = dtau_faces[c].const_array(mfi);
                const auto th = q.theta.const_array(mfi);
                const auto rt = R.theta.const_array(mfi);
                const auto at = At.const_array(mfi);
                const auto z = q.z_phys_nd.const_array(mfi);
                const auto m = has_mask ? mask.const_array(mfi) : A.const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const amrex::Real keep = has_mask ? amrex::Real(1) - m(i, j, k) : amrex::Real(1);
                    if (keep == 0) { return; }
                    // Cell-centre heights of the two cells (the mean of eight nodes).
                    auto zc = [&] (int kk) {
                        return amrex::Real(0.125) * (z(i, j, kk) + z(i + 1, j, kk) + z(i, j + 1, kk) + z(i + 1, j + 1, kk)
                                                     + z(i, j, kk + 1) + z(i + 1, j, kk + 1) + z(i, j + 1, kk + 1)
                                                     + z(i + 1, j + 1, kk + 1));
                    };
                    const amrex::Real gamma = (th(i, j, k) - th(i, j, k - 1)) / (zc(k) - zc(k - 1));
                    // theta's effective step at the face: 1 / (1/dtau + a_P,theta), a_P,theta the face mean.
                    const amrex::Real ap_theta = amrex::Real(0.5) * (at(i, j, k - 1, coef_P) + at(i, j, k, coef_P));
                    const amrex::Real tau = amrex::Real(1) / (amrex::Real(1) / dt(i, j, k) + ap_theta);
                    const amrex::Real stiff = gt * amrex::max(gamma, amrex::Real(0)) * tau;
                    p(i, j, k) += stiff;
                    dd(i, j, k) += stiff;
                    r(i, j, k) -= gt * tau * amrex::Real(0.5) * (rt(i, j, k - 1) + rt(i, j, k));
                });
            }
        }
        amrex::MultiFab delta;
        *last[c] = line_sweeps(sweep, loc, A, *aps[c], rhs, sweeps, mesh, geom, delta);
        *dmax[c] = delta.norm0();
        amrex::MultiFab::Add(*comps[c], delta, 0, 0, 1, 0);
    }
    return out;
}

} // namespace ablsteady
