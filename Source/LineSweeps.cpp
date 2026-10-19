#include "LineSweeps.H"

#include "ImplicitOperator.H"
#include "LineSolve_K.H"

#include <AMReX_iMultiFab.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_GpuLaunch.H>

#include <algorithm>
#include <string>

namespace ablsteady {

namespace {

amrex::IntVect loc_type (Loc loc)
{
    switch (loc) {
    case Loc::xface: return amrex::IntVect(1, 0, 0);
    case Loc::yface: return amrex::IntVect(0, 1, 0);
    case Loc::zface: return amrex::IntVect(0, 0, 1);
    default: return amrex::IntVect(0);
    }
}

using CArr = amrex::Array4<amrex::Real const>;
using MArr = amrex::Array4<amrex::Real>;

// One pass over the columns of S (in place in d; c, e the Thomas scratch;
// a, dg, r and c, e must cover S, d must cover S grown by one in x and y).
// One kernel per ordering and direction: no option branch inside a device
// lambda.

void pass_sequential_forward (const amrex::Box& S, CArr const& a, CArr const& dg, CArr const& r, MArr const& d,
                              MArr const& c, MArr const& e)
{
    const int ilo = S.smallEnd(0), ihi = S.bigEnd(0), jlo = S.smallEnd(1), jhi = S.bigEnd(1);
    const int klo = S.smallEnd(2), khi = S.bigEnd(2);
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
    {
        for (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) { solve_column(a, dg, r, d, c, e, d, i, j, klo, khi); }
        }
    });
}

void pass_sequential_backward (const amrex::Box& S, CArr const& a, CArr const& dg, CArr const& r, MArr const& d,
                               MArr const& c, MArr const& e)
{
    const int ilo = S.smallEnd(0), ihi = S.bigEnd(0), jlo = S.smallEnd(1), jhi = S.bigEnd(1);
    const int klo = S.smallEnd(2), khi = S.bigEnd(2);
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
    {
        for (int j = jhi; j >= jlo; --j) {
            for (int i = ihi; i >= ilo; --i) { solve_column(a, dg, r, d, c, e, d, i, j, klo, khi); }
        }
    });
}

// The columns of the wavefront (i - ilo) + (j - jlo) = s of S, one per
// thread: m-th column at i = ilo + amin + m.
void plane (int s, const amrex::Box& S, CArr const& a, CArr const& dg, CArr const& r, MArr const& d, MArr const& c,
            MArr const& e)
{
    const int ilo = S.smallEnd(0), jlo = S.smallEnd(1), klo = S.smallEnd(2), khi = S.bigEnd(2);
    const int nx = S.length(0), ny = S.length(1);
    const int amin = std::max(0, s - (ny - 1)), amax = std::min(s, nx - 1);
    amrex::ParallelFor(amax - amin + 1, [=] AMREX_GPU_DEVICE (int m) noexcept
    {
        const int ia = amin + m;
        solve_column(a, dg, r, d, c, e, d, ilo + ia, jlo + s - ia, klo, khi);
    });
}

void pass_diagonal_forward (const amrex::Box& S, CArr const& a, CArr const& dg, CArr const& r, MArr const& d,
                            MArr const& c, MArr const& e)
{
    const int last = S.length(0) + S.length(1) - 2;
    for (int s = 0; s <= last; ++s) { plane(s, S, a, dg, r, d, c, e); }
}

void pass_diagonal_backward (const amrex::Box& S, CArr const& a, CArr const& dg, CArr const& r, MArr const& d,
                             MArr const& c, MArr const& e)
{
    const int last = S.length(0) + S.length(1) - 2;
    for (int s = last; s >= 0; --s) { plane(s, S, a, dg, r, d, c, e); }
}

using PassFn = void (*) (const amrex::Box&, CArr const&, CArr const&, CArr const&, MArr const&, MArr const&, MArr const&);

PassFn pass_of (LineOrder order, bool forward)
{
    if (order == LineOrder::diagonal) { return forward ? pass_diagonal_forward : pass_diagonal_backward; }
    return forward ? pass_sequential_forward : pass_sequential_backward;
}

void abort_on (const std::vector<std::string>& err)
{
    if (err.empty()) { return; }
    std::string all = "ablsteady: zline_lusgs: invalid options:";
    for (const auto& e : err) { all += "\n  " + e; }
    amrex::Abort(all);
}

// Sub-phase 2.1's path (steady.line_ghosts = lagged): the previous pass's
// values across box faces and periodic seams, every valid point swept.
amrex::Real lusgs_lagged (const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs, int sweeps,
                          LineOrder order, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
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
            // copies of a shared face differ within a pass and OverrideSync
            // keeps the owner's at the end).
            delta.FillBoundary(geom.periodicity());
            const PassFn pass = pass_of(order, forward);
            for (amrex::MFIter mfi(delta, false); mfi.isValid(); ++mfi) {
                pass(mfi.validbox(), A.const_array(mfi), ap.const_array(mfi), rhs.const_array(mfi), delta.array(mfi),
                     cprime.array(mfi), dprime.array(mfi));
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

// Sub-phase 2.2 (steady.line_ghosts = zero): box-local sweeps over the owned
// region grown by the overlap, delta = 0 outside it, no communication
// between the passes.
amrex::Real lusgs_zero (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                        int sweeps, LineOrder order, int overlap, const amrex::Geometry& geom, amrex::MultiFab& delta)
{
    const amrex::BoxArray& ba = A.boxArray();
    const amrex::DistributionMapping& dm = A.DistributionMap();
    // The system on the grown region: copies with overlap ghost cells and
    // the periodic images (outside a non-periodic domain face the ghosts are
    // 0 and never read: the swept region is clipped to the domain).
    amrex::MultiFab Ag, apg, rg;
    const amrex::MultiFab* As = &A;
    const amrex::MultiFab* aps = &ap;
    const amrex::MultiFab* rs = &rhs;
    if (overlap > 0) {
        const amrex::IntVect ng(overlap, overlap, 0);
        Ag.define(ba, dm, A.nComp(), ng);
        apg.define(ba, dm, 1, ng);
        rg.define(ba, dm, 1, ng);
        for (auto* p : {&Ag, &apg, &rg}) { p->setVal(0.0); }
        amrex::MultiFab::Copy(Ag, A, 0, 0, A.nComp(), 0);
        amrex::MultiFab::Copy(apg, ap, 0, 0, 1, 0);
        amrex::MultiFab::Copy(rg, rhs, 0, 0, 1, 0);
        for (auto* p : {&Ag, &apg, &rg}) { p->FillBoundary(geom.periodicity()); }
        As = &Ag;
        aps = &apg;
        rs = &rg;
    }
    // The box-local iterate: zero outside the swept region (the ghost delta).
    amrex::MultiFab work(ba, dm, 1, amrex::IntVect(overlap + 1, overlap + 1, 1));
    work.setVal(0.0);
    const amrex::IntVect ns(overlap, overlap, 0);
    amrex::MultiFab cprime(ba, dm, 1, ns), dprime(ba, dm, 1, ns), before(ba, dm, 1, 0);
    amrex::Real change = 0;
    for (int s = 0; s < sweeps; ++s) {
        for (const bool forward : {true, false}) {
            if (s == sweeps - 1 && !forward) { amrex::MultiFab::Copy(before, work, 0, 0, 1, 0); }
            const PassFn pass = pass_of(order, forward);
            for (amrex::MFIter mfi(work, false); mfi.isValid(); ++mfi) {
                const amrex::Box S = line_sweep_box(loc, line_owned_box(loc, mfi.validbox(), geom), overlap, geom);
                pass(S, As->const_array(mfi), aps->const_array(mfi), rs->const_array(mfi), work.array(mfi),
                     cprime.array(mfi), dprime.array(mfi));
            }
        }
    }
    // The last backward pass's change over the owned points (an overlap
    // point or a non-owned face copy is another box's unknown).
    amrex::MultiFab::Subtract(before, work, 0, 0, 1, 0);
    for (amrex::MFIter mfi(before, false); mfi.isValid(); ++mfi) {
        const amrex::Box owned = line_owned_box(loc, mfi.validbox(), geom);
        const auto b = before.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (!owned.contains(i, j, k)) { b(i, j, k) = 0; }
        });
    }
    change = before.norm0();
    // The owned values into delta; every non-owned copy of a face takes its
    // owner's value (the owner mask of line_owned_box, not AMReX's default).
    amrex::MultiFab::Copy(delta, work, 0, 0, 1, 0);
    if (loc == Loc::xface || loc == Loc::yface) {
        amrex::iMultiFab mask(ba, dm, 1, 0);
        for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
            const amrex::Box owned = line_owned_box(loc, mfi.validbox(), geom);
            const auto m = mask.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                m(i, j, k) = owned.contains(i, j, k) ? 1 : 0;
            });
        }
        amrex::OverrideSync(delta, mask, geom.periodicity());
    }
    return change;
}

} // namespace

const char* sweep_name (SweepKind kind) noexcept { return kind == SweepKind::lusgs ? "lusgs" : "zjacobi"; }
const char* line_ghosts_name (LineGhosts g) noexcept { return g == LineGhosts::zero ? "zero" : "lagged"; }
const char* line_order_name (LineOrder o) noexcept { return o == LineOrder::diagonal ? "diagonal" : "sequential"; }

std::vector<std::string> check_line_options (const LineSweepOptions& opt)
{
    std::vector<std::string> err;
    if (opt.overlap < 0 || opt.overlap > max_line_overlap) {
        err.emplace_back("steady.line_overlap = " + std::to_string(opt.overlap) + ": must be an integer >= 0 and <= "
                         + std::to_string(max_line_overlap));
    } else if (opt.overlap > 0 && opt.ghosts != LineGhosts::zero) {
        err.emplace_back("steady.line_overlap = " + std::to_string(opt.overlap)
                         + ": needs steady.line_ghosts = zero (the lagged ghosts sweep the valid box alone)");
    }
    return err;
}

std::string check_line_overlap_layout (const amrex::BoxArray& ba, int overlap)
{
    int narrowest = 0;
    for (int b = 0; b < ba.size(); ++b) {
        const amrex::Box bx = amrex::enclosedCells(ba[b]);
        const int w = std::min(bx.length(0), bx.length(1));
        narrowest = (b == 0) ? w : std::min(narrowest, w);
    }
    if (ba.size() > 0 && overlap > narrowest) {
        return "steady.line_overlap = " + std::to_string(overlap) + ": must be <= the narrowest box's width in x and y ("
               + std::to_string(narrowest) + " cells; mesh.max_grid_size)";
    }
    return {};
}

amrex::Box line_owned_box (Loc loc, const amrex::Box& vbox, const amrex::Geometry& geom)
{
    amrex::Box owned = vbox;
    const amrex::Box& dom = geom.Domain();
    for (int d = 0; d < 2; ++d) {
        if (!vbox.type(d)) { continue; }   // cell-centred in d: the box's own
        // the high face is the next box's low face, unless it is the domain's
        if (vbox.bigEnd(d) != dom.bigEnd(d) + 1) { owned.growHi(d, -1); }
        // on a periodic axis the domain's low face is the image of the high one
        if (geom.isPeriodic(d) && vbox.smallEnd(d) == dom.smallEnd(d)) { owned.growLo(d, -1); }
    }
    amrex::ignore_unused(loc);
    return owned;
}

amrex::Box line_sweep_box (Loc loc, const amrex::Box& owned, int overlap, const amrex::Geometry& geom)
{
    amrex::Box S = owned;
    const amrex::Box dom = amrex::convert(geom.Domain(), loc_type(loc));
    for (int d = 0; d < 2; ++d) {
        S.grow(d, overlap);
        if (!geom.isPeriodic(d)) {
            S.setSmall(d, std::max(S.smallEnd(d), dom.smallEnd(d)));
            S.setBig(d, std::min(S.bigEnd(d), dom.bigEnd(d)));
        }
    }
    return S;
}

void check_line_system (const char* op, Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap,
                        const amrex::MultiFab& rhs, int sweeps, const Mesh& mesh, const amrex::Geometry& geom,
                        amrex::MultiFab& delta)
{
    if (sweeps < 1) { amrex::Abort(std::string("ablsteady: ") + op + ": sweeps must be >= 1 (got " + std::to_string(sweeps) + ")"); }
    check_operator_field(op, "A", A, loc, mesh, 0);
    {
        // Every box must span the column, or the line solve would be a block
        // Jacobi in z (make_box_array never splits z; Mesh::define would take it).
        const amrex::Box dom = amrex::convert(geom.Domain(), loc_type(loc));
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

amrex::Real zline_lusgs (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                         int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta,
                         const LineSweepOptions& opt)
{
    check_line_system("zline_lusgs", loc, A, ap, rhs, sweeps, mesh, geom, delta);
    std::vector<std::string> err = check_line_options(opt);
    if (err.empty()) {
        const std::string lay = check_line_overlap_layout(A.boxArray(), opt.overlap);
        if (!lay.empty()) { err.push_back(lay); }
    }
    abort_on(err);
    if (opt.ghosts == LineGhosts::zero) { return lusgs_zero(loc, A, ap, rhs, sweeps, opt.order, opt.overlap, geom, delta); }
    return lusgs_lagged(A, ap, rhs, sweeps, opt.order, geom, delta);
}

} // namespace ablsteady
