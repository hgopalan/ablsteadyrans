// Sub-phase 2.2: the box-local line LU-SGS sweeps (LineSweeps.H).
// UT-2.2.1 (every C-grid face is swept by exactly one box, on 1, 2 and 4
// boxes, against a host replay of the box-local sweeps), UT-2.2.2 (the
// overlaps 0, 2 and 4 reach the same converged state; a wider overlap needs
// fewer outer iterations; an overlap that covers the domain reproduces the
// one-box sweep bit for bit), UT-2.2.3 (the converged state does not depend
// on the layout: 1, 4 and 16 boxes; the iterates do), UT-2.2.4 (the
// diagonal wavefront ordering equals the sequential sweep bit for bit),
// UT-2.2.6 (beyond the plan: the alternating ordering) and the refusals.
#include "ImplicitOperator.H"
#include "LineSweeps.H"
#include "Predictor.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ablsteady::LineGhosts;
using ablsteady::LineOrder;
using ablsteady::LineSweepOptions;
using ablsteady::Loc;
using ablsteady::SweepKind;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

const Loc kLocs[] = {Loc::cell, Loc::xface, Loc::yface, Loc::zface};

const char* loc_name (Loc loc)
{
    return loc == Loc::cell ? "cell" : loc == Loc::xface ? "xface" : loc == Loc::yface ? "yface" : "zface";
}

amrex::IntVect type_of (Loc loc)
{
    return loc == Loc::xface ? amrex::IntVect(1, 0, 0) : loc == Loc::yface ? amrex::IntVect(0, 1, 0)
         : loc == Loc::zface ? amrex::IntVect(0, 0, 1) : amrex::IntVect(0);
}

// A flat mesh of nx x ny x nz cells in boxes of at most mx x my columns.
PCase box_mesh (int nx, int ny, int nz, bool periodic, int mx, int my)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = 0;
    const auto f = periodic ? ablsteady::PoissonFace::periodic : ablsteady::PoissonFace::neumann;
    PCase c = ablsteady::ptest::make_case(nx, ny, nz, f, f, ms, 0, true);
    c.in.mesh.max_grid_size = {mx, my};
    c.make([ms] (double x, double y) { return ms.h(x, y); });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

LineSweepOptions lusgs (LineGhosts g, int overlap = 0, LineOrder order = LineOrder::sequential)
{
    LineSweepOptions o(SweepKind::lusgs);
    o.ghosts = g;
    o.overlap = overlap;
    o.order = order;
    return o;
}

// A random M-matrix system at loc: the six neighbour coefficients in [0.1,
// 1] with the x ones scaled by xscale (a strong coupling along x, as the
// upwind advection at a large CFL), the diagonal their sum plus extra, a
// random right-hand side; the copies of a shared face carry one system.
void random_system (const PCase& c, Loc loc, unsigned seed, double extra, double xscale, amrex::MultiFab& A,
                    amrex::MultiFab& ap, amrex::MultiFab& rhs)
{
    const amrex::BoxArray ba = amrex::convert(c.ba, type_of(loc));
    A.define(ba, c.dm, ablsteady::n_implicit_coefs, 0);
    ap.define(ba, c.dm, 1, 0);
    rhs.define(ba, c.dm, 1, 0);
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        // seeded by box index, so the system is the same on any number of ranks
        std::mt19937 gen(seed + 7919u * static_cast<unsigned>(mfi.index()));
        std::uniform_real_distribution<double> dist(0.1, 1.0), rd(-1.0, 1.0);
        amrex::FArrayBox fa(b, ablsteady::n_implicit_coefs, amrex::The_Pinned_Arena()), fp(b, 1, amrex::The_Pinned_Arena()),
                         fr(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            double sum = 0;
            for (int n = 1; n < ablsteady::n_implicit_coefs; ++n) {
                double v = dist(gen);
                if (n == ablsteady::coef_XM || n == ablsteady::coef_XP) { v *= xscale; }
                fa(iv, n) = static_cast<amrex::Real>(v);
                sum += v;
            }
            fa(iv, ablsteady::coef_P) = static_cast<amrex::Real>(sum);
            fp(iv) = static_cast<amrex::Real>(sum + extra);
            fr(iv) = static_cast<amrex::Real>(rd(gen));
        });
        A[mfi].copy<amrex::RunOn::Device>(fa, b);
        ap[mfi].copy<amrex::RunOn::Device>(fp, b);
        rhs[mfi].copy<amrex::RunOn::Device>(fr, b);
        amrex::Gpu::streamSynchronize();
    }
    A.OverrideSync(c.geom.periodicity());
    ap.OverrideSync(c.geom.periodicity());
    rhs.OverrideSync(c.geom.periodicity());
}

// rhs - N x on the valid points (N the full system: ap on the diagonal, -a_nb
// off it; periodic images through the geometry, nothing outside a
// non-periodic face or the column).
void residual_of (const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs, const amrex::MultiFab& x,
                  const amrex::Geometry& geom, amrex::MultiFab& res)
{
    amrex::MultiFab xg(x.boxArray(), x.DistributionMap(), 1, 1);
    xg.setVal(0.0);
    amrex::MultiFab::Copy(xg, x, 0, 0, 1, 0);
    xg.FillBoundary(geom.periodicity());
    res.define(x.boxArray(), x.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(res, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
        const auto a = A.const_array(mfi);
        const auto p = ap.const_array(mfi);
        const auto r = rhs.const_array(mfi);
        const auto d = xg.const_array(mfi);
        const auto o = res.array(mfi);
        amrex::ParallelFor(vb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real nb = a(i, j, k, ablsteady::coef_XM) * d(i - 1, j, k) + a(i, j, k, ablsteady::coef_XP) * d(i + 1, j, k)
                             + a(i, j, k, ablsteady::coef_YM) * d(i, j - 1, k) + a(i, j, k, ablsteady::coef_YP) * d(i, j + 1, k);
            if (k > klo) { nb += a(i, j, k, ablsteady::coef_ZM) * d(i, j, k - 1); }
            if (k < khi) { nb += a(i, j, k, ablsteady::coef_ZP) * d(i, j, k + 1); }
            o(i, j, k) = r(i, j, k) - (p(i, j, k) * d(i, j, k) - nb);
        });
    }
}

// The outer defect correction x <- x + B (rhs - N x), B the line solve with
// opt and sweeps, until |rhs - N x| <= rtol |rhs|. Returns the iterations.
int outer_solve (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                 const LineSweepOptions& opt, int sweeps, const PCase& c, double rtol, amrex::MultiFab& x)
{
    x.define(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
    x.setVal(0.0);
    const double r0 = rhs.norm0();
    for (int it = 0; it < 2000; ++it) {
        amrex::MultiFab res;
        residual_of(A, ap, rhs, x, c.geom, res);
        if (res.norm0() <= rtol * r0) { return it; }
        amrex::MultiFab delta;
        ablsteady::line_sweeps(opt, loc, A, ap, res, sweeps, c.mesh, c.geom, delta);
        amrex::MultiFab::Add(x, delta, 0, 0, 1, 0);
    }
    return -1;
}

amrex::MultiFab gather (const amrex::MultiFab& mf)
{
    amrex::BoxArray one(mf.boxArray().minimalBox());
    const amrex::DistributionMapping dm1(amrex::Vector<int>{0});
    amrex::MultiFab g(one, dm1, mf.nComp(), 0);
    g.ParallelCopy(mf, 0, 0, mf.nComp(), 0, 0);
    return g;
}

amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

double max_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

// Host replay (rank 0) of the box-local sweeps (ghosts = zero), written from
// the rule of the header, not from the code's helpers: a box owns the
// points of its valid box except, in a nodal x or y direction, its high face
// unless that is the domain's high face, and on a periodic axis the domain's
// low face (the image of the last box's high face); it sweeps the owned
// region grown by n_o (clipped to the domain on a non-periodic axis) with
// zero outside it, j outer and i inner, and keeps the owned points. Along a
// periodic axis one cell wide both neighbours are the point itself.
// Returns max |delta - replay| over every valid point of delta (a non-owned
// copy compared with its owner's value) and the largest |replay|.
std::pair<double, double> replay_error (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap,
                                        const amrex::MultiFab& rhs, const amrex::MultiFab& delta, int sweeps, int n_o,
                                        const PCase& c)
{
    const amrex::MultiFab gA = gather(A), gp = gather(ap), gr = gather(rhs), gd = gather(delta);
    double worst = 0, scale = 0;
    for (amrex::MFIter mfi(gd, false); mfi.isValid(); ++mfi) {   // rank 0
        const auto a = host_copy(gA, mfi), p = host_copy(gp, mfi), r = host_copy(gr, mfi), got = host_copy(gd, mfi);
        const amrex::Box gb = mfi.validbox();
        const amrex::Box& cdom = c.geom.Domain();
        const amrex::IntVect t = type_of(loc);
        auto wrap = [&] (amrex::IntVect iv) {
            for (int d = 0; d < 2; ++d) {
                if (!c.geom.isPeriodic(d)) { continue; }
                const int L = cdom.length(d), lo = cdom.smallEnd(d);
                iv[d] = lo + (((iv[d] - lo) % L) + L) % L;
            }
            return iv;
        };
        const bool selfx = c.geom.isPeriodic(0) && cdom.length(0) == 1, selfy = c.geom.isPeriodic(1) && cdom.length(1) == 1;
        amrex::FArrayBox out(gb, 1, amrex::The_Pinned_Arena());
        out.setVal<amrex::RunOn::Host>(0.0);
        const amrex::BoxArray ba = amrex::convert(c.ba, t);
        for (int b = 0; b < ba.size(); ++b) {
            const amrex::Box vb = ba[b];
            amrex::Box own = vb;
            for (int d = 0; d < 2; ++d) {
                if (t[d] == 0) { continue; }
                if (vb.bigEnd(d) != cdom.bigEnd(d) + 1) { own.setBig(d, vb.bigEnd(d) - 1); }
                if (c.geom.isPeriodic(d) && vb.smallEnd(d) == cdom.smallEnd(d)) { own.setSmall(d, vb.smallEnd(d) + 1); }
            }
            amrex::Box S = own;
            for (int d = 0; d < 2; ++d) {
                if (cdom.length(d) == 1) { continue; }   // no overlap along a one-cell axis (a 2-D deck)
                S.grow(d, n_o);
                if (!c.geom.isPeriodic(d)) {
                    S.setSmall(d, std::max(S.smallEnd(d), gb.smallEnd(d)));
                    S.setBig(d, std::min(S.bigEnd(d), gb.bigEnd(d)));
                }
            }
            amrex::Box Sg = S;
            Sg.grow(0, 1).grow(1, 1);
            amrex::FArrayBox W(Sg, 1, amrex::The_Pinned_Arena()), cc(S, 1, amrex::The_Pinned_Arena()), ee(S, 1, amrex::The_Pinned_Arena());
            W.setVal<amrex::RunOn::Host>(0.0);
            const int klo = S.smallEnd(2), khi = S.bigEnd(2);
            auto column = [&] (int i, int j) {
                for (int k = klo; k <= khi; ++k) {
                    const amrex::IntVect iv(i, j, k), w = wrap(iv);
                    // a periodic one-cell axis: both neighbours are this point (implicit)
                    amrex::Real pw = p(w), xn = 0, yn = 0;
                    if (selfx) { pw -= a(w, ablsteady::coef_XM) + a(w, ablsteady::coef_XP); }
                    else { xn = a(w, ablsteady::coef_XM) * W(amrex::IntVect(i - 1, j, k)) + a(w, ablsteady::coef_XP) * W(amrex::IntVect(i + 1, j, k)); }
                    if (selfy) { pw -= a(w, ablsteady::coef_YM) + a(w, ablsteady::coef_YP); }
                    else { yn = a(w, ablsteady::coef_YM) * W(amrex::IntVect(i, j - 1, k)) + a(w, ablsteady::coef_YP) * W(amrex::IntVect(i, j + 1, k)); }
                    const amrex::Real rk = r(w) + xn + yn;
                    const amrex::Real lower = -a(w, ablsteady::coef_ZM), upper = -a(w, ablsteady::coef_ZP);
                    if (k == klo) {
                        cc(iv) = upper / pw;
                        ee(iv) = rk / pw;
                    } else {
                        const amrex::Real m = pw - lower * cc(amrex::IntVect(i, j, k - 1));
                        cc(iv) = upper / m;
                        ee(iv) = (rk - lower * ee(amrex::IntVect(i, j, k - 1))) / m;
                    }
                }
                W(amrex::IntVect(i, j, khi)) = ee(amrex::IntVect(i, j, khi));
                for (int k = khi - 1; k >= klo; --k) {
                    W(amrex::IntVect(i, j, k)) = ee(amrex::IntVect(i, j, k)) - cc(amrex::IntVect(i, j, k)) * W(amrex::IntVect(i, j, k + 1));
                }
            };
            for (int s = 0; s < sweeps; ++s) {
                for (int j = S.smallEnd(1); j <= S.bigEnd(1); ++j) { for (int i = S.smallEnd(0); i <= S.bigEnd(0); ++i) { column(i, j); } }
                for (int j = S.bigEnd(1); j >= S.smallEnd(1); --j) { for (int i = S.bigEnd(0); i >= S.smallEnd(0); --i) { column(i, j); } }
            }
            amrex::LoopOnCpu(own, [&] (int i, int j, int k) { out(wrap(amrex::IntVect(i, j, k))) = W(amrex::IntVect(i, j, k)); });
        }
        amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            worst = std::max(worst, std::abs(double(got(iv)) - double(out(wrap(iv)))));
            scale = std::max(scale, std::abs(double(out(wrap(iv)))));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceRealMax(scale);
    return {worst, scale};
}

} // namespace

// UT-2.2.1. The owned regions (line_owned_box, which the sweep uses) cover
// every unknown of the domain exactly once, periodic images identified, for
// cells and the three face locations on 1, 2 and 4 boxes (and an uneven
// 7 x 6 split into 4 x 3 columns at most, and a 2-D 8 x 1 mesh in two
// boxes), periodic and not. The swept
// result equals the host replay of the box-local sweeps of the header's
// rule (one and two sweeps, n_o = 0 and 2): a face swept by two boxes, or
// a ghost delta taken from a neighbour, changes the result next to the box
// face. On one non-periodic box the box-local sweep is sub-phase 2.1's.
TEST (LineSweeps, EveryFaceIsOwnedAndSweptOnce)
{
    struct Layout { int nx, ny, mx, my, boxes; };
    const Layout layouts[] = {{8, 8, 8, 8, 1}, {8, 8, 4, 8, 2}, {8, 8, 4, 4, 4}, {7, 6, 4, 3, 4}, {8, 1, 4, 1, 2}};
    for (const bool periodic : {false, true}) {
        for (const Layout& L : layouts) {
            PCase c = box_mesh(L.nx, L.ny, 4, periodic, L.mx, L.my);
            ASSERT_EQ(c.ba.size(), L.boxes);
            for (const Loc loc : kLocs) {
                const amrex::IntVect t = type_of(loc);
                const amrex::Box dom = amrex::convert(c.geom.Domain(), t);
                // the number of owners of each unknown (a periodic image is the same unknown)
                std::vector<int> count(static_cast<std::size_t>(dom.numPts()), 0);
                const amrex::BoxArray ba = amrex::convert(c.ba, t);
                for (int b = 0; b < ba.size(); ++b) {
                    const amrex::Box own = ablsteady::line_owned_box(loc, ba[b], c.geom);
                    EXPECT_TRUE(ba[b].contains(own)) << loc_name(loc);
                    amrex::LoopOnCpu(own, [&] (int i, int j, int k) {
                        amrex::IntVect iv(i, j, k);
                        for (int d = 0; d < 2; ++d) {
                            if (periodic && t[d] == 1 && iv[d] == c.geom.Domain().bigEnd(d) + 1) { iv[d] = c.geom.Domain().smallEnd(d); }
                        }
                        ++count[static_cast<std::size_t>(dom.index(iv))];
                    });
                }
                int bad = 0;
                amrex::LoopOnCpu(dom, [&] (int i, int j, int k) {
                    for (int d = 0; d < 2; ++d) {   // the image of a seam face is not a separate unknown
                        if (periodic && t[d] == 1 && (d == 0 ? i : j) == c.geom.Domain().bigEnd(d) + 1) { return; }
                    }
                    if (count[static_cast<std::size_t>(dom.index(amrex::IntVect(i, j, k)))] != 1) { ++bad; }
                });
                // the images carry no owner of their own
                amrex::LoopOnCpu(dom, [&] (int i, int j, int k) {
                    for (int d = 0; d < 2; ++d) {
                        if (periodic && t[d] == 1 && (d == 0 ? i : j) == c.geom.Domain().bigEnd(d) + 1
                            && count[static_cast<std::size_t>(dom.index(amrex::IntVect(i, j, k)))] != 0) { ++bad; }
                    }
                });
                EXPECT_EQ(bad, 0) << loc_name(loc) << " periodic " << periodic << " boxes " << L.boxes;
                // the swept result against the replay
                amrex::MultiFab A, ap, rhs;
                random_system(c, loc, 11u + static_cast<unsigned>(L.boxes), 0.3, 3.0, A, ap, rhs);
                for (const int n_o : {0, 2}) {
                    if (!ablsteady::check_line_overlap_layout(c.ba, n_o).empty()) { continue; }
                    for (const int sweeps : {1, 2}) {
                        amrex::MultiFab delta;
                        ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, delta, lusgs(LineGhosts::zero, n_o));
                        const auto [err, scale] = replay_error(loc, A, ap, rhs, delta, sweeps, n_o, c);
                        EXPECT_GT(scale, 0.1);
                        EXPECT_LT(err, tol(1e-13, 1e-5) * scale) << loc_name(loc) << " periodic " << periodic << " boxes " << L.boxes
                                                                 << " n_o " << n_o << " sweeps " << sweeps;
                    }
                }
                if (!periodic && L.boxes == 1) {
                    // one non-periodic box: nothing lies across a box face, and
                    // the box owns every point: zero and lagged are one sweep
                    amrex::MultiFab dz, dl;
                    ablsteady::zline_lusgs(loc, A, ap, rhs, 2, c.mesh, c.geom, dz, lusgs(LineGhosts::zero));
                    ablsteady::zline_lusgs(loc, A, ap, rhs, 2, c.mesh, c.geom, dl, lusgs(LineGhosts::lagged));
                    EXPECT_EQ(max_diff(dz, dl), 0) << loc_name(loc);
                }
            }
        }
    }
}

// UT-2.2.2. On a 16 x 8 mesh in 4 x 2 boxes of 4 x 4 columns with a strong
// x coupling: the outer defect correction with the box-local sweeps (two per
// step) reaches |rhs - N x| <= 1e-12 |rhs| for n_o = 0, 2 and 4, the three
// states agree to 1e-10 (relative), and the overlap cuts the outer
// iterations (strictly from 0 to 2 and from 2 to 4). Then an overlap that
// covers the domain from every box (two boxes of four cells, n_o = 4,
// non-periodic) reproduces the one-box sweep bit for bit.
TEST (LineSweeps, OverlapsReachTheSameState)
{
    for (const bool periodic : {false, true}) {
        PCase c = box_mesh(16, 8, 4, periodic, 4, 4);
        ASSERT_EQ(c.ba.size(), 8);
        for (const Loc loc : {Loc::cell, Loc::xface}) {
            amrex::MultiFab A, ap, rhs;
            random_system(c, loc, 21u, 0.05, 4.0, A, ap, rhs);
            amrex::MultiFab x0;
            int prev = -1;
            for (const int n_o : {0, 2, 4}) {
                amrex::MultiFab x;
                const int its = outer_solve(loc, A, ap, rhs, lusgs(LineGhosts::zero, n_o), 2, c, tol(1e-12, 1e-5), x);
                ASSERT_GT(its, 0) << loc_name(loc) << " n_o " << n_o;
                amrex::Print() << "UT-2.2.2 " << loc_name(loc) << " periodic " << periodic << " n_o " << n_o << ": " << its
                               << " outer iterations\n";
                if (n_o == 0) {
                    x0.define(x.boxArray(), x.DistributionMap(), 1, 0);
                    amrex::MultiFab::Copy(x0, x, 0, 0, 1, 0);
                } else {
                    EXPECT_LT(max_diff(x, x0), tol(1e-10, 1e-4) * static_cast<double>(x0.norm0())) << loc_name(loc) << " n_o " << n_o;
                    EXPECT_LT(its, prev) << loc_name(loc) << " n_o " << n_o;
                }
                prev = its;
            }
        }
    }
    {
        PCase two = box_mesh(8, 4, 4, false, 4, 8), one = box_mesh(8, 4, 4, false, 8, 8);
        ASSERT_EQ(two.ba.size(), 2);
        ASSERT_EQ(one.ba.size(), 1);
        for (const Loc loc : {Loc::cell, Loc::yface, Loc::zface}) {   // owned regions of four in x
            amrex::MultiFab A1, ap1, r1, A2, ap2, r2;
            random_system(one, loc, 31u, 0.3, 3.0, A1, ap1, r1);
            const amrex::BoxArray ba2 = amrex::convert(two.ba, type_of(loc));
            A2.define(ba2, two.dm, ablsteady::n_implicit_coefs, 0);
            ap2.define(ba2, two.dm, 1, 0);
            r2.define(ba2, two.dm, 1, 0);
            A2.ParallelCopy(A1, 0, 0, ablsteady::n_implicit_coefs);
            ap2.ParallelCopy(ap1);
            r2.ParallelCopy(r1);
            for (const auto order : {LineOrder::sequential, LineOrder::diagonal}) {
                amrex::MultiFab d1, d2;
                ablsteady::zline_lusgs(loc, A1, ap1, r1, 3, one.mesh, one.geom, d1, lusgs(LineGhosts::zero, 0, order));
                ablsteady::zline_lusgs(loc, A2, ap2, r2, 3, two.mesh, two.geom, d2, lusgs(LineGhosts::zero, 4, order));
                amrex::MultiFab g2(d1.boxArray(), d1.DistributionMap(), 1, 0);
                g2.ParallelCopy(d2);
                EXPECT_EQ(max_diff(d1, g2), 0) << loc_name(loc);
                amrex::MultiFab d0;   // without the overlap the box face shows
                ablsteady::zline_lusgs(loc, A2, ap2, r2, 3, two.mesh, two.geom, d0, lusgs(LineGhosts::zero, 0, order));
                amrex::MultiFab g0(d1.boxArray(), d1.DistributionMap(), 1, 0);
                g0.ParallelCopy(d0);
                EXPECT_GT(max_diff(d1, g0), 1e-3 * static_cast<double>(d1.norm0())) << loc_name(loc);
            }
        }
    }
}

// UT-2.2.3. A 16 x 16 mesh in 1, 4 and 16 boxes (periodic and not): the
// box-local sweeps depend on the box shape, so one step differs between the
// layouts (pinned below: the test is not of bit identity), but the outer
// iteration reaches the same converged state on every layout, to 1e-10
// (relative), for n_o = 0 and 2; the lagged ghosts of sub-phase 2.1 too.
TEST (LineSweeps, ConvergedStateDoesNotDependOnTheLayout)
{
    for (const bool periodic : {false, true}) {
        PCase c1 = box_mesh(16, 16, 4, periodic, 16, 16);
        ASSERT_EQ(c1.ba.size(), 1);
        for (const Loc loc : {Loc::cell, Loc::xface, Loc::yface}) {
            amrex::MultiFab A1, ap1, r1;
            random_system(c1, loc, 41u, 0.1, 2.0, A1, ap1, r1);
            for (const LineSweepOptions& opt : {lusgs(LineGhosts::zero, 0), lusgs(LineGhosts::zero, 2), lusgs(LineGhosts::lagged)}) {
                amrex::MultiFab xref, step1;
                ASSERT_GT(outer_solve(loc, A1, ap1, r1, opt, 2, c1, tol(1e-12, 1e-5), xref), 0);
                ablsteady::line_sweeps(opt, loc, A1, ap1, r1, 2, c1.mesh, c1.geom, step1);
                for (const int m : {8, 4}) {
                    PCase cm = box_mesh(16, 16, 4, periodic, m, m);
                    ASSERT_EQ(cm.ba.size(), (16 / m) * (16 / m));
                    const amrex::BoxArray ba = amrex::convert(cm.ba, type_of(loc));
                    amrex::MultiFab A(ba, cm.dm, ablsteady::n_implicit_coefs, 0), ap(ba, cm.dm, 1, 0), rhs(ba, cm.dm, 1, 0);
                    A.ParallelCopy(A1, 0, 0, ablsteady::n_implicit_coefs);
                    ap.ParallelCopy(ap1);
                    rhs.ParallelCopy(r1);
                    amrex::MultiFab x, step;
                    const int its = outer_solve(loc, A, ap, rhs, opt, 2, cm, tol(1e-12, 1e-5), x);
                    ASSERT_GT(its, 0);
                    amrex::Print() << "UT-2.2.3 " << loc_name(loc) << " periodic " << periodic << " ghosts "
                                   << ablsteady::line_ghosts_name(opt.ghosts) << " n_o " << opt.overlap << " boxes "
                                   << cm.ba.size() << ": " << its << " outer iterations\n";
                    amrex::MultiFab g(xref.boxArray(), xref.DistributionMap(), 1, 0);
                    g.ParallelCopy(x);
                    EXPECT_LT(max_diff(g, xref), tol(1e-10, 1e-4) * static_cast<double>(xref.norm0())) << loc_name(loc) << " " << m;
                    ablsteady::line_sweeps(opt, loc, A, ap, rhs, 2, cm.mesh, cm.geom, step);
                    amrex::MultiFab gs(step1.boxArray(), step1.DistributionMap(), 1, 0);
                    gs.ParallelCopy(step);
                    EXPECT_GT(max_diff(gs, step1), 1e-4 * static_cast<double>(step1.norm0())) << loc_name(loc) << " " << m;
                }
            }
        }
    }
}

// UT-2.2.4. The diagonal wavefront ordering gives the sequential sweep's
// result bit for bit: every location, 1, 2 and 4 boxes (an uneven 7 x 6
// split too), periodic and not, one and three sweeps, the lagged ghosts and
// the zero ghosts with n_o = 0 and 2; the reported change too. The launch
// counter shows each ordering reaches its own kernels (one launch per box
// and pass against one per wavefront plane), so the test cannot pass by
// running the sequential kernels twice.
TEST (LineSweeps, DiagonalOrderingEqualsTheSequentialSweep)
{
    struct Layout { int nx, ny, mx, my; };
    const Layout layouts[] = {{8, 8, 8, 8}, {8, 8, 4, 8}, {8, 8, 4, 4}, {7, 6, 4, 3}};
    for (const bool periodic : {false, true}) {
        for (const Layout& L : layouts) {
            PCase c = box_mesh(L.nx, L.ny, 4, periodic, L.mx, L.my);
            for (const Loc loc : kLocs) {
                amrex::MultiFab A, ap, rhs;
                random_system(c, loc, 51u, 0.2, 3.0, A, ap, rhs);
                for (const LineSweepOptions& base : {lusgs(LineGhosts::lagged), lusgs(LineGhosts::zero, 0), lusgs(LineGhosts::zero, 2)}) {
                    for (const int sweeps : {1, 3}) {
                        LineSweepOptions diag = base;
                        diag.order = LineOrder::diagonal;
                        amrex::MultiFab ds, dd;
                        const auto n0 = ablsteady::line_sweep_launches();
                        const amrex::Real cs = ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, ds, base);
                        const auto n1 = ablsteady::line_sweep_launches();
                        const amrex::Real cd = ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, dd, diag);
                        const auto n2 = ablsteady::line_sweep_launches();
                        // the kernels each ordering reaches: one launch per box and
                        // pass (sequential), one per plane (diagonal)
                        std::uint64_t planes = 0, boxes = 0;
                        for (amrex::MFIter mfi(ds, false); mfi.isValid(); ++mfi) {
                            const amrex::Box S = base.ghosts == LineGhosts::lagged
                                                 ? mfi.validbox()
                                                 : ablsteady::line_sweep_box(loc, ablsteady::line_owned_box(loc, mfi.validbox(), c.geom),
                                                                             base.overlap, c.geom);
                            planes += static_cast<std::uint64_t>(S.length(0) + S.length(1) - 1);
                            ++boxes;
                        }
                        EXPECT_EQ(n1 - n0, 2u * static_cast<unsigned>(sweeps) * boxes) << loc_name(loc);
                        EXPECT_EQ(n2 - n1, 2u * static_cast<unsigned>(sweeps) * planes) << loc_name(loc);
                        if (boxes > 0) { EXPECT_GT(n2 - n1, n1 - n0); }   // counts are per rank
                        EXPECT_GT(ds.norm0(), 0.1);
                        EXPECT_EQ(max_diff(ds, dd), 0) << loc_name(loc) << " periodic " << periodic << " boxes " << c.ba.size()
                                                       << " ghosts " << ablsteady::line_ghosts_name(base.ghosts) << " n_o "
                                                       << base.overlap << " sweeps " << sweeps;
                        EXPECT_EQ(cs, cd);
                    }
                }
            }
        }
    }
}

// The refusals: an overlap outside [0, max_line_overlap], an overlap with
// the lagged ghosts, an overlap wider than the narrowest box; each message
// names steady.line_overlap. No overlap along a one-cell axis (a 2-D deck).
// The names of the options.
// Beyond the plan (UT-2.2.6): the alternating ordering (steady.line_alternate).
// (a) On a system coupled only to the (i - 1) and (j + 1) neighbours (the
// upwind couplings of u > 0, v < 0) the j-reversed forward pass is an exact
// solve, so two alternating sweeps reach the solution to round-off; two
// sweeps of the ascending pairs do not. (b) The diagonal ordering equals the
// sequential one bit for bit with the j-reversed pairs too (the wavefronts
// i - j = const), at the launch counts of UT-2.2.4. (c) One sweep is not
// changed (only the odd sweeps are reversed); two sweeps are.
TEST (LineSweeps, AlternatingOrderingSweepsBothDiagonals)
{
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    {   // (a)
        PCase c = box_mesh(8, 8, 4, false, 8, 8);
        amrex::MultiFab A, ap, rhs;
        random_system(c, Loc::cell, 71u, 0.2, 3.0, A, ap, rhs);
        A.setVal(0.0, ablsteady::coef_XP, 1);
        A.setVal(0.0, ablsteady::coef_YM, 1);
        LineSweepOptions alt = lusgs(LineGhosts::lagged);
        alt.alternate = 1;
        amrex::MultiFab d0, d1, r0, r1;
        ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d0, lusgs(LineGhosts::lagged));
        ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d1, alt);
        residual_of(A, ap, rhs, d0, c.geom, r0);
        residual_of(A, ap, rhs, d1, c.geom, r1);
        const amrex::Real scale = rhs.norm0();
        EXPECT_LE(r1.norm0(), 1000 * eps * scale) << "alternating";
        EXPECT_GT(r0.norm0(), amrex::Real(1e-3) * scale) << "ascending pairs";
    }
    struct Layout { int nx, ny, mx, my; };
    const Layout layouts[] = {{8, 8, 8, 8}, {8, 8, 4, 4}, {7, 6, 4, 3}};
    for (const bool periodic : {false, true}) {
        for (const Layout& L : layouts) {
            PCase c = box_mesh(L.nx, L.ny, 4, periodic, L.mx, L.my);
            for (const Loc loc : kLocs) {
                amrex::MultiFab A, ap, rhs;
                random_system(c, loc, 73u, 0.2, 3.0, A, ap, rhs);
                for (const LineSweepOptions& plain : {lusgs(LineGhosts::lagged), lusgs(LineGhosts::zero, 2)}) {
                    LineSweepOptions seq = plain;
                    seq.alternate = 1;
                    LineSweepOptions diag = seq;
                    diag.order = LineOrder::diagonal;
                    for (const int sweeps : {1, 2, 3}) {
                        amrex::MultiFab dp, ds, dd;
                        ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, dp, plain);
                        const auto n0 = ablsteady::line_sweep_launches();
                        const amrex::Real cs = ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, ds, seq);
                        const auto n1 = ablsteady::line_sweep_launches();
                        const amrex::Real cd = ablsteady::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, dd, diag);
                        const auto n2 = ablsteady::line_sweep_launches();
                        std::uint64_t planes = 0, boxes = 0;
                        for (amrex::MFIter mfi(ds, false); mfi.isValid(); ++mfi) {
                            const amrex::Box S = plain.ghosts == LineGhosts::lagged
                                                 ? mfi.validbox()
                                                 : ablsteady::line_sweep_box(loc, ablsteady::line_owned_box(loc, mfi.validbox(), c.geom),
                                                                             plain.overlap, c.geom);
                            planes += static_cast<std::uint64_t>(S.length(0) + S.length(1) - 1);
                            ++boxes;
                        }
                        EXPECT_EQ(n1 - n0, 2u * static_cast<unsigned>(sweeps) * boxes) << loc_name(loc);
                        EXPECT_EQ(n2 - n1, 2u * static_cast<unsigned>(sweeps) * planes) << loc_name(loc);
                        EXPECT_EQ(max_diff(ds, dd), 0) << loc_name(loc) << " periodic " << periodic << " boxes " << c.ba.size()
                                                       << " ghosts " << ablsteady::line_ghosts_name(plain.ghosts) << " sweeps " << sweeps;
                        EXPECT_EQ(cs, cd);
                        const double gap = max_diff(dp, ds);
                        if (sweeps == 1) {
                            EXPECT_EQ(gap, 0) << loc_name(loc);
                        } else {
                            EXPECT_GT(gap, 1e-6) << loc_name(loc) << " sweeps " << sweeps;
                        }
                    }
                }
            }
        }
    }
    // the range check names its key
    LineSweepOptions bad = lusgs(LineGhosts::lagged);
    for (const int v : {-1, 2}) {
        bad.alternate = v;
        const auto err = ablsteady::check_line_options(bad);
        ASSERT_EQ(err.size(), 1u) << v;
        EXPECT_NE(err.front().find("steady.line_alternate"), std::string::npos) << err.front();
    }
}

TEST (LineSweeps, RefusalsAndNames)
{
    PCase c = box_mesh(8, 8, 4, true, 4, 4);
    amrex::MultiFab A, ap, rhs;
    random_system(c, Loc::cell, 61u, 0.3, 1.0, A, ap, rhs);
    for (const LineSweepOptions& bad : {lusgs(LineGhosts::zero, -1), lusgs(LineGhosts::zero, ablsteady::max_line_overlap + 1),
                                        lusgs(LineGhosts::lagged, 2), lusgs(LineGhosts::zero, 5)}) {
        amrex::MultiFab delta;
        try {
            ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 1, c.mesh, c.geom, delta, bad);
            ADD_FAILURE() << "accepted overlap " << bad.overlap;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("steady.line_overlap"), std::string::npos) << e.what();
        }
    }
    amrex::MultiFab delta;   // the widest the 4 x 4 boxes allow
    EXPECT_NO_THROW(ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 1, c.mesh, c.geom, delta, lusgs(LineGhosts::zero, 4)));
    // a 2-D deck (one cell in y): the overlap applies in x alone
    {
        PCase c2 = box_mesh(8, 1, 4, true, 4, 1);
        EXPECT_EQ(ablsteady::line_overlap_width(3, c2.geom.Domain()), amrex::IntVect(3, 0, 0));
        EXPECT_EQ(ablsteady::line_overlap_width(3, c.geom.Domain()), amrex::IntVect(3, 3, 0));
        EXPECT_TRUE(ablsteady::check_line_overlap_layout(c2.ba, 4).empty());
        EXPECT_FALSE(ablsteady::check_line_overlap_layout(c2.ba, 5).empty());
        amrex::MultiFab A2, ap2, r2, d2;
        random_system(c2, Loc::xface, 62u, 0.3, 1.0, A2, ap2, r2);
        EXPECT_NO_THROW(ablsteady::zline_lusgs(Loc::xface, A2, ap2, r2, 1, c2.mesh, c2.geom, d2, lusgs(LineGhosts::zero, 4)));
    }
    EXPECT_STREQ(ablsteady::line_ghosts_name(LineGhosts::lagged), "lagged");
    EXPECT_STREQ(ablsteady::line_ghosts_name(LineGhosts::zero), "zero");
    EXPECT_STREQ(ablsteady::line_order_name(LineOrder::sequential), "sequential");
    EXPECT_STREQ(ablsteady::line_order_name(LineOrder::diagonal), "diagonal");
    // the dispatcher routes the options: zero ghosts through line_sweeps
    amrex::MultiFab d1, d2;
    ablsteady::line_sweeps(lusgs(LineGhosts::zero, 2), Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d1);
    ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d2, lusgs(LineGhosts::zero, 2));
    EXPECT_EQ(max_diff(d1, d2), 0);
    ablsteady::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d2, lusgs(LineGhosts::lagged));
    EXPECT_GT(max_diff(d1, d2), 0);
}
