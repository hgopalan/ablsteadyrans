// Transfer operators between the levels of the grid
// hierarchy (Transfer.H).
//   restriction conserves the face mass flux (coarse = sum of fine);
//            the coarse divergence is the J-weighted mean of the fine one
//   prolongation + projection gives a divergence-free field
//   prolongation of bilinear fields is exact (flat terrain), up to
//            the domain faces and across a periodic seam
//   prolonged k >= k_floor
//   also: a gathered level gives the transfers of the distributed
//            one bit for bit; restrict(prolong(c)) on flat ground is the
//            closed-form smoother S_x S_y; the target's node heights come
//            from the hierarchy; a coarse direction one cell wide gives a
//            constant; a wrong level pair aborts
// The pinned values are computed on the host from closed forms, in double
// precision, from host copies of the fields gathered on every rank.
#include "BoundaryConditions.H"
#include "Forcing.H"
#include "Inputs.H"
#include "MassFlux.H"
#include "Mesh.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"
#include "Transfer.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double tol (double dbl, double flt) { return kDouble ? dbl : flt; }

using gust::FaceType;

// The run's level 0 (as App.cpp lays it out) and the hierarchy built on it.
struct TCase
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab z;
    gust::Mesh mesh;
    std::unique_ptr<gust::MeshHierarchy> h;

    void make (int gather_from = 0)
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::hierarchy_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        z.define(amrex::convert(ba, amrex::IntVect(1)), dm, 1, gust::ngrow_nodal);
        z.setVal(0.0);
        mesh.define(ba, dm);
        mesh.build(in, geom, z);
        h = std::make_unique<gust::MeshHierarchy>(in, geom, z, mesh, gather_from);
    }
};

// 32 x 24 x 10 cells in boxes of at most 12: 6 boxes of unequal size on
// level 0 (12, 12, 8 by 12, 12; uneven on 4 ranks), three levels (level 2:
// 8 x 6 cells). periodic_x / periodic_y choose the
// periodic direction; hill: a round Agnesi hill (STF) instead of flat ground.
// The caller calls make(): the hierarchy keeps references to level 0.
TCase case_inputs (bool periodic_x, bool periodic_y, bool hill)
{
    TCase c;
    c.in.mesh.n_cell = amrex::IntVect(32, 24, 10);
    c.in.mesh.prob_lo = {0, 0, 0};
    c.in.mesh.prob_hi = {3200, 2400, 1000};
    c.in.mesh.is_periodic = {periodic_x ? 1 : 0, periodic_y ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {12, 12};
    c.in.mesh.dz0 = 40;
    c.in.mesh.hierarchy_levels = 3;
    if (hill) {
        c.in.terrain.type = "agnesi";
        c.in.terrain.hill_shape = "round";
        c.in.terrain.hill_height = 150;
        c.in.terrain.hill_half_width = 300;
        c.in.terrain.smoothing = 1;
    }
    c.in.terrain.hierarchy_passes = 1;
    return c;
}

// The valid points of a MultiFab on every rank, by global index (a point
// shared by boxes holds their mean, which is its one value).
struct HostField
{
    amrex::Box dom;
    std::vector<double> v;
    [[nodiscard]] double operator() (int i, int j, int k) const
    {
        const amrex::IntVect l = dom.smallEnd();
        const amrex::IntVect n = dom.length();
        return v[(static_cast<std::size_t>(k - l[2]) * n[1] + (j - l[1])) * n[0] + (i - l[0])];
    }
};

HostField gather (const amrex::MultiFab& mf)
{
    HostField h;
    h.dom = mf.boxArray().minimalBox();
    const auto npts = static_cast<std::size_t>(h.dom.numPts());
    std::vector<double> sum(npts, 0.0), cnt(npts, 0.0);
    const amrex::IntVect l = h.dom.smallEnd();
    const amrex::IntVect n = h.dom.length();
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(mf[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const auto m = (static_cast<std::size_t>(k - l[2]) * n[1] + (j - l[1])) * n[0] + (i - l[0]);
            sum[m] += static_cast<double>(f(amrex::IntVect(i, j, k)));
            cnt[m] += 1;
        });
    }
    const auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Sum(sum.data(), static_cast<int>(npts), comm);
    amrex::ParallelAllReduce::Sum(cnt.data(), static_cast<int>(npts), comm);
    h.v.resize(npts);
    for (std::size_t m = 0; m < npts; ++m) { h.v[m] = sum[m] / cnt[m]; }
    return h;
}

// The horizontal position of point (i, j) of mf's type on geometry g.
double xpos (const amrex::Geometry& g, const amrex::MultiFab& mf, int i, int d)
{
    const double off = mf.ixType().nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// Fill every point of mf (ghosts included) with fn(x, y, k). Along a
// periodic direction the index is wrapped into the domain first, so the
// field is periodic bit for bit (the two copies of a nodal seam point
// agree, as OverrideSync keeps them in the solver).
void fill_xy (amrex::MultiFab& mf, const amrex::Geometry& g, const std::function<double(double, double, int)>& fn)
{
    const amrex::Box dom = g.Domain();
    auto wrap = [&] (int i, int d) {
        if (!g.isPeriodic(d)) { return i; }
        const int n = dom.length(d);
        return ((i - dom.smallEnd(d)) % n + n) % n + dom.smallEnd(d);
    };
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(xpos(g, mf, wrap(i, 0), 0), xpos(g, mf, wrap(j, 1), 1), k));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

// Fill every point of mf (ghosts included) with fn(i, j, k) of its index.
void fill_ijk (amrex::MultiFab& mf, const std::function<double(int, int, int)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

void define_level (const TCase& c, int l, gust::State& q)
{
    q.define(c.h->boxArray(l), c.h->DistributionMap(l));
    amrex::MultiFab::Copy(q.z_phys_nd, c.h->z_phys_nd(l), 0, 0, 1, gust::ngrow_nodal);
}

// A smooth, non-uniform state on level l (ghosts included): a wind with
// shear and a horizontal wave, w from a small wave, theta, k and p'.
void smooth_state (const TCase& c, int l, gust::State& q)
{
    const amrex::Geometry& g = c.h->geom(l);
    const double Lx = 3200, Ly = 2400, pi = 3.14159265358979;
    fill_xy(q.u, g, [=] (double x, double y, int k) { return 8.0 + 0.3 * k + 1.5 * std::sin(2 * pi * x / Lx) * std::cos(2 * pi * y / Ly); });
    fill_xy(q.v, g, [=] (double x, double y, int k) { return 1.0 - 0.1 * k + 0.8 * std::cos(2 * pi * x / Lx + 0.3) * std::sin(4 * pi * y / Ly); });
    fill_xy(q.w, g, [=] (double x, double y, int k) { return 0.2 * std::sin(2 * pi * x / Lx) * std::sin(2 * pi * y / Ly + 0.1 * k); });
    fill_xy(q.theta, g, [=] (double x, double y, int k) { return 300.0 + 0.5 * k + 0.2 * std::cos(2 * pi * (x + y) / Lx); });
    fill_xy(q.tke, g, [=] (double x, double y, int k) { return 0.5 * std::exp(-0.2 * k) + 0.1 * std::sin(2 * pi * x / Lx) * std::sin(2 * pi * y / Ly) + 0.2; });
    fill_xy(q.p, g, [=] (double x, double y, int k) { return 3.0 * std::cos(2 * pi * x / Lx) + 2.0 * std::sin(2 * pi * y / Ly) - 0.1 * k; });
}

double max_abs (const std::vector<double>& v)
{
    double m = 0;
    for (const double x : v) { m = std::max(m, std::abs(x)); }
    return m;
}

} // namespace

// ---------------------------------------------------------------------------
// Restriction conserves the face mass flux: every coarse face
// carries the sum of the fluxes of the fine faces it covers,
//   a_x,c u_c r = sum a_x,f u_f,  a_y,c v_c r = sum a_y,f v_f,
//   a_z,c Omega_c r^2 = sum a_z,f Omega_f,
// to 1e-12 of the flux scale, on flat ground and over a hill (whose coarse
// areas are not sums of the fine ones), for r = 2 and 4. So the coarse
// divergence is the J-weighted mean of the fine one, r^2 J_c D_c = sum J_f
// D_f, and the cell fields are the detJ-weighted means.
TEST (Transfer, RestrictionConservesFaceFlux)
{
    const double rel = tol(1e-12, 2e-5);
    for (const bool hill : {false, true}) {
        TCase c = case_inputs(false, true, hill);
        c.make();
        const std::pair<int, int> pairs[] = {{0, 1}, {1, 2}, {0, 2}};
        for (const auto& [lf, lc] : pairs) {
            const std::string where = std::string(hill ? "hill" : "flat") + ", " + std::to_string(lf) + " -> " + std::to_string(lc);
            const int r = 1 << (lc - lf);
            gust::State qf, qc;
            define_level(c, lf, qf);
            define_level(c, lc, qc);
            smooth_state(c, lf, qf);
            amrex::MultiFab om_c;
            gust::restrict_state(*c.h, lf, qf, lc, qc, &om_c);

            const gust::Mesh& mf = c.h->mesh(lf);
            const gust::Mesh& mc = c.h->mesh(lc);
            amrex::MultiFab om_f(qf.w.boxArray(), qf.w.DistributionMap(), 1, 0);
            gust::omega_from_w(qf.u, qf.v, qf.w, c.h->z_phys_nd(lf), mf, c.h->geom(lf), om_f, 0);
            const HostField uf = gather(qf.u), vf = gather(qf.v), of = gather(om_f);
            const HostField axf = gather(mf.ax), ayf = gather(mf.ay), azf = gather(mf.az), Jf = gather(mf.detJ);
            const HostField uc = gather(qc.u), vc = gather(qc.v), oc = gather(om_c);
            const HostField axc = gather(mc.ax), ayc = gather(mc.ay), azc = gather(mc.az), Jc = gather(mc.detJ);
            const HostField tf = gather(qf.theta), tc = gather(qc.theta), kf = gather(qf.tke), kc = gather(qc.tke);
            const HostField pf = gather(qf.p), pc = gather(qc.p);

            // x-faces, y-faces, z-faces.
            double ex = 0, sx = 0, ey = 0, sy = 0, ez = 0, sz = 0;
            const amrex::Box& cd = c.h->geom(lc).Domain();
            const int NX = cd.length(0), NY = cd.length(1), NZ = cd.length(2);
            bool areas_differ = false;
            for (int k = 0; k < NZ; ++k) {
                for (int J = 0; J < NY; ++J) {
                    for (int I = 0; I <= NX; ++I) {
                        double s = 0, sa = 0;
                        for (int b = 0; b < r; ++b) {
                            s += axf(r * I, r * J + b, k) * uf(r * I, r * J + b, k);
                            sa += axf(r * I, r * J + b, k);
                        }
                        ex = std::max(ex, std::abs(axc(I, J, k) * uc(I, J, k) * r - s));
                        sx = std::max(sx, std::abs(s));
                        if (std::abs(sa - r * axc(I, J, k)) > 1e-3 * sa) { areas_differ = true; }
                    }
                }
                for (int J = 0; J <= NY; ++J) {
                    for (int I = 0; I < NX; ++I) {
                        double s = 0;
                        for (int a = 0; a < r; ++a) { s += ayf(r * I + a, r * J, k) * vf(r * I + a, r * J, k); }
                        ey = std::max(ey, std::abs(ayc(I, J, k) * vc(I, J, k) * r - s));
                        sy = std::max(sy, std::abs(s));
                    }
                }
            }
            for (int k = 0; k <= NZ; ++k) {
                for (int J = 0; J < NY; ++J) {
                    for (int I = 0; I < NX; ++I) {
                        double s = 0;
                        for (int b = 0; b < r; ++b) {
                            for (int a = 0; a < r; ++a) { s += azf(r * I + a, r * J + b, k) * of(r * I + a, r * J + b, k); }
                        }
                        ez = std::max(ez, std::abs(azc(I, J, k) * oc(I, J, k) * r * r - s));
                        sz = std::max(sz, std::abs(s));
                    }
                }
            }
            EXPECT_LT(ex, rel * sx) << "x-face flux, " << where;
            EXPECT_LT(ey, rel * sy) << "y-face flux, " << where;
            EXPECT_LT(ez, rel * sz) << "z-face flux, " << where;
            EXPECT_GT(sz, 1e-3 * sx) << "the vertical flux is not vacuous, " << where;
            // Over the hill the coarse areas are not the sums of the fine ones
            // (the coarse ground is smoothed): the area weight matters.
            EXPECT_EQ(areas_differ, hill) << where;

            // Cell fields: q_c sum J_f = sum J_f q_f.
            double ec = 0, sc = 0;
            bool J_varies = false;
            const std::pair<const HostField*, const HostField*> cells[3] = {{&tf, &tc}, {&kf, &kc}, {&pf, &pc}};
            for (const auto& [f, co] : cells) {
                for (int k = 0; k < NZ; ++k) {
                    for (int J = 0; J < NY; ++J) {
                        for (int I = 0; I < NX; ++I) {
                            double sj = 0, sjq = 0, jmin = 1e30, jmax = 0;
                            for (int b = 0; b < r; ++b) {
                                for (int a = 0; a < r; ++a) {
                                    const double w = Jf(r * I + a, r * J + b, k);
                                    sj += w;
                                    sjq += w * (*f)(r * I + a, r * J + b, k);
                                    jmin = std::min(jmin, w);
                                    jmax = std::max(jmax, w);
                                }
                            }
                            ec = std::max(ec, std::abs((*co)(I, J, k) - sjq / sj));
                            sc = std::max(sc, std::abs(sjq / sj));
                            if (jmax - jmin > 1e-3 * jmax) { J_varies = true; }
                        }
                    }
                }
            }
            EXPECT_LT(ec, tol(4e-15, 4e-7) * sc) << "detJ-weighted cell means, " << where;
            EXPECT_EQ(J_varies, hill) << where;

            // The divergence: r^2 J_c D_c = sum J_f D_f.
            amrex::MultiFab Df(qf.p.boxArray(), qf.p.DistributionMap(), 1, 0), Dc(qc.p.boxArray(), qc.p.DistributionMap(), 1, 0);
            gust::velocity_divergence(qf.u, qf.v, qf.w, c.h->z_phys_nd(lf), mf, c.h->geom(lf), Df);
            gust::face_divergence(mc, c.h->geom(lc), {&qc.u, &qc.v, &om_c}, Dc);
            const HostField df = gather(Df), dc = gather(Dc);
            double ed = 0, sd = 0;
            for (int k = 0; k < NZ; ++k) {
                for (int J = 0; J < NY; ++J) {
                    for (int I = 0; I < NX; ++I) {
                        double s = 0;
                        for (int b = 0; b < r; ++b) {
                            for (int a = 0; a < r; ++a) { s += Jf(r * I + a, r * J + b, k) * df(r * I + a, r * J + b, k); }
                        }
                        ed = std::max(ed, std::abs(r * r * Jc(I, J, k) * dc(I, J, k) - s));
                        sd = std::max(sd, std::abs(s));
                    }
                }
            }
            // The flux scale per cell: |u| / dx.
            const double dscale = sx / static_cast<double>(c.h->geom(lc).CellSize(0)) * r;
            EXPECT_LT(ed, rel * dscale) << "divergence, " << where;
            EXPECT_GT(sd, 1e-3 * dscale) << "the divergence is not vacuous, " << where;
        }
    }
}

// ---------------------------------------------------------------------------
// Prolongation of a field bilinear in x and y is exact on flat
// terrain, at every fine point: up to a non-periodic domain face (linear
// extrapolation in the half cell next to it) and across a periodic seam. In
// the periodic direction the field is a tent of period L with its two kinks
// on coarse points of its own type (cell centres or faces), so it is linear
// between any two neighbouring coarse points, across the seam too.
TEST (Transfer, ProlongationOfBilinearFieldsIsExact)
{
    const double rel = tol(1e-13, 2e-6);
    for (const bool px : {false, true}) {
        TCase c = case_inputs(px, !px, false);
        c.make();
        const std::pair<int, int> pairs[] = {{0, 1}, {1, 2}, {0, 2}};
        for (const auto& [lf, lc] : pairs) {
            const amrex::Geometry& gc = c.h->geom(lc);
            gust::State qf, qc;
            define_level(c, lf, qf);
            define_level(c, lc, qc);
            amrex::MultiFab* cf[6] = {&qc.u, &qc.v, &qc.w, &qc.theta, &qc.tke, &qc.p};
            amrex::MultiFab* ff[6] = {&qf.u, &qf.v, &qf.w, &qf.theta, &qf.tke, &qf.p};
            const char* names[6] = {"u", "v", "w", "theta", "k", "p"};
            std::vector<std::function<double(double, double, int)>> fn(6);
            for (int n = 0; n < 6; ++n) {
                // The periodic direction pd, its period and the kinks on the
                // coarse points of this field's type.
                const int pd = px ? 0 : 1;
                const double L = static_cast<double>(gc.ProbLength(pd));
                const double D = static_cast<double>(gc.CellSize(pd));
                const bool nodal = cf[n]->ixType().nodeCentered(pd);
                const int N = gc.Domain().length(pd);
                const double ya = (N / 4 + (nodal ? 0.0 : 0.5)) * D;
                const double a = 1.0 + n, b = 0.01 * (n + 1), cc = 0.02 - 0.003 * n, e = 1e-5 * (n + 2), g = 0.7 * n - 1;
                fn[n] = [=] (double x, double y, int k) {
                    const double s = px ? x : y;   // periodic coordinate
                    const double t = px ? y : x;   // the other one
                    double m = std::fmod(s - ya, L);
                    if (m < 0) { m += L; }
                    const double tent = std::abs(m - L / 2);
                    return (a + b * t) * (1.0 + cc * tent) + e * t * tent + g * k;
                };
                fill_xy(*cf[n], gc, fn[n]);
                // Coarse ghost cells are not read: a large value would show.
                cf[n]->setBndry(amrex::Real(1e30));
                ff[n]->setVal(amrex::Real(-1e30));
            }
            gust::prolong_fields(*c.h, lc, qc, lf, qf);
            const amrex::Geometry& gf = c.h->geom(lf);
            for (int n = 0; n < 6; ++n) {
                const HostField h = gather(*ff[n]);
                double err = 0, scale = 0;
                const amrex::Box& d = h.dom;
                for (int k = d.smallEnd(2); k <= d.bigEnd(2); ++k) {
                    for (int j = d.smallEnd(1); j <= d.bigEnd(1); ++j) {
                        for (int i = d.smallEnd(0); i <= d.bigEnd(0); ++i) {
                            const double ex = fn[n](xpos(gf, *ff[n], i, 0), xpos(gf, *ff[n], j, 1), k);
                            err = std::max(err, std::abs(h(i, j, k) - ex));
                            scale = std::max(scale, std::abs(ex));
                        }
                    }
                }
                EXPECT_LT(err, rel * scale) << names[n] << ", periodic " << (px ? "x" : "y") << ", " << lc << " -> " << lf;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// restrict(prolong(c)) on flat ground. Bilinear
// interpolation is not conservative, so the round trip is not the identity:
// for r = 2 it is the smoother S_x S_y with, along a cell-centred direction,
//   (S c)_I = (c_{I-1} + 6 c_I + c_{I+1}) / 8
// in the interior and across a periodic seam, and (S c)_I = c_I in the
// first and last cells next to a non-periodic face (the extrapolated and
// the interpolated half cells average to c_I there), and S = 1 along a
// nodal direction (u in x, v in y: every coarse face is a fine face).
// So theta, k, p', w -> S_x S_y; u -> S_y; v -> S_x. Pinned against a
// non-smooth coarse field.
namespace {
double smooth1 (const std::function<double(int)>& c, int I, int N, bool periodic)
{
    if (!periodic && (I == 0 || I == N - 1)) { return c(I); }
    const int im = periodic ? (I - 1 + N) % N : I - 1;
    const int ip = periodic ? (I + 1) % N : I + 1;
    return (c(im) + 6 * c(I) + c(ip)) / 8;
}
} // namespace

TEST (Transfer, RoundTripIsTheClosedFormSmoother)
{
    const double rel = tol(2e-14, 2e-6);
    for (const bool px : {false, true}) {
        TCase c = case_inputs(px, !px, false);
        c.make();
        const int lf = 0, lc = 1;
        gust::State qc, qf, back;
        define_level(c, lc, qc);
        define_level(c, lf, qf);
        define_level(c, lc, back);
        amrex::MultiFab* cf[6] = {&qc.u, &qc.v, &qc.w, &qc.theta, &qc.tke, &qc.p};
        amrex::MultiFab* bf[6] = {&back.u, &back.v, &back.w, &back.theta, &back.tke, &back.p};
        const char* names[6] = {"u", "v", "w", "theta", "k", "p"};
        for (int n = 0; n < 6; ++n) {
            fill_ijk(*cf[n], [n] (int i, int j, int k) {
                return 1.0 + static_cast<double>(gust::test::pattern(i + 40, j + 40, k + 40, n)) * 1e-5;
            });
        }
        // The periodic seam faces are one face: give them one value.
        qc.u.OverrideSync(c.h->geom(lc).periodicity());
        qc.v.OverrideSync(c.h->geom(lc).periodicity());
        gust::prolong_fields(*c.h, lc, qc, lf, qf);
        qf.u.FillBoundary(c.h->geom(lf).periodicity());
        qf.v.FillBoundary(c.h->geom(lf).periodicity());
        gust::restrict_state(*c.h, lf, qf, lc, back);
        const amrex::Box& cd = c.h->geom(lc).Domain();
        const int NX = cd.length(0), NY = cd.length(1);
        for (int n = 0; n < 6; ++n) {
            const HostField hc = gather(*cf[n]), hb = gather(*bf[n]);
            const bool nx = cf[n]->ixType().nodeCentered(0), ny = cf[n]->ixType().nodeCentered(1);
            double err = 0, change = 0, scale = 0;
            const amrex::Box& d = hb.dom;
            for (int k = d.smallEnd(2); k <= d.bigEnd(2); ++k) {
                for (int J = d.smallEnd(1); J <= d.bigEnd(1); ++J) {
                    for (int I = d.smallEnd(0); I <= d.bigEnd(0); ++I) {
                        // S_y of the S_x-smoothed column values.
                        auto sx = [&] (int jj) {
                            if (nx) { return hc(I, jj, k); }
                            return smooth1([&] (int ii) { return hc(ii, jj, k); }, I, NX, px);
                        };
                        const double ex = ny ? sx(J) : smooth1(sx, J, NY, !px);
                        err = std::max(err, std::abs(hb(I, J, k) - ex));
                        change = std::max(change, std::abs(ex - hc(I, J, k)));
                        scale = std::max(scale, std::abs(ex));
                    }
                }
            }
            EXPECT_LT(err, rel * scale) << names[n] << ", periodic " << (px ? "x" : "y");
            EXPECT_GT(change, 1e-7) << "not vacuous: " << names[n];
        }
    }
}

// ---------------------------------------------------------------------------
// The prolonged k is at least the floor: a coarse k that dips
// below the floor (and below 0) gives fine values equal to the bilinear
// interpolant where it is above the floor and to the floor exactly where it
// is not; the clipped count is the number of such cells.
TEST (Transfer, ProlongedTkeIsAboveTheFloor)
{
    TCase c = case_inputs(true, false, false);
    c.make();
    const int lf = 0, lc = 1;
    const amrex::Real floor = amrex::Real(1e-4);
    gust::State qc, qf, interp;
    define_level(c, lc, qc);
    define_level(c, lf, qf);
    define_level(c, lf, interp);
    smooth_state(c, lc, qc);
    const double pi = 3.14159265358979;
    fill_xy(qc.tke, c.h->geom(lc), [=] (double x, double y, int k) {
        return 1e-4 * (0.6 + 1.2 * std::sin(2 * pi * x / 3200 + 0.4 * k) * std::cos(2 * pi * y / 2400));
    });
    // Periodic x, walls in y, slip ground and top: no wall model needed.
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::slip_wall, FaceType::slip_wall, FaceType::periodic, FaceType::slip_wall,
                FaceType::slip_wall};
    const gust::BoundaryConditions bc(bin, c.h->geom(lf), c.h->mesh(lf), gust::CoriolisOptions{});
    gust::PoissonOptions popt;
    auto solver = c.h->make_poisson(lf, popt);
    const gust::ProlongResult res = gust::prolong_state(*c.h, lc, qc, lf, qf, bc, nullptr,
                                                                   gust::BuoyancyOptions{}, *solver, floor);
    gust::prolong_fields(*c.h, lc, qc, lf, interp);
    const HostField k = gather(qf.tke), ki = gather(interp.tke);
    long below = 0;
    double err = 0, kmin = 1e30, kimin = 1e30;
    for (std::size_t m = 0; m < k.v.size(); ++m) {
        const double expect = ki.v[m] < floor ? static_cast<double>(floor) : ki.v[m];
        if (ki.v[m] < floor) { ++below; }
        err = std::max(err, std::abs(k.v[m] - expect));
        kmin = std::min(kmin, k.v[m]);
        kimin = std::min(kimin, ki.v[m]);
    }
    EXPECT_GE(kmin, static_cast<double>(floor));
    EXPECT_EQ(err, 0.0);
    EXPECT_LT(kimin, 0.0) << "not vacuous: the interpolant dips below 0";
    EXPECT_GT(below, 100);
    EXPECT_EQ(res.tke_clipped, below);
    // The floor's range check aborts, naming turb.tke_floor.
    for (const amrex::Real bad : {amrex::Real(-1e-6), std::numeric_limits<amrex::Real>::quiet_NaN(),
                                  std::numeric_limits<amrex::Real>::infinity()}) {
        try {
            (void)gust::prolong_state(*c.h, lc, qc, lf, qf, bc, nullptr, gust::BuoyancyOptions{}, *solver, bad);
            ADD_FAILURE() << "no abort for tke_floor = " << bad;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("turb.tke_floor"), std::string::npos) << e.what();
        }
    }
}

// ---------------------------------------------------------------------------
// Prolongation + projection gives a divergence-free field: the
// divergence left is the Poisson solver's residual, below its tolerance in
// the norm its back end stops on (gmres_mlmg: the J-weighted 2-norm, the
// default; gmres_fft: the plain 2-norm). The Dirichlet normal velocities of
// the fine level's boundary conditions (inflow, walls, ground, top) keep
// the values fill_ghosts gave them, bit for bit. Flat ground and a hill,
// inflow / outflow in x, periodic y, the surface layer, a slip top.
namespace {
double norm2_weighted (const amrex::MultiFab& div, const amrex::MultiFab& J, bool weighted)
{
    amrex::MultiFab t(div.boxArray(), div.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(t, div, 0, 0, 1, 0);
    if (weighted) { amrex::MultiFab::Multiply(t, J, 0, 0, 1, 0); }
    return static_cast<double>(t.norm2());
}
} // namespace

TEST (Transfer, ProlongationAndProjectionIsDivergenceFree)
{
    for (const bool hill : {false, true}) {
        for (const std::string backend : {"gmres_mlmg", "gmres_fft"}) {
            TCase c = case_inputs(false, true, hill);
            c.make();
            const int lf = 0, lc = 1;
            const std::string where = std::string(hill ? "hill" : "flat") + ", " + backend;
            gust::State qc, qf, ref;
            define_level(c, lc, qc);
            define_level(c, lf, qf);
            define_level(c, lf, ref);
            smooth_state(c, lc, qc);
            gust::BoundaryInputs bin;
            bin.face = {FaceType::inflow, FaceType::periodic, FaceType::surface_layer, FaceType::outflow,
                        FaceType::periodic, FaceType::slip_wall};
            const gust::BoundaryConditions bc(bin, c.h->geom(lf), c.h->mesh(lf), gust::CoriolisOptions{});
            const gust::SurfaceLayerOptions so;
            const gust::BuoyancyOptions buoy;
            gust::PoissonOptions popt;
            popt.backend = backend;
            popt.reltol = amrex::Real(tol(1e-9, 1e-4));
            auto solver = c.h->make_poisson(lf, popt);
            const gust::ProlongResult res = gust::prolong_state(*c.h, lc, qc, lf, qf, bc, &so, buoy, *solver,
                                                                           amrex::Real(1e-6));
            // The interpolation and the boundary fill alone.
            gust::prolong_fields(*c.h, lc, qc, lf, ref);
            bc.fill_ghosts(ref, &so, buoy);

            const gust::Mesh& m = c.h->mesh(lf);
            const amrex::Geometry& g = c.h->geom(lf);
            amrex::MultiFab D0(qf.p.boxArray(), qf.p.DistributionMap(), 1, 0), D1(qf.p.boxArray(), qf.p.DistributionMap(), 1, 0);
            gust::velocity_divergence(ref.u, ref.v, ref.w, ref.z_phys_nd, m, g, D0);
            gust::velocity_divergence(qf.u, qf.v, qf.w, qf.z_phys_nd, m, g, D1);
            const bool weighted = (backend != "gmres_fft");
            const double n0 = norm2_weighted(D0, m.detJ, weighted), n1 = norm2_weighted(D1, m.detJ, weighted);
            EXPECT_GT(res.projection.divergence_before, 1e-4) << "not vacuous, " << where;
            EXPECT_GT(n0, 1e-3) << where;
            EXPECT_LT(n1, 1.01 * static_cast<double>(popt.reltol) * n0) << where;
            // The plain max norm the pressure correction reports, from the same field.
            EXPECT_EQ(static_cast<double>(res.projection.divergence_after), static_cast<double>(D1.norm0())) << where;

            // The Dirichlet face points keep their boundary values; p' is the
            // interpolated one; theta and k are untouched.
            const amrex::Box& dom = g.Domain();
            const HostField u = gather(qf.u), ur = gather(ref.u), w = gather(qf.w), wr = gather(ref.w);
            const HostField om = [&] {
                amrex::MultiFab o(qf.w.boxArray(), qf.w.DistributionMap(), 1, 0);
                gust::omega_from_w(qf.u, qf.v, qf.w, qf.z_phys_nd, m, g, o, 0);
                return gather(o);
            }();
            double eu = 0, eo = 0, du = 0;
            for (int k = 0; k < dom.length(2); ++k) {
                for (int j = 0; j < dom.length(1); ++j) {
                    eu = std::max(eu, std::abs(u(0, j, k) - ur(0, j, k)));   // inflow face
                    for (int i = 1; i < dom.length(0); ++i) { du = std::max(du, std::abs(u(i, j, k) - ur(i, j, k))); }
                }
            }
            for (int j = 0; j < dom.length(1); ++j) {
                for (int i = 0; i < dom.length(0); ++i) {
                    eo = std::max({eo, std::abs(om(i, j, 0)), std::abs(om(i, j, dom.length(2)))});
                }
            }
            EXPECT_EQ(eu, 0.0) << "inflow u, " << where;
            EXPECT_LT(eo, tol(1e-12, 1e-5) * max_abs(u.v)) << "no flux through the ground and the top, " << where;
            EXPECT_GT(du, 1e-3) << "the projection moved the interior faces, " << where;
            EXPECT_EQ(gust::test::count_bit_differences(qf.p, ref.p, 0), 0) << where;
            EXPECT_EQ(gust::test::count_bit_differences(qf.theta, ref.theta, 0), 0) << where;
        }
    }
}

// ---------------------------------------------------------------------------
// MeshHierarchy.GatheredLevelEqualsDistributed for the transfers: with the
// coarse levels gathered onto one box (gather_from = 1) the restriction, the
// interpolation and, into the distributed level 0, the whole prolongation give
// the distributed hierarchy's results bit for bit, over the hill.
TEST (Transfer, GatheredLevelEqualsDistributed)
{
    TCase a = case_inputs(false, true, true);
    a.make();
    TCase b = case_inputs(false, true, true);
    b.make(1);
    ASSERT_TRUE(b.h->gathered(1));
    ASSERT_GT(a.h->boxArray(1).size(), 1);
    ASSERT_NE(a.h->boxArray(0)[0].numPts(), a.h->boxArray(0)[a.h->boxArray(0).size() - 1].numPts()) << "unequal boxes";
    // A gathered field copied back onto the distributed layout.
    auto back = [] (const amrex::MultiFab& src, const amrex::MultiFab& like) {
        amrex::MultiFab out(like.boxArray(), like.DistributionMap(), 1, 0);
        out.setVal(0.0);
        out.ParallelCopy(src, 0, 0, 1);
        return out;
    };
    auto diffs = [&] (const gust::State& qa, const gust::State& qb) {
        long n = 0;
        const amrex::MultiFab* fa[6] = {&qa.u, &qa.v, &qa.w, &qa.theta, &qa.tke, &qa.p};
        const amrex::MultiFab* fb[6] = {&qb.u, &qb.v, &qb.w, &qb.theta, &qb.tke, &qb.p};
        for (int f = 0; f < 6; ++f) { n += gust::test::count_bit_differences(back(*fb[f], *fa[f]), *fa[f], 0); }
        return n;
    };
    const std::pair<int, int> pairs[] = {{0, 1}, {1, 2}, {0, 2}};
    for (const auto& [lf, lc] : pairs) {
        const std::string where = std::to_string(lf) + " <-> " + std::to_string(lc);
        {
            gust::State fa, fb, ca, cb;
            define_level(a, lf, fa);
            define_level(b, lf, fb);
            define_level(a, lc, ca);
            define_level(b, lc, cb);
            smooth_state(a, lf, fa);
            smooth_state(b, lf, fb);
            gust::restrict_state(*a.h, lf, fa, lc, ca);
            gust::restrict_state(*b.h, lf, fb, lc, cb);
            EXPECT_EQ(diffs(ca, cb), 0) << "restriction, " << where;
        }
        {
            gust::State fa, fb, ca, cb;
            define_level(a, lf, fa);
            define_level(b, lf, fb);
            define_level(a, lc, ca);
            define_level(b, lc, cb);
            smooth_state(a, lc, ca);
            smooth_state(b, lc, cb);
            gust::prolong_fields(*a.h, lc, ca, lf, fa);
            gust::prolong_fields(*b.h, lc, cb, lf, fb);
            EXPECT_EQ(diffs(fa, fb), 0) << "interpolation, " << where;
        }
    }
    // The whole prolongation into level 0 (distributed in both).
    for (const int lc : {1, 2}) {
        gust::State fa, fb, ca, cb;
        define_level(a, 0, fa);
        define_level(b, 0, fb);
        define_level(a, lc, ca);
        define_level(b, lc, cb);
        smooth_state(a, lc, ca);
        smooth_state(b, lc, cb);
        gust::BoundaryInputs bin;
        bin.face = {FaceType::inflow, FaceType::periodic, FaceType::surface_layer, FaceType::outflow, FaceType::periodic,
                    FaceType::slip_wall};
        const gust::SurfaceLayerOptions so;
        const gust::BuoyancyOptions buoy;
        const gust::BoundaryConditions bca(bin, a.geom, a.mesh, gust::CoriolisOptions{});
        const gust::BoundaryConditions bcb(bin, b.geom, b.mesh, gust::CoriolisOptions{});
        gust::PoissonOptions popt;
        auto sa = a.h->make_poisson(0, popt);
        auto sb = b.h->make_poisson(0, popt);
        gust::prolong_state(*a.h, lc, ca, 0, fa, bca, &so, buoy, *sa, amrex::Real(1e-6));
        gust::prolong_state(*b.h, lc, cb, 0, fb, bcb, &so, buoy, *sb, amrex::Real(1e-6));
        EXPECT_EQ(diffs(fa, fb), 0) << "prolong_state, " << lc << " -> 0";
    }
}

// ---------------------------------------------------------------------------
// restrict_state and prolong_fields set the target State's
// node heights (ghost nodes included) to its level's, bit for bit, whatever
// they held before. Over the hill, where the levels' grounds differ.
TEST (Transfer, TargetNodeHeightsComeFromTheHierarchy)
{
    TCase c = case_inputs(false, true, true);
    c.make();
    const int ng = gust::ngrow_nodal;
    const std::pair<int, int> pairs[] = {{0, 1}, {1, 2}, {0, 2}};
    for (const auto& [lf, lc] : pairs) {
        const std::string where = std::to_string(lf) + " <-> " + std::to_string(lc);
        {
            gust::State qf, qc;
            define_level(c, lf, qf);
            qc.define(c.h->boxArray(lc), c.h->DistributionMap(lc));
            qc.z_phys_nd.setVal(-1.0);
            smooth_state(c, lf, qf);
            gust::restrict_state(*c.h, lf, qf, lc, qc);
            EXPECT_EQ(gust::test::count_bit_differences(qc.z_phys_nd, c.h->z_phys_nd(lc), ng), 0)
                << "restrict_state, " << where;
        }
        {
            gust::State qf, qc;
            define_level(c, lc, qc);
            qf.define(c.h->boxArray(lf), c.h->DistributionMap(lf));
            qf.z_phys_nd.setVal(-1.0);
            smooth_state(c, lc, qc);
            gust::prolong_fields(*c.h, lc, qc, lf, qf);
            EXPECT_EQ(gust::test::count_bit_differences(qf.z_phys_nd, c.h->z_phys_nd(lf), ng), 0)
                << "prolong_fields, " << where;
        }
    }
}

// ---------------------------------------------------------------------------
// A non-periodic direction one coarse cell wide. 32 x 4
// cells, periodic x, walls in y, three levels: level 2 is one cell wide in
// y, so the field is constant along y there and the interpolation reads
// that cell only (no second point to extrapolate with). A tent in x,
// constant in y, is reproduced at every fine point of each field, from
// level 2 to levels 1 and 0. The coarse image's ghost cells outside the
// domain are never filled: under amrex.init_snan (the unit-test flags) a
// read of them gives NaN, even with a zero weight.
TEST (Transfer, OneCellDirectionIsConstant)
{
    const double rel = tol(1e-13, 2e-6);
    TCase c;
    c.in.mesh.n_cell = amrex::IntVect(32, 4, 10);
    c.in.mesh.prob_lo = {0, 0, 0};
    c.in.mesh.prob_hi = {3200, 400, 1000};
    c.in.mesh.is_periodic = {1, 0, 0};
    c.in.mesh.max_grid_size = {12, 12};
    c.in.mesh.dz0 = 40;
    c.in.mesh.hierarchy_levels = 3;
    c.in.terrain.hierarchy_passes = 1;
    c.make();
    const int lc = 2;
    ASSERT_EQ(c.h->geom(lc).Domain().length(1), 1);
    for (const int lf : {0, 1}) {
        const amrex::Geometry& gc = c.h->geom(lc);
        gust::State qf, qc;
        define_level(c, lf, qf);
        define_level(c, lc, qc);
        amrex::MultiFab* cf[6] = {&qc.u, &qc.v, &qc.w, &qc.theta, &qc.tke, &qc.p};
        amrex::MultiFab* ff[6] = {&qf.u, &qf.v, &qf.w, &qf.theta, &qf.tke, &qf.p};
        const char* names[6] = {"u", "v", "w", "theta", "k", "p"};
        std::vector<std::function<double(double, double, int)>> fn(6);
        for (int n = 0; n < 6; ++n) {
            const double L = static_cast<double>(gc.ProbLength(0));
            const double D = static_cast<double>(gc.CellSize(0));
            const bool nodal = cf[n]->ixType().nodeCentered(0);
            const double xa = (gc.Domain().length(0) / 4 + (nodal ? 0.0 : 0.5)) * D;
            const double a = 1.0 + n, cc = 0.02 - 0.003 * n, g = 0.7 * n - 1;
            fn[n] = [=] (double x, double, int k) {
                double m = std::fmod(x - xa, L);
                if (m < 0) { m += L; }
                return a * (1.0 + cc * std::abs(m - L / 2)) + g * k;
            };
            fill_xy(*cf[n], gc, fn[n]);
            cf[n]->setBndry(amrex::Real(1e30));
            ff[n]->setVal(amrex::Real(-1e30));
        }
        gust::prolong_fields(*c.h, lc, qc, lf, qf);
        const amrex::Geometry& gf = c.h->geom(lf);
        for (int n = 0; n < 6; ++n) {
            const HostField h = gather(*ff[n]);
            double err = 0, scale = 0;
            const amrex::Box& d = h.dom;
            for (int k = d.smallEnd(2); k <= d.bigEnd(2); ++k) {
                for (int j = d.smallEnd(1); j <= d.bigEnd(1); ++j) {
                    for (int i = d.smallEnd(0); i <= d.bigEnd(0); ++i) {
                        const double ex = fn[n](xpos(gf, *ff[n], i, 0), xpos(gf, *ff[n], j, 1), k);
                        const double e = std::abs(h(i, j, k) - ex);
                        err = std::isfinite(e) ? std::max(err, e) : std::numeric_limits<double>::infinity();
                        scale = std::max(scale, std::abs(ex));
                    }
                }
            }
            EXPECT_LT(err, rel * scale) << names[n] << ", " << lc << " -> " << lf;
        }
    }
}

// ---------------------------------------------------------------------------
// A level pair that is not fine < coarse aborts, naming the
// rule, in restrict_state and prolong_fields.
TEST (Transfer, WrongLevelPairsAbort)
{
    TCase c = case_inputs(false, true, false);
    c.make();
    const std::pair<int, int> bad[] = {{1, 1}, {2, 1}, {-1, 1}, {0, 3}};
    for (const auto& [lf, lc] : bad) {
        const std::string where = std::to_string(lf) + ", " + std::to_string(lc);
        const int sf = std::clamp(lf, 0, 2), sc = std::clamp(lc, 0, 2);
        gust::State qf, qc;
        define_level(c, sf, qf);
        define_level(c, sc, qc);
        smooth_state(c, sf, qf);
        smooth_state(c, sc, qc);
        try {
            gust::restrict_state(*c.h, lf, qf, lc, qc);
            ADD_FAILURE() << "restrict_state: no abort for " << where;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("the levels must satisfy"), std::string::npos) << e.what();
        }
        try {
            gust::prolong_fields(*c.h, lc, qc, lf, qf);
            ADD_FAILURE() << "prolong_fields: no abort for " << where;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("the levels must satisfy"), std::string::npos) << e.what();
        }
    }
}
