// The first-order implicit-side operator (sub-phase 1.5): UT-1.5.4 (diagonal
// dominance), and its consistency with the explicit operators (beyond the
// plan): first-order upwind convection in the bounded form is A1 exactly, and
// on a flat mesh the diffusion operators have no cross terms, so their
// principal part is all of them.
#include "Convection.H"
#include "Diffusion.H"
#include "ImplicitOperator.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace {

using gust::Limiter;
using gust::Loc;
using gust::PoissonFace;
using gust::ptest::PCase;
using gust::trtest::tol;

const Loc kLocs[] = {Loc::cell, Loc::xface, Loc::yface, Loc::zface};

const char* loc_name (Loc loc)
{
    return loc == Loc::cell ? "cell" : loc == Loc::xface ? "u" : loc == Loc::yface ? "v" : "w";
}

amrex::IntVect itype (Loc loc)
{
    const int t = static_cast<int>(loc);
    return {t & 1, (t >> 1) & 1, (t >> 2) & 1};
}

// (32, 24, 16) cells over the cosine hill of height H0 (0: flat), geometric layers.
PCase make_mesh (double H0, int max_grid = 16)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = gust::ptest::make_case(32, 24, 16, PoissonFace::periodic, PoissonFace::periodic, ms, 0, true);
    if (max_grid == c.in.mesh.max_grid_size[0]) { return c; }
    PCase d;
    d.in = c.in;
    d.in.mesh.max_grid_size = {max_grid, max_grid};
    d.make([ms] (double x, double y) { return ms.h(x, y); });
    return d;
}

void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, double zero_fraction,
                  const amrex::Geometry& geom)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi), u01(0, 1);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const double v = dist(gen);
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(u01(gen) < zero_fraction ? 0.0 : v);
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    mf.FillBoundary(geom.periodicity());
}

struct Fields
{
    amrex::Array<amrex::MultiFab, 3> F;
    amrex::MultiFab nu;
};

Fields random_fields (const PCase& c, double fmax, double numax, unsigned seed)
{
    Fields f;
    for (int d = 0; d < 3; ++d) {
        f.F[d].define(amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d)), c.dm, 1, 1);
        fill_random(f.F[d], -fmax, fmax, seed + d, 0.0, c.geom);
    }
    f.nu.define(c.ba, c.dm, 1, 1);
    fill_random(f.nu, 0, numax, seed + 5, numax > 0 ? 0.2 : 1.0, c.geom);
    return f;
}

// UT-1.5.4. For random transport fluxes (+-20 m/s, not divergence-free) and a
// random diffusivity (0 to 50 m^2/s, a fifth of the cells 0) on the hill mesh,
// every row of A1 has non-negative neighbour weights and a_P equal to their
// sum: weakly diagonally dominant with non-positive off-diagonals.
TEST(ImplicitOperator, DiagonallyDominantForRandomWindAndNu)
{
    PCase c = make_mesh(120);
    Fields f = random_fields(c, 20, 50, 31);
    for (Loc loc : kLocs) {
        amrex::MultiFab A;
        gust::implicit_coefficients(loc, {&f.F[0], &f.F[1], &f.F[2]}, f.nu, c.z, c.mesh, c.geom, A);
        amrex::Long bad_sign = 0, bad_sum = 0, rows = 0;
        double max_rel = 0;
        for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            amrex::FArrayBox fa(vb, gust::n_implicit_coefs, amrex::The_Pinned_Arena());
            fa.copy<amrex::RunOn::Device>(A[mfi], vb);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                double sum = 0;
                for (int n = gust::coef_XM; n <= gust::coef_ZP; ++n) {
                    const double a = fa(iv, n);
                    if (!(a >= 0)) { ++bad_sign; }
                    sum += a;
                }
                const double ap = fa(iv, gust::coef_P);
                const double rel = std::abs(ap - sum) / std::max(ap, 1e-300);
                max_rel = std::max(max_rel, rel);
                if (!(ap > 0) || rel > tol(1e-14, 1e-6)) { ++bad_sum; }
                ++rows;
            });
        }
        amrex::ParallelDescriptor::ReduceLongSum(bad_sign);
        amrex::ParallelDescriptor::ReduceLongSum(bad_sum);
        amrex::ParallelDescriptor::ReduceLongSum(rows);
        EXPECT_GT(rows, 0);
        EXPECT_EQ(bad_sign, 0) << loc_name(loc) << ": negative neighbour weights";
        EXPECT_EQ(bad_sum, 0) << loc_name(loc) << ": a_P is not the sum of the neighbour weights (max rel " << max_rel << ")";
    }
}

// Max |a - b| over the valid points, relative to max |a|.
double rel_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::LinComb(d, 1.0, a, 0, -1.0, b, 0, 0, 1, 0);
    return static_cast<double>(d.norm0()) / static_cast<double>(a.norm0());
}

// Convection in A1 is the bounded first-order upwind operator: for random q
// and transport on the hill mesh, convection(upwind, bounded) equals A1 q
// (diffusivity 0) at every point, for all four locations.
TEST(ImplicitOperator, ConvectionPartIsBoundedUpwind)
{
    PCase c = make_mesh(120);
    Fields f = random_fields(c, 20, 0, 41);
    for (Loc loc : kLocs) {
        const amrex::BoxArray ba = amrex::convert(c.ba, itype(loc));
        amrex::MultiFab q(ba, c.dm, 1, 2), Cq(ba, c.dm, 1, 0), Aq(ba, c.dm, 1, 0), A;
        fill_random(q, -1, 2, 51, 0.0, c.geom);
        gust::ConvectionOptions opt;
        opt.limiter = Limiter::upwind;
        opt.bounded = true;
        gust::convection(loc, q, {&f.F[0], &f.F[1], &f.F[2]}, opt, c.mesh, c.geom, Cq);
        gust::implicit_coefficients(loc, {&f.F[0], &f.F[1], &f.F[2]}, f.nu, c.z, c.mesh, c.geom, A);
        gust::apply_implicit(loc, A, q, c.mesh, Aq);
        EXPECT_LE(rel_diff(Cq, Aq), tol(1e-12, 1e-5)) << loc_name(loc);
    }
}

// On a flat (stretched) mesh the metric cross terms vanish. Then minus the
// scalar diffusion equals A1 q (transport 0) for random q and gamma; and the
// momentum diffusion of a field with one non-zero component u_c equals -A1 u_c
// for that component: tau_cj = nu (d_j u_c + d_c u_j) reduces to
// nu (1 + delta_cj) d_j u_c, the principal part with the own-direction factor.
TEST(ImplicitOperator, DiffusionPartIsTheFlatOperator)
{
    PCase c = make_mesh(0);
    Fields f = random_fields(c, 0, 50, 61);
    for (int d = 0; d < 3; ++d) { f.F[d].setVal(0.0); }

    amrex::MultiFab q(c.ba, c.dm, 1, 1), Dq(c.ba, c.dm, 1, 0), Aq(c.ba, c.dm, 1, 0), A;
    fill_random(q, -1, 2, 71, 0.0, c.geom);
    gust::scalar_diffusion(q, f.nu, c.z, c.mesh, c.geom, Dq);
    gust::implicit_coefficients(Loc::cell, {&f.F[0], &f.F[1], &f.F[2]}, f.nu, c.z, c.mesh, c.geom, A);
    gust::apply_implicit(Loc::cell, A, q, c.mesh, Aq);
    Dq.mult(-1.0);
    EXPECT_LE(rel_diff(Dq, Aq), tol(1e-12, 1e-5)) << "scalar";

    for (int cc = 0; cc < 3; ++cc) {
        const Loc loc = gust::velocity_loc(cc);
        amrex::Array<amrex::MultiFab, 3> vel, out;
        for (int d = 0; d < 3; ++d) {
            const amrex::BoxArray ba = amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d));
            vel[d].define(ba, c.dm, 1, 2);
            out[d].define(ba, c.dm, 1, 0);
            vel[d].setVal(0.0);
        }
        fill_random(vel[cc], -1, 2, 81 + cc, 0.0, c.geom);
        gust::momentum_diffusion({&vel[0], &vel[1], &vel[2]}, f.nu, c.z, c.mesh, c.geom, {&out[0], &out[1], &out[2]});
        amrex::MultiFab Ac, Au(out[cc].boxArray(), c.dm, 1, 0);
        gust::implicit_coefficients(loc, {&f.F[0], &f.F[1], &f.F[2]}, f.nu, c.z, c.mesh, c.geom, Ac);
        gust::apply_implicit(loc, Ac, vel[cc], c.mesh, Au);
        out[cc].mult(-1.0);
        EXPECT_LE(rel_diff(out[cc], Au), tol(1e-12, 1e-5)) << loc_name(loc);
    }
}

// Fill mf with f(index along d), the same in every column: random values per
// index, ghosts included.
void fill_along (amrex::MultiFab& mf, int d, unsigned seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1, 2);
    std::vector<double> f(512);
    for (auto& v : f) { v = dist(gen); }
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox fab(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const int n = (d == 0 ? i : d == 1 ? j : k) + 8;
            fab(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(f[static_cast<std::size_t>(n)]);
        });
        mf[mfi].copy<amrex::RunOn::Device>(fab, b);
        amrex::Gpu::streamSynchronize();
    }
}

// Every weight of A1 is the explicit operator's flux for a field that varies
// along one index direction only, over the hill: for q = f(i_d) (one velocity
// component at a time, the others zero) the cross-metric parts of the face
// fluxes in direction d vanish, so
//   G_d(low face of P) = a_{d-}(P) J(P) dxi_d (q(P) - q(P - e_d)).
// This pins the own-direction factors 1 + delta_cd and the z-face factor
// |n|^2 + n_c^2 on sloping faces, which the flat test cannot see.
TEST(ImplicitOperator, WeightsAreOneDimensionalFluxes)
{
    PCase c = make_mesh(120);
    Fields f = random_fields(c, 0, 50, 91);
    for (int d = 0; d < 3; ++d) { f.F[d].setVal(0.0); }
    const auto dx = c.geom.CellSizeArray();
    const int minus[3] = {gust::coef_XM, gust::coef_YM, gust::coef_ZM};
    for (int comp = -1; comp < 3; ++comp) {
        const Loc loc = comp < 0 ? Loc::cell : gust::velocity_loc(comp);
        const int t = static_cast<int>(loc);
        amrex::MultiFab A;
        gust::implicit_coefficients(loc, {&f.F[0], &f.F[1], &f.F[2]}, f.nu, c.z, c.mesh, c.geom, A);
        for (int d = 0; d < 3; ++d) {
            amrex::Array<amrex::MultiFab, 3> G;
            amrex::MultiFab q;
            if (comp < 0) {
                q.define(c.ba, c.dm, 1, 1);
                fill_along(q, d, 100 + d);
                gust::scalar_diffusive_fluxes(q, f.nu, c.z, c.mesh, c.geom, G);
            } else {
                amrex::Array<amrex::MultiFab, 3> vel;
                amrex::Array<amrex::Array<amrex::MultiFab, 3>, 3> GG;
                for (int e = 0; e < 3; ++e) {
                    vel[e].define(amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(e)), c.dm, 1, 2);
                    vel[e].setVal(0.0);
                }
                fill_along(vel[comp], d, 200 + 10 * comp + d);
                gust::momentum_diffusive_fluxes({&vel[0], &vel[1], &vel[2]}, f.nu, c.z, c.mesh, c.geom, GG);
                for (int e = 0; e < 3; ++e) { G[e] = std::move(GG[comp][e]); }
                q = std::move(vel[comp]);
            }
            double max_diff = 0, max_flux = 0;
            for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                const amrex::Box qb = amrex::grow(vb, 1);
                const amrex::Box gb = amrex::grow(G[d][mfi].box(), 0);
                amrex::FArrayBox fa(vb, gust::n_implicit_coefs, amrex::The_Pinned_Arena()), fq(qb, 1, amrex::The_Pinned_Arena()),
                    fg(gb, 1, amrex::The_Pinned_Arena()), fj(c.mesh.detJ[mfi].box(), 1, amrex::The_Pinned_Arena());
                fa.copy<amrex::RunOn::Device>(A[mfi], vb);
                fq.copy<amrex::RunOn::Device>(q[mfi], qb);
                fg.copy<amrex::RunOn::Device>(G[d][mfi], gb);
                fj.copy<amrex::RunOn::Device>(c.mesh.detJ[mfi], fj.box());
                amrex::Gpu::streamSynchronize();
                amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                    const amrex::IntVect iv(i, j, k);
                    amrex::IntVect e(0);
                    e[d] = 1;
                    // The low flux location: the same index where loc is cell-centred in d, one below where nodal.
                    const amrex::IntVect lo = gust::stencil::nodal(t, d) ? iv - e : iv;
                    double J = 0;
                    int n = 0;
                    for (int cc = 0; cc < 3; ++cc) {
                        if (gust::stencil::nodal(t, cc)) {
                            amrex::IntVect o(0);
                            o[cc] = 1;
                            J += fj(iv - o) + fj(iv);
                            n += 2;
                        }
                    }
                    J = (n == 0) ? static_cast<double>(fj(iv)) : J / n;
                    const double pred = fa(iv, minus[d]) * J * dx[d] * (static_cast<double>(fq(iv)) - static_cast<double>(fq(iv - e)));
                    max_diff = std::max(max_diff, std::abs(pred - static_cast<double>(fg(lo))));
                    max_flux = std::max(max_flux, std::abs(static_cast<double>(fg(lo))));
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(max_diff);
            amrex::ParallelDescriptor::ReduceRealMax(max_flux);
            EXPECT_GT(max_flux, 0.0);
            EXPECT_LE(max_diff, tol(1e-12, 1e-5) * max_flux) << loc_name(loc) << " direction " << d;
        }
    }
}

// The operators give the same answer on any box layout: 16^2-column boxes
// and 8^2-column boxes (4x as many box faces) for MUSCL convection, momentum
// and scalar diffusion and A1, all on the manufactured fields. Agreement is
// to round-off, not bit for bit: where a compiler contracts multiply-adds
// (GCC on arm64), a point inside a vectorised loop and the same point in a
// loop remainder round differently, and box edges move the remainders.
TEST(ImplicitOperator, IndependentOfBoxLayout)
{
    auto run = [] (int max_grid, amrex::Vector<amrex::MultiFab>& res) {
        PCase c = make_mesh(120, max_grid);
        const gust::trtest::Flow fl = gust::trtest::make_flow(c, 120);
        amrex::Array<amrex::MultiFab, 3> F, vel;
        gust::trtest::fill_transport(F, c, fl);
        amrex::MultiFab nu(c.ba, c.dm, 1, 1), q(c.ba, c.dm, 1, 2);
        gust::trtest::fill_loc(nu, Loc::cell, c, [&] (double x, double y, double z) { return fl.nu(x, y, z); });
        gust::trtest::fill_loc(q, Loc::cell, c, [&] (double x, double y, double z) { return fl.q(x, y, z); });
        for (int d = 0; d < 3; ++d) {
            vel[d].define(amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d)), c.dm, 1, 2);
            gust::trtest::fill_loc(vel[d], gust::velocity_loc(d), c,
                                       [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
        }
        res.clear();
        res.resize(8);
        gust::ConvectionOptions opt;
        for (int d = 0; d < 3; ++d) {
            res[d].define(vel[d].boxArray(), c.dm, 1, 0);
            gust::convection(gust::velocity_loc(d), vel[d], {&F[0], &F[1], &F[2]}, opt, c.mesh, c.geom, res[d]);
            res[3 + d].define(vel[d].boxArray(), c.dm, 1, 0);
        }
        gust::momentum_diffusion({&vel[0], &vel[1], &vel[2]}, nu, c.z, c.mesh, c.geom, {&res[3], &res[4], &res[5]});
        res[6].define(c.ba, c.dm, 1, 0);
        gust::scalar_diffusion(q, nu, c.z, c.mesh, c.geom, res[6]);
        amrex::MultiFab A;
        gust::implicit_coefficients(Loc::xface, {&F[0], &F[1], &F[2]}, nu, c.z, c.mesh, c.geom, A);
        res[7].define(A.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        amrex::MultiFab::Copy(res[7], A, 0, 0, gust::n_implicit_coefs, 0);
    };
    amrex::Vector<amrex::MultiFab> a, b;
    run(16, a);
    run(8, b);
    for (int n = 0; n < 8; ++n) {
        amrex::MultiFab bb(a[n].boxArray(), a[n].DistributionMap(), a[n].nComp(), 0);
        bb.ParallelCopy(b[n], 0, 0, a[n].nComp());
        amrex::MultiFab::Subtract(bb, a[n], 0, 0, a[n].nComp(), 0);
        for (int comp = 0; comp < a[n].nComp(); ++comp) {
            EXPECT_LE(static_cast<double>(bb.norm0(comp)), tol(1e-13, 1e-5) * static_cast<double>(a[n].norm0(comp)))
                << "result " << n << " comp " << comp;
        }
        // Within one layout, the two boxes that share a face point (or its
        // periodic image) agree there.
        for (auto* r : {&a[n], &b[n]}) {
            if (r->ixType().cellCentered()) { continue; }
            amrex::MultiFab synced(r->boxArray(), r->DistributionMap(), r->nComp(), 0);
            amrex::MultiFab::Copy(synced, *r, 0, 0, r->nComp(), 0);
            synced.OverrideSync(amrex::Periodicity(amrex::IntVect(32, 24, 0)));
            amrex::MultiFab::Subtract(synced, *r, 0, 0, r->nComp(), 0);
            for (int comp = 0; comp < r->nComp(); ++comp) {
                EXPECT_LE(static_cast<double>(synced.norm0(comp)), tol(1e-13, 1e-5) * static_cast<double>(r->norm0(comp)))
                    << "shared points, result " << n;
            }
        }
    }
}

} // namespace
