// Sub-phase 2.3: the horizontal-only grid hierarchy (MeshHierarchy.H).
//   UT-2.3.1 z-levels bit-identical on every level
//   UT-2.3.2 ground smoothing passes p (2^l - 1); constant terrain unchanged
//   UT-2.3.3 free-stream preservation on every level
//   UT-2.3.4 flat terrain: every level gives the same column (Poisson and SIMPLEC)
//   UT-2.3.5 invalid layouts abort, naming the key
//   UT-2.3.6 a level gathered onto one box equals the distributed level
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Inputs.H"
#include "Mesh.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "Simplec.H"
#include "SmoothSurface.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double tol (double dbl, double flt) { return kDouble ? dbl : flt; }

// The run's level 0 (as App.cpp lays it out) and the hierarchy built on it.
struct HCase
{
    ablsteady::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab z;
    ablsteady::Mesh mesh;
    std::unique_ptr<ablsteady::MeshHierarchy> h;

    void layout ()
    {
        geom = ablsteady::make_geometry(in.mesh);
        ba = ablsteady::hierarchy_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        z.define(amrex::convert(ba, amrex::IntVect(1)), dm, 1, ablsteady::ngrow_nodal);
        z.setVal(0.0);
        mesh.define(ba, dm);
    }
    // Level 0 from the terrain inputs, then the hierarchy.
    void make (int gather_from = 0)
    {
        layout();
        mesh.build(in, geom, z);
        h = std::make_unique<ablsteady::MeshHierarchy>(in, geom, z, mesh, gather_from);
    }
};

// A round hill on a stretched mesh, non-periodic in x and periodic in y, cut
// into several boxes on every level (32 x 16 cells, boxes of 8: 8 boxes).
HCase hill3d (int smoothing, int levels = 3)
{
    HCase c;
    c.in.mesh.n_cell = amrex::IntVect(32, 16, 12);
    c.in.mesh.prob_lo = {0, 0, 0};
    c.in.mesh.prob_hi = {3200, 1600, 1200};
    c.in.mesh.is_periodic = {0, 1, 0};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.dz0 = 40;
    c.in.mesh.hierarchy_levels = levels;
    c.in.terrain.type = "agnesi";
    c.in.terrain.hill_shape = "round";
    c.in.terrain.hill_height = 150;
    c.in.terrain.hill_half_width = 400;
    c.in.terrain.smoothing = smoothing;
    c.in.terrain.hierarchy_passes = 1;
    return c;
}

// The k = 0 slab of z (ground), every node of every fab, ghosts included,
// as host values keyed by (i, j) of the fab: max |a - b| over the slab of
// two MultiFabs on one layout, and whether they are bit-identical there.
struct SlabDiff { double max_abs{0}; long bit_diffs{0}; };
SlabDiff ground_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    SlabDiff d;
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.fabbox();
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        amrex::FArrayBox ha(bx, 1, amrex::The_Pinned_Arena()), hb(bx, 1, amrex::The_Pinned_Arena());
        ha.copy<amrex::RunOn::Device>(a[mfi], bx);
        hb.copy<amrex::RunOn::Device>(b[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            d.max_abs = std::max(d.max_abs, std::abs(static_cast<double>(ha(iv) - hb(iv))));
            if (std::memcmp(&ha(iv), &hb(iv), sizeof(amrex::Real)) != 0) { ++d.bit_diffs; }
        });
    }
    amrex::ParallelAllReduce::Max(d.max_abs, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(d.bit_diffs, amrex::ParallelDescriptor::Communicator());
    return d;
}

// The ground of level l rebuilt independently on the host: level 0's ground
// at node (r I, r J), then n passes of h + beta (sum of the four neighbours -
// 4 h), neighbours wrapped in a periodic direction and clamped to the edge
// node at a non-periodic face, in double precision. Returns the nodal
// (nx + 1) x (ny + 1) array, i fastest.
std::vector<double> host_ground (const HCase& c, int l, int n)
{
    const int r = 1 << l;
    const int NX = c.in.mesh.n_cell[0], NY = c.in.mesh.n_cell[1];
    const int nx = NX / r, ny = NY / r;
    // Level 0's ground on every rank: the mean over the boxes that hold each
    // node (a node shared by boxes has one value, so the mean is that value).
    std::vector<double> g0((NX + 1) * (NY + 1), 0.0);
    std::vector<double> local((NX + 1) * (NY + 1), 0.0), nloc((NX + 1) * (NY + 1), 0.0);
    for (amrex::MFIter mfi(c.z, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(c.z[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
            const int m = j * (NX + 1) + i;
            local[m] += static_cast<double>(f(amrex::IntVect(i, j, 0)));
            nloc[m] += 1;
        });
    }
    amrex::ParallelAllReduce::Sum(local.data(), static_cast<int>(local.size()), amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(nloc.data(), static_cast<int>(nloc.size()), amrex::ParallelDescriptor::Communicator());
    for (std::size_t m = 0; m < g0.size(); ++m) { g0[m] = local[m] / nloc[m]; }

    std::vector<double> h((nx + 1) * (ny + 1));
    for (int J = 0; J <= ny; ++J) {
        for (int I = 0; I <= nx; ++I) { h[J * (nx + 1) + I] = g0[(r * J) * (NX + 1) + r * I]; }
    }
    const bool px = c.in.mesh.is_periodic[0] != 0, py = c.in.mesh.is_periodic[1] != 0;
    // Periodic: node n is node 0, so the neighbours of node 0 are 1 and n - 1.
    auto at = [&] (const std::vector<double>& a, int I, int J) {
        if (px) { I = (I < 0) ? I + nx : (I > nx ? I - nx : I); } else { I = std::max(0, std::min(nx, I)); }
        if (py) { J = (J < 0) ? J + ny : (J > ny ? J - ny : J); } else { J = std::max(0, std::min(ny, J)); }
        return a[J * (nx + 1) + I];
    };
    const double beta = ablsteady::kHierarchyBeta;
    for (int p = 0; p < n; ++p) {
        const std::vector<double> o = h;
        for (int J = 0; J <= ny; ++J) {
            for (int I = 0; I <= nx; ++I) {
                h[J * (nx + 1) + I] = at(o, I, J) + beta * (at(o, I + 1, J) + at(o, I - 1, J) + at(o, I, J + 1)
                                                            + at(o, I, J - 1) - 4 * at(o, I, J));
            }
        }
    }
    return h;
}

// max |ground of level l - ref| over the valid nodes.
double max_ground_error (const amrex::MultiFab& z, const std::vector<double>& ref, int nx)
{
    double e = 0;
    for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.setSmall(2, 0);
        bx.setBig(2, 0);
        amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(z[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
            e = std::max(e, std::abs(static_cast<double>(f(amrex::IntVect(i, j, 0))) - ref[j * (nx + 1) + i]));
        });
    }
    amrex::ParallelAllReduce::Max(e, amrex::ParallelDescriptor::Communicator());
    return e;
}

// The column of mf at (i, j) = (0, 0), every k of its valid range, on every rank.
std::vector<double> column0 (const amrex::MultiFab& mf)
{
    const amrex::Box nb = amrex::convert(mf.boxArray().minimalBox(), mf.ixType());
    const int klo = nb.smallEnd(2), n = nb.length(2);
    std::vector<double> col(n, 0.0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (vb.smallEnd(0) != 0 || vb.smallEnd(1) != 0) { continue; }
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(mf[mfi], vb);
        amrex::Gpu::streamSynchronize();
        for (int k = klo; k < klo + n; ++k) { col[k - klo] = static_cast<double>(f(amrex::IntVect(0, 0, k))); }
    }
    amrex::ParallelAllReduce::Sum(col.data(), n, amrex::ParallelDescriptor::Communicator());
    return col;
}

// max over the valid points of |mf(i, j, k) - col(k)|.
double max_dev_from_column (const amrex::MultiFab& mf, const std::vector<double>& col)
{
    const int klo = amrex::convert(mf.boxArray().minimalBox(), mf.ixType()).smallEnd(2);
    double e = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(mf[mfi], vb);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            e = std::max(e, std::abs(static_cast<double>(f(amrex::IntVect(i, j, k))) - col[k - klo]));
        });
    }
    amrex::ParallelAllReduce::Max(e, amrex::ParallelDescriptor::Communicator());
    return e;
}

double max_abs (const std::vector<double>& v)
{
    double m = 0;
    for (const double x : v) { m = std::max(m, std::abs(x)); }
    return m;
}

// Fill every point of mf (ghosts included) with prof(z) of its k: z the
// cell-centre height (nodes zl) for a cell or x/y-face field, the node
// height for a z-face field; ghost rows take the nearest valid row's value.
template <class Fn>
void fill_profile (amrex::MultiFab& mf, const amrex::Vector<amrex::Real>& zl, const Fn& prof)
{
    const int nz = static_cast<int>(zl.size()) - 1;
    const bool zface = mf.ixType().nodeCentered(2);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            double zz = 0;
            if (zface) {
                zz = static_cast<double>(zl[std::max(0, std::min(nz, k))]);
            } else {
                const int kc = std::max(0, std::min(nz - 1, k));
                zz = 0.5 * static_cast<double>(zl[kc] + zl[kc + 1]);
            }
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(prof(zz));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

// A conventionally neutral column: the k-L closure, the surface layer,
// Coriolis with a geostrophic wind (10, 0) m/s, a capping inversion from 600
// to 700 m with a lapse rate above, buoyancy, a slip top.
struct ColumnModel
{
    ablsteady::CoriolisOptions cor;
    ablsteady::ReferenceProfile ref;
    ablsteady::BoundaryInputs bin;
    ablsteady::SurfaceLayerOptions surf;
    ablsteady::ClosureOptions copt;
    ablsteady::SimplecOptions sopt;
    ablsteady::PoissonOptions popt;
    ColumnModel ()
    {
        cor.enabled = true;
        cor.latitude = 45;
        cor.geostrophic = true;
        cor.ug = 10;
        cor.vg = 0;
        ref.z = {0, 600, 700, 3000};
        ref.theta = {300, 300, 305, 312};
        bin.face = {ablsteady::FaceType::periodic, ablsteady::FaceType::periodic, ablsteady::FaceType::surface_layer,
                    ablsteady::FaceType::periodic, ablsteady::FaceType::periodic, ablsteady::FaceType::slip_wall};
        bin.theta_ref = ref;
        bin.zhi_theta_grad = amrex::Real(0.003);
        popt.backend = "gmres_fft";   // exact column symmetry on flat ground (see test (a) below)
    }
};

double u_prof (double z) { return 8.0 * std::log1p(z / 0.1) / std::log1p(1500.0 / 0.1); }
double v_prof (double z) { return -1.5 * std::sin(3.14159265358979 * z / 1500.0); }
double zero_prof (double) { return 0.0; }
double theta_prof (double z)
{
    return z < 600 ? 300.0 : (z < 700 ? 300.0 + 5.0 * (z - 600) / 100 : 305.0 + 7.0 * (z - 700) / 2300);
}
double k_prof (double z) { return 0.4 * std::exp(-z / 400.0) + 1e-3; }

// On level l of c: the column state, then `iterations` SIMPLEC iterations
// with a uniform pseudo-time step of 2 s, then the residual R of the result
// (ghosts filled by the boundary conditions first).
void column_level (const HCase& c, int l, const ColumnModel& cm, int iterations, ablsteady::State& q,
                   ablsteady::ResidualFields& R)
{
    const amrex::Geometry& g = c.h->geom(l);
    const ablsteady::Mesh& m = c.h->mesh(l);
    const amrex::BoxArray& ba = c.h->boxArray(l);
    const amrex::DistributionMapping& dm = c.h->DistributionMap(l);
    q.define(ba, dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.h->z_phys_nd(l), 0, 0, 1, ablsteady::ngrow_nodal);
    const auto& zl = m.zlevels();
    fill_profile(q.u, zl, u_prof);
    fill_profile(q.v, zl, v_prof);
    fill_profile(q.w, zl, zero_prof);
    fill_profile(q.p, zl, zero_prof);
    fill_profile(q.theta, zl, theta_prof);
    fill_profile(q.tke, zl, k_prof);

    ablsteady::Closure closure(cm.copt);
    amrex::MultiFab dist;
    ablsteady::height_above_ground(q.z_phys_nd, dist);
    const ablsteady::BoundaryConditions bc(cm.bin, g, m, cm.cor);
    ablsteady::ResidualInputs res_in;
    res_in.closure = &closure;
    res_in.wall_distance = &dist;
    res_in.theta_ref = cm.ref;
    res_in.coriolis = cm.cor;
    res_in.surface = &cm.surf;
    res_in.bc = &bc;
    if (iterations > 0) {
        ablsteady::SimplecStage stage(cm.sopt, 2, res_in, &closure, cm.popt, m, g, c.h->z_phys_nd(l), cm.copt.tke_floor, ba, dm);
        stage.init(q);
        amrex::MultiFab dtau(ba, dm, 1, 0);
        dtau.setVal(2.0);
        for (int it = 0; it < iterations; ++it) { stage.iterate(q, dtau); }
    }
    bc.fill_ghosts(q, &cm.surf, res_in.buoyancy, false, nullptr);
    ablsteady::residual(q, res_in, m, g, R);
}

} // namespace

// ---------------------------------------------------------------------------
// UT-2.3.1. Every level has the z-levels of level 0, bit for bit (a stretched
// grid whose ratio is solved by Newton), the same nz and z extent, and on
// flat terrain its node heights are those levels exactly in every column.
TEST (MeshHierarchy, ZLevelsAreIdenticalOnEveryLevel)
{
    for (const bool flat : {true, false}) {
        HCase c = hill3d(1);
        if (flat) { c.in.terrain = ablsteady::TerrainInputs{}; c.in.terrain.hierarchy_passes = 2; }
        c.make();
        const auto& z0 = c.h->mesh(0).zlevels();
        ASSERT_EQ(c.h->num_levels(), 3);
        for (int l = 1; l < c.h->num_levels(); ++l) {
            const auto& zl = c.h->mesh(l).zlevels();
            ASSERT_EQ(zl.size(), z0.size());
            EXPECT_EQ(std::memcmp(zl.data(), z0.data(), z0.size() * sizeof(amrex::Real)), 0) << "level " << l;
            EXPECT_EQ(c.h->geom(l).Domain().length(2), c.geom.Domain().length(2));
            EXPECT_EQ(c.h->geom(l).ProbLo(2), c.geom.ProbLo(2));
            EXPECT_EQ(c.h->geom(l).ProbHi(2), c.geom.ProbHi(2));
            EXPECT_EQ(c.h->geom(l).Domain().length(0) * (1 << l), c.geom.Domain().length(0));
            EXPECT_EQ(c.h->geom(l).Domain().length(1) * (1 << l), c.geom.Domain().length(1));
            if (flat) {
                // Every valid node: z(i, j, k) == z_k.
                long nbad = 0;
                const auto& zn = c.h->z_phys_nd(l);
                for (amrex::MFIter mfi(zn, false); mfi.isValid(); ++mfi) {
                    const amrex::Box vb = mfi.validbox();
                    amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
                    f.copy<amrex::RunOn::Device>(zn[mfi], vb);
                    amrex::Gpu::streamSynchronize();
                    amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                        if (std::memcmp(&f(amrex::IntVect(i, j, k)), &z0[k], sizeof(amrex::Real)) != 0) { ++nbad; }
                    });
                }
                amrex::ParallelAllReduce::Sum(nbad, amrex::ParallelDescriptor::Communicator());
                EXPECT_EQ(nbad, 0) << "level " << l;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// UT-2.3.2. The pass counts are p (2^l - 1): 0 on level 0, pinned for p = 0
// to 3 and l = 0 to 4. Each coarse ground equals an independent host
// rebuild (injection, then that many passes of the STF stencil with beta =
// 1/8) and differs from one pass more or fewer; level 0's ground is the
// input's bit for bit; a constant ground is unchanged on every level, bit
// for bit, ghost nodes included.
TEST (MeshHierarchy, GroundSmoothingPasses)
{
    const int expect[4][5] = {{0, 0, 0, 0, 0}, {0, 1, 3, 7, 15}, {0, 2, 6, 14, 30}, {0, 3, 9, 21, 45}};
    for (int p = 0; p < 4; ++p) {
        for (int l = 0; l < 5; ++l) { EXPECT_EQ(ablsteady::ground_passes(l, p), expect[p][l]) << "p " << p << " l " << l; }
    }

    for (const int p : {0, 1, 2}) {
        HCase c = hill3d(1);
        c.in.terrain.hierarchy_passes = p;
        c.make();
        // Level 0 is the run's own mesh: not smoothed.
        EXPECT_EQ(&c.h->z_phys_nd(0), &c.z);
        EXPECT_EQ(c.h->ground_passes(0), 0);
        for (int l = 1; l < c.h->num_levels(); ++l) {
            const int n = c.h->ground_passes(l);
            EXPECT_EQ(n, p * ((1 << l) - 1));
            const int nx = c.in.mesh.n_cell[0] >> l;
            const double H = 150;
            const double e = max_ground_error(c.h->z_phys_nd(l), host_ground(c, l, n), nx);
            EXPECT_LT(e, tol(1e-12, 1e-4) * H) << "p " << p << " level " << l;
            // One pass more or fewer is visibly different (the count is pinned).
            if (n > 0) {
                EXPECT_GT(max_ground_error(c.h->z_phys_nd(l), host_ground(c, l, n - 1), nx), 1e-3) << "level " << l;
            }
            EXPECT_GT(max_ground_error(c.h->z_phys_nd(l), host_ground(c, l, n + 1), nx), 1e-3) << "level " << l;
        }
    }

    // Constant terrain: 37.5 m everywhere (exact in binary), non-periodic in x.
    HCase c = hill3d(1);
    c.in.terrain = ablsteady::TerrainInputs{};
    c.in.terrain.hierarchy_passes = 3;
    c.layout();
    c.z.setVal(37.5);
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    c.h = std::make_unique<ablsteady::MeshHierarchy>(c.in, c.geom, c.z, c.mesh);
    for (int l = 0; l < c.h->num_levels(); ++l) {
        amrex::MultiFab ref(c.h->z_phys_nd(l).boxArray(), c.h->z_phys_nd(l).DistributionMap(), 1, ablsteady::ngrow_nodal);
        ref.setVal(37.5);
        const SlabDiff d = ground_diff(c.h->z_phys_nd(l), ref);
        EXPECT_EQ(d.bit_diffs, 0) << "level " << l << ", max diff " << d.max_abs;
    }
}

// ---------------------------------------------------------------------------
// UT-2.3.3. Free-stream preservation on every level, for BTF, STF and
// Sullivan: the metric divergence of a uniform wind (10, 5) m/s is below
// 1e-12 of (|U| + |V|) / dx of the level.
TEST (MeshHierarchy, FreeStreamOnEveryLevel)
{
    for (const int sm : {0, 1, 2}) {
        HCase c = hill3d(sm);
        c.make();
        for (int l = 0; l < c.h->num_levels(); ++l) {
            const double scale = 15.0 * c.h->geom(l).InvCellSize(0);
            const double div = ablsteady::free_stream_divergence(c.h->mesh(l), c.h->geom(l), c.h->z_phys_nd(l), 10, 5);
            EXPECT_LT(div, tol(1e-12, 1e-4) * scale) << "smoothing " << sm << ", level " << l;
        }
        // The coarse grounds are not flat (the test is not vacuous).
        for (int l = 1; l < c.h->num_levels(); ++l) {
            EXPECT_GT(c.h->mesh(l).quality().max_ground_height, 50) << "level " << l;
        }
    }
}

// ---------------------------------------------------------------------------
// UT-2.3.4. Flat terrain: every level gives the same column.
//  (a) A Poisson problem whose right-hand side depends on z alone, solved
//      with each level's own solver: the same phi column on every level.
//  (b) The steady equations of a conventionally neutral column (the k-L
//      closure, the surface layer, Coriolis with a geostrophic wind, a
//      capping inversion with buoyancy): the full residual of the same
//      horizontally uniform state has the same u, v, w, p, theta and k
//      columns on every level, so every level has the same column solutions.
//  (c) Ten SIMPLEC iterations from that state keep every column of a level
//      alike, on every level.
// The iterates of (c) differ between levels: the pseudo-time step and the
// implicit diagonal hold horizontal terms (|u|/dx, nu/dx^2) that cancel only
// at the fixed point. So the steady equations (b) are compared, not runs to
// steady; a column with a capping inversion has no strict steady state to
// compare (the inversion keeps rising).
// Relative to the largest value: < 1e-12.
TEST (MeshHierarchy, FlatTerrainColumnsAgree)
{
    HCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 8, 24);
    c.in.mesh.prob_hi = {3200, 1600, 1500};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {8, 4};
    c.in.mesh.dz0 = 15;
    c.in.mesh.hierarchy_levels = 3;
    c.in.terrain.hierarchy_passes = 2;
    c.make();
    const int L = c.h->num_levels();
    const double rel = tol(1e-12, 1e-5);

    // (a) Poisson: rhs = cos(pi z / H) (Neumann ground and top, periodic sides).
    {
        std::vector<double> col_ref;
        for (int l = 0; l < L; ++l) {
            // Back end A: its FFT preconditioner is exact on flat ground, so the
            // columns agree to round-off; MLMG's red-black smoother breaks the
            // column symmetry at its tolerance (this test is of the hierarchy).
            ablsteady::PoissonOptions popt;
            popt.backend = "gmres_fft";
            auto solver = c.h->make_poisson(l, popt);
            amrex::MultiFab rhs(c.h->boxArray(l), c.h->DistributionMap(l), 1, 0);
            amrex::MultiFab phi(c.h->boxArray(l), c.h->DistributionMap(l), 1, 1);
            phi.setVal(0.0);
            fill_profile(rhs, c.h->mesh(l).zlevels(), [] (double z) { return std::cos(3.14159265358979 * z / 1500.0); });
            ablsteady::PoissonBC bc;
            bc.face = {ablsteady::PoissonFace::periodic, ablsteady::PoissonFace::periodic, ablsteady::PoissonFace::neumann,
                       ablsteady::PoissonFace::periodic, ablsteady::PoissonFace::periodic, ablsteady::PoissonFace::neumann};
            solver->solve(rhs, phi, ablsteady::FaceCoeffs{}, bc);
            const std::vector<double> col = column0(phi);
            if (l == 0) { col_ref = col; ASSERT_GT(max_abs(col), 1.0); }
            EXPECT_LT(max_dev_from_column(phi, col_ref), rel * max_abs(col_ref)) << "Poisson, level " << l;
        }
    }

    // (b) The steady equations: the full residual of the column state on
    // every level, (c) ten SIMPLEC iterations on every level.
    ColumnModel cm;
    std::vector<std::vector<double>> r_ref(6), q_ref(6);
    for (int l = 0; l < L; ++l) {
        ablsteady::State q;
        ablsteady::ResidualFields R;
        column_level(c, l, cm, 0, q, R);
        const amrex::MultiFab* rf[6] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
        const char* names[6] = {"u", "v", "w", "p", "theta", "k"};
        for (int n = 0; n < 6; ++n) {
            const std::vector<double> col = column0(*rf[n]);
            if (l == 0) { r_ref[n] = col; }
            // R_w and R_p (the divergence) are zero to round-off on a
            // horizontally uniform state with w = 0: scaled by R_u's size.
            const double scale = (n == 2 || n == 3) ? max_abs(r_ref[0]) : max_abs(r_ref[n]);
            EXPECT_LT(max_dev_from_column(*rf[n], r_ref[n]), rel * scale) << "R_" << names[n] << ", level " << l;
        }
        // Not vacuous: the state is not steady (Coriolis, the wall, the closure).
        if (l == 0) {
            for (const int n : {0, 1, 4, 5}) { EXPECT_GT(max_abs(r_ref[n]), 1e-6) << "R_" << names[n]; }
        }

        // (c) The stage: every column of a level stays alike (the iterates of
        // different levels differ in pseudo-time, through the horizontal
        // terms of each cell's implicit diagonal, not at the fixed point).
        ablsteady::State q10;
        ablsteady::ResidualFields R10;
        column_level(c, l, cm, 10, q10, R10);
        const amrex::MultiFab* f[6] = {&q10.u, &q10.v, &q10.w, &q10.p, &q10.theta, &q10.tke};
        for (int n = 0; n < 6; ++n) {
            const std::vector<double> col = column0(*f[n]);
            const double scale = (n == 2 || n == 3) ? 10.0 : max_abs(col);
            EXPECT_LT(max_dev_from_column(*f[n], col), rel * scale) << names[n] << " after 10 iterations, level " << l;
        }
    }
}

// ---------------------------------------------------------------------------
// UT-2.3.5. Invalid layouts and inputs abort, naming the key; valid ones pass.
TEST (MeshHierarchy, InvalidLayoutsAbort)
{
    auto msgs = [] (const ablsteady::MeshInputs& m, const ablsteady::TerrainInputs& t) {
        return ablsteady::check_hierarchy_inputs(m, t);
    };
    ablsteady::MeshInputs m;
    m.n_cell = amrex::IntVect(32, 16, 8);
    m.max_grid_size = {8, 8};
    ablsteady::TerrainInputs t;
    m.hierarchy_levels = 4;
    EXPECT_TRUE(msgs(m, t).empty());
    m.hierarchy_levels = 1;
    EXPECT_TRUE(msgs(m, t).empty());

    struct Case { amrex::IntVect n; int mgs; int L; int p; bool given; const char* expect; };
    const Case cases[] = {
        {{4, 4, 8}, 8, 3, 1, false, "mesh.hierarchy_levels = 3 leaves the coarsest level 1 x 1 cells (mesh.n_cell 4 4)"},
        {{30, 16, 8}, 8, 3, 1, false, "mesh.hierarchy_levels = 3 needs mesh.n_cell x and y divisible by 2^2 = 4 (got 30 16)"},
        {{32, 18, 8}, 8, 3, 1, false, "mesh.hierarchy_levels = 3 needs mesh.n_cell x and y divisible by 2^2 = 4 (got 32 18)"},
        {{32, 16, 8}, 64, 6, 1, false, "mesh.hierarchy_levels = 6 needs mesh.n_cell x and y divisible by 2^5 = 32 (got 32 16)"},
        {{32, 16, 8}, 2, 3, 1, false, "mesh.max_grid_size must be a multiple of 2^(mesh.hierarchy_levels - 1) = 4 in x and y (got 2 2)"},
        {{32, 16, 8}, 6, 3, 1, false, "mesh.max_grid_size must be a multiple of 2^(mesh.hierarchy_levels - 1) = 4 in x and y (got 6 6)"},
        {{32, 16, 8}, 8, 0, 1, false, "mesh.hierarchy_levels must be an integer from 1 to 10 (got 0)"},
        {{32, 16, 8}, 8, 11, 1, false, "mesh.hierarchy_levels must be an integer from 1 to 10 (got 11)"},
        {{32, 16, 8}, 8, 3, -1, false, "terrain.hierarchy_passes must be an integer from 0 to 64 (got -1)"},
        {{32, 16, 8}, 8, 3, 65, false, "terrain.hierarchy_passes must be an integer from 0 to 64 (got 65)"},
        {{32, 16, 8}, 8, 1, 2, true, "terrain.hierarchy_passes has no effect unless mesh.hierarchy_levels > 1"},
    };
    for (const Case& k : cases) {
        ablsteady::MeshInputs mm = m;
        mm.n_cell = k.n;
        mm.max_grid_size = {k.mgs, k.mgs};
        mm.hierarchy_levels = k.L;
        ablsteady::TerrainInputs tt;
        tt.hierarchy_passes = k.p;
        if (k.given) { tt.given.emplace_back("hierarchy_passes"); }
        const auto e = msgs(mm, tt);
        ASSERT_EQ(e.size(), 1u) << k.expect;
        EXPECT_EQ(e[0].rfind(k.expect, 0), 0u) << e[0];
        // validate() reports it too.
        ablsteady::Inputs in;
        in.mesh = mm;
        in.mesh.prob_hi = {1000, 1000, 1000};
        in.terrain = tt;
        const auto all = ablsteady::validate(in);
        bool found = false;
        for (const auto& s : all) { found = found || s.rfind(k.expect, 0) == 0; }
        EXPECT_TRUE(found) << k.expect;
    }

    // hierarchy_box_array: boxes of at most max_grid_size that every level
    // coarsens exactly, covering the domain, each spanning the full column.
    m.n_cell = amrex::IntVect(40, 24, 6);
    m.max_grid_size = {12, 8};
    m.hierarchy_levels = 3;
    const amrex::BoxArray ba = ablsteady::hierarchy_box_array(m);
    EXPECT_TRUE(ba.coarsenable(amrex::IntVect(4, 4, 1)));
    EXPECT_EQ(ba.numPts(), 40L * 24 * 6);
    EXPECT_TRUE(ba.minimalBox() == amrex::Box(amrex::IntVect(0), amrex::IntVect(39, 23, 5)));
    for (int b = 0; b < ba.size(); ++b) {
        EXPECT_LE(ba[b].length(0), 12);
        EXPECT_LE(ba[b].length(1), 8);
        EXPECT_EQ(ba[b].length(2), 6);
    }

    // A level 0 laid out without hierarchy_box_array (boxes of 6, not a
    // multiple of 4) is refused, naming mesh.max_grid_size.
    HCase c = hill3d(1);
    c.in.mesh.n_cell = amrex::IntVect(24, 16, 12);
    c.layout();
    ablsteady::MeshInputs six = c.in.mesh;
    six.max_grid_size = {6, 8};
    c.ba = ablsteady::make_box_array(six);
    c.dm = amrex::DistributionMapping(c.ba);
    c.z.define(amrex::convert(c.ba, amrex::IntVect(1)), c.dm, 1, ablsteady::ngrow_nodal);
    c.z.setVal(0.0);
    c.mesh.define(c.ba, c.dm);
    c.mesh.build(c.in, c.geom, c.z);
    try {
        ablsteady::MeshHierarchy h(c.in, c.geom, c.z, c.mesh);
        ADD_FAILURE() << "a layout that does not coarsen was accepted";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("mesh.max_grid_size: the level-0 boxes do not coarsen"), std::string::npos)
            << e.what();
    }
    // Inputs that do not divide abort in the constructor too.
    HCase d = hill3d(1);
    d.make();
    ablsteady::Inputs bad = d.in;
    bad.mesh.hierarchy_levels = 6;
    EXPECT_THROW(ablsteady::MeshHierarchy(bad, d.geom, d.z, d.mesh), std::runtime_error);
}

// ---------------------------------------------------------------------------
// UT-2.3.6. A coarse level gathered onto one box with ParallelCopy equals the
// distributed level bit for bit: node heights (ghost nodes included), detJ,
// the face areas, z_cc and the quality report; on 1, 2 and 4 ranks.
TEST (MeshHierarchy, GatheredLevelEqualsDistributed)
{
    for (const int sm : {0, 1}) {
        HCase a = hill3d(sm);
        a.make();
        HCase b = hill3d(sm);
        b.make(1);
        for (int l = 1; l < a.h->num_levels(); ++l) {
            ASSERT_TRUE(b.h->gathered(l));
            ASSERT_FALSE(a.h->gathered(l));
            EXPECT_EQ(b.h->boxArray(l).size(), 1);
            EXPECT_GT(a.h->boxArray(l).size(), 1);
            // The gathered fields copied back onto the distributed layout.
            auto back = [&] (const amrex::MultiFab& src, const amrex::MultiFab& like, int ng) {
                amrex::MultiFab out(like.boxArray(), like.DistributionMap(), 1, ng);
                out.setVal(0.0);
                // No periodic shift (as the gather): every node at its own index.
                out.ParallelCopy(src, 0, 0, 1, amrex::IntVect(ng), amrex::IntVect(ng));
                return out;
            };
            const amrex::MultiFab& za = a.h->z_phys_nd(l);
            EXPECT_EQ(ablsteady::test::count_bit_differences(back(b.h->z_phys_nd(l), za, 1), za, 1), 0)
                << "z_phys_nd, smoothing " << sm << ", level " << l;
            const ablsteady::Mesh& ma = a.h->mesh(l);
            const ablsteady::Mesh& mb = b.h->mesh(l);
            const std::pair<const amrex::MultiFab*, const amrex::MultiFab*> fields[] = {
                {&ma.detJ, &mb.detJ}, {&ma.ax, &mb.ax}, {&ma.ay, &mb.ay}, {&ma.az, &mb.az}, {&ma.z_cc, &mb.z_cc}};
            const char* names[] = {"detJ", "ax", "ay", "az", "z_cc"};
            for (int n = 0; n < 5; ++n) {
                EXPECT_EQ(ablsteady::test::count_bit_differences(back(*fields[n].second, *fields[n].first, 0),
                                                                 *fields[n].first, 0), 0)
                    << names[n] << ", smoothing " << sm << ", level " << l;
            }
            const ablsteady::MeshQuality& qa = ma.quality();
            const ablsteady::MeshQuality& qb = mb.quality();
            EXPECT_EQ(qa.min_thickness_ratio, qb.min_thickness_ratio);
            EXPECT_EQ(qa.max_aspect_ratio, qb.max_aspect_ratio);
            EXPECT_EQ(qa.max_surface_slope_deg, qb.max_surface_slope_deg);
            EXPECT_EQ(qa.stf_capped_levels, qb.stf_capped_levels);
            EXPECT_EQ(ma.quality_report(), mb.quality_report());
        }
    }
}
