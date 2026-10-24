// Tests of the Newton–Krylov core added after #19. The
// operator on terrain under the wall law (where one ghost fill is not a fixed
// point) against a central difference of q -> R(fill(q)), and Newton steps
// there; the rule that skips an update whose Krylov residual is above |b| or
// whose |b| is not finite, alone and inside a step; the memory budget per rank
// on several boxes; the inner product with every box on rank 0; the step's
// diagnostics; and a non-uniform dtau on the faces. test_Newton.cpp holds the
// core's own tests.
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Inputs.H"
#include "Newton.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using gust::FaceType;
using gust::Loc;
using gust::NewtonKrylov;
using gust::NewtonOperator;
using gust::NewtonOptions;
using gust::NewtonVector;
using gust::ptest::PCase;
using gust::trtest::tol;

constexpr double kPi = 3.14159265358979323846;
constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr Loc kLoc[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};

amrex::MultiFab& field (gust::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

const amrex::MultiFab& rfield (const gust::ResidualFields& R, int f)
{
    const amrex::MultiFab* p[6] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

// A smooth state with positive u, v and a perturbation of every field (the
// fields of test_Newton.cpp's QuadraticCase).
double base_value (int f, double x, double y, double z, double Lx, double Ly, double H)
{
    const double ax = 2 * kPi * x / Lx, ay = 2 * kPi * y / Ly, az = kPi * z / H;
    switch (f) {
    case 0: return 8 + 2 * std::sin(ax + 0.3) + z / H;
    case 1: return 3 + std::cos(ay) + 0.5 * std::sin(az);
    case 2: return 0.5 * std::sin(az) * (1 + 0.2 * std::cos(ax));
    case 3: return 5 * std::sin(ax) * std::cos(ay);
    case 4: return 300 + 0.01 * z + 0.5 * std::sin(ax);
    default: return 0.4 + 0.1 * std::cos(ay);
    }
}

double pert_value (int f, double x, double y, double z, double Lx, double Ly, double H)
{
    const double ax = 2 * kPi * x / Lx, ay = 2 * kPi * y / Ly, az = kPi * z / H;
    switch (f) {
    case 0: return std::cos(ay + z / H);
    case 1: return 0.5 * std::sin(ax);
    case 2: return 0.1 * std::sin(az) * std::sin(ay);
    case 3: return 2 * std::cos(ay);
    case 4: return 0.3 * std::cos(ax) * std::sin(az);
    default: return 0.05 * std::sin(ax);
    }
}

// A 9 x 6 x 8 mesh on 900 x 600 x 400 m over a Gaussian hill of height
// `hill` (30 m: slopes up to 0.12), periodic y, x periodic or inflow /
// outflow, a slip top and the ground `ground`; fixed diffusivities of
// 2 m^2/s, upwind convection, no buoyancy, z0 = 0.1 m.
struct HillCase
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::SurfaceLayerOptions sl;
    gust::ResidualInputs in;
    gust::State q, v;
    static constexpr double Lx = 900, Ly = 600, H = 400;
    HillCase (int max_grid, FaceType ground, bool perx, double hill)
    {
        c.in.mesh.n_cell = amrex::IntVect(9, 6, 8);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {perx ? 1 : 0, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([=] (double x, double) { return hill * std::exp(-(x - 450) * (x - 450) / (150.0 * 150.0)); });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        gust::trtest::smooth_vertical_ghost_nodes(c);
        sl.z0 = amrex::Real(0.1);
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(2.0);
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = 0;
        in.momentum.limiter = gust::Limiter::upwind;
        in.theta.limiter = gust::Limiter::upwind;
        in.surface = (ground == FaceType::surface_layer) ? &sl : nullptr;
        const FaceType xlo = perx ? FaceType::periodic : FaceType::inflow;
        const FaceType xhi = perx ? FaceType::periodic : FaceType::outflow;
        bin.face = {xlo, FaceType::periodic, ground, xhi, FaceType::periodic, FaceType::slip_wall};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        for (gust::State* s : {&q, &v}) {
            s->define(c.ba, c.dm);
            amrex::MultiFab::Copy(s->z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        }
        for (int f = 0; f < 6; ++f) {
            gust::trtest::fill_loc(field(q, f), kLoc[f], c,
                                        [&] (double x, double y, double z) { return base_value(f, x, y, z, Lx, Ly, H); });
            gust::trtest::fill_loc(field(v, f), kLoc[f], c,
                                        [&] (double x, double y, double z) { return pert_value(f, x, y, z, Lx, Ly, H); });
        }
        bc->fill_ghosts(q, in.surface, in.buoyancy);
        bc->fill_ghosts(v, in.surface, in.buoyancy, true);
    }
};

// to = fill(a + s b) (every field, ghosts included).
void combine_fill (gust::State& to, gust::State& a, gust::State& b, double s, const HillCase& hc)
{
    to.define(hc.c.ba, hc.c.dm);
    amrex::MultiFab::Copy(to.z_phys_nd, a.z_phys_nd, 0, 0, 1, a.z_phys_nd.nGrow());
    for (int f = 0; f < 6; ++f) {
        auto& t = field(to, f);
        amrex::MultiFab::LinComb(t, 1, field(a, f), 0, static_cast<amrex::Real>(s), field(b, f), 0, 0, 1, t.nGrow());
    }
    hc.bc->fill_ghosts(to, hc.in.surface, hc.in.buoyancy);
}

// x = S^-1 v on the valid points.
void to_vector (gust::State& v, const NewtonOperator& op, NewtonVector& x)
{
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab::Copy(x.f[f], field(v, f), 0, 0, 1, 0);
        x.f[f].mult(amrex::Real(1) / op.unknown_scales()[f], 0, 1, 0);
    }
}

double rel_diff (const NewtonOperator& op, const NewtonVector& a, const NewtonVector& b)
{
    NewtonVector d = op.makeVecRHS();
    op.linComb(d, 1, a, -1, b);
    return static_cast<double>(op.norm2(d)) / static_cast<double>(op.norm2(b));
}

// The exact J v of q -> R(fill(q)) to O(h^2): T^-1 (R(fill(q + h v)) -
// R(fill(q - h v))) / 2h, plus the identity of the Dirichlet rows.
void central_jv (HillCase& hc, const NewtonOperator& op, const NewtonVector& x, NewtonVector& E)
{
    const double h = 1e-4;
    gust::State qp, qm;
    combine_fill(qp, hc.q, hc.v, h, hc);
    combine_fill(qm, hc.q, hc.v, -h, hc);
    gust::ResidualFields Rp, Rm;
    gust::residual(qp, hc.in, hc.c.mesh, hc.c.geom, Rp);
    gust::residual(qm, hc.in, hc.c.mesh, hc.c.geom, Rm);
    for (int f = 0; f < 6; ++f) {
        const auto s = static_cast<amrex::Real>(0.5 / h) / op.row_scales()[f];
        amrex::MultiFab::LinComb(E.f[f], s, rfield(Rp, f), 0, -s, rfield(Rm, f), 0, 0, 1, 0);
    }
    const gust::Equation eqs[3] = {gust::Equation::u, gust::Equation::v, gust::Equation::w};
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab dir;
        hc.in.bc->dirichlet_rows(eqs[d], kLoc[d], hc.c.mesh, dir);
        amrex::MultiFab::AddProduct(E.f[d], dir, 0, x.f[d], 0, 0, 1, 0);
    }
}

// R(fill(q)) on a copy of q: the base point NewtonOperator needs.
void residual_of_refill (const HillCase& hc, gust::State& q, gust::ResidualFields& R)
{
    gust::State qf;
    gust::State zero;
    zero.define(hc.c.ba, hc.c.dm);
    combine_fill(qf, q, zero, 0, hc);
    gust::residual(qf, hc.in, hc.c.mesh, hc.c.geom, R);
    amrex::Gpu::streamSynchronize();
}

} // namespace

// ---------------------------------------------------------------------------
// On a 30 m hill with bc.zlo = surface_layer, NewtonOperator::apply matches
// the central difference of q -> R(fill(q)) to 1e-5 (2.5e-8 to 4.6e-8
// measured) on one box and on six uneven boxes, with x periodic and with
// inflow / outflow. The ghost fill is a fixed point
// (BoundaryConditions.TerrainFillIsAFixedPoint), so the base point R(q) of a
// filled q is R(fill(q)) and gives the same operator bit for bit. Before the
// fill was a fixed point the base point R(q) was off by more than 1 (322 and
// 568 measured). Flat ground under the wall law and the hill with a slip
// ground too. The wall law's ghost gradient depends on q, so the operator
// must also refill q + eps v (a refill = false mutant fails here).
TEST(Newton, TerrainWallLawOperatorMatchesCentralDifference)
{
    struct Config { FaceType ground; bool perx; double hill; };
    const Config configs[] = {
        {FaceType::surface_layer, true, 30},
        {FaceType::surface_layer, false, 30},
        {FaceType::surface_layer, false, 0},
        {FaceType::slip_wall, false, 30},
    };
    amrex::Print() << "Terrain J v: ground, periodic x, hill, boxes, error with R(fill(q)), with R(q)\n";
    for (const auto& cf : configs) {
        for (const int max_grid : {64, 4}) {
            HillCase hc(max_grid, cf.ground, cf.perx, cf.hill);
            gust::ResidualFields Rf, Rq;
            residual_of_refill(hc, hc.q, Rf);
            gust::residual(hc.q, hc.in, hc.c.mesh, hc.c.geom, Rq);
            NewtonOperator op(hc.q, Rf, hc.in, hc.c.mesh, hc.c.geom, nullptr);
            NewtonOperator opq(hc.q, Rq, hc.in, hc.c.mesh, hc.c.geom, nullptr);
            NewtonVector x = op.makeVecLHS(), Ax = op.makeVecRHS(), Aq = op.makeVecRHS(), E = op.makeVecRHS();
            to_vector(hc.v, op, x);
            op.apply(Ax, x);
            opq.apply(Aq, x);
            central_jv(hc, op, x, E);
            const double e = rel_diff(op, Ax, E), eq = rel_diff(op, Aq, E);
            amrex::Print() << std::setprecision(3) << "  " << static_cast<int>(cf.ground) << " " << cf.perx << " "
                           << cf.hill << " " << hc.c.ba.size() << "  " << e << "  " << eq << "\n";
            EXPECT_LT(e, tol(1e-5, 1e-1)) << static_cast<int>(cf.ground) << " " << cf.perx << " " << cf.hill << " " << max_grid;
            EXPECT_EQ(eq, e) << static_cast<int>(cf.ground) << " " << cf.perx << " " << cf.hill << " " << max_grid;
        }
    }
}

// ---------------------------------------------------------------------------
// The ghost fill is a fixed point: fill(fill(q)) = fill(q) on every field and
// every ghost cell, bit for bit, and so R(fill(fill(q))) = R(fill(q)). On the
// 30 m hill under the wall law, periodic and inflow/outflow in x, on one box and
// on six (box seams in x and y). Before the wall model took the ground cell's
// normal velocity from the no-flux rule on the cell's own u and v, it read the
// ground w face, which the fill then rewrites, and a second fill moved the
// ghosts by up to 1e-4 and R by about 1e-5. HillCase's q is filled once.
TEST(BoundaryConditions, TerrainFillIsAFixedPoint)
{
    for (const bool perx : {true, false}) {
        for (const int max_grid : {64, 4}) {
            HillCase hc(max_grid, FaceType::surface_layer, perx, 30);
            gust::State q2;
            gust::State zero;
            zero.define(hc.c.ba, hc.c.dm);
            combine_fill(q2, hc.q, zero, 0, hc);
            for (int f = 0; f < 6; ++f) {
                const amrex::MultiFab& a = field(hc.q, f);
                const amrex::MultiFab& b = field(q2, f);
                amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, a.nGrowVect());
                amrex::MultiFab::Copy(d, b, 0, 0, 1, d.nGrowVect());
                amrex::MultiFab::Subtract(d, a, 0, 0, 1, d.nGrowVect());
                EXPECT_EQ(d.norminf(0, d.nGrow()), 0) << "field " << f << " periodic x " << perx << " boxes "
                                                      << hc.c.ba.size();
            }
            gust::ResidualFields R1, R2;
            gust::residual(hc.q, hc.in, hc.c.mesh, hc.c.geom, R1);
            gust::residual(q2, hc.in, hc.c.mesh, hc.c.geom, R2);
            for (int f = 0; f < 6; ++f) {
                amrex::MultiFab d(rfield(R1, f).boxArray(), rfield(R1, f).DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(d, rfield(R2, f), 0, 0, 1, 0);
                amrex::MultiFab::Subtract(d, rfield(R1, f), 0, 0, 1, 0);
                EXPECT_EQ(d.norminf(0, 0), 0) << "residual " << f << " periodic x " << perx << " boxes " << hc.c.ba.size();
            }
            amrex::Gpu::streamSynchronize();
        }
    }
}

// ---------------------------------------------------------------------------
// Newton steps on the hill under the wall law (inflow / outflow, dtau = 20 s,
// the default restart 30 and rtol 1e-2): every GMRES solve meets its
// tolerance, every update is applied and |R| decreases at every step
// (11.3, 1.13, 0.537, 0.303, 0.211 measured). Before the base point was
// R(fill(q)), the first solve stopped at the 2000-iteration cap with a
// Krylov residual 1.9e44 times |b|. One box and the uneven boxes agree to
// 1e-6 after the steps (4.0e-8 measured: the two layouts add their sums in
// different orders, and four inexact solves carry that round-off).
TEST(Newton, TerrainWallLawStepsReduceResidual)
{
    std::vector<std::vector<double>> hist;
    NewtonVector qs[2];
    HillCase ref(64, FaceType::surface_layer, false, 30);
    for (int L = 0; L < 2; ++L) {
        HillCase hc(L == 0 ? 64 : 4, FaceType::surface_layer, false, 30);
        amrex::MultiFab dtau(hc.c.ba, hc.c.dm, 1, 0);
        dtau.setVal(20.0);
        NewtonKrylov nk(NewtonOptions{}, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        std::vector<double> r;
        for (int n = 0; n < 4; ++n) {
            const auto info = nk.step(hc.q, &dtau);
            r.push_back(static_cast<double>(info.residual_norm));
            EXPECT_EQ(info.gmres_status, 0) << L << " " << n;
            EXPECT_TRUE(info.applied) << L << " " << n;
            EXPECT_LE(static_cast<double>(info.krylov_relative), 1e-2 * (1 + 1e-6)) << L << " " << n;
        }
        gust::ResidualFields R;
        residual_of_refill(hc, hc.q, R);
        NewtonOperator op(hc.q, R, nk.residual_inputs(), hc.c.mesh, hc.c.geom, nullptr);
        NewtonVector b = op.makeVecRHS();
        op.rhs(R, b);
        r.push_back(static_cast<double>(op.norm2(b)));
        hist.push_back(r);
        qs[L].define(ref.c.ba, ref.c.dm);
        for (int f = 0; f < 6; ++f) {
            qs[L].f[f].setVal(0.0);
            qs[L].f[f].ParallelCopy(field(hc.q, f), 0, 0, 1, 0, 0, hc.c.geom.periodicity());
        }
    }
    amrex::Print() << "Terrain Newton |R|:";
    for (const double x : hist[1]) { amrex::Print() << " " << std::setprecision(3) << x; }
    amrex::Print() << "\n";
    for (const auto& r : hist) {
        for (std::size_t n = 0; n + 1 < r.size(); ++n) { EXPECT_LT(r[n + 1], r[n]) << n; }
        EXPECT_LT(r.back(), 0.05 * r.front());
    }
    gust::ResidualFields R;
    NewtonOperator opr(ref.q, R, ref.in, ref.c.mesh, ref.c.geom, nullptr);
    EXPECT_LT(rel_diff(opr, qs[1], qs[0]), tol(1e-6, 1e-1));
}

// ---------------------------------------------------------------------------
// The rule that skips an update: applied when the Krylov residual is at most
// |b| (the edge 1 included), skipped above it and for NaN and inf, and NaN
// raises no invalid flag (StepSkipsUnusableUpdate runs it inside a step).
TEST(Newton, UpdateSkippedAboveRhsNorm)
{
    using gust::newton_update_usable;
    const amrex::Real one = 1;
    EXPECT_TRUE(newton_update_usable(0));
    EXPECT_TRUE(newton_update_usable(amrex::Real(1e-2)));
    EXPECT_TRUE(newton_update_usable(one));
    EXPECT_FALSE(newton_update_usable(std::nextafter(one, amrex::Real(2))));
    EXPECT_FALSE(newton_update_usable(amrex::Real(1.9e30)));
    EXPECT_FALSE(newton_update_usable(std::numeric_limits<amrex::Real>::infinity()));
    EXPECT_FALSE(newton_update_usable(std::numeric_limits<amrex::Real>::quiet_NaN()));
    // NaN must not raise the invalid flag: an ordered comparison does with
    // gcc (x86 and arm64; AppleClang uses a quiet compare), and the CI's
    // amrex.fpe_trap_invalid=1 turns that into an abort. The flag is tested
    // directly, so a gcc build catches it where traps are off (macOS).
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    std::feclearexcept(FE_INVALID);
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    const volatile bool usable = newton_update_usable(nan);
    const bool raised = std::fetestexcept(FE_INVALID) != 0;
    std::feclearexcept(FE_INVALID);
    amrex::setFPExcept(prev_fpe);
    EXPECT_FALSE(usable);
    EXPECT_FALSE(raised);
}

// ---------------------------------------------------------------------------
// NewtonKrylov::step skips the update: (a) a hook that returns NaN makes the
// Krylov residual NaN, and (b) a NaN in the state makes |R| NaN; in both the
// step reports applied = false, prints a line saying so (at verbose 0) and
// the valid points of q are unchanged (a
// flat case with a slip ground, whose ghost fill is idempotent: exactly the
// same bits). The invalid trap is off around both: NaN is the input here.
namespace {
class NaNHook : public gust::NewtonPreconditioner
{
public:
    void apply (NewtonVector& dq, const NewtonVector& /*r*/) override
    {
        for (auto& f : dq.f) { f.setVal(std::numeric_limits<amrex::Real>::quiet_NaN()); }
    }
};

// max |a - b| over the six fields' valid points; NaN where both are NaN counts as equal.
double valid_diff (gust::State& a, gust::State& b)
{
    double d = 0;
    for (int f = 0; f < 6; ++f) {
        for (amrex::MFIter mfi(field(a, f), false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            amrex::FArrayBox fa(vb, 1, amrex::The_Pinned_Arena()), fb(vb, 1, amrex::The_Pinned_Arena());
            fa.copy<amrex::RunOn::Device>(field(a, f)[mfi], vb);
            fb.copy<amrex::RunOn::Device>(field(b, f)[mfi], vb);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                const double x = fa(amrex::IntVect(i, j, k)), y = fb(amrex::IntVect(i, j, k));
                if (std::isnan(x) && std::isnan(y)) { return; }
                d = std::max(d, std::isnan(x - y) ? std::numeric_limits<double>::infinity() : std::abs(x - y));
            });
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(d);
    return d;
}

void copy_state (gust::State& to, gust::State& from, const HillCase& hc)
{
    gust::State zero;
    zero.define(hc.c.ba, hc.c.dm);
    combine_fill(to, from, zero, 0, hc);
}
} // namespace

TEST(Newton, StepSkipsUnusableUpdate)
{
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    NewtonOptions opt;
    opt.gmres_restart = 5;
    {
        HillCase hc(4, FaceType::slip_wall, true, 0);
        gust::State before;
        copy_state(before, hc.q, hc);
        NewtonKrylov nk(opt, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        NaNHook nan;
        testing::internal::CaptureStdout();
        const auto info = nk.step(hc.q, nullptr, &nan);
        const std::string out = testing::internal::GetCapturedStdout();
        if (amrex::ParallelDescriptor::IOProcessor()) {
            EXPECT_NE(out.find("the update is not applied"), std::string::npos) << out;   // printed at verbose 0 too
        }
        EXPECT_TRUE(std::isfinite(static_cast<double>(info.residual_norm)));
        EXPECT_FALSE(std::isfinite(static_cast<double>(info.krylov_relative)));
        EXPECT_FALSE(info.applied);
        EXPECT_EQ(valid_diff(hc.q, before), 0);
    }
    {
        HillCase hc(4, FaceType::slip_wall, true, 0);
        // One NaN theta point in box 0 (on whichever rank owns it; |R| is
        // all-reduced, so every rank sees NaN).
        for (amrex::MFIter mfi(hc.q.theta, false); mfi.isValid(); ++mfi) {
            if (mfi.index() != 0) { continue; }
            const amrex::Box vb = mfi.validbox();
            const amrex::Box one(vb.smallEnd(), vb.smallEnd());
            hc.q.theta[mfi].setVal<amrex::RunOn::Device>(std::numeric_limits<amrex::Real>::quiet_NaN(), one);
        }
        amrex::Gpu::streamSynchronize();
        gust::State before;
        copy_state(before, hc.q, hc);
        NewtonKrylov nk(opt, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        const auto info = nk.step(hc.q, nullptr);
        EXPECT_FALSE(std::isfinite(static_cast<double>(info.residual_norm)));
        EXPECT_FALSE(info.applied);
        EXPECT_EQ(valid_diff(hc.q, before), 0);
    }
    amrex::setFPExcept(prev_fpe);
}

// ---------------------------------------------------------------------------
// The memory budget applies to the rank with the most points: one
// 512 x 512 x 128 box per rank with restart 20 (43 vectors, 69.4 GB per box
// in double) fits the 70 GB default on every rank count, while the total
// over 2 or more ranks exceeds it. A box layout whose boxes all sit on rank 0
// gives that rank the total. On one rank both are the same layout.
TEST(Newton, KrylovBudgetIsPerRank)
{
    const int np = amrex::ParallelDescriptor::NProcs();
    amrex::BoxList bl;
    for (int n = 0; n < np; ++n) {
        bl.push_back(amrex::Box(amrex::IntVect(512 * n, 0, 0), amrex::IntVect(512 * n + 511, 511, 127)));
    }
    const amrex::BoxArray ba(bl);
    amrex::Vector<int> spread(np), packed(np, 0);
    for (int n = 0; n < np; ++n) { spread[n] = n; }
    const amrex::DistributionMapping dms(spread), dmp(packed);
    const double per = (3 * 512.0 * 512 * 128 + 513.0 * 512 * 128 + 512.0 * 513 * 128 + 512.0 * 512 * 129)
                       * sizeof(amrex::Real);
    const auto ms = gust::krylov_memory(20, ba, dms);
    const auto mp = gust::krylov_memory(20, ba, dmp);
    EXPECT_EQ(ms.vectors, 43);
    EXPECT_EQ(ms.total_bytes, 43 * per * np);
    EXPECT_EQ(ms.max_rank_bytes, 43 * per);
    EXPECT_EQ(mp.max_rank_bytes, 43 * per * np);
    EXPECT_EQ(mp.max_rank, 0);

    if constexpr (kDouble) {
        // Any small case gives the residual inputs; only the layout is checked.
        HillCase hc(64, FaceType::slip_wall, true, 0);
        NewtonOptions opt;
        opt.gmres_restart = 20;
        EXPECT_NO_THROW(NewtonKrylov(opt, hc.in, hc.c.mesh, hc.c.geom, ba, dms, 0));
        if (np > 1) { EXPECT_THROW(NewtonKrylov(opt, hc.in, hc.c.mesh, hc.c.geom, ba, dmp, 0), std::runtime_error); }
    }
}

// ---------------------------------------------------------------------------
// dotProduct with every box on rank 0 against the distributed layout of the
// same 15 boxes: the same bits. Moving every box one rank on would not test
// this on 2 ranks (each rank's group moves whole, and an all-reduce of two
// partial sums is the same in either order); putting every box on one rank
// changes the grouping. Rough values make the regrouped partial sums round
// apart (a rank-local sum fails here on 2 and 4 ranks). On one rank the two
// layouts are the same.
TEST(Newton, InnerProductSameBitsWithAllBoxesOnRankZero)
{
    HillCase a(2, FaceType::slip_wall, true, 30);
    // Rough values (a hash of the position, as rough as a Krylov vector), so
    // that partial sums grouped differently round apart.
    for (gust::State* s : {&a.q, &a.v}) {
        const double seed = (s == &a.q) ? 0.0 : 1.7;
        for (int f = 0; f < 6; ++f) {
            gust::trtest::fill_loc(field(*s, f), kLoc[f], a.c, [&] (double x, double y, double z) {
                const double h = std::sin(12.9898 * x + 78.233 * y + 37.719 * z + 4.1 * f + seed) * 43758.5453;
                return h - std::floor(h) - 0.5;
            });
        }
    }
    gust::ResidualFields R;
    NewtonOperator op(a.q, R, a.in, a.c.mesh, a.c.geom, nullptr);
    NewtonVector x = op.makeVecLHS(), y = op.makeVecLHS();
    to_vector(a.q, op, x);
    to_vector(a.v, op, y);
    const amrex::Real dot = op.dotProduct(x, y);
    ASSERT_GE(a.c.ba.size(), 15);

    const amrex::DistributionMapping dm0(amrex::Vector<int>(a.c.ba.size(), 0));
    HillCase b(2, FaceType::slip_wall, true, 30);
    PCase& c0 = b.c;
    ASSERT_EQ(c0.ba, a.c.ba);
    c0.dm = dm0;
    c0.z.define(a.c.z.boxArray(), dm0, 1, a.c.z.nGrow());
    c0.z.ParallelCopy(a.c.z, 0, 0, 1, a.c.z.nGrow(), a.c.z.nGrow());
    c0.mesh.define(c0.ba, dm0);
    c0.mesh.build_from_ground(c0.in, c0.geom, c0.z);
    gust::State q0;
    q0.define(c0.ba, dm0);
    amrex::MultiFab::Copy(q0.z_phys_nd, c0.z, 0, 0, 1, c0.z.nGrow());
    NewtonOperator op0(q0, R, a.in, c0.mesh, c0.geom, nullptr);
    NewtonVector x0 = op0.makeVecLHS(), y0 = op0.makeVecLHS();
    for (int f = 0; f < 6; ++f) {
        x0.f[f].ParallelCopy(x.f[f]);
        y0.f[f].ParallelCopy(y.f[f]);
    }
    EXPECT_EQ(op0.dotProduct(x0, y0), dot);
    EXPECT_NE(dot, 0);
}

// ---------------------------------------------------------------------------
// The step's diagnostics: krylov_relative = krylov_residual / residual_norm,
// above 0 and at most the tolerance once GMRES converges; eps_min and
// eps_max the smallest and largest jvp step of the operator. jvp's eps is
// h (1 + mean |q~|) / rms(v~), so x and 10 x give eps and eps / 10 whatever
// the order of the calls.
TEST(Newton, StepDiagnosticsAndEpsRange)
{
    HillCase hc(4, FaceType::surface_layer, false, 30);
    {
        gust::ResidualFields R;
        residual_of_refill(hc, hc.q, R);
        NewtonOperator op(hc.q, R, hc.in, hc.c.mesh, hc.c.geom, nullptr);
        NewtonVector x = op.makeVecLHS(), x10 = op.makeVecLHS(), Ax = op.makeVecRHS();
        to_vector(hc.v, op, x);
        op.linComb(x10, 10, x, 0, x);
        op.apply(Ax, x10);
        const amrex::Real e10 = op.eps_max();
        op.apply(Ax, x);
        EXPECT_EQ(op.jvp_calls(), 2);
        EXPECT_NEAR(static_cast<double>(op.eps_max() / op.eps_min()), 10, tol(1e-10, 1e-4));
        EXPECT_EQ(op.eps_min(), e10);
    }
    amrex::MultiFab dtau(hc.c.ba, hc.c.dm, 1, 0);
    dtau.setVal(20.0);
    NewtonOptions opt;
    NewtonKrylov nk(opt, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
    const auto info = nk.step(hc.q, &dtau);
    EXPECT_EQ(info.gmres_status, 0);
    EXPECT_GT(info.residual_norm, 0);
    EXPECT_GT(info.krylov_residual, 0);
    EXPECT_NEAR(static_cast<double>(info.krylov_relative),
                static_cast<double>(info.krylov_residual / info.residual_norm), tol(1e-14, 1e-6));
    EXPECT_LE(static_cast<double>(info.krylov_relative), static_cast<double>(opt.krylov_rtol) * (1 + 1e-6));
    EXPECT_GT(info.eps_min, 0);
    EXPECT_LE(info.eps_min, info.eps_max);
    EXPECT_GE(info.jvp_calls, info.gmres_iterations);   // one per iteration, plus the residual of each cycle
}

// ---------------------------------------------------------------------------
// A non-uniform dtau (10 + i + 2 j s in cell (i, j, k)): apply with it minus
// apply without is x / dtau_f on the momentum faces, with dtau_f the mean of
// the two cells (the periodic image across the x seam), x / dtau on theta
// and k, and zero on continuity and on the Dirichlet w faces.
TEST(Newton, NonUniformDtauOnFaces)
{
    HillCase hc(4, FaceType::slip_wall, true, 30);
    const int nx = 9, ny = 6;
    amrex::MultiFab dtau(hc.c.ba, hc.c.dm, 1, 0);
    gust::trtest::fill_loc(dtau, Loc::cell, hc.c, [&] (double x, double y, double) {
        const auto i = static_cast<int>(std::floor(x / (HillCase::Lx / nx)));
        const auto j = static_cast<int>(std::floor(y / (HillCase::Ly / ny)));
        return 10.0 + i + 2.0 * j;
    });
    // dtau on the faces, written from the cell indices (i - 1 of i = 0 is the
    // periodic image nx - 1 in x and ny - 1 in y).
    auto cell = [&] (int i, int j) { return 10.0 + ((i + nx) % nx) + 2.0 * ((j + ny) % ny); };
    const double dx = HillCase::Lx / nx, dy = HillCase::Ly / ny;
    gust::ResidualFields R;
    gust::residual(hc.q, hc.in, hc.c.mesh, hc.c.geom, R);
    NewtonOperator op0(hc.q, R, hc.in, hc.c.mesh, hc.c.geom, nullptr);
    NewtonOperator opd(hc.q, R, hc.in, hc.c.mesh, hc.c.geom, &dtau);
    NewtonVector x = op0.makeVecLHS(), A0 = op0.makeVecRHS(), Ad = op0.makeVecRHS(), D = op0.makeVecRHS(),
                 E = op0.makeVecRHS();
    to_vector(hc.v, op0, x);
    op0.apply(A0, x);
    opd.apply(Ad, x);
    op0.linComb(D, 1, Ad, -1, A0);
    // The expected x / dtau, point by point.
    gust::State inv;
    inv.define(hc.c.ba, hc.c.dm);
    gust::trtest::fill_loc(inv.u, Loc::xface, hc.c, [&] (double x0, double y0, double) {
        const auto i = static_cast<int>(std::lround(x0 / dx)), j = static_cast<int>(std::floor(y0 / dy));
        return 1.0 / (0.5 * (cell(i - 1, j) + cell(i, j)));
    });
    gust::trtest::fill_loc(inv.v, Loc::yface, hc.c, [&] (double x0, double y0, double) {
        const auto i = static_cast<int>(std::floor(x0 / dx)), j = static_cast<int>(std::lround(y0 / dy));
        return 1.0 / (0.5 * (cell(i, j - 1) + cell(i, j)));
    });
    gust::trtest::fill_loc(inv.w, Loc::zface, hc.c, [&] (double x0, double y0, double) {
        return 1.0 / cell(static_cast<int>(std::floor(x0 / dx)), static_cast<int>(std::floor(y0 / dy)));
    });
    for (amrex::MultiFab* mf : {&inv.theta, &inv.tke}) {
        gust::trtest::fill_loc(*mf, Loc::cell, hc.c, [&] (double x0, double y0, double) {
            return 1.0 / cell(static_cast<int>(std::floor(x0 / dx)), static_cast<int>(std::floor(y0 / dy)));
        });
    }
    // The Dirichlet w faces (ground and top) have no pseudo-time term.
    amrex::MultiFab dir;
    hc.in.bc->dirichlet_rows(gust::Equation::w, Loc::zface, hc.c.mesh, dir);
    for (int f = 0; f < 6; ++f) {
        E.f[f].setVal(0.0);
        if (f == gust::nf_p) { continue; }
        amrex::MultiFab::AddProduct(E.f[f], field(inv, f), 0, x.f[f], 0, 0, 1, 0);
    }
    amrex::MultiFab keep(dir.boxArray(), dir.DistributionMap(), 1, 0);
    keep.setVal(1.0);
    amrex::MultiFab::Subtract(keep, dir, 0, 0, 1, 0);
    amrex::MultiFab::Multiply(E.f[gust::nf_w], keep, 0, 0, 1, 0);
    EXPECT_LT(rel_diff(op0, D, E), tol(1e-12, 1e-5));
}

// ---------------------------------------------------------------------------
// The fixed point in more configurations (from the implementation review of
// audit A2): a hill varying in x and y, inflow / outflow in x and y (corners
// outside two non-periodic faces), the faces reversed, heat flux of both signs
// (the Newton solve for u*), bc.surface_zref above the first centre (the k > 0
// interval, and periodic with some columns in each), a geostrophic top; one box,
// 3 x 2 boxes split evenly and unevenly; the fill and the homogeneous fill of the
// perturbation. With the old ground-cell wind every configuration with k = 0
// fails (ghosts off by 0.002 to 0.055).
namespace {
struct Cfg { const char* name; bool perx, pery, heat; double qflux, zref; bool geo, reverse; double xc = 450; };

struct HillXY
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::SurfaceLayerOptions sl;
    gust::ResidualInputs in;
    gust::State q, v;
    static constexpr double Lx = 900, Ly = 600, H = 400;
    HillXY (int max_grid, const Cfg& cf)
    {
        c.in.mesh.n_cell = amrex::IntVect(9, 6, 8);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {cf.perx ? 1 : 0, cf.pery ? 1 : 0, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        const double xc = cf.xc;
        c.make([=] (double x, double y) { return 30 * std::exp(-((x - xc) * (x - xc) + 0.5 * (y - 300) * (y - 300)) / (150.0 * 150.0)); });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        gust::trtest::smooth_vertical_ghost_nodes(c);
        sl.z0 = amrex::Real(0.1);
        sl.zref = amrex::Real(cf.zref);
        if (cf.heat) { sl.stability = gust::SurfaceLayerOptions::Stability::heat_flux; sl.heat_flux = amrex::Real(cf.qflux); }
        for (auto* mf : {&nu, &gt, &gk}) { mf->define(c.ba, c.dm, 1, 1); mf->setVal(2.0); }
        in.nu = &nu; in.gamma_theta = &gt; in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = cf.heat ? amrex::Real(9.81) : amrex::Real(0);
        in.momentum.limiter = gust::Limiter::upwind;
        in.theta.limiter = gust::Limiter::upwind;
        in.surface = &sl;
        if (cf.geo) { in.coriolis.enabled = true; in.coriolis.geostrophic = true; in.coriolis.latitude = 45; in.coriolis.ug = 9; in.coriolis.vg = 2; }
        FaceType xlo = cf.perx ? FaceType::periodic : (cf.reverse ? FaceType::outflow : FaceType::inflow);
        FaceType xhi = cf.perx ? FaceType::periodic : (cf.reverse ? FaceType::inflow : FaceType::outflow);
        FaceType ylo = cf.pery ? FaceType::periodic : FaceType::inflow;
        FaceType yhi = cf.pery ? FaceType::periodic : FaceType::outflow;
        bin.face = {xlo, ylo, FaceType::surface_layer, xhi, yhi, cf.geo ? FaceType::geostrophic : FaceType::slip_wall};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        for (gust::State* s : {&q, &v}) {
            s->define(c.ba, c.dm);
            amrex::MultiFab::Copy(s->z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        }
        for (int f = 0; f < 6; ++f) {
            gust::trtest::fill_loc(field(q, f), kLoc[f], c, [&] (double x, double y, double z) { return base_value(f, x, y, z, Lx, Ly, H); });
            gust::trtest::fill_loc(field(v, f), kLoc[f], c, [&] (double x, double y, double z) { return pert_value(f, x, y, z, Lx, Ly, H); });
        }
        bc->fill_ghosts(q, in.surface, in.buoyancy);
        bc->fill_ghosts(v, in.surface, in.buoyancy, true);
    }
};

double max_ghost_diff (gust::State& a, gust::State& b, int f)
{
    const amrex::MultiFab& x = field(a, f);
    amrex::MultiFab d(x.boxArray(), x.DistributionMap(), 1, x.nGrowVect());
    amrex::MultiFab::Copy(d, field(b, f), 0, 0, 1, d.nGrowVect());
    amrex::MultiFab::Subtract(d, x, 0, 0, 1, d.nGrowVect());
    return static_cast<double>(d.norminf(0, d.nGrow()));
}

void copy_all (gust::State& to, gust::State& from)
{
    to.define(from.p.boxArray(), from.p.DistributionMap());
    amrex::MultiFab::Copy(to.z_phys_nd, from.z_phys_nd, 0, 0, 1, from.z_phys_nd.nGrow());
    for (int f = 0; f < 6; ++f) { amrex::MultiFab::Copy(field(to, f), field(from, f), 0, 0, 1, field(to, f).nGrowVect()); }
}
} // namespace

TEST(BoundaryConditions, TerrainFillIsAFixedPointInMoreConfigurations)
{
    const Cfg cfgs[] = {
        {"inout_x_inout_y", false, false, false, 0, 0, false, false},
        {"reverse_x", false, true, false, 0, 0, false, true},
        {"heat_unstable", false, true, true, 0.1, 0, false, false},
        {"heat_stable", false, true, true, -0.02, 0, false, false},
        {"zref_100", false, true, false, 0, 100, false, false},
        {"zref_100_heat_xy", false, false, true, 0.05, 100, false, false},
        {"geostrophic_top", true, true, false, 0, 0, true, false},
        {"all_periodic_zref60", true, true, false, 0, 60, false, false},
    };
    for (const auto& cf : cfgs) {
        for (const int mg : {64, 4, 3}) {
            HillXY hc(mg, cf);
            gust::State q2, v2;
            copy_all(q2, hc.q);
            hc.bc->fill_ghosts(q2, hc.in.surface, hc.in.buoyancy);
            copy_all(v2, hc.v);
            hc.bc->fill_ghosts(v2, hc.in.surface, hc.in.buoyancy, true);
            for (int f = 0; f < 6; ++f) {
                const double dq = max_ghost_diff(hc.q, q2, f), dv = max_ghost_diff(hc.v, v2, f);
                EXPECT_EQ(dq, 0) << cf.name << " " << mg << " " << f;
                EXPECT_EQ(dv, 0) << cf.name << " " << mg << " " << f;
            }
            amrex::Gpu::streamSynchronize();
        }
    }
}

// ---------------------------------------------------------------------------
// Every copy of a ghost cell holds the same value (audit G3). After the fill
// on 3 x 2 boxes, split evenly (3 cells) and unevenly (4 + 2 in y), each ghost
// cell equals the value one box over the whole domain gives the same cell:
// its own column inside the domain (the periodic image in a periodic
// direction) or, outside a non-periodic face, the one box's ghost there. The
// configurations of TerrainFillIsAFixedPointInMoreConfigurations; the fill and
// the homogeneous fill of the perturbation. Before, the outermost lateral
// ghost layer differed by up to 5.0 m/s for u and 0.24 for the perturbation's w.
namespace {
//! The largest |a - ref| over the cells of a's fabs the reference reached
//! (marked cells were not reached).
double max_diff_reached (const amrex::MultiFab& a, const amrex::MultiFab& ref, amrex::Real marker)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, a.nGrowVect());
    for (amrex::MFIter mfi(d, false); mfi.isValid(); ++mfi) {
        const auto x = a.const_array(mfi), r = ref.const_array(mfi);
        const auto e = d.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            e(i, j, k) = (r(i, j, k) == marker) ? amrex::Real(0) : std::abs(x(i, j, k) - r(i, j, k));
        });
    }
    return static_cast<double>(d.norminf(0, d.nGrow()));
}

double max_ghost_copy_diff (gust::State& a, gust::State& one, int f, const amrex::Geometry& geom)
{
    const amrex::MultiFab& x = field(a, f);
    const amrex::MultiFab& s = field(one, f);
    const amrex::Real marker = std::numeric_limits<amrex::Real>::max();
    amrex::MultiFab ref(x.boxArray(), x.DistributionMap(), 1, x.nGrowVect());
    ref.setVal(marker);
    const int ng = x.nGrow();
    const amrex::IntVect src_ng(geom.isPeriodic(0) ? 0 : ng, geom.isPeriodic(1) ? 0 : ng, ng);
    ref.ParallelCopy(s, 0, 0, 1, src_ng, x.nGrowVect(), geom.periodicity());
    return max_diff_reached(x, ref, marker);
}
} // namespace

TEST(BoundaryConditions, GhostCopiesAgreeAcrossBoxes)
{
    const Cfg cfgs[] = {
        {"inout_x_inout_y", false, false, false, 0, 0, false, false},
        {"reverse_x", false, true, false, 0, 0, false, true},
        {"heat_unstable", false, true, true, 0.1, 0, false, false},
        {"heat_stable", false, true, true, -0.02, 0, false, false},
        {"zref_100", false, true, false, 0, 100, false, false},
        {"zref_100_heat_xy", false, false, true, 0.05, 100, false, false},
        {"geostrophic_top", true, true, false, 0, 0, true, false},
        {"all_periodic_zref60", true, true, false, 0, 60, false, false},
    };
    for (const auto& cf : cfgs) {
        HillXY one(64, cf);
        for (const int mg : {4, 3}) {
            HillXY hc(mg, cf);
            for (int f = 0; f < 6; ++f) {
                const double dq = max_ghost_copy_diff(hc.q, one.q, f, hc.c.geom);
                const double dv = max_ghost_copy_diff(hc.v, one.v, f, hc.c.geom);
                EXPECT_EQ(dq, 0) << cf.name << " " << mg << " " << f;
                EXPECT_EQ(dv, 0) << cf.name << " " << mg << " " << f;
            }
            amrex::Gpu::streamSynchronize();
        }
    }
}

// ---------------------------------------------------------------------------
// The residual does not depend on the box layout where an inflow face crosses
// terrain (audit A3, from Pressel's review). The hill sits on the inflow face
// (and on the inflow / outflow y faces), nu_t comes from the closure with the
// realizability cap C_r k / S (turb.realizability = 0.309), and R on 3 x 2
// boxes, split evenly and unevenly, equals R on one box to 1e-12 of its norm
// per field (0 measured). The closure evaluates nu_t in the first lateral
// ghost layer from the strain there, which reads the second layer: the inflow
// face points of u. Before every ghost copy held one value those came from
// the clamped heights, and R differed by up to 1.6e-5 of its norm (w).
TEST(BoundaryConditions, ResidualWithClosureAtInflowSeamsIsLayoutFree)
{
    const Cfg cfgs[] = {
        {"hill_on_inflow_x", false, true, false, 0, 0, false, false, 0},
        {"hill_on_inflow_xy", false, false, false, 0, 0, false, false, 0},
        {"hill_on_inflow_reverse", false, true, false, 0, 0, false, true, 900},
    };
    gust::ClosureOptions opt;
    opt.nu_molecular = amrex::Real(1.5e-5);
    opt.realizability = amrex::Real(0.309);   // nu_t capped by C_r k / S: the strain enters nu_t
    const gust::Closure closure(opt);
    auto run = [&] (HillXY& h, gust::ResidualFields& R, amrex::MultiFab& dist) {
        h.in.momentum.limiter = gust::Limiter::van_leer;
        h.in.theta.limiter = gust::Limiter::van_leer;
        h.in.nu = nullptr;
        h.in.gamma_theta = nullptr;
        h.in.gamma_tke = nullptr;
        h.in.closure = &closure;
        gust::height_above_ground(h.q.z_phys_nd, dist);
        h.in.wall_distance = &dist;
        gust::residual(h.q, h.in, h.c.mesh, h.c.geom, R);
    };
    for (const auto& cf : cfgs) {
        HillXY one(64, cf);
        gust::ResidualFields R1;
        amrex::MultiFab d1;
        run(one, R1, d1);
        for (const int mg : {4, 3}) {
            HillXY hc(mg, cf);
            gust::ResidualFields R;
            amrex::MultiFab dd;
            run(hc, R, dd);
            for (int f = 0; f < 6; ++f) {
                const amrex::MultiFab& a = rfield(R1, f);
                amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
                d.setVal(0.0);
                d.ParallelCopy(rfield(R, f), 0, 0, 1, 0, 0);
                amrex::MultiFab::Subtract(d, a, 0, 0, 1, 0);
                const double n = static_cast<double>(a.norminf(0, 0)), e = static_cast<double>(d.norminf(0, 0));
                EXPECT_LE(e, tol(1e-12, 1e-4) * n) << cf.name << " " << mg << " " << f;
            }
            amrex::Gpu::streamSynchronize();
        }
    }
}
