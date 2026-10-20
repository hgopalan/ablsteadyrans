#include "Closure.H"

#include "CVFluxes.H"
#include "Stencil_K.H"

#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace ablsteady {

namespace {

const char* const kKeys[] = {"Cmu0", "Cb", "Rt_crit", "Rt_min", "max_geom_lscale", "tke_floor", "sigma_k",
                             "nu_molecular", "alpha_molecular", "kappa", "nut_relax", "frozen_switch", "face_mean",
                             "wall_distance"};

//! Encode and decode the switches of one point.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int encode (const Switches& s) noexcept
{
    return s.branch | (s.smooth_Rt_g ? 4 : 0) | (s.smooth_Rt ? 8 : 0) | (s.floored ? 16 : 0);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Switches decode (int code) noexcept
{
    Switches s;
    s.branch = code & 3;
    s.smooth_Rt_g = (code & 4) != 0;
    s.smooth_Rt = (code & 8) != 0;
    s.floored = (code & 16) != 0;
    return s;
}

//! The valid box grown by one cell in x and y only (the closure's cells).
amrex::Box lateral_grown (const amrex::Box& vb) { return amrex::grow(vb, amrex::IntVect(1, 1, 0)); }

//! The vertical ghost layers of a cell field take the adjacent cell's values
//! (zero gradient), over the whole lateral extent of each box.
void copy_vertical_ghosts (amrex::MultiFab& mf)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox(), fb = mfi.fabbox();
        const auto a = mf.array(mfi);
        const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
        amrex::Box slab = fb;
        slab.setRange(2, klo - 1, 1);
        amrex::ParallelFor(slab, mf.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            a(i, j, k, n) = a(i, j, klo, n);
        });
        slab.setRange(2, khi + 1, 1);
        amrex::ParallelFor(slab, mf.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            a(i, j, k, n) = a(i, j, khi, n);
        });
    }
}

//! The options after the checks (every message, each naming its key); the
//! constants are formed from checked options only (1/Cb^2 of a bad Cb is not
//! taken).
const ClosureOptions& checked (const ClosureOptions& opt)
{
    const auto err = check_closure_options(opt);
    if (!err.empty()) {
        std::string all = "ablsteady: invalid closure options:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

} // namespace

std::vector<std::string> check_closure_options (const ClosureOptions& opt)
{
    std::vector<std::string> err;
    // isnan first, and the ordered comparison only afterwards (inside the
    // predicate): gcc's ordered comparisons trap on NaN under fpe_trap_invalid.
    auto bad = [&] (amrex::Real x, const char* key, const char* range, auto ok) {
        if (std::isnan(x) || !std::isfinite(x) || !ok()) {
            err.emplace_back(std::string("turb.") + key + " = " + std::to_string(x) + ": must be " + range);
        }
    };
    bad(opt.Cmu0, "Cmu0", "finite and > 0", [&] { return opt.Cmu0 > 0; });
    bad(opt.Cb, "Cb", "finite and > 0", [&] { return opt.Cb > 0; });
    bad(opt.Rt_crit, "Rt_crit", "finite and <= 0", [&] { return opt.Rt_crit <= 0; });
    bad(opt.Rt_min, "Rt_min", "finite, < turb.Rt_crit and > -3.6 (the stability functions have poles near -3.61)",
        [&] { return !std::isnan(opt.Rt_crit) && opt.Rt_min < opt.Rt_crit && opt.Rt_min > amrex::Real(-3.6); });
    bad(opt.max_geom_lscale, "max_geom_lscale", "finite and > 0 (m)", [&] { return opt.max_geom_lscale > 0; });
    bad(opt.tke_floor, "tke_floor", "0 (machine epsilon) or >= machine epsilon (m^2/s^2)",
        [&] { return opt.tke_floor == 0 || opt.tke_floor >= std::numeric_limits<amrex::Real>::epsilon(); });
    bad(opt.sigma_k, "sigma_k", "finite and > 0", [&] { return opt.sigma_k > 0; });
    bad(opt.nu_molecular, "nu_molecular", "finite and >= 0 (m^2/s)", [&] { return opt.nu_molecular >= 0; });
    bad(opt.alpha_molecular, "alpha_molecular", "finite and >= 0 (m^2/s)", [&] { return opt.alpha_molecular >= 0; });
    bad(opt.kappa, "kappa", "finite and > 0", [&] { return opt.kappa > 0; });
    bad(opt.nut_relax, "nut_relax", "finite, > 0 and <= 1", [&] { return opt.nut_relax > 0 && opt.nut_relax <= 1; });
    return err;
}

ClosureOptions read_closure_options ()
{
    ClosureOptions opt;
    amrex::ParmParse pp("turb");
    for (const char* key : kKeys) {
        if (pp.countname(key) > 1) { amrex::Abort(std::string("ablsteady: turb.") + key + " is given more than once"); }
        if (pp.contains(key) && pp.countval(key) != 1) {
            amrex::Abort(std::string("ablsteady: turb.") + key + " takes exactly one value");
        }
    }
    pp.query("Cmu0", opt.Cmu0);
    pp.query("Cb", opt.Cb);
    pp.query("Rt_crit", opt.Rt_crit);
    pp.query("Rt_min", opt.Rt_min);
    pp.query("max_geom_lscale", opt.max_geom_lscale);
    pp.query("tke_floor", opt.tke_floor);
    pp.query("sigma_k", opt.sigma_k);
    pp.query("nu_molecular", opt.nu_molecular);
    pp.query("alpha_molecular", opt.alpha_molecular);
    pp.query("kappa", opt.kappa);
    pp.query("nut_relax", opt.nut_relax);
    int frozen = 0;
    pp.query("frozen_switch", frozen);
    if (frozen != 0 && frozen != 1) { amrex::Abort("ablsteady: turb.frozen_switch must be 0 or 1"); }
    opt.frozen_switch = (frozen == 1);
    std::string fm = "arithmetic", wd = "terrain_height";
    pp.query("face_mean", fm);
    pp.query("wall_distance", wd);
    if (fm == "arithmetic") { opt.face_mean = FaceMean::arithmetic; }
    else if (fm == "logarithmic") { opt.face_mean = FaceMean::logarithmic; }
    else { amrex::Abort("ablsteady: turb.face_mean = " + fm + ": must be arithmetic or logarithmic"); }
    if (wd == "terrain_height") { opt.wall_distance = WallDistanceKind::terrain_height; }
    else if (wd == "tucker") { opt.wall_distance = WallDistanceKind::tucker; }
    else { amrex::Abort("ablsteady: turb.wall_distance = " + wd + ": must be terrain_height or tucker"); }
    return checked(opt);
}

ClosureConstants closure_constants (const ClosureOptions& opt)
{
    ClosureConstants c;
    c.kappa = opt.kappa;
    c.Cmu0 = opt.Cmu0;
    c.Cmu0_pow3 = opt.Cmu0 * opt.Cmu0 * opt.Cmu0;
    c.inv_Cb_sq = amrex::Real(1) / (opt.Cb * opt.Cb);
    c.Rt_crit = opt.Rt_crit;
    c.Rt_min = opt.Rt_min;
    c.l_g_max = opt.max_geom_lscale;
    // ERF: the floor is erf.tke_floor if set, otherwise machine epsilon.
    c.tke_floor = amrex::max(opt.tke_floor, std::numeric_limits<amrex::Real>::epsilon());
    c.n2_eps = std::numeric_limits<amrex::Real>::epsilon();
    return c;
}

void ClosureFields::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    for (auto* mf : {&nu, &gamma_theta, &gamma_tke, &nut}) { mf->define(ba, dm, 1, 1); }
    for (auto* mf : {&length, &Rt, &Ps, &Pb, &eps, &dSdk, &dSdk_exact, &dnutdk}) { mf->define(ba, dm, 1, 0); }
}

void height_above_ground (const amrex::MultiFab& z_nd, amrex::MultiFab& dist)
{
    const amrex::BoxArray cba = amrex::convert(z_nd.boxArray(), amrex::IntVect::TheZeroVector());
    if (dist.empty()) { dist.define(cba, z_nd.DistributionMap(), 1, 1); }
    if (dist.boxArray() != cba || dist.DistributionMap() != z_nd.DistributionMap() || dist.nGrow() < 1) {
        amrex::Abort("ablsteady: height_above_ground: dist must be cell-centred on the nodes' layout with one ghost cell");
    }
    // Every box spans the column (State.cpp's layout): the ground is each
    // box's first layer and the vertical ghosts are outside the domain.
    const amrex::Box dom = cba.minimalBox();
    for (int n = 0; n < cba.size(); ++n) {
        if (cba[n].smallEnd(2) != dom.smallEnd(2) || cba[n].bigEnd(2) != dom.bigEnd(2)) {
            amrex::Abort("ablsteady: height_above_ground: every box must span the whole z column");
        }
    }
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = lateral_grown(mfi.validbox());
        const int klo = mfi.validbox().smallEnd(2);   // the ground layer (boxes span the column)
        const auto d = dist.array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real zc = 0, zg = 0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    zg += z(i + a, j + b, klo);
                    zc += z(i + a, j + b, k) + z(i + a, j + b, k + 1);
                }
            }
            d(i, j, k) = amrex::Real(0.125) * zc - amrex::Real(0.25) * zg;
        });
    }
    copy_vertical_ghosts(dist);
}

void stratification_and_strain (const State& q, const BuoyancyOptions& buoyancy, const amrex::Geometry& geom,
                                amrex::MultiFab& N2, amrex::MultiFab& S2)
{
    using namespace stencil;
    const auto dxinv = geom.InvCellSizeArray();
    const amrex::Real g_over_theta0 = buoyancy.gravity / buoyancy.theta0;
    if (N2.empty()) { N2.define(q.p.boxArray(), q.p.DistributionMap(), 1, 1); }
    if (S2.empty()) { S2.define(q.p.boxArray(), q.p.DistributionMap(), 1, 1); }
    for (amrex::MFIter mfi(N2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = lateral_grown(mfi.validbox());
        const auto n2 = N2.array(mfi);
        const auto s2 = S2.array(mfi);
        const auto th = q.theta.const_array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Metric m = metric_at<CC>(z, i, j, k, dxinv);
            const auto gt = grad_at<CC, CC>(th, i, j, k, m, dxinv);
            n2(i, j, k) = g_over_theta0 * gt[2];
            const auto gu = grad_at<XF, CC>(u, i, j, k, m, dxinv);
            const auto gv = grad_at<YF, CC>(v, i, j, k, m, dxinv);
            const auto gw = grad_at<ZF, CC>(w, i, j, k, m, dxinv);
            // 2 S_ij S_ij = 1/2 sum_ij (d_j u_i + d_i u_j)^2
            const amrex::Real G[3][3] = {{gu[0], gu[1], gu[2]}, {gv[0], gv[1], gv[2]}, {gw[0], gw[1], gw[2]}};
            amrex::Real s = 0;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const amrex::Real e = G[a][b] + G[b][a];
                    s += e * e;
                }
            }
            s2(i, j, k) = amrex::Real(0.5) * s;
        });
    }
    copy_vertical_ghosts(N2);
    copy_vertical_ghosts(S2);
}

Closure::Closure (const ClosureOptions& opt) : m_opt(checked(opt)), m_c(closure_constants(m_opt)) {}

void Closure::evaluate (const State& q, const amrex::MultiFab& dist, amrex::Real z0, const BuoyancyOptions& buoyancy,
                        const Mesh& mesh, const amrex::Geometry& geom, ClosureFields& f) const
{
    const char* op = "Closure::evaluate";
    check_operator_field(op, "tke", q.tke, Loc::cell, mesh, 2);
    check_operator_field(op, "theta", q.theta, Loc::cell, mesh, 2);
    check_operator_field(op, "dist", dist, Loc::cell, mesh, 1);
    if (std::isnan(z0) || !(z0 >= 0 && std::isfinite(z0))) {
        amrex::Abort(std::string("ablsteady: ") + op + ": z0 = " + std::to_string(z0) + " must be finite and >= 0");
    }
    for (const auto& e : check_buoyancy_options(buoyancy)) { amrex::Abort(std::string("ablsteady: ") + op + ": " + e); }
    // l_g = kappa (d + z0) must be positive everywhere the kernels run (1/l_g^2).
    if (!(dist.min(0, 1) + z0 > 0)) {
        amrex::Abort(std::string("ablsteady: ") + op + ": the wall distance plus z0 must be > 0 at every cell (min "
                     + std::to_string(dist.min(0, 1)) + " + " + std::to_string(z0) + ")");
    }
    if (f.empty()) { f.define(q.p.boxArray(), q.p.DistributionMap()); }
    check_operator_field(op, "nu", f.nu, Loc::cell, mesh, 1);
    if (has_history()) { check_operator_field(op, "history", m_nut_prev, Loc::cell, mesh, 1); }
    if (frozen() && (m_switches.boxArray() != q.p.boxArray() || m_switches.DistributionMap() != q.p.DistributionMap()
                     || m_switches.nGrow() < 1)) {
        amrex::Abort(std::string("ablsteady: ") + op + ": the frozen switches are not on this layout");
    }

    amrex::MultiFab N2, S2;
    stratification_and_strain(q, buoyancy, geom, N2, S2);
    const ClosureConstants c = m_c;
    const amrex::Real nu_mol = m_opt.nu_molecular, alpha_mol = m_opt.alpha_molecular, inv_sigma = 1 / m_opt.sigma_k;
    const bool has = has_history();
    const amrex::Real w = has ? m_opt.nut_relax : amrex::Real(1);
    const bool use_frozen = frozen();
    for (amrex::MFIter mfi(f.nu, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = lateral_grown(mfi.validbox());
        const amrex::Box vb = mfi.validbox();
        const auto kk = q.tke.const_array(mfi);
        const auto n2 = N2.const_array(mfi);
        const auto s2 = S2.const_array(mfi);
        const auto d = dist.const_array(mfi);
        // A valid, finite array in every case: the history when there is one,
        // else k itself with the weight 1 - w = 0 (0 times a finite value is
        // 0 exactly; a fresh field could hold signalling NaNs).
        const auto prev = has ? m_nut_prev.const_array(mfi) : q.tke.const_array(mfi);
        const auto sw = use_frozen ? m_switches.const_array(mfi) : amrex::Array4<int const>{};
        const auto nu = f.nu.array(mfi), gth = f.gamma_theta.array(mfi), gk = f.gamma_tke.array(mfi), nut = f.nut.array(mfi);
        const auto Lf = f.length.array(mfi), Rtf = f.Rt.array(mfi), Psf = f.Ps.array(mfi), Pbf = f.Pb.array(mfi),
                   epsf = f.eps.array(mfi), dsf = f.dSdk.array(mfi), dse = f.dSdk_exact.array(mfi), dnf = f.dnutdk.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real l_g = closure::geom_length(c.kappa * (d(i, j, k) + z0), c.l_g_max);
            const Switches s = use_frozen ? decode(sw(i, j, k)) : closure::switches_for(kk(i, j, k), n2(i, j, k), l_g, c);
            const auto p = closure::closure_point<amrex::Real>(kk(i, j, k), n2(i, j, k), s2(i, j, k), l_g, c, s);
            const amrex::Real ratio = p.cmu_prime / p.cmu;   // 1/Pr_t(Rt)
            const amrex::Real nt = w * p.nut + (1 - w) * prev(i, j, k);
            nut(i, j, k) = nt;
            nu(i, j, k) = nu_mol + nt;
            gth(i, j, k) = alpha_mol + ratio * nt;
            gk(i, j, k) = nu_mol + nt * inv_sigma;
            if (vb.contains(amrex::IntVect(i, j, k))) {
                // The sources with the relaxed nu_t: P_s = nu_t S^2, P_b = -nu_t' N^2.
                const amrex::Real Ps = nt * s2(i, j, k), Pb = -ratio * nt * n2(i, j, k);
                Lf(i, j, k) = p.L;
                Rtf(i, j, k) = p.Rt;
                Psf(i, j, k) = Ps;
                Pbf(i, j, k) = Pb;
                epsf(i, j, k) = p.eps;
                closure::Point<amrex::Real> pr = p;   // the relaxed sources for the diagonal
                pr.Ps = Ps;
                pr.Pb = Pb;
                dsf(i, j, k) = closure::patankar_diagonal(pr);
                // The exact derivative on dual numbers: nu_t relaxed by w, eps as is.
                const auto pd = closure::closure_point<Dual>(Dual{kk(i, j, k), 1}, n2(i, j, k), s2(i, j, k), l_g, c, s);
                const amrex::Real dnut = w * pd.nut.d;
                dnf(i, j, k) = dnut;
                dse(i, j, k) = dnut * s2(i, j, k) - ratio * dnut * n2(i, j, k) - pd.eps.d;
            }
        });
    }
    for (auto* mf : {&f.nu, &f.gamma_theta, &f.gamma_tke, &f.nut}) { copy_vertical_ghosts(*mf); }
}

void Closure::accept (const ClosureFields& f)
{
    if (f.empty()) { amrex::Abort("ablsteady: Closure::accept: the fields are not defined"); }
    if (m_nut_prev.empty() || m_nut_prev.boxArray() != f.nut.boxArray() || m_nut_prev.DistributionMap() != f.nut.DistributionMap()) {
        m_nut_prev.clear();
        m_nut_prev.define(f.nut.boxArray(), f.nut.DistributionMap(), 1, 1);
    }
    amrex::MultiFab::Copy(m_nut_prev, f.nut, 0, 0, 1, 1);
}

void Closure::write_checkpoint (const std::string& dir) const
{
    if (!has_history()) { return; }
    amrex::VisMF::Write(m_nut_prev, dir + "/Level_0/closure_nut");
}

void Closure::read_checkpoint (const std::string& dir, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                               const amrex::Geometry& geom)
{
    const std::string name = dir + "/Level_0/closure_nut";
    // One decision for every rank (the file system may show the directory late).
    int present = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(name + "_H") ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&present, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (present == 0) {
        m_nut_prev.clear();
        return;
    }
    amrex::MultiFab src;
    amrex::VisMF::Read(src, name);   // defined on the checkpoint's own layout, ghosts included
    m_nut_prev.clear();
    m_nut_prev.define(ba, dm, 1, 1);
    // Valid cells and the ghost ring (VisMF wrote both): a ghost outside the
    // domain is the same value in every box that holds that index (all
    // computed from the state's filled ghosts), so copying ghosts from a
    // checkpoint on other boxes restores it exactly. Ghosts first, valid
    // data last, as the state's restart does.
    m_nut_prev.ParallelCopy(src, 0, 0, 1, amrex::IntVect(1), amrex::IntVect(1), geom.periodicity());
    m_nut_prev.ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(0), geom.periodicity());
    (void)geom;
}

void Closure::freeze_switches (const State& q, const amrex::MultiFab& dist, amrex::Real z0, const BuoyancyOptions& buoyancy,
                               const Mesh& mesh, const amrex::Geometry& geom)
{
    check_operator_field("Closure::freeze_switches", "dist", dist, Loc::cell, mesh, 1);
    amrex::MultiFab N2, S2;
    stratification_and_strain(q, buoyancy, geom, N2, S2);
    m_switches.clear();
    m_switches.define(q.p.boxArray(), q.p.DistributionMap(), 1, 1);
    m_switches.setVal(0);
    const ClosureConstants c = m_c;
    for (amrex::MFIter mfi(m_switches, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = lateral_grown(mfi.validbox());
        const auto sw = m_switches.array(mfi);
        const auto kk = q.tke.const_array(mfi);
        const auto n2 = N2.const_array(mfi);
        const auto d = dist.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real l_g = closure::geom_length(c.kappa * (d(i, j, k) + z0), c.l_g_max);
            sw(i, j, k) = encode(closure::switches_for(kk(i, j, k), n2(i, j, k), l_g, c));
        });
    }
}

} // namespace ablsteady
