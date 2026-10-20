#include "SurfaceLayer.H"

#include "CVFluxes.H"
#include "ERF_MOSTUtils.H"
#include "ImplicitOperator.H"
#include "Stencil_K.H"

#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>
#include <string>

namespace ablsteady {

namespace {

constexpr amrex::Real kWsMin = amrex::Real(0.1);        // ERF's minimum wind speed
constexpr amrex::Real kUstarMin = amrex::Real(0.001);   // ERF's u* floor in L
constexpr amrex::Real kWallKUstarMin = amrex::Real(0.01);
constexpr int kMaxNewton = 50;

amrex::BoxArray ground_slab (const amrex::BoxArray& ba)
{
    amrex::BoxList bl;
    for (int n = 0; n < ba.size(); ++n) {
        amrex::Box b = ba[n];
        b.setBig(2, b.smallEnd(2));
        bl.push_back(b);
    }
    return amrex::BoxArray(bl);
}

void check_slab (const char* op, const SurfaceState& s, const amrex::MultiFab& like)
{
    if (s.ustar.boxArray() != ground_slab(amrex::convert(like.boxArray(), amrex::IntVect(0)))
        || s.ustar.DistributionMap() != like.DistributionMap()) {
        amrex::Abort(std::string("ablsteady: ") + op + ": the surface state is not on this layout");
    }
}

// Vertical height of the centre of cell (i, j, k) above the column's ground
// (the mean of its four ground nodes, at the domain's low z, klo).
struct Column
{
    amrex::Array4<amrex::Real const> z;
    int i, j, klo;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real ground () const noexcept
    {
        return amrex::Real(0.25) * ((z(i, j, klo) + z(i + 1, j, klo)) + (z(i, j + 1, klo) + z(i + 1, j + 1, klo)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int k) const noexcept
    {
        return stencil::interp<stencil::NODE, stencil::CC>(z, i, j, klo + k) - ground();
    }
};

// The wind parallel to the ground at the centre of cell (i, j, k).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tangential_wind (amrex::Array4<amrex::Real const> const& u, amrex::Array4<amrex::Real const> const& v,
                      amrex::Array4<amrex::Real const> const& w, int i, int j, int k, const amrex::Real n[3],
                      amrex::Real ut[3]) noexcept
{
    const amrex::Real half = amrex::Real(0.5);
    const amrex::Real U[3] = {half * (u(i, j, k) + u(i + 1, j, k)), half * (v(i, j, k) + v(i, j + 1, k)),
                              half * (w(i, j, k) + w(i, j, k + 1))};
    const amrex::Real un = U[0] * n[0] + U[1] * n[1] + U[2] * n[2];
    for (int d = 0; d < 3; ++d) { ut[d] = U[d] - un * n[d]; }
}

// u - kappa U / den(u): the MOST relation for u*, den = max(ln((d+z0)/z0) - psi_m(d/L), m),
// m the floor of bc.surface_denominator_floor.
struct MostRelation
{
    amrex::Real kappa, U, d, lz, m, gq;
    amrex::Real ustar_min;   //!< kUstarMin, held by value: a namespace constant is not usable in device code
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real inv_L (amrex::Real u) const noexcept
    {
        const amrex::Real uf = amrex::max(u, ustar_min);
        return -kappa * gq / (uf * uf * uf);
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real u) const noexcept
    {
        const similarity_funs sfuns;
        const amrex::Real zeta = amrex::max(amrex::min(d * inv_L(u), amrex::Real(100)), -amrex::Real(100));
        return u - kappa * U / amrex::max(lz - sfuns.calc_psi_m(zeta), m);
    }
};

template <bool HEAT, bool FIRST_CENTRE, bool FLOOR_ONE>
void columns (const State& q, const SurfaceLayerOptions& opt, const BuoyancyOptions& buoy,
              const amrex::Geometry& geom, SurfaceState& s)
{
    const auto dxinv = geom.InvCellSizeArray();
    const int nz = geom.Domain().length(2);
    const int klo = geom.Domain().smallEnd(2);
    const amrex::Real zref_in = opt.zref, z0 = opt.z0, kappa = opt.kappa;
    const amrex::Real gq = (buoy.gravity / buoy.theta0) * opt.heat_flux;
    const amrex::Real inv_cmu2 = amrex::Real(1) / (opt.Cmu0 * opt.Cmu0);
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real tol = 64 * eps;
    const amrex::Real hrel = std::sqrt(eps);
    // The floors as locals: nvcc cannot use a namespace-scope constant in
    // device code (it is odr-used through amrex::max's references).
    const amrex::Real ws_min_c = kWsMin, ustar_min_c = kUstarMin;
    for (amrex::MFIter mfi(s.ustar, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(amrex::IntVect(1, 1, 0));
        const auto z = q.z_phys_nd.const_array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto us = s.ustar.array(mfi);
        const auto sp = s.speed.array(mfi);
        const auto zr_a = s.zref.array(mfi);
        const auto d0_a = s.d0.array(mfi);
        const auto il = s.inv_L.array(mfi);
        const auto ta = s.tau.array(mfi);
        const auto dg = s.drag.array(mfi);
        const auto kw = s.kwall.array(mfi);
        const auto ar = s.area.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            // Captured here, outside the constexpr-if blocks below: nvcc refuses
            // an extended device lambda whose first capture of a variable is
            // inside an if constexpr.
            const amrex::Real gq_c = gq, hrel_c = hrel, tol_c = tol, ustar_min_l = ustar_min_c;
            const Column h{z, i, j, klo};
            const amrex::Real zx = stencil::deriv<stencil::NODE, stencil::ZF, 0>(z, i, j, klo) * dxinv[0];
            const amrex::Real zy = stencil::deriv<stencil::NODE, stencil::ZF, 1>(z, i, j, klo) * dxinv[1];
            const amrex::Real area = std::sqrt(1 + zx * zx + zy * zy);
            const amrex::Real n[3] = {-zx / area, -zy / area, 1 / area};
            const amrex::Real h0 = h(0);
            const amrex::Real zr = FIRST_CENTRE ? h0 : zref_in;
            // Clamped: the host has refused a z_ref outside every column.
            const int k = amrex::max(0, amrex::min(nz - 2, reference_interval(h, nz, zr)));
            amrex::Real ua[3], ub[3];
            tangential_wind(u, v, w, i, j, klo + k, n, ua);
            tangential_wind(u, v, w, i, j, klo + k + 1, n, ub);
            // Interpolation in ln(d + z0), d = h / |A|.
            const amrex::Real ia = amrex::Real(1) / area;
            const amrex::Real la = std::log(h(k) * ia + z0);
            const amrex::Real f = (std::log(zr * ia + z0) - la) / (std::log(h(k + 1) * ia + z0) - la);
            amrex::Real ut[3];
            for (int c = 0; c < 3; ++c) { ut[c] = ua[c] + f * (ub[c] - ua[c]); }
            const amrex::Real mag = std::sqrt(ut[0] * ut[0] + ut[1] * ut[1] + ut[2] * ut[2]);
            const amrex::Real U = amrex::max(mag, ws_min_c);
            const amrex::Real d = zr * ia;
            const amrex::Real lz = std::log((d + z0) / z0);
            // The denominator's floor: min(lz, 1) leaves the neutral log law
            // whole (lz > 0 for every d > 0); 1 cuts it below d = (e - 1) z0.
            amrex::Real m = amrex::min(lz, amrex::Real(1));
            if constexpr (FLOOR_ONE) { m = amrex::Real(1); }
            amrex::Real ustar = kappa * U / amrex::max(lz, m);
            amrex::Real invL = 0;
            if constexpr (HEAT) {
                const MostRelation F{kappa, U, d, lz, m, gq_c, ustar_min_l};
                bool converged = false;
                for (int it = 0; it < kMaxNewton && !converged; ++it) {
                    const amrex::Real fu = F(ustar);
                    const amrex::Real du = hrel_c * ustar;
                    const amrex::Real fp = (F(ustar + du) - fu) / du;
                    amrex::Real next = (fp > 0) ? ustar - fu / fp : ustar - amrex::Real(0.5) * fu;
                    next = amrex::max(amrex::Real(0.5) * ustar, amrex::min(amrex::Real(2) * ustar, next));
                    converged = std::abs(next - ustar) <= tol_c * next;
                    ustar = next;
                }
                if (!converged) { ustar = std::numeric_limits<amrex::Real>::quiet_NaN(); }
                invL = F.inv_L(ustar);
            }
            const amrex::Real d0 = h0 * ia;
            us(i, j, 0) = ustar;
            sp(i, j, 0) = U;
            zr_a(i, j, 0) = zr;
            d0_a(i, j, 0) = d0;
            il(i, j, 0) = invL;
            ar(i, j, 0) = area;
            const amrex::Real c = ustar * ustar / U * area;
            dg(i, j, 0) = c;
            for (int m = 0; m < 3; ++m) { ta(i, j, 0, m) = c * ut[m]; }
            amrex::Real kwall = ustar * ustar * inv_cmu2;
            if constexpr (HEAT) {
                // Axell and Liungman 2001, Eq. 16, as ERF: only an upward
                // (destabilising) flux raises k_w.
                kwall = inv_cmu2 * std::pow(ustar * ustar * ustar + kappa * amrex::max(gq_c, amrex::Real(0)) * d0,
                                            amrex::Real(2) / amrex::Real(3));
            }
            kw(i, j, 0) = kwall;
        });
    }
}

} // namespace

std::vector<std::string> check_surface_layer_options (const SurfaceLayerOptions& opt)
{
    std::vector<std::string> err;
    // isnan first: gcc's ordered comparisons trap on NaN.
    auto positive = [&] (amrex::Real x, const char* name) {
        if (std::isnan(x) || !(x > 0 && std::isfinite(x))) {
            err.emplace_back(std::string("surface layer ") + name + " = " + std::to_string(x) + ": must be finite and > 0");
        }
    };
    positive(opt.z0, "z0");
    positive(opt.Cmu0, "Cmu0");
    positive(opt.kappa, "kappa");
    if (!std::isfinite(opt.zref)) {
        err.emplace_back("surface layer zref = " + std::to_string(opt.zref) + ": must be finite (<= 0: first cell centre)");
    }
    if (!std::isfinite(opt.heat_flux)) {
        err.emplace_back("surface layer heat_flux = " + std::to_string(opt.heat_flux) + ": must be finite");
    } else if (opt.stability == SurfaceLayerOptions::Stability::neutral && opt.heat_flux != 0) {
        err.emplace_back("surface layer heat_flux is set but the stability is neutral; it would be ignored");
    }
    return err;
}

SurfaceLayerOptions read_surface_layer_options ()
{
    // The ground is a boundary condition: its keys live in bc.surface_*.
    SurfaceLayerOptions opt;
    amrex::ParmParse pp("bc");
    for (const char* key : {"surface_z0", "surface_zref", "surface_stability", "surface_heat_flux", "surface_wall_k",
                            "surface_denominator_floor"}) {
        if (pp.contains(key) && pp.countval(key) != 1) {
            amrex::Abort(std::string("ablsteady: bc.") + key + " takes exactly one value");
        }
    }
    pp.query("surface_z0", opt.z0);
    pp.query("surface_zref", opt.zref);
    pp.query("surface_heat_flux", opt.heat_flux);
    std::string s = "neutral", w = "dirichlet", f = "log_law";
    pp.query("surface_stability", s);
    pp.query("surface_wall_k", w);
    pp.query("surface_denominator_floor", f);
    if (s == "neutral") { opt.stability = SurfaceLayerOptions::Stability::neutral; }
    else if (s == "heat_flux") { opt.stability = SurfaceLayerOptions::Stability::heat_flux; }
    else { amrex::Abort("ablsteady: bc.surface_stability = " + s + ": must be neutral or heat_flux"); }
    if (w == "dirichlet") { opt.wall_k = SurfaceLayerOptions::WallK::dirichlet; }
    else if (w == "source") { opt.wall_k = SurfaceLayerOptions::WallK::source; }
    else { amrex::Abort("ablsteady: bc.surface_wall_k = " + w + ": must be dirichlet or source"); }
    if (f == "log_law") { opt.denominator_floor = SurfaceLayerOptions::DenominatorFloor::log_law; }
    else if (f == "one") { opt.denominator_floor = SurfaceLayerOptions::DenominatorFloor::one; }
    else { amrex::Abort("ablsteady: bc.surface_denominator_floor = " + f + ": must be log_law or one"); }
    const auto err = check_surface_layer_options(opt);
    if (!err.empty()) {
        // The checks say "surface layer z0"; the key is bc.surface_z0.
        std::string all = "ablsteady: invalid inputs:";
        const std::string from = "surface layer ";
        for (std::string msg : err) {
            if (msg.rfind(from, 0) == 0) { msg = "bc.surface_" + msg.substr(from.size()); }
            all += "\n  " + msg;
        }
        amrex::Abort(all);
    }
    return opt;
}

void compute_surface_state (const State& q, const SurfaceLayerOptions& opt, const BuoyancyOptions& buoyancy,
                            const Mesh& mesh, const amrex::Geometry& geom, SurfaceState& s)
{
    const char* op = "compute_surface_state";
    auto err = check_surface_layer_options(opt);
    for (const auto& e : check_buoyancy_options(buoyancy)) { err.push_back(e); }
    if (!err.empty()) { amrex::Abort(std::string("ablsteady: ") + op + ": " + err.front()); }
    check_operator_field(op, "u", q.u, Loc::xface, mesh, 1);
    check_operator_field(op, "v", q.v, Loc::yface, mesh, 1);
    check_operator_field(op, "w", q.w, Loc::zface, mesh, 1);
    check_nodal_heights(op, q.z_phys_nd, mesh);
    const amrex::Box& dom = geom.Domain();
    const int nz = dom.length(2), klo = dom.smallEnd(2);
    if (nz < 2) { amrex::Abort("ablsteady: compute_surface_state: the mesh needs at least two layers"); }
    const amrex::BoxArray& cba = mesh.detJ.boxArray();
    for (int n = 0; n < cba.size(); ++n) {
        if (cba[n].smallEnd(2) != klo || cba[n].bigEnd(2) != dom.bigEnd(2)) {
            amrex::Abort("ablsteady: compute_surface_state: every box must span the whole z column");
        }
    }

    // A given z_ref must lie on every valid column, between the highest
    // first-centre height and the lowest last-centre height.
    if (opt.zref > 0) {
        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMin> ops;
        amrex::ReduceData<amrex::Real, amrex::Real> data(ops);
        for (amrex::MFIter mfi(mesh.detJ, false); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.validbox();
            bx.setBig(2, bx.smallEnd(2));
            const auto z = q.z_phys_nd.const_array(mfi);
            ops.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<amrex::Real, amrex::Real> {
                const Column h{z, i, j, klo};
                return {h(0), h(nz - 1)};
            });
        }
        const auto r = data.value();
        amrex::Real lo = amrex::get<0>(r), hi = amrex::get<1>(r);
        amrex::ParallelAllReduce::Max(lo, amrex::ParallelDescriptor::Communicator());
        amrex::ParallelAllReduce::Min(hi, amrex::ParallelDescriptor::Communicator());
        if (!(opt.zref >= lo && opt.zref <= hi)) {
            amrex::Abort("ablsteady: surface layer zref = " + std::to_string(opt.zref)
                         + " m must lie between the first and last cell-centre heights above the ground of every column ("
                         + std::to_string(lo) + " to " + std::to_string(hi) + " m)");
        }
    }

    const amrex::BoxArray ba2 = ground_slab(cba);
    const amrex::DistributionMapping& dm = mesh.detJ.DistributionMap();
    if (s.ustar.empty() || s.ustar.boxArray() != ba2 || s.ustar.DistributionMap() != dm) {
        const amrex::IntVect ng(1, 1, 0);
        for (auto* mf : {&s.ustar, &s.speed, &s.zref, &s.d0, &s.inv_L, &s.drag, &s.kwall, &s.area}) {
            mf->clear();
            mf->define(ba2, dm, 1, ng);
        }
        s.tau.clear();
        s.tau.define(ba2, dm, 3, ng);
    }
    const bool heat = (opt.stability == SurfaceLayerOptions::Stability::heat_flux);
    const bool first = !(opt.zref > 0);
    if (opt.denominator_floor == SurfaceLayerOptions::DenominatorFloor::one) {
        if (heat && first) { columns<true, true, true>(q, opt, buoyancy, geom, s); }
        else if (heat) { columns<true, false, true>(q, opt, buoyancy, geom, s); }
        else if (first) { columns<false, true, true>(q, opt, buoyancy, geom, s); }
        else { columns<false, false, true>(q, opt, buoyancy, geom, s); }
    } else {
        if (heat && first) { columns<true, true, false>(q, opt, buoyancy, geom, s); }
        else if (heat) { columns<true, false, false>(q, opt, buoyancy, geom, s); }
        else if (first) { columns<false, true, false>(q, opt, buoyancy, geom, s); }
        else { columns<false, false, false>(q, opt, buoyancy, geom, s); }
    }

    // Every column the ground fluxes read (ghost columns included) needs a
    // finite, positive u* and a finite k_w. NaN is tested with != first.
    amrex::ReduceOps<amrex::ReduceOpSum> ops;
    amrex::ReduceData<int> data(ops);
    for (amrex::MFIter mfi(s.ustar, false); mfi.isValid(); ++mfi) {
        const auto us = s.ustar.const_array(mfi);
        const auto kw = s.kwall.const_array(mfi);
        ops.eval(mfi.growntilebox(amrex::IntVect(1, 1, 0)), data,
                 [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<int> {
            const amrex::Real a = us(i, j, k), b = kw(i, j, k);
            if (a != a || b != b) { return {1}; }
            return {(a > 0 && a < std::numeric_limits<amrex::Real>::max()) ? 0 : 1};
        });
    }
    int bad = amrex::get<0>(data.value());
    amrex::ParallelAllReduce::Sum(bad, amrex::ParallelDescriptor::Communicator());
    if (bad > 0) {
        amrex::Abort("ablsteady: compute_surface_state: no positive, finite friction velocity at "
                     + std::to_string(bad) + " columns (surface layer heat_flux = " + std::to_string(opt.heat_flux)
                     + ", z0 = " + std::to_string(opt.z0) + ")");
    }
}

void set_ground_fluxes (const SurfaceState& s, const SurfaceLayerOptions& opt, amrex::MultiFab& Gu,
                        amrex::MultiFab& Gv, amrex::MultiFab& Gtheta, amrex::MultiFab& Gtke)
{
    check_slab("set_ground_fluxes", s, Gtheta);
    const amrex::Real qs = (opt.stability == SurfaceLayerOptions::Stability::heat_flux) ? opt.heat_flux : amrex::Real(0);
    const amrex::Real half = amrex::Real(0.5);
    for (amrex::MFIter mfi(Gtheta, false); mfi.isValid(); ++mfi) {
        amrex::Box cb = amrex::enclosedCells(mfi.validbox());
        cb.setBig(2, cb.smallEnd(2));
        const int k0 = cb.smallEnd(2);
        const auto ta = s.tau.const_array(mfi);
        const auto ar = s.area.const_array(mfi);
        const auto gu = Gu.array(mfi);
        const auto gv = Gv.array(mfi);
        const auto gt = Gtheta.array(mfi);
        const auto gk = Gtke.array(mfi);
        amrex::ParallelFor(amrex::surroundingNodes(cb, 0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            gu(i, j, k0) = half * (ta(i - 1, j, 0, 0) + ta(i, j, 0, 0));
        });
        amrex::ParallelFor(amrex::surroundingNodes(cb, 1), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            gv(i, j, k0) = half * (ta(i, j - 1, 0, 1) + ta(i, j, 0, 1));
        });
        amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            gt(i, j, k0) = -qs * ar(i, j, 0);
            gk(i, j, k0) = 0;
        });
    }
}

void apply_wall_k (const SurfaceState& s, const SurfaceLayerOptions& opt, const amrex::MultiFab& tke,
                   amrex::MultiFab& Rk)
{
    check_slab("apply_wall_k", s, Rk);
    const bool dirichlet = (opt.wall_k == SurfaceLayerOptions::WallK::dirichlet);
    for (amrex::MFIter mfi(Rk, false); mfi.isValid(); ++mfi) {
        amrex::Box cb = mfi.validbox();
        cb.setBig(2, cb.smallEnd(2));
        const int k0 = cb.smallEnd(2);
        const auto r = Rk.array(mfi);
        const auto t = tke.const_array(mfi);
        const auto us = s.ustar.const_array(mfi);
        const auto d0 = s.d0.const_array(mfi);
        const auto kw = s.kwall.const_array(mfi);
        const amrex::Real us_min = kWallKUstarMin;   // by value for device code
        if (dirichlet) {
            amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                r(i, j, k0) = amrex::max(us(i, j, 0), us_min) / d0(i, j, 0) * (t(i, j, k0) - kw(i, j, 0));
            });
        } else {
            amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                r(i, j, k0) += amrex::max(us(i, j, 0), us_min) / d0(i, j, 0) * (t(i, j, k0) - kw(i, j, 0));
            });
        }
    }
}

void apply_wall_to_implicit (const SurfaceState& s, Loc loc, const Mesh& mesh, const amrex::Geometry& geom,
                             amrex::MultiFab& A)
{
    check_slab("apply_wall_to_implicit", s, A);
    if (loc == Loc::zface) { amrex::Abort("ablsteady: apply_wall_to_implicit: w has no wall row (a boundary row)"); }
    const amrex::Real dzinv = geom.InvCellSize(2);
    const amrex::Real half = amrex::Real(0.5);
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        amrex::Box gb = mfi.validbox();
        gb.setBig(2, gb.smallEnd(2));
        const int k0 = gb.smallEnd(2);
        const auto a = A.array(mfi);
        const auto dg = s.drag.const_array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        // The ground neighbour's weight leaves the row.
        amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            a(i, j, k0, coef_P) -= a(i, j, k0, coef_ZM);
            a(i, j, k0, coef_ZM) = 0;
        });
        if (loc == Loc::xface) {
            amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                a(i, j, k0, coef_P) += half * (dg(i - 1, j, 0) + dg(i, j, 0)) * dzinv
                                     / stencil::interp<stencil::CC, stencil::XF>(J, i, j, k0);
            });
        } else if (loc == Loc::yface) {
            amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                a(i, j, k0, coef_P) += half * (dg(i, j - 1, 0) + dg(i, j, 0)) * dzinv
                                     / stencil::interp<stencil::CC, stencil::YF>(J, i, j, k0);
            });
        }
    }
}

void apply_wall_k_to_implicit (const SurfaceState& s, const SurfaceLayerOptions& opt, amrex::MultiFab& A)
{
    check_slab("apply_wall_k_to_implicit", s, A);
    if (A.nComp() != n_implicit_coefs || !A.ixType().cellCentered()) {
        amrex::Abort("ablsteady: apply_wall_k_to_implicit: A must be the cell-centred implicit coefficients");
    }
    const bool dirichlet = (opt.wall_k == SurfaceLayerOptions::WallK::dirichlet);
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        amrex::Box cb = mfi.validbox();
        cb.setBig(2, cb.smallEnd(2));
        const int k0 = cb.smallEnd(2);
        const auto a = A.array(mfi);
        const auto us = s.ustar.const_array(mfi);
        const auto d0 = s.d0.const_array(mfi);
        const amrex::Real us_min = kWallKUstarMin;   // by value for device code
        if (dirichlet) {
            amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k0, n) = 0; }
                a(i, j, k0, coef_P) = amrex::max(us(i, j, 0), us_min) / d0(i, j, 0);
            });
        } else {
            amrex::ParallelFor(cb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                a(i, j, k0, coef_P) += amrex::max(us(i, j, 0), us_min) / d0(i, j, 0);
            });
        }
    }
}

} // namespace ablsteady
