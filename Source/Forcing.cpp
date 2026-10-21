#include "Forcing.H"

#include "CVFluxes.H"
#include "Stencil_K.H"

#include <AMReX_GpuContainers.H>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace gust {

// NaN is tested with isnan before any ordered comparison: gcc's ordered
// comparisons raise the invalid-operation trap on a NaN.

std::vector<std::string> check_buoyancy_options (const BuoyancyOptions& opt)
{
    std::vector<std::string> err;
    if (std::isnan(opt.gravity) || !(opt.gravity >= 0 && std::isfinite(opt.gravity))) {
        err.emplace_back("buoyancy gravity = " + std::to_string(opt.gravity) + ": must be finite and >= 0");
    }
    if (std::isnan(opt.theta0) || !(opt.theta0 > 0 && std::isfinite(opt.theta0))) {
        err.emplace_back("buoyancy theta0 = " + std::to_string(opt.theta0) + ": must be finite and > 0");
    }
    return err;
}

std::vector<std::string> check_reference_profile (const ReferenceProfile& prof)
{
    std::vector<std::string> err;
    if (prof.z.empty() || prof.z.size() != prof.theta.size()) {
        err.emplace_back("reference profile: z and theta need the same number of points, at least one");
        return err;
    }
    for (int n = 0; n < static_cast<int>(prof.z.size()); ++n) {
        const amrex::Real z = prof.z[n], t = prof.theta[n];
        if (!std::isfinite(z) || std::isnan(t) || !(t > 0 && std::isfinite(t))) {
            err.emplace_back("reference profile point " + std::to_string(n) + ": z must be finite and theta finite and > 0");
            return err;
        }
        if (n > 0 && !(z > prof.z[n - 1])) {
            err.emplace_back("reference profile: heights must increase strictly (point " + std::to_string(n) + ")");
            return err;
        }
    }
    return err;
}

std::vector<std::string> check_coriolis_options (const CoriolisOptions& opt)
{
    std::vector<std::string> err;
    if (std::isnan(opt.latitude) || !(opt.latitude >= -90 && opt.latitude <= 90)) {
        err.emplace_back("coriolis latitude = " + std::to_string(opt.latitude) + ": must be in [-90, 90] degrees");
    }
    if (std::isnan(opt.rotation_period) || !(opt.rotation_period > 0 && std::isfinite(opt.rotation_period))) {
        err.emplace_back("coriolis rotation_period = " + std::to_string(opt.rotation_period)
                         + ": must be finite and > 0");
    }
    if (!std::isfinite(opt.ug) || !std::isfinite(opt.vg)) {
        err.emplace_back("coriolis ug, vg = " + std::to_string(opt.ug) + ", " + std::to_string(opt.vg)
                         + ": must be finite");
    }
    if (opt.geostrophic && !opt.enabled) {
        err.emplace_back("coriolis geostrophic forcing needs Coriolis enabled (it is f times the geostrophic wind)");
    } else if (opt.geostrophic && opt.latitude == 0) {
        err.emplace_back("coriolis geostrophic forcing at latitude 0: f = 0, so it would not drive the flow");
    }
    if (!opt.geostrophic && (opt.ug != 0 || opt.vg != 0)) {
        err.emplace_back("coriolis ug, vg are set but geostrophic forcing is off; they would be ignored");
    }
    return err;
}

amrex::Real coriolis_parameter (const CoriolisOptions& opt)
{
    if (!opt.enabled) { return 0; }
    const double pi = 3.14159265358979323846;
    return static_cast<amrex::Real>(4 * pi * std::sin(static_cast<double>(opt.latitude) * pi / 180)
                                    / static_cast<double>(opt.rotation_period));
}

void fill_theta_ref (const ReferenceProfile& prof, const amrex::MultiFab& z_nd, amrex::MultiFab& theta_ref)
{
    const auto err = check_reference_profile(prof);
    if (!err.empty()) { amrex::Abort("gust: fill_theta_ref: " + err.front()); }
    if (theta_ref.boxArray() != amrex::convert(z_nd.boxArray(), amrex::IntVect(0))
        || theta_ref.DistributionMap() != z_nd.DistributionMap() || theta_ref.nGrow() > z_nd.nGrow()) {
        amrex::Abort("gust: fill_theta_ref: theta_ref must be cell-centred on z_nd's cells, ghosts <= z_nd's");
    }
    const int n = static_cast<int>(prof.z.size());
    amrex::Gpu::DeviceVector<amrex::Real> zs(n), ts(n);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, prof.z.begin(), prof.z.end(), zs.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, prof.theta.begin(), prof.theta.end(), ts.begin());
    const amrex::Real* zp = zs.data();
    const amrex::Real* tp = ts.data();
    for (amrex::MFIter mfi(theta_ref, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.fabbox();
        const auto t = theta_ref.array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real zc = stencil::interp<stencil::NODE, stencil::CC>(z, i, j, k);
            amrex::Real v = tp[0];
            if (zc >= zp[n - 1]) {
                v = tp[n - 1];
            } else {
                for (int m = 1; m < n; ++m) {
                    if (zc < zp[m]) {
                        if (zc > zp[m - 1]) {
                            const amrex::Real w = (zc - zp[m - 1]) / (zp[m] - zp[m - 1]);
                            v = tp[m - 1] + w * (tp[m] - tp[m - 1]);
                        } else {
                            v = tp[m - 1];
                        }
                        break;
                    }
                }
            }
            t(i, j, k) = v;
        });
    }
    amrex::Gpu::streamSynchronize();
}

void subtract_buoyancy (const amrex::MultiFab& theta, const amrex::MultiFab& theta_ref, const BuoyancyOptions& opt,
                        const Mesh& mesh, amrex::MultiFab& Rw)
{
    const auto err = check_buoyancy_options(opt);
    if (!err.empty()) { amrex::Abort("gust: subtract_buoyancy: " + err.front()); }
    check_operator_field("subtract_buoyancy", "theta", theta, Loc::cell, mesh, 1);
    check_operator_field("subtract_buoyancy", "theta_ref", theta_ref, Loc::cell, mesh, 1);
    check_operator_field("subtract_buoyancy", "Rw", Rw, Loc::zface, mesh, 0);
    const amrex::Real gt = opt.gravity / opt.theta0;
    for (amrex::MFIter mfi(Rw, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto r = Rw.array(mfi);
        const auto t = theta.const_array(mfi);
        const auto tr = theta_ref.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real half = amrex::Real(0.5);
            r(i, j, k) -= gt * half * ((t(i, j, k - 1) - tr(i, j, k - 1)) + (t(i, j, k) - tr(i, j, k)));
        });
    }
}

void subtract_coriolis (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                        const CoriolisOptions& opt, const Mesh& mesh, amrex::MultiFab& Ru, amrex::MultiFab& Rv,
                        amrex::MultiFab& Rw)
{
    const auto err = check_coriolis_options(opt);
    if (!err.empty()) { amrex::Abort("gust: subtract_coriolis: " + err.front()); }
    if (!opt.enabled) { return; }
    const char* op = "subtract_coriolis";
    check_operator_field(op, "u", u, Loc::xface, mesh, 1);
    check_operator_field(op, "v", v, Loc::yface, mesh, 1);
    check_operator_field(op, "w", w, Loc::zface, mesh, 1);
    check_operator_field(op, "Ru", Ru, Loc::xface, mesh, 0);
    check_operator_field(op, "Rv", Rv, Loc::yface, mesh, 0);
    check_operator_field(op, "Rw", Rw, Loc::zface, mesh, 0);
    const double pi = 3.14159265358979323846;
    const double phi = static_cast<double>(opt.latitude) * pi / 180;
    const double fc = 4 * pi / static_cast<double>(opt.rotation_period);
    const auto fs = static_cast<amrex::Real>(fc * std::sin(phi));
    const auto fcos = static_cast<amrex::Real>(opt.three_d ? fc * std::cos(phi) : 0.0);
    const amrex::Real ug = opt.geostrophic ? opt.ug : amrex::Real(0);
    const amrex::Real vg = opt.geostrophic ? opt.vg : amrex::Real(0);
    using stencil::XF;
    using stencil::YF;
    using stencil::ZF;
    for (amrex::MFIter mfi(Ru, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto r = Ru.array(mfi);
        const auto vv = v.const_array(mfi);
        const auto ww = w.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) -= fs * (stencil::interp<YF, XF>(vv, i, j, k) - vg) - fcos * stencil::interp<ZF, XF>(ww, i, j, k);
        });
    }
    for (amrex::MFIter mfi(Rv, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto r = Rv.array(mfi);
        const auto uu = u.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) += fs * (stencil::interp<XF, YF>(uu, i, j, k) - ug);
        });
    }
    if (opt.three_d) {
        for (amrex::MFIter mfi(Rw, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto r = Rw.array(mfi);
            const auto uu = u.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                r(i, j, k) -= fcos * stencil::interp<XF, ZF>(uu, i, j, k);
            });
        }
    }
}

std::vector<std::string> check_forcing_options (const ForcingOptions& opt)
{
    std::vector<std::string> err;
    for (const auto& e : check_buoyancy_options(opt.buoyancy)) { err.push_back("forcing." + e.substr(e.find(' ') + 1)); }
    for (const auto& e : check_coriolis_options(opt.coriolis)) { err.push_back("forcing." + e.substr(e.find(' ') + 1)); }
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(opt.pressure_gradient[d])) {
            err.emplace_back("forcing.pressure_gradient = " + std::to_string(opt.pressure_gradient[d])
                             + " (component " + std::to_string(d) + "): must be finite (m/s^2)");
            break;
        }
    }
    return err;
}

namespace {

const char* const kForcingKeys[] = {"gravity", "theta0", "coriolis", "latitude", "rotation_period", "coriolis_3d",
                                    "geostrophic", "ug", "vg", "pressure_gradient"};

} // namespace

ForcingOptions read_forcing_options ()
{
    ForcingOptions opt;
    amrex::ParmParse pp("forcing");
    std::vector<std::string> err;
    for (const char* key : kForcingKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("forcing.") + key + " is given more than once"); }
        const int n = (std::string(key) == "pressure_gradient") ? 3 : 1;
        if (pp.countval(key) != n) {
            err.emplace_back(std::string("forcing.") + key + " takes exactly " + std::to_string(n) + " value" + (n > 1 ? "s" : ""));
        }
    }
    if (err.empty()) {
        pp.query("gravity", opt.buoyancy.gravity);
        pp.query("theta0", opt.buoyancy.theta0);
        auto flag = [&] (const char* key, bool& b) {
            int v = b ? 1 : 0;
            if (pp.query(key, v)) {
                if (v != 0 && v != 1) { err.emplace_back(std::string("forcing.") + key + " must be 0 or 1"); }
                b = (v == 1);
            }
        };
        flag("coriolis", opt.coriolis.enabled);
        flag("coriolis_3d", opt.coriolis.three_d);
        flag("geostrophic", opt.coriolis.geostrophic);
        pp.query("latitude", opt.coriolis.latitude);
        pp.query("rotation_period", opt.coriolis.rotation_period);
        pp.query("ug", opt.coriolis.ug);
        pp.query("vg", opt.coriolis.vg);
        amrex::Vector<amrex::Real> pg;
        if (pp.queryarr("pressure_gradient", pg) && pg.size() == 3) {
            for (int d = 0; d < 3; ++d) { opt.pressure_gradient[d] = pg[d]; }
        }
        for (const auto& e : check_forcing_options(opt)) { err.push_back(e); }
        auto given = [&] (const char* k) { return std::find(opt.given.begin(), opt.given.end(), k) != opt.given.end(); };
        if (opt.buoyancy.gravity == 0 && given("theta0")) {
            err.emplace_back("forcing.theta0 has no effect with forcing.gravity = 0");
        }
        if (!opt.coriolis.enabled) {
            for (const char* k : {"latitude", "rotation_period", "coriolis_3d", "geostrophic", "ug", "vg"}) {
                if (given(k)) { err.emplace_back(std::string("forcing.") + k + " has no effect without forcing.coriolis = 1"); }
            }
        } else if (!opt.coriolis.geostrophic) {
            for (const char* k : {"ug", "vg"}) {
                if (given(k)) { err.emplace_back(std::string("forcing.") + k + " has no effect without forcing.geostrophic = 1"); }
            }
        }
    }
    if (!err.empty()) {
        std::string all = "gust: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

void subtract_pressure_gradient (const amrex::Array<amrex::Real, 3>& F, const Mesh& mesh, amrex::MultiFab& Ru,
                                 amrex::MultiFab& Rv, amrex::MultiFab& Rw)
{
    const char* op = "subtract_pressure_gradient";
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(F[d])) { amrex::Abort(std::string("gust: ") + op + ": the force must be finite"); }
    }
    check_operator_field(op, "Ru", Ru, Loc::xface, mesh, 0);
    check_operator_field(op, "Rv", Rv, Loc::yface, mesh, 0);
    check_operator_field(op, "Rw", Rw, Loc::zface, mesh, 0);
    amrex::MultiFab* R[3] = {&Ru, &Rv, &Rw};
    for (int d = 0; d < 3; ++d) {
        if (F[d] != 0) { R[d]->plus(-F[d], 0, 1, 0); }
    }
}

} // namespace gust
