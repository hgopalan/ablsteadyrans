#include "InflowProfile.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

namespace ablsteady {

std::vector<std::string> check_inflow_profile_options (const InflowProfileOptions& opt)
{
    std::vector<std::string> err;
    // isnan first: gcc's ordered comparisons trap on NaN.
    auto positive = [&] (amrex::Real x, const char* name) {
        if (std::isnan(x) || !(x > 0 && std::isfinite(x))) {
            err.emplace_back(std::string("bc.inflow_") + name + " = " + std::to_string(x) + ": must be finite and > 0");
        }
    };
    using K = InflowProfileOptions::Kind;
    positive(opt.z0, "z0");   // every kind: the roughness of the equilibrium k and nu_t
    if (opt.kind == K::log_law) {
        positive(opt.ustar, "ustar");
    } else if (opt.kind == K::power_law) {
        if (std::isnan(opt.speed_ref) || !(opt.speed_ref >= 0 && std::isfinite(opt.speed_ref))) {
            err.emplace_back("bc.inflow_speed = " + std::to_string(opt.speed_ref) + ": must be finite and >= 0");
        }
        positive(opt.z_ref, "zref");
        positive(opt.alpha, "alpha");
    } else if (opt.file.empty()) {
        err.emplace_back("bc.inflow_file must be given for bc.inflow_profile = file");
    }
    if (!std::isfinite(opt.direction)) {
        err.emplace_back("bc.inflow_direction = " + std::to_string(opt.direction) + ": must be finite (degrees)");
    }
    positive(opt.Cmu0, "Cmu0");
    positive(opt.kappa, "kappa");
    return err;
}

InflowProfile::InflowProfile (const InflowProfileOptions& opt, const ReferenceProfile& theta_ref)
    : m_opt(opt)
{
    auto err = check_inflow_profile_options(opt);
    for (const auto& e : check_reference_profile(theta_ref)) { err.push_back(e); }
    if (!err.empty()) { amrex::Abort("ablsteady: InflowProfile: " + err.front()); }
    m_hrz = theta_ref.z;
    m_hrt = theta_ref.theta;
    m_rz.resize(theta_ref.z.size());
    m_rt.resize(theta_ref.theta.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, theta_ref.z.begin(), theta_ref.z.end(), m_rz.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, theta_ref.theta.begin(), theta_ref.theta.end(), m_rt.begin());
    if (opt.kind != InflowProfileOptions::Kind::file) { return; }

    // The file, read on every rank (a small text table).
    std::ifstream is(opt.file);
    if (!is) { amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": cannot be opened"); }
    amrex::Vector<amrex::Real> z, u, v, t, k;
    std::string line;
    int ncol = 0, lineno = 0;
    while (std::getline(is, line)) {
        ++lineno;
        const auto hash = line.find('#');
        if (hash != std::string::npos) { line.erase(hash); }
        std::istringstream ss(line);
        std::vector<double> vals;
        double x;
        while (ss >> x) { vals.push_back(x); }
        if (vals.empty()) { continue; }
        if (vals.size() != 4 && vals.size() != 5) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": line " + std::to_string(lineno)
                         + " must hold z u v theta [k] (4 or 5 values)");
        }
        if (ncol == 0) { ncol = static_cast<int>(vals.size()); }
        if (static_cast<int>(vals.size()) != ncol) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": line " + std::to_string(lineno)
                         + " has a different number of columns than the first row");
        }
        for (double y : vals) {
            if (!std::isfinite(y)) {
                amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": line " + std::to_string(lineno)
                             + " holds a value that is not finite");
            }
        }
        if (!z.empty() && !(vals[0] > static_cast<double>(z.back()))) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": heights must increase strictly (line "
                         + std::to_string(lineno) + ")");
        }
        if (z.empty() && vals[0] != 0) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": the first height must be 0 (the ground)");
        }
        if (vals[3] <= 0) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": theta must be > 0 (line "
                         + std::to_string(lineno) + ")");
        }
        if (ncol == 5 && vals[4] < 0) {
            amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": k must be >= 0 (line "
                         + std::to_string(lineno) + ")");
        }
        z.push_back(static_cast<amrex::Real>(vals[0]));
        u.push_back(static_cast<amrex::Real>(vals[1]));
        v.push_back(static_cast<amrex::Real>(vals[2]));
        t.push_back(static_cast<amrex::Real>(vals[3]));
        k.push_back(static_cast<amrex::Real>(ncol == 5 ? vals[4] : 0.0));
    }
    if (z.size() < 2) {
        amrex::Abort("ablsteady: bc.inflow_file = " + opt.file + ": needs at least two rows");
    }
    m_has_k = (ncol == 5);
    m_hz = z; m_hu = u; m_hv = v; m_ht = t; m_hk = k;
    for (auto* p : {&m_z, &m_u, &m_v, &m_t, &m_k}) { p->resize(z.size()); }
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, z.begin(), z.end(), m_z.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u.begin(), u.end(), m_u.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, v.begin(), v.end(), m_v.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, t.begin(), t.end(), m_t.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, k.begin(), k.end(), m_k.begin());
}

InflowProfile::View InflowProfile::view () const
{
    return view_with(m_z.data(), m_u.data(), m_v.data(), m_t.data(), m_k.data(), m_rz.data(), m_rt.data());
}

InflowProfile::View InflowProfile::host_view () const
{
    return view_with(m_hz.data(), m_hu.data(), m_hv.data(), m_ht.data(), m_hk.data(), m_hrz.data(), m_hrt.data());
}

InflowProfile::View InflowProfile::view_with (const amrex::Real* tz, const amrex::Real* tu, const amrex::Real* tv,
                                              const amrex::Real* tt, const amrex::Real* tk, const amrex::Real* rz,
                                              const amrex::Real* rt) const
{
    View p{};
    using K = InflowProfileOptions::Kind;
    p.kind = (m_opt.kind == K::log_law) ? 0 : (m_opt.kind == K::power_law) ? 1 : 2;
    p.ustar = m_opt.ustar;
    p.z0 = m_opt.z0;
    p.speed_ref = m_opt.speed_ref;
    p.z_ref = m_opt.z_ref;
    p.alpha = m_opt.alpha;
    const double a = static_cast<double>(m_opt.direction) * 3.14159265358979323846 / 180;
    p.cosd = static_cast<amrex::Real>(std::cos(a));
    p.sind = static_cast<amrex::Real>(std::sin(a));
    p.Cmu0 = m_opt.Cmu0;
    p.kappa = m_opt.kappa;
    p.tz = tz;
    p.tu = tu;
    p.tv = tv;
    p.tt = tt;
    p.tk = tk;
    p.n = static_cast<int>(m_hz.size());
    p.has_k = m_has_k;
    p.rz = rz;
    p.rt = rt;
    p.rn = static_cast<int>(m_hrz.size());
    return p;
}

} // namespace ablsteady
