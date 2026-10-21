#include "SteadyCriterion.H"

#include <AMReX.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace gust {

namespace {

std::string fmt (amrex::Real x)
{
    std::ostringstream os;
    os.precision(3);
    os << x;
    return os.str();
}

} // namespace

// Every check tests std::isfinite first: an ordered comparison (<, >) with a
// NaN raises FE_INVALID, which the GPU-CI flag amrex.fpe_trap_invalid=1 traps.
std::vector<std::string> validate (const SteadyParams& p)
{
    std::vector<std::string> err;
    if (!std::isfinite(p.residual_tol) || !(p.residual_tol > 0)) {
        err.emplace_back("steady.residual_tol must be finite and > 0");
    }
    // Above max_exponent10 orders, 10^-orders underflows to 0 and only an
    // exact zero residual would pass.
    if (!std::isfinite(p.residual_drop_orders) || !(p.residual_drop_orders >= 0)
        || !(p.residual_drop_orders <= std::numeric_limits<amrex::Real>::max_exponent10)) {
        err.emplace_back("steady.residual_drop_orders must be >= 0 and <= "
                         + std::to_string(std::numeric_limits<amrex::Real>::max_exponent10));
    }
    if (p.residual_ref_iterations < 1) {
        err.emplace_back("steady.residual_ref_iterations must be >= 1");
    }
    if (!std::isfinite(p.monitor_rel_tol) || !(p.monitor_rel_tol > 0)) {
        err.emplace_back("steady.monitor_rel_tol must be finite and > 0");
    }
    if (p.monitor_window < 2) {
        err.emplace_back("steady.monitor_window must be >= 2");
    }
    return err;
}

SteadyCriterion::SteadyCriterion (const SteadyParams& p, std::vector<std::string> residual_names,
                                  std::vector<MonitorSpec> monitors)
    : m_p(p), m_res_names(std::move(residual_names)), m_mon(std::move(monitors))
{
    auto err = validate(m_p);
    if (m_res_names.empty()) { err.emplace_back("SteadyCriterion needs at least one residual"); }
    for (const auto& m : m_mon) {
        if (!std::isfinite(m.scale) || !(m.scale > 0)) {
            err.emplace_back("monitor '" + m.name + "': scale must be finite and > 0");
        }
    }
    if (!err.empty()) {
        std::string all = "gust: invalid steady criterion:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    m_h.res_ref.assign(m_res_names.size(), amrex::Real(0));
    m_h.res_last.assign(m_res_names.size(), amrex::Real(0));
    m_h.mon_hist.assign(m_mon.size(), {});
}

void SteadyCriterion::restore (const SteadyHistory& h)
{
    const auto keep = static_cast<std::size_t>(m_p.monitor_window) + 1;
    bool ok = h.count >= 0 && h.res_ref.size() == m_res_names.size()
              && h.res_last.size() == m_res_names.size() && h.mon_hist.size() == m_mon.size();
    for (const auto& mh : h.mon_hist) { ok = ok && mh.size() <= keep; }
    if (!ok) {
        amrex::Abort("gust: SteadyCriterion::restore: the saved history does not match "
                     "this criterion's residuals, monitors or steady.monitor_window");
    }
    m_h = h;
}

void SteadyCriterion::seed_reference (const std::vector<amrex::Real>& ref)
{
    bool ok = m_h.count == 0 && ref.size() == m_res_names.size();
    for (const auto r : ref) { ok = ok && std::isfinite(r) && r >= 0; }
    if (!ok) {
        amrex::Abort("gust: SteadyCriterion::seed_reference: needs one finite reference >= 0 per residual, "
                     "before the first iteration");
    }
    m_h.res_ref = ref;
}

void SteadyCriterion::add (const std::vector<amrex::Real>& residuals,
                           const std::vector<amrex::Real>& monitors)
{
    if (residuals.size() != m_res_names.size() || monitors.size() != m_mon.size()) {
        amrex::Abort("gust: SteadyCriterion::add got " + std::to_string(residuals.size())
                     + " residuals and " + std::to_string(monitors.size()) + " monitors, expected "
                     + std::to_string(m_res_names.size()) + " and " + std::to_string(m_mon.size()));
    }

    for (std::size_t q = 0; q < residuals.size(); ++q) {
        const amrex::Real r = residuals[q];
        if (!std::isfinite(r)) {
            if (!m_h.nonfinite) { m_h.nonfinite_what = "residual '" + m_res_names[q] + "'"; }
            m_h.nonfinite = true;
            continue;
        }
        if (m_h.count < m_p.residual_ref_iterations) {
            m_h.res_ref[q] = std::max(m_h.res_ref[q], std::abs(r));
        }
        m_h.res_last[q] = std::abs(r);
    }

    const auto keep = static_cast<std::size_t>(m_p.monitor_window) + 1;
    for (std::size_t m = 0; m < monitors.size(); ++m) {
        const amrex::Real v = monitors[m];
        if (!std::isfinite(v)) {
            if (!m_h.nonfinite) { m_h.nonfinite_what = "monitor '" + m_mon[m].name + "'"; }
            m_h.nonfinite = true;
            continue;
        }
        auto& h = m_h.mon_hist[m];
        h.push_back(v);
        if (h.size() > keep) { h.erase(h.begin()); }
    }
    ++m_h.count;
}

SteadyVerdict SteadyCriterion::verdict () const
{
    SteadyVerdict v;
    v.iteration = m_h.count - 1;
    if (m_h.count == 0) {
        v.reasons.emplace_back("no iterations recorded");
        return v;
    }
    if (m_h.nonfinite) {
        v.status = SteadyStatus::Diverged;
        v.reasons.emplace_back(m_h.nonfinite_what + " became NaN or infinite");
        return v;
    }

    // Residuals: below the absolute tolerance AND dropped by the required
    // orders below the reference (the largest of the first few values).
    const amrex::Real drop = std::pow(amrex::Real(10), -m_p.residual_drop_orders);
    for (std::size_t q = 0; q < m_res_names.size(); ++q) {
        const amrex::Real r = m_h.res_last[q];
        // Floored at the tolerance: a residual that starts at or near zero
        // is not held to a drop below round-off.
        const amrex::Real ref = std::max(m_h.res_ref[q], m_p.residual_tol);
        if (!(r < m_p.residual_tol)) {
            v.reasons.emplace_back("residual '" + m_res_names[q] + "' = " + fmt(r)
                                   + " is not below steady.residual_tol = " + fmt(m_p.residual_tol));
        }
        if (!(r <= ref * drop)) {
            const amrex::Real orders = (r > 0) ? std::log10(ref / r) : amrex::Real(0);
            v.reasons.emplace_back("residual '" + m_res_names[q] + "' dropped " + fmt(orders)
                                   + " orders, less than steady.residual_drop_orders = "
                                   + fmt(m_p.residual_drop_orders));
        }
    }

    // Monitors: a full window, and every value in it within the tolerance of
    // the latest value.
    const auto need = static_cast<std::size_t>(m_p.monitor_window) + 1;
    for (std::size_t m = 0; m < m_mon.size(); ++m) {
        const auto& h = m_h.mon_hist[m];
        if (h.size() < need) {
            v.reasons.emplace_back("monitor '" + m_mon[m].name + "': " + std::to_string(h.size())
                                   + " of the " + std::to_string(need)
                                   + " values steady.monitor_window needs");
            continue;
        }
        const amrex::Real last = h.back();
        amrex::Real change = 0;
        for (const auto x : h) { change = std::max(change, std::abs(x - last)); }
        const amrex::Real ref = std::max(std::abs(last), m_mon[m].scale);
        const amrex::Real tol = m_p.monitor_rel_tol * ref;
        if (!(change <= tol)) {
            v.reasons.emplace_back("monitor '" + m_mon[m].name + "' changed by " + fmt(change / ref)
                                   + " of max(|value|, scale) over the last "
                                   + std::to_string(m_p.monitor_window)
                                   + " iterations; steady.monitor_rel_tol = " + fmt(m_p.monitor_rel_tol));
            continue;
        }
        // A small change per window does not bound the distance to the
        // limit when the monitor converges slowly. Aitken's estimate from
        // the changes over the window's two halves, d1 and d2 = rho d1:
        // remaining = d2 rho / (1 - rho), valid for -1 < rho < 1.
        const std::size_t W = h.size() - 1;
        const std::size_t half = W / 2;
        const amrex::Real d1 = h[W - half] - h[W - 2 * half];
        const amrex::Real d2 = h[W] - h[W - half];
        // Changes far below the tolerance (round-off) say nothing about rho.
        const amrex::Real noise = amrex::Real(1.0e-3) * tol;
        if (std::max(std::abs(d1), std::abs(d2)) <= noise) { continue; }
        const amrex::Real rho = (d1 != 0) ? d2 / d1 : std::numeric_limits<amrex::Real>::infinity();
        if (!(rho > -1 && rho < 1)) {
            v.reasons.emplace_back("monitor '" + m_mon[m].name
                                   + "' is not converging over the window (half-window change ratio "
                                   + fmt(rho) + ")");
            continue;
        }
        const amrex::Real remaining = std::abs(d2 * rho / (1 - rho));
        if (!(remaining <= tol)) {
            v.reasons.emplace_back("monitor '" + m_mon[m].name + "' is an extrapolated "
                                   + fmt(remaining / ref)
                                   + " of max(|value|, scale) from its limit; steady.monitor_rel_tol = "
                                   + fmt(m_p.monitor_rel_tol));
        }
    }

    v.status = v.reasons.empty() ? SteadyStatus::Steady : SteadyStatus::NotSteady;
    return v;
}

} // namespace gust
