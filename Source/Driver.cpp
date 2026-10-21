#include "Driver.H"

#include "PseudoTime.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

namespace ablsteady {

namespace {

const char* const kKeys[] = {"residual_tol", "residual_drop_orders", "residual_ref_iterations", "monitor_rel_tol",
                             "monitor_window", "max_iterations", "cfl0", "cfl_max", "ser_exponent", "ser_max_growth",
                             "switch_window", "switch_residual", "switch_back_growth", "zline_sweeps", "sweep", "line_ghosts", "line_overlap",
                             "line_order", "line_alternate", "machine", "nodes", "stage",
                             "newton_preconditioner", "switch_monitor_tol", "newton_fallback_failures",
                             "newton_cfl_max"};

// 2: the reference hold; 3: the residual scale height (versions 1 and 2 still read)
constexpr int kDriverVersion = 3;
const char* const kDriverMagic = "ablsteady-driver";

std::string hexfloat (double x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", x);
    return buf;
}

std::string fmt (double x)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", x);
    return buf;
}

void history_error (const std::string& dir, const std::string& what)   // aborts
{
    amrex::Abort("ablsteady: " + dir + "/Driver: " + what);
}

} // namespace

// ---------------------------------------------------------------------------
std::vector<std::string> check_steady_options (const SteadyOptions& opt)
{
    std::vector<std::string> err = validate(opt.criterion);
    auto bad = [&] (amrex::Real x, const char* key, const char* range, const auto& ok) {
        // isnan first: gcc traps on an ordered compare with a NaN under the invalid trap
        if (std::isnan(x) || !ok()) {
            err.emplace_back(std::string("steady.") + key + " = " + std::to_string(double(x)) + ": must be " + range);
        }
    };
    if (opt.max_iterations < 1) { err.emplace_back("steady.max_iterations must be >= 1"); }
    bad(opt.cfl0, "cfl0", "finite and > 0", [&] { return std::isfinite(opt.cfl0) && opt.cfl0 > 0; });
    bad(opt.cfl_max, "cfl_max", "finite and >= steady.cfl0",
        [&] { return std::isfinite(opt.cfl_max) && !std::isnan(opt.cfl0) && opt.cfl_max >= opt.cfl0; });
    // Not given: the default of the sweep (default_cfl_max); say so, since the
    // message names a key the deck does not have.
    if (!err.empty() && err.back().rfind("steady.cfl_max = ", 0) == 0
        && std::find(opt.given.begin(), opt.given.end(), "cfl_max") == opt.given.end()) {
        err.back() += std::string(" (the default with steady.sweep = ") + sweep_name(opt.sweep) + "; give steady.cfl_max)";
    }
    bad(opt.ser_exponent, "ser_exponent", "finite and >= 0", [&] { return std::isfinite(opt.ser_exponent) && opt.ser_exponent >= 0; });
    bad(opt.ser_max_growth, "ser_max_growth", "finite and >= 1",
        [&] { return std::isfinite(opt.ser_max_growth) && opt.ser_max_growth >= 1; });
    if (opt.switch_window < 1) { err.emplace_back("steady.switch_window must be >= 1"); }
    bad(opt.switch_residual, "switch_residual", "> 0 and <= 1",
        [&] { return opt.switch_residual > 0 && opt.switch_residual <= 1; });
    bad(opt.switch_back_growth, "switch_back_growth", "finite and > 1",
        [&] { return std::isfinite(opt.switch_back_growth) && opt.switch_back_growth > 1; });
    bad(opt.switch_monitor_tol, "switch_monitor_tol", "0 (off), or > 0 and < 1",
        [&] { return opt.switch_monitor_tol >= 0 && opt.switch_monitor_tol < 1; });
    if (opt.newton_fallback_failures < 1) { err.emplace_back("steady.newton_fallback_failures must be >= 1"); }
    // The Newton stage's cap (only that stage reads it; with one stage it aborts as a key with no effect, App.cpp).
    if (opt.stage == StageKind::simplec_newton) {
        bad(opt.newton_cfl_max, "newton_cfl_max", "finite and >= steady.cfl_max",
            [&] { return std::isfinite(opt.newton_cfl_max) && !std::isnan(opt.cfl_max) && opt.newton_cfl_max >= opt.cfl_max; });
    }
    if (opt.nodes < 1) { err.emplace_back("steady.nodes must be >= 1"); }
    if (opt.sweep == SweepKind::lusgs) {
        if (opt.zline_sweeps < 1 || opt.zline_sweeps > 100) {
            err.emplace_back("steady.zline_sweeps must be >= 1 and <= 100 with steady.sweep = lusgs");
        }
    } else if (opt.zline_sweeps < 2 || opt.zline_sweeps > 100) {
        err.emplace_back("steady.zline_sweeps must be >= 2 and <= 100 (one lagged sweep is unstable at large horizontal CFL)");
    }
    for (auto& e : check_line_options(line_sweep_options(opt))) { err.push_back(std::move(e)); }
    // Only the odd sweeps are reversed, so with one sweep the option would be
    // read and ignored.
    if (opt.sweep == SweepKind::lusgs && opt.line_alternate == 1 && opt.zline_sweeps < 2) {
        err.emplace_back("steady.line_alternate = 1 needs steady.zline_sweeps >= 2 (only the odd sweeps are reversed)");
    }
    if (opt.sweep != SweepKind::lusgs) {
        for (const char* key : {"line_ghosts", "line_overlap", "line_order", "line_alternate"}) {
            if (std::find(opt.given.begin(), opt.given.end(), key) != opt.given.end()) {
                err.emplace_back(std::string("steady.") + key + " has no effect unless steady.sweep = lusgs");
            }
        }
    }
    return err;
}

SweepKind default_sweep () noexcept
{
#ifdef AMREX_USE_GPU
    return SweepKind::zjacobi;
#else
    return SweepKind::lusgs;
#endif
}

amrex::Real default_cfl_max (SweepKind s) noexcept
{
    return s == SweepKind::lusgs ? amrex::Real(8) : amrex::Real(10);
}

LineSweepOptions line_sweep_options (const SteadyOptions& opt)
{
    LineSweepOptions o(opt.sweep);
    o.ghosts = opt.line_ghosts;
    o.overlap = opt.line_overlap;
    o.order = opt.line_order;
    o.alternate = opt.line_alternate;
    return o;
}

SteadyOptions read_steady_options ()
{
    SteadyOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("steady");
    for (const char* key : kKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("steady.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { err.emplace_back(std::string("steady.") + key + " takes exactly one value"); }
    }
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    pp.query("residual_tol", opt.criterion.residual_tol);
    pp.query("residual_drop_orders", opt.criterion.residual_drop_orders);
    pp.query("residual_ref_iterations", opt.criterion.residual_ref_iterations);
    pp.query("monitor_rel_tol", opt.criterion.monitor_rel_tol);
    pp.query("monitor_window", opt.criterion.monitor_window);
    pp.query("max_iterations", opt.max_iterations);
    pp.query("cfl0", opt.cfl0);
    pp.query("cfl_max", opt.cfl_max);
    pp.query("ser_exponent", opt.ser_exponent);
    pp.query("ser_max_growth", opt.ser_max_growth);
    pp.query("switch_window", opt.switch_window);
    pp.query("switch_residual", opt.switch_residual);
    pp.query("switch_back_growth", opt.switch_back_growth);
    pp.query("switch_monitor_tol", opt.switch_monitor_tol);
    pp.query("newton_fallback_failures", opt.newton_fallback_failures);
    pp.query("newton_cfl_max", opt.newton_cfl_max);
    pp.query("nodes", opt.nodes);
    pp.query("zline_sweeps", opt.zline_sweeps);
    std::string sweep;
    opt.sweep = default_sweep();
    if (pp.query("sweep", sweep)) {
        if (sweep == "zjacobi") { opt.sweep = SweepKind::zjacobi; }
        else if (sweep == "lusgs") { opt.sweep = SweepKind::lusgs; }
        else { amrex::Abort("ablsteady: steady.sweep = " + sweep + ": must be zjacobi or lusgs"); }
    }
    if (!pp.contains("cfl_max")) { opt.cfl_max = default_cfl_max(opt.sweep); }
    std::string pc;
    if (pp.query("newton_preconditioner", pc)) {
        if (pc == "identity") { opt.newton_preconditioner = NewtonPreconditionerKind::identity; }
        else if (pc == "simplec") { opt.newton_preconditioner = NewtonPreconditionerKind::simplec; }
        else { amrex::Abort("ablsteady: steady.newton_preconditioner = " + pc + ": must be identity or simplec"); }
    }
    std::string word;
    if (pp.query("line_ghosts", word)) {
        if (word == "lagged") { opt.line_ghosts = LineGhosts::lagged; }
        else if (word == "zero") { opt.line_ghosts = LineGhosts::zero; }
        else { amrex::Abort("ablsteady: steady.line_ghosts = " + word + ": must be lagged or zero"); }
    }
    pp.query("line_overlap", opt.line_overlap);
    if (pp.query("line_order", word)) {
        if (word == "sequential") { opt.line_order = LineOrder::sequential; }
        else if (word == "diagonal") { opt.line_order = LineOrder::diagonal; }
        else { amrex::Abort("ablsteady: steady.line_order = " + word + ": must be sequential or diagonal"); }
    }
    pp.query("line_alternate", opt.line_alternate);
    std::string stage;
    if (pp.query("stage", stage)) {
        if (stage == "none") { opt.stage = StageKind::none; }
        else if (stage == "simplec") { opt.stage = StageKind::simplec; }
        else if (stage == "simplec_newton") { opt.stage = StageKind::simplec_newton; }
        else { amrex::Abort("ablsteady: steady.stage = " + stage + ": must be none, simplec or simplec_newton"); }
    }
    std::string machine;
    if (pp.query("machine", machine) && !parse_machine(machine, opt.machine)) {
        amrex::Abort("ablsteady: steady.machine = " + machine + ": must be local, kestrel-cpu or kestrel-gpu");
    }
    err = check_steady_options(opt);
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

// ---------------------------------------------------------------------------
namespace {

// The first stage, after the checks that must pass before the criterion can
// be built from it (so that their messages, not the criterion's, name the
// fault). Every constructor argument that touches the stages goes through
// this call: the order in which arguments are evaluated is unspecified (gcc
// evaluates right to left), so no argument may dereference a stage on its own.
const Stage& first_stage (const SteadyOptions& opt, const std::vector<std::unique_ptr<Stage>>& stages)
{
    for (const auto& e : check_steady_options(opt)) { amrex::Abort("ablsteady: Driver: " + e); }
    if (stages.empty()) { amrex::Abort("ablsteady: Driver: at least one stage is needed"); }
    for (const auto& s : stages) {
        if (s == nullptr) { amrex::Abort("ablsteady: Driver: a stage is null"); }
    }
    return *stages.front();
}

} // namespace

Driver::Driver (const SteadyOptions& opt, std::vector<std::unique_ptr<Stage>> stages, const Mesh& mesh,
                const amrex::Geometry& geom, const amrex::MultiFab* nu)
    : m_opt(opt),
      m_stages(std::move(stages)),
      m_mesh(mesh),
      m_geom(geom),
      m_nu(nu),
      m_criterion(opt.criterion, first_stage(opt, m_stages).residual_names(), first_stage(opt, m_stages).monitor_specs())
{
    const auto& names = m_stages.front()->residual_names();
    const auto& specs = m_stages.front()->monitor_specs();
    for (const auto& s : m_stages) {
        bool same = s->residual_names() == names && s->monitor_specs().size() == specs.size();
        for (std::size_t m = 0; same && m < specs.size(); ++m) {
            same = s->monitor_specs()[m].name == specs[m].name && s->monitor_specs()[m].scale == specs[m].scale;
        }
        if (!same) {
            amrex::Abort("ablsteady: Driver: stage " + s->name() + " does not report the residuals and monitors of stage "
                         + m_stages.front()->name());
        }
    }
    m_h.accounts.resize(m_stages.size());
    m_h.cfl = opt.cfl0;
    m_h.scale_height = residual_scale_height(m_mesh, m_geom);
}

void Driver::seed_reference (const std::vector<amrex::Real>& ref)
{
    if (!m_h.records.empty() || m_h.criterion.count != 0) {
        amrex::Abort("ablsteady: Driver::seed_reference: only before the first iteration");
    }
    m_criterion.seed_reference(ref);
    m_h.criterion = m_criterion.history();
}

void Driver::release_reference_hold ()
{
    if (!m_h.records.empty() || m_h.criterion.count != 0) {
        amrex::Abort("ablsteady: Driver::release_reference_hold: only before the first iteration");
    }
    bool seeded = false;
    for (const amrex::Real x : m_h.criterion.res_ref) { seeded = seeded || x > 0; }
    if (!seeded) {
        // Unseeded, the reference is still 0 and r_n would be measured
        // against residual_tol alone.
        amrex::Abort("ablsteady: Driver::release_reference_hold: the reference is not seeded (seed_reference first)");
    }
    m_h.hold = 0;
}

void Driver::start_cfl (amrex::Real cfl)
{
    if (!m_h.records.empty() || m_h.criterion.count != 0) {
        amrex::Abort("ablsteady: Driver::start_cfl: only before the first iteration");
    }
    if (hold_iterations() > 0) { amrex::Abort("ablsteady: Driver::start_cfl: the reference hold is on"); }
    if (std::isnan(cfl) || !std::isfinite(cfl) || !(cfl > 0) || cfl > m_opt.cfl_max) {
        amrex::Abort("ablsteady: Driver::start_cfl: " + fmt(double(cfl)) + " is not finite, > 0 and <= steady.cfl_max = "
                     + fmt(double(m_opt.cfl_max)));
    }
    // The SER law goes below cfl0 when r_n > 1 (a residual above its
    // reference), so a coarse level can end below it; as restore does, the
    // start is brought up to cfl0.
    m_h.cfl = amrex::max(cfl, m_opt.cfl0);
}

int Driver::hold_iterations () const
{
    return m_h.hold < 0 ? m_opt.criterion.residual_ref_iterations : m_h.hold;
}

amrex::Real Driver::relative_residual () const
{
    return m_h.records.empty() ? amrex::Real(0) : m_h.records.back().relative_residual;
}

amrex::Real Driver::next_cfl (amrex::Real r) const
{
    // The robust start: cfl0 while the criterion's reference is still being
    // taken (the first residual_ref_iterations, the hold; none on a level
    // whose reference is seeded and released), then the SER law capped by
    // cfl_max and by the growth per iteration (ser_cfl, PseudoTime.H).
    if (m_criterion.iterations() <= hold_iterations()) { return m_opt.cfl0; }
    return ser_cfl(m_h.cfl, r, m_opt.cfl0, m_opt.cfl_max, m_opt.ser_exponent, m_opt.ser_max_growth);
}

void Driver::apply_switch_rules (int iteration, amrex::Real r, const std::vector<amrex::Real>& residuals, State& q)
{
    const int W = m_opt.switch_window;
    m_h.below_count = (r <= m_opt.switch_residual) ? m_h.below_count + 1 : 0;
    int to = m_h.stage;
    std::string reason;
    if (m_h.stage > 0 && m_stages[m_h.stage]->fallback_requested()) {
        // The stage's own request (the Newton stage's failed steps, or its
        // hand-back for the final pass): at once, whatever the hysteresis.
        to = m_h.stage - 1;
        reason = m_stages[m_h.stage]->fallback_reason();
    } else if (m_h.last_switch >= 0 && iteration - m_h.last_switch < W) {
        return;   // hysteresis
    } else if (m_h.below_count >= W && m_h.stage + 1 < static_cast<int>(m_stages.size())
               && m_stages[m_h.stage + 1]->accepts_switch() && monitors_settled()) {
        to = m_h.stage + 1;
        reason = "relative residual " + fmt(double(r)) + " <= steady.switch_residual " + fmt(double(m_opt.switch_residual))
                 + " for " + std::to_string(W) + " iterations";
        if (m_opt.switch_monitor_tol > 0) {
            reason += ", monitors within steady.switch_monitor_tol " + fmt(double(m_opt.switch_monitor_tol));
        }
    } else if (m_h.stage > 0 && static_cast<int>(m_h.records.size()) > W && iteration - W > m_h.last_switch) {
        // the record W iterations back belongs to the current stage
        const amrex::Real before = m_h.records[m_h.records.size() - 1 - W].relative_residual;
        if (r >= m_opt.switch_back_growth * before) {
            to = m_h.stage - 1;
            reason = "relative residual grew from " + fmt(double(before)) + " to " + fmt(double(r)) + " over "
                     + std::to_string(W) + " iterations (steady.switch_back_growth " + fmt(double(m_opt.switch_back_growth)) + ")";
        }
    }
    if (to == m_h.stage) { return; }
    SwitchEvent e;
    e.iteration = iteration;
    e.from = m_h.stage;
    e.to = to;
    e.reason = reason;
    e.residuals = residuals;
    m_h.switches.push_back(e);
    std::string state;
    const auto& names = m_stages.front()->residual_names();
    for (std::size_t i = 0; i < residuals.size(); ++i) {
        state += (i ? ", " : "") + names[i] + " " + fmt(double(residuals[i]));
    }
    amrex::Print() << "ablsteady: iteration " << iteration << ": switching from stage " << m_stages[e.from]->name()
                   << " to " << m_stages[to]->name() << ": " << reason << " (scaled residuals: " << state << ")\n";
    m_stages[e.from]->leave(q);
    m_h.stage = to;
    m_h.below_count = 0;
    m_h.last_switch = iteration;
    m_stages[to]->init(q);
}

bool Driver::monitors_settled () const
{
    if (m_opt.switch_monitor_tol <= 0) { return true; }
    const int W = m_opt.switch_window;
    const auto& rec = m_h.records;
    if (static_cast<int>(rec.size()) < W + 1) { return false; }
    const auto& specs = m_stages.front()->monitor_specs();
    const auto& now = rec.back().monitors;
    const auto& then = rec[rec.size() - 1 - static_cast<std::size_t>(W)].monitors;
    for (std::size_t m = 0; m < now.size(); ++m) {
        const amrex::Real scale = amrex::max(std::abs(now[m]), specs[m].scale);
        if (!(std::abs(now[m] - then[m]) <= m_opt.switch_monitor_tol * scale)) { return false; }
    }
    return true;
}

SteadyVerdict Driver::step (State& q)
{
    Stage& stage = *m_stages[m_h.stage];
    if (m_need_init) {
        stage.init(q);
        m_need_init = false;
    }
    const double t0 = amrex::second();
    local_dtau(q, m_nu, m_mesh, m_geom, m_h.cfl, m_dtau);
    stage.iterate(q, m_dtau);
    const std::vector<amrex::Real> res = stage.residual_norms();
    const std::vector<amrex::Real> mon = stage.monitors();
    m_criterion.add(res, mon);
    // The iteration's wall time is the slowest rank's, so every rank holds
    // the same history (the checkpoint is rank 0's, read by all).
    double wall = amrex::second() - t0;
    amrex::ParallelDescriptor::ReduceRealMax(wall);

    // r_n: every residual over its reference (the criterion's, floored at
    // the tolerance as the criterion does); 1 during the reference hold.
    // Without the hold (a seeded level) the first residual_ref_iterations
    // values still raise the reference, so r_n <= 1 while they do.
    amrex::Real r = 1;
    const SteadyHistory& ch = m_criterion.history();
    if (ch.nonfinite) {
        r = relative_residual();   // Diverged: no arithmetic on the NaN or inf (the traps)
    } else if (m_criterion.iterations() > hold_iterations()) {
        r = 0;
        for (std::size_t i = 0; i < res.size(); ++i) {
            const amrex::Real ref = amrex::max(ch.res_ref[i], m_opt.criterion.residual_tol);
            r = amrex::max(r, std::abs(res[i]) / ref);
        }
    }
    IterationRecord rec;
    rec.iteration = static_cast<int>(m_h.records.size());
    rec.stage = m_h.stage;
    rec.cfl = m_h.cfl;
    rec.relative_residual = r;
    rec.wall_seconds = wall;
    rec.residuals = res;
    rec.monitors = mon;
    m_h.records.push_back(rec);
    StageAccount& acc = m_h.accounts[m_h.stage];
    acc.iterations += 1;
    acc.wall_seconds += wall;
    acc.au += allocation_units(m_opt.machine, m_opt.nodes, wall);
    m_h.criterion = ch;

    const SteadyVerdict v = m_criterion.verdict();
    if (v.status == SteadyStatus::NotSteady) {
        m_h.cfl = next_cfl(r);
        apply_switch_rules(rec.iteration, r, res, q);
    }
    return v;
}

DriverResult Driver::run (State& q)
{
    DriverResult out;
    out.verdict = m_criterion.verdict();
    while (static_cast<int>(m_h.records.size()) < m_opt.max_iterations) {
        out.verdict = step(q);
        ++out.iterations;
        if (out.verdict.status != SteadyStatus::NotSteady) { break; }
    }
    return out;
}

void Driver::restore (const DriverHistory& h_in, State& q)
{
    // A run with more stages may continue a history with fewer (a
    // steady.stage = simplec checkpoint finished with simplec_newton): the
    // missing stages start with empty accounts. The reverse loses a stage.
    DriverHistory h_pad;
    const bool pad = h_in.accounts.size() < m_stages.size() && h_in.stage < static_cast<int>(h_in.accounts.size());
    if (pad) {
        h_pad = h_in;
        h_pad.accounts.resize(m_stages.size());
    }
    const DriverHistory& h = pad ? h_pad : h_in;
    if (h.accounts.size() != m_stages.size()) {
        amrex::Abort("ablsteady: Driver::restore: the history has " + std::to_string(h.accounts.size()) + " stages, this driver "
                     + std::to_string(m_stages.size()) + " (a steady.stage = simplec_newton checkpoint cannot continue with "
                     "steady.stage = simplec)");
    }
    if (h.stage < 0 || h.stage >= static_cast<int>(m_stages.size())) {
        amrex::Abort("ablsteady: Driver::restore: the history's stage " + std::to_string(h.stage) + " does not exist");
    }
    for (const auto& r : h.records) {
        if (r.residuals.size() != m_stages.front()->residual_names().size()
            || r.monitors.size() != m_stages.front()->monitor_specs().size()) {
            amrex::Abort("ablsteady: Driver::restore: a record's residual or monitor count does not match the stages");
        }
    }
    if (std::isnan(h.cfl) || !std::isfinite(h.cfl) || !(h.cfl > 0)) {
        amrex::Abort("ablsteady: Driver::restore: the history's CFL is not finite and > 0");
    }
    if (h.below_count < 0 || h.last_switch < -1 || h.last_switch >= static_cast<int>(h.records.size())
        || h.criterion.count != static_cast<int>(h.records.size())) {
        amrex::Abort("ablsteady: Driver::restore: the history's switch state or iteration count is inconsistent");
    }
    if (h.hold != -1 && h.hold != 0) {
        amrex::Abort("ablsteady: Driver::restore: the history's reference hold " + std::to_string(h.hold) + " is not -1 or 0");
    }
    if (std::isnan(h.scale_height) || !std::isfinite(h.scale_height) || h.scale_height < 0) {
        amrex::Abort("ablsteady: Driver::restore: the history's residual scale height is not finite and >= 0");
    }
    // The residual norms scale with H (rms / (U^2/H), ...). A history written
    // with another H (format 1 and 2: prob_hi z - prob_lo z, which on a
    // stretched mesh is not the depth of the layers) is brought to this H, so
    // the reference, the drop and steady.residual_tol continue in one scale.
    // r_n, a ratio of two norms, is kept as recorded.
    DriverHistory hs = h;
    const amrex::Real H_old = (h.scale_height > 0) ? h.scale_height : m_geom.ProbHi(2) - m_geom.ProbLo(2);
    const amrex::Real H_new = residual_scale_height(m_mesh, m_geom);
    if (H_old != H_new) {
        const amrex::Real f = H_new / H_old;
        auto rescale = [f] (std::vector<amrex::Real>& v) { for (auto& x : v) { x *= f; } };
        rescale(hs.criterion.res_ref);
        rescale(hs.criterion.res_last);
        for (auto& r : hs.records) { rescale(r.residuals); }
        for (auto& e : hs.switches) { rescale(e.residuals); }
        amrex::Print() << "ablsteady: restart: the history's residuals were scaled with H = " << fmt(double(H_old))
                       << " m; rescaled by " << fmt(double(f)) << " to H = " << fmt(double(H_new)) << " m\n";
    }
    hs.scale_height = H_new;
    m_criterion.restore(hs.criterion);   // checks its own sizes
    m_h = hs;
    // The options may have changed across the restart: the CFL continues
    // within the new bounds.
    m_h.cfl = amrex::min(amrex::max(m_h.cfl, m_opt.cfl0), m_opt.cfl_max);
    m_need_init = true;
    (void)q;
}

// ---------------------------------------------------------------------------
namespace {

void write_reals (std::ostream& os, const std::vector<amrex::Real>& v)
{
    os << v.size();
    for (const amrex::Real x : v) { os << " " << hexfloat(static_cast<double>(x)); }
}

// Hex floats are read through strtod: libstdc++'s stream extraction does not
// accept them (StateIO.cpp and Reference.cpp do the same).
double read_double (std::istream& is, const std::string& dir, const char* what)
{
    std::string tok;
    is >> tok;
    char* end = nullptr;
    const double x = tok.empty() ? 0.0 : std::strtod(tok.c_str(), &end);
    if (!is || tok.empty() || *end != '\0') { history_error(dir, std::string("corrupt: bad value for ") + what); }
    return x;
}

int read_int (std::istream& is, const std::string& dir, const char* what)
{
    int i = 0;
    if (!(is >> i)) { history_error(dir, std::string("corrupt: bad value for ") + what); }
    return i;
}

//! A count of items still to be read: non-negative and no more than the
//! buffer could hold (each item takes at least two bytes), so a corrupt
//! count cannot allocate gigabytes before the end of the file is noticed.
int read_count (std::istream& is, const std::string& dir, const char* what, std::size_t bytes)
{
    const int n = read_int(is, dir, what);
    if (n < 0 || static_cast<std::size_t>(n) > bytes / 2) { history_error(dir, std::string("corrupt: count out of range for ") + what); }
    return n;
}

std::vector<amrex::Real> read_reals (std::istream& is, const std::string& dir, const char* what, std::size_t bytes)
{
    const int n = read_count(is, dir, what, bytes);
    std::vector<amrex::Real> v(static_cast<std::size_t>(n));
    for (auto& x : v) { x = static_cast<amrex::Real>(read_double(is, dir, what)); }
    return v;
}

void expect (std::istream& is, const std::string& dir, const char* word)
{
    std::string got;
    if (!(is >> got) || got != word) { history_error(dir, std::string("corrupt: expected '") + word + "'"); }
}

// A string that may hold spaces: its length, then the bytes after one space.
void write_string (std::ostream& os, const std::string& s) { os << s.size() << " " << s; }

std::string read_string (std::istream& is, const std::string& dir, const char* what, std::size_t bytes)
{
    const int n = read_int(is, dir, what);
    if (n < 0 || static_cast<std::size_t>(n) > bytes) { history_error(dir, std::string("corrupt: length out of range for ") + what); }
    is.get();   // the separating space
    std::string s(static_cast<std::size_t>(n), '\0');
    if (n > 0 && !is.read(&s[0], static_cast<std::streamsize>(n))) {
        history_error(dir, std::string("corrupt: short string for ") + what);
    }
    return s;
}

} // namespace

void write_driver_history (const std::string& dir, const DriverHistory& h)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        amrex::ParallelDescriptor::Barrier();
        return;
    }
    const std::string name = dir + "/Driver";
    std::ofstream os(name);
    if (!os.good()) { amrex::FileOpenFailed(name); }
    os << kDriverMagic << " " << kDriverVersion << "\n";
    os << "state " << h.stage << " " << hexfloat(static_cast<double>(h.cfl)) << " " << h.below_count << " " << h.last_switch << "\n";
    os << "hold " << h.hold << "\n";
    os << "scale_height " << hexfloat(static_cast<double>(h.scale_height)) << "\n";
    os << "accounts " << h.accounts.size() << "\n";
    for (const auto& a : h.accounts) {
        os << a.iterations << " " << hexfloat(a.wall_seconds) << " " << hexfloat(a.au) << "\n";
    }
    os << "records " << h.records.size() << "\n";
    for (const auto& r : h.records) {
        os << r.iteration << " " << r.stage << " " << hexfloat(static_cast<double>(r.cfl)) << " "
           << hexfloat(static_cast<double>(r.relative_residual)) << " " << hexfloat(r.wall_seconds) << " ";
        write_reals(os, r.residuals);
        os << " ";
        write_reals(os, r.monitors);
        os << "\n";
    }
    os << "switches " << h.switches.size() << "\n";
    for (const auto& s : h.switches) {
        os << s.iteration << " " << s.from << " " << s.to << " ";
        write_reals(os, s.residuals);
        os << " ";
        write_string(os, s.reason);
        os << "\n";
    }
    const SteadyHistory& c = h.criterion;
    os << "criterion " << c.count << " " << (c.nonfinite ? 1 : 0) << " ";
    write_string(os, c.nonfinite_what);
    os << "\n";
    write_reals(os, c.res_ref);
    os << "\n";
    write_reals(os, c.res_last);
    os << "\n";
    os << c.mon_hist.size() << "\n";
    for (const auto& m : c.mon_hist) {
        write_reals(os, m);
        os << "\n";
    }
    if (!os.good()) { history_error(dir, "could not be written"); }
    os.close();
    amrex::ParallelDescriptor::Barrier();
}

DriverHistory read_driver_history (const std::string& dir)
{
    // Every rank must learn of a missing file before the broadcast, or the
    // other ranks wait in it forever while rank 0 aborts.
    int exists = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(dir + "/Driver") ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&exists, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (exists == 0) { history_error(dir, "does not exist"); }
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(dir + "/Driver", buf);
    buf.push_back('\0');
    const std::size_t bytes = buf.size();
    std::istringstream is(buf.data());
    std::string magic;
    int version = 0;
    is >> magic >> version;
    if (!is || magic != kDriverMagic || version < 1 || version > kDriverVersion) {
        history_error(dir, "not an ablsteady driver history of format version 1 to " + std::to_string(kDriverVersion));
    }
    DriverHistory h;
    expect(is, dir, "state");
    h.stage = read_int(is, dir, "stage");
    h.cfl = static_cast<amrex::Real>(read_double(is, dir, "cfl"));
    h.below_count = read_int(is, dir, "below_count");
    h.last_switch = read_int(is, dir, "last_switch");
    if (version >= 2) {
        expect(is, dir, "hold");
        h.hold = read_int(is, dir, "hold");
        if (h.hold != -1 && h.hold != 0) { history_error(dir, "corrupt: hold must be -1 or 0"); }
    }   // version 1 had no line: the hold of steady.residual_ref_iterations
    if (version >= 3) {
        expect(is, dir, "scale_height");
        h.scale_height = static_cast<amrex::Real>(read_double(is, dir, "scale_height"));
        if (!std::isfinite(h.scale_height) || !(h.scale_height > 0)) { history_error(dir, "corrupt: scale_height must be > 0"); }
    }   // versions 1 and 2: 0, prob_hi z - prob_lo z (Driver::restore)
    expect(is, dir, "accounts");
    const int na = read_count(is, dir, "the account count", bytes);
    if (na < 1) { history_error(dir, "corrupt: no accounts"); }
    h.accounts.resize(na);
    for (auto& a : h.accounts) {
        a.iterations = read_int(is, dir, "account iterations");
        a.wall_seconds = read_double(is, dir, "account wall_seconds");
        a.au = read_double(is, dir, "account au");
    }
    expect(is, dir, "records");
    const int nr = read_count(is, dir, "the record count", bytes);
    h.records.resize(nr);
    for (auto& r : h.records) {
        r.iteration = read_int(is, dir, "record iteration");
        r.stage = read_int(is, dir, "record stage");
        r.cfl = static_cast<amrex::Real>(read_double(is, dir, "record cfl"));
        r.relative_residual = static_cast<amrex::Real>(read_double(is, dir, "record relative_residual"));
        r.wall_seconds = read_double(is, dir, "record wall_seconds");
        r.residuals = read_reals(is, dir, "record residuals", bytes);
        r.monitors = read_reals(is, dir, "record monitors", bytes);
    }
    expect(is, dir, "switches");
    const int ns = read_count(is, dir, "the switch count", bytes);
    h.switches.resize(ns);
    for (auto& s : h.switches) {
        s.iteration = read_int(is, dir, "switch iteration");
        s.from = read_int(is, dir, "switch from");
        s.to = read_int(is, dir, "switch to");
        s.residuals = read_reals(is, dir, "switch residuals", bytes);
        s.reason = read_string(is, dir, "switch reason", bytes);
    }
    expect(is, dir, "criterion");
    h.criterion.count = read_int(is, dir, "criterion count");
    h.criterion.nonfinite = read_int(is, dir, "criterion nonfinite") != 0;
    h.criterion.nonfinite_what = read_string(is, dir, "criterion nonfinite_what", bytes);
    h.criterion.res_ref = read_reals(is, dir, "criterion res_ref", bytes);
    h.criterion.res_last = read_reals(is, dir, "criterion res_last", bytes);
    const int nm = read_count(is, dir, "the monitor count", bytes);
    h.criterion.mon_hist.resize(nm);
    for (auto& m : h.criterion.mon_hist) { m = read_reals(is, dir, "criterion mon_hist", bytes); }
    for (const auto& r : h.records) {
        if (r.stage < 0 || r.stage >= na) { history_error(dir, "corrupt: a record's stage does not exist"); }
    }
    if (h.stage < 0 || h.stage >= na) { history_error(dir, "corrupt: the current stage does not exist"); }
    return h;
}

} // namespace ablsteady
