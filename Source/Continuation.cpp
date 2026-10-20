#include "Continuation.H"

#include "BoundaryConditions.H"
#include "PseudoTime.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace ablsteady {

namespace {

constexpr int kContinuationVersion = 1;
const char* const kContinuationMagic = "ablsteady-continuation";

std::string num (double x)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<amrex::Real>::digits10) << x;
    return s.str();
}

std::string hexfloat (double x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", x);
    return buf;
}

std::string sci (double x)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4e", x);
    return buf;
}

void state_error (const std::string& dir, const std::string& what)   // aborts
{
    amrex::Abort("ablsteady: " + dir + "/Continuation: " + what);
}

amrex::MultiFab* solution_field (State& q, int f)
{
    amrex::MultiFab* p[n_newton_fields] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return p[f];
}

//! to = from on every point of the six fields and the heights, ghosts included.
void copy_state (State& to, State& from)
{
    auto t = to.fields();
    const auto f = from.fields();
    for (std::size_t n = 0; n < t.size(); ++n) {
        amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect());
    }
}

//! k = max(k, floor) on the valid points.
void clip_floor (amrex::MultiFab& k, amrex::Real kmin)
{
    for (amrex::MFIter mfi(k, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto a = k.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            a(i, j, kk) = amrex::max(a(i, j, kk), kmin);
        });
    }
}

const char* outcome_name (AttemptOutcome o)
{
    switch (o) {
    case AttemptOutcome::accepted: return "accepted";
    case AttemptOutcome::rejected: return "rejected (no lambda >= steady.newton_lambda_min met the Armijo condition)";
    default: return "skipped (the Newton update is not usable)";
    }
}

} // namespace

// ---------------------------------------------------------------------------
ContinuationOptions continuation_options (const SteadyOptions& steady)
{
    ContinuationOptions o;
    o.cfl0 = steady.cfl0;
    o.cfl_max = steady.cfl_max;
    o.ser_exponent = steady.ser_exponent;
    o.ser_max_growth = steady.ser_max_growth;
    return o;
}

std::vector<std::string> check_continuation_options (const ContinuationOptions& opt)
{
    // The SER keys: the driver's own checks (the other fields at their defaults).
    SteadyOptions s;
    s.cfl0 = opt.cfl0;
    s.cfl_max = opt.cfl_max;
    s.ser_exponent = opt.ser_exponent;
    s.ser_max_growth = opt.ser_max_growth;
    std::vector<std::string> err = check_steady_options(s);
    // isnan first: an ordered compare with NaN raises the invalid flag (gcc).
    auto open01 = [&] (amrex::Real x, const char* key, double hi) {
        const auto d = static_cast<double>(x);
        if (std::isnan(d) || !(d > 0) || !(d < hi)) {
            err.emplace_back(std::string("steady.") + key + " = " + num(d) + ": must be > 0 and < " + num(hi));
        }
    };
    open01(opt.armijo_c, "newton_armijo_c", 0.5);
    open01(opt.backtrack, "newton_backtrack", 1);
    open01(opt.lambda_min, "newton_lambda_min", 1);
    open01(opt.cfl_cut, "newton_cfl_cut", 1);
    if (opt.max_retries < 1) {
        err.emplace_back("steady.newton_max_retries = " + std::to_string(opt.max_retries) + ": must be an integer >= 1");
    }
    if (opt.dtau_inf_at_cap != 0 && opt.dtau_inf_at_cap != 1) {
        err.emplace_back("steady.newton_dtau_inf_at_cap = " + std::to_string(opt.dtau_inf_at_cap) + ": must be 0 or 1");
    }
    return err;
}

ContinuationOptions read_continuation_options (const SteadyOptions& steady)
{
    ContinuationOptions opt = continuation_options(steady);
    std::vector<std::string> err;
    amrex::ParmParse pp("steady");
    for (const char* key : kContinuationInputKeys) {
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
    pp.query("newton_armijo_c", opt.armijo_c);
    pp.query("newton_backtrack", opt.backtrack);
    pp.query("newton_lambda_min", opt.lambda_min);
    pp.query("newton_cfl_cut", opt.cfl_cut);
    pp.query("newton_max_retries", opt.max_retries);
    pp.query("newton_dtau_inf_at_cap", opt.dtau_inf_at_cap);
    err = check_continuation_options(opt);
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

// ---------------------------------------------------------------------------
// isfinite first: NaN fails the comparison anyway, but an ordered compare
// with NaN raises the invalid flag (the trap of amrex.fpe_trap_invalid).
bool armijo_accepts (amrex::Real r_trial, amrex::Real r, amrex::Real lambda, amrex::Real c)
{
    return std::isfinite(static_cast<double>(r_trial)) && r_trial <= (amrex::Real(1) - c * lambda) * r;
}

LineSearchResult armijo_line_search (const std::function<amrex::Real (amrex::Real)>& phi, amrex::Real r,
                                     const ContinuationOptions& opt)
{
    LineSearchResult ls;
    ls.lambda = 1;
    while (true) {
        ls.residual = phi(ls.lambda);
        if (armijo_accepts(ls.residual, r, ls.lambda, opt.armijo_c)) {
            ls.accepted = true;
            return ls;
        }
        const amrex::Real next = opt.backtrack * ls.lambda;
        if (!(next >= opt.lambda_min)) { return ls; }
        ls.lambda = next;
        ++ls.backtracks;
    }
}

ContinuationState continuation_after_accept (const ContinuationOptions& opt, ContinuationState st, amrex::Real r_new)
{
    st.steps += 1;
    st.r_last = r_new;
    const amrex::Real r = (st.r0 > 0) ? r_new / st.r0 : amrex::Real(0);
    const amrex::Real from = st.infinite ? opt.cfl_max : st.cfl;
    st.cfl = ser_cfl(from, r, opt.cfl0, opt.cfl_max, opt.ser_exponent, opt.ser_max_growth);
    st.infinite = opt.dtau_inf_at_cap == 1 && st.cfl >= opt.cfl_max;
    return st;
}

ContinuationState continuation_after_failure (const ContinuationOptions& opt, ContinuationState st)
{
    st.cfl = opt.cfl_cut * (st.infinite ? opt.cfl_max : st.cfl);
    st.infinite = false;
    return st;
}

// ---------------------------------------------------------------------------
void write_continuation_state (const std::string& dir, const ContinuationState& st)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        amrex::ParallelDescriptor::Barrier();
        return;
    }
    const std::string name = dir + "/Continuation";
    std::ofstream os(name);
    if (!os.good()) { amrex::FileOpenFailed(name); }
    os << kContinuationMagic << " " << kContinuationVersion << "\n";
    os << "state " << st.steps << " " << hexfloat(static_cast<double>(st.cfl)) << " " << (st.infinite ? 1 : 0) << " "
       << hexfloat(static_cast<double>(st.r0)) << " " << hexfloat(static_cast<double>(st.r_last)) << "\n";
    if (!os.good()) { state_error(dir, "could not be written"); }
    os.close();
    amrex::ParallelDescriptor::Barrier();
}

ContinuationState read_continuation_state (const std::string& dir)
{
    // Every rank learns of a missing file before the broadcast.
    int exists = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(dir + "/Continuation") ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&exists, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (exists == 0) { state_error(dir, "does not exist"); }
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(dir + "/Continuation", buf);
    buf.push_back('\0');
    std::istringstream is(buf.data());
    std::string magic, word;
    int version = 0;
    is >> magic >> version;
    if (!is || magic != kContinuationMagic || version != kContinuationVersion) {
        state_error(dir, "not an ablsteady continuation state of format version " + std::to_string(kContinuationVersion));
    }
    // Hex floats through strtod: libstdc++'s stream extraction does not take them.
    auto real = [&] (const char* what) {
        std::string tok;
        is >> tok;
        char* end = nullptr;
        const double x = tok.empty() ? 0.0 : std::strtod(tok.c_str(), &end);
        if (!is || tok.empty() || *end != '\0') { state_error(dir, std::string("corrupt: bad value for ") + what); }
        return static_cast<amrex::Real>(x);
    };
    ContinuationState st;
    int inf = 0;
    if (!(is >> word) || word != "state" || !(is >> st.steps)) { state_error(dir, "corrupt: bad state line"); }
    st.cfl = real("cfl");
    if (!(is >> inf) || (inf != 0 && inf != 1)) { state_error(dir, "corrupt: bad value for infinite"); }
    st.infinite = inf == 1;
    st.r0 = real("r0");
    st.r_last = real("r_last");
    return st;
}

// ---------------------------------------------------------------------------
Continuation::Continuation (const ContinuationOptions& copt, const NewtonOptions& nopt, const ResidualInputs& in,
                            const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab* nu,
                            amrex::Real tke_floor, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                            int verbose)
    : m_opt(copt),
      m_nk(nopt, in, mesh, geom, ba, dm, verbose > 1 ? 1 : 0),
      m_mesh(mesh),
      m_geom(geom),
      m_nu(nu),
      m_tke_floor(tke_floor),
      m_verbose(verbose)
{
    const auto err = check_continuation_options(copt);
    if (!err.empty()) {
        std::string all = "ablsteady: Continuation: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    if (std::isnan(static_cast<double>(tke_floor)) || !std::isfinite(static_cast<double>(tke_floor)) || !(tke_floor >= 0)) {
        amrex::Abort("ablsteady: Continuation: tke_floor must be finite and >= 0 (turb.tke_floor)");
    }
    m_st.cfl = copt.cfl0;
    m_dq.define(ba, dm);
    m_qb.define(ba, dm);
    m_qt.define(ba, dm);
    m_dtau.define(ba, dm, 1, 0);
}

void Continuation::restore (const ContinuationState& st)
{
    auto finite_nonneg = [] (amrex::Real x) {
        const auto d = static_cast<double>(x);
        return std::isfinite(d) && d >= 0;
    };
    if (st.steps < 0 || !finite_nonneg(st.cfl) || !(st.cfl > 0) || !finite_nonneg(st.r0) || !finite_nonneg(st.r_last)) {
        amrex::Abort("ablsteady: Continuation::restore: the state is inconsistent (steps >= 0, CFL finite and > 0, |R| finite and >= 0)");
    }
    m_st = st;
    // The options may have changed across the restart.
    m_st.cfl = amrex::min(m_st.cfl, m_opt.cfl_max);
}

void Continuation::log (const ContinuationRecord& rec) const
{
    if (m_verbose <= 0 && rec.outcome == AttemptOutcome::accepted) { return; }
    amrex::Print() << "ablsteady: PTC step " << rec.step << " attempt " << rec.attempt << ": |R| "
                   << sci(static_cast<double>(rec.residual)) << " -> " << sci(static_cast<double>(rec.residual_new))
                   << ", CFL " << (rec.infinite ? std::string("inf") : sci(static_cast<double>(rec.cfl)))
                   << ", lambda " << rec.lambda << ", " << rec.backtracks << " backtracks, GMRES "
                   << rec.gmres_iterations << " iterations (status " << rec.gmres_status << "), Krylov relative "
                   << sci(static_cast<double>(rec.krylov_relative)) << ": " << outcome_name(rec.outcome) << "\n";
}

ContinuationStatus Continuation::step (State& q, NewtonPreconditioner* pc)
{
    const ResidualInputs& in = m_nk.residual_inputs();
    for (int attempt = 0; attempt <= m_opt.max_retries; ++attempt) {
        ContinuationRecord rec;
        rec.step = m_st.steps;
        rec.attempt = attempt;
        rec.cfl = m_st.cfl;
        rec.infinite = m_st.infinite;
        // Every attempt starts from the caller's q (a copy: on terrain the
        // fill moves the ground w face, so q itself is left untouched until a
        // step is accepted).
        copy_state(m_qb, q);
        in.bc->fill_ghosts(m_qb, in.surface, in.buoyancy);
        const amrex::MultiFab* dtau = nullptr;
        if (!m_st.infinite) {
            local_dtau(m_qb, m_nu, m_mesh, m_geom, m_st.cfl, m_dtau);
            dtau = &m_dtau;
        }
        const NewtonStepInfo info = m_nk.direction(m_qb, dtau, m_dq, pc);
        rec.residual = info.residual_norm;
        rec.residual_new = info.residual_norm;
        rec.gmres_iterations = info.gmres_iterations;
        rec.gmres_status = info.gmres_status;
        rec.krylov_relative = info.krylov_relative;
        const bool finite_r = std::isfinite(static_cast<double>(info.residual_norm));
        if (m_st.r0 == 0 && finite_r) { m_st.r0 = info.residual_norm; }
        if (finite_r && info.residual_norm == 0) {
            // A solution: nothing to do (dq = 0).
            rec.outcome = AttemptOutcome::accepted;
            m_records.push_back(rec);
            log(rec);
            m_st = continuation_after_accept(m_opt, m_st, 0);
            return ContinuationStatus::accepted;
        }
        if (info.applied) {
            const amrex::Real r = info.residual_norm;
            auto phi = [&] (amrex::Real lambda) {
                copy_state(m_qt, m_qb);
                for (int f = 0; f < n_newton_fields; ++f) {
                    amrex::MultiFab::Saxpy(*solution_field(m_qt, f), lambda, *solution_field(m_dq, f), 0, 0, 1, 0);
                }
                clip_floor(m_qt.tke, m_tke_floor);
                return m_nk.residual_norm(m_qt);   // fills m_qt
            };
            const LineSearchResult ls = armijo_line_search(phi, r, m_opt);
            rec.lambda = ls.lambda;
            rec.backtracks = ls.backtracks;
            rec.residual_new = ls.residual;
            if (ls.accepted) {
                copy_state(q, m_qt);
                amrex::Gpu::streamSynchronize();
                rec.outcome = AttemptOutcome::accepted;
                m_records.push_back(rec);
                log(rec);
                m_st = continuation_after_accept(m_opt, m_st, ls.residual);
                return ContinuationStatus::accepted;
            }
            rec.outcome = AttemptOutcome::rejected;
        } else {
            rec.outcome = AttemptOutcome::skipped;
        }
        m_records.push_back(rec);
        log(rec);
        m_st = continuation_after_failure(m_opt, m_st);
    }
    // q is as it was given.
    amrex::Print() << "ablsteady: PTC step " << m_st.steps << ": failed after " << m_opt.max_retries
                   << " retries (steady.newton_max_retries); the next CFL is " << sci(static_cast<double>(m_st.cfl)) << "\n";
    return ContinuationStatus::failed;
}

ContinuationRunResult Continuation::run (State& q, int max_steps, amrex::Real rel_tol, NewtonPreconditioner* pc)
{
    ContinuationRunResult out;
    for (int n = 0; n < max_steps; ++n) {
        if (m_st.steps > 0 && m_st.r_last <= rel_tol * m_st.r0) {
            out.converged = true;
            return out;
        }
        out.status = step(q, pc);
        if (out.status == ContinuationStatus::failed) { return out; }
        ++out.steps;
    }
    out.converged = m_st.steps > 0 && m_st.r_last <= rel_tol * m_st.r0;
    return out;
}

} // namespace ablsteady
