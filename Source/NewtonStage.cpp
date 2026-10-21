#include "NewtonStage.H"

#include "BoundaryConditions.H"
#include "PseudoTime.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <fstream>
#include <sstream>

namespace ablsteady {

namespace {

//! The continuation of the stage: the SER law starts at the cap
//! (steady.cfl_max), not at steady.cfl0, so once the first step reaches it
//! the continuation stays at infinite dtau (SIMPLEC's ramp is done).
ContinuationOptions at_cap (ContinuationOptions o)
{
    o.cfl0 = o.cfl_max;
    return o;
}

const char* const kMagic = "ablsteady-newton-stage";
constexpr int kVersion = 1;

void state_error (const std::string& dir, const std::string& what)
{
    amrex::Abort("ablsteady: " + dir + "/NewtonStage: " + what);
}

bool exists_everywhere (const std::string& name)
{
    int e = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(name) ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&e, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return e == 1;
}

std::string sci (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", x);
    return b;
}

} // namespace

NewtonStageOptions read_newton_stage_options (const SteadyOptions& sopt, const ClosureOptions& copt)
{
    NewtonStageOptions o;
    o.newton = read_newton_options();
    o.continuation = read_continuation_options(sopt);
    o.smoothness = smoothness_options(copt);
    return o;
}

NewtonStage::NewtonStage (const SteadyOptions& sopt, const NewtonStageOptions& opt, const SimplecOptions& simplec,
                          const PoissonOptions& popt, const ResidualInputs& in, Closure* closure, const Mesh& mesh,
                          const amrex::Geometry& geom, const amrex::MultiFab& z_nd, amrex::Real tke_floor, amrex::Real norm_height, const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm, NewtonPreconditioner* hook, int verbose)
    : m_sopt(sopt), m_opt(opt), m_in(in), m_nin(newton_residual_inputs(in, opt.smoothness)), m_closure(closure),
      m_mesh(mesh), m_geom(geom), m_height(norm_height),
      m_cont(at_cap(opt.continuation), opt.newton, m_nin, mesh, geom, in.closure != nullptr ? &m_cf.nu : in.nu, tke_floor, ba,
             dm, verbose),
      m_pc(hook), m_smooth(opt.smoothness), m_verbose(verbose)
{
    if (in.bc == nullptr) { amrex::Abort("ablsteady: NewtonStage: the residual inputs need the boundary conditions (bc)"); }
    if (in.closure != closure) {
        amrex::Abort("ablsteady: NewtonStage: closure must be the object the residual inputs point to (or both null)");
    }
    if (std::isnan(norm_height) || !(norm_height > 0 && std::isfinite(norm_height))) {
        amrex::Abort("ablsteady: NewtonStage: the norm height must be finite and > 0");
    }
    if (m_pc == nullptr) {
        m_own_pc = make_newton_preconditioner(sopt.newton_preconditioner, simplec, sopt.zline_sweeps, m_nin, popt, mesh, geom,
                                              z_nd, sopt.cfl_max, ba, dm, line_sweep_options(sopt));
        m_pc = m_own_pc.get();
    }
    m_R.define(ba, dm);
    m_cf.define(ba, dm);
    m_names = stage_residual_names();
    m_specs = stage_monitor_specs(in);
    m_res.assign(m_names.size(), 0);
    m_mon.assign(m_specs.size(), 0);
}

void NewtonStage::init (State& q)
{
    m_request = false;
    m_reason.clear();
    m_failures = m_restored ? m_failures : 0;
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy);
    if (!m_smooth.active()) {
        m_smooth.enter(q, m_closure, m_in, m_mesh, m_geom);
        // A fresh continuation at the SIMPLEC cap: with newton_dtau_inf_at_cap
        // its first accepted step goes to infinite dtau (plain Newton).
        ContinuationState st;
        st.cfl = m_opt.continuation.cfl_max;
        m_cont.restore(st);
    }
    m_restored = false;
    if (m_closure != nullptr && m_closure->frozen()) {
        const amrex::iMultiFab& sw = m_closure->frozen_switches();
        m_saved.clear();
        m_saved.define(sw.boxArray(), sw.DistributionMap(), 1, sw.nGrow());
        amrex::iMultiFab::Copy(m_saved, sw, 0, 0, 1, sw.nGrow());
    }
    run_norms(q);   // the closure's fields at q: the continuation's diffusivity
    if (m_verbose > 0) {
        amrex::Print() << "ablsteady: Newton stage: " << (m_nin.momentum.limiter != m_in.momentum.limiter ? "smooth limiter, " : "")
                       << (m_closure != nullptr && m_closure->frozen() ? "k floor frozen, " : "") << "continuation from CFL "
                       << m_cont.state().cfl << (m_cont.state().infinite ? " (infinite)" : "") << ", "
                       << newton_preconditioner_name(m_sopt.newton_preconditioner) << " preconditioner\n";
    }
}

void NewtonStage::iterate (State& q, const amrex::MultiFab& /*dtau*/)
{
    m_request = false;
    const ContinuationState& before = m_cont.state();
    const amrex::Real drop = std::pow(amrex::Real(10), -m_sopt.criterion.residual_drop_orders);
    if (before.steps > 0 && before.r0 > 0 && before.r_last <= drop * before.r0) {   // after an accepted step
        // Converged on the stage's equations (the run's own: differs() was
        // false): no more steps, which would only fail at round-off; the
        // driver's criterion fills its monitor window at the cost of one
        // residual evaluation per iteration.
        m_status = ContinuationStatus::accepted;
        m_converged = true;
        run_norms(q);
        return;
    }
    m_converged = false;
    m_status = m_cont.step(q, m_pc);
    m_failures = (m_status == ContinuationStatus::failed) ? m_failures + 1 : 0;
    run_norms(q);
    if (m_failures >= m_sopt.newton_fallback_failures) {
        m_request = true;
        m_reason = std::to_string(m_failures) + " consecutive failed Newton steps (steady.newton_fallback_failures "
                   + std::to_string(m_sopt.newton_fallback_failures) + ")";
        return;
    }
    const ContinuationState& st = m_cont.state();
    if (m_status == ContinuationStatus::accepted && st.steps > 0 && st.r0 > 0 && st.r_last <= drop * st.r0 && differs(q)) {
        m_request = true;
        m_handed_back = true;
        char b[64];
        std::snprintf(b, sizeof(b), "%g", static_cast<double>(m_sopt.criterion.residual_drop_orders));
        m_reason = "Newton converged on its own equations (|R| " + sci(double(st.r_last)) + ", " + b
                   + " orders below |R_0| " + sci(double(st.r0))
                   + "), which differ from the run's (limiter or k floor): the final SIMPLEC pass";
    }
}

void NewtonStage::leave (State& q)
{
    m_request = false;
    m_failures = 0;
    if (!m_smooth.active()) { return; }
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy);
    const amrex::Long moved = m_smooth.leave(q, m_closure, m_in, m_mesh, m_geom);
    m_saved.clear();
    if (m_verbose > 0) {
        amrex::Print() << "ablsteady: Newton stage left after " << m_cont.state().steps << " steps; " << moved
                       << " k-floor decisions moved\n";
    }
}

bool NewtonStage::differs (const State& q) const
{
    if (m_nin.momentum.limiter != m_in.momentum.limiter || m_nin.theta.limiter != m_in.theta.limiter) { return true; }
    if (m_closure == nullptr || !m_closure->frozen()) { return false; }
    const amrex::Real z0 = (m_in.surface != nullptr) ? m_in.surface->z0 : amrex::Real(0);
    return m_closure->changed_switches(q, *m_in.wall_distance, z0, m_in.buoyancy, m_mesh, m_geom) > 0;
}

void NewtonStage::run_norms (State& q)
{
    // The accepted point's ghosts are filled (Continuation.H); the wall state
    // for the monitors.
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy, false, m_in.surface != nullptr ? &m_ss : nullptr);
    const bool frozen = (m_closure != nullptr && m_closure->frozen());
    if (frozen) { m_closure->thaw(); }   // the run's own (live) decisions
    residual(q, m_in, m_mesh, m_geom, m_R, m_in.closure != nullptr ? &m_cf : nullptr);
    if (frozen) { m_closure->restore_switches(m_saved, FreezeScope::floor, q.p.boxArray(), q.p.DistributionMap(), m_geom); }
    stage_norms(m_R, q, m_in.surface != nullptr ? &m_ss : nullptr, m_in.jvp_scales, m_height, m_geom, m_res, m_mon);
}

void NewtonStage::write_checkpoint (const std::string& dir) const
{
    m_smooth.write_checkpoint(dir, m_closure);
    if (m_smooth.active()) { write_continuation_state(dir, m_cont.state()); }
    int ok = 1;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(dir + "/NewtonStage");
        os << kMagic << " " << kVersion << "\n";
        os << "failures " << m_failures << " handed_back " << (m_handed_back ? 1 : 0) << "\n";
        ok = os.good() ? 1 : 0;
    }
    amrex::ParallelDescriptor::Bcast(&ok, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (ok == 0) { state_error(dir, "could not be written"); }
    amrex::ParallelDescriptor::Barrier();
}

void NewtonStage::read_checkpoint (const std::string& dir, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                                   const amrex::Geometry& geom)
{
    m_smooth.read_checkpoint(dir, m_closure, ba, dm, geom);
    if (!exists_everywhere(dir + "/NewtonStage")) {
        // A checkpoint from a run that never entered the stage (or older).
        m_failures = 0;
        m_handed_back = false;
        return;
    }
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(dir + "/NewtonStage", buf);
    buf.push_back('\0');
    std::istringstream is(buf.data());
    std::string magic, w1, w2;
    int version = 0, failures = -1, back = -1;
    is >> magic >> version;
    if (!is || magic != kMagic || version != kVersion) {
        state_error(dir, "not an ablsteady Newton stage state of format version " + std::to_string(kVersion));
    }
    if (!(is >> w1 >> failures >> w2 >> back) || w1 != "failures" || w2 != "handed_back" || failures < 0
        || (back != 0 && back != 1)) {
        state_error(dir, "corrupt: bad state line");
    }
    m_failures = failures;
    m_handed_back = (back == 1);
    if (m_smooth.active()) {
        m_cont.restore(read_continuation_state(dir));
        m_restored = true;
    }
}

} // namespace ablsteady
