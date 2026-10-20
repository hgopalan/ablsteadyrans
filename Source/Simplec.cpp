#include "Simplec.H"

#include "CVFluxes.H"
#include "Convection.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace ablsteady {

// NaN is tested with isnan before any ordered comparison (gcc's ordered
// comparisons raise the invalid-operation trap on a NaN).

std::vector<std::string> check_simplec_options (const SimplecOptions& opt)
{
    std::vector<std::string> err;
    const amrex::Real eps10 = 10 * std::numeric_limits<amrex::Real>::epsilon();
    if (std::isnan(opt.pressure_reltol) || !(opt.pressure_reltol >= eps10 && opt.pressure_reltol < 1)) {
        err.emplace_back("simplec.pressure_reltol = " + std::to_string(opt.pressure_reltol)
                         + ": must be >= 10 machine epsilon and < 1 (the Poisson solver's own bound)");
    }
    if (std::isnan(opt.alpha_p) || !(opt.alpha_p > 0 && opt.alpha_p <= 1)) {
        err.emplace_back("simplec.alpha_p = " + std::to_string(opt.alpha_p) + ": must be > 0 and <= 1");
    }
    if (std::isnan(opt.rotation_relax) || !(opt.rotation_relax > 0 && opt.rotation_relax <= 1)) {
        err.emplace_back("simplec.rotation_relax = " + std::to_string(opt.rotation_relax) + ": must be > 0 and <= 1");
    }
    // The floor: caps of 4 and 8 burst on the Askervein hill with a 1 m first cell (Predictor.H).
    if (std::isnan(opt.w_vnn_cap) || !(opt.w_vnn_cap == 0 || (opt.w_vnn_cap >= w_vnn_cap_min && std::isfinite(opt.w_vnn_cap)))) {
        err.emplace_back("simplec.w_vnn_cap = " + std::to_string(opt.w_vnn_cap) + ": must be 0 (off), or finite and >= "
                         + std::to_string(static_cast<int>(w_vnn_cap_min)));
    }
    // Below 1 the bound itself is crossed (Predictor.H): no margin.
    if (std::isnan(opt.w_vnn_cap_slope_margin)
        || !(opt.w_vnn_cap_slope_margin == 0 || (opt.w_vnn_cap_slope_margin >= 1 && std::isfinite(opt.w_vnn_cap_slope_margin)))) {
        err.emplace_back("simplec.w_vnn_cap_slope_margin = " + std::to_string(opt.w_vnn_cap_slope_margin)
                         + ": must be 0 (off), or finite and >= 1");
    }
    return err;
}

SimplecOptions read_simplec_options ()
{
    SimplecOptions opt;
    amrex::ParmParse pp("simplec");
    std::vector<std::string> err;
    for (const char* key : {"pressure_reltol", "alpha_p", "rotation_relax", "w_vnn_cap", "w_vnn_cap_slope_margin"}) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("simplec.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { err.emplace_back(std::string("simplec.") + key + " takes exactly one value"); }
    }
    if (err.empty()) {
        pp.query("pressure_reltol", opt.pressure_reltol);
        pp.query("alpha_p", opt.alpha_p);
        pp.query("rotation_relax", opt.rotation_relax);
        pp.query("w_vnn_cap", opt.w_vnn_cap);
        pp.query("w_vnn_cap_slope_margin", opt.w_vnn_cap_slope_margin);
        err = check_simplec_options(opt);
        // Read but ignored without the cap: refuse it.
        const bool margin_given = std::find(opt.given.begin(), opt.given.end(), "w_vnn_cap_slope_margin") != opt.given.end();
        if (margin_given && opt.w_vnn_cap == 0) {
            err.emplace_back("simplec.w_vnn_cap_slope_margin has no effect unless simplec.w_vnn_cap > 0; remove it");
        }
    }
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

namespace {

//! beta = 1/d where d > 0 (keep = 1), beta = 1 and keep = 0 where d = 0.
//! Aborts on a negative or non-finite d.
void face_coefficients (const char* what, const amrex::MultiFab& d, amrex::MultiFab& beta, amrex::MultiFab& keep)
{
    if (beta.empty()) {
        beta.define(d.boxArray(), d.DistributionMap(), 1, 0);
        keep.define(d.boxArray(), d.DistributionMap(), 1, 0);
    }
    {
        const amrex::Real lo = d.min(0);   // collective: every rank aborts together
        if (std::isnan(lo) || !(lo >= 0) || !std::isfinite(d.max(0))) {
            amrex::Abort(std::string("ablsteady: pressure_correction: the SIMPLEC denominator ") + what
                         + " must be finite and >= 0 everywhere");
        }
    }
    for (amrex::MFIter mfi(beta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto b = beta.array(mfi);
        const auto kp = keep.array(mfi);
        const auto dd = d.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool on = dd(i, j, k) > 0;
            b(i, j, k) = on ? amrex::Real(1) / dd(i, j, k) : amrex::Real(1);
            kp(i, j, k) = on ? amrex::Real(1) : amrex::Real(0);
        });
    }
}

amrex::Real max_abs_masked (const amrex::MultiFab& f, const amrex::MultiFab& keep)
{
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real> data(op);
    for (amrex::MFIter mfi(f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = f.const_array(mfi);
        const auto kp = keep.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
        {
            return {std::abs(a(i, j, k) * kp(i, j, k))};
        });
    }
    amrex::Real m = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceRealMax(m);
    return m;
}

} // namespace

PressureCorrectionResult pressure_correction (State& q, const PredictorFields& d, PoissonSolver& solver, const PoissonBC& bc,
                                              amrex::Real alpha_p, const Mesh& mesh, const amrex::Geometry& geom,
                                              SimplecWorkspace& w)
{
    const char* op = "pressure_correction";
    if (std::isnan(alpha_p) || !(alpha_p > 0 && alpha_p <= 1)) {
        amrex::Abort("ablsteady: pressure_correction: alpha_p must be > 0 and <= 1");
    }
    check_operator_field(op, "d_u", d.d_u, Loc::xface, mesh, 0);
    check_operator_field(op, "d_v", d.d_v, Loc::yface, mesh, 0);
    check_operator_field(op, "d_w", d.d_w, Loc::zface, mesh, 0);
    check_operator_field(op, "u", q.u, Loc::xface, mesh, 1);
    check_operator_field(op, "v", q.v, Loc::yface, mesh, 1);
    check_operator_field(op, "w", q.w, Loc::zface, mesh, 1);
    check_operator_field(op, "p", q.p, Loc::cell, mesh, 1);
    const amrex::MultiFab& z = q.z_phys_nd;
    const amrex::BoxArray& ba = q.p.boxArray();
    const amrex::DistributionMapping& dm = q.p.DistributionMap();
    PressureCorrectionResult res;

    if (w.rhs.empty()) {
        w.rhs.define(ba, dm, 1, 0);
        w.div.define(ba, dm, 1, 0);
        w.phi.define(ba, dm, 1, 1);
        w.omega.define(mesh.az.boxArray(), dm, 1, 0);
        const amrex::MultiFab* areas[3] = {&mesh.ax, &mesh.ay, &mesh.az};
        for (int c = 0; c < 3; ++c) { w.F[c].define(areas[c]->boxArray(), dm, 1, 0); }
    }
    face_coefficients("d_u", d.d_u, w.beta_x, w.keep_x);
    face_coefficients("d_v", d.d_v, w.beta_y, w.keep_y);
    face_coefficients("d_w", d.d_w, w.beta_z, w.keep_z);

    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, w.rhs);
    res.divergence_before = w.rhs.norm0();
    // The scale of the divergence, U / h, for the net-source check.
    amrex::Real speed = amrex::max(q.u.norm0(0, 0, false), q.v.norm0(0, 0, false), q.w.norm0(0, 0, false));
    amrex::ParallelDescriptor::ReduceRealMax(speed);
    // The thinnest layer, J dzeta (J > 0 in the domain).
    const amrex::Real dz_min = mesh.detJ.min(0) * geom.CellSize(2);
    const amrex::Real inv_h = geom.InvCellSize(0) + geom.InvCellSize(1) + amrex::Real(1) / dz_min;

    FaceCoeffs coeffs;
    coeffs.x = &w.beta_x;
    coeffs.y = &w.beta_y;
    coeffs.z = &w.beta_z;
    coeffs.precond_ratio = 1;
    w.phi.setVal(0.0);
    res.poisson = solver.solve(w.rhs, w.phi, coeffs, bc, &w.F);
    if (res.poisson.singular && std::abs(res.poisson.removed_mean) > amrex::Real(1.0e-8) * speed * inv_h) {
        amrex::Abort("ablsteady: pressure_correction: the predicted velocity has a net source the boundaries cannot carry: "
                     "J-weighted mean divergence " + std::to_string(double(res.poisson.removed_mean)) + " 1/s against the scale "
                     + std::to_string(double(speed * inv_h)) + " 1/s (a boundary condition or the predictor is inconsistent)");
    }

    // p += alpha_p p'; u, v, Omega += the masked fluxes; w follows.
    amrex::MultiFab::Saxpy(q.p, alpha_p, w.phi, 0, 0, 1, 0);
    res.max_dp = alpha_p * w.phi.norm0();
    res.max_du = max_abs_masked(w.F[0], w.keep_x);
    res.max_dv = max_abs_masked(w.F[1], w.keep_y);
    res.max_dw = max_abs_masked(w.F[2], w.keep_z);
    for (int c = 0; c < 3; ++c) {
        const amrex::MultiFab& keep = (c == 0) ? w.keep_x : (c == 1) ? w.keep_y : w.keep_z;
        amrex::MultiFab::Multiply(w.F[c], keep, 0, 0, 1, 0);
    }
    omega_from_w(q.u, q.v, q.w, z, mesh, geom, w.omega, 0);
    amrex::MultiFab::Add(q.u, w.F[0], 0, 0, 1, 0);
    amrex::MultiFab::Add(q.v, w.F[1], 0, 0, 1, 0);
    amrex::MultiFab::Add(w.omega, w.F[2], 0, 0, 1, 0);
    // The two copies of a periodic seam face are one unknown (the fluxes on
    // them come from identical operands; this keeps them identical to the
    // bit). Then the periodic images of u and v for w_from_omega and the
    // divergence; the boundary conditions refill every ghost before the next
    // residual.
    q.u.OverrideSync(geom.periodicity());
    q.v.OverrideSync(geom.periodicity());
    w.omega.OverrideSync(geom.periodicity());
    q.u.FillBoundary(geom.periodicity());
    q.v.FillBoundary(geom.periodicity());
    w_from_omega(q.u, q.v, w.omega, z, mesh, geom, q.w);
    q.w.OverrideSync(geom.periodicity());
    q.w.FillBoundary(geom.periodicity());
    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, w.div);
    res.divergence_after = w.div.norm0();
    return res;
}

ScalarSweepResult scalar_sweep (Equation eq, amrex::MultiFab& field, const amrex::MultiFab& R_eq,
                                const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                                const ResidualInputs& in, const amrex::MultiFab& dtau, const Mesh& mesh,
                                const amrex::Geometry& geom, int sweeps, amrex::Real floor, const SurfaceState* surface,
                                const ClosureFields* fields, SimplecWorkspace& w, const LineSweepOptions& sweep)
{
    const char* op = "scalar_sweep";
    if (eq != Equation::theta && eq != Equation::tke) { amrex::Abort("ablsteady: scalar_sweep: eq must be theta or tke"); }
    check_operator_field(op, "field", field, Loc::cell, mesh, 1);
    check_operator_field(op, "R", R_eq, Loc::cell, mesh, 0);
    check_operator_field(op, "dtau", dtau, Loc::cell, mesh, 0);
    if (std::isnan(floor) || !(floor >= 0)) { amrex::Abort("ablsteady: scalar_sweep: the floor must be >= 0"); }
    {
        const amrex::Real lo = dtau.min(0);
        if (std::isnan(lo) || !(lo > 0) || !std::isfinite(dtau.max(0))) {
            amrex::Abort("ablsteady: scalar_sweep: dtau must be finite and > 0 everywhere");
        }
    }
    ScalarSweepResult res;
    implicit_rows(F, z_nd, in, eq, mesh, geom, w.A, surface, fields);
    const amrex::BoxArray& ba = field.boxArray();
    const amrex::DistributionMapping& dm = field.DistributionMap();
    if (w.ap.empty()) {
        w.ap.define(ba, dm, 1, 0);
        w.delta.define(ba, dm, 1, 1);
    }
    amrex::MultiFab rhs(ba, dm, 1, 0);
    amrex::MultiFab::Copy(rhs, R_eq, 0, 0, 1, 0);
    rhs.mult(amrex::Real(-1), 0, 1, 0);
    for (amrex::MFIter mfi(w.ap, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto p = w.ap.array(mfi);
        const auto a = w.A.const_array(mfi);
        const auto dt = dtau.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            p(i, j, k) = amrex::Real(1) / dt(i, j, k) + a(i, j, k, coef_P);
        });
    }
    res.last_sweep = line_sweeps(sweep, Loc::cell, w.A, w.ap, rhs, sweeps, mesh, geom, w.delta);
    res.max_dq = w.delta.norm0();
    amrex::MultiFab::Add(field, w.delta, 0, 0, 1, 0);
    if (eq == Equation::tke) {
        amrex::ReduceOps<amrex::ReduceOpSum> rop;
        amrex::ReduceData<amrex::Long> data(rop);
        for (amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto f = field.array(mfi);
            rop.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
            {
                const bool low = f(i, j, k) < floor;
                if (low) { f(i, j, k) = floor; }
                return {low ? amrex::Long(1) : amrex::Long(0)};
            });
        }
        res.clipped = amrex::get<0>(data.value(rop));
        amrex::ParallelDescriptor::ReduceLongSum(res.clipped);
    }
    return res;
}

namespace {

PoissonOptions with_reltol (PoissonOptions popt, amrex::Real reltol)
{
    popt.reltol = reltol;
    return popt;
}

} // namespace

SimplecStage::SimplecStage (const SimplecOptions& opt, int sweeps, const ResidualInputs& in, Closure* closure,
                            const PoissonOptions& popt, const Mesh& mesh, const amrex::Geometry& geom,
                            const amrex::MultiFab& z_nd, amrex::Real tke_floor, const amrex::BoxArray& ba,
                            const amrex::DistributionMapping& dm, const LineSweepOptions& sweep)
    : m_opt(opt), m_sweeps(sweeps), m_sweep(sweep), m_in(in), m_closure(closure), m_mesh(mesh), m_geom(geom), m_z(z_nd),
      m_tke_floor(tke_floor), m_popt(with_reltol(popt, opt.pressure_reltol)),
      m_pbc(in.bc != nullptr ? poisson_bc_from_faces(in.bc->inputs(), geom) : PoissonBC{}),
      m_solver(mesh, geom, z_nd, m_popt)
{
    for (const auto& e : check_simplec_options(opt)) { amrex::Abort("ablsteady: SimplecStage: " + e); }
    if (sweeps < 1) { amrex::Abort("ablsteady: SimplecStage: sweeps must be >= 1 (steady.zline_sweeps)"); }
    if (in.bc == nullptr) { amrex::Abort("ablsteady: SimplecStage: the residual inputs need the boundary conditions (bc)"); }
    if (in.closure != closure) {
        amrex::Abort("ablsteady: SimplecStage: closure must be the object the residual inputs point to (or both null)");
    }
    if (in.closure == nullptr && in.nu == nullptr) {
        amrex::Abort("ablsteady: SimplecStage: without a closure the residual inputs need the diffusivities");
    }
    if (std::isnan(tke_floor) || !(tke_floor >= 0)) { amrex::Abort("ablsteady: SimplecStage: tke_floor must be >= 0 (turb.tke_floor)"); }
    for (const amrex::Real sc : {in.jvp_scales.velocity, in.jvp_scales.theta, in.jvp_scales.tke}) {
        if (std::isnan(sc) || !(sc > 0 && std::isfinite(sc))) {
            amrex::Abort("ablsteady: SimplecStage: the residual scales (jvp_scales) must be finite and > 0");
        }
    }
    m_R.define(ba, dm);
    m_cf.define(ba, dm);
    m_cf.nu.setVal(0.0);
    m_height = geom.ProbHi(2) - geom.ProbLo(2);
    m_names = {"u", "v", "w", "p", "theta", "tke"};
    const JvpScales& s = in.jvp_scales;
    m_specs = {{"max_speed", s.velocity}, {"mean_tke", s.tke}};
    if (in.surface != nullptr) { m_specs.push_back({"mean_ustar", amrex::Real(0.1)}); }
    m_res.assign(m_names.size(), 0);
    m_mon.assign(m_specs.size(), 0);
}

void SimplecStage::init (State&)
{
    // The stage's state lives in q and the closure's history (both
    // checkpointed by the application). The driver's step limit reads nu()
    // before the first residual of the next iteration: the diffusivity of
    // the last accepted closure fields, nu_molecular + the history, so a
    // restart continues bit for bit (zero before the first iteration of a
    // fresh start, as the run it would continue).
    if (m_closure != nullptr && m_closure->has_history()) {
        const amrex::MultiFab& h = m_closure->history();
        amrex::MultiFab::Copy(m_cf.nut, h, 0, 0, 1, amrex::min(h.nGrow(), m_cf.nut.nGrow()));
        amrex::MultiFab::Copy(m_cf.nu, h, 0, 0, 1, amrex::min(h.nGrow(), m_cf.nu.nGrow()));
        m_cf.nu.plus(m_closure->options().nu_molecular, 0, 1, m_cf.nu.nGrow());
    }
}

void SimplecStage::compute_norms (const State& q)
{
    const JvpScales& s = m_in.jvp_scales;
    const amrex::Real U = s.velocity, H = m_height;
    const amrex::MultiFab* R[6] = {&m_R.u, &m_R.v, &m_R.w, &m_R.p, &m_R.theta, &m_R.tke};
    const amrex::Real scale[6] = {U * U / H, U * U / H, U * U / H, U / H, U * s.theta / H, U * s.tke / H};
    for (int n = 0; n < 6; ++n) {
        // Dot, not norm2 (which asserts cell-centred data); collective. The
        // face fields count a box-seam face twice and both periodic copies,
        // as numPts does: a consistent weighted RMS.
        const amrex::Real n2 = std::sqrt(amrex::MultiFab::Dot(*R[n], 0, 1, 0));
        const auto npts = static_cast<amrex::Real>(R[n]->boxArray().numPts());
        m_res[static_cast<std::size_t>(n)] = n2 / std::sqrt(npts) / scale[n];
    }
    amrex::Real speed = amrex::max(q.u.norm0(0, 0, false), q.v.norm0(0, 0, false), q.w.norm0(0, 0, false));
    amrex::ParallelDescriptor::ReduceRealMax(speed);
    m_mon[0] = speed;
    m_mon[1] = q.tke.sum(0, false) / static_cast<amrex::Real>(q.tke.boxArray().numPts());
    if (m_in.surface != nullptr) {
        const auto& dom = m_geom.Domain();
        m_mon[2] = m_ss.ustar.sum(0, false) / static_cast<amrex::Real>(dom.length(0) * dom.length(1));
    }
}

void SimplecStage::iterate (State& q, const amrex::MultiFab& dtau)
{
    check_operator_field("SimplecStage::iterate", "dtau", dtau, Loc::cell, m_mesh, 0);
    // 1. Ghosts and the wall state at q.
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy, false, m_in.surface != nullptr ? &m_ss : nullptr);
    // 2. The residual (the closure's fields at q).
    residual(q, m_in, m_mesh, m_geom, m_R, m_in.closure != nullptr ? &m_cf : nullptr);
    compute_norms(q);
    // 3. The momentum predictor.
    transport_fluxes(q.u, q.v, q.w, m_z, m_mesh, m_geom, m_w.transport);
    const amrex::Array<const amrex::MultiFab*, 3> F{&m_w.transport[0], &m_w.transport[1], &m_w.transport[2]};
    face_dtau(dtau, m_mesh, m_geom, m_w.dtau_faces);
    const SurfaceState* ss = (m_in.surface != nullptr) ? &m_ss : nullptr;
    const ClosureFields* cf = (m_in.closure != nullptr) ? &m_cf : nullptr;
    m_info.predictor = momentum_predictor(q, m_R, F, m_in, m_w.dtau_faces, m_mesh, m_geom, m_sweeps, m_pf, ss, cf, m_opt.rotation_relax,
                                          m_sweep, m_opt.w_vnn_cap, m_opt.w_vnn_cap_slope_margin);
    // 4. The pressure correction on u* (its ghosts refilled first).
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy, false, nullptr);
    m_info.pressure = pressure_correction(q, m_pf, m_solver, m_pbc, m_opt.alpha_p, m_mesh, m_geom, m_w);
    // 5. theta and k. theta's convection is re-evaluated with the corrected
    // velocity's transport (the second half of the point-implicit buoyancy
    // pair, Predictor.H); its diffusion, k's residual and the closure's
    // sources stay those of the start of the iteration. The implicit rows of
    // both use the corrected transport.
    if (m_w.conv0.empty()) {
        m_w.conv0.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
        m_w.conv1.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
    }
    convection(Loc::cell, q.theta, F, m_in.theta, m_mesh, m_geom, m_w.conv0);
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy, false, nullptr);
    transport_fluxes(q.u, q.v, q.w, m_z, m_mesh, m_geom, m_w.transport);
    convection(Loc::cell, q.theta, F, m_in.theta, m_mesh, m_geom, m_w.conv1);
    amrex::MultiFab::Subtract(m_w.conv1, m_w.conv0, 0, 0, 1, 0);
    amrex::MultiFab::Add(m_R.theta, m_w.conv1, 0, 0, 1, 0);
    m_info.theta = scalar_sweep(Equation::theta, q.theta, m_R.theta, F, m_z, m_in, dtau, m_mesh, m_geom, m_sweeps, 0, ss, cf, m_w,
                                m_sweep);
    m_info.tke = scalar_sweep(Equation::tke, q.tke, m_R.tke, F, m_z, m_in, dtau, m_mesh, m_geom, m_sweeps, m_tke_floor, ss, cf, m_w,
                              m_sweep);
    // 6. The nu_t history.
    if (m_closure != nullptr) { m_closure->accept(m_cf); }
}

} // namespace ablsteady
