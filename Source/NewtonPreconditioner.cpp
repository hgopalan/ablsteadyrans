#include "NewtonPreconditioner.H"

#include "Convection.H"
#include "MassFlux.H"
#include "PseudoTime.H"

#include <AMReX_Print.H>

#include <cmath>
#include <string>

namespace ablsteady {

namespace {

//! The State field of each NewtonField.
amrex::MultiFab& state_field (State& q, int f)
{
    amrex::MultiFab* p[n_newton_fields] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

amrex::MultiFab& residual_field (ResidualFields& R, int f)
{
    amrex::MultiFab* p[n_newton_fields] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

//! to = from on every point (ghost cells and node heights included).
void copy_state (State& to, const State& from)
{
    auto t = to.fields();
    const auto f = const_cast<State&>(from).fields();
    for (std::size_t n = 0; n < t.size(); ++n) { amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect()); }
}

//! y = r on the points with mask = 1 (the Dirichlet rows), y elsewhere.
void dirichlet_rows_take_r (amrex::MultiFab& y, const amrex::MultiFab& r, const amrex::MultiFab& mask)
{
    for (amrex::MFIter mfi(y, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto ya = y.array(mfi);
        const auto ra = r.const_array(mfi);
        const auto ma = mask.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ya(i, j, k) = ma(i, j, k) * ra(i, j, k) + (amrex::Real(1) - ma(i, j, k)) * ya(i, j, k);
        });
    }
}

//! d = a_P on the faces with d > 0 (the SIMPLE denominator); the faces
//! with d = 0 (the Dirichlet rows: no correction) keep it.
void diagonal_denominators (amrex::MultiFab& d, const amrex::MultiFab& ap)
{
    for (amrex::MFIter mfi(d, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto da = d.array(mfi);
        const auto aa = ap.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (da(i, j, k) > 0) { da(i, j, k) = aa(i, j, k); }
        });
    }
}

PoissonOptions with_reltol (PoissonOptions popt, amrex::Real reltol)
{
    popt.reltol = reltol;
    return popt;
}

} // namespace

SimplecPreconditioner::SimplecPreconditioner (const SimplecOptions& opt, int sweeps, const ResidualInputs& in,
                                              const PoissonOptions& popt, const Mesh& mesh, const amrex::Geometry& geom,
                                              const amrex::MultiFab& z_nd, amrex::Real inf_cfl, const amrex::BoxArray& ba,
                                              const amrex::DistributionMapping& dm, const LineSweepOptions& sweep,
                                              PressureDenominator denominator)
    : m_opt(opt), m_sweeps(sweeps), m_sweep(sweep), m_denominator(denominator), m_in(in), m_mesh(mesh), m_geom(geom), m_z(z_nd), m_inf_cfl(inf_cfl),
      m_popt(with_reltol(popt, opt.pressure_reltol)),
      m_pbc(in.bc != nullptr ? poisson_bc_from_faces(in.bc->inputs(), geom) : PoissonBC{}),
      m_solver(mesh, geom, z_nd, m_popt)
{
    const std::string who = "ablsteady: SimplecPreconditioner: ";
    for (const auto& e : check_simplec_options(opt)) { amrex::Abort(who + e); }
    if (sweeps < 1) { amrex::Abort(who + "sweeps must be >= 1 (steady.zline_sweeps)"); }
    if (in.bc == nullptr) { amrex::Abort(who + "the residual inputs need the boundary conditions (bc)"); }
    if (in.closure == nullptr && in.nu == nullptr) {
        amrex::Abort(who + "without a closure the residual inputs need the diffusivities");
    }
    for (const auto& e : check_line_options(sweep)) { amrex::Abort(who + e); }
    // isnan first: an ordered comparison with NaN raises the invalid flag.
    if (std::isnan(inf_cfl) || !(std::isfinite(inf_cfl) && inf_cfl > 0)) {
        amrex::Abort(who + "inf_cfl (steady.cfl_max) must be finite and > 0");
    }
    for (const auto& e : check_poisson_bc(m_pbc, geom)) { amrex::Abort(who + e); }
    m_singular = true;
    for (const PoissonFace f : m_pbc.face) { m_singular = m_singular && f != PoissonFace::dirichlet; }
    for (State* s : {&m_q0, &m_d, &m_qn}) {
        s->define(ba, dm);
        amrex::MultiFab::Copy(s->z_phys_nd, z_nd, 0, 0, 1, s->z_phys_nd.nGrow());
    }
    const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
    const Loc locs[3] = {Loc::xface, Loc::yface, Loc::zface};
    for (int d = 0; d < 3; ++d) { in.bc->dirichlet_rows(eqs[d], locs[d], mesh, m_dir[d]); }
    m_R.define(ba, dm);
    if (in.closure != nullptr) { m_cf.define(ba, dm); }
    m_dtau.define(ba, dm, 1, 0);
    m_conv0.define(ba, dm, 1, 0);
}

void SimplecPreconditioner::setup (const State& q, const amrex::MultiFab* dtau)
{
    if (q.p.boxArray() != m_q0.p.boxArray() || q.p.DistributionMap() != m_q0.p.DistributionMap()) {
        amrex::Abort("ablsteady: SimplecPreconditioner::setup: q is not on the layout given to the constructor");
    }
    ++m_setups;
    // The linearization point with the stage's fill and wall state (step 1
    // of SimplecStage::iterate), and the closure's fields at it (step 2).
    copy_state(m_q0, q);
    SurfaceState* ss = (m_in.surface != nullptr) ? &m_ss : nullptr;
    m_in.bc->fill_ghosts(m_q0, m_in.surface, m_in.buoyancy, false, ss);
    if (m_in.closure != nullptr) { residual(m_q0, m_in, m_mesh, m_geom, m_R, &m_cf); }
    m_used_inf = (dtau == nullptr);
    if (dtau != nullptr) {
        check_operator_field("SimplecPreconditioner::setup", "dtau", *dtau, Loc::cell, m_mesh, 0);
        amrex::MultiFab::Copy(m_dtau, *dtau, 0, 0, 1, 0);
    } else {
        const amrex::MultiFab* nu = (m_in.closure != nullptr) ? &m_cf.nu : m_in.nu;
        local_dtau(m_q0, nu, m_mesh, m_geom, m_inf_cfl, m_dtau);
    }
    {
        // Once per setup: the bound of the sweeps (scalar_sweep checks it too).
        const amrex::Real lo = m_dtau.min(0);
        if (std::isnan(lo) || !(lo > 0) || !std::isfinite(m_dtau.max(0))) {
            amrex::Abort("ablsteady: SimplecPreconditioner::setup: dtau must be finite and > 0 everywhere");
        }
    }
    transport_fluxes(m_q0.u, m_q0.v, m_q0.w, m_z, m_mesh, m_geom, m_F0);
    face_dtau(m_dtau, m_mesh, m_geom, m_w.dtau_faces);
    const amrex::Array<const amrex::MultiFab*, 3> F0{&m_F0[0], &m_F0[1], &m_F0[2]};
    convection(Loc::cell, m_q0.theta, F0, m_in.theta, m_mesh, m_geom, m_conv0);
    m_ready = true;
}

void SimplecPreconditioner::apply (NewtonVector& dq, const NewtonVector& r)
{
    if (!m_ready) { amrex::Abort("ablsteady: SimplecPreconditioner::apply: setup has not been called"); }
    ++m_applies;
    const SurfaceState* ss = (m_in.surface != nullptr) ? &m_ss : nullptr;
    const ClosureFields* cf = (m_in.closure != nullptr) ? &m_cf : nullptr;

    // -r in place of the residual (the stage's R(q) for r = -R(q)).
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab& R = residual_field(m_R, f);
        amrex::MultiFab::Copy(R, r.f[f], 0, 0, 1, 0);
        R.mult(amrex::Real(-1), 0, 1, 0);
    }
    // 1. The momentum predictor on the increment: du* from zero. Its
    // point-implicit buoyancy reads theta's vertical gradient, so theta is
    // q's for this step.
    for (int f = 0; f < n_newton_fields; ++f) { state_field(m_d, f).setVal(0.0); }
    amrex::MultiFab::Copy(m_d.theta, m_q0.theta, 0, 0, 1, m_d.theta.nGrow());
    const amrex::Array<const amrex::MultiFab*, 3> F0{&m_F0[0], &m_F0[1], &m_F0[2]};
    m_info.predictor = momentum_predictor(m_d, m_R, F0, m_in, m_w.dtau_faces, m_mesh, m_geom, m_sweeps, m_pf, ss, cf,
                                          m_opt.rotation_relax, m_sweep, m_opt.w_vnn_cap);

    // 2. The pressure correction of du*, with r_p as the source: the
    // homogeneous fill gives du*'s ghost cells (theta back to 0 first).
    m_d.theta.setVal(0.0);
    m_in.bc->fill_ghosts(m_d, m_in.surface, m_in.buoyancy, true);
    if (m_denominator == PressureDenominator::diagonal) {
        diagonal_denominators(m_pf.d_u, m_pf.ap_u);
        diagonal_denominators(m_pf.d_v, m_pf.ap_v);
        diagonal_denominators(m_pf.d_w, m_pf.ap_w);
    }
    m_info.pressure = pressure_correction(m_d, m_pf, m_solver, m_pbc, m_opt.alpha_p, m_mesh, m_geom, m_w, &r.f[nf_p], true);

    // 3. theta's convection with the corrected velocity u + du (the stage's
    // fill of the corrected state), minus that at q.
    copy_state(m_qn, m_q0);
    amrex::MultiFab::Add(m_qn.u, m_d.u, 0, 0, 1, 0);
    amrex::MultiFab::Add(m_qn.v, m_d.v, 0, 0, 1, 0);
    amrex::MultiFab::Add(m_qn.w, m_d.w, 0, 0, 1, 0);
    m_in.bc->fill_ghosts(m_qn, m_in.surface, m_in.buoyancy, false, nullptr);
    transport_fluxes(m_qn.u, m_qn.v, m_qn.w, m_z, m_mesh, m_geom, m_w.transport);
    const amrex::Array<const amrex::MultiFab*, 3> F1{&m_w.transport[0], &m_w.transport[1], &m_w.transport[2]};
    if (m_w.conv1.empty()) { m_w.conv1.define(m_conv0.boxArray(), m_conv0.DistributionMap(), 1, 0); }
    convection(Loc::cell, m_qn.theta, F1, m_in.theta, m_mesh, m_geom, m_w.conv1);
    amrex::MultiFab::Subtract(m_w.conv1, m_conv0, 0, 0, 1, 0);
    amrex::MultiFab::Add(m_R.theta, m_w.conv1, 0, 0, 1, 0);

    // 4. The theta and k steps on the increments (zero, ghosts included),
    // with the corrected transport; no floor on dk.
    m_info.theta = scalar_sweep(Equation::theta, m_d.theta, m_R.theta, F1, m_z, m_in, m_dtau, m_mesh, m_geom, m_sweeps, 0, ss, cf,
                                m_w, m_sweep);
    m_info.tke = scalar_sweep(Equation::tke, m_d.tke, m_R.tke, F1, m_z, m_in, m_dtau, m_mesh, m_geom, m_sweeps, 0, ss, cf, m_w,
                              m_sweep, false);

    // 5. dq: the increment's valid points; every shared face takes its
    // owner's copy (the contract of the hook, Newton.H).
    for (int f = 0; f < n_newton_fields; ++f) { amrex::MultiFab::Copy(dq.f[f], state_field(m_d, f), 0, 0, 1, 0); }
    // The Dirichlet rows of u, v, w are identity rows of the Newton system
    // (Newton.H; the operator's fill sets those points of dq from the rest):
    // dq = r there.
    for (int f = 0; f < 3; ++f) { dirichlet_rows_take_r(dq.f[f], r.f[f], m_dir[f]); }
    for (int f = 0; f < 3; ++f) { dq.f[f].OverrideSync(m_geom.periodicity()); }
}

std::unique_ptr<NewtonPreconditioner> make_newton_preconditioner (NewtonPreconditionerKind kind, const SimplecOptions& opt,
                                                                  int sweeps, const ResidualInputs& in,
                                                                  const PoissonOptions& popt, const Mesh& mesh,
                                                                  const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                                                                  amrex::Real inf_cfl, const amrex::BoxArray& ba,
                                                                  const amrex::DistributionMapping& dm,
                                                                  const LineSweepOptions& sweep)
{
    if (kind == NewtonPreconditionerKind::identity) { return nullptr; }
    return std::make_unique<SimplecPreconditioner>(opt, sweeps, in, popt, mesh, geom, z_nd, inf_cfl, ba, dm, sweep);
}

} // namespace ablsteady
