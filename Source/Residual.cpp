#include "Residual.H"

#include "Diffusion.H"
#include "ImplicitOperator.H"
#include "BoundaryConditions.H"
#include "MassFlux.H"

#include <AMReX_iMultiFab.H>

#include <cmath>
#include <limits>
#include <string>

namespace ablsteady {

void ResidualFields::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    using amrex::IntVect;
    u.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    v.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    w.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    p.define(ba, dm, 1, 0);
    theta.define(ba, dm, 1, 0);
    tke.define(ba, dm, 1, 0);
}

amrex::Vector<ResidualFields::Field> ResidualFields::fields ()
{
    return {{"u", &u}, {"v", &v}, {"w", &w}, {"p", &p}, {"theta", &theta}, {"tke", &tke}};
}

amrex::Vector<const amrex::MultiFab*> ResidualFields::fields () const
{
    return {&u, &v, &w, &p, &theta, &tke};
}

namespace {

void check_diffusivities (const char* op, const ResidualInputs& in)
{
    std::string err;
    if (in.closure != nullptr) {
        if (in.nu != nullptr || in.gamma_theta != nullptr || in.gamma_tke != nullptr) {
            err = "with a closure, nu, gamma_theta and gamma_tke must be null (the closure computes them)";
        } else if (in.wall_distance == nullptr) {
            err = "a closure needs wall_distance";
        }
    } else if (in.nu == nullptr) { err = "nu is not set"; }
    else if (in.gamma_theta == nullptr) { err = "gamma_theta is not set"; }
    else if (in.gamma_tke == nullptr) { err = "gamma_tke is not set"; }
    if (!err.empty()) { amrex::Abort(std::string("ablsteady: ") + op + ": ResidualInputs: " + err); }
}

//! The face average in use.
FaceMean face_mean_of (const ResidualInputs& in)
{
    return (in.closure != nullptr) ? in.closure->options().face_mean : in.face_mean;
}

//! Evaluate the closure at q into f (the surface layer's z0, or 0).
void evaluate_closure (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
                       ClosureFields& f)
{
    const amrex::Real z0 = (in.surface != nullptr) ? in.surface->z0 : amrex::Real(0);
    in.closure->evaluate(q, *in.wall_distance, z0, in.buoyancy, mesh, geom, f);
}

void check_inputs (const char* op, const ResidualInputs& in)
{
    check_diffusivities(op, in);
    std::string err;
    std::vector<std::string> all;
    for (const auto& e : check_reference_profile(in.theta_ref)) { all.push_back(e); }
    const auto& js = in.jvp_scales;
    for (amrex::Real x : {js.velocity, js.theta, js.tke, js.pressure}) {
        if (std::isnan(x) || !(x > 0 && std::isfinite(x))) {
            all.emplace_back("jvp_scales: every scale must be finite and > 0");
            break;
        }
    }
    for (const auto& e : check_convection_options(in.momentum)) { all.push_back("momentum " + e); }
    for (const auto& e : check_convection_options(in.theta)) { all.push_back("theta " + e); }
    for (const auto& e : check_buoyancy_options(in.buoyancy)) { all.push_back(e); }
    for (const auto& e : check_coriolis_options(in.coriolis)) { all.push_back(e); }
    if (in.surface != nullptr) {
        for (const auto& e : check_surface_layer_options(*in.surface)) { all.push_back(e); }
    }
    if (in.bc != nullptr) {
        const bool ground_wall = (in.bc->inputs().face[2] == FaceType::surface_layer);
        if (ground_wall != (in.surface != nullptr)) {
            all.emplace_back("bc.zlo = surface_layer needs ResidualInputs::surface, and other ground types must not have it");
        }
    }
    if (err.empty() && !all.empty()) { err = all.front(); }
    if (!err.empty()) { amrex::Abort(std::string("ablsteady: ") + op + ": ResidualInputs: " + err); }
}

// Over the owned points of q (shared faces once): sum of |(q - offset) / s|,
// sum of (v / s)^2 and the number of points, added to sums.
template <bool OFFSET>
void scaled_sums (const amrex::MultiFab& q, const amrex::MultiFab& v, const amrex::MultiFab* offset, amrex::Real s,
                  const amrex::Periodicity& per, double sums[3])
{
    const auto mask = q.OwnerMask(per);
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> ops;
    amrex::ReduceData<double, double, double> data(ops);
    const double sinv = 1.0 / static_cast<double>(s);
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const auto qa = q.const_array(mfi);
        const auto va = v.const_array(mfi);
        const auto m = mask->const_array(mfi);
        const auto o = OFFSET ? offset->const_array(mfi) : qa;
        ops.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<double, double, double>
        {
            if (m(i, j, k) == 0) { return {0.0, 0.0, 0.0}; }
            const double qq = OFFSET ? static_cast<double>(qa(i, j, k) - o(i, j, k)) : static_cast<double>(qa(i, j, k));
            const double vv = static_cast<double>(va(i, j, k)) * sinv;
            return {std::abs(qq) * sinv, vv * vv, 1.0};
        });
    }
    const auto r = data.value();
    sums[0] += amrex::get<0>(r);
    sums[1] += amrex::get<1>(r);
    sums[2] += amrex::get<2>(r);
}

amrex::Vector<const amrex::MultiFab*> solution_fields (const State& q)
{
    return {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
}

} // namespace

void residual (const State& q, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
               ResidualFields& R, ClosureFields* fields)
{
    check_inputs("residual", in);
    if (R.u.empty()) { R.define(q.p.boxArray(), q.p.DistributionMap()); }
    const amrex::MultiFab& z = q.z_phys_nd;
    // The diffusivities: given, or the closure's at q.
    ClosureFields local_fields;
    ClosureFields* cf = (fields != nullptr) ? fields : &local_fields;
    const amrex::MultiFab* nu = in.nu;
    const amrex::MultiFab* gamma_theta = in.gamma_theta;
    const amrex::MultiFab* gamma_tke = in.gamma_tke;
    if (in.closure != nullptr) {
        evaluate_closure(q, in, mesh, geom, *cf);
        nu = &cf->nu;
        gamma_theta = &cf->gamma_theta;
        gamma_tke = &cf->gamma_tke;
    }
    const FaceMean mean = face_mean_of(in);
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    transport_fluxes(q.u, q.v, q.w, z, mesh, geom, F);
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};

    // Momentum: C - D + G p - S.
    const amrex::Array<const amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    const amrex::Array<amrex::MultiFab*, 3> Rm{&R.u, &R.v, &R.w};
    amrex::Array<amrex::MultiFab, 3> tmp;
    for (int c = 0; c < 3; ++c) {
        convection(velocity_loc(c), *vel[c], Fp, in.momentum, mesh, geom, *Rm[c]);
        tmp[c].define(Rm[c]->boxArray(), Rm[c]->DistributionMap(), 1, 0);
    }
    // Diffusion as fluxes, so the wall model can set the ground fluxes.
    amrex::Array<amrex::Array<amrex::MultiFab, 3>, 3> Gm;
    amrex::Array<amrex::MultiFab, 3> Gt, Gk;
    momentum_diffusive_fluxes(vel, *nu, z, mesh, geom, Gm, mean);
    scalar_diffusive_fluxes(q.theta, *gamma_theta, z, mesh, geom, Gt, mean);
    scalar_diffusive_fluxes(q.tke, *gamma_tke, z, mesh, geom, Gk, mean);
    SurfaceState ss;
    if (in.surface != nullptr) {
        compute_surface_state(q, *in.surface, in.buoyancy, mesh, geom, ss);
        set_ground_fluxes(ss, *in.surface, Gm[0][2], Gm[1][2], Gt[2], Gk[2]);
    }
    for (int c = 0; c < 3; ++c) {
        flux_divergence(velocity_loc(c), Gm[c], mesh, geom, tmp[c]);
        amrex::MultiFab::Subtract(*Rm[c], tmp[c], 0, 0, 1, 0);
    }
    pressure_gradient(q.p, z, mesh, geom, {&tmp[0], &tmp[1], &tmp[2]});
    for (int c = 0; c < 3; ++c) { amrex::MultiFab::Add(*Rm[c], tmp[c], 0, 0, 1, 0); }
    amrex::MultiFab theta_ref(q.p.boxArray(), q.p.DistributionMap(), 1, 1);
    fill_theta_ref(in.theta_ref, z, theta_ref);
    subtract_buoyancy(q.theta, theta_ref, in.buoyancy, mesh, R.w);
    subtract_coriolis(q.u, q.v, q.w, in.coriolis, mesh, R.u, R.v, R.w);
    subtract_pressure_gradient(in.pressure_gradient, mesh, R.u, R.v, R.w);

    // Continuity.
    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, R.p);

    // Scalars: C - D.
    amrex::MultiFab d(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
    convection(Loc::cell, q.theta, Fp, in.theta, mesh, geom, R.theta);
    flux_divergence(Loc::cell, Gt, mesh, geom, d);
    amrex::MultiFab::Subtract(R.theta, d, 0, 0, 1, 0);
    ConvectionOptions kopt;
    kopt.limiter = Limiter::upwind;
    kopt.bounded = true;
    convection(Loc::cell, q.tke, Fp, kopt, mesh, geom, R.tke);
    flux_divergence(Loc::cell, Gk, mesh, geom, d);
    amrex::MultiFab::Subtract(R.tke, d, 0, 0, 1, 0);
    if (in.closure != nullptr) {
        // R_k -= S_k = P_s + P_b - eps (the wall-k rows below replace the first cell's).
        amrex::MultiFab::Subtract(R.tke, cf->Ps, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(R.tke, cf->Pb, 0, 0, 1, 0);
        amrex::MultiFab::Add(R.tke, cf->eps, 0, 0, 1, 0);
    }
    if (in.surface != nullptr) { apply_wall_k(ss, *in.surface, q.tke, R.tke); }
    if (in.bc != nullptr) {
        in.bc->add_damping(q, R);
        in.bc->apply_rows(q, R);
    }

    // A face shared by two boxes (or a periodic image) gets one value: the two
    // copies can differ in the last bits where a compiler contracts
    // multiply-adds or the mesh's shared nodes differ, and a per-box update
    // would let them drift apart.
    for (amrex::MultiFab* mf : {&R.u, &R.v, &R.w}) { mf->OverrideSync(geom.periodicity()); }
}

void line_jacobian (const State& q, const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh,
                    const amrex::Geometry& geom, amrex::MultiFab& tri)
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, mesh, geom, F);
    SurfaceState ss;
    const bool wall = (in.surface != nullptr && eq != Equation::w);
    if (wall) { compute_surface_state(q, *in.surface, in.buoyancy, mesh, geom, ss); }
    ClosureFields cf;
    if (in.closure != nullptr) { evaluate_closure(q, in, mesh, geom, cf); }
    line_jacobian({&F[0], &F[1], &F[2]}, q.z_phys_nd, in, eq, dir, mesh, geom, tri, wall ? &ss : nullptr,
                  in.closure != nullptr ? &cf : nullptr);
}

void implicit_rows (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& A, const SurfaceState* surface, const ClosureFields* fields)
{
    check_diffusivities("implicit_rows", in);
    const bool closure = (in.closure != nullptr);
    if (closure && (fields == nullptr || fields->empty())) {
        amrex::Abort("ablsteady: implicit_rows: with a closure the rows need its fields at the state (Closure::evaluate)");
    }
    Loc loc = Loc::cell;
    const amrex::MultiFab* g = closure ? &fields->nu : in.nu;
    switch (eq) {
    case Equation::u: loc = Loc::xface; break;
    case Equation::v: loc = Loc::yface; break;
    case Equation::w: loc = Loc::zface; break;
    case Equation::theta: g = closure ? &fields->gamma_theta : in.gamma_theta; break;
    case Equation::tke: g = closure ? &fields->gamma_tke : in.gamma_tke; break;
    }
    implicit_coefficients(loc, F, *g, z_nd, mesh, geom, A, face_mean_of(in));
    if (closure && eq == Equation::tke) {
        // Point-implicit sources: the Patankar diagonal (eps + max(-P_b, 0))/k.
        amrex::MultiFab::Add(A, fields->dSdk, 0, coef_P, 1, 0);
    }
    const bool wall = (in.surface != nullptr && eq != Equation::w);
    if (wall) {
        if (surface == nullptr) {
            amrex::Abort("ablsteady: implicit_rows: with a wall model the first cells' rows need the surface state");
        }
        apply_wall_to_implicit(*surface, loc, mesh, geom, A);
        if (eq == Equation::tke) { apply_wall_k_to_implicit(*surface, *in.surface, A); }
    }
    if (in.bc != nullptr) { in.bc->apply_to_implicit(eq, loc, z_nd, A); }
}

void line_jacobian (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                    const ResidualInputs& in, Equation eq, int dir, const Mesh& mesh, const amrex::Geometry& geom,
                    amrex::MultiFab& tri, const SurfaceState* surface, const ClosureFields* fields)
{
    if (dir < 0 || dir > 2) { amrex::Abort("ablsteady: line_jacobian: dir must be 0, 1 or 2"); }
    amrex::MultiFab A;
    implicit_rows(F, z_nd, in, eq, mesh, geom, A, surface, fields);
    const Loc loc = (eq == Equation::u) ? Loc::xface : (eq == Equation::v) ? Loc::yface : (eq == Equation::w) ? Loc::zface : Loc::cell;
    if (tri.empty()) { tri.define(A.boxArray(), A.DistributionMap(), n_line_coefs, 0); }
    check_operator_field("line_jacobian", "tri", tri, loc, mesh, 0);
    if (tri.nComp() != n_line_coefs) { amrex::Abort("ablsteady: line_jacobian: tri needs 3 components"); }
    amrex::MultiFab::Copy(tri, A, coef_XM + 2 * dir, line_lower, 1, 0);
    amrex::MultiFab::Copy(tri, A, coef_P, line_diag, 1, 0);
    amrex::MultiFab::Copy(tri, A, coef_XP + 2 * dir, line_upper, 1, 0);
    tri.mult(amrex::Real(-1), line_lower, 1, 0);
    tri.mult(amrex::Real(-1), line_upper, 1, 0);
}

amrex::Real jvp (const State& q, const State& v, const ResidualInputs& in, const Mesh& mesh,
                 const amrex::Geometry& geom, ResidualFields& Jv, const ResidualFields* R0)
{
    if (Jv.u.empty()) { Jv.define(q.p.boxArray(), q.p.DistributionMap()); }
    const auto qf = solution_fields(q);
    const auto vf = solution_fields(v);
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        if (vf[n]->boxArray() != qf[n]->boxArray() || vf[n]->DistributionMap() != qf[n]->DistributionMap()
            || vf[n]->nGrow() < qf[n]->nGrow()) {
            amrex::Abort("ablsteady: jvp: v must be a state on q's layout with q's ghost cells");
        }
    }
    check_inputs("jvp", in);
    const amrex::Periodicity per = geom.periodicity();
    amrex::MultiFab theta_ref(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
    fill_theta_ref(in.theta_ref, q.z_phys_nd, theta_ref);
    const auto& js = in.jvp_scales;
    const amrex::Real scale[6] = {js.velocity, js.velocity, js.velocity, js.pressure, js.theta, js.tke};
    double sums[3] = {0, 0, 0};
    for (int n = 0; n < 6; ++n) {
        if (n == 4) {
            scaled_sums<true>(*qf[n], *vf[n], &theta_ref, scale[n], per, sums);
        } else {
            scaled_sums<false>(*qf[n], *vf[n], nullptr, scale[n], per, sums);
        }
    }
    amrex::ParallelAllReduce::Sum(sums, 3, amrex::ParallelDescriptor::Communicator());
    if (sums[1] == 0) {
        for (auto& f : Jv.fields()) { f.mf->setVal(0.0); }
        return 0;
    }
    const double mean_q = sums[0] / sums[2];
    const double rms_v = std::sqrt(sums[1] / sums[2]);
    const double eps = std::sqrt(static_cast<double>(std::numeric_limits<amrex::Real>::epsilon())) * (1 + mean_q) / rms_v;

    ResidualFields base;
    if (R0 == nullptr) {
        residual(q, in, mesh, geom, base);
        R0 = &base;
    }
    State qp;
    qp.define(q.p.boxArray(), q.p.DistributionMap());
    amrex::MultiFab::Copy(qp.z_phys_nd, q.z_phys_nd, 0, 0, 1, q.z_phys_nd.nGrow());
    const amrex::Vector<amrex::MultiFab*> pf{&qp.u, &qp.v, &qp.w, &qp.p, &qp.theta, &qp.tke};
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        amrex::MultiFab::LinComb(*pf[n], 1, *qf[n], 0, static_cast<amrex::Real>(eps), *vf[n], 0, 0, 1, qf[n]->nGrow());
    }
    residual(qp, in, mesh, geom, Jv);
    auto jf = Jv.fields();
    const auto rf = R0->fields();
    for (int n = 0; n < static_cast<int>(jf.size()); ++n) {
        amrex::MultiFab::Subtract(*jf[n].mf, *rf[n], 0, 0, 1, 0);
        jf[n].mf->mult(static_cast<amrex::Real>(1 / eps), 0, 1, 0);
    }
    return static_cast<amrex::Real>(eps);
}

} // namespace ablsteady
