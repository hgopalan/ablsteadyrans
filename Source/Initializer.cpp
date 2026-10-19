#include "Initializer.H"

#include "CVFluxes.H"
#include "Heights_K.H"
#include "MassFlux.H"
#include "Stencil_K.H"

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace ablsteady {

namespace {

const char* const kKeys[] = {"first_guess", "z_agl_min", "alpha_h", "alpha_v", "k_decay_height", "obrien", "projection",
                             "column_tol", "column_max_iterations"};

//! A profile component at every point (ghosts included) of a field of type
//! T, at the point's height above its column's ground floored at h_min.
//! comp: 0 u, 1 v, 2 theta, 3 k.
template <int T>
void fill_from_profile (amrex::MultiFab& mf, int comp, const InflowProfile::View& prof, amrex::Real h_min,
                        amrex::Real zi, const amrex::MultiFab& z_nd, int klo)
{
    // 1/z_i formed once, finite for z_i = 0 too (the decay is then off), so no
    // division sits in a branch a compiler may speculate under the traps.
    const amrex::Real zinv = amrex::Real(1) / amrex::max(zi, std::numeric_limits<amrex::Real>::min());
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real h = amrex::max(height_above_ground<T>(nz, i, j, k, klo), h_min);
            amrex::Real u, v, t, kk, nut;
            prof.eval(h, u, v, t, kk, nut);
            amrex::Real val = (comp == 0) ? u : (comp == 1) ? v : (comp == 2) ? t : kk;
            if (comp == 3 && zi > 0) {
                // k decays above z_i (the argument is 0 below it: exp(0) = 1)
                const amrex::Real arg = amrex::max(h - zi, amrex::Real(0)) * zinv;
                val *= std::exp(-arg);
            }
            a(i, j, k) = val;
        });
    }
}

} // namespace

InitOptions read_init_options ()
{
    InitOptions opt;
    amrex::ParmParse pp("init");
    for (const char* key : kKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { amrex::Abort(std::string("ablsteady: init.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { amrex::Abort(std::string("ablsteady: init.") + key + " takes exactly one value"); }
    }
    std::string fg = "none";
    pp.query("first_guess", fg);
    if (fg == "none") { opt.first_guess = InitOptions::FirstGuess::none; }
    else if (fg == "profile") { opt.first_guess = InitOptions::FirstGuess::profile; }
    else if (fg == "column") { opt.first_guess = InitOptions::FirstGuess::column; }
    else { amrex::Abort("ablsteady: init.first_guess = " + fg + ": must be none, profile or column"); }
    int proj = 1;
    pp.query("projection", proj);
    if (proj != 0 && proj != 1) { amrex::Abort("ablsteady: init.projection must be 0 or 1"); }
    opt.projection = (proj == 1);
    pp.query("column_tol", opt.column_tol);
    pp.query("column_max_iterations", opt.column_max_iterations);
    pp.query("z_agl_min", opt.z_agl_min);
    pp.query("alpha_h", opt.alpha_h);
    pp.query("alpha_v", opt.alpha_v);
    pp.query("k_decay_height", opt.k_decay_height);
    int ob = 0;
    pp.query("obrien", ob);
    if (ob != 0 && ob != 1) { amrex::Abort("ablsteady: init.obrien must be 0 or 1"); }
    opt.obrien = (ob == 1);
    const auto err = check_init_options(opt);
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

std::vector<std::string> check_init_options (const InitOptions& opt)
{
    std::vector<std::string> err;
    // isnan first, the ordered comparison inside the predicate (gcc traps on
    // NaN compares under fpe_trap_invalid).
    auto bad = [&] (amrex::Real x, const char* key, const char* range, auto ok) {
        if (std::isnan(x) || !std::isfinite(x) || !ok()) {
            err.emplace_back(std::string("init.") + key + " = " + std::to_string(x) + ": must be " + range);
        }
    };
    bad(opt.z_agl_min, "z_agl_min", "finite and >= 0 (m)", [&] { return opt.z_agl_min >= 0; });
    bad(opt.alpha_h, "alpha_h", "finite and > 0", [&] { return opt.alpha_h > 0; });
    bad(opt.alpha_v, "alpha_v", "finite and > 0", [&] { return opt.alpha_v > 0; });
    bad(opt.k_decay_height, "k_decay_height", "finite and >= 0 (m; 0: no decay)", [&] { return opt.k_decay_height >= 0; });
    bad(opt.column_tol, "column_tol", "finite and > 0", [&] { return opt.column_tol > 0; });
    if (opt.column_max_iterations < 1) { err.emplace_back("init.column_max_iterations must be >= 1"); }
    auto given = [&] (const char* k) { for (const auto& g : opt.given) { if (g == k) { return true; } } return false; };
    if (opt.first_guess == InitOptions::FirstGuess::none) {
        for (const auto& k : opt.given) {
            if (k != "first_guess") { err.emplace_back("init." + k + " has no effect without init.first_guess = profile or column"); }
        }
    } else {
        if (opt.first_guess != InitOptions::FirstGuess::column) {
            for (const char* k : {"column_tol", "column_max_iterations"}) {
                if (given(k)) { err.emplace_back(std::string("init.") + k + " has no effect without init.first_guess = column"); }
            }
        } else if (given("k_decay_height")) {
            err.emplace_back("init.k_decay_height would decay the column's converged k again: not with init.first_guess = column");
        }
        if (!opt.projection) {
            for (const char* k : {"alpha_h", "alpha_v"}) {
                if (given(k)) { err.emplace_back(std::string("init.") + k + " has no effect with init.projection = 0"); }
            }
        }
    }
    return err;
}

PoissonBC initializer_poisson_bc (const BoundaryInputs& in, const amrex::Geometry& geom)
{
    PoissonBC bc;
    for (int o = 0; o < 6; ++o) {
        const amrex::Orientation ori(o % 3, o < 3 ? amrex::Orientation::low : amrex::Orientation::high);
        if (geom.isPeriodic(ori.coordDir())) {
            bc.face[o] = PoissonFace::periodic;
        } else if (in.face[o] == FaceType::inflow || in.face[o] == FaceType::outflow) {
            bc.face[o] = PoissonFace::dirichlet;
        } else {
            bc.face[o] = PoissonFace::neumann;
        }
    }
    return bc;
}

amrex::Real obrien_omega (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& z_nd, const Mesh& mesh,
                          const amrex::Geometry& geom, amrex::MultiFab& omega)
{
    const char* op = "obrien_omega";
    check_operator_field(op, "u", u, Loc::xface, mesh, 0);
    check_operator_field(op, "v", v, Loc::yface, mesh, 0);
    check_nodal_heights(op, z_nd, mesh);
    const amrex::BoxArray zba = amrex::convert(mesh.detJ.boxArray(), amrex::IntVect(0, 0, 1));
    if (omega.empty()) { omega.define(zba, mesh.detJ.DistributionMap(), 1, 0); }
    check_operator_field(op, "omega", omega, Loc::zface, mesh, 0);
    // The horizontal part of face_divergence: (1/J) [d_xi (a_x u) + d_eta (a_y v)].
    amrex::MultiFab zero(zba, mesh.detJ.DistributionMap(), 1, 0), dh(mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    zero.setVal(0.0);
    face_divergence(mesh, geom, {&u, &v, &zero}, dh);
    // Column integration on the host (boxes span the column; the work is a
    // vanishing fraction of a solve): a_z Omega_{k+1/2} = a_z Omega_{k-1/2}
    // - J_k Dh_k dzeta from Omega = 0 at the ground (face_divergence divides
    // the vertical difference by dzeta), then the top residual removed with
    // (height above the ground / column depth)^2.
    amrex::Real worst = 0;
    const amrex::Real dzeta = geom.CellSize(2);
    const int klo = geom.Domain().smallEnd(2), nz = geom.Domain().length(2);
    for (amrex::MFIter mfi(omega, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox om(vb, 1, amrex::The_Pinned_Arena());
        const amrex::Box cb = mesh.detJ[mfi].box();
        amrex::FArrayBox J(cb, 1, amrex::The_Pinned_Arena());
        amrex::Box dcb = amrex::enclosedCells(vb);
        amrex::FArrayBox Dc(dcb, 1, amrex::The_Pinned_Arena());
        const amrex::Box azb = mesh.az[mfi].box();
        amrex::FArrayBox az(azb, 1, amrex::The_Pinned_Arena());
        const amrex::Box zb = z_nd[mfi].box();
        amrex::FArrayBox zn(zb, 1, amrex::The_Pinned_Arena());
        if (vb.smallEnd(2) != klo || vb.length(2) != nz + 1) {
            amrex::Abort("ablsteady: obrien_omega: every box must span the whole column (mesh.max_grid_size is horizontal)");
        }
        J.copy<amrex::RunOn::Device>(mesh.detJ[mfi], cb);
        Dc.copy<amrex::RunOn::Device>(dh[mfi], dcb);
        az.copy<amrex::RunOn::Device>(mesh.az[mfi], azb);
        zn.copy<amrex::RunOn::Device>(z_nd[mfi], zb);
        amrex::Gpu::streamSynchronize();
        amrex::Box cols = vb;
        cols.setRange(2, klo, 1);
        amrex::LoopOnCpu(cols, [&] (int i, int j, int) {
            // the node-mean height of the column's z faces
            auto zf = [&] (int k) {
                return amrex::Real(0.25) * (zn(amrex::IntVect(i, j, k)) + zn(amrex::IntVect(i + 1, j, k))
                                            + zn(amrex::IntVect(i, j + 1, k)) + zn(amrex::IntVect(i + 1, j + 1, k)));
            };
            om(amrex::IntVect(i, j, klo)) = 0;
            for (int k = klo; k < klo + nz; ++k) {
                const amrex::IntVect c(i, j, k), fl(i, j, k), fh(i, j, k + 1);
                om(fh) = (az(fl) * om(fl) - J(c) * Dc(c) * dzeta) / az(fh);
            }
            const amrex::Real E = om(amrex::IntVect(i, j, klo + nz));
            worst = std::max(worst, static_cast<amrex::Real>(std::abs(double(E))));
            // O'Brien's weight: quadratic in the height above the ground over
            // the column's depth (on a stretched mesh not the face index).
            const amrex::Real zg = zf(klo), depth = zf(klo + nz) - zg;
            for (int k = klo + 1; k <= klo + nz; ++k) {
                const amrex::Real frac = (zf(k) - zg) / depth;
                om(amrex::IntVect(i, j, k)) -= frac * frac * E;
            }
            om(amrex::IntVect(i, j, klo + nz)) = 0;   // frac = 1 there: exactly
        });
        omega[mfi].copy<amrex::RunOn::Device>(om, vb);
        amrex::Gpu::streamSynchronize();
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    return worst;
}

amrex::Real first_guess_state (const InitOptions& opt, const BoundaryInputs& bc, const Mesh& mesh,
                               const amrex::Geometry& geom, State& q, const CoriolisOptions& coriolis,
                               const InflowProfile* column)
{
    {
        const auto err = check_init_options(opt);
        if (!err.empty()) { amrex::Abort("ablsteady: first_guess_state: " + err.front()); }
    }
    if (opt.first_guess == InitOptions::FirstGuess::none) { return 0; }
    const char* op = "first_guess_state";
    check_operator_field(op, "u", q.u, Loc::xface, mesh, 2);
    check_operator_field(op, "theta", q.theta, Loc::cell, mesh, 2);
    check_nodal_heights(op, q.z_phys_nd, mesh);
    const BoundaryConditions bcs(bc, geom, mesh, coriolis);   // checks bc for this geometry (a geostrophic top needs the run's Coriolis)
    if (opt.first_guess == InitOptions::FirstGuess::column && column == nullptr) {
        amrex::Abort("ablsteady: first_guess_state: init.first_guess = column needs the column's profile");
    }
    const InflowProfile::View prof = (column != nullptr) ? column->view() : bcs.profile().view();
    const int klo = geom.Domain().smallEnd(2);
    const amrex::MultiFab& z = q.z_phys_nd;

    // 1. The first guess at every point of u and v (ghosts included); w = 0
    // except on the ground faces, which carry no flux (Omega_0 = 0 there): the
    // correction's ground flux is zero, so the constraint u.n = 0 at the
    // ground enters through the first guess (Sherman's problem minimises
    // |u - u_0| over the interior faces with the ground flux fixed).
    fill_from_profile<stencil::XF>(q.u, 0, prof, opt.z_agl_min, 0, z, klo);
    fill_from_profile<stencil::YF>(q.v, 1, prof, opt.z_agl_min, 0, z, klo);
    q.w.setVal(0.0);
    {
        amrex::MultiFab omega0(mesh.az.boxArray(), mesh.detJ.DistributionMap(), 1, 0);
        omega_from_w(q.u, q.v, q.w, z, mesh, geom, omega0, 0);
        for (amrex::MFIter mfi(omega0, false); mfi.isValid(); ++mfi) {
            amrex::Box ground = mfi.validbox();
            ground.setRange(2, klo, 1);
            omega0[mfi].template setVal<amrex::RunOn::Device>(amrex::Real(0), ground, 0, 1);
        }
        w_from_omega(q.u, q.v, omega0, z, mesh, geom, q.w);
    }
    // 2. theta from the sounding (absolute height; the reference itself, so
    // one ghost cell of the state is filled and the second stays 0 until the
    // boundary fill), k from the profile with the decay, p' = 0.
    fill_from_profile<stencil::CC>(q.tke, 3, prof, opt.z_agl_min, opt.k_decay_height, z, klo);
    {
        amrex::MultiFab theta_ref(q.p.boxArray(), q.p.DistributionMap(), 1, 1);
        fill_theta_ref(bc.theta_ref, z, theta_ref);
        q.theta.setVal(0.0);
        amrex::MultiFab::Copy(q.theta, theta_ref, 0, 0, 1, 1);
    }
    // p' = -int_z^top g (theta - theta_ref)/theta_0 dz is 0: theta is the
    // reference itself. (The integration comes with the first guess whose
    // theta differs from the sounding.)
    q.p.setVal(0.0);

    // 3. O'Brien: w0 from the column continuity with the top residual removed.
    amrex::Real residual = 0;
    if (opt.obrien) {
        amrex::MultiFab omega0;
        residual = obrien_omega(q.u, q.v, z, mesh, geom, omega0);
        w_from_omega(q.u, q.v, omega0, z, mesh, geom, q.w);
    }
    amrex::Gpu::streamSynchronize();   // the profile's device tables may go out of scope in the caller
    return residual;
}

PoissonOptions initializer_poisson_options (const InitOptions& opt, const PoissonOptions& popt)
{
    PoissonOptions p = popt;
    if (opt.alpha_v != opt.alpha_h && p.backend != "gmres_fft") { p.backend = "gmres_fft"; }
    return p;
}

InitResult initialize_state (const InitOptions& opt, const BoundaryInputs& bc, const PoissonOptions& popt, const Mesh& mesh,
                             const amrex::Geometry& geom, State& q, amrex::MultiFab* lambda, const CoriolisOptions& coriolis,
                             const InflowProfile* column)
{
    InitResult res;
    res.obrien_top_residual = first_guess_state(opt, bc, mesh, geom, q, coriolis, column);
    if (opt.first_guess == InitOptions::FirstGuess::none) { return res; }
    const amrex::MultiFab& z = q.z_phys_nd;
    if (!opt.projection) {
        // The first guess as it is: the stage's own pressure correction removes
        // the divergence (reported here as before = after).
        amrex::MultiFab div(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
        velocity_divergence(q.u, q.v, q.w, z, mesh, geom, div);
        res.divergence_before = res.divergence_after = div.norm2();
        return res;
    }

    // 4. The correction. The transport fluxes are Omega-based, so the
    // correction is applied to (u, v, Omega) and w follows.
    amrex::MultiFab rhs(q.p.boxArray(), q.p.DistributionMap(), 1, 0), div(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, rhs);
    res.divergence_before = rhs.norm2();
    amrex::MultiFab local_lambda;
    amrex::MultiFab* lam = (lambda != nullptr) ? lambda : &local_lambda;
    if (lam->empty()) { lam->define(q.p.boxArray(), q.p.DistributionMap(), 1, 1); }
    lam->setVal(0.0);
    const PoissonBC pbc = initializer_poisson_bc(bc, geom);
    const PoissonOptions ip = initializer_poisson_options(opt, popt);
    res.poisson_backend = ip.backend;
    res.backend_fallback = (ip.backend != popt.backend);
    PoissonSolver solver(mesh, geom, z, ip);
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    {
        const amrex::MultiFab* areas[3] = {&mesh.ax, &mesh.ay, &mesh.az};
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), mesh.detJ.DistributionMap(), 1, 0); }
    }
    FaceCoeffs coeffs;
    amrex::MultiFab bx, by, bz, bzc;
    const bool unit = (opt.alpha_h == 1 && opt.alpha_v == 1);
    if (!unit) {
        // Sigma = diag(alpha_h^2, alpha_h^2, alpha_v^2): x, y and the cross part
        // of the z faces carry alpha_h^2, the z part alpha_v^2 (Poisson.H).
        const amrex::DistributionMapping& dm = mesh.detJ.DistributionMap();
        bx.define(mesh.ax.boxArray(), dm, 1, 0);
        by.define(mesh.ay.boxArray(), dm, 1, 0);
        bz.define(mesh.az.boxArray(), dm, 1, 0);
        bzc.define(mesh.az.boxArray(), dm, 1, 0);
        const amrex::Real ah2 = opt.alpha_h * opt.alpha_h, av2 = opt.alpha_v * opt.alpha_v;
        bx.setVal(ah2);
        by.setVal(ah2);
        bz.setVal(av2);
        bzc.setVal(ah2);
        coeffs.x = &bx;
        coeffs.y = &by;
        coeffs.z = &bz;
        // With alpha_v = alpha_h the cross part takes the same coefficient as
        // the z part, so the split is not needed (and back end B can solve it).
        if (opt.alpha_v != opt.alpha_h) { coeffs.z_cross = &bzc; }
        coeffs.precond_ratio = av2 / ah2;
    }
    res.poisson = solver.solve(rhs, *lam, coeffs, pbc, &F);
    amrex::MultiFab omega(mesh.az.boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    omega_from_w(q.u, q.v, q.w, z, mesh, geom, omega, 0);
    amrex::MultiFab::Add(q.u, F[0], 0, 0, 1, 0);
    amrex::MultiFab::Add(q.v, F[1], 0, 0, 1, 0);
    amrex::MultiFab::Add(omega, F[2], 0, 0, 1, 0);
    // u and v with one ghost cell for w_from_omega and the divergence: the
    // periodic images; outside the domain the first-guess values stay.
    q.u.FillBoundary(geom.periodicity());
    q.v.FillBoundary(geom.periodicity());
    w_from_omega(q.u, q.v, omega, z, mesh, geom, q.w);
    velocity_divergence(q.u, q.v, q.w, z, mesh, geom, div);
    res.divergence_after = div.norm2();
    return res;
}

} // namespace ablsteady
