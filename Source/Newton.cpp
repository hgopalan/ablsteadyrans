#include "Newton.H"

#include "BoundaryConditions.H"
#include "Predictor.H"

#include <AMReX_GMRES.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace gust {

namespace {


//! x with the digits of amrex::Real: a float input prints as it was given.
std::string num (double x)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<amrex::Real>::digits10) << x;
    return s.str();
}

//! x with every digit of amrex::Real: a bound printed this way is accepted
//! when given back.
std::string num_exact (double x)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << x;
    return s.str();
}

//! The State field of each NewtonField.
amrex::MultiFab& state_field (State& q, int f)
{
    amrex::MultiFab* p[n_newton_fields] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

const amrex::MultiFab& residual_field (const ResidualFields& R, int f)
{
    const amrex::MultiFab* p[n_newton_fields] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

//! y = (inv_dtau x_s + Jv) / t + dir x on the valid points, the kernels of
//! NewtonOperator::apply: with or without the pseudo-time term, with or
//! without Dirichlet rows (the face velocities).
void combine (amrex::MultiFab& y, const amrex::MultiFab& x, const amrex::MultiFab& dq, const amrex::MultiFab& Jv,
              const amrex::MultiFab* inv_dtau, const amrex::MultiFab* dir, amrex::Real tinv)
{
    for (amrex::MFIter mfi(y, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto ya = y.array(mfi);
        const auto xa = x.const_array(mfi);
        const auto da = dq.const_array(mfi);
        const auto ja = Jv.const_array(mfi);
        if (inv_dtau != nullptr && dir != nullptr) {
            const auto ia = inv_dtau->const_array(mfi);
            const auto ma = dir->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ya(i, j, k) = (ia(i, j, k) * da(i, j, k) + ja(i, j, k)) * tinv + ma(i, j, k) * xa(i, j, k);
            });
        } else if (inv_dtau != nullptr) {
            const auto ia = inv_dtau->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ya(i, j, k) = (ia(i, j, k) * da(i, j, k) + ja(i, j, k)) * tinv;
            });
        } else if (dir != nullptr) {
            const auto ma = dir->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ya(i, j, k) = ja(i, j, k) * tinv + ma(i, j, k) * xa(i, j, k);
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ya(i, j, k) = ja(i, j, k) * tinv;
            });
        }
    }
}

//! inv = (1 - dir) / dtau (dir null: no Dirichlet rows).
void inverse_dtau (amrex::MultiFab& inv, const amrex::MultiFab& dtau, const amrex::MultiFab* dir)
{
    for (amrex::MFIter mfi(inv, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto ia = inv.array(mfi);
        const auto da = dtau.const_array(mfi);
        if (dir != nullptr) {
            const auto ma = dir->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ia(i, j, k) = (amrex::Real(1) - ma(i, j, k)) / da(i, j, k);
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ia(i, j, k) = amrex::Real(1) / da(i, j, k);
            });
        }
    }
}

//! The sum over the owned points of a box of a * b, in double: one thread
//! per (i, j) column adds its k values in order, the host adds the columns
//! in storage order. The same bits for the same box whatever the rank.
double box_dot (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::iMultiFab* owner,
                const amrex::MFIter& mfi)
{
    // mfi iterates over the cells; the field's own (face) box of the same index.
    const amrex::Box bx = amrex::convert(mfi.validbox(), a.ixType());
    const int klo = bx.smallEnd(2);
    const int khi = bx.bigEnd(2);
    amrex::Box cols = bx;
    cols.setBig(2, klo);
    const auto npts = static_cast<std::size_t>(cols.numPts());
    amrex::Gpu::DeviceVector<double> colsum(npts);
    double* cs = colsum.data();
    const auto lo = amrex::lbound(cols);
    const auto len = amrex::length(cols);
    const auto aa = a.const_array(mfi);
    const auto ba = b.const_array(mfi);
    if (owner != nullptr) {
        const auto m = owner->const_array(mfi);
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            double s = 0;
            for (int k = klo; k <= khi; ++k) {
                if (m(i, j, k) != 0) { s += static_cast<double>(aa(i, j, k)) * static_cast<double>(ba(i, j, k)); }
            }
            cs[(i - lo.x) + (j - lo.y) * len.x + (k0 - lo.z) * len.x * len.y] = s;
        });
    } else {
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            double s = 0;
            for (int k = klo; k <= khi; ++k) { s += static_cast<double>(aa(i, j, k)) * static_cast<double>(ba(i, j, k)); }
            cs[(i - lo.x) + (j - lo.y) * len.x + (k0 - lo.z) * len.x * len.y] = s;
        });
    }
    amrex::Vector<double> h(npts);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, colsum.begin(), colsum.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    double s = 0;
    for (const double x : h) { s += x; }
    return s;
}

} // namespace

std::vector<std::string> check_newton_options (const NewtonOptions& opt)
{
    std::vector<std::string> err;
    if (opt.gmres_restart < 5 || opt.gmres_restart > 200) {
        err.emplace_back("steady.newton_gmres_restart = " + std::to_string(opt.gmres_restart)
                         + ": must be an integer >= 5 and <= 200");
    }
    if (opt.gmres_max_iters < 1) {
        err.emplace_back("steady.newton_gmres_max_iters = " + std::to_string(opt.gmres_max_iters)
                         + ": must be an integer >= 1");
    }
    // isnan first: NaN fails every ordered comparison.
    // The lower bound, machine epsilon to the power 3/4 (1.8e-12 in double,
    // 6.4e-6 in single): there J v keeps about three digits (relative error
    // near 1e-3). Below it the round-off grows as 1/h: at machine epsilon J v
    // has no correct digits and Newton stalls or diverges on the column of
    // Newton.ColumnConvergesQuadratically, and further down q + eps v rounds
    // back to q until J v is 0 and GMRES reports a converged zero step.
    const auto e = static_cast<double>(opt.fd_eps_scale);
    const double e_min = newton_fd_eps_scale_min();
    if (std::isnan(e) || !(e >= e_min) || !(e <= 1.0e-2)) {
        err.emplace_back("steady.newton_fd_eps_scale = " + num(e) + ": must be >= " + num_exact(e_min)
                         + " (machine epsilon^(3/4)) and <= 1e-2");
    }
    const auto r = static_cast<double>(opt.krylov_rtol);
    if (std::isnan(r) || !(r > 0) || !(r < 1)) {
        err.emplace_back("steady.newton_krylov_rtol = " + num(r) + ": must be > 0 and < 1");
    }
    const auto b = static_cast<double>(opt.mem_budget_gb);
    if (std::isnan(b) || !(b > 0) || !(b <= 80)) {
        err.emplace_back("steady.newton_mem_budget_gb = " + num(b) + ": must be > 0 and <= 80 (GB of one 80 GB H100)");
    }
    return err;
}

NewtonOptions read_newton_options ()
{
    NewtonOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("steady");
    for (const char* key : kNewtonInputKeys) {
        if (!pp.contains(key)) { continue; }
        opt.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("steady.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { err.emplace_back(std::string("steady.") + key + " takes exactly one value"); }
    }
    if (!err.empty()) {
        std::string all = "gust: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    pp.query("newton_gmres_restart", opt.gmres_restart);
    pp.query("newton_gmres_max_iters", opt.gmres_max_iters);
    pp.query("newton_fd_eps_scale", opt.fd_eps_scale);
    pp.query("newton_krylov_rtol", opt.krylov_rtol);
    pp.query("newton_mem_budget_gb", opt.mem_budget_gb);
    err = check_newton_options(opt);
    if (!err.empty()) {
        std::string all = "gust: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return opt;
}

void NewtonVector::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    using amrex::IntVect;
    f[nf_u].define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    f[nf_v].define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    f[nf_w].define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    f[nf_p].define(ba, dm, 1, 0);
    f[nf_theta].define(ba, dm, 1, 0);
    f[nf_tke].define(ba, dm, 1, 0);
}

NewtonOperator::NewtonOperator (const State& q, const ResidualFields& R0, const ResidualInputs& in, const Mesh& mesh,
                                const amrex::Geometry& geom, const amrex::MultiFab* dtau, NewtonPreconditioner* pc)
    : m_q(q), m_R0(R0), m_in(in), m_mesh(mesh), m_geom(geom), m_pc(pc),
      m_ba(q.p.boxArray()), m_dm(q.p.DistributionMap())
{
    if (in.bc == nullptr) {
        amrex::Abort("gust: NewtonOperator: ResidualInputs::bc is null (the perturbation's homogeneous ghost fill needs the boundary conditions)");
    }
    const auto& js = in.jvp_scales;
    m_s = {js.velocity, js.velocity, js.velocity, js.pressure, js.theta, js.tke};
    m_t = {js.velocity, js.velocity, js.velocity, amrex::Real(1), js.theta, js.tke};

    NewtonVector tmp;
    tmp.define(m_ba, m_dm);
    for (int d = 0; d < 3; ++d) { m_owner[d] = tmp.f[d].OwnerMask(geom.periodicity()); }

    // The Dirichlet rows of the face velocities (identity rows, no pseudo-time term).
    const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
    const Loc locs[3] = {Loc::xface, Loc::yface, Loc::zface};
    for (int d = 0; d < 3; ++d) { in.bc->dirichlet_rows(eqs[d], locs[d], mesh, m_dir[d]); }

    if (dtau != nullptr) {
        m_has_dtau = true;
        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> dtf;
        face_dtau(*dtau, mesh, geom, dtf);
        for (int d = 0; d < 3; ++d) {
            m_inv_dtau[d].define(tmp.f[d].boxArray(), m_dm, 1, 0);
            inverse_dtau(m_inv_dtau[d], dtf[d], &m_dir[d]);
        }
        for (const int f : {nf_theta, nf_tke}) {
            m_inv_dtau[f].define(m_ba, m_dm, 1, 0);
            inverse_dtau(m_inv_dtau[f], *dtau, nullptr);
        }
        amrex::Gpu::streamSynchronize();   // dtf goes out of scope
    }
    m_dq.define(m_ba, m_dm);
    amrex::MultiFab::Copy(m_dq.z_phys_nd, q.z_phys_nd, 0, 0, 1, q.z_phys_nd.nGrow());
    m_Jv.define(m_ba, m_dm);
    if (m_pc != nullptr) {
        m_pr.define(m_ba, m_dm);
        m_pz.define(m_ba, m_dm);
    }
}

void NewtonOperator::apply (NewtonVector& Ax, const NewtonVector& x)
{
    // dq = S x on the valid points, then its ghost cells (and Dirichlet face
    // points) from the homogeneous fill, before jvp reads them (its step
    // reads the valid points). jvp then refills q + eps dq in full, so the
    // ghost fill's own dependence on q (the wall law's ground gradient) is
    // differentiated too.
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab& d = state_field(m_dq, f);
        amrex::MultiFab::Copy(d, x.f[f], 0, 0, 1, 0);
        d.mult(m_s[f], 0, 1, 0);
    }
    m_in.bc->fill_ghosts(m_dq, m_in.surface, m_in.buoyancy, true);
    const amrex::Real eps = jvp(m_q, m_dq, m_in, m_mesh, m_geom, m_Jv, &m_R0, true);
    if (eps > 0) {
        m_eps_min = (m_eps_min > 0) ? std::min(m_eps_min, eps) : eps;
        m_eps_max = std::max(m_eps_max, eps);
    }
    ++m_jvp_calls;
    for (int f = 0; f < n_newton_fields; ++f) {
        const amrex::MultiFab* inv = (m_has_dtau && !m_inv_dtau[f].empty()) ? &m_inv_dtau[f] : nullptr;
        const amrex::MultiFab* dir = (f < 3) ? &m_dir[f] : nullptr;
        combine(Ax.f[f], x.f[f], state_field(m_dq, f), residual_field(m_Jv, f), inv, dir, amrex::Real(1) / m_t[f]);
    }
}

void NewtonOperator::precond (NewtonVector& z, const NewtonVector& r)
{
    if (m_pc == nullptr) {
        assign(z, r);
        return;
    }
    // z = S^-1 M T r: the hook works in physical units.
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab::Copy(m_pr.f[f], r.f[f], 0, 0, 1, 0);
        m_pr.f[f].mult(m_t[f], 0, 1, 0);
    }
    m_pc->apply(m_pz, m_pr);
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab::Copy(z.f[f], m_pz.f[f], 0, 0, 1, 0);
        z.f[f].mult(amrex::Real(1) / m_s[f], 0, 1, 0);
    }
}

NewtonVector NewtonOperator::makeVecRHS () const
{
    NewtonVector v;
    v.define(m_ba, m_dm);
    return v;
}

NewtonVector NewtonOperator::makeVecLHS () const
{
    return makeVecRHS();
}

void NewtonOperator::assign (NewtonVector& lhs, const NewtonVector& rhs) const
{
    for (int f = 0; f < n_newton_fields; ++f) { amrex::MultiFab::Copy(lhs.f[f], rhs.f[f], 0, 0, 1, 0); }
}

void NewtonOperator::increment (NewtonVector& lhs, const NewtonVector& rhs, RT a) const
{
    for (int f = 0; f < n_newton_fields; ++f) { amrex::MultiFab::Saxpy(lhs.f[f], a, rhs.f[f], 0, 0, 1, 0); }
}

void NewtonOperator::linComb (NewtonVector& lhs, RT a, const NewtonVector& rhs_a, RT b,
                              const NewtonVector& rhs_b) const
{
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab::LinComb(lhs.f[f], a, rhs_a.f[f], 0, b, rhs_b.f[f], 0, 0, 1, 0);
    }
}

void NewtonOperator::scale (NewtonVector& v, RT fac) const
{
    for (int f = 0; f < n_newton_fields; ++f) { v.f[f].mult(fac, 0, 1, 0); }
}

void NewtonOperator::setToZero (NewtonVector& v) const
{
    for (int f = 0; f < n_newton_fields; ++f) { v.f[f].setVal(0.0); }
}

NewtonOperator::RT NewtonOperator::dotProduct (const NewtonVector& a, const NewtonVector& b) const
{
    // One slot per box (every field of the box, in field order), all-reduced
    // (one contributor per slot: exact), added in box order.
    const int nboxes = static_cast<int>(m_ba.size());
    amrex::Vector<double> box_sum(nboxes, 0.0);
    for (amrex::MFIter mfi(a.f[nf_p], false); mfi.isValid(); ++mfi) {
        double s = 0;
        for (int f = 0; f < n_newton_fields; ++f) {
            s += box_dot(a.f[f], b.f[f], (f < 3) ? m_owner[f].get() : nullptr, mfi);
        }
        box_sum[mfi.index()] = s;
    }
    amrex::ParallelDescriptor::ReduceRealSum(box_sum.dataPtr(), nboxes);
    double total = 0;
    for (const double x : box_sum) { total += x; }
    return static_cast<RT>(total);
}

NewtonOperator::RT NewtonOperator::norm2 (const NewtonVector& v) const
{
    return std::sqrt(dotProduct(v, v));
}

void NewtonOperator::rhs (const ResidualFields& R, NewtonVector& b) const
{
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab::Copy(b.f[f], residual_field(R, f), 0, 0, 1, 0);
        b.f[f].mult(-amrex::Real(1) / m_t[f], 0, 1, 0);
    }
}

void NewtonOperator::add_update (State& dq, const NewtonVector& x) const
{
    for (int f = 0; f < n_newton_fields; ++f) {
        amrex::MultiFab::Saxpy(state_field(dq, f), m_s[f], x.f[f], 0, 0, 1, 0);
    }
}

KrylovMemory krylov_memory (int restart, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    KrylovMemory m;
    m.vectors = 2 * restart + 3;
    m.arnoldi_vectors = restart + 1;
    amrex::Vector<double> rank_bytes(amrex::ParallelDescriptor::NProcs(), 0.0);
    for (int i = 0; i < static_cast<int>(ba.size()); ++i) {
        const amrex::Box b = ba[i];
        double pts = 3.0 * static_cast<double>(b.numPts());
        for (int d = 0; d < 3; ++d) { pts += static_cast<double>(amrex::surroundingNodes(b, d).numPts()); }
        const double bytes = pts * static_cast<double>(sizeof(amrex::Real));
        m.vector_bytes += bytes;
        rank_bytes[dm[i]] += bytes;
    }
    const auto it = std::max_element(rank_bytes.begin(), rank_bytes.end());
    m.max_rank = static_cast<int>(it - rank_bytes.begin());
    m.total_bytes = m.vector_bytes * m.vectors;
    m.max_rank_bytes = *it * m.vectors;
    return m;
}

// isfinite first: an ordered comparison with NaN raises the invalid flag
// (and so the trap of amrex.fpe_trap_invalid); isfinite does not.
bool newton_update_usable (amrex::Real krylov_relative)
{
    return std::isfinite(static_cast<double>(krylov_relative)) && krylov_relative <= 1;
}

const char* newton_preconditioner_name (NewtonPreconditionerKind kind) noexcept
{
    return (kind == NewtonPreconditionerKind::simplec) ? "simplec" : "identity";
}

double newton_fd_eps_scale_min ()
{
    return std::pow(static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()), 0.75);
}

NewtonKrylov::NewtonKrylov (const NewtonOptions& opt, const ResidualInputs& in, const Mesh& mesh,
                            const amrex::Geometry& geom, const amrex::BoxArray& ba,
                            const amrex::DistributionMapping& dm, int verbose)
    : m_opt(opt), m_in(in), m_mesh(mesh), m_geom(geom), m_verbose(verbose)
{
    const auto err = check_newton_options(opt);
    if (!err.empty()) {
        std::string all = "gust: NewtonKrylov: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    if (in.bc == nullptr) {
        amrex::Abort("gust: NewtonKrylov: ResidualInputs::bc is null (the perturbation's homogeneous ghost fill needs the boundary conditions)");
    }
    m_in.jvp_scales.step = opt.fd_eps_scale;
    m_mem = krylov_memory(opt.gmres_restart, ba, dm);
    const double gb = 1.0e9;
    if (m_verbose > 0) {
        amrex::Print() << "gust: Newton-Krylov: FGMRES(" << opt.gmres_restart << ") Krylov memory estimate "
                       << m_mem.vectors << " vectors (" << m_mem.arnoldi_vectors << " Arnoldi + "
                       << opt.gmres_restart << " preconditioned + 2 work) x " << m_mem.vector_bytes / gb
                       << " GB = " << m_mem.total_bytes / gb << " GB in total, " << m_mem.max_rank_bytes / gb
                       << " GB on rank " << m_mem.max_rank << " (budget " << opt.mem_budget_gb << " GB per rank)\n";
    }
    if (m_mem.max_rank_bytes > static_cast<double>(opt.mem_budget_gb) * gb) {
        amrex::Abort("gust: steady.newton_mem_budget_gb = " + num(static_cast<double>(opt.mem_budget_gb))
                     + ": the Krylov memory estimate of rank " + std::to_string(m_mem.max_rank) + ", "
                     + num(m_mem.max_rank_bytes / gb) + " GB (steady.newton_gmres_restart = "
                     + std::to_string(opt.gmres_restart) + ": " + std::to_string(m_mem.vectors)
                     + " vectors), exceeds it; lower steady.newton_gmres_restart or use more ranks");
    }
}

void NewtonKrylov::base_residual (const State& q, ResidualFields& R0) const
{
    // R0 = R(fill(q)), the base point of every jvp (which refills q + eps dq).
    // On terrain with bc.zlo = surface_layer one fill is not a fixed point
    // (the surface state reads the ground w face, which the no-flux rule then
    // resets from the u, v ghosts), so differencing against R(q) would add
    // (R(fill(q)) - R(q)) / eps to every J v. A second fill of a copy gives
    // the same function on both sides of the difference. Where the fill is
    // idempotent (flat ground, the walls) R0 is R(q) bit for bit.
    State qf;
    qf.define(q.p.boxArray(), q.p.DistributionMap());
    auto to = qf.fields();
    const auto from = const_cast<State&>(q).fields();
    for (std::size_t n = 0; n < to.size(); ++n) {
        amrex::MultiFab::Copy(*to[n].mf, *from[n].mf, 0, 0, 1, to[n].mf->nGrowVect());
    }
    m_in.bc->fill_ghosts(qf, m_in.surface, m_in.buoyancy);
    residual(qf, m_in, m_mesh, m_geom, R0);
    amrex::Gpu::streamSynchronize();   // qf goes out of scope
}

NewtonStepInfo NewtonKrylov::solve (State& q, const amrex::MultiFab* dtau, NewtonPreconditioner* pc,
                                    const std::function<void (const NewtonOperator&, const NewtonVector&)>& use)
{
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy);
    ResidualFields R0;
    base_residual(q, R0);

    NewtonStepInfo info;
    {
        NewtonOperator op(q, R0, m_in, m_mesh, m_geom, dtau, pc);
        if (pc != nullptr) { pc->setup(q, dtau); }
        NewtonVector b = op.makeVecRHS();
        NewtonVector x = op.makeVecLHS();
        op.rhs(R0, b);
        info.residual_norm = op.norm2(b);

        amrex::GMRES<NewtonVector, NewtonOperator> gmres;
        gmres.setRestartLength(m_opt.gmres_restart);
        gmres.setMaxIters(m_opt.gmres_max_iters);
        gmres.setFlexible(true);
        gmres.define(op);
        gmres.solve(x, b, m_opt.krylov_rtol, amrex::Real(0));

        info.gmres_iterations = gmres.getNumIters();
        info.gmres_status = gmres.getStatus();
        info.krylov_residual = (info.residual_norm > 0) ? gmres.getResidualNorm() : amrex::Real(0);
        info.krylov_relative = (info.residual_norm > 0) ? info.krylov_residual / info.residual_norm : amrex::Real(0);
        info.eps_min = op.eps_min();
        info.eps_max = op.eps_max();
        info.jvp_calls = op.jvp_calls();

        // A Krylov residual above |b| (or not finite), or a non-finite |b|,
        // means x does worse than no step: q's valid points are left as they
        // are. The same values on every rank (the GMRES norms come from the
        // all-reduced inner products), so every rank takes the same branch.
        // A solve stopped by steady.newton_gmres_max_iters (status 1) goes
        // through the same rule.
        info.applied = std::isfinite(static_cast<double>(info.residual_norm))
                       && newton_update_usable(info.krylov_relative);
        if (info.applied) { use(op, x); }
        amrex::Gpu::streamSynchronize();   // the Krylov vectors go out of scope
    }
    return info;
}

void NewtonKrylov::report (const NewtonStepInfo& info) const
{
    if (m_verbose > 0 || !info.applied) {
        amrex::Print() << "gust: Newton step: |R| = " << info.residual_norm << ", GMRES "
                       << info.gmres_iterations << " iterations (status " << info.gmres_status
                       << "), Krylov residual " << info.krylov_residual << " (" << info.krylov_relative
                       << " relative), " << info.jvp_calls << " jvp calls, eps " << info.eps_min << " .. "
                       << info.eps_max << (info.applied ? "" : "; the update is not applied (|R| or the Krylov residual not finite, or the Krylov residual above |R|)")
                       << "\n";
    }
}

NewtonStepInfo NewtonKrylov::step (State& q, const amrex::MultiFab* dtau, NewtonPreconditioner* pc)
{
    NewtonStepInfo info = solve(q, dtau, pc, [&] (const NewtonOperator& op, const NewtonVector& x) { op.add_update(q, x); });
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy);
    report(info);
    m_last = info;
    return info;
}

NewtonStepInfo NewtonKrylov::direction (State& q, const amrex::MultiFab* dtau, State& dq, NewtonPreconditioner* pc)
{
    if (dq.p.boxArray() != q.p.boxArray() || dq.p.DistributionMap() != q.p.DistributionMap()) {
        amrex::Abort("gust: NewtonKrylov::direction: dq is not on q's layout");
    }
    for (int f = 0; f < n_newton_fields; ++f) { state_field(dq, f).setVal(0.0); }
    NewtonStepInfo info = solve(q, dtau, pc, [&] (const NewtonOperator& op, const NewtonVector& x) { op.add_update(dq, x); });
    report(info);
    m_last = info;
    return info;
}

amrex::Real NewtonKrylov::residual_norm (State& q)
{
    m_in.bc->fill_ghosts(q, m_in.surface, m_in.buoyancy);
    ResidualFields R0;
    base_residual(q, R0);
    NewtonOperator op(q, R0, m_in, m_mesh, m_geom, nullptr);
    NewtonVector b = op.makeVecRHS();
    op.rhs(R0, b);
    const amrex::Real r = op.norm2(b);
    amrex::Gpu::streamSynchronize();   // b goes out of scope
    return r;
}

} // namespace gust
