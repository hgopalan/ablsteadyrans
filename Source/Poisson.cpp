#include "Poisson.H"
#include "PoissonCross_K.H"

#include <ERF_TerrainPoisson.H>

#include <AMReX_GMRES.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>
#include <sstream>

namespace ablsteady {

namespace {

bool is_dirichlet_anywhere (const PoissonBC& bc)
{
    for (const auto f : bc.face) {
        if (f == PoissonFace::dirichlet) { return true; }
    }
    return false;
}

// J-weighted mean of the valid cells of mf.
amrex::Real weighted_mean (const amrex::MultiFab& mf, const amrex::MultiFab& J)
{
    return amrex::MultiFab::Dot(J, 0, mf, 0, 1, 0) / J.sum(0);
}

[[noreturn]] void abort_all (const std::vector<std::string>& err)
{
    std::string all = "ablsteady: invalid inputs:";
    for (const auto& e : err) { all += "\n  " + e; }
    amrex::Abort(all);
    std::abort();
}

// True if every valid value of mf is finite and > 0.
bool all_positive (const amrex::MultiFab& mf)
{
    const amrex::Real lo = mf.min(0);   // NaN compares false below
    bool ok = std::isfinite(lo) && lo > 0 && std::isfinite(mf.max(0));
    if (ok && mf.contains_nan(0, 1, 0)) { ok = false; }
    return ok;
}

} // namespace

// The linear operator for amrex::GMRES: ERF's TerrainPoisson, with
//  * the z_cross correction (beta_z - z_cross) cross on the z-faces, and
//  * the preconditioner applied to the residual divided by the cell's
//    horizontal coefficient (a constant factor would not change GMRES; a
//    varying one, as SIMPLEC's 1/a_P, it captures to first order).
// Everything else is TerrainPoisson's.
class PoissonSolver::Op
{
public:
    using RT = amrex::Real;   // amrex::GMRES's scalar type

    Op (TerrainPoisson& tp, const amrex::Geometry& geom, const Mesh& mesh, const amrex::MultiFab& z,
        const amrex::MultiFab* dcross, const amrex::MultiFab* beta_cell, const PoissonBC& bc)
        : m_tp(tp), m_geom(geom), m_mesh(mesh), m_z(z), m_dcross(dcross), m_beta_cell(beta_cell), m_bc(bc)
    {}

    void apply (amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
    {
        m_tp.apply(lhs, rhs);   // fills rhs's ghosts for the boundary conditions
        if (m_dcross != nullptr) { add_cross(lhs, rhs); }
    }

    // A phi += (1/J) d_zeta (a_z (beta_z - z_cross) cross) on the faces whose
    // flux the operator keeps (not the ground, not a Neumann top).
    void add_cross (amrex::MultiFab& lhs, amrex::MultiFab const& x)
    {
        const auto dxinv = m_geom.InvCellSizeArray();
        const int ktop = m_geom.Domain().bigEnd(2) + 1;
        const bool top_open = m_bc.face[amrex::Orientation(2, amrex::Orientation::high)] != PoissonFace::neumann;
        for (amrex::MFIter mfi(lhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto y = lhs.array(mfi);
            const auto xa = x.const_array(mfi);
            const auto zp = m_z.const_array(mfi);
            const auto az = m_mesh.az.const_array(mfi);
            const auto J = m_mesh.detJ.const_array(mfi);
            const auto d = m_dcross->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                auto face = [&] (int kf) -> amrex::Real {
                    if (kf == 0 || (kf == ktop && !top_open)) { return amrex::Real(0); }
                    return az(i, j, kf) * d(i, j, kf) * z_face_cross(i, j, kf, xa, zp, dxinv[0], dxinv[1]);
                };
                y(i, j, k) += (face(k + 1) - face(k)) * dxinv[2] / J(i, j, k);
            });
        }
    }

    // The z_cross part of the fluxes: F_z -= (beta_z - z_cross) cross.
    void add_cross_flux (amrex::MultiFab& fz, amrex::MultiFab const& x)
    {
        const auto dxinv = m_geom.InvCellSizeArray();
        const int ktop = m_geom.Domain().bigEnd(2) + 1;
        const bool top_open = m_bc.face[amrex::Orientation(2, amrex::Orientation::high)] != PoissonFace::neumann;
        for (amrex::MFIter mfi(fz, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto f = fz.array(mfi);
            const auto xa = x.const_array(mfi);
            const auto zp = m_z.const_array(mfi);
            const auto d = m_dcross->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == 0 || (k == ktop && !top_open)) { return; }
                f(i, j, k) -= d(i, j, k) * z_face_cross(i, j, k, xa, zp, dxinv[0], dxinv[1]);
            });
        }
    }

    void precond (amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
    {
        if (m_beta_cell == nullptr) {
            m_tp.precond(lhs, rhs);
            return;
        }
        amrex::MultiFab scaled(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(scaled, rhs, 0, 0, 1, 0);
        amrex::MultiFab::Divide(scaled, *m_beta_cell, 0, 0, 1, 0);
        m_tp.precond(lhs, scaled);
    }

    void assign (amrex::MultiFab& lhs, amrex::MultiFab const& rhs) { m_tp.assign(lhs, rhs); }
    void scale (amrex::MultiFab& lhs, amrex::Real fac) { m_tp.scale(lhs, fac); }
    amrex::Real dotProduct (amrex::MultiFab const& a, amrex::MultiFab const& b) { return m_tp.dotProduct(a, b); }
    void increment (amrex::MultiFab& lhs, amrex::MultiFab const& rhs, amrex::Real a) { m_tp.increment(lhs, rhs, a); }
    void linComb (amrex::MultiFab& lhs, amrex::Real a, amrex::MultiFab const& ra, amrex::Real b,
                  amrex::MultiFab const& rb)
    {
        m_tp.linComb(lhs, a, ra, b, rb);
    }
    amrex::MultiFab makeVecRHS () { return m_tp.makeVecRHS(); }
    amrex::MultiFab makeVecLHS () { return m_tp.makeVecLHS(); }
    amrex::Real norm2 (amrex::MultiFab const& v) { return m_tp.norm2(v); }
    void setToZero (amrex::MultiFab& v) { m_tp.setToZero(v); }

private:
    TerrainPoisson& m_tp;
    const amrex::Geometry& m_geom;
    const Mesh& m_mesh;
    const amrex::MultiFab& m_z;
    const amrex::MultiFab* m_dcross;
    const amrex::MultiFab* m_beta_cell;
    PoissonBC m_bc;
};

std::vector<std::string> check_poisson_bc (const PoissonBC& bc, const amrex::Geometry& geom)
{
    std::vector<std::string> err;
    const char* name[3] = {"x", "y", "z"};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (const auto side : {amrex::Orientation::low, amrex::Orientation::high}) {
            const PoissonFace f = bc.face[amrex::Orientation(d, side)];
            const std::string where = std::string(name[d]) + (side == amrex::Orientation::low ? "lo" : "hi");
            if (geom.isPeriodic(d) != (f == PoissonFace::periodic)) {
                err.emplace_back("Poisson BC " + where + ": periodic exactly when the domain is periodic in "
                                 + name[d]);
            }
        }
    }
    if (bc.face[amrex::Orientation(2, amrex::Orientation::low)] != PoissonFace::neumann) {
        err.emplace_back("Poisson BC zlo: the ground is a wall (Neumann)");
    }
    return err;
}

std::vector<std::string> check_poisson_options (const PoissonOptions& opt)
{
    std::vector<std::string> err;
    // GMRES cannot reduce the residual much below the rounding level.
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    if (!std::isfinite(opt.reltol) || !(opt.reltol >= 10 * eps) || !(opt.reltol < 1)) {
        std::ostringstream os;
        os << "solver.poisson_reltol must be >= " << 10 * eps << " (10 machine epsilon) and < 1";
        err.emplace_back(os.str());
    }
    if (!std::isfinite(opt.abstol) || !(opt.abstol >= 0)) {
        err.emplace_back("solver.poisson_abstol must be finite and >= 0");
    }
    if (opt.max_iter < 1) { err.emplace_back("solver.poisson_max_iter must be >= 1"); }
    if (opt.restart < 1) { err.emplace_back("solver.poisson_restart must be >= 1"); }
    if (opt.precond != "fft" && opt.precond != "none") {
        err.emplace_back("solver.poisson_precond must be fft or none (got " + opt.precond + ")");
    }
    if (opt.verbose < 0) { err.emplace_back("solver.poisson_verbose must be >= 0"); }
    return err;
}

PoissonOptions read_poisson_options ()
{
    PoissonOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("solver");
    for (const char* key : {"poisson_reltol", "poisson_abstol", "poisson_max_iter", "poisson_restart",
                            "poisson_precond", "poisson_verbose"}) {
        if (pp.countname(key) > 1) { err.emplace_back(std::string("solver.") + key + " is given more than once"); }
        if (pp.contains(key) && pp.countval(key) != 1) {
            err.emplace_back(std::string("solver.") + key + " takes exactly one value");
        }
    }
    pp.query("poisson_reltol", opt.reltol);
    pp.query("poisson_abstol", opt.abstol);
    pp.query("poisson_max_iter", opt.max_iter);
    pp.query("poisson_restart", opt.restart);
    pp.query("poisson_precond", opt.precond);
    pp.query("poisson_verbose", opt.verbose);
    auto range = check_poisson_options(opt);
    err.insert(err.end(), range.begin(), range.end());
    if (!err.empty()) { abort_all(err); }
    return opt;
}

PoissonSolver::PoissonSolver (const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd,
                              const PoissonOptions& opt)
    : m_mesh(mesh), m_geom(geom), m_z(z_phys_nd), m_opt(opt)
{
    auto err = check_poisson_options(opt);
    if (!err.empty()) { abort_all(err); }
    const auto& zl = mesh.zlevels();
    for (int k = 0; k + 1 < static_cast<int>(zl.size()); ++k) { m_dz.push_back(zl[k + 1] - zl[k]); }
}

PoissonSolver::~PoissonSolver () = default;

void PoissonSolver::check_layout (const amrex::MultiFab& mf, const char* what) const
{
    const amrex::BoxArray cc = amrex::convert(mf.boxArray(), amrex::IntVect(0));
    if (cc != m_mesh.detJ.boxArray() || mf.DistributionMap() != m_mesh.detJ.DistributionMap()) {
        amrex::Abort(std::string("ablsteady: PoissonSolver: ") + what
                     + " must be on the mesh's BoxArray and DistributionMapping");
    }
}

void PoissonSolver::prepare (const FaceCoeffs& coeffs, const PoissonBC& bc)
{
    auto err = check_poisson_bc(bc, m_geom);
    if (!err.empty()) { abort_all(err); }
    if (!std::isfinite(coeffs.precond_ratio) || !(coeffs.precond_ratio > 0)) {
        amrex::Abort("ablsteady: FaceCoeffs::precond_ratio must be finite and > 0");
    }
    const amrex::MultiFab* b[4] = {coeffs.x, coeffs.y, coeffs.z, coeffs.z_cross};
    const char* bname[4] = {"FaceCoeffs::x", "FaceCoeffs::y", "FaceCoeffs::z", "FaceCoeffs::z_cross"};
    const amrex::MultiFab* a[3] = {&m_mesh.ax, &m_mesh.ay, &m_mesh.az};
    for (int d = 0; d < 4; ++d) {
        if (b[d] == nullptr) { continue; }
        const amrex::MultiFab& ref = *a[d < 3 ? d : 2];
        if (b[d]->boxArray() != ref.boxArray() || b[d]->DistributionMap() != ref.DistributionMap()) {
            amrex::Abort(std::string("ablsteady: ") + bname[d] + " must be a face MultiFab on the mesh's box array");
        }
        if (!all_positive(*b[d])) {
            amrex::Abort(std::string("ablsteady: ") + bname[d] + " must be finite and > 0 everywhere");
        }
    }

    // A new mesh layout invalidates everything cached.
    const bool new_layout = !m_ax.ok() || m_op_ba != m_mesh.detJ.boxArray()
                            || m_ax.DistributionMap() != m_mesh.ax.DistributionMap();
    if (new_layout) {
        m_ax.define(m_mesh.ax.boxArray(), m_mesh.ax.DistributionMap(), 1, 0);
        m_ay.define(m_mesh.ay.boxArray(), m_mesh.ay.DistributionMap(), 1, 0);
        m_az.define(m_mesh.az.boxArray(), m_mesh.az.DistributionMap(), 1, 0);
        m_op.reset();
    }

    // Premultiplied areas a beta (beta = 1: an exact copy).
    amrex::MultiFab* ab[3] = {&m_ax, &m_ay, &m_az};
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab::Copy(*ab[d], *a[d], 0, 0, 1, 0);
        if (b[d] != nullptr) { amrex::MultiFab::Multiply(*ab[d], *b[d], 0, 0, 1, 0); }
    }

    // z_cross: the z-face correction beta_z - z_cross.
    m_split = coeffs.z_cross != nullptr;
    if (m_split) {
        if (!m_dcross.ok() || new_layout) { m_dcross.define(m_mesh.az.boxArray(), m_mesh.az.DistributionMap(), 1, 0); }
        if (coeffs.z != nullptr) {
            amrex::MultiFab::Copy(m_dcross, *coeffs.z, 0, 0, 1, 0);
        } else {
            m_dcross.setVal(1.0);
        }
        amrex::MultiFab::Subtract(m_dcross, *coeffs.z_cross, 0, 0, 1, 0);
    }

    // Cell horizontal coefficient for the preconditioner: the mean of the
    // cell's four x- and y-face values.
    m_scaled_precond = coeffs.x != nullptr || coeffs.y != nullptr;
    if (m_scaled_precond) {
        if (!m_beta_cell.ok() || new_layout) {
            m_beta_cell.define(m_mesh.detJ.boxArray(), m_mesh.detJ.DistributionMap(), 1, 0);
        }
        for (amrex::MFIter mfi(m_beta_cell, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto c = m_beta_cell.array(mfi);
            const bool hx = coeffs.x != nullptr, hy = coeffs.y != nullptr;
            const auto bxa = hx ? coeffs.x->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            const auto bya = hy ? coeffs.y->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real sx = hx ? bxa(i, j, k) + bxa(i + 1, j, k) : amrex::Real(2);
                const amrex::Real sy = hy ? bya(i, j, k) + bya(i, j + 1, k) : amrex::Real(2);
                c(i, j, k) = amrex::Real(0.25) * (sx + sy);
            });
        }
    }

    if (m_op && m_op_bc == bc && m_op_ratio == coeffs.precond_ratio) { return; }

    // ERF's boundary names: Outflow/Open map to a Dirichlet (odd) face, any
    // other name to a Neumann (even) face; periodic comes from the geometry.
    amrex::Array<std::string, 2 * AMREX_SPACEDIM> names;
    for (int f = 0; f < 2 * AMREX_SPACEDIM; ++f) {
        names[f] = (bc.face[f] == PoissonFace::dirichlet) ? "Open"
                 : (bc.face[f] == PoissonFace::periodic)  ? "Periodic" : "SlipWall";
    }
    // Preconditioner: the flat Laplacian with dz' = dz sqrt(beta_h / beta_v).
    const amrex::Real scale = std::sqrt(amrex::Real(1) / coeffs.precond_ratio);
    amrex::Vector<amrex::Real> dz_h(m_dz.size());
    for (int k = 0; k < static_cast<int>(m_dz.size()); ++k) { dz_h[k] = m_dz[k] * scale; }
    amrex::Gpu::DeviceVector<amrex::Real> dz_d(dz_h.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dz_h.begin(), dz_h.end(), dz_d.begin());

    m_op.reset();
    m_op = std::make_unique<TerrainPoisson>(m_geom, m_geom, m_mesh.detJ.boxArray(), m_mesh.detJ.DistributionMap(),
                                            names, dz_d, m_ax, m_ay, m_az, m_mesh.detJ, &m_z, false);
    amrex::Gpu::streamSynchronize();   // dz_d is copied asynchronously
    m_op_bc = bc;
    m_op_ratio = coeffs.precond_ratio;
    m_op_ba = m_mesh.detJ.boxArray();
}

void PoissonSolver::apply (amrex::MultiFab& Aphi, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                           const PoissonBC& bc)
{
    if (phi.nGrowVect().min() < 1) { amrex::Abort("ablsteady: PoissonSolver::apply needs phi with a ghost cell"); }
    check_layout(phi, "phi");
    check_layout(Aphi, "Aphi");
    prepare(coeffs, bc);
    Op op(*m_op, m_geom, m_mesh, m_z, m_split ? &m_dcross : nullptr, nullptr, bc);
    op.apply(Aphi, phi);
}

PoissonResult PoissonSolver::solve (amrex::MultiFab& rhs, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                                    const PoissonBC& bc, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* fluxes)
{
    if (phi.nGrowVect().min() < 1) { amrex::Abort("ablsteady: PoissonSolver::solve needs phi with a ghost cell"); }
    check_layout(rhs, "rhs");
    check_layout(phi, "phi");
    if (fluxes != nullptr) {
        const amrex::MultiFab* a[3] = {&m_mesh.ax, &m_mesh.ay, &m_mesh.az};
        for (int d = 0; d < 3; ++d) {
            const auto& f = (*fluxes)[d];
            // TerrainPoisson::getFluxes writes every point of its fluxes.
            if (f.nGrowVect().max() != 0 || f.boxArray() != a[d]->boxArray()
                || f.DistributionMap() != a[d]->DistributionMap()) {
                amrex::Abort("ablsteady: PoissonSolver::solve: fluxes must be face MultiFabs on the mesh's box "
                             "array, without ghost cells");
            }
        }
    }
    prepare(coeffs, bc);
    Op op(*m_op, m_geom, m_mesh, m_z, m_split ? &m_dcross : nullptr, m_scaled_precond ? &m_beta_cell : nullptr,
          bc);

    PoissonResult res;
    res.singular = !is_dirichlet_anywhere(bc);
    if (res.singular) {
        // Only the mean-free part of rhs is in the range of A (J is its left
        // null vector).
        res.removed_mean = weighted_mean(rhs, m_mesh.detJ);
        rhs.plus(-res.removed_mean, 0, 1, 0);
    }

    amrex::GMRES<amrex::MultiFab, Op> gmres;
    gmres.define(op);
    gmres.setVerbose(m_opt.verbose);
    gmres.setRestartLength(m_opt.restart);
    gmres.setMaxIters(m_opt.max_iter);
    // The preconditioner is a fixed linear map, so standard (not flexible)
    // GMRES applies; flexible GMRES would keep restart more vectors.
    gmres.setFlexible(false);
    m_op->usePrecond(m_opt.precond == "fft");
    gmres.solve(phi, rhs, m_opt.reltol, m_opt.abstol, m_opt.max_iter);

    res.iterations = gmres.getNumIters();
    res.initial_residual = gmres.getInitialResidualNorm();
    if (res.singular) {
        phi.plus(-weighted_mean(phi, m_mesh.detJ), 0, 1, 0);
    }
    // The true residual, not GMRES's estimate.
    amrex::MultiFab r(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
    op.apply(r, phi);
    amrex::MultiFab::Xpay(r, amrex::Real(-1), rhs, 0, 0, 1, 0);
    res.final_residual = r.norm2();

    if (gmres.getStatus() != 0) {
        std::ostringstream os;
        os << "ablsteady: the Poisson solve did not converge in solver.poisson_max_iter = " << m_opt.max_iter
           << " iterations (residual " << res.final_residual << " of " << res.initial_residual
           << "; solver.poisson_reltol = " << m_opt.reltol << ", solver.poisson_abstol = " << m_opt.abstol << ")";
        amrex::Abort(os.str());
    }

    m_op->apply_bcs(phi);
    if (fluxes != nullptr) {
        m_op->getFluxes(phi, *fluxes);
        const amrex::MultiFab* b[3] = {coeffs.x, coeffs.y, coeffs.z};
        for (int d = 0; d < 3; ++d) {
            if (b[d] != nullptr) { amrex::MultiFab::Multiply((*fluxes)[d], *b[d], 0, 0, 1, 0); }
        }
        if (m_split) { op.add_cross_flux((*fluxes)[2], phi); }
    }
    if (m_opt.verbose > 0) {
        amrex::Print() << "ablsteady: Poisson solve: " << res.iterations << " iterations, residual "
                       << res.final_residual << " of " << res.initial_residual;
        if (res.singular) { amrex::Print() << ", removed mean " << res.removed_mean; }
        amrex::Print() << "\n";
    }
    return res;
}

void face_divergence (const Mesh& mesh, const amrex::Geometry& geom,
                      const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, amrex::MultiFab& div)
{
    const amrex::MultiFab* a[3] = {&mesh.ax, &mesh.ay, &mesh.az};
    for (int d = 0; d < 3; ++d) {
        if (F[d]->boxArray() != a[d]->boxArray() || F[d]->DistributionMap() != a[d]->DistributionMap()) {
            amrex::Abort("ablsteady: face_divergence: F must be face MultiFabs on the mesh's box array");
        }
    }
    if (div.boxArray() != mesh.detJ.boxArray() || div.DistributionMap() != mesh.detJ.DistributionMap()) {
        amrex::Abort("ablsteady: face_divergence: div must be on the mesh's box array");
    }
    const auto dxinv = geom.InvCellSizeArray();
    for (amrex::MFIter mfi(div, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto d = div.array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        const auto ax = mesh.ax.const_array(mfi);
        const auto ay = mesh.ay.const_array(mfi);
        const auto az = mesh.az.const_array(mfi);
        const auto fx = F[0]->const_array(mfi);
        const auto fy = F[1]->const_array(mfi);
        const auto fz = F[2]->const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i, j, k) = ((ax(i + 1, j, k) * fx(i + 1, j, k) - ax(i, j, k) * fx(i, j, k)) * dxinv[0]
                        + (ay(i, j + 1, k) * fy(i, j + 1, k) - ay(i, j, k) * fy(i, j, k)) * dxinv[1]
                        + (az(i, j, k + 1) * fz(i, j, k + 1) - az(i, j, k) * fz(i, j, k)) * dxinv[2])
                       / J(i, j, k);
        });
    }
}

} // namespace ablsteady
