#include "Poisson.H"
#include "PoissonCross_K.H"

#include <ERF_TerrainPoisson.H>

#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTerrainPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <sstream>

namespace gust {

namespace {

bool op_dir_ground (const PoissonBC& bc)
{
    return bc.face[amrex::Orientation(2, amrex::Orientation::low)] == PoissonFace::dirichlet;
}

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
    std::string all = "gust: invalid inputs:";
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

// The cell's horizontal coefficient: the mean of its four x- and y-face
// values (1 where a coefficient is not given). A free function: CUDA does not
// allow device lambdas in private member functions such as prepare().
void fill_beta_cell (amrex::MultiFab& cell, const amrex::MultiFab* bx_face, const amrex::MultiFab* by_face)
{
    for (amrex::MFIter mfi(cell, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto c = cell.array(mfi);
        const bool hx = bx_face != nullptr, hy = by_face != nullptr;
        const auto bxa = hx ? bx_face->const_array(mfi) : amrex::Array4<amrex::Real const>{};
        const auto bya = hy ? by_face->const_array(mfi) : amrex::Array4<amrex::Real const>{};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real sx = hx ? bxa(i, j, k) + bxa(i + 1, j, k) : amrex::Real(2);
            const amrex::Real sy = hy ? bya(i, j, k) + bya(i, j + 1, k) : amrex::Real(2);
            c(i, j, k) = amrex::Real(0.25) * (sx + sy);
        });
    }
}

// Fill phi's ghost cells for bc: periodic images, then even (Neumann) or odd
// (Dirichlet) copies across the non-periodic faces, x first, then y over the
// x ghosts, then z over both (TerrainPoisson::apply_bcs's order, with the
// ground odd when it is Dirichlet). A free function: CUDA does not allow
// device lambdas in private member functions.
void fill_bc_ghosts (amrex::MultiFab& phi, const amrex::Geometry& geom, const PoissonBC& bc)
{
    phi.FillBoundary(geom.periodicity());
    const amrex::Box& dom = geom.Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) { continue; }
        for (const auto side : {amrex::Orientation::low, amrex::Orientation::high}) {
            const PoissonFace f = bc.face[amrex::Orientation(d, side)];
            const amrex::Real sign = (f == PoissonFace::dirichlet) ? amrex::Real(-1) : amrex::Real(1);
            const bool lo = side == amrex::Orientation::low;
            const int edge = lo ? dom.smallEnd(d) : dom.bigEnd(d);
            for (amrex::MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                if ((lo && vb.smallEnd(d) != edge) || (!lo && vb.bigEnd(d) != edge)) { continue; }
                // The ghost slab next to the face, grown over the other
                // directions: fully over those already filled (e < d), and
                // over the later ones (e > d) only where the neighbours hold
                // values now (inside the domain or across a periodic face),
                // so every corner ends up filled.
                amrex::Box gb = vb;
                for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                    if (e == d) { continue; }
                    gb.grow(e, 1);
                    if (e > d && !geom.isPeriodic(e)) {
                        gb.setSmall(e, std::max(gb.smallEnd(e), dom.smallEnd(e)));
                        gb.setBig(e, std::min(gb.bigEnd(e), dom.bigEnd(e)));
                    }
                }
                gb.setSmall(d, lo ? edge - 1 : edge + 1);
                gb.setBig(d, lo ? edge - 1 : edge + 1);
                gb &= mfi.fabbox();
                const auto p = phi.array(mfi);
                const int s = lo ? 1 : -1;
                const int di = (d == 0) ? s : 0, dj = (d == 1) ? s : 0, dk = (d == 2) ? s : 0;
                amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    p(i, j, k) = sign * p(i + di, j + dj, k + dk);
                });
            }
        }
    }
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
        const amrex::MultiFab* dcross, const amrex::MultiFab* beta_cell, const PoissonBC& bc,
        const amrex::MultiFab& ax, const amrex::MultiFab& ay, const amrex::MultiFab& az)
        : m_tp(tp), m_geom(geom), m_mesh(mesh), m_z(z), m_dcross(dcross), m_beta_cell(beta_cell), m_bc(bc),
          m_ax(ax), m_ay(ay), m_az(az),
          m_dir_ground(bc.face[amrex::Orientation(2, amrex::Orientation::low)] == PoissonFace::dirichlet)
    {}

    void apply (amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
    {
        if (m_dir_ground) {
            apply_dirichlet_ground(lhs, rhs);
        } else {
            m_tp.apply(lhs, rhs);   // fills rhs's ghosts for the boundary conditions
        }
        if (m_dcross != nullptr) { add_cross(lhs, rhs); }
    }

    // A phi with phi = 0 on the ground: ERF's kernel with an odd ghost below
    // the ground (so the bottom row's cross terms see the Dirichlet value),
    // plus the ground-face flux ERF's kernel leaves out, from ERF's own
    // one-sided Dirichlet formula (terrpoisson_flux_zlo_dir).
    void apply_dirichlet_ground (amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
    {
        auto& x = const_cast<amrex::MultiFab&>(rhs);
        fill_bc_ghosts(x, m_geom, m_bc);
        const auto dxinv = m_geom.InvCellSizeArray();
        const int zf = m_tp.zeroFluxFaces();
        const amrex::Dim3 dlo = amrex::lbound(m_geom.Domain());
        const amrex::Dim3 dhi = amrex::ubound(m_geom.Domain());
        const int k0 = m_geom.Domain().smallEnd(2);
        for (amrex::MFIter mfi(lhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto y = lhs.array(mfi);
            const auto xa = x.const_array(mfi);
            const auto ax = m_ax.const_array(mfi);
            const auto ay = m_ay.const_array(mfi);
            const auto az = m_az.const_array(mfi);
            const auto J = m_mesh.detJ.const_array(mfi);
            const auto zp = m_z.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                terrpoisson_adotx(i, j, k, y, xa, ax, ay, az, J, zp, dxinv[0], dxinv[1], dxinv[2], dlo, dhi, zf);
                if (k == k0) {
                    y(i, j, k) += az(i, j, k) * terrpoisson_flux_zlo_dir(i, j, k, xa, zp, dxinv[0], dxinv[1])
                                * dxinv[2] / J(i, j, k);
                }
            });
        }
    }

    // The fluxes with phi = 0 on the ground (as apply_dirichlet_ground).
    void fluxes_dirichlet_ground (amrex::MultiFab& phi, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fl)
    {
        fill_bc_ghosts(phi, m_geom, m_bc);
        const auto dxinv = m_geom.InvCellSizeArray();
        const int zf = m_tp.zeroFluxFaces();
        const amrex::Dim3 dlo = amrex::lbound(m_geom.Domain());
        const amrex::Dim3 dhi = amrex::ubound(m_geom.Domain());
        const int k0 = m_geom.Domain().smallEnd(2);
        for (amrex::MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
            const auto xa = phi.const_array(mfi);
            const auto zp = m_z.const_array(mfi);
            const auto fx = fl[0].array(mfi);
            const auto fy = fl[1].array(mfi);
            const auto fz = fl[2].array(mfi);
            amrex::ParallelFor(fl[0][mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const bool zero = ((zf & tp_face_xlo) && i == dlo.x) || ((zf & tp_face_xhi) && i == dhi.x + 1);
                fx(i, j, k) = zero ? amrex::Real(0) : terrpoisson_flux_x(i, j, k, xa, zp, dxinv[0]);
            });
            amrex::ParallelFor(fl[1][mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const bool zero = ((zf & tp_face_ylo) && j == dlo.y) || ((zf & tp_face_yhi) && j == dhi.y + 1);
                fy(i, j, k) = zero ? amrex::Real(0) : terrpoisson_flux_y(i, j, k, xa, zp, dxinv[1]);
            });
            amrex::ParallelFor(fl[2][mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if ((zf & tp_face_zhi) && k == dhi.z + 1) {
                    fz(i, j, k) = amrex::Real(0);
                } else if (k == k0) {
                    fz(i, j, k) = terrpoisson_flux_zlo_dir(i, j, k, xa, zp, dxinv[0], dxinv[1]);
                } else {
                    fz(i, j, k) = terrpoisson_flux_z(i, j, k, xa, zp, dxinv[0], dxinv[1]);
                }
            });
        }
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
    const amrex::MultiFab& m_ax;
    const amrex::MultiFab& m_ay;
    const amrex::MultiFab& m_az;
    bool m_dir_ground;
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
    if (opt.backend != "gmres_fft" && opt.backend != "mlmg" && opt.backend != "gmres_mlmg") {
        err.emplace_back("solver.poisson_backend must be gmres_fft, mlmg or gmres_mlmg (got " + opt.backend + ")");
    }
    auto given = [&opt] (const char* k) { return std::find(opt.given.begin(), opt.given.end(), k) != opt.given.end(); };
    if (given("poisson_precond") && opt.backend != "gmres_fft") {
        err.emplace_back("solver.poisson_precond has no effect with solver.poisson_backend = " + opt.backend
                         + " (back end A only)"
                         + (given("poisson_backend") ? std::string()
                                                     : std::string("; gmres_mlmg is the default: remove the key, or set "
                                                                   "solver.poisson_backend = gmres_fft")));
    }
    if (given("poisson_restart") && opt.backend == "mlmg") {
        err.emplace_back("solver.poisson_restart has no effect with solver.poisson_backend = mlmg (no GMRES)");
    }
    return err;
}

PoissonOptions read_poisson_options ()
{
    PoissonOptions opt;
    std::vector<std::string> err;
    amrex::ParmParse pp("solver");
    for (const char* key : {"poisson_reltol", "poisson_abstol", "poisson_max_iter", "poisson_restart",
                            "poisson_precond", "poisson_verbose", "poisson_backend"}) {
        if (pp.countname(key) > 1) { err.emplace_back(std::string("solver.") + key + " is given more than once"); }
        if (pp.contains(key) && pp.countval(key) != 1) {
            err.emplace_back(std::string("solver.") + key + " takes exactly one value");
        }
        if (pp.contains(key)) { opt.given.emplace_back(key); }
    }
    pp.query("poisson_reltol", opt.reltol);
    pp.query("poisson_abstol", opt.abstol);
    pp.query("poisson_max_iter", opt.max_iter);
    pp.query("poisson_restart", opt.restart);
    pp.query("poisson_precond", opt.precond);
    pp.query("poisson_verbose", opt.verbose);
    pp.query("poisson_backend", opt.backend);
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
        amrex::Abort(std::string("gust: PoissonSolver: ") + what
                     + " must be on the mesh's BoxArray and DistributionMapping");
    }
    // The ghost fills assume boxes that span the z column (make_box_array).
    for (int b = 0; b < cc.size(); ++b) {
        if (cc[b].smallEnd(2) != m_geom.Domain().smallEnd(2) || cc[b].bigEnd(2) != m_geom.Domain().bigEnd(2)) {
            amrex::Abort("gust: PoissonSolver: every box must span the full z column");
        }
    }
}

void PoissonSolver::prepare (const FaceCoeffs& coeffs, const PoissonBC& bc)
{
    auto err = check_poisson_bc(bc, m_geom);
    if (!err.empty()) { abort_all(err); }
    if (!std::isfinite(coeffs.precond_ratio) || !(coeffs.precond_ratio > 0)) {
        amrex::Abort("gust: FaceCoeffs::precond_ratio must be finite and > 0");
    }
    const amrex::MultiFab* b[4] = {coeffs.x, coeffs.y, coeffs.z, coeffs.z_cross};
    const char* bname[4] = {"FaceCoeffs::x", "FaceCoeffs::y", "FaceCoeffs::z", "FaceCoeffs::z_cross"};
    const amrex::MultiFab* a[3] = {&m_mesh.ax, &m_mesh.ay, &m_mesh.az};
    for (int d = 0; d < 4; ++d) {
        if (b[d] == nullptr) { continue; }
        const amrex::MultiFab& ref = *a[d < 3 ? d : 2];
        if (b[d]->boxArray() != ref.boxArray() || b[d]->DistributionMap() != ref.DistributionMap()) {
            amrex::Abort(std::string("gust: ") + bname[d] + " must be a face MultiFab on the mesh's box array");
        }
        if (!all_positive(*b[d])) {
            amrex::Abort(std::string("gust: ") + bname[d] + " must be finite and > 0 everywhere");
        }
    }

    // A new mesh layout invalidates everything cached.
    const bool new_layout = !m_ax.ok() || m_layout_ba != m_mesh.detJ.boxArray()
                            || m_layout_dm != m_mesh.detJ.DistributionMap();
    if (new_layout) {
        m_layout_ba = m_mesh.detJ.boxArray();
        m_layout_dm = m_mesh.detJ.DistributionMap();
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

    // z_cross: the z-face correction beta_z - z_cross (back end A, Neumann
    // ground; its kernel has no Dirichlet ground face).
    m_split = coeffs.z_cross != nullptr;
    if (m_split && m_opt.backend != "gmres_fft") {
        amrex::Abort("gust: FaceCoeffs::z_cross needs solver.poisson_backend = gmres_fft "
                     "(MLTerrainPoisson has no split coefficient)");
    }
    if (m_split && bc.face[amrex::Orientation(2, amrex::Orientation::low)] == PoissonFace::dirichlet) {
        amrex::Abort("gust: FaceCoeffs::z_cross with a Dirichlet ground is not supported");
    }
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
        fill_beta_cell(m_beta_cell, coeffs.x, coeffs.y);
    }

    if (m_opt.backend != "gmres_fft") { return; }   // back end B builds its own operator
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
    if (phi.nGrowVect().min() < 1) { amrex::Abort("gust: PoissonSolver::apply needs phi with a ghost cell"); }
    check_layout(phi, "phi");
    check_layout(Aphi, "Aphi");
    if (m_opt.backend != "gmres_fft") {
        amrex::Abort("gust: PoissonSolver::apply is back end A's operator (solver.poisson_backend = gmres_fft)");
    }
    prepare(coeffs, bc);
    Op op(*m_op, m_geom, m_mesh, m_z, m_split ? &m_dcross : nullptr, nullptr, bc, m_ax, m_ay, m_az);
    op.apply(Aphi, phi);
}

void PoissonSolver::fill_ghosts (amrex::MultiFab& phi, const PoissonBC& bc)
{
    if (phi.nGrowVect().min() < 1) { amrex::Abort("gust: PoissonSolver::fill_ghosts needs phi with a ghost cell"); }
    check_layout(phi, "phi");
    auto err = check_poisson_bc(bc, m_geom);
    if (!err.empty()) { abort_all(err); }
    fill_bc_ghosts(phi, m_geom, bc);
}

PoissonResult PoissonSolver::solve (amrex::MultiFab& rhs, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                                    const PoissonBC& bc, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* fluxes)
{
    if (phi.nGrowVect().min() < 1) { amrex::Abort("gust: PoissonSolver::solve needs phi with a ghost cell"); }
    check_layout(rhs, "rhs");
    check_layout(phi, "phi");
    if (fluxes != nullptr) {
        const amrex::MultiFab* a[3] = {&m_mesh.ax, &m_mesh.ay, &m_mesh.az};
        for (int d = 0; d < 3; ++d) {
            const auto& f = (*fluxes)[d];
            // TerrainPoisson::getFluxes writes every point of its fluxes.
            if (f.nGrowVect().max() != 0 || f.boxArray() != a[d]->boxArray()
                || f.DistributionMap() != a[d]->DistributionMap()) {
                amrex::Abort("gust: PoissonSolver::solve: fluxes must be face MultiFabs on the mesh's box "
                             "array, without ghost cells");
            }
        }
    }
    prepare(coeffs, bc);

    PoissonResult res;
    res.singular = !is_dirichlet_anywhere(bc);
    if (res.singular) {
        // Only the mean-free part of rhs is in the range of A (J is its left
        // null vector).
        res.removed_mean = weighted_mean(rhs, m_mesh.detJ);
        rhs.plus(-res.removed_mean, 0, 1, 0);
    }
    if (m_opt.backend != "gmres_fft") { return solve_mlmg(rhs, phi, coeffs, bc, fluxes, res); }

    Op op(*m_op, m_geom, m_mesh, m_z, m_split ? &m_dcross : nullptr, m_scaled_precond ? &m_beta_cell : nullptr,
          bc, m_ax, m_ay, m_az);

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
        os << "gust: the Poisson solve did not converge in solver.poisson_max_iter = " << m_opt.max_iter
           << " iterations (residual " << res.final_residual << " of " << res.initial_residual
           << "; solver.poisson_reltol = " << m_opt.reltol << ", solver.poisson_abstol = " << m_opt.abstol << ")";
        amrex::Abort(os.str());
    }

    fill_bc_ghosts(phi, m_geom, bc);
    if (fluxes != nullptr) {
        if (op_dir_ground(bc)) {
            op.fluxes_dirichlet_ground(phi, *fluxes);
        } else {
            m_op->getFluxes(phi, *fluxes);
        }
        const amrex::MultiFab* b[3] = {coeffs.x, coeffs.y, coeffs.z};
        for (int d = 0; d < 3; ++d) {
            if (b[d] != nullptr) { amrex::MultiFab::Multiply((*fluxes)[d], *b[d], 0, 0, 1, 0); }
        }
        if (m_split) { op.add_cross_flux((*fluxes)[2], phi); }
    }
    if (m_opt.verbose > 0) {
        amrex::Print() << "gust: Poisson solve: " << res.iterations << " iterations, residual "
                       << res.final_residual << " of " << res.initial_residual;
        if (res.singular) { amrex::Print() << ", removed mean " << res.removed_mean; }
        amrex::Print() << "\n";
    }
    return res;
}

PoissonResult PoissonSolver::solve_mlmg (amrex::MultiFab& rhs, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                                         const PoissonBC& bc,
                                         amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* fluxes, PoissonResult res)
{
    const amrex::BoxArray& ba = m_mesh.detJ.boxArray();
    const amrex::DistributionMapping& dm = m_mesh.detJ.DistributionMap();
    if (!m_ml || !(m_ml_bc == bc) || m_ml_ba != ba || m_ml_dm != dm) {
        amrex::LPInfo info;
        // A periodic one-cell direction (a 2-D case) is hidden: nothing can
        // vary along it. (A non-periodic one-cell direction keeps its
        // terrain slopes, as back end A does.)
        for (int d = 0; d < 2; ++d) {
            if (m_geom.Domain().length(d) == 1 && m_geom.isPeriodic(d)) { info.setHiddenDirection(d); }
        }
        m_ml = std::make_unique<amrex::MLTerrainPoisson>(amrex::Vector<amrex::Geometry>{m_geom},
                                                         amrex::Vector<amrex::BoxArray>{ba},
                                                         amrex::Vector<amrex::DistributionMapping>{dm}, info);
        amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lo{}, hi{};
        auto type = [] (PoissonFace f) {
            return f == PoissonFace::periodic  ? amrex::LinOpBCType::Periodic
                 : f == PoissonFace::dirichlet ? amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
        };
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lo[d] = type(bc.face[amrex::Orientation(d, amrex::Orientation::low)]);
            hi[d] = type(bc.face[amrex::Orientation(d, amrex::Orientation::high)]);
        }
        m_ml->setDomainBC(lo, hi);
        m_ml->setLevelBC(0, nullptr);
        m_ml->setZPhys(0, m_z);
        m_ml->setDetJ(0, m_mesh.detJ);
        m_ml_bc = bc;
        m_ml_ba = ba;
        m_ml_dm = dm;
    }
    m_ml->setAreas(0, {&m_ax, &m_ay, &m_az});   // premultiplied a beta

    amrex::MLMG mlmg(*m_ml);
    mlmg.setMaxIter(m_opt.max_iter);
    mlmg.setVerbose(m_opt.verbose);
    mlmg.setThrowException(true);
    res.initial_residual = rhs.norm2();
    phi.setVal(0.0);
    bool converged = true;
    std::string why;
    // MLTerrainPoisson's CPU red-black smoother (AMReX
    // AMReX_MLTerrainPoisson.cpp, the fused loop's sa += ta) adds whole vector
    // registers and stores only the active colour; the inactive lanes are
    // uninitialized and can raise FE_INVALID with no effect on the result.
    // The invalid trap is off during the solve; the result is checked for
    // finite values below instead.
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    try {
        if (m_opt.backend == "mlmg") {
            mlmg.solve({&phi}, {&rhs}, m_opt.reltol, m_opt.abstol);
            res.iterations = mlmg.getNumIters();
        } else {
            amrex::GMRESMLMG gmres(mlmg);
            gmres.setMaxIters(m_opt.max_iter);
            gmres.setVerbose(m_opt.verbose);
            gmres.getGMRES().setRestartLength(m_opt.restart);
            gmres.solve(phi, rhs, m_opt.reltol, m_opt.abstol);
            res.iterations = gmres.getNumIters();
            converged = gmres.getGMRES().getStatus() == 0;
        }
    } catch (const std::exception& e) {
        converged = false;
        why = e.what();
        res.iterations = m_opt.max_iter;
    }
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev_fpe);
    if (phi.contains_nan(0, 1, 0) || phi.contains_inf(0, 1, 0)) {
        amrex::Abort("gust: the Poisson solve (solver.poisson_backend = " + m_opt.backend
                     + ") returned non-finite values");
    }
    if (res.singular) { phi.plus(-weighted_mean(phi, m_mesh.detJ), 0, 1, 0); }

    // The true residual: MLTerrainPoisson applies J A phi.
    amrex::MultiFab r(ba, dm, 1, 0);
    mlmg.apply({&r}, {&phi});
    amrex::MultiFab::Divide(r, m_mesh.detJ, 0, 0, 1, 0);
    amrex::MultiFab::Xpay(r, amrex::Real(-1), rhs, 0, 0, 1, 0);
    res.final_residual = r.norm2();
    if (!converged) {
        std::ostringstream os;
        os << "gust: the Poisson solve (solver.poisson_backend = " << m_opt.backend
           << ") did not converge in solver.poisson_max_iter = " << m_opt.max_iter << " iterations (residual "
           << res.final_residual << " of " << res.initial_residual << "; solver.poisson_reltol = " << m_opt.reltol
           << ", solver.poisson_abstol = " << m_opt.abstol << ")";
        if (!why.empty()) { os << ": " << why; }
        amrex::Abort(os.str());
    }

    fill_bc_ghosts(phi, m_geom, bc);
    if (fluxes != nullptr) {
        m_ml->getFluxes({{&(*fluxes)[0], &(*fluxes)[1], &(*fluxes)[2]}}, {&phi}, amrex::MLLinOp::Location::FaceCenter);
        const amrex::MultiFab* b[3] = {coeffs.x, coeffs.y, coeffs.z};
        for (int d = 0; d < 3; ++d) {
            if (b[d] != nullptr) { amrex::MultiFab::Multiply((*fluxes)[d], *b[d], 0, 0, 1, 0); }
        }
    }
    if (m_opt.verbose > 0) {
        amrex::Print() << "gust: Poisson solve (" << m_opt.backend << "): " << res.iterations
                       << " iterations, residual " << res.final_residual << " of " << res.initial_residual << "\n";
    }
    return res;
}

void face_divergence (const Mesh& mesh, const amrex::Geometry& geom,
                      const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, amrex::MultiFab& div)
{
    const amrex::MultiFab* a[3] = {&mesh.ax, &mesh.ay, &mesh.az};
    for (int d = 0; d < 3; ++d) {
        if (F[d]->boxArray() != a[d]->boxArray() || F[d]->DistributionMap() != a[d]->DistributionMap()) {
            amrex::Abort("gust: face_divergence: F must be face MultiFabs on the mesh's box array");
        }
    }
    if (div.boxArray() != mesh.detJ.boxArray() || div.DistributionMap() != mesh.detJ.DistributionMap()) {
        amrex::Abort("gust: face_divergence: div must be on the mesh's box array");
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

} // namespace gust
