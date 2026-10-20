#include "Transfer.H"

#include "CVFluxes.H"
#include "MassFlux.H"
#include "Residual.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <string>

namespace ablsteady {

namespace {

//! Where a fine point sits among the coarse points along one direction:
//! i_c = clamp(floor(m / den), imin, imax) and t = (m - den i_c) / den, with
//! m = 2 (i - flo) + 1 - r, den = 2r for a cell-centred direction (the fine
//! centre at (i_c + 1/2 + t) coarse cells) and m = i - flo, den = r with
//! floor((m - 1) / den) for a nodal one. The clamp bounds give the linear
//! extrapolation at a non-periodic face (t < 0 or t > 1); a periodic
//! direction reads the periodic images at i_c = -1 and N. lin = 0 and
//! step = 0 for a non-periodic cell-centred direction one coarse cell wide
//! (constant there: no second point to extrapolate with).
struct Axis
{
    int flo{0}, clo{0};
    int r{2};
    int imin{0}, imax{0};
    int step{1};
    amrex::Real lin{1};
    amrex::Real inv{0};   // 1/den, a power of two
};

Axis make_axis (const amrex::Geometry& gf, const amrex::Geometry& gc, int d, int r, bool nodal)
{
    Axis a;
    a.flo = gf.Domain().smallEnd(d);
    a.clo = gc.Domain().smallEnd(d);
    a.r = r;
    const int N = gc.Domain().length(d);
    if (nodal) {
        a.imin = 0;
        a.imax = N - 1;
        a.inv = amrex::Real(1) / amrex::Real(r);
    } else {
        a.inv = amrex::Real(1) / amrex::Real(2 * r);
        if (gc.isPeriodic(d)) {
            a.imin = -1;
            a.imax = N - 1;
        } else if (N >= 2) {
            a.imin = 0;
            a.imax = N - 2;
        } else {
            a.imin = 0;
            a.imax = 0;
            a.step = 0;
            a.lin = 0;
        }
    }
    return a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int floor_div (int a, int b) noexcept
{
    return (a >= 0) ? a / b : -((-a + b - 1) / b);
}

template <bool Nodal>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void locate (const Axis& a, int i, int& ic, amrex::Real& t) noexcept
{
    const int den = Nodal ? a.r : 2 * a.r;
    const int m = Nodal ? (i - a.flo) : (2 * (i - a.flo) + 1 - a.r);
    int c = floor_div(Nodal ? m - 1 : m, den);
    c = amrex::min(amrex::max(c, a.imin), a.imax);
    t = a.lin * amrex::Real(m - den * c) * a.inv;
    ic = c + a.clo;
}

//! fine = the bilinear interpolation of the coarse image c (on fine's
//! BoxArray coarsened, fine's DistributionMapping, one ghost filled) at
//! fine's valid points of type T (stencil bits: nodal in x, y).
template <int T>
void interpolate (const amrex::MultiFab& c, amrex::MultiFab& fine, const Axis& ax, const Axis& ay)
{
    constexpr bool nx = stencil::nodal(T, 0);
    constexpr bool ny = stencil::nodal(T, 1);
    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto f = fine.array(mfi);
        const auto q = c.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int I = 0, J = 0;
            amrex::Real tx = 0, ty = 0;
            locate<nx>(ax, i, I, tx);
            locate<ny>(ay, j, J, ty);
            const int I1 = I + ax.step, J1 = J + ay.step;
            const amrex::Real one(1);
            f(i, j, k) = (one - ty) * ((one - tx) * q(I, J, k) + tx * q(I1, J, k))
                       + ty * ((one - tx) * q(I, J1, k) + tx * q(I1, J1, k));
        });
    }
}

void check_levels (const char* op, const MeshHierarchy& h, int lf, int lc)
{
    if (!(lf >= 0 && lf < lc && lc < h.num_levels())) {
        amrex::Abort(std::string("ablsteady: ") + op + ": the levels must satisfy 0 <= fine (" + std::to_string(lf)
                     + ") < coarse (" + std::to_string(lc) + ") < " + std::to_string(h.num_levels()));
    }
    const amrex::IntVect ratio(hierarchy_ratio(lc - lf), hierarchy_ratio(lc - lf), 1);
    if (!h.boxArray(lf).coarsenable(ratio)) {
        amrex::Abort(std::string("ablsteady: ") + op + ": level " + std::to_string(lf) + "'s boxes do not coarsen by "
                     + std::to_string(ratio[0]));
    }
}

void check_state (const char* op, const char* what, const State& q, const MeshHierarchy& h, int l)
{
    const Mesh& m = h.mesh(l);
    const std::string w(what);
    check_operator_field(op, (w + ".u").c_str(), q.u, Loc::xface, m, 0);
    check_operator_field(op, (w + ".v").c_str(), q.v, Loc::yface, m, 0);
    check_operator_field(op, (w + ".w").c_str(), q.w, Loc::zface, m, 0);
    check_operator_field(op, (w + ".p").c_str(), q.p, Loc::cell, m, 0);
    check_operator_field(op, (w + ".theta").c_str(), q.theta, Loc::cell, m, 0);
    check_operator_field(op, (w + ".tke").c_str(), q.tke, Loc::cell, m, 0);
}

//! The coarse image of fine's layout: fine's BoxArray coarsened by
//! (r, r, 1), at type `like`'s, on fine's DistributionMapping.
amrex::MultiFab image (const MeshHierarchy& h, int lf, int r, const amrex::MultiFab& like, int ngrow)
{
    const amrex::BoxArray ba = amrex::convert(amrex::coarsen(h.boxArray(lf), amrex::IntVect(r, r, 1)), like.ixType());
    return amrex::MultiFab(ba, h.DistributionMap(lf), 1, ngrow);
}

//! c = sum J q / sum J over the r x r fine cells of each coarse image cell.
void restrict_cells (const amrex::MultiFab& J, const amrex::MultiFab& q, int r, amrex::MultiFab& c)
{
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto jj = J.const_array(mfi);
        const auto qq = q.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int Jc, int k) noexcept
        {
            amrex::Real sjq = 0, sj = 0;
            for (int b = 0; b < r; ++b) {
                for (int a = 0; a < r; ++a) {
                    const amrex::Real w = jj(r * I + a, r * Jc + b, k);
                    sjq += w * qq(r * I + a, r * Jc + b, k);
                    sj += w;
                }
            }
            o(I, Jc, k) = sjq / sj;
        });
    }
}

//! c = (1/r) sum a q over the r fine x-faces (D = 0) or y-faces (D = 1) of
//! each coarse image face: the coarse flux per coarse unit width.
template <int D>
void restrict_faces (const amrex::MultiFab& area, const amrex::MultiFab& q, int r, amrex::MultiFab& c)
{
    const amrex::Real rinv = amrex::Real(1) / amrex::Real(r);
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto aa = area.const_array(mfi);
        const auto qq = q.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int J, int k) noexcept
        {
            amrex::Real s = 0;
            for (int b = 0; b < r; ++b) {
                const int i = (D == 0) ? r * I : r * I + b;
                const int j = (D == 0) ? r * J + b : r * J;
                s += aa(i, j, k) * qq(i, j, k);
            }
            o(I, J, k) = s * rinv;
        });
    }
}

//! c = (1/r^2) sum a_z Omega over the r x r fine z-faces of each coarse
//! image z-face.
void restrict_zfaces (const amrex::MultiFab& az, const amrex::MultiFab& omega, int r, amrex::MultiFab& c)
{
    const amrex::Real r2inv = amrex::Real(1) / amrex::Real(r * r);
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto aa = az.const_array(mfi);
        const auto qq = omega.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int J, int k) noexcept
        {
            amrex::Real s = 0;
            for (int b = 0; b < r; ++b) {
                for (int a = 0; a < r; ++a) {
                    s += aa(r * I + a, r * J + b, k) * qq(r * I + a, r * J + b, k);
                }
            }
            o(I, J, k) = s * r2inv;
        });
    }
}

//! The rows below the ground and above the top of a face field (one ghost
//! row each) take the nearest valid row.
void zero_gradient_rows (amrex::MultiFab& mf, const amrex::Box& domain)
{
    const int klo = domain.smallEnd(2), khi = domain.bigEnd(2);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const auto a = mf.array(mfi);
        if (vb.smallEnd(2) == klo) {
            amrex::Box row = vb;
            row.setSmall(2, klo - 1);
            row.setBig(2, klo - 1);
            amrex::ParallelFor(row, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = a(i, j, k + 1); });
        }
        if (vb.bigEnd(2) == khi) {
            amrex::Box row = vb;
            row.setSmall(2, khi + 1);
            row.setBig(2, khi + 1);
            amrex::ParallelFor(row, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = a(i, j, k - 1); });
        }
    }
}

} // namespace

void restrict_state (const MeshHierarchy& h, int lf, const State& fine, int lc, State& coarse, amrex::MultiFab* omega)
{
    const char* op = "restrict_state";
    check_levels(op, h, lf, lc);
    check_state(op, "fine", fine, h, lf);
    check_state(op, "coarse", coarse, h, lc);
    check_operator_field(op, "fine.u", fine.u, Loc::xface, h.mesh(lf), 1);
    check_operator_field(op, "fine.v", fine.v, Loc::yface, h.mesh(lf), 1);
    const int r = hierarchy_ratio(lc - lf);
    amrex::MultiFab::Copy(coarse.z_phys_nd, h.z_phys_nd(lc), 0, 0, 1, ngrow_nodal);
    const Mesh& mf = h.mesh(lf);
    const Mesh& mc = h.mesh(lc);
    const amrex::Geometry& gc = h.geom(lc);

    // Cell centres: the detJ-weighted mean, computed on the coarse image
    // and copied to the coarse layout (valid cells; no periodic shift).
    {
        const std::pair<const amrex::MultiFab*, amrex::MultiFab*> cells[3] = {
            {&fine.theta, &coarse.theta}, {&fine.tke, &coarse.tke}, {&fine.p, &coarse.p}};
        amrex::MultiFab img = image(h, lf, r, fine.p, 0);
        for (const auto& [src, dst] : cells) {
            restrict_cells(mf.detJ, *src, r, img);
            dst->ParallelCopy(img, 0, 0, 1);
        }
    }
    // u, v: the coarse flux per unit width, then divided by the coarse area.
    {
        amrex::MultiFab img = image(h, lf, r, fine.u, 0);
        restrict_faces<0>(mf.ax, fine.u, r, img);
        coarse.u.ParallelCopy(img, 0, 0, 1);
        amrex::MultiFab::Divide(coarse.u, mc.ax, 0, 0, 1, 0);
    }
    {
        amrex::MultiFab img = image(h, lf, r, fine.v, 0);
        restrict_faces<1>(mf.ay, fine.v, r, img);
        coarse.v.ParallelCopy(img, 0, 0, 1);
        amrex::MultiFab::Divide(coarse.v, mc.ay, 0, 0, 1, 0);
    }
    // w: the vertical flux Omega, restricted as u and v; then w from it on
    // the coarse mesh, the slope term reading u and v with zero-gradient rows
    // below the ground and above the top (Transfer.H).
    {
        amrex::MultiFab om_f(fine.w.boxArray(), fine.w.DistributionMap(), 1, 0);
        omega_from_w(fine.u, fine.v, fine.w, h.z_phys_nd(lf), mf, h.geom(lf), om_f, 0);
        amrex::MultiFab img = image(h, lf, r, fine.w, 0);
        restrict_zfaces(mf.az, om_f, r, img);
        amrex::MultiFab local;
        amrex::MultiFab& om_c = (omega != nullptr) ? *omega : local;
        if (om_c.empty()) { om_c.define(coarse.w.boxArray(), coarse.w.DistributionMap(), 1, 0); }
        check_operator_field(op, "omega", om_c, Loc::zface, mc, 0);
        om_c.ParallelCopy(img, 0, 0, 1);
        amrex::MultiFab::Divide(om_c, mc.az, 0, 0, 1, 0);

        amrex::MultiFab u1(coarse.u.boxArray(), coarse.u.DistributionMap(), 1, 1);
        amrex::MultiFab v1(coarse.v.boxArray(), coarse.v.DistributionMap(), 1, 1);
        amrex::MultiFab::Copy(u1, coarse.u, 0, 0, 1, 0);
        amrex::MultiFab::Copy(v1, coarse.v, 0, 0, 1, 0);
        u1.FillBoundary(gc.periodicity());
        v1.FillBoundary(gc.periodicity());
        zero_gradient_rows(u1, gc.Domain());
        zero_gradient_rows(v1, gc.Domain());
        w_from_omega(u1, v1, om_c, h.z_phys_nd(lc), mc, gc, coarse.w);
    }
}

void prolong_fields (const MeshHierarchy& h, int lc, const State& coarse, int lf, State& fine)
{
    const char* op = "prolong_fields";
    check_levels(op, h, lf, lc);
    check_state(op, "coarse", coarse, h, lc);
    check_state(op, "fine", fine, h, lf);
    const int r = hierarchy_ratio(lc - lf);
    const amrex::Geometry& gf = h.geom(lf);
    const amrex::Geometry& gc = h.geom(lc);
    const Axis cx = make_axis(gf, gc, 0, r, false), cy = make_axis(gf, gc, 1, r, false);
    const Axis nx = make_axis(gf, gc, 0, r, true), ny = make_axis(gf, gc, 1, r, true);
    amrex::MultiFab::Copy(fine.z_phys_nd, h.z_phys_nd(lf), 0, 0, 1, ngrow_nodal);

    // The coarse values on the coarse image of the fine layout, with one
    // ghost cell: the neighbouring boxes and the periodic images.
    auto img = [&] (const amrex::MultiFab& src) {
        amrex::MultiFab c = image(h, lf, r, src, 1);
        c.ParallelCopy(src, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), gc.periodicity());
        return c;
    };
    interpolate<stencil::CC>(img(coarse.theta), fine.theta, cx, cy);
    interpolate<stencil::CC>(img(coarse.tke), fine.tke, cx, cy);
    interpolate<stencil::CC>(img(coarse.p), fine.p, cx, cy);
    interpolate<stencil::XF>(img(coarse.u), fine.u, nx, cy);
    interpolate<stencil::YF>(img(coarse.v), fine.v, cx, ny);
    interpolate<stencil::ZF>(img(coarse.w), fine.w, cx, cy);
}

ProlongResult prolong_state (const MeshHierarchy& h, int lc, const State& coarse, int lf, State& fine,
                             const BoundaryConditions& bc, const SurfaceLayerOptions* surface,
                             const BuoyancyOptions& buoyancy, PoissonSolver& solver, amrex::Real tke_floor)
{
    if (!std::isfinite(tke_floor) || !(tke_floor >= 0)) {
        amrex::Abort("ablsteady: prolong_state: tke_floor must be finite and >= 0 (turb.tke_floor)");
    }
    ProlongResult res;
    prolong_fields(h, lc, coarse, lf, fine);
    const Mesh& mesh = h.mesh(lf);
    const amrex::Geometry& geom = h.geom(lf);

    // k >= the floor (before the boundary fill, so the ghosts see it).
    {
        amrex::ReduceOps<amrex::ReduceOpSum> rop;
        amrex::ReduceData<amrex::Long> data(rop);
        for (amrex::MFIter mfi(fine.tke, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto f = fine.tke.array(mfi);
            rop.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
            {
                const bool low = f(i, j, k) < tke_floor;
                if (low) { f(i, j, k) = tke_floor; }
                return {low ? amrex::Long(1) : amrex::Long(0)};
            });
        }
        res.tke_clipped = amrex::get<0>(data.value(rop));
        amrex::ParallelDescriptor::ReduceLongSum(res.tke_clipped);
    }

    bc.fill_ghosts(fine, surface, buoyancy);

    // The projection: SIMPLEC's correction with d = 1 off the Dirichlet rows.
    PredictorFields d;
    {
        const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
        amrex::MultiFab* dd[3] = {&d.d_u, &d.d_v, &d.d_w};
        for (int c = 0; c < 3; ++c) {
            bc.dirichlet_rows(eqs[c], velocity_loc(c), mesh, *dd[c]);
            // d = 1 - mask
            dd[c]->mult(amrex::Real(-1), 0, 1, 0);
            dd[c]->plus(amrex::Real(1), 0, 1, 0);
        }
    }
    amrex::MultiFab p_keep(fine.p.boxArray(), fine.p.DistributionMap(), 1, fine.p.nGrow());
    amrex::MultiFab::Copy(p_keep, fine.p, 0, 0, 1, fine.p.nGrow());
    SimplecWorkspace w;
    res.projection = pressure_correction(fine, d, solver, poisson_bc_from_faces(bc.inputs(), geom), amrex::Real(1), mesh,
                                         geom, w);
    amrex::MultiFab::Copy(fine.p, p_keep, 0, 0, 1, fine.p.nGrow());
    return res;
}

} // namespace ablsteady
