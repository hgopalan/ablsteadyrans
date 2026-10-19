#include "Mesh.H"

#include "TerrainSurface.H"
#include "ZLevels.H"

#include <ERF_IndexDefines.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_Utils.H>

#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>
#include <cstdio>
#include <sstream>

namespace ablsteady {

namespace {

constexpr amrex::Real kRadToDeg = amrex::Real(57.29577951308232);

// ERF's boundary types for the terrain construction: only symmetry would
// change it (mirrored ghost nodes), and ablsteady has none yet.
amrex::GpuArray<ERF_BC, AMREX_SPACEDIM * 2> terrain_bc (const amrex::Geometry& geom)
{
    amrex::GpuArray<ERF_BC, AMREX_SPACEDIM * 2> bc{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const ERF_BC side = geom.isPeriodic(d) ? ERF_BC::periodic : ERF_BC::outflow;
        bc[amrex::Orientation(d, amrex::Orientation::low)] = side;
        bc[amrex::Orientation(d, amrex::Orientation::high)] = side;
    }
    return bc;
}

//! Make every ghost node of z_phys_nd stand for the node it represents, at
//! every k. After ERF's construction and FillBoundary two kinds were left:
//!  - ghost columns outside a non-periodic face whose periodic image in the
//!    other direction is a ghost too (corners such as (-1, -1) with x walls
//!    and periodic y): ERF clamps them in both directions, to the column of
//!    (0, 0) instead of (0, ny - 1);
//!  - the rows below the ground and above the top of ghost columns: ERF
//!    extrapolates them before FillBoundary, from the clamped columns, and
//!    FillBoundary does not reach them (they are outside every valid box).
//! Here, after FillBoundary: (1) a node outside a non-periodic face takes
//! the height of the edge node of its own row (clamp only in the
//! non-periodic direction; the edge column's periodic ghosts are already
//! filled); (2) every column's rows below the ground and above the top are
//! recomputed from that column with ERF's own expressions, 2 z_0 - z_1 and
//! z_n + (k - n)(z_n - z_{n-1}), so a valid column keeps its values when
//! nothing above changed it, and a ghost column gets exactly its image's.
void fill_ghost_columns (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom)
{
    const amrex::Box& domain = geom.Domain();
    const int domlo_x = domain.smallEnd(0);
    const int domhi_x = domain.bigEnd(0) + 1;   // nodal
    const int domlo_y = domain.smallEnd(1);
    const int domhi_y = domain.bigEnd(1) + 1;
    const int k0 = domain.smallEnd(2);
    const int kn = domain.bigEnd(2) + 1;
    const bool per_x = geom.isPeriodic(0);
    const bool per_y = geom.isPeriodic(1);
    for (amrex::MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const amrex::Box gbx = mfi.fabbox();
        const auto z = z_phys_nd.array(mfi);
        amrex::Box lat = gbx;
        lat.setSmall(2, k0);
        lat.setBig(2, kn);
        amrex::ParallelFor(lat, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ic = per_x ? i : amrex::max(domlo_x, amrex::min(domhi_x, i));
            const int jc = per_y ? j : amrex::max(domlo_y, amrex::min(domhi_y, j));
            if (ic != i || jc != j) { z(i, j, k) = z(ic, jc, k); }
        });
        amrex::Box cols = gbx;
        cols.setSmall(2, k0);
        cols.setBig(2, k0);
        const int klo = gbx.smallEnd(2);
        const int khi = gbx.bigEnd(2);
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            constexpr amrex::Real two = 2;
            if (klo <= k0 - 1) { z(i, j, k0 - 1) = two * z(i, j, k0) - z(i, j, k0 + 1); }
            for (int k = kn + 1; k <= khi; ++k) {
                z(i, j, k) = z(i, j, kn) + (k - kn) * (z(i, j, kn) - z(i, j, kn - 1));
            }
        });
    }
}

} // namespace

void Mesh::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    using amrex::IntVect;
    detJ.define(ba, dm, 1, 1);
    ax.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 1);
    ay.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 1);
    az.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 1);
    z_cc.define(ba, dm, 1, 1);
    for (auto* mf : {&detJ, &ax, &ay, &az, &z_cc}) { mf->setVal(0.0); }
}

void Mesh::build (const Inputs& in, const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd)
{
    fill_terrain_surface(z_phys_nd, geom, in.terrain);
    build_from_ground(in, geom, z_phys_nd);
}

void Mesh::build_from_ground (const Inputs& in, const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd)
{
    m_zlevels = make_zlevels(in.mesh, geom);
    m_smoothing = in.terrain.smoothing;
    auto bc = terrain_bc(geom);
    int n_capped = 0;
    make_terrain_fitted_coords(0, geom, z_phys_nd, m_zlevels, bc, in.terrain.smoothing,
                               in.terrain.gamma_m, &n_capped);
    m_quality.stf_capped_levels = n_capped;
    // A node shared by two boxes (or a periodic image) gets one value. Each
    // box computes its own copy, and where a compiler contracts multiply-adds
    // (gcc) the copies differ in the last bits (1e-13 m on a BTF hill); the
    // operators read each box's copy, so the boxes that share a face would
    // compute different values there. The ghost nodes are refilled to match.
    z_phys_nd.OverrideSync(geom.periodicity());
    z_phys_nd.FillBoundary(geom.periodicity());
    fill_ghost_columns(z_phys_nd, geom);

    make_metrics(geom, z_phys_nd);
    make_quality(geom, z_phys_nd);

    const auto& q = m_quality;
    if (!(q.max_ground_height < m_zlevels.back())) {
        amrex::Abort("ablsteady: the terrain (terrain.type = " + in.terrain.type + ") reaches "
                     + std::to_string(q.max_ground_height) + " m, at or above the top level "
                     + std::to_string(m_zlevels.back()) + " m; raise mesh.prob_hi z");
    }
    if (!(q.min_thickness_ratio > 0) || !(q.min_nodal_thickness_ratio > 0)) {
        amrex::Abort("ablsteady: the terrain-following mesh has inverted or empty cells (min layer thickness "
                     "ratio " + std::to_string(q.min_nodal_thickness_ratio) + "); raise mesh.prob_hi z or "
                     "use terrain.smoothing = 1 (STF)");
    }
}

void Mesh::check_restart (const Inputs& in, const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd)
{
    amrex::MultiFab rebuilt(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, z_phys_nd.nGrowVect());
    rebuilt.setVal(0.0);
    build(in, geom, rebuilt);

    // Valid nodes only; compared exactly, with no arithmetic on the values.
    long ndiff = 0;
    for (amrex::MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto a = z_phys_nd.const_array(mfi);
        const auto b = rebuilt.const_array(mfi);
        amrex::ReduceOps<amrex::ReduceOpSum> op;
        amrex::ReduceData<long> data(op);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<long> {
            return {(a(i, j, k) == b(i, j, k)) ? 0L : 1L};
        });
        ndiff += amrex::get<0>(data.value());
    }
    amrex::ParallelAllReduce::Sum(ndiff, amrex::ParallelDescriptor::Communicator());
    if (ndiff != 0) {
        amrex::Abort("ablsteady: io.restart: the checkpoint's node heights differ from the ones the "
                     "mesh and terrain inputs give at " + std::to_string(ndiff)
                     + " nodes; restart with the inputs the checkpoint was written with");
    }
    // The ghost nodes too: a checkpoint written before fill_ghost_columns
    // holds the old corner and below-ground ghost heights (valid nodes are
    // equal, checked above).
    amrex::MultiFab::Copy(z_phys_nd, rebuilt, 0, 0, 1, z_phys_nd.nGrowVect());
}

void Mesh::make_metrics (const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd)
{
    make_J(geom, z_phys_nd, detJ);
    make_areas(geom, z_phys_nd, ax, ay, az);
    make_zcc(geom, z_phys_nd, z_cc);
}

void Mesh::make_quality (const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd)
{
    const auto dx = geom.CellSizeArray();
    const amrex::Real dz = dx[2];
    const amrex::Real hmax_dx = std::max(dx[0], dx[1]);
    const int k0 = geom.Domain().smallEnd(2);
    const int nz = geom.Domain().length(2);

    amrex::Gpu::DeviceVector<amrex::Real> dzlev(nz);
    amrex::Vector<amrex::Real> dzlev_h(nz);
    for (int k = 0; k < nz; ++k) { dzlev_h[k] = m_zlevels[k + 1] - m_zlevels[k]; }
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dzlev_h.begin(), dzlev_h.end(), dzlev.begin());
    const amrex::Real* dzl = dzlev.data();

    // Per cell: thickness ratio and aspect ratio. Per node: slopes of the
    // coordinate surfaces (and of the ground at k = 0), the nodal layer
    // thickness ratio, the ground height and the top-node displacement. Per
    // first-layer cell: centre height above the mean of its four ground nodes.
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpMin> cell_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> cell_data(cell_ops);
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax,
                     amrex::ReduceOpMax> node_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> node_data(node_ops);
    const amrex::Real ztop_level = m_zlevels[nz];
    const amrex::Real lowest = std::numeric_limits<amrex::Real>::lowest();
    const amrex::Real big = std::numeric_limits<amrex::Real>::max();

    for (amrex::MFIter mfi(detJ, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto J = detJ.const_array(mfi);
        const auto zc = z_cc.const_array(mfi);
        const auto z = z_phys_nd.const_array(mfi);
        cell_ops.eval(bx, cell_data, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real>
        {
            const amrex::Real thick = J(i, j, k) * dz;
            const amrex::Real ratio = thick / dzl[k - k0];
            const amrex::Real aspect = hmax_dx / thick;
            amrex::Real first = big;
            if (k == k0) {
                const amrex::Real ground = amrex::Real(0.25) * (z(i, j, k) + z(i + 1, j, k)
                                                                + z(i, j + 1, k) + z(i + 1, j + 1, k));
                first = zc(i, j, k) - ground;
            }
            return {ratio, aspect, first};
        });

        // Every node of the box; slopes only where the x and y neighbours are
        // in the box too.
        const amrex::Box nbx = amrex::convert(bx, amrex::IntVect(1));
        const int ihi = nbx.bigEnd(0);
        const int jhi = nbx.bigEnd(1);
        node_ops.eval(nbx, node_data, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>
        {
            amrex::Real s = 0;
            if (i < ihi && j < jhi) {
                const amrex::Real sx = (z(i + 1, j, k) - z(i, j, k)) / dx[0];
                const amrex::Real sy = (z(i, j + 1, k) - z(i, j, k)) / dx[1];
                s = std::sqrt(sx * sx + sy * sy);
            }
            const bool top = (k == k0 + nz);
            const amrex::Real thick = top ? big : (z(i, j, k + 1) - z(i, j, k)) / dzl[k - k0];
            return {s, (k == k0) ? s : amrex::Real(0), thick, (k == k0) ? z(i, j, k) : lowest,
                    top ? std::abs(z(i, j, k) - ztop_level) : amrex::Real(0)};
        });
    }
    auto cv = cell_data.value();
    auto nv = node_data.value();
    amrex::Real min_ratio = amrex::get<0>(cv);
    amrex::Real max_aspect = amrex::get<1>(cv);
    amrex::Real min_first = amrex::get<2>(cv);
    amrex::Real max_slope = amrex::get<0>(nv);
    amrex::Real max_ground = amrex::get<1>(nv);
    amrex::Real min_nodal = amrex::get<2>(nv);
    amrex::Real max_height = amrex::get<3>(nv);
    amrex::Real max_top = amrex::get<4>(nv);
    const auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Min(min_ratio, comm);
    amrex::ParallelAllReduce::Max(max_aspect, comm);
    amrex::ParallelAllReduce::Min(min_first, comm);
    amrex::ParallelAllReduce::Max(max_slope, comm);
    amrex::ParallelAllReduce::Max(max_ground, comm);
    amrex::ParallelAllReduce::Min(min_nodal, comm);
    amrex::ParallelAllReduce::Max(max_height, comm);
    amrex::ParallelAllReduce::Max(max_top, comm);

    m_quality.min_thickness_ratio = min_ratio;
    m_quality.min_nodal_thickness_ratio = min_nodal;
    m_quality.max_ground_height = max_height;
    m_quality.max_top_displacement = max_top;
    m_quality.top_level = ztop_level;
    m_quality.max_aspect_ratio = max_aspect;
    m_quality.min_first_centre_height = min_first;
    m_quality.max_surface_slope_deg = std::atan(max_slope) * kRadToDeg;
    m_quality.max_ground_slope_deg = std::atan(max_ground) * kRadToDeg;
    m_quality.nz = nz;
}

std::string Mesh::quality_report () const
{
    const auto& q = m_quality;
    std::ostringstream os;
    char line[160];
    auto put = [&] (const char* name, double v) {
        std::snprintf(line, sizeof(line), "  %-42s = %.9e\n", name, v);
        os << line;
    };
    os << "ablsteady: mesh quality\n";
    put("min cell thickness / flat layer thickness", q.min_thickness_ratio);
    put("min nodal layer thickness / flat thickness", q.min_nodal_thickness_ratio);
    put("max ground height (m)", q.max_ground_height);
    put("max ground slope (deg)", q.max_ground_slope_deg);
    put("max coordinate-surface slope (deg)", q.max_surface_slope_deg);
    put("max aspect ratio max(dx,dy) / thickness", q.max_aspect_ratio);
    put("min first-cell centre above ground (m)", q.min_first_centre_height);
    put("max top-node displacement (m)", q.max_top_displacement);
    std::snprintf(line, sizeof(line), "  %-42s = %d of %d\n", "STF levels capped by terrain.gamma_m",
                  q.stf_capped_levels, q.nz);
    os << line;
    if (q.stf_capped_levels > 0) {
        os << "WARNING: the STF smoothing stopped early at " << q.stf_capped_levels
           << " levels (terrain.gamma_m); the terrain is steep for this vertical grid\n";
    }
    if (q.max_top_displacement > amrex::Real(1.e-6) * q.top_level) {
        os << "WARNING: the terrain imprint reaches the top: top nodes are up to " << q.max_top_displacement
           << " m from the flat top level; raise mesh.prob_hi z\n";
    }
    return os.str();
}

amrex::Real free_stream_divergence (const Mesh& mesh, const amrex::Geometry& geom,
                                    const amrex::MultiFab& z_phys_nd, amrex::Real U, amrex::Real V)
{
    const auto dxinv = geom.InvCellSizeArray();
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real> data(op);
    for (amrex::MFIter mfi(mesh.detJ, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto J = mesh.detJ.const_array(mfi);
        const auto Ax = mesh.ax.const_array(mfi);
        const auto Ay = mesh.ay.const_array(mfi);
        const auto z = z_phys_nd.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
        {
            const amrex::Real zx_lo = Compute_h_xi_AtKface(i, j, k, dxinv, z);
            const amrex::Real zx_hi = Compute_h_xi_AtKface(i, j, k + 1, dxinv, z);
            const amrex::Real zy_lo = Compute_h_eta_AtKface(i, j, k, dxinv, z);
            const amrex::Real zy_hi = Compute_h_eta_AtKface(i, j, k + 1, dxinv, z);
            const amrex::Real div = U * (Ax(i + 1, j, k) - Ax(i, j, k)) * dxinv[0]
                                  + V * (Ay(i, j + 1, k) - Ay(i, j, k)) * dxinv[1]
                                  - (U * (zx_hi - zx_lo) + V * (zy_hi - zy_lo)) * dxinv[2];
            return {std::abs(div / J(i, j, k))};
        });
    }
    amrex::Real m = amrex::get<0>(data.value());
    amrex::ParallelAllReduce::Max(m, amrex::ParallelDescriptor::Communicator());
    return m;
}

} // namespace ablsteady
