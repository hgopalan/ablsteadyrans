#include "WallDistance.H"

#include <ERF_TerrainMetrics.H>

#include <AMReX_Print.H>

#include <cmath>

namespace gust {

namespace {

// d from phi (ghosts filled for the Tucker BCs: odd below the ground). A
// free function: CUDA does not allow device lambdas in private members.
void distance_from_phi (const amrex::MultiFab& phi, const amrex::MultiFab& z, const amrex::Geometry& geom,
                        amrex::MultiFab& dist)
{
    const auto dxinv = geom.InvCellSizeArray();
    const int k0 = geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(dist, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto p = phi.const_array(mfi);
        const auto zp = z.const_array(mfi);
        const auto d = dist.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // ERF's cell-centred gradient in physical space (chain rule for a
            // mesh deformed in z only).
            const amrex::Real half = amrex::Real(0.5);
            const amrex::Real dpdxi = half * (p(i + 1, j, k) - p(i - 1, j, k)) * dxinv[0];
            const amrex::Real dpdeta = half * (p(i, j + 1, k) - p(i, j - 1, k)) * dxinv[1];
            const amrex::Real dpdzeta = half * (p(i, j, k + 1) - p(i, j, k - 1)) * dxinv[2];
            const amrex::Real h_zeta = Compute_h_zeta_AtCellCenter(i, j, k, dxinv, zp);
            const amrex::Real h_xi = Compute_h_xi_AtCellCenter(i, j, k, dxinv, zp);
            const amrex::Real h_eta = Compute_h_eta_AtCellCenter(i, j, k, dxinv, zp);
            amrex::Real dpdz = dpdzeta / h_zeta;
            if (k == k0) {
                // In the first cell, the quadratic along the column through
                // phi = 0 on the ground and the two lowest cell centres, at
                // heights a and b above the ground:
                // dphi/dz = (phi_0 b (b - 2a) + phi_1 a^2) / (a b (b - a)).
                // On uniform levels this is (phi_0 + phi_1/3)/dz; it cancels
                // the dz^2/8 offset of the discrete flat-ground solution better
                // than the centred difference with the odd ghost, and keeps that
                // gain on stretched levels.
                auto column = [&] (int kk) {
                    return amrex::Real(0.25) * (zp(i, j, kk) + zp(i + 1, j, kk) + zp(i, j + 1, kk) + zp(i + 1, j + 1, kk));
                };
                const amrex::Real zw = column(k);
                const amrex::Real a = half * (column(k) + column(k + 1)) - zw;
                const amrex::Real b = half * (column(k + 1) + column(k + 2)) - zw;
                dpdz = (p(i, j, k) * b * (b - amrex::Real(2) * a) + p(i, j, k + 1) * a * a) / (a * b * (b - a));
            }
            const amrex::Real dpdx = dpdxi - h_xi * dpdz;
            const amrex::Real dpdy = dpdeta - h_eta * dpdz;
            const amrex::Real g2 = dpdx * dpdx + dpdy * dpdy + dpdz * dpdz;
            // Tucker 2003, Eq. 2
            d(i, j, k) = -std::sqrt(g2) + std::sqrt(g2 + amrex::Real(2) * p(i, j, k));
        });
    }
}

} // namespace

PoissonBC tucker_bc (const amrex::Geometry& geom)
{
    PoissonBC bc;
    for (int d = 0; d < 2; ++d) {
        const PoissonFace f = geom.isPeriodic(d) ? PoissonFace::periodic : PoissonFace::neumann;
        bc.face[amrex::Orientation(d, amrex::Orientation::low)] = f;
        bc.face[amrex::Orientation(d, amrex::Orientation::high)] = f;
    }
    bc.face[amrex::Orientation(2, amrex::Orientation::low)] = PoissonFace::dirichlet;
    bc.face[amrex::Orientation(2, amrex::Orientation::high)] = PoissonFace::neumann;
    return bc;
}

PoissonResult tucker_wall_distance (const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd,
                                    amrex::MultiFab& dist, const PoissonOptions& opt)
{
    if (dist.boxArray() != mesh.detJ.boxArray() || dist.DistributionMap() != mesh.detJ.DistributionMap()) {
        amrex::Abort("gust: tucker_wall_distance: dist must be on the mesh's BoxArray and DistributionMapping");
    }
    const PoissonBC bc = tucker_bc(geom);
    amrex::MultiFab rhs(mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    amrex::MultiFab phi(mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), 1, 1);
    rhs.setVal(-1.0);
    phi.setVal(0.0);
    PoissonSolver solver(mesh, geom, z_phys_nd, opt);
    const PoissonResult res = solver.solve(rhs, phi, {}, bc);
    solver.fill_ghosts(phi, bc);
    distance_from_phi(phi, z_phys_nd, geom, dist);
    return res;
}

} // namespace gust
