#include <ERF_TerrainMetrics.H>
#include <ERF_Utils.H>
#include <AMReX_ParmParse.H>
#include <ERF_Constants.H>
#include <ERF_Interpolation_1D.H>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

using namespace amrex;

void
validate_flat_terrain (int lev,
                       const MultiFab& z_phys_nd,
                       const Vector<Real>& z_levels_h)
{
    // The comparison is against heights that have been round-tripped through z_phys_nd, so the
    // threshold has to track the working precision: a hard-wired 1e-10 is far below single-
    // precision round-off at domain-top magnitudes and would abort on a genuinely flat mesh.
    // Keep the original (double-precision) threshold, floored at a few hundred ulps.
    const Real scale = std::max(Real(1.0), std::abs(z_levels_h.back()));
    const Real rel_tol = std::max(Real(1.0e-10), Real(100.0) * std::numeric_limits<Real>::epsilon());
    const Real tolerance = rel_tol * scale;

    MultiFab expected(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, 0);
    Gpu::DeviceVector<Real> z_levels_d(z_levels_h.size());
    Gpu::copy(Gpu::hostToDevice, z_levels_h.begin(), z_levels_h.end(), z_levels_d.begin());
    const Real* z_levels = z_levels_d.data();

    for (MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& expected_arr = expected.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            expected_arr(i,j,k) = z_levels[k];
        });
    }

    MultiFab error(expected.boxArray(), expected.DistributionMap(), 1, 0);
    MultiFab::Copy(error, expected, 0, 0, 1, 0);
    MultiFab::Subtract(error, z_phys_nd, 0, 0, 1, 0);
    const Real max_error = error.norminf(0, 0, false);

    if (max_error > tolerance) {
        Abort("erf.flat_terrain = true, but z_phys_nd is not horizontally "
              "flat or does not match the configured vertical spacing at "
              "level " + std::to_string(lev) + ".");
    }
}

/**
 * Define a default z_phys so we have it even if a completely regular grid
 * This will be over-written if we use z_levels, or grid stretching, or terrain-fitted grids
 */
void
init_default_zphys (int /*lev*/,
                    const Geometry& geom,
                    MultiFab& z_phys_nd,
                    MultiFab& z_phys_cc,
                    Real z_offset)
{
    const auto& dx = geom.CellSize();
    Real dz = dx[2];

    for (MFIter mfi(z_phys_nd,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        const Array4< Real> z_nd_arr = z_phys_nd.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            z_nd_arr(i,j,k) = k * dz - z_offset;
        });
    }

    for (MFIter mfi(z_phys_cc,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        const Array4< Real> z_cc_arr = z_phys_cc.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            z_cc_arr(i,j,k) = (k + myhalf) * dz - z_offset;
        });
    }
}

/**
 * Computation of the terrain grid from BTF, STF, or Sullivan TF model
 */

void
make_terrain_fitted_coords (int lev,
                            const Geometry& geom,
                            MultiFab& z_phys_nd,
                            Vector<Real> const& z_levels_h,
                            GpuArray<ERF_BC, AMREX_SPACEDIM*2>& phys_bc_type,
                            int terrain_smoothing,
                            Real terrain_gamma_m,
                            int* n_capped)
{
    const Box& domain = geom.Domain();

    int domlo_z = domain.smallEnd(2);
    int domhi_z = domain.bigEnd(2) + 1;

    // Just in case ...
    z_phys_nd.setDomainBndry(bogus_large_value,0,1,geom);

    // ****************************************************************************

    if (lev == 0) {
        BoxArray ba(z_phys_nd.boxArray());
        bool all_boxes_touch_bottom = true;
        for (int i = 0; i < ba.size(); i++) {
            if (ba[i].smallEnd(2) != domlo_z) {
                all_boxes_touch_bottom = false;
            }
        }

        if (all_boxes_touch_bottom) {
            init_which_terrain_grid(lev, geom, z_phys_nd, z_levels_h, terrain_smoothing, terrain_gamma_m, n_capped);
        } else {

            BoxArray ba_new(domain);
            ChopGrids2D(ba_new, domain, ParallelDescriptor::NProcs());

            DistributionMapping dm_new(ba_new);
            ba_new.surroundingNodes();

            MultiFab z_phys_nd_new(ba_new, dm_new, 1, z_phys_nd.nGrowVect());

            z_phys_nd_new.ParallelCopy(z_phys_nd,0,0,1,z_phys_nd.nGrowVect(),z_phys_nd.nGrowVect());

            init_which_terrain_grid(lev, geom, z_phys_nd_new, z_levels_h, terrain_smoothing, terrain_gamma_m, n_capped);

            z_phys_nd.ParallelCopy(z_phys_nd_new,0,0,1,z_phys_nd.nGrowVect(),z_phys_nd.nGrowVect());
        }
    } else {
        // ablsteady builds every grid level as its own full-domain mesh (lev = 0);
        // ERF's refined-level path is not copied.
        Abort("make_terrain_fitted_coords: only lev = 0 is supported");
    }

    //
    // Fill ghost layers and corners (including periodic) -- no matter what level
    //
    z_phys_nd.FillBoundary(geom.periodicity());

    if (phys_bc_type[Orientation(0,Orientation::low )] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(0,Orientation::high)] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(1,Orientation::low )] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(1,Orientation::high)] == ERF_BC::symmetry) {

        const auto& dom_lo = lbound(convert(domain,IntVect(1,1,1)));
        const auto& dom_hi = ubound(convert(domain,IntVect(1,1,1)));

        //
        // Mirror the whole ghost region, not just its first layer: the base state carries
        // ComputeGhostCells+1 layers and the state stencils reach several cells past the
        // boundary, so mirroring only the first would leave the rest of the ghost region
        // holding whatever the terrain construction extrapolated into it.
        //
        const IntVect& ngv = z_phys_nd.nGrowVect();

        for (MFIter mfi(z_phys_nd,true); mfi.isValid(); ++mfi) {
            //
            // Whether this tile reaches a domain boundary is a question about its valid
            // nodes, so ask it of the tilebox.  The grown tilebox begins a ghost width short
            // of the boundary and so never equals it -- asking it of the grown box is what
            // kept this mirroring from ever running.  The grown box is still what we sweep,
            // so that the ghost rows in the other two directions are covered as well.
            //
            const Box& vbx = mfi.tilebox();
            const Box& bx  = mfi.growntilebox();
            const Array4< Real> z_nd_arr = z_phys_nd.array(mfi);
            if (phys_bc_type[Orientation(0,Orientation::low)] == ERF_BC::symmetry && vbx.smallEnd(0) == dom_lo.x) {
                const int ng = ngv[0];
                ParallelFor(makeSlab(bx,0,1), [=] AMREX_GPU_DEVICE (int , int j, int k)
                {
                    for (int n = 1; n <= ng; ++n) {
                        z_nd_arr(dom_lo.x-n,j,k) = z_nd_arr(dom_lo.x+n,j,k);
                    }
                });
            }
            if (phys_bc_type[Orientation(0,Orientation::high)] == ERF_BC::symmetry && vbx.bigEnd(0) == dom_hi.x) {
                const int ng = ngv[0];
                ParallelFor(makeSlab(bx,0,1), [=] AMREX_GPU_DEVICE (int , int j, int k)
                {
                    for (int n = 1; n <= ng; ++n) {
                        z_nd_arr(dom_hi.x+n,j,k) = z_nd_arr(dom_hi.x-n,j,k);
                    }
                });
            }
            if (phys_bc_type[Orientation(1,Orientation::low)] == ERF_BC::symmetry && vbx.smallEnd(1) == dom_lo.y) {
                const int ng = ngv[1];
                ParallelFor(makeSlab(bx,1,1), [=] AMREX_GPU_DEVICE (int i, int  , int k)
                {
                    for (int n = 1; n <= ng; ++n) {
                        z_nd_arr(i,dom_lo.y-n,k) = z_nd_arr(i,dom_lo.y+n,k);
                    }
                });
            }
            if (phys_bc_type[Orientation(1,Orientation::high)] == ERF_BC::symmetry && vbx.bigEnd(1) == dom_hi.y) {
                const int ng = ngv[1];
                ParallelFor(makeSlab(bx,1,1), [=] AMREX_GPU_DEVICE (int i, int  , int k)
                {
                    for (int n = 1; n <= ng; ++n) {
                        z_nd_arr(i,dom_hi.y+n,k) = z_nd_arr(i,dom_hi.y-n,k);
                    }
                });
            }
        }
    }

    //********************************************************************************
    // Populate domain boundary cells in z-direction
    //********************************************************************************
    int ngrow = z_phys_nd.nGrow();

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Only set values above top of domain if this box reaches that far
        Box nd_bx = mfi.tilebox();

        // Note that domhi_z is already nodal in the z-direction
        if (nd_bx.bigEnd(2) >= domhi_z) {
            // Grown box with no z range
            Box bx_zhi = mfi.growntilebox(ngrow);
            bx_zhi.setSmall(2,domhi_z+1);
            Array4<Real> const& z_arr = z_phys_nd.array(mfi);

            // Extrapolate top layer
            ParallelFor(bx_zhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                z_arr(i,j,k) = z_arr(i,j,domhi_z)
                             + (k-domhi_z) * (z_arr(i,j,domhi_z) - z_arr(i,j,domhi_z-1));
            });
        }
    }
} // make_terrain_fitted_coords

/**
 * Initialize the terrain grid using the selected model (BTF, STF, or Sullivan TF).
 *
 * @param[in] lev Level in the AMR hierarchy.
 * @param[in] geom Geometry defining the domain.
 * @param[out] z_phys_nd MultiFab to be filled with nodal physical heights.
 * @param[in] z_levels_h Reference z-levels for grid stretching.
 */
void
init_which_terrain_grid (int lev,
                         Geometry const& geom,
                         MultiFab& z_phys_nd,
                         Vector<Real> const& z_levels_h,
                         int terrain_smoothing,
                         Real terrain_gamma_m,
                         int* n_capped_out)
{
    // ablsteady: the method and gamma_m come from the caller (terrain.smoothing,
    // terrain.gamma_m) instead of ParmParse, and only lev = 0 is built.
    if (n_capped_out != nullptr) { *n_capped_out = 0; }

    if (lev > 0) {
        Abort("init_which_terrain_grid: only lev = 0 is supported");
    }

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();
    IntVect ngrowVect = z_phys_nd.nGrowVect();

    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2); int domhi_z = domain.bigEnd(2) + 1;

    int imin = domlo_x; // if (geom.isPeriodic(0)) imin -= z_phys_nd.nGrowVect()[0];
    int jmin = domlo_y; // if (geom.isPeriodic(1)) jmin -= z_phys_nd.nGrowVect()[1];

    int imax = domhi_x; // if (geom.isPeriodic(0)) imax += z_phys_nd.nGrowVect()[0];
    int jmax = domhi_y; // if (geom.isPeriodic(1)) jmax += z_phys_nd.nGrowVect()[1];

    int nz = static_cast<int>(z_levels_h.size());
    Real z_top = z_levels_h[nz-1];

    Gpu::DeviceVector<Real> z_levels_d;
    z_levels_d.resize(nz);
    Gpu::copy(Gpu::hostToDevice, z_levels_h.begin(), z_levels_h.end(), z_levels_d.begin());

    switch(terrain_smoothing) {
    case 0: // BTF Method
    {
        for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Note that this box is nodal because it is based on z_phys_nd
            const Box&  bx = mfi.validbox();

            int k0 = bx.smallEnd()[2];

            // Grown box with corrected ghost cells at top
            Box gbx = mfi.growntilebox(ngrowVect);

            if (bx.smallEnd(2) == domlo_z) {
                gbx.setSmall(2,domlo_z);
            } else {
                gbx.growLo(2,-1);
            }

            // Note that we don't overwrite the values at the high end of the box
            // regardless of whether the box reaches the top of the domain or not.
            // In the case of lev > 0, this ensures that the fine nodes at the top
            // of a fine box are those that are interpolated from the coarse grid
            if (bx.bigEnd(2) == domhi_z) {
                gbx.setBig(2,domhi_z);
            } else {
                gbx.growHi(2,-1);
                if (gbx.bigEnd(2) > domhi_z) {
                    gbx.setBig(2,domhi_z);
                }
            }

            Array4<Real> const& z_arr = z_phys_nd.array(mfi);
            auto const&         z_lev = z_levels_d.data();

            //
            // Vertical grid stretching using BTF model from p2163 of Klemp2011
            // z_levels are only defined from k = dom_lo to dom_hi (nodal)
            //
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // ablsteady: the ground layer k0 is not rewritten. ERF writes
                // z_sfc z_top / z_top there while other threads read z_arr(ii,jj,k0),
                // a race on GPU; the value differs from z_sfc by at most 1 ulp.
                if (k == k0) { return; }
                int ii = amrex::max(amrex::min(i,imax),imin);
                int jj = amrex::max(amrex::min(j,jmax),jmin);

                //
                // Start with flat z_lev set either with uniform cell size or specified z_levels
                // If k0 = 0 then z_arr at k0 has already been filled from the terrain data
                // If k0 > 0 then z_arr at k0 has already been filled from interpolation
                //
                Real z         = z_lev[k];
                Real z_sfc     = z_arr(ii,jj,k0);
                Real z_lev_sfc = z_lev[k0];

                z_arr(i,j,k) = ( (z_sfc - z_lev_sfc) * z_top +
                                 (z_top - z_sfc    ) * z     ) / (z_top - z_lev_sfc);
            });
        } // mfi

        z_phys_nd.FillBoundary(geom.periodicity());

        for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Note that this box is nodal because it is based on z_phys_nd
            const Box&  bx = mfi.validbox();
            Box gbx = mfi.growntilebox(ngrowVect);

            int k0 = bx.smallEnd()[2];

            if (k0 == 0)
            {
                Array4<Real> const& z_arr = z_phys_nd.array(mfi);

                // Fill lateral boundaries below the bottom surface
                ParallelFor(makeSlab(gbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    z_arr(i,j,-1) = two*z_arr(i,j,0) - z_arr(i,j,1);
                });
            }
        } // mfi
        break;
    } // case 0

    case 1: // STF Method
    {
        // Get MultiFab spanning domain with 1 level of ghost cells
        MultiFab h_mf(    z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, ngrow+1);
        MultiFab h_mf_old(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, ngrow+1);

        // Save max height for smoothing
        Real h_m;

        // Create 2D MF without allocation
        MultiFab mf2d;
        {
            BoxList bl2d = h_mf.boxArray().boxList();
            for (auto& b : bl2d) { b.setRange(2,b.smallEnd(2)); }
            BoxArray ba2d(std::move(bl2d));
            mf2d = MultiFab(ba2d, h_mf.DistributionMap(), 1, ngrow, MFInfo().SetAlloc(false));
        }

        // Get MultiArray4s from the multifabs
        MultiArray4<Real> const& ma_h_s     = h_mf.arrays();
        MultiArray4<Real> const& ma_h_s_old = h_mf_old.arrays();
        MultiArray4<Real> const& ma_z_phys  = z_phys_nd.arrays();

        // Bottom boundary
        int k0 = domlo_z;

        // Get max value
        h_m = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, mf2d, IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> GpuTuple<Real>
                {
                  // Get Array4s
                  const auto & h     = ma_h_s[box_no];
                  const auto & z_arr = ma_z_phys[box_no];

                  int ii = amrex::max(amrex::min(i,imax),imin);
                  int jj = amrex::max(amrex::min(j,jmax),jmin);

                  // Fill the lateral boundaries
                  z_arr(i,j,k0) = z_arr(ii,jj,k0);

                  // Populate h with terrain
                  h(i,j,k0) = z_arr(i,j,k0);

                  // Return height for max
                  return { z_arr(i,j,k0) };
                });
        amrex::ParallelDescriptor::ReduceRealMax(h_m);

        if (h_m < std::numeric_limits<Real>::epsilon()) h_m = Real(1e-16);

        // Fill ghost cells (neglects domain boundary if not periodic)
        h_mf.FillBoundary(geom.periodicity());

        // Make h_mf copy for old values
        MultiFab::Copy(h_mf_old, h_mf,0,0,1,h_mf_old.nGrow());

        // Minimum allowed fractional grid spacing
        Real gamma_m = terrain_gamma_m;
        Real z_H     = Real(2.44)*h_m/(1-gamma_m); // Klemp2011 Eqn. 11

        // Levels at which the minimum-thickness criterion forced us to stop smoothing
        // early, reported once below: a terrain that needs this at many levels is at the
        // edge of what STF can represent on this vertical grid.
        int n_capped = 0;

        // Populate h_mf at k>0 with h_s, solving in ordered 2D slices
        for (int k = domlo_z+1; k <= domhi_z; k++) // skip terrain level
        {
            auto const& z_lev_h = z_levels_h.data();

            Real zz       = z_lev_h[k];
            Real zz_minus = z_lev_h[k-1];

            // Hybrid attenuation profile, Klemp2011 Eqn. 9
            Real A;
            Real foo = std::cos((PI/2)*(zz/z_H));
            if(zz < z_H) { A = foo*foo*foo*foo*foo*foo; } // A controls rate of return to atm
            else         { A = 0; }
            Real foo_minus = std::cos((PI/2)*(zz_minus/z_H));
            Real A_minus;
            if(zz_minus < z_H) { A_minus = foo_minus*foo_minus*foo_minus*foo_minus*foo_minus*foo_minus; } // A controls rate of return to atm
            else               { A_minus = 0; }

            unsigned maxIter = 50; // M_k in paper
            unsigned iter    = 0;
            Real threshold   = gamma_m;
            Real diff        = bogus_large_value;
            while (iter < maxIter)
            {

                diff = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, mf2d, IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> GpuTuple<Real>
                {
                    const auto & h_s     = ma_h_s[box_no];
                    const auto & h_s_old = ma_h_s_old[box_no];

                    Real beta_k = Real(0.2)*std::min(zz/(2*h_m),Real(1.0)); //smoothing coefficient (Eqn. 8)

                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                    int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                    if (iter == 0) {
                        h_s(i,j,k) = h_s_old(i,j,k-1) + beta_k*(h_s_old(ii+1,jj  ,k-1)
                                                              + h_s_old(ii-1,jj  ,k-1)
                                                              + h_s_old(ii  ,jj+1,k-1)
                                                              + h_s_old(ii  ,jj-1,k-1) - 4*h_s_old(ii,jj,k-1));
                    }
                    else {
                        h_s(i,j,k) = h_s_old(i,j,k  ) + beta_k*(h_s_old(ii+1,jj  ,k  )
                                                              + h_s_old(ii-1,jj  ,k  )
                                                              + h_s_old(ii  ,jj+1,k  )
                                                              + h_s_old(ii  ,jj-1,k  ) - 4*h_s_old(ii,jj,k  ));
                    }

                    // Minimum vertical grid spacing condition (Klemp2011 Eqn. 7)
                    return { (zz + A * h_s(i,j,k) - (zz_minus + A_minus * h_s(i,j,k-1))) / (zz - zz_minus) };

                }); //ParReduce

                ParallelDescriptor::ReduceRealMin(diff);

                //
                // Each pass smooths h_s at THIS level while the level below is already
                // fixed, so over a peak it always thins the layer: the passes have to
                // stop before Eqn. 7 is violated, not after.  Measuring only once a pass
                // has been committed is not enough -- nothing then stops the pass that
                // breaks the criterion from being the one we keep, which is how a column
                // over steep terrain ended up at 0.21 of its nominal layer thickness
                // against gamma_m = 0.5.  Reject that pass and restore what preceded it.
                //
                // At iter 0 there is no earlier pass at this level to fall back to, so
                // carry the level below up unsmoothed instead.  That leaves only the
                // A-decay term, whose worst case is exactly what z_H = 2.44 h_m/(1-gamma_m)
                // was chosen to bound: max|dA/dz| = 2.439/z_H, so h_m max|dA/dz| dz equals
                // (1-gamma_m) dz and the thickness ratio lands on gamma_m.
                //
                if (diff < threshold) {
                    ++n_capped;
                    if (iter == 0) {
                        for ( MFIter mfi(h_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
                        {
                            Box xybx = mfi.growntilebox(ngrow);
                            xybx.setRange(2,0);
                            Array4<Real> const& h_s = h_mf.array(mfi);
                            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                            {
                                h_s(i,j,k) = h_s(i,j,k-1);
                            });
                        }
                        MultiFab::Copy(h_mf_old, h_mf,0,0,1,h_mf_old.nGrow());
                        h_mf_old.FillBoundary(geom.periodicity());
                    } else {
                        // h_mf_old still holds the last pass that satisfied the criterion
                        MultiFab::Copy(h_mf, h_mf_old,0,0,1,h_mf.nGrow());
                    }
                    break;
                }

                MultiFab::Copy(h_mf_old, h_mf,0,0,1,h_mf_old.nGrow());

                iter++;

                // Fill ghost points
                h_mf_old.FillBoundary(geom.periodicity());

            } // while

            auto const& z_lev_d = z_levels_d.data();

            // Populate z_phys_nd by solving z_arr(i,j,k) = z + A*h_s(i,j,k)
            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Grown box with no z range
                Box xybx = mfi.growntilebox(ngrow);
                xybx.setRange(2,0);

                Array4<Real> const& h_s   = h_mf_old.array(mfi);
                Array4<Real> const& z_arr = z_phys_nd.array(mfi);

                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

                    // Location of nodes
                    Real z = z_lev_d[k];

                    // STF model from p2164 of Klemp2011 (Eqn. 4)
                    z_arr(i,j,k) = z + A*h_s(i,j,k);

                    // Fill below the bottom surface
                    if (k == 1) {
                        z_arr(i,j,k0-1) = two*z_arr(i,j,k0) - z_arr(i,j,k);
                    }
                });
            } // mfi
        } // k

            if (n_capped > 0) {
                amrex::Print() << "NOTE: the smoothed terrain-following grid at level " << lev
                               << " hit its minimum layer-thickness floor (terrain.gamma_m = "
                               << gamma_m << ") at " << n_capped << " of "
                               << (domhi_z - domlo_z) << " levels, so the terrain smoothing was "
                               << "capped there.  Lowering terrain.gamma_m, coarsening the "
                               << "vertical grid or using terrain.smoothing = 0 would give "
                               << "this terrain more room.\n";
            }
            if (n_capped_out != nullptr) { *n_capped_out = n_capped; }

            Gpu::streamSynchronize();

            break;
        } // case 1

        case 2: // Sullivan TF Method
        {
            int k0    = 0;

            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Grown box with corrected ghost cells at top
                Box gbx = mfi.growntilebox(ngrow);
                gbx.setRange(2,domlo_z,domhi_z+1);

                Array4<Real> const& z_arr = z_phys_nd.array(mfi);
                auto const&         z_lev = z_levels_d.data();

                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Vertical grid stretching
                    Real z =  z_lev[k];

                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    // Fill levels using model from Sullivan et. al. 2014
                    int omega = 3; //Used to adjust how rapidly grid lines level out. omega=1 is BTF!
                    z_arr(i,j,k) = z + static_cast<Real>(std::pow((one - (z/z_top)),omega) * z_arr(ii,jj,k0));

                    // Fill lateral boundaries and below the bottom surface
                    if (k == k0) {
                        z_arr(i,j,k0  ) = z_arr(ii,jj,k0);
                    }
                    if (k == 1) {
                        z_arr(i,j,k0-1) = two*z_arr(ii,jj,k0) - z_arr(i,j,k);
                    }
                });
            } // mfi
            break;
        } // case 2

        case 3: // Debugging Test Method -- applies Sullivan TF starting at k = 1 so that domain does not change size
        {
            int k0    = 0;

            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Grown box with corrected ghost cells at top
                Box gbx = mfi.growntilebox(ngrow);
                gbx.setRange(2,domlo_z,domhi_z+1);

                Array4<Real> const& z_arr = z_phys_nd.array(mfi);
                auto const&         z_lev = z_levels_d.data();

                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Vertical grid stretching
                    Real z =  z_lev[k];

                    int ii = amrex::max(amrex::min(i,imax),imin);
                    int jj = amrex::max(amrex::min(j,jmax),jmin);

                    // Fill values outside the lateral boundaries and below the bottom surface (necessary if init_type = "WRFInput")
                    if (k == k0+1) {
                        z_arr(i,j,k) = z + z_arr(ii,jj,k0);
                    } else {
                        // Fill levels using model from Sullivan et. al. 2014
                        int omega = 3; //Used to adjust how rapidly grid lines level out. omega=1 is BTF!
                        z_arr(i,j,k) = z + static_cast<Real>(std::pow((one - (z/z_top)),omega) * z_arr(ii,jj,k0));
                    }
                });
                gbx.setBig(2,0);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    z_arr(i,j,k  ) = zero;
                    z_arr(i,j,k-1) = -z_arr(i,j,k+1);
                });
            } // mfi
            break;
        } // case 3
    } //switch
}

/**
 * Computation of detJ at cell-center
 */
void
make_J (const Geometry& geom,
        MultiFab& z_phys_nd,
        MultiFab& detJ_cc)
{
    const auto *dx = geom.CellSize();
    Real dzInv = one/dx[2];

    // Domain valid box (z_nd is nodal)
    const Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow= detJ_cc.nGrow();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(detJ_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ngrow);
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }
        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > detJ = detJ_cc.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           detJ(i, j, k) = Real(.25) * dzInv * (
                   z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
                  -z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i+1,j+1,k  ) );
        });
    }
    detJ_cc.FillBoundary(geom.periodicity());
}

/**
 * Computation of area fractions on faces
 */
void
make_areas (const Geometry& geom,
            MultiFab& z_phys_nd,
            MultiFab& ax,
            MultiFab& ay,
            MultiFab& az)
{
    const auto* dx = geom.CellSize();
    Real dzInv = one/dx[2];

    // Domain valid box (z_nd is nodal)
    const Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // The z-faces are always full when using terrain-fitted coordinates
    az.setVal(1.0);

    //
    // x-areas
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(ax, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ax.nGrow());
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }

        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > ax_arr = ax.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               ax_arr(i, j, k) = Real(.5) * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i,j+1,k+1) - z_nd(i,j,k) - z_nd(i,j+1,k));
        });
    }

    //
    // y-areas
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(ay, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ay.nGrow());
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }

        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > ay_arr = ay.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               ay_arr(i, j, k) = Real(.5) * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i+1,j,k+1) - z_nd(i,j,k) - z_nd(i+1,j,k));
        });
    }

    ax.FillBoundary(geom.periodicity());
    ay.FillBoundary(geom.periodicity());
    az.FillBoundary(geom.periodicity());
}

/**
 * Computation of z_phys at cell-center
 */
void
make_zcc (const Geometry& geom,
          MultiFab& z_phys_nd,
          MultiFab& z_phys_cc)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(z_phys_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox();
        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > z_cc = z_phys_cc.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           z_cc(i, j, k) = Real(.125) * ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
                                         +z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) );
       });
    }
    z_phys_cc.FillBoundary(geom.periodicity());
}


/**
 * Computation min dz at cell-center
 */
Real
get_dzmin_terrain (MultiFab& z_phys_nd)
{
    auto const& ma_z_nd_arr = z_phys_nd.const_arrays();
    GpuTuple<Real> min = ParReduce(TypeList<ReduceOpMin>{},
                                   TypeList<Real>{},
                                   z_phys_nd, IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real>
                                   {
                                       amrex::Real dz = Compute_Z_AtWFace(i,j,k+1,ma_z_nd_arr[box_no]) -
                                                        Compute_Z_AtWFace(i,j,k  ,ma_z_nd_arr[box_no]);
                                       return { dz };
                                   });
    Real r = (get<0>(min) + std::numeric_limits<amrex::Real>::epsilon());
    ParallelDescriptor::ReduceRealMin(r);
    return r;
}
