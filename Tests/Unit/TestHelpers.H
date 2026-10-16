#ifndef ABLSTEADY_TEST_HELPERS_H_
#define ABLSTEADY_TEST_HELPERS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>

#include <cstring>
#include <limits>

namespace ablsteady::test {

//! Number of points in the region grown by ng where a and b differ in any
//! bit (a and b must share BoxArray and DistributionMapping). Points are
//! compared with memcmp, so a NaN equals the same NaN and -0 differs from +0.
inline long count_bit_differences (const amrex::MultiFab& a, const amrex::MultiFab& b, int ng)
{
    long ndiff = 0;
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(ng);
        amrex::FArrayBox ha(bx, a.nComp(), amrex::The_Pinned_Arena());
        amrex::FArrayBox hb(bx, b.nComp(), amrex::The_Pinned_Arena());
        ha.copy<amrex::RunOn::Device>(a[mfi], bx);
        hb.copy<amrex::RunOn::Device>(b[mfi], bx);
        amrex::Gpu::streamSynchronize();
        const auto n = static_cast<std::size_t>(bx.numPts()) * a.nComp();
        for (std::size_t m = 0; m < n; ++m) {
            if (std::memcmp(ha.dataPtr() + m, hb.dataPtr() + m, sizeof(amrex::Real)) != 0) {
                ++ndiff;
            }
        }
    }
    amrex::ParallelAllReduce::Sum(ndiff, amrex::ParallelDescriptor::Communicator());
    return ndiff;
}

//! A value mant 2^e with mant < 2^20 and e in [-40, -1] (magnitudes up to 2^19),
//! built from integers only (no libm), so it is the same on every platform
//! and sums of it round differently in different orders.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pattern (int i, int j, int k, int c = 0)
{
    auto h = static_cast<unsigned int>(i) * 73856093u
           ^ static_cast<unsigned int>(j) * 19349663u
           ^ static_cast<unsigned int>(k) * 83492791u
           ^ static_cast<unsigned int>(c) * 2654435761u;
    h ^= h >> 13;
    h *= 0x5bd1e995u;
    h ^= h >> 15;
    const int mant = static_cast<int>(h & 0xFFFFFu);
    const int expo = static_cast<int>((h >> 20) % 40u) - 40;
    double v = static_cast<double>(mant);
    for (int e = 0; e > expo; --e) { v *= 0.5; }    // exact: v = mant 2^expo
    return static_cast<amrex::Real>((h >> 31) ? -v : v);
}

//! Fill component 0 of mf, ghost cells included, with pattern(i, j, k, c).
//! A free function: device lambdas must not sit inside TEST bodies.
inline void fill_pattern (amrex::MultiFab& mf, int c)
{
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox();
        const auto a = mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) = pattern(i, j, k, c);
        });
    }
}

//! Set the valid point iv of component 0 to a quiet NaN (no arithmetic on it,
//! so no floating-point trap).
inline void set_one_nan (amrex::MultiFab& mf, const amrex::IntVect& iv)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (!mfi.validbox().contains(iv)) { continue; }
        const auto a = mf.array(mfi);
        const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::ParallelFor(amrex::Box(iv, iv), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) = nan;
        });
    }
}

} // namespace ablsteady::test

#endif
