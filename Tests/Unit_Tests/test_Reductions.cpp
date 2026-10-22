// Deterministic sum and max, identical to the last bit on 1, 2 and
// 4 ranks. CTest runs this file's tests at each rank count; every run must
// reproduce the same pinned bits.
#include "Reductions.H"
#include "State.H"
#include "TestHelpers.H"

#include <gtest/gtest.h>

#include <cstdint>
#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>

namespace {

using gust::test::fill_pattern;
using gust::test::pattern;

// 24 x 20 x 16 cells cut 8 x 8 in x/y: 3 x 3 boxes (the y boxes are 8, 8, 4).
gust::MeshInputs reduction_mesh ()
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(24, 20, 16);
    m.prob_hi = {1, 1, 1};
    m.max_grid_size = {8, 8};
    return m;
}

std::uint64_t bits (double x)
{
    std::uint64_t b = 0;
    std::memcpy(&b, &x, sizeof(b));
    return b;
}

// -0x1.f9a008f38d09cp+19, pinned from the 1-rank run (macOS arm64, AppleClang);
// the 2- and 4-rank runs, and Linux GCC/Clang, must reproduce it. The addends
// are dyadic and the only operations are additions, so no FMA contraction or
// libm difference can change it.
constexpr std::uint64_t kSumBits = 0xc12f9a008f38d09cULL;

TEST (Reductions, SumIsIdenticalOnEveryRankCount)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) {
        GTEST_SKIP() << "pinned value is for double precision";
    } else {
        const auto ba = gust::make_box_array(reduction_mesh());
        ASSERT_EQ(ba.size(), 9);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab mf(ba, dm, 1, 1);
        fill_pattern(mf, 3);

        const double s = gust::deterministic_sum(mf, 0);
        EXPECT_EQ(bits(s), kSumBits) << std::hexfloat << s;

        // The same value with a different DistributionMapping of the same boxes.
        amrex::Vector<int> pmap(ba.size());
        for (int b = 0; b < ba.size(); ++b) {
            pmap[b] = (ba.size() - 1 - b) % amrex::ParallelDescriptor::NProcs();
        }
        const amrex::DistributionMapping dm_rev(pmap);
        amrex::MultiFab mf_rev(ba, dm_rev, 1, 1);
        fill_pattern(mf_rev, 3);
        EXPECT_EQ(bits(gust::deterministic_sum(mf_rev, 0)), bits(s));

        // Accuracy: within 16 eps sum|x| of the exact sum, pinned from Python's
        // math.fsum of the same addends (-0x1.f9a008f38d094p+19).
        double abs_sum = 0;
        for (int k = 0; k < 16; ++k) {
            for (int j = 0; j < 20; ++j) {
                for (int i = 0; i < 24; ++i) { abs_sum += std::abs(static_cast<double>(pattern(i, j, k, 3))); }
            }
        }
        const double exact = -0x1.f9a008f38d094p+19;
        EXPECT_LE(std::abs(s - exact), 16 * std::numeric_limits<double>::epsilon() * abs_sum);
    }
}

TEST (Reductions, SumIgnoresGhostCells)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab a(ba, dm, 1, 1), b(ba, dm, 1, 1);
    fill_pattern(a, 3);
    fill_pattern(b, 3);
    b.setBndry(1.0e30);
    EXPECT_EQ(gust::deterministic_sum(a, 0), gust::deterministic_sum(b, 0));
}

TEST (Reductions, MaxEqualsTheExactMaximum)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 1);
    fill_pattern(mf, 5);
    mf.setBndry(1.0e30);   // ghost cells must not count

    amrex::Real ref = std::numeric_limits<amrex::Real>::lowest();
    for (int k = 0; k < 16; ++k) {
        for (int j = 0; j < 20; ++j) {
            for (int i = 0; i < 24; ++i) { ref = std::max(ref, pattern(i, j, k, 5)); }
        }
    }
    EXPECT_EQ(gust::global_max(mf, 0), ref);
}

TEST (Reductions, MaxReportsNaNAndPositiveZero)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 1);

    mf.setVal(-0.0);
    const amrex::Real z = gust::global_max(mf, 0);
    EXPECT_EQ(z, 0);
    EXPECT_FALSE(std::signbit(z));

    mf.setVal(-1.0);
    gust::test::set_one_nan(mf, amrex::IntVect(17, 3, 9));
    EXPECT_TRUE(std::isnan(gust::global_max(mf, 0)));
}

TEST (Reductions, SumRejectsFaceData)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab u(amrex::convert(ba, amrex::IntVect(1, 0, 0)), dm, 1, 0);
    u.setVal(1.0);
    EXPECT_THROW(gust::deterministic_sum(u, 0), std::runtime_error);
}

} // namespace
