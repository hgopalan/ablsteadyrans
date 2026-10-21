#ifndef GUST_REDUCTIONS_H_
#define GUST_REDUCTIONS_H_

#include <AMReX_MultiFab.H>

namespace gust {

//! Sum of component comp (cell-centred data only) over the valid region, identical to the last bit for
//! any number of MPI ranks and any DistributionMapping of a fixed BoxArray.
//! A different BoxArray (e.g. after a restart onto another layout) can change
//! the last bits; anything compared bitwise across runs must keep the layout.
//!
//! Each column (i, j) of a box is summed in k order by one GPU thread; the
//! column sums are added in storage order (i fastest, then j) on the host;
//! each box sum lands in its own slot of a vector of length nboxes, which is
//! all-reduced (every slot has one nonzero contributor, so adding the zeros
//! is exact) and then added in box-index order. The result does depend on
//! the BoxArray: a different box layout adds in a different order.
amrex::Real deterministic_sum (const amrex::MultiFab& mf, int comp);

//! Maximum of component comp over the valid region, for any data type. NaN
//! if any valid value is NaN; zero is returned as +0. With those two cases
//! fixed, the maximum does not depend on the order of evaluation.
amrex::Real global_max (const amrex::MultiFab& mf, int comp);

} // namespace gust

#endif
