#ifndef GUST_CLOSUREWAKE_H_
#define GUST_CLOSUREWAKE_H_

// The realizability cap of the k-L closure for separated and wake flows:
//
//   nu_t <= cap,  cap = C_r k / S   (turb.realizability = C_r > 0; S = sqrt(2 S_ij S_ij)),
//
// the form of the SST limiter (a_1 = 0.31). In neutral flow nu_t S / k = Cmu0^2 sqrt(P / eps), so C_r = 0.309 (just
// below Cmu0^2 = 0.30936 for the default Cmu0) trims nu_t by 0-3 % in an equilibrium surface layer and limits it
// strongly where the strain far exceeds equilibrium: the separating shear layer behind a crest
// (docs/turbulence_model.rst). Without the option the cap is the
// largest Real, and min(nu_t, cap) is nu_t exactly (bit-identical to the plain closure).

#include "Closure.H"

#include <AMReX_MultiFab.H>

namespace gust {

//! The cap on the closure's cells (the valid cells and one lateral ghost ring; defined here if empty) from
//! k (q.tke) and S2 = S^2.
void nut_cap_field (const State& q, const amrex::MultiFab& S2, const ClosureOptions& opt, const ClosureConstants& c,
                    amrex::MultiFab& cap);

} // namespace gust

#endif
