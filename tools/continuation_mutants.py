#!/usr/bin/env python3
"""
The one-line mutants of the pseudo-transient continuation
(Source/Continuation.cpp, the SER law in Source/PseudoTime.cpp and the GMRES
limit in Source/Newton.cpp).

Run against Continuation.*, Newton.* and Driver.* (the driver also calls
ser_cfl) by the shared runner tools/mutants.py (two ranks, the unit-test
flags). The table is in docs/testing.rst.

    python3 tools/continuation_mutants.py BUILD_DIR [C1 C4 ...]
"""
import glob
import os
import re
import sys

import mutants

C = "Source/Continuation.cpp"
P = "Source/PseudoTime.cpp"
N = "Source/Newton.cpp"
M = [
 # the line search
 ("C1", "Armijo without the sufficient decrease (any decrease accepted)", C,
  "r_trial <= (amrex::Real(1) - c * lambda) * r;", "r_trial <= r;"),
 ("C2", "Armijo compares a NaN trial (no isfinite guard)", C,
  "return std::isfinite(static_cast<double>(r_trial)) && r_trial <=", "return r_trial <="),
 ("C3", "backtracking factor fixed at 1/2", C,
  "const amrex::Real next = opt.backtrack * ls.lambda;", "const amrex::Real next = amrex::Real(0.5) * ls.lambda;"),
 ("C4", "no lower bound on lambda", C,
  "if (!(next >= opt.lambda_min)) { break; }", "if (!(next > amrex::Real(0))) { break; }"),
 ("C5", "backtracks not counted", C, "++ls.backtracks;", ";"),
 ("C6", "full-step residual not kept (the growth bound never passes)", C,
  "if (ls.lambda == 1) { r_full = ls.residual; }", ";"),
 ("C7", "growth bound strict at the bound", C,
  "r_full <= growth * r) {", "r_full < growth * r) {"),
 ("C8", "growth-bound step not re-evaluated (last trial kept)", C,
  "        ls.lambda = 1;\n        ls.residual = phi(ls.lambda);\n        ls.accepted = true;",
  "        ls.lambda = 1;\n        ls.accepted = true;"),
 ("C9", "growth-bound flag not recorded", C, "ls.by_growth = true;", ";"),
 # the SER update and the CFL cut
 ("C10", "SER restarts from the stored CFL after infinite dtau", C,
  "const amrex::Real from = st.infinite ? opt.cfl_max : st.cfl;", "const amrex::Real from = st.cfl;"),
 ("C11", "never switches to infinite dtau", C,
  "st.infinite = opt.dtau_inf_at_cap == 1 && st.cfl >= opt.cfl_max;", "st.infinite = false;"),
 ("C12", "last residual not stored", C, "st.r_last = r_new;", ";"),
 ("C13", "accepted steps not counted", C, "st.steps += 1;", ";"),
 ("C14", "no CFL cut after a failed attempt", C,
  "opt.cfl_cut * (st.infinite ? opt.cfl_max : st.cfl)", "(st.infinite ? opt.cfl_max : st.cfl)"),
 ("C15", "no floor on the cut CFL", C, "std::numeric_limits<amrex::Real>::min());", "amrex::Real(0));"),
 ("C16", "a failed attempt keeps infinite dtau", C,
  "    st.infinite = false;\n    return st;", "    return st;"),
 ("C17", "growth bound also at infinite dtau", C,
  "return st.infinite ? amrex::Real(1) : opt.ptc_growth;", "return opt.ptc_growth;"),
 # restart
 ("C18", "checkpoint drops the infinite flag", C, "(st.infinite ? 1 : 0)", "0"),
 ("C19", "reader ignores r_last", C, "st.r_last = real(\"r_last\");", "(void)real(\"r_last\");"),
 ("C20", "restore keeps a CFL above a lowered cfl_max", C,
  "m_st.cfl = amrex::min(m_st.cfl, m_opt.cfl_max);", ";"),
 # reference residual and step
 ("C21", "set_reference_residual ignores its argument", C, "    m_st.r0 = r0;\n", "\n"),
 ("C22", "set_reference_residual accepts 0", C,
  "|| !(d > 0)) {\n        amrex::Abort(\"gust: Continuation::set_reference_residual",
  "|| !(d >= 0)) {\n        amrex::Abort(\"gust: Continuation::set_reference_residual"),
 ("C23", "no layout check in step", C,
  "if (q.p.boxArray() != m_qb.p.boxArray() || q.p.DistributionMap() != m_qb.p.DistributionMap()) {",
  "if (false) {"),
 ("C24", "|R_0| never taken from the first step", C,
  "if (m_st.r0 == 0 && finite_r) { m_st.r0 = info.residual_norm; }", ";"),
 ("C25", "trial k not clipped at the floor", C, "clip_floor(m_qt.tke, m_tke_floor);", ";"),
 ("C26", "accepted step not copied back to q", C, "copy_state(q, m_qt);", ";"),
 ("C27", "finite dtau even when infinite", C, "if (!m_st.infinite) {", "if (true) {"),
 ("C28", "one retry fewer", C,
  "attempt <= m_opt.max_retries; ++attempt", "attempt < m_opt.max_retries; ++attempt"),
 ("C29", "run() never stops at the tolerance inside the loop", C,
  "if (m_st.steps > 0 && m_st.r_last <= rel_tol * m_st.r0) {", "if (false) {"),
 # inputs
 ("C30", "newton_armijo_c accepted up to 1", C,
  "open01(opt.armijo_c, \"newton_armijo_c\", 0.5);", "open01(opt.armijo_c, \"newton_armijo_c\", 1);"),
 ("C31", "newton_max_retries = 0 accepted", C, "if (opt.max_retries < 1) {", "if (opt.max_retries < 0) {"),
 ("C32", "newton_ptc_growth below 1 accepted", C,
  "if (std::isnan(d) || !std::isfinite(d) || !(d >= 1)) {", "if (std::isnan(d) || !std::isfinite(d) || !(d >= 0)) {"),
 ("C33", "newton_ptc_growth not read", C, "pp.query(\"newton_ptc_growth\", opt.ptc_growth);", ";"),
 ("C34", "newton_dtau_inf_at_cap not read", C, "pp.query(\"newton_dtau_inf_at_cap\", opt.dtau_inf_at_cap);", ";"),
 ("C35", "a repeated key accepted", C,
  "if (pp.countname(key) > 1) { err.emplace_back", "if (false) { err.emplace_back"),
 ("C36", "tke_floor not checked", C,
  "if (std::isnan(static_cast<double>(tke_floor)) || !std::isfinite(static_cast<double>(tke_floor)) || !(tke_floor >= 0)) {",
  "if (false) {"),
 # the SER law (also the driver's) and the GMRES limit
 ("C37", "SER without the growth limit", P,
  "const amrex::Real capped = amrex::min(cfl_max, growth * cfl);", "const amrex::Real capped = cfl_max;"),
 ("C38", "SER exponent sign flipped", P,
  "return amrex::min(capped, cfl0 * std::pow(r, -p));", "return amrex::min(capped, cfl0 * std::pow(r, p));"),
 ("C39", "GMRES iteration limit not applied", N, "gmres.setMaxIters(m_opt.gmres_max_iters);", ";"),
 ("C40", "newton_gmres_max_iters = 0 accepted", N,
  "if (opt.gmres_max_iters < 1) {", "if (opt.gmres_max_iters < 0) {"),
]

# C2 (no isfinite guard) is told apart only through the invalid flag that
# ArmijoOnKnownQuadratic checks after comparing a NaN. GCC compiles the
# ordered compare to a signalling instruction (comisd on x86-64, fcmpe on
# aarch64) and kills C2; Clang and AppleClang emit a quiet compare
# (ucomisd, fcmp), so no test can tell C2 from the original there.
def compiler_id(build):
    """CMAKE_CXX_COMPILER_ID of the build (CMake records it per version)."""
    for path in glob.glob(os.path.join(build, "CMakeFiles", "*", "CMakeCXXCompiler.cmake")):
        with open(path) as fh:
            m = re.search(r'^set\(CMAKE_CXX_COMPILER_ID "([^"]*)"\)', fh.read(), re.M)
        if m:
            return m.group(1)
    return ""

EQUIVALENT = {}
if len(sys.argv) > 1 and "Clang" in compiler_id(sys.argv[1]):
    EQUIVALENT["C2"] = "quiet NaN compare with Clang"

if __name__ == "__main__":
    mutants.run(M, "Continuation.*:Newton.*:Driver.*", EQUIVALENT, __doc__)
