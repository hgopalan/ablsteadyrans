#!/usr/bin/env python3
"""
Sub-phase 3.5: the one-line mutants of the SIMPLEC <-> Newton switching
(Source/NewtonStage.cpp and the switch rules in Source/Driver.cpp).

Run against Driver.* and NewtonStage.* by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/newton_stage_mutants.py BUILD_DIR [N1 N4 ...]
"""
import mutants

D = "Source/Driver.cpp"
N = "Source/NewtonStage.cpp"
M = [
 ("N1", "the fallback waits for the hysteresis window", D,
  "    if (m_h.stage > 0 && m_stages[m_h.stage]->fallback_requested()) {",
  "    if (m_h.stage > 0 && m_stages[m_h.stage]->fallback_requested() && !(m_h.last_switch >= 0 && iteration - m_h.last_switch < W)) {"),
 ("N2", "the forward switch ignores the monitors", D,
  "               && m_stages[m_h.stage + 1]->accepts_switch() && monitors_settled()) {",
  "               && m_stages[m_h.stage + 1]->accepts_switch()) {"),
 ("N3", "the forward switch ignores accepts_switch", D,
  "               && m_stages[m_h.stage + 1]->accepts_switch() && monitors_settled()) {",
  "               && monitors_settled()) {"),
 ("N4", "leave not called on a switch", D,
  "    m_stages[e.from]->leave(q);", ""),
 ("N5", "the monitor test against the window's first record off by one", D,
  "    const auto& then = rec[rec.size() - 1 - static_cast<std::size_t>(W)].monitors;",
  "    const auto& then = rec[rec.size() - static_cast<std::size_t>(W)].monitors;"),
 ("N6", "the fallback after one failure whatever the key", N,
  "    if (m_failures >= m_sopt.newton_fallback_failures) {", "    if (m_failures >= 1) {"),
 ("N7", "a failed step does not count", N,
  "    m_failures = (m_status == ContinuationStatus::failed) ? m_failures + 1 : 0;",
  "    m_failures = 0;"),
 ("N8", "converged before any accepted step", N,
  "    if (before.steps > 0 && before.r0 > 0 && before.r_last <= drop * before.r0) {",
  "    if (before.r0 > 0 && before.r_last <= drop * before.r0) {"),
 ("N9", "no hand-back for the final pass", N,
  "        m_handed_back = true;\n        char b[64];", "        char b[64];"),
 ("N10", "the continuation's start CFL not reset on entry", N,
  "        st.cfl = m_sopt.cfl_max;", ""),
 ("N11", "the continuation state not restored on restart", N,
  "        m_cont.restore(read_continuation_state(dir));", ""),
 ("N12", "no limit on failed visits", N,
  "        if (m_failed_visits >= kMaxFailedVisits) {", "        if (false) {"),
 ("N13", "no idle hand-back", N,
  "        if (++m_idle > m_sopt.criterion.monitor_window) {", "        if (false) {"),
 ("N14", "the Newton cap not used (the SIMPLEC cap, plain Newton after one step)", N,
  "    o.cfl0 = sopt.newton_cfl_max;\n    o.cfl_max = sopt.newton_cfl_max;", "    o.cfl0 = sopt.cfl_max;\n    o.cfl_max = sopt.cfl_max;"),
 ("N15", "the preconditioner of an infinite dtau at the SIMPLEC cap", N,
  "z_nd, sopt.newton_cfl_max, ba, dm, line_sweep_options(sopt));", "z_nd, sopt.cfl_max, ba, dm, line_sweep_options(sopt));"),
]

if __name__ == "__main__":
    mutants.run(M, "Driver.*:NewtonStage.*", {}, __doc__)
