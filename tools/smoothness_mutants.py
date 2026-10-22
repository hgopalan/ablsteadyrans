#!/usr/bin/env python3
"""
The one-line mutants of the smoothness handling
(Source/Smoothness.cpp) and of the closure's floor freeze and stratification
gate (Source/Closure.cpp).

Run against Smoothness.* and Closure.* by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/smoothness_mutants.py BUILD_DIR [S1 S4 ...]
"""
import mutants

C = "Source/Closure.cpp"
S = "Source/Smoothness.cpp"
M = [
 ("S1", "the floor scope holds every decision", C,
  "return scope == FreezeScope::floor ? kFloorBit : kAllBits;", "return kAllBits;"),
 ("S2", "changed_switches compares every bit", C,
  "return {(now & fm) != (sw(i, j, k) & fm) ?", "return {now != sw(i, j, k) ?"),
 ("S3", "the floor scope decides from the live state (frozen bits ignored)", C,
  "        floor_held_decisions(q.tke, N2, dist, z0, m_c, m_switches, out);", "        live_decisions(q.tke, N2, dist, z0, m_c, out);"),
 ("S11", "a held cell decided at its own k, not the floor (negative nu_t)", C,
  "const amrex::Real k_in = ((sw(i, j, k) & kFloorBit) != 0) ? c.tke_floor : kk(i, j, k);",
  "const amrex::Real k_in = kk(i, j, k);"),
 ("S12", "the floor scope holds every decision in evaluate", C,
  "    } else if (m_scope == FreezeScope::all) {", "    } else if (true) {"),
 ("S13", "P_b from the gated N^2", C,
  "Pb = -ratio * nt * n2r(i, j, k);", "Pb = -ratio * nt * n2(i, j, k);"),
 ("S4", "the gate with S^2 not squared", C,
  "b = s2(i, j, k) * s2(i, j, k);", "b = s2(i, j, k);"),
 ("S5", "the gate not applied", C,
  "if (m_opt.ri_gate > 0) { gate_stratification(N2, S2, m_opt.ri_gate); }", ""),
 ("S6", "turb.ri_gate not read", C,
  '    pp.query("ri_gate", opt.ri_gate);\n', ""),
 ("S7", "no limiter swap", S,
  "return (l == Limiter::van_leer) ? Limiter::van_albada : l;", "return l;"),
 ("S8", "leave does not thaw", S,
  "        closure->thaw();\n    }\n    m_active = false;", "    }\n    m_active = false;"),
 ("S9", "enter freezes every decision", S,
  "geom, FreezeScope::floor);", "geom, FreezeScope::all);"),
 ("S10", "the checkpoint's scope ignored on read", S,
  'scope == "all" ? FreezeScope::all : FreezeScope::floor', "FreezeScope::all"),
]

if __name__ == "__main__":
    # S5 removes the gate from Closure::stratification, which only freeze_switches and changed_switches use (evaluate
    # gates N^2 itself): with the floor scope the held bit does not depend on N^2.
    mutants.run(M, "Smoothness.*:Closure.*", {"S5": "the floor bit does not depend on N^2"}, __doc__)
