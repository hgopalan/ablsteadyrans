#!/usr/bin/env python3
"""
Sub-phase 2.2: the one-line mutants of the box-local line sweeps
(Source/LineSweeps.cpp) and of their input reader (Source/Driver.cpp).

Each mutant replaces one line (the old text must occur exactly once), rebuilds
unit_tests, runs LineSweeps.*, Driver.* and Predictor.Lusgs* on two ranks,
restores the file and prints one row: id | change | killed or SURVIVED |
failing tests. The table is in docs/testing.rst.

    python3 tools/line_sweep_mutants.py BUILD_DIR [M1 M4 ...]
"""
import os
import re
import subprocess
import sys

M = [
 ("M1", "overlap ignored (swept region not grown)", "Source/LineSweeps.cpp",
  "line_sweep_box(loc, line_owned_box(loc, mfi.validbox(), geom), overlap, geom);",
  "line_sweep_box(loc, line_owned_box(loc, mfi.validbox(), geom), 0, geom);"),
 ("M2", "face updated twice (box keeps its high face)", "Source/LineSweeps.cpp",
  "if (vbox.bigEnd(d) != dom.bigEnd(d) + 1) { owned.growHi(d, -1); }",
  "if (false) { owned.growHi(d, -1); }"),
 ("M3", "periodic seam updated twice (first box keeps the domain's low face)", "Source/LineSweeps.cpp",
  "if (geom.isPeriodic(d) && vbox.smallEnd(d) == dom.smallEnd(d)) { owned.growLo(d, -1); }",
  "if (false) { owned.growLo(d, -1); }"),
 ("M4", "forward wavefront skips its last plane", "Source/LineSweeps.cpp",
  "for (int s = 0; s <= last; ++s) { plane(s, S, a, dg, r, d, c, e); }",
  "for (int s = 0; s < last; ++s) { plane(s, S, a, dg, r, d, c, e); }"),
 ("M5", "backward wavefront skips plane 0", "Source/LineSweeps.cpp",
  "for (int s = last; s >= 0; --s) { plane(s, S, a, dg, r, d, c, e); }",
  "for (int s = last; s > 0; --s) { plane(s, S, a, dg, r, d, c, e); }"),
 ("M6", "ghost delta not zeroed (neighbours' values exchanged between passes)", "Source/LineSweeps.cpp",
  "            const PassFn pass = pass_of(order, forward, alternate && s % 2 == 1);\n            for (amrex::MFIter mfi(work, false);",
  "            const PassFn pass = pass_of(order, forward, alternate && s % 2 == 1); work.FillBoundary(geom.periodicity());\n            for (amrex::MFIter mfi(work, false);"),
 ("M7", "non-owned face copies keep their own value (owner mask all 1)", "Source/LineSweeps.cpp",
  "m(i, j, k) = owned.contains(i, j, k) ? 1 : 0;", "m(i, j, k) = 1;"),
 ("M8", "one-cell periodic axis not folded into the diagonal", "Source/LineSweeps.cpp",
  "if (fold[d] == 0) { continue; }", "if (true) { continue; }"),
 ("M9", "steady.line_overlap not read", "Source/Driver.cpp",
  'pp.query("line_overlap", opt.line_overlap);', ';'),
 ("M10", "steady.line_ghosts not read", "Source/Driver.cpp",
  'if (pp.query("line_ghosts", word)) {', 'if (false && pp.query("line_ghosts", word)) {'),
 ("M11", "steady.line_order not read", "Source/Driver.cpp",
  'if (pp.query("line_order", word)) {', 'if (false && pp.query("line_order", word)) {'),
 ("M12", "line_overlap above 8 accepted", "Source/LineSweeps.cpp",
  "if (opt.overlap < 0 || opt.overlap > max_line_overlap) {", "if (opt.overlap < 0 || opt.overlap > max_line_overlap + 1) {"),
 ("M13", "negative line_overlap accepted", "Source/LineSweeps.cpp",
  "if (opt.overlap < 0 || opt.overlap > max_line_overlap) {", "if (opt.overlap < -1 || opt.overlap > max_line_overlap) {"),
 ("M14", "overlap with lagged ghosts accepted", "Source/LineSweeps.cpp",
  "} else if (opt.overlap > 0 && opt.ghosts != LineGhosts::zero) {", "} else if (false) {"),
 ("M15", "overlap wider than the narrowest box accepted", "Source/LineSweeps.cpp",
  "narrowest < std::numeric_limits<int>::max() && overlap > narrowest) {", "narrowest < std::numeric_limits<int>::max() && overlap > narrowest + 1) {"),
 ("M16", "line key given with zjacobi accepted", "Source/Driver.cpp",
  "if (opt.sweep != SweepKind::lusgs) {\n        for (const char* key", "if (false) {\n        for (const char* key"),
 # beyond the plan: the alternating ordering (steady.line_alternate)
 ("M17", "odd sweeps not reversed (alternation ignored)", "Source/LineSweeps.cpp",
  "    if (jrev) {", "    if (false) {"),
 ("M18", "reversed forward wavefront skips its last plane", "Source/LineSweeps.cpp",
  "for (int s = 0; s <= last; ++s) { plane_jrev(s, S, a, dg, r, d, c, e); }",
  "for (int s = 0; s < last; ++s) { plane_jrev(s, S, a, dg, r, d, c, e); }"),
 ("M19", "reversed sequential forward pass walks j ascending", "Source/LineSweeps.cpp",
  "        for (int j = jhi; j >= jlo; --j) {\n            for (int i = ilo; i <= ihi; ++i)",
  "        for (int j = jlo; j <= jhi; ++j) {\n            for (int i = ilo; i <= ihi; ++i)"),
 ("M20", "steady.line_alternate not read", "Source/Driver.cpp",
  'pp.query("line_alternate", opt.line_alternate);', ';'),
 ("M21", "line_alternate = 2 accepted", "Source/LineSweeps.cpp",
  "if (opt.alternate != 0 && opt.alternate != 1) {", "if (opt.alternate < 0 || opt.alternate > 2) {"),
 ("M22", "line_alternate given with zjacobi accepted", "Source/Driver.cpp",
  '{"line_ghosts", "line_overlap", "line_order", "line_alternate"}', '{"line_ghosts", "line_overlap", "line_order"}'),
]


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    build = os.path.abspath(sys.argv[1])
    only = sys.argv[2:]
    survived = 0
    for mid, what, f, old, new in M:
        if only and mid not in only:
            continue
        path = os.path.join(root, f)
        with open(path) as fh:
            src = fh.read()
        if src.count(old) != 1:
            sys.exit(f"{mid}: the original line occurs {src.count(old)} times in {f}")
        with open(path, "w") as fh:
            fh.write(src.replace(old, new))
        try:
            b = subprocess.run(["cmake", "--build", build, "--target", "unit_tests"], capture_output=True, text=True)
            failed = []
            if b.returncode:
                res = "build error"
            else:
                t = subprocess.run(["mpirun", "--oversubscribe", "-np", "2", os.path.join(build, "unit_tests"),
                                    "--gtest_filter=LineSweeps.*:Driver.*:Predictor.Lusgs*",
                                    "amrex.the_arena_init_size=0"], capture_output=True, text=True, timeout=900)
                failed = sorted(set(re.findall(r"\[  FAILED  \] (\S+\.\S+)", t.stdout)))
                res = "killed" if failed or t.returncode else "SURVIVED"
        finally:
            with open(path, "w") as fh:
                fh.write(src)
        survived += res == "SURVIVED"
        print(f"{mid} | {what} | {res} | {', '.join(failed)}", flush=True)
    subprocess.run(["cmake", "--build", build, "--target", "unit_tests"], capture_output=True)
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
