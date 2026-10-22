#!/usr/bin/env python3
"""Compare the steady states of two runs of one deck (rt_steady_hill3d_neutral_matches).

    check_line_match.py DIR_A DIR_B PREFIX [--rel_tol 1e-8] [--mean_free p[,...]]

Each directory must hold exactly one plotfile PREFIX<step> (the run's final
state; the per-level plotfiles of a grid hierarchy, PREFIX_level<l>_<step>,
are not counted; the test's setup step empties the directories first, so a stale
plotfile cannot be compared). u, v, w, p and tke must agree to rel_tol of
their largest magnitude (the infinity norm, as amrex_fcompare). theta is
compared after subtracting each run's mean: with no theta condition and no
gravity its additive constant is free and drifts at round-off, so only its
variation is part of the steady state. --mean_free adds fields compared the
same way (p where no face fixes its level and the runs start from different
states, as grid sequencing does).
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

FIELDS = ["u", "v", "w", "p", "tke", "theta"]


def the_plotfile(directory, prefix):
    # PREFIX followed by the step only: a run with a grid hierarchy also
    # writes PREFIX_level<l>_<step> (MeshHierarchy::write_plotfiles).
    found = sorted(p for p in glob.glob(os.path.join(directory, prefix + "*"))
                   if os.path.isdir(p) and os.path.basename(p)[len(prefix):].isdigit())
    if len(found) != 1:
        raise SystemExit("check_line_match: expected one %s* in %s, found %d" % (prefix, directory, len(found)))
    return found[0]


def flat(field):
    return [v for plane in field for col in plane for v in col]


def main(argv):
    rel_tol, mean_free = 1e-8, {"theta"}
    if len(argv) < 4:
        raise SystemExit(__doc__)
    i = 4
    while i < len(argv):
        if argv[i] == "--rel_tol" and i + 1 < len(argv):
            rel_tol = float(argv[i + 1])
        elif argv[i] == "--mean_free" and i + 1 < len(argv) and set(argv[i + 1].split(",")) <= set(FIELDS):
            mean_free |= set(argv[i + 1].split(","))
        else:
            raise SystemExit(__doc__)
        i += 2
    a_plt, b_plt = the_plotfile(argv[1], argv[3]), the_plotfile(argv[2], argv[3])
    _, a = erf_plotfile.read_fields(a_plt, FIELDS)
    _, b = erf_plotfile.read_fields(b_plt, FIELDS)
    ok = True
    for f in FIELDS:
        x, y = flat(a[f]), flat(b[f])
        scale = max(abs(v) for v in x)
        if f in mean_free:
            mx, my = sum(x) / len(x), sum(y) / len(y)
            x = [v - mx for v in x]
            y = [v - my for v in y]
        err = max(abs(p - q) for p, q in zip(x, y))
        good = err <= rel_tol * scale
        ok = ok and good
        label = f + " - mean" if f in mean_free else f
        print("check_line_match: %-12s max |a - b| = %.3e, %.1e x max |a| = %.3e: %s"
              % (label, err, rel_tol, rel_tol * scale, "ok" if good else "FAIL"))
    print("check_line_match: %s vs %s: %s" % (a_plt, b_plt, "match" if ok else "MISMATCH"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
