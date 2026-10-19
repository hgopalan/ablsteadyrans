#!/usr/bin/env python3
"""The converged-column first guess removes the inertial spin-up (sub-phase 2.1b).

    check_column_guess.py PROFILE_CHK COLUMN_CHK

Reads the Driver histories of two steady runs of the same deck, one from the
inflow profile (init.first_guess = profile) and one from the converged column
(init.first_guess = column), and requires that at iteration 300 (after both
CFL ramps) the column run's u and v residuals are at most half of the profile
run's (the Ekman balance is in the first guess, so the inertial spin-up is
gone; the first iteration is not compared: over a hill the lifted spiral's
flank imbalance makes the column seed's first residual larger) and that the
column run reaches the verdict in fewer iterations. Either argument may be a
glob prefix: the newest matching checkpoint is used.
"""
import glob
import os
import sys


def read_history(path):
    with open(os.path.join(path, "Driver"), encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    i = 0
    while i < len(lines) and not lines[i].startswith("records "):
        i += 1
    n = int(lines[i].split()[1])
    out = []
    for ln in lines[i + 1:i + 1 + n]:
        parts = ln.split()
        nres = int(parts[5])
        out.append([float.fromhex(x) for x in parts[6:6 + nres]])
    return out


def resolve(arg):
    if os.path.isdir(arg):
        return arg
    hits = sorted(glob.glob(arg + "*"), key=os.path.getmtime)
    hits = [h for h in hits if os.path.isdir(h)]
    if not hits:
        raise SystemExit("check_column_guess: no checkpoint matches " + arg)
    return hits[-1]


def main(argv):
    if len(argv) != 3:
        raise SystemExit(__doc__)
    prof, col = read_history(resolve(argv[1])), read_history(resolve(argv[2]))
    ok = True
    print("check_column_guess: profile %d iterations, column %d iterations" % (len(prof), len(col)))
    at = min(300, len(prof), len(col)) - 1
    for name, idx in (("u", 0), ("v", 1)):
        ratio = col[at][idx] / prof[at][idx] if prof[at][idx] > 0 else float("inf")
        passed = ratio <= 0.5
        ok = ok and passed
        print("  residual %s at iteration %d: column / profile = %.3e (<= 0.5): %s" % (name, at + 1, ratio, "pass" if passed else "FAIL"))
    passed = len(col) < len(prof)
    ok = ok and passed
    print("  iterations to the verdict: column %d < profile %d: %s" % (len(col), len(prof), "pass" if passed else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
