#!/usr/bin/env python3
"""Inflow table (bc.inflow_profile = file) from an ERF Dirichlet inflow file.

    erf_inflow_table.py ERF_FILE --theta THETA [--out TABLE]
    erf_inflow_table.py ERF_FILE --theta THETA --check TABLE

ERF_FILE holds rows of "z u v w" (ERF's xlo.dirichlet_file, for example
extern_erf/inflow_powerlaw.txt). The table has rows of "z u v theta": a z = 0
row first (u = v = 0), then ERF's rows with w dropped and theta = THETA (a
neutral case). k is left out, so ablsteady sets the equilibrium k of its wall
law at the inflow (docs/boundary_conditions.rst).

Prints the wind direction (degrees, meteorological: where it blows from) and
speed of the top row. With --check TABLE (regression test of the committed
table) the exit code is 0 only if TABLE's rows equal the rows this tool would
write, to 1e-9.
"""
import argparse
import math
import sys


def read_erf(path):
    rows = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            s = line.split("#", 1)[0].split()
            if not s:
                continue
            if len(s) != 4:
                sys.exit(f"{path}:{n}: expected 4 columns (z u v w), got {len(s)}")
            rows.append([float(x) for x in s])
    if len(rows) < 2:
        sys.exit(f"{path}: needs at least 2 rows")
    for a, b in zip(rows, rows[1:]):
        if not b[0] > a[0]:
            sys.exit(f"{path}: z must increase (z = {a[0]} then {b[0]})")
    if rows[0][0] <= 0.0:
        sys.exit(f"{path}: the first z must be > 0 (the table adds the z = 0 row)")
    return rows


def table(rows, theta):
    return [[0.0, 0.0, 0.0, theta]] + [[z, u, v, theta] for z, u, v, _ in rows]


def read_table(path):
    out = []
    with open(path) as f:
        for line in f:
            s = line.split("#", 1)[0].split()
            if s:
                out.append([float(x) for x in s])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("erf_file")
    ap.add_argument("--theta", type=float, required=True, help="potential temperature of every row (K, > 0)")
    ap.add_argument("--out")
    ap.add_argument("--check")
    a = ap.parse_args()
    if not (a.theta > 0.0 and math.isfinite(a.theta)):
        sys.exit("--theta must be finite and > 0")
    t = table(read_erf(a.erf_file), a.theta)
    z, u, v, _ = t[-1]
    wdir = math.degrees(math.atan2(-u, -v)) % 360.0
    print(f"top row z = {z:g} m: speed {math.hypot(u, v):.4f} m/s from {wdir:.2f} deg")
    if a.check:
        got = read_table(a.check)
        if len(got) != len(t):
            print(f"FAIL: {a.check} has {len(got)} rows, expected {len(t)}")
            return 1
        for i, (g, e) in enumerate(zip(got, t)):
            if len(g) != 4 or any(abs(x - y) > 1e-9 * max(1.0, abs(y)) for x, y in zip(g, e)):
                print(f"FAIL: {a.check} row {i + 1}: {g} != {e}")
                return 1
        print(f"PASS: {a.check} matches {a.erf_file}")
        return 0
    if a.out:
        with open(a.out, "w") as f:
            f.write(f"# z u v theta: {a.erf_file.split('/')[-1]} (ERF, z u v w) with w dropped, "
                    f"theta = {a.theta:g} K, a z = 0 row added (tools/erf_inflow_table.py)\n")
            for r in t:
                f.write(f"{r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {r[3]:.6f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
