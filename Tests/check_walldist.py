#!/usr/bin/env python3
"""Run wall_distance and apply ERF's wall-distance checks (RT-1.11.2).

    check_walldist.py --hmax H --half-width L -- COMMAND [ARGS ...]

Runs COMMAND (Tests/Tools/wall_distance.cpp on a 2-D Witch of Agnesi ridge,
x-centred in the domain) and checks the plotfile it wrote against the exact
distance to the ridge, with the tolerances of ERF's
Exec/CanonicalTests/Canonical_RANS/Neutral_Hill_2D/check_hill.py (ERF 92410d57):
  * max relative error <= 0.15, mean relative error <= 0.03,
  * max absolute error where the exact distance is below 100 m <= 0.2 vertical cells,
  * every value finite.
The exact distance is ERF's scan over the foot point with a golden-section
refinement, ported unchanged. Prints one row per check; the exit code is the
verdict. The field itself is compared with amrex_fcompare in a separate test.
"""
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)


def ridge(x, hmax, L):
    return hmax / (1.0 + (x / L) ** 2)


def exact_distance(x, z, hmax, L):
    """ERF's check_hill.py exact_distance."""
    best_s, best_d = x, abs(z - ridge(x, hmax, L))
    span = 4.0 * L
    n = 200
    for m in range(-n, n + 1):
        s = x + span * m / n
        d = math.hypot(x - s, z - ridge(s, hmax, L))
        if d < best_d:
            best_s, best_d = s, d
    a, b = best_s - span / n, best_s + span / n
    gr = 0.5 * (math.sqrt(5.0) - 1.0)
    c, d_ = b - gr * (b - a), a + gr * (b - a)
    for _ in range(60):
        fc = math.hypot(x - c, z - ridge(c, hmax, L))
        fd = math.hypot(x - d_, z - ridge(d_, hmax, L))
        if fc < fd:
            b, d_ = d_, c
            c = b - gr * (b - a)
        else:
            a, c = c, d_
            d_ = a + gr * (b - a)
    s = 0.5 * (a + b)
    return math.hypot(x - s, z - ridge(s, hmax, L))


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    try:
        hmax = float(opts[opts.index("--hmax") + 1])
        L = float(opts[opts.index("--half-width") + 1])
    except (ValueError, IndexError):
        print(__doc__)
        return 2
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    print(proc.stdout)
    if proc.returncode != 0:
        print(f"check_walldist: the run exited {proc.returncode}")
        return 1
    plt = re.findall(r"wall_distance: wrote (\S+)", proc.stdout)
    if not plt:
        print("check_walldist: the run wrote no plotfile")
        return 1
    hdr, f = erf_plotfile.read_fields(plt[0], ["walldist", "z"])
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    dx, dz = hdr["dx"][0], hdr["dx"][2]
    x_lo = hdr["prob_lo"][0]
    xcen = 0.5 * (hdr["prob_lo"][0] + hdr["prob_hi"][0])
    max_rel = sum_rel = max_abs_near = 0.0
    n = 0
    finite = True
    for i in range(nx):
        x = x_lo + (i + 0.5) * dx - xcen
        for j in range(ny):
            for k in range(nz):
                d = f["walldist"][i][j][k]
                finite &= math.isfinite(d)
                de = exact_distance(x, f["z"][i][j][k], hmax, L)
                err = d - de
                if de < 100.0:
                    max_abs_near = max(max_abs_near, abs(err))
                rel = abs(err) / de
                max_rel = max(max_rel, rel)
                sum_rel += rel
                n += 1
    rows = [
        ("all values finite", finite, True, finite),
        ("max rel err walldist vs exact ridge distance", max_rel, "<= 0.15", max_rel <= 0.15),
        ("mean rel err walldist vs exact ridge distance", sum_rel / n, "<= 0.03", sum_rel / n <= 0.03),
        ("max abs err walldist, d < 100 m [cells]", max_abs_near / dz, "<= 0.2", max_abs_near / dz <= 0.2),
    ]
    ok = True
    print(f"{'check':48s} {'measured':>22} {'target':>10} pass")
    for name, meas, target, passed in rows:
        print(f"{name:48s} {meas!s:>22} {target!s:>10} {'pass' if passed else 'FAIL'}")
        ok &= bool(passed)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
