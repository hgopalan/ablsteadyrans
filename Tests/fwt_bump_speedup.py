#!/usr/bin/env python3
"""Write the FastWindTerrain section of the RT-1.12.2 gold (Tests/Gold/init_bump_trend.json).

    fwt_bump_speedup.py GOLD.json SHA ALPHA=PLOTFILE [ALPHA=PLOTFILE ...]

Each PLOTFILE is a FastWindTerrain plotfile (grid.output_format = plt) of
regtests/phase6_solve_correction/inputs_bump run with poisson.alpha_v = ALPHA
(alpha_h = 1). The hilltop column is the one with the highest terrain_z;
the speed-up at a height h above the ground is |(u, v)| / |(u0, v0)| of the
fluid cells there (u0 is FastWindTerrain's own inflow profile by height
above the ground), interpolated linearly to the gold's heights. SHA is the
FastWindTerrain commit the plotfiles came from (external/VERSIONS). The
other sections of GOLD.json are kept; the "fwt" section is replaced.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def speedup(plotfile, heights):
    _, d = erf_plotfile.read_fields(plotfile, ["z_cc", "terrain_z", "u", "v", "u0", "v0"])
    nx, ny = len(d["u"]), len(d["u"][0])
    zt, i, j = max((d["terrain_z"][i][j][0], i, j) for i in range(nx) for j in range(ny))
    col = []
    for k in range(len(d["u"][i][j])):
        s0 = math.hypot(d["u0"][i][j][k], d["v0"][i][j][k])
        if s0 == 0:   # a solid cell
            continue
        col.append((d["z_cc"][i][j][k] - zt, math.hypot(d["u"][i][j][k], d["v"][i][j][k]) / s0))
    col.sort()
    out = {}
    for hs in heights:
        for (h0, r0), (h1, r1) in zip(col, col[1:]):
            if h0 <= hs <= h1:
                out[str(hs)] = r0 + (r1 - r0) * (hs - h0) / (h1 - h0)
                break
        else:
            raise SystemExit(f"fwt_bump_speedup: height {hs} m is outside the hilltop column of {plotfile}")
    return (i, j, zt), out


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    gold, sha, cases = argv[0], argv[1], argv[2:]
    with open(gold) as fh:
        g = json.load(fh)
    fwt = {"sha": sha, "speedup": {}, "hilltop": {}}
    for case in cases:
        alpha, plotfile = case.split("=", 1)
        top, s = speedup(plotfile, g["heights"])
        fwt["speedup"][alpha] = s
        fwt["hilltop"][alpha] = {"i": top[0], "j": top[1], "terrain_z": top[2]}
        print(alpha, top, {k: round(v, 4) for k, v in s.items()})
    g["fwt"] = fwt
    with open(gold, "w") as fh:
        json.dump(g, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"fwt_bump_speedup: wrote the fwt section of {gold}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
