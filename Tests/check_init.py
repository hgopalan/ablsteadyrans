#!/usr/bin/env python3
"""Run ablsteady with the mass-consistent initializer and check the result (RT-1.12.1, RT-1.12.2).

    check_init.py askervein GOLD.json [--write] -- COMMAND [ARGS ...]
    check_init.py trend GOLD.json [--write] -- COMMAND [ARGS ...]

askervein (RT-1.12.1) runs COMMAND once and checks
  * the Poisson iterations, at most 10 % above the gold count;
  * the divergence after the correction, at most 10 x solver.poisson_reltol
    (as the run prints it) times the divergence before;
  * the divergence before, equal to the gold to a relative 1e-9 (the first
    guess and the mesh are fixed, so it only changes if they do);
  * from the plotfile: max |w| and the largest horizontal speed equal to the
    gold to a relative 1e-6, theta equal to the sounding everywhere and p'
    exactly 0 (theta is the reference itself), k positive everywhere.
The plotfile itself is compared with amrex_fcompare in a separate test.

trend (RT-1.12.2) runs COMMAND once per alpha_v of the gold (appending
init.alpha_v) and compares the hilltop speed-up with
FastWindTerrain's on the same hill, profile and alpha set (the gold's "fwt"
section, written by Tests/fwt_bump_speedup.py from FastWindTerrain's
plotfiles). The speed-up is |u_h| / S_0(h) in the hilltop column, S_0 the
inflow power law at the height h above the ground, interpolated to the
gold's heights. Checks:
  * every speed-up above 1 (the hill speeds the flow up) in both codes;
  * the ranking of the alpha_v values by speed-up equals FastWindTerrain's
    at the gold's ranking heights (not where two values nearly tie);
  * |ablsteady - FastWindTerrain| within the gold's tolerance per height;
  * ablsteady's own values equal to the gold to a relative 1e-6 and its
    iterations at most 10 % above the gold's (regression).
--write stores the measured ablsteady values as the new gold, keeping the
"fwt" section (review the diff before committing it). The exit code is the
verdict.
"""
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"check_init: the run exited {proc.returncode}")
    return proc.stdout


def measure(out):
    m = re.search(r"ablsteady: initializer: (\d+) Poisson iterations \(reltol (\S+)\), divergence (\S+) before, (\S+) after", out)
    w = re.search(r"ablsteady: wrote (\S+)", out)
    if not m or not w:
        raise SystemExit("check_init: the run did not print the initializer report and the plotfile name")
    return {"iterations": int(m.group(1)), "reltol": float(m.group(2)), "div_before": float(m.group(3)),
            "div_after": float(m.group(4)), "plotfile": w.group(1)}


def rows_common(m, g):
    return [
        ("iterations", m["iterations"], f"<= {1.1 * g['iterations']:.1f}", m["iterations"] <= 1.1 * g["iterations"]),
        ("divergence after / before", m["div_after"] / m["div_before"], f"<= {10 * m['reltol']:g}",
         m["div_after"] <= 10 * m["reltol"] * m["div_before"]),
        ("divergence before", m["div_before"], f"{g['div_before']:.9e} (rel 1e-9)",
         abs(m["div_before"] - g["div_before"]) <= 1e-9 * abs(g["div_before"])),
    ]


def report(rows):
    ok = True
    print(f"{'check':34s} {'measured':>18} {'target':>34} pass")
    for name, meas, target, passed in rows:
        print(f"{name:34s} {meas!s:>18} {target:>34} {'pass' if passed else 'FAIL'}")
        ok &= passed
    return ok


# --- RT-1.12.1 -------------------------------------------------------------
def askervein_fields(plotfile):
    _, d = erf_plotfile.read_fields(plotfile, ["u", "v", "w", "theta", "p", "tke"])
    wmax = speed = 0.0
    theta_lo, theta_hi, p_abs, k_min = math.inf, -math.inf, 0.0, math.inf
    for ui, vi, wi, ti, pi, ki in zip(d["u"], d["v"], d["w"], d["theta"], d["p"], d["tke"]):
        for uj, vj, wj, tj, pj, kj in zip(ui, vi, wi, ti, pi, ki):
            for u, v, w, t, p, k in zip(uj, vj, wj, tj, pj, kj):
                wmax = max(wmax, abs(w))
                speed = max(speed, math.hypot(u, v))
                theta_lo, theta_hi = min(theta_lo, t), max(theta_hi, t)
                p_abs = max(p_abs, abs(p))
                k_min = min(k_min, k)
    return {"wmax": wmax, "speed_max": speed, "theta_min": theta_lo, "theta_max": theta_hi, "p_absmax": p_abs,
            "k_min": k_min}


def askervein(gold, write, cmd):
    m = measure(run(cmd))
    m.update(askervein_fields(m["plotfile"]))
    del m["plotfile"]
    if write:
        with open(gold, "w") as fh:
            json.dump(m, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"check_init: wrote {gold}")
        return 0
    with open(gold) as fh:
        g = json.load(fh)
    rows = rows_common(m, g) + [
        ("max |w|", m["wmax"], f"{g['wmax']:.9e} (rel 1e-6)", abs(m["wmax"] - g["wmax"]) <= 1e-6 * g["wmax"]),
        ("max horizontal speed", m["speed_max"], f"{g['speed_max']:.9e} (rel 1e-6)",
         abs(m["speed_max"] - g["speed_max"]) <= 1e-6 * g["speed_max"]),
        ("theta", (m["theta_min"], m["theta_max"]), f"= {g['theta_min']} everywhere",
         m["theta_min"] == g["theta_min"] and m["theta_max"] == g["theta_min"]),
        ("max |p'|", m["p_absmax"], "= 0", m["p_absmax"] == 0.0),
        ("min k", m["k_min"], "> 0", m["k_min"] > 0),
    ]
    return 0 if report(rows) else 1


# --- RT-1.12.2 -------------------------------------------------------------
def speedup_column(plotfile, mesh, heights):
    """Speed-up |u_h| / S_0(h) at the gold's heights above the ground in the hilltop column."""
    _, d = erf_plotfile.read_fields(plotfile, ["u", "v", "z"])
    i, j = mesh["column"]
    dz0, r, nz = mesh["dz0"], mesh["stretch_ratio"], mesh["n_layers"]
    ztop = dz0 * (r ** nz - 1) / (r - 1)   # the physical top: where the layers end

    def zflat(k):   # the flat mesh's cell-centre heights
        return dz0 * (r ** k - 1) / (r - 1) + 0.5 * dz0 * r ** k

    # The basic terrain-following map: z = zg + z_flat (1 - zg / ztop), so the
    # ground height of the column follows from the first cell centre.
    zg = (d["z"][i][j][0] - zflat(0)) / (1 - zflat(0) / ztop)
    col = []
    for k in range(len(d["z"][i][j])):
        h = d["z"][i][j][k] - zg
        s0 = mesh["speed_ref"] * (h / mesh["z_ref"]) ** mesh["exponent"]
        col.append((h, math.hypot(d["u"][i][j][k], d["v"][i][j][k]) / s0))
    out = {}
    for hs in heights:
        for (h0, r0), (h1, r1) in zip(col, col[1:]):
            if h0 <= hs <= h1:
                out[str(hs)] = r0 + (r1 - r0) * (hs - h0) / (h1 - h0)
                break
        else:
            raise SystemExit(f"check_init: height {hs} m is outside the hilltop column")
    return zg, out


def ranking(table, alphas, h):
    return sorted(alphas, key=lambda a: table[a][h])


def trend(gold, write, cmd):
    with open(gold) as fh:
        g = json.load(fh)
    fwt, mesh, heights, alphas = g["fwt"], g["mesh"], g["heights"], g["alphas"]
    abl = {}
    for a in alphas:
        m = measure(run(cmd + [f"init.alpha_v={a}"]))   # each run overwrites the deck's plotfile; read it at once
        zg, s = speedup_column(m["plotfile"], mesh, heights)
        abl[a] = {"iterations": m["iterations"], "reltol": m["reltol"], "div_before": m["div_before"],
                  "div_after": m["div_after"], "ground": zg, "speedup": s}
    if write:
        g["ablsteady"] = abl
        with open(gold, "w") as fh:
            json.dump(g, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"check_init: wrote the ablsteady section of {gold}")
        return 0
    rows = []
    print(f"{'alpha_v':>8} " + " ".join(f"{'abl ' + str(h) + ' m':>12} {'fwt ' + str(h) + ' m':>12}" for h in heights))
    for a in alphas:
        print(f"{a:>8} " + " ".join(f"{abl[a]['speedup'][str(h)]:12.4f} {fwt['speedup'][a][str(h)]:12.4f}" for h in heights))
    for a in alphas:
        m, ga = abl[a], g["ablsteady"][a]
        rows += rows_common(m, ga)
        for h in heights:
            hs = str(h)
            rows.append((f"alpha_v {a}: speed-up at {h} m", m["speedup"][hs], f"{ga['speedup'][hs]:.9e} (rel 1e-6)",
                         abs(m["speedup"][hs] - ga["speedup"][hs]) <= 1e-6 * ga["speedup"][hs]))
            rows.append((f"alpha_v {a}: vs FastWindTerrain at {h} m", m["speedup"][hs] - fwt["speedup"][a][hs],
                         f"|.| <= {g['tolerance'][hs]}", abs(m["speedup"][hs] - fwt["speedup"][a][hs]) <= g["tolerance"][hs]))
    for h in heights:
        hs = str(h)
        if h in g.get("ranking_heights", heights):
            ra, rf = ranking({a: abl[a]["speedup"] for a in alphas}, alphas, hs), ranking(fwt["speedup"], alphas, hs)
            rows.append((f"ranking of alpha_v at {h} m", " < ".join(ra), " < ".join(rf), ra == rf))
        rows.append((f"speed-up > 1 at {h} m (both codes)", min(min(abl[a]["speedup"][hs] for a in alphas),
                                                                   min(fwt["speedup"][a][hs] for a in alphas)), "> 1",
                     all(abl[a]["speedup"][hs] > 1 and fwt["speedup"][a][hs] > 1 for a in alphas)))
    return 0 if report(rows) else 1


def main(argv):
    if len(argv) < 3 or "--" not in argv or argv[0] not in ("askervein", "trend"):
        print(__doc__)
        return 2
    sep = argv.index("--")
    mode, gold = argv[0], argv[1]
    write = "--write" in argv[2:sep]
    cmd = argv[sep + 1:]
    if not cmd:
        print(__doc__)
        return 2
    return askervein(gold, write, cmd) if mode == "askervein" else trend(gold, write, cmd)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
