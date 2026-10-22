#!/usr/bin/env python3
"""The 2-D ridge of Ishihara & Zhou (2025, Data in Brief 63, 112260; CC BY 4.0) at full scale (tunnel x1000):
z = h cos^2(pi (x - xc) / 2 L) for |x - xc| <= L, h = 40 m, L = 100 m (maximum slope about 32 degrees),
xc = 800 m in a 2400 m domain under an 896 m slip ceiling. Smooth (z0 0.01 m, u* 0.21 m/s; mmc1, mmc3) and
rough (z0 0.3 m, u* 0.32 m/s; mmc2, mmc4) surfaces. The setup follows H. Gopalan's Kynema-SGF scripts
(kynema/kynema-sgf#2031, #2036); the LES values are those of Liu, Diao & Ishihara (2019, Renewable Energy 136,
968) digitized there.

    ishihara_zhou.py prepare DATA OUTDIR
        From the Data-in-Brief files in DATA (mmc1-mmc4.xls): the terrain (ridge.xyz), the inflow tables
        (inflow_smooth.txt, inflow_rough.txt: z u v theta k, k = (su^2 + sv^2 + sw^2)/2; the log law with the
        paper's u* and z0 below the lowest point, the free stream held above the highest) and the stations
        (stations_smooth.csv, stations_rough.csv: x, z, U, W, k, all at full scale, x from the crest).
    ishihara_zhou.py compare PLOTFILE --surface smooth|rough --outdir OUTDIR [--label NAME] [--out CSV]
        The run's u along the near-surface line (5 m smooth, 7 m rough, above the local ground), the
        reversed-flow extent there (separation and reattachment, m from the crest), and U at every station
        point with the rms against the tunnel; the tunnel and LES reattachment for reference.
"""
import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))   # askervein_lines.py, tiff_to_xyz.py
import askervein_lines as al  # noqa: E402

H, L, XC, LX, LZ = 40.0, 100.0, 800.0, 2400.0, 896.0
KAPPA, THETA = 0.41, 300.0
SURF = {"smooth": dict(approach="mmc1.xls", ridge="mmc3.xls", z0=0.01, ustar=0.21, zline=5.0),
        "rough": dict(approach="mmc2.xls", ridge="mmc4.xls", z0=0.3, ustar=0.32, zline=7.0)}
# Reversed flow at the near-surface line (m from the crest): tunnel ranges from the Data-in-Brief stations,
# LES from Liu, Diao & Ishihara (2019) Fig. 11 (as digitized in kynema-sgf#2036).
REF = {"smooth": dict(tunnel_reattach=(150.0, 200.0), les=(45.0, 180.0)),
       "rough": dict(tunnel_reattach=(250.0, 300.0), les=(37.0, 239.0))}


def ground(x):
    s = abs(x - XC)
    return H * math.cos(0.5 * math.pi * s / L) ** 2 if s <= L else 0.0


def read_xls(path):
    import xlrd
    sh = xlrd.open_workbook(path).sheet_by_index(0)
    return [[sh.cell_value(r, c) for c in range(sh.ncols)] for r in range(sh.nrows) if isinstance(sh.cell_value(r, 0), float)]


def prepare(a):
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "ridge.xyz"), "w") as f:
        f.write("# Ishihara-Zhou 2-D ridge at full scale: z = 40 cos^2(pi (x - 800) / 200) for |x - 800| <= 100 m\n")
        for i in range(int(LX / 2) + 1):
            x = 2.0 * i
            for y in (-16.0, 0.0, 16.0, 32.0, 48.0):
                f.write(f"{x:.1f} {y:.1f} {ground(x):.5f}\n")
    for name, s in SURF.items():
        rows = sorted(read_xls(os.path.join(a.data, s["approach"])))
        prof = [(r[0], r[1], 0.5 * (r[2] ** 2 + r[3] ** 2 + r[4] ** 2)) for r in rows]   # mm -> m: x1000 / 1000
        zmin, umin, kmin = prof[0]
        below = [(z, min(max(s["ustar"] / KAPPA * math.log(z / s["z0"]), 0.0), umin) if z > s["z0"] else 0.0, kmin)
                 for z in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0) if z < zmin]
        top = prof[-1]
        table = [(0.0, 0.0, kmin)] + below + prof + [(LZ + 10.0, top[1], top[2])]
        with open(os.path.join(a.outdir, f"inflow_{name}.txt"), "w") as f:
            f.write(f"# z u v theta k: Ishihara & Zhou (2025) {s['approach']} at full scale (CC BY 4.0); log law with "
                    f"u* {s['ustar']} m/s, z0 {s['z0']} m below {zmin} m\n")
            f.writelines(f"{z:.3f} {u:.5f} 0.0 {THETA:.1f} {k:.6f}\n" for z, u, k in table)
        with open(os.path.join(a.outdir, f"stations_{name}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x", "z", "U", "W", "k"])
            for r in read_xls(os.path.join(a.data, s["ridge"])):
                w.writerow([r[0], r[1], r[2], r[3], 0.5 * (r[4] ** 2 + r[5] ** 2 + r[6] ** 2)])
    print(f"{a.outdir}: ridge.xyz, inflow_smooth.txt, inflow_rough.txt, stations_smooth.csv, stations_rough.csv")


def column_u(cols, x):
    """u and the absolute heights of the cell-centre column nearest y = middle at x (linear in x)."""
    h = cols.hdr
    lo, dx = h["prob_lo"], h["dx"]
    ny = h["hi"][1] - h["lo"][1] + 1
    fx = (x - lo[0]) / dx[0] - 0.5
    i0 = int(math.floor(fx))
    t = fx - i0
    u0, _, z0 = cols.column(i0, ny // 2)
    u1, _, z1 = cols.column(i0 + 1, ny // 2)
    return [(1 - t) * a + t * b for a, b in zip(u0, u1)], [(1 - t) * a + t * b for a, b in zip(z0, z1)]


def at_height(us, zs, z):
    if z <= zs[0]:
        return us[0]
    for k in range(len(zs) - 1):
        if zs[k] <= z <= zs[k + 1]:
            t = (z - zs[k]) / (zs[k + 1] - zs[k])
            return (1 - t) * us[k] + t * us[k + 1]
    return us[-1]


def compare(a):
    s = SURF[a.surface]
    cols = al.Columns(a.plotfile)
    line = []
    for i in range(0, 451):
        xr = -100.0 + i          # from 100 m before the crest to 350 m behind it
        x = XC + xr
        us, zs = column_u(cols, x)
        g = ground(x)
        # Heights above the local ground; below the first centre the wall law from the first centre.
        zrel = [z - g for z in zs]
        zt = s["zline"]
        if zt < zrel[0]:
            u = us[0] * math.log((zt + s["z0"]) / s["z0"]) / math.log((zrel[0] + s["z0"]) / s["z0"])
        else:
            u = at_height(us, zrel, zt)
        line.append((xr, u))
    rev = [x for x, u in line if u < 0]
    sep, rea = (min(rev), max(rev)) if rev else (None, None)
    st = list(csv.DictReader(open(os.path.join(a.outdir, f"stations_{a.surface}.csv"))))
    err, out = [], []
    for r in st:
        xr, z, um = float(r["x"]), float(r["z"]), float(r["U"])
        us, zs = column_u(cols, XC + xr)
        um_model = at_height(us, zs, z)
        err.append(um_model - um)
        out.append((xr, z, um, um_model))
    rms = math.sqrt(sum(e * e for e in err) / len(err))
    ref = REF[a.surface]
    print(f"{a.label or a.plotfile} ({a.surface}): reversed flow at {s['zline']:g} m: "
          + (f"{sep:.0f} to {rea:.0f} m" if rev else "none")
          + f"; tunnel reattachment {ref['tunnel_reattach'][0]:.0f}-{ref['tunnel_reattach'][1]:.0f} m, LES {ref['les'][0]:.0f}"
          f" to {ref['les'][1]:.0f} m; U rms vs tunnel {rms:.3f} m/s over {len(err)} points")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["kind", "x", "z", "U_tunnel", "U_model"])
            for xr, u in line:
                w.writerow(["line", xr, s["zline"], "", f"{u:.4f}"])
            for xr, z, um, mm in out:
                w.writerow(["station", xr, z, um, f"{mm:.4f}"])


GUST, TUNNEL, MWT, LES = "#2ca02c", "black", "#1f77b4", "grey"


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s = SURF[a.surface]
    cols = al.Columns(a.plotfile)
    xs, us = [], []
    for i in range(0, 551):
        xr = -150.0 + i
        u, zs = column_u(cols, XC + xr)
        zrel = [z - ground(XC + xr) for z in zs]
        zt = s["zline"]
        us.append(u[0] * math.log((zt + s["z0"]) / s["z0"]) / math.log((zrel[0] + s["z0"]) / s["z0"])
                  if zt < zrel[0] else at_height(u, zrel, zt))
        xs.append(xr)
    st = list(csv.DictReader(open(os.path.join(a.outdir, f"stations_{a.surface}.csv"))))
    xst = sorted({float(r["x"]) for r in st})
    fig = plt.figure(figsize=(13, 5.5))
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(xs, us, color=GUST, label="GUST - CFD")
    line = [(float(r["x"]), float(r["U"])) for r in st if abs(float(r["z"]) - s["zline"]) < 1e-6]
    ax.plot([p[0] for p in line], [p[1] for p in line], "o", color=TUNNEL, label=f"tunnel (z = {s['zline']:g} m)")
    ref = REF[a.surface]
    ax.axvspan(*ref["tunnel_reattach"], color=TUNNEL, alpha=0.12, label="tunnel reattachment")
    ax.axvline(ref["les"][1], color=LES, ls=":", label="LES reattachment")
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_xlabel("x from the crest (m)")
    ax.set_ylabel(f"u at {s['zline']:g} m above the ground (m/s)")
    ax.set_title(f"Ishihara-Zhou ridge, {a.surface}: near-surface line")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax2 = fig.add_subplot(1, 3, (2, 3))
    shift = 0.0
    for x in xst:
        pts = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
        u, zs = column_u(cols, XC + x)
        g = ground(XC + x)
        zz = [z for z in zs if z - g <= 200]
        uu = u[:len(zz)]
        sc = 1 / 10.0                      # 10 m/s per station spacing
        ax2.plot([x + uu[k] * sc * 50 for k in range(len(zz))], zz, color=GUST, lw=1.2,
                 label="GUST - CFD" if x == xst[0] else None)
        ax2.plot([x + p[1] * sc * 50 for p in pts], [p[0] for p in pts], "o", ms=3, color=TUNNEL,   # z: absolute
                 label="tunnel" if x == xst[0] else None)
        ax2.axvline(x, color="0.85", lw=0.6)
    gx = [XC - 200 + 2 * i for i in range(301)]
    ax2.fill_between([x - XC for x in gx], [ground(x) for x in gx], 0, color="0.75")
    ax2.set_xlabel("x from the crest (m); profiles: U, 50 m per 10 m/s")
    ax2.set_ylabel("z (m)")
    ax2.set_ylim(0, 200)
    ax2.set_title("U profiles at the tunnel stations")
    ax2.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("data")
    p.add_argument("outdir")
    c = sub.add_parser("compare")
    c.add_argument("plotfile")
    c.add_argument("--surface", choices=SURF, required=True)
    c.add_argument("--outdir", required=True)
    c.add_argument("--label")
    c.add_argument("--out")
    q = sub.add_parser("plot")
    q.add_argument("plotfile")
    q.add_argument("--surface", choices=SURF, required=True)
    q.add_argument("--outdir", required=True)
    q.add_argument("--out", required=True)
    a = ap.parse_args()
    {"prepare": prepare, "compare": compare, "plot": plot}[a.cmd](a)


if __name__ == "__main__":
    main()
