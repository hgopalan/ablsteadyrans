#!/usr/bin/env python3
"""The RUSHIL two-dimensional hills (Khurshudyan, Snyder & Nekrasov 1981, EPA-600/4-81-067; ERCOFTAC Classic
Collection case 69, Application Challenge AC5-05) at the full scale of Manning, Hancock & Whiting (2010, Wind
Engineering 34, 477): tunnel x1000, H = 117 m, a/H = 3, 5, 8 (Hills 3, 5, 8), z0 = 0.157 m, the flat-floor
profile as inflow. The hill shape (Khurshudyan et al. 1981), for |xi| <= a:

    x = xi/2 [1 + a^2 / (xi^2 + m^2 (a^2 - xi^2))],
    z = m/2 sqrt(a^2 - xi^2) [1 - a^2 / (xi^2 + m^2 (a^2 - xi^2))],   H = a (m - 1/m) / 2,

which gives maximum slopes 0.492, 0.295 and 0.185 (Manning et al.'s Table 2.1: 0.49, 0.29).

    rushil.py prepare DATA OUTDIR
        From the ERCOFTAC AC5-05 files in DATA: hill3.xyz, hill5.xyz, hill8.xyz (the terrain, with the summit at
        x = XC), inflow.txt (z u v theta k from AC5-05_Flatrefcase.dat; k = 0.75 (u'^2 + w'^2), v'^2 taken as the
        mean of u'^2 and w'^2) and stations_h3.csv, stations_h5.csv, stations_h8.csv (x from the summit, z above
        the surface, U, all at full scale).
    rushil.py compare PLOTFILE --hill 3|5|8 --outdir OUTDIR [--label NAME] [--out CSV]
        The reversed-flow extent at Z = 0.1 H above the ground (in units of a, from the summit), the speed-up
        U / U_in(0.8 H) on the planes Z = 0.1 H ... 0.6 H (Manning et al.'s Figs. 4.3, 4.4), and U at every
        station point with its rms against the tunnel.
    rushil.py plot PLOTFILE --hill 3|5|8 --outdir OUTDIR [--manning CSV] --out PNG
        Top: speed-up on the planes Z = 0.1 H ... 0.6 H. Bottom: U profiles at the tunnel stations. With --manning
        (a digitized file of Manning et al.'s Figs. 4.3-4.6, Hills 3 and 5) the Meteodyn WT curves and their tunnel
        points are added.
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))   # askervein_lines.py, tiff_to_xyz.py
import askervein_lines as al  # noqa: E402

SCALE, H, XC = 1000.0, 117.0, 4680.0   # summit at 40 H from the inflow
THETA = 300.0
HILLS = {3: 3.0, 5: 5.0, 8: 8.0}
# Reversed flow at Z = 0.1 H in units of a from the summit: tunnel (Manning et al.'s Table 2.2, from the
# RUSHIL reports) and Meteodyn WT (their Fig. 4.3 at Z = 0.1 H).
REF = {3: dict(tunnel=(0.27, 2.17), mwt=(0.8, 1.8)), 5: dict(tunnel=None, mwt=None), 8: dict(tunnel=None, mwt=None)}


def shape(a):
    """(x, z) points of the hill of half-length a (m), from xi = -a to a."""
    q = 2 * H / a
    m = (q + math.sqrt(q * q + 4)) / 2
    pts = []
    n = 4000
    for i in range(-n, n + 1):
        xi = a * i / n
        d = xi * xi + m * m * (a * a - xi * xi)
        pts.append((0.5 * xi * (1 + a * a / d), 0.5 * m * math.sqrt(max(a * a - xi * xi, 0.0)) * (1 - a * a / d)))
    return pts


def ground_fn(a):
    pts = shape(a)
    xs = [p[0] for p in pts]

    def g(x):
        r = x - XC
        if r <= xs[0] or r >= xs[-1]:
            return 0.0
        lo, hi = 0, len(xs) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= r:
                lo = mid
            else:
                hi = mid
        t = (r - xs[lo]) / (xs[hi] - xs[lo]) if xs[hi] > xs[lo] else 0.0
        return pts[lo][1] + t * (pts[hi][1] - pts[lo][1])
    return g


def read_dat(path):
    rows, x = [], None
    for line in open(path, errors="replace"):
        m = re.search(r"x\s*=\s*([-+]?\d*\.?\d+)\s*a?\s*=\s*([-+]?\d*\.?\d+)\s*m", line)
        if m:
            x = float(m.group(2))
        p = line.split()
        if p and not line.lstrip().startswith(("#", "%")):
            try:
                rows.append([float(v) for v in p])
            except ValueError:
                pass
    return x, rows


def prepare(a):
    os.makedirs(a.outdir, exist_ok=True)
    for n, ar in HILLS.items():
        g = ground_fn(ar * H)
        with open(os.path.join(a.outdir, f"hill{n}.xyz"), "w") as f:
            f.write(f"# RUSHIL Hill {n} at full scale (H = 117 m, a = {ar * H:.0f} m), summit at x = {XC:.0f} m\n")
            for i in range(0, int(2 * XC / 5) + 1):
                x = 5.0 * i
                for y in (-40.0, 0.0, 40.0, 80.0, 120.0):
                    f.write(f"{x:.1f} {y:.1f} {g(x):.4f}\n")
        with open(os.path.join(a.outdir, f"stations_h{n}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x", "z", "U"])
            for p in sorted(glob.glob(os.path.join(a.data, f"AC5-05_h{n}[0-9][0-9]s.dat"))):
                x, rows = read_dat(p)
                if x is None:
                    continue
                for r in rows:
                    w.writerow([x * SCALE, r[1] * SCALE, r[2]])
    _, rows = read_dat(os.path.join(a.data, "AC5-05_Flatrefcase.dat"))
    prof = sorted((r[1] * SCALE, r[2], 0.75 * (r[4] ** 2 + r[5] ** 2)) for r in rows if len(r) >= 6)
    z0, ustar, kappa = 0.157, 0.178, 0.41
    zmin, umin, kmin = prof[0]
    below = [(z, min(ustar / kappa * math.log(z / z0), umin), kmin) for z in (0.5, 1.0, 2.0, 3.0, 5.0) if z0 < z < zmin]
    top = prof[-1]
    with open(os.path.join(a.outdir, "inflow.txt"), "w") as f:
        f.write("# z u v theta k: RUSHIL flat-floor reference profile (ERCOFTAC AC5-05_Flatrefcase.dat) at full scale; "
                "k = 0.75 (u'^2 + w'^2); log law with u* 0.178 m/s, z0 0.157 m below the lowest point\n")
        table = [(0.0, 0.0, kmin)] + below + prof + [(1700.0, top[1], top[2])]
        f.writelines(f"{z:.3f} {u:.5f} 0.0 {THETA:.1f} {k:.6f}\n" for z, u, k in table)
    print(f"{a.outdir}: hill3/5/8.xyz, stations_h3/5/8.csv ({len(prof)}-point inflow profile)")


def column(cols, x):
    h = cols.hdr
    lo, dx = h["prob_lo"], h["dx"]
    ny = h["hi"][1] - h["lo"][1] + 1
    fx = (x - lo[0]) / dx[0] - 0.5
    i0 = int(math.floor(fx))
    t = fx - i0
    u0, _, z0 = cols.column(i0, ny // 2)
    u1, _, z1 = cols.column(i0 + 1, ny // 2)
    return [(1 - t) * p + t * q for p, q in zip(u0, u1)], [(1 - t) * p + t * q for p, q in zip(z0, z1)]


def u_at(us, zrel, zt, z0=0.157):
    if zt < zrel[0]:
        return us[0] * math.log((zt + z0) / z0) / math.log((zrel[0] + z0) / z0)
    for k in range(len(zrel) - 1):
        if zrel[k] <= zt <= zrel[k + 1]:
            t = (zt - zrel[k]) / (zrel[k + 1] - zrel[k])
            return (1 - t) * us[k] + t * us[k + 1]
    return us[-1]


def compare(a):
    hl = HILLS[a.hill] * H
    g = ground_fn(hl)
    cols = al.Columns(a.plotfile)
    us, zs = column(cols, 200.0)                       # the inflow column, flat ground
    uref = u_at(us, zs, 0.8 * H)
    h = cols.hdr
    xlo = h["prob_lo"][0] + 0.5 * h["dx"][0]
    xhi = h["prob_lo"][0] + (h["hi"][0] - h["lo"][0] + 0.5) * h["dx"][0]
    planes = {}
    for zf in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        line = []
        for i in range(-400, 601):
            xa = i / 100.0
            x = XC + xa * hl
            if not xlo <= x < xhi:
                continue
            us, zs = column(cols, x)
            gr = g(x)
            line.append((xa, u_at(us, [z - gr for z in zs], zf * H) / uref))
        planes[zf] = line
    rev = [xa for xa, s in planes[0.1] if s < 0 and xa > 0]
    sep, rea = (min(rev), max(rev)) if rev else (None, None)
    st = list(csv.DictReader(open(os.path.join(a.outdir, f"stations_h{a.hill}.csv"))))
    err, rows, skipped = [], [], 0
    for r in st:
        x, zr, um = float(r["x"]), float(r["z"]), float(r["U"])
        if not xlo <= XC + x < xhi:     # stations beyond the domain (Hill 3's far wake reaches 15 a)
            skipped += 1
            continue
        us, zs = column(cols, XC + x)
        gr = g(XC + x)
        mm = u_at(us, [z - gr for z in zs], zr)
        err.append(mm - um)
        rows.append((x / hl, zr / H, um, mm))
    rms = math.sqrt(sum(e * e for e in err) / len(err))
    ref = REF[a.hill]
    msg = f"{a.label or a.plotfile} (Hill {a.hill}): reversed flow at 0.1 H: " + (f"{sep:.2f}a to {rea:.2f}a" if rev else "none")
    if ref["tunnel"]:
        msg += f" (tunnel {ref['tunnel'][0]}a-{ref['tunnel'][1]}a, Meteodyn WT {ref['mwt'][0]}a-{ref['mwt'][1]}a)"
    print(msg + f"; U rms vs tunnel {rms:.3f} m/s over {len(err)} points ({skipped} beyond the domain); U_ref(0.8H) {uref:.3f} m/s")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["kind", "zf_or_xa", "x_or_z", "tunnel", "model"])
            for zf, line in planes.items():
                for xa, s in line:
                    w.writerow(["plane", zf, xa, "", f"{s:.4f}"])
            for xa, zh, um, mm in rows:
                w.writerow(["station", f"{xa:.3f}", f"{zh:.4f}", um, f"{mm:.4f}"])


GUST, TUNNEL, MWT, LES = "#2ca02c", "black", "#1f77b4", "grey"


def manning_series(path, fig, panel):
    out = {}
    for r in csv.DictReader(l for l in open(path) if l[0] != "#"):
        if r["figure"] == fig and r["panel"] == panel:
            out.setdefault(r["series"], []).append((float(r["x"]), float(r["y"])))
    return out


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    hl = HILLS[a.hill] * H
    g = ground_fn(hl)
    cols = al.Columns(a.plotfile)
    us, zs = column(cols, 200.0)
    uref = u_at(us, zs, 0.8 * H)
    figs = {3: ("fig4_3", "fig4_5"), 5: ("fig4_4", "fig4_6")}.get(a.hill)
    fig, axs = plt.subplots(3, 4, figsize=(15, 10))
    planes = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    st = list(csv.DictReader(open(os.path.join(a.outdir, f"stations_h{a.hill}.csv"))))
    for n, zf in enumerate(planes):
        ax = axs.flat[n]
        xa, sp = [], []
        for i in range(-500, 601):
            x = XC + i / 100.0 * hl
            if not 50 < x < 9300:
                continue
            u, z = column(cols, x)
            xa.append(i / 100.0)
            sp.append(u_at(u, [zz - g(x) for zz in z], zf * H) / uref)
        ax.plot(xa, sp, color=GUST, label="GUST - CFD")
        if a.manning and figs:
            ser = manning_series(a.manning, figs[0], f"Z{zf:g}h")
            if "MWT" in ser:
                p = sorted(ser["MWT"])
                ax.plot([q[0] for q in p], [q[1] for q in p], "--", color=MWT, label="Meteodyn WT (Manning et al. 2010)")
            if "measured" in ser:
                ax.plot([q[0] for q in ser["measured"]], [q[1] for q in ser["measured"]], "o", ms=4, color=TUNNEL,
                        label="tunnel (Manning et al. 2010)")
        else:   # the tunnel from the ERCOFTAC stations, interpolated to the plane
            pts = []
            for x in sorted({float(r["x"]) for r in st}):
                prof = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
                zt = zf * H
                for k in range(len(prof) - 1):
                    if prof[k][0] <= zt <= prof[k + 1][0]:
                        t = (zt - prof[k][0]) / (prof[k + 1][0] - prof[k][0])
                        pts.append((x / hl, ((1 - t) * prof[k][1] + t * prof[k + 1][1]) / uref))
                        break
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=4, color=TUNNEL, label="tunnel (ERCOFTAC)")
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xlim(-5, 6)
        ax.set_title(f"Z = {zf:g} H", fontsize=10)
        ax.grid(alpha=0.3)
        if n == 0:
            ax.legend(fontsize=7)
        if n % 4 == 0:
            ax.set_ylabel("U / U_in(0.8 H)")
    # Profiles at the stations.
    stations = [x for x in sorted({float(r["x"]) for r in st}) if -1.0 * hl <= x <= 3.2 * hl][:6]
    for n, x in enumerate(stations):
        ax = axs.flat[6 + n]
        prof = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
        u, z = column(cols, XC + x)
        zr = [zz - g(XC + x) for zz in z]
        k = [i for i, v in enumerate(zr) if v <= 1.2 * H]
        ax.plot([u[i] / uref for i in k], [zr[i] / H for i in k], color=GUST, label="GUST - CFD")
        ax.plot([p[1] / uref for p in prof], [p[0] / H for p in prof], "o", ms=3, color=TUNNEL,
                label="tunnel (ERCOFTAC)")
        if a.manning and figs:
            lab = f"X{x / hl:.2g}a".replace("X0a", "X0").replace("X-0a", "X0")
            ser = manning_series(a.manning, figs[1], lab)
            if "MWT" in ser:
                p = sorted(ser["MWT"], key=lambda q: q[1])
                ax.plot([q[0] for q in p], [q[1] for q in p], "--", color=MWT, label="Meteodyn WT")
        ax.axvline(0, color="0.7", lw=0.8)
        ax.set_ylim(0, 1.2)
        ax.set_title(f"X = {x / hl:.2f} a", fontsize=10)
        ax.grid(alpha=0.3)
        if n == 0:
            ax.legend(fontsize=7)
        if (6 + n) % 4 == 0:
            ax.set_ylabel("Z / H")
        ax.set_xlabel("U / U_in(0.8 H)")
    for n in range(6 + len(stations), 12):
        axs.flat[n].axis("off")
    fig.suptitle(f"RUSHIL Hill {a.hill} (a / H = {a.hill}), full scale H = 117 m", fontsize=12)
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
    c.add_argument("--hill", type=int, choices=HILLS, required=True)
    c.add_argument("--outdir", required=True)
    c.add_argument("--label")
    c.add_argument("--out")
    q = sub.add_parser("plot")
    q.add_argument("plotfile")
    q.add_argument("--hill", type=int, choices=HILLS, required=True)
    q.add_argument("--outdir", required=True)
    q.add_argument("--manning", help="digitized Manning et al. (2010) curves (figure,panel,series,x,y); not distributed")
    q.add_argument("--out", required=True)
    a = ap.parse_args()
    {"prepare": prepare, "compare": compare, "plot": plot}[a.cmd](a)


if __name__ == "__main__":
    main()
