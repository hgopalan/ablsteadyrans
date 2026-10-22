#!/usr/bin/env python3
"""Speed-up on the transects of Wagenbrenner et al. (2019, Atmosphere 10, 672):
Askervein (lines A, AA, B), Bolund (cases 1, 3, 4: lines B, A, B) and Big
Southern Butte (TSW, TWSW; observations not yet available).

    windninja2019.py prepare WNTOOLS OUTDIR
        From a checkout of github.com/nwagenbrenner/windninja-tools: the
        terrain of Bolund and Big Southern Butte as terrain.type = xyz files
        (bolund.xyz, bsb.xyz) and the observations as one CSV (masts.csv).
    windninja2019.py sample PLOTFILE --site SITE --case CASE --masts masts.csv
            (--terrain XYZ | --erf_terrain FILE) [--out CSV]
        The model's horizontal speed at each mast's height above the ground,
        and along each line every 5 m (Bolund) or 25 m (Askervein), as
        speed-up (S - U_ref) / U_ref.
    windninja2019.py plot CSV [CSV ...] --masts masts.csv --site SITE --case CASE --out PNG
        One panel per line: the observations (points) and each CSV (a curve,
        labelled by its file name).

Speed-up follows the paper: (S - U_ref) / U_ref with U_ref the input wind of
the case (Askervein 8.9 m/s at 10 m; Bolund 10.9, 8.7 and 7.6 m/s at 5 m;
Big Southern Butte 8.3 m/s at 3 m), for the observations and the model alike.
The ground under each column is the terrain file's (the same inverse-distance
weighting as GUST's terrain.type = xyz: 6 neighbours, power 2; or ERF's
grid, bilinear); the vertical interpolation is that of askervein_lines.py.
"""
import argparse
import csv
import math
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))   # askervein_lines.py, tiff_to_xyz.py
import askervein_lines as al  # noqa: E402

# site, case: reference speed (m/s), mast height above the ground (m), lines
# as name: (x0, y0, x1, y1) end points, d measured from the origin of the line's frame.
CASES = {
    ("askervein", "210"): dict(uref=8.9, agl=10.0, step=25.0, paper=("fig8", {"A": "LineA", "AA": "LineAA", "B": "LineB"})),
    ("bolund", "270"): dict(uref=10.9, agl=5.0, step=5.0, lines={"B": (-200.0, 0.0, 200.0, 0.0)},
                            paper=("fig12", {"B": "Case1"})),
    # Case 3: the paper's axis is the mast's x, not the distance along line A.
    ("bolund", "239"): dict(uref=8.7, agl=5.0, step=5.0, lines={"A": (-200.0, -113.2, 200.0, 115.5)}, d_is_x=True,
                            paper=("fig12", {"A": "Case3"})),
    ("bolund", "90"): dict(uref=7.6, agl=5.0, step=5.0, lines={"B": (-200.0, 0.0, 200.0, 0.0)},
                           paper=("fig12", {"B": "Case4"})),
    ("bsb", "222"): dict(uref=8.3, agl=3.0, step=50.0, paper=("fig15", {"TSW": "TSW", "TWSW": "TWSW"})),
}
# The paper's curves drawn by plot --paper (Wagenbrenner et al. 2019, digitized): the recommended CFD
# setting (KELU: k-epsilon, linear upwind), the mass-consistent solver (COM) and the LES.
PAPER_SERIES = {"KELU": dict(color="tab:brown", ls="-.", label="WindNinja CFD KELU (paper)"),
                "COM": dict(color="k", ls="--", label="WindNinja COM (paper)")}   # the paper's LES is not plotted
OURS = ["tab:blue", "tab:green", "tab:purple", "tab:red", "tab:cyan"]  # our runs; the paper's curves use PAPER_SERIES
BSB_CORNER = (326705.0, 4795665.0)  # big_butte.asc's xllcorner, yllcorner (UTM 12N)
BOLUND_LINES = {"270": ("B", ["M7", "M6", "M3", "M8"]), "90": ("B", ["M7", "M6", "M3", "M8"]),
                "239": ("A", ["M0", "M1", "M2", "M3", "M4"])}
FIELDS = ["site", "case", "line", "mast", "d", "x", "y", "agl", "speed", "speedup"]


# ---------------------------------------------------------------------------
def prepare(a):
    src = a.wntools
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    # Askervein: points.asc, blocks A, AA, B of "x y speed d" (Taylor and Teunissen; the frame of extern_erf/askervein.txt).
    blocks, cur = [], []
    for line in open(os.path.join(src, "data", "observed_data", "points.asc")):
        if line.strip():
            cur.append([float(v) for v in line.split()])
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)
    assert len(blocks) == 3, f"points.asc: {len(blocks)} blocks, expected A, AA, B"
    uref = CASES[("askervein", "210")]["uref"]
    for name, blk in zip(("A", "AA", "B"), blocks):
        for x, y, s, d in blk:
            if s > 0:
                rows.append(dict(site="askervein", case="210", line=name, mast=f"{name}{int(d)}", d=d, x=x, y=y,
                                 agl=10.0, speed=s, speedup=(s - uref) / uref))
    # Bolund: the 5 m masts of the blind comparison (Bechmann et al. 2011), speed = (vel/u*) u*.
    with zipfile.ZipFile(os.path.join(src, "data", "bolund", "Bolund_Measurements.zip")) as z:
        for case, (lname, masts) in BOLUND_LINES.items():
            text = z.read(f"Measurements/Dir_{case}.dat").decode().splitlines()
            x0, y0, x1, y1 = CASES[("bolund", case)]["lines"][lname]
            ux, uy = (x1 - x0), (y1 - y0)
            n = math.hypot(ux, uy)
            uref = CASES[("bolund", case)]["uref"]
            for t in text[1:]:
                f = t.split()
                mid = f[0]
                if not mid.endswith("Z05S") or mid[:-4] not in masts:
                    continue
                x, y, ustar, vel = float(f[3]), float(f[4]), float(f[7]), float(f[8])
                s = vel * ustar
                d = x if CASES[("bolund", case)].get("d_is_x") else (x * ux + y * uy) / n
                rows.append(dict(site="bolund", case=case, line=lname, mast=mid, d=d, x=x, y=y,
                                 agl=5.0, speed=s, speedup=(s - uref) / uref))
    # Big Southern Butte: the measured 10-min means of bsb_obs_1700.csv (every site with data), in the frame of
    # bsb.xyz (x, y from big_butte.asc's lower-left corner). The sensors on TSW and TWSW get their distance along
    # the line (the projection on Figure 14a's digitized transect, distance 0 at its south-west end); the others
    # are on no line ("-": in the statistics, not the transect plots).
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Tests", "Real_Terrain", "WindNinja2019")
    nocomment = lambda path: (l for l in open(path) if l[0] != "#")  # noqa: E731
    ends = {}
    for r in csv.DictReader(nocomment(os.path.join(here, "bsb_transects_fig14.csv"))):
        if r["kind"] in ("start", "end"):
            ends.setdefault(r["transect"], {})[r["kind"]] = (float(r["easting"]), float(r["northing"]))
    uref = CASES[("bsb", "222")]["uref"]
    for r in csv.DictReader(nocomment(os.path.join(here, "bsb_obs_1700.csv"))):
        e, n, sp = float(r["easting"]), float(r["northing"]), float(r["speed"])
        line, d = "-", 0.0
        for name in ("TSW", "TWSW"):
            if r["plot"].startswith(name) and r["plot"][len(name):].split("_")[0].isdigit():
                (e0, n0), (e1, n1) = ends[name]["start"], ends[name]["end"]
                L = math.hypot(e1 - e0, n1 - n0)
                line, d = name, ((e - e0) * (e1 - e0) + (n - n0) * (n1 - n0)) / L
        rows.append(dict(site="bsb", case="222", line=line, mast=r["plot"], d=d, x=e - BSB_CORNER[0],
                         y=n - BSB_CORNER[1], agl=3.0, speed=sp, speedup=sp / uref - 1))
    with open(os.path.join(a.outdir, "masts.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"{a.outdir}/masts.csv: {len(rows)} observations")
    # Terrain.
    import tiff_to_xyz
    tiff_to_xyz.main([os.path.join(src, "data", "bolund", "bolund_hill_extended.tif"),
                      os.path.join(a.outdir, "bolund.xyz"), "4", "0"])
    asc_to_xyz(os.path.join(src, "test_cases", "big_butte.asc"), os.path.join(a.outdir, "bsb.xyz"))


def read_asc(src):
    """An ESRI ASCII grid: (z rows north first, xll, yll, cellsize, nodata)."""
    import numpy as np
    hdr = {}
    with open(src) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    z = np.loadtxt(src, skiprows=6)
    return z, hdr["xllcorner"], hdr["yllcorner"], hdr["cellsize"], hdr.get("nodata_value", -9999.0)


def smooth_asc(src, out, sigma, cap=None, soft=100.0):
    """big_butte.asc smoothed by a Gaussian of standard deviation sigma (m) and, with cap (m ASL), limited by
    the soft cap z -> cap - soft ln(1 + exp((cap - z)/soft)) (z below the cap unchanged, above it -> cap):
    the form that reproduces WindNinja's CFD mesh terrain of the Butte (Figure 14b,c; sigma 120 m, cap 2050 m,
    soft 150 m: 25 m and 24 m RMS along TSW and TWSW, against 41-45 m for the best Gaussian alone, which widens
    the mountain). As xyz in bsb.xyz's frame and with bsb.xyz's z shift (the unsmoothed minimum)."""
    import numpy as np
    from scipy.ndimage import gaussian_filter
    z, _, _, h, nod = read_asc(src)
    zmin = z[z != nod].min()
    zs = gaussian_filter(z, sigma / h, mode="nearest") if sigma > 0 else z
    if cap is not None:
        zs = cap - soft * np.logaddexp(0.0, (cap - zs) / soft)
    ny, nx = z.shape
    with open(out, "w") as f:
        f.write(f"# {src} smoothed by a Gaussian of sigma {sigma} m, soft cap {cap} m ASL (softness {soft} m); "
                f"x, y from the lower-left corner; z - {zmin}\n")
        for j in range(ny):
            y = (ny - j - 0.5) * h
            for i in range(nx):
                f.write(f"{(i + 0.5) * h:.1f} {y:.1f} {zs[j, i] - zmin:.2f}\n")
    print(f"{out}: sigma {sigma} m, max {zs.max():.1f} m ASL (unsmoothed {z.max():.1f})")


def fit_smoothing(a):
    """The Gaussian width that best reproduces WindNinja's CFD terrain along TSW and TWSW (Figure 14b,c):
    the RMS of (smoothed DEM - CFD line) along both transects for each sigma."""
    import numpy as np
    from scipy.ndimage import gaussian_filter, map_coordinates
    z, xll, yll, h, nod = read_asc(a.asc)
    ny = z.shape[0]
    ends = {}
    for r in csv.DictReader(l for l in open(a.transects) if l[0] != "#"):
        if r["kind"] in ("start", "end"):
            ends.setdefault(r["transect"], {})[r["kind"]] = (float(r["easting"]), float(r["northing"]))
    prof = [r for r in csv.DictReader(l for l in open(a.profiles) if l[0] != "#")]
    print("sigma (m) | RMS TSW | RMS TWSW | max TSW | max TWSW  (smoothed DEM minus CFD line, m)")
    for sigma in [float(v) for v in a.sigmas.split(",")]:
        zs = gaussian_filter(z, sigma / h, mode="nearest") if sigma > 0 else z
        out = []
        for name in ("TSW", "TWSW"):
            (e0, n0), (e1, n1) = ends[name]["start"], ends[name]["end"]
            L = math.hypot(e1 - e0, n1 - n0)
            pts = [(float(r["x"]), float(r["elev"])) for r in prof if r["panel"] == name and r["series"] == "CFD"]
            d = []
            for x, el in pts:
                e, n = e0 + (e1 - e0) * x / L, n0 + (n1 - n0) * x / L
                col, row = (e - xll) / h - 0.5, ny - (n - yll) / h - 0.5
                d.append(float(map_coordinates(zs, [[row], [col]], order=1)[0]) - el)
            d = np.array(d)
            out.append((math.sqrt((d ** 2).mean()), d.max(), d.min()))
        print(f"{sigma:8.0f} | {out[0][0]:7.1f} | {out[1][0]:8.1f} | {out[0][1]:+7.1f} | {out[1][1]:+8.1f}")


def asc_to_xyz(src, out):
    """An ESRI ASCII grid as xyz points: x, y from the lower-left corner (local metres), z minus its minimum."""
    hdr, vals = {}, []
    with open(src) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[0].lower() in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "nodata_value"):
                hdr[p[0].lower()] = float(p[1])
            else:
                vals.extend(float(v) for v in p)
    nx, ny, h = int(hdr["ncols"]), int(hdr["nrows"]), hdr["cellsize"]
    nod = hdr.get("nodata_value", -9999.0)
    good = [v for v in vals if v != nod]
    zmin = min(good)
    with open(out, "w") as f:
        f.write(f"# {src}: {nx} x {ny} cells of {h} m; x, y from the corner ({hdr['xllcorner']}, {hdr['yllcorner']}); "
                f"z - {zmin}\n")
        for j in range(ny):
            y = (ny - j - 0.5) * h
            for i in range(nx):
                z = vals[j * nx + i]
                if z != nod:
                    f.write(f"{(i + 0.5) * h:.1f} {y:.1f} {z - zmin:.2f}\n")
    print(f"{out}: {len(good)} points, {nx * h:.0f} x {ny * h:.0f} m, relief {max(good) - zmin:.1f} m")


# ---------------------------------------------------------------------------
class Idw:
    """GUST's terrain.type = xyz interpolation: inverse distance, k neighbours, power p."""

    def __init__(self, path, k=6, p=2.0):
        import numpy as np
        from scipy.spatial import cKDTree
        pts = []
        for line in open(path):
            line = line.split("#")[0].replace(",", " ").split()
            if len(line) == 3:
                try:
                    pts.append([float(v) for v in line])
                except ValueError:
                    pass
        self.p = np.array(pts)
        self.tree = cKDTree(self.p[:, :2])
        self.k, self.pw = k, p

    def __call__(self, x, y):
        d, i = self.tree.query([x, y], k=self.k)
        if d[0] == 0:
            return float(self.p[i[0], 2])
        w = 1.0 / d ** self.pw
        return float((w * self.p[i, 2]).sum() / w.sum())


def read_masts(path, site, case):
    return [r for r in csv.DictReader(open(path)) if r["site"] == site and r["case"] == case]


def sample(a):
    c = CASES[(a.site, a.case)]
    cols = al.Columns(a.plotfile)
    if a.terrain:
        idw = Idw(a.terrain)
        hx, hy = cols.hdr["dx"][0] / 2, cols.hdr["dx"][1] / 2

        def ground(x, y):
            # The mesh's ground under the column (x, y is its centre): the mean of the four corner nodes'
            # heights, as GUST places them (IDW at the nodes). IDW at the centre itself differs by
            # more than half a 1 m first layer on steep 30 m terrain.
            return 0.25 * sum(idw(x + sx * hx, y + sy * hy) for sx in (-1, 1) for sy in (-1, 1))
    else:
        grid = al.read_terrain(a.erf_terrain)
        ground = lambda x, y: al.ground(grid, x, y)  # noqa: E731
    z0 = a.z0
    out = []
    masts = read_masts(a.masts, a.site, a.case)
    for r in masts:
        x, y = float(r["x"]), float(r["y"])
        try:
            s = al.speed_at(cols, ground, x, y, c["agl"], z0)
        except SystemExit:  # a site outside the plotfile's cell centres (far stations of the Butte network)
            continue
        out.append(dict(site=a.site, case=a.case, line=r["line"], mast=r["mast"], d=float(r["d"]), x=x, y=y,
                        agl=c["agl"], speed=s, speedup=(s - c["uref"]) / c["uref"]))
    # Dense points along each line: the given end points, or the masts' own line (first to last mast, extended 20 %).
    lines = dict(c.get("lines", {}))
    for name in sorted({r["line"] for r in masts}):
        if name in lines:
            continue
        m = [r for r in masts if r["line"] == name]
        m.sort(key=lambda r: float(r["d"]))
        (xa, ya, da), (xb, yb, db) = [(float(r["x"]), float(r["y"]), float(r["d"])) for r in (m[0], m[-1])]
        L = db - da
        lines[name] = (xa - 0.2 * (xb - xa), ya - 0.2 * (yb - ya), xb + 0.2 * (xb - xa), yb + 0.2 * (yb - ya),
                       da - 0.2 * L, db + 0.2 * L)
    for name, seg in lines.items():
        x0, y0, x1, y1 = seg[:4]
        n = math.hypot(x1 - x0, y1 - y0)
        if len(seg) == 6:  # the masts' own line: d runs between the given values
            d0, d1 = seg[4], seg[5]
        else:  # a line through the frame's origin: d is the signed distance from it
            d0 = (x0 * (x1 - x0) + y0 * (y1 - y0)) / n
            d1 = d0 + n
        steps = int(n // c["step"])
        for k in range(steps + 1):
            t = k / steps
            x, y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
            try:
                s = al.speed_at(cols, ground, x, y, c["agl"], z0)
            except SystemExit:  # outside the plotfile's cell centres
                continue
            out.append(dict(site=a.site, case=a.case, line=name, mast="", d=x if c.get("d_is_x") else d0 + t * (d1 - d0), x=x, y=y,
                            agl=c["agl"], speed=s, speedup=(s - c["uref"]) / c["uref"]))
    f = open(a.out, "w", newline="") if a.out else sys.stdout
    w = csv.DictWriter(f, fieldnames=FIELDS)
    w.writeheader()
    for r in out:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    if a.out:
        f.close()
        print(f"{a.out}: {len(out)} rows")


def legend_label(path):
    """GUST - Fast (the mass-consistent run) or GUST - CFD (the RANS run), from the sample's file name
    <case>[_stable][_100k]_<mc|mc_stable|rans_<limiter>>.csv, with the variant in brackets."""
    name = os.path.splitext(os.path.basename(path))[0]
    extra = []
    if "_rans__" in name:  # a sweep variant (sweep.sh): GUST - CFD with one change
        return f"GUST - CFD ({name.split('_rans__')[1]})"
    if "_mc" in name:
        label = "GUST - Fast"
        if name.endswith("_mc_stable"):
            extra.append("slightly stable")
    elif "_rans_" in name:
        label = "GUST - CFD"
        if "_stable_" in name:
            extra.append("slightly stable")
        lim = name.split("_rans_")[1]
        if lim not in ("van_leer", "recommended"):   # recommended: the recommended GUST - CFD setting
            extra.append(lim)
    else:
        return name
    if "_100k_" in name:
        extra.append("100 K cells")
    return label + (f" ({', '.join(extra)})" if extra else "")


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    obs = read_masts(a.masts, a.site, a.case)
    names = [n for n in dict.fromkeys(r["line"] for r in obs) if n != "-"]  # masts.csv's order; "-": no line
    fig, axes = plt.subplots(1, len(names), figsize=(4.2 * len(names), 3.6), squeeze=False)
    for ax, name in zip(axes[0], names):
        for path in a.csv:
            rows = [r for r in csv.DictReader(open(path)) if r["line"] == name and not r["mast"]]
            rows.sort(key=lambda r: float(r["d"]))
            ax.plot([float(r["d"]) for r in rows], [float(r["speedup"]) for r in rows], lw=1.6,
                    color=OURS[a.csv.index(path)] if len(a.csv) <= len(OURS) else plt.cm.tab20(a.csv.index(path) % 20),
                    label=legend_label(path))
        if a.paper:
            pfig, panels = CASES[(a.site, a.case)]["paper"]
            for series, style in PAPER_SERIES.items():
                pts = sorted((float(r["x"]), float(r["y"])) for r in csv.DictReader(l for l in open(a.paper) if l[0] != "#")
                             if r["figure"] == pfig and r["panel"] == panels.get(name) and r["series"] == series)
                if pts:
                    # Break the curve where the digitized points leave a gap (another curve covered it in
                    # the figure): more than 5 % of the curve's x range between neighbours.
                    gap = 0.05 * (pts[-1][0] - pts[0][0])
                    xs, ys = [pts[0][0]], [pts[0][1]]
                    for (x0, _), (x1, y1) in zip(pts, pts[1:]):
                        if x1 - x0 > gap:
                            xs.append(float("nan"))
                            ys.append(float("nan"))
                        xs.append(x1)
                        ys.append(y1)
                    ax.plot(xs, ys, lw=1, **style)
        o = [r for r in obs if r["line"] == name]
        ax.plot([float(r["d"]) for r in o], [float(r["speedup"]) for r in o], "ko", label="observed")
        ax.axhline(0, color="0.8", lw=0.8)
        ax.set_title(f"{a.site} {a.case}: line {name}")
        ax.set_xlabel("distance (m)")
        ax.set_ylabel("speed-up")
        ax.grid(alpha=0.3)
    axes[0][-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(a.out)


# Mast classes for the error statistics, as the paper's windward (w), ridge-top (r) and lee (l) groups.
# Askervein: line B runs along the ridge (r); on lines A and AA the masts within 100 m of the hilltop are r,
# upwind of them w, downwind l. Bolund: by the case's wind. Big Southern Butte: TSW7 (the summit) is r;
# on each line the sensors before the summit are w, after it l (TWSW has none on the summit).
BOLUND_CLASS = {"270": {"M7": "w", "M6": "r", "M3": "r", "M8": "l"},
                "239": {"M0": "w", "M1": "w", "M2": "r", "M3": "r", "M4": "l"},
                "90": {"M8": "w", "M3": "r", "M6": "r", "M7": "l"}}
_BSB = {}


def BSB_GROUND(x, y):
    """The ground of big_butte.asc at (x, y) in bsb.xyz's frame, the summit height and position (all m)."""
    if not _BSB:
        import numpy as np
        src = os.environ.get("BSB_ASC", "big_butte.asc")   # windninja-tools' test_cases/big_butte.asc
        if not os.path.exists(src):
            sys.exit(f"{src}: not found; set BSB_ASC to windninja-tools' test_cases/big_butte.asc (Butte mast classes)")
        z, _, _, h, _ = read_asc(src)
        j, i = np.unravel_index(np.argmax(z), z.shape)
        _BSB.update(z=z, h=h, zs=float(z.max()), xs=(i + 0.5) * h, ys=(z.shape[0] - j - 0.5) * h)
    z, h = _BSB["z"], _BSB["h"]
    i = min(max(int(x / h), 0), z.shape[1] - 1)
    j = min(max(z.shape[0] - 1 - int(y / h), 0), z.shape[0] - 1)
    return float(z[j, i]), _BSB["zs"], (_BSB["xs"], _BSB["ys"])


def mast_class(r):
    d = float(r["d"])
    if r["site"] == "askervein":
        return "r" if r["line"] == "B" or abs(d) <= 100 else ("w" if d < 0 else "l")
    if r["site"] == "bolund":
        return BOLUND_CLASS[r["case"]][r["mast"][:2]]
    # Big Southern Butte: ridge-top within 150 m of the summit's height (bsb.xyz), else windward or lee by the
    # side of the summit along the wind (from 222 degrees, towards 48 degrees counter-clockwise from x).
    zg, zs, (xs, ys) = BSB_GROUND(float(r["x"]), float(r["y"]))
    if zg >= zs - 150:
        return "r"
    along = (float(r["x"]) - xs) * math.cos(math.radians(48)) + (float(r["y"]) - ys) * math.sin(math.radians(48))
    return "w" if along < 0 else "l"


def paper_at(paper_rows, site, case, line, d):
    """The paper's digitized curve at distance d (linear between the two nearest points), or None outside it."""
    pfig, panels = CASES[(site, case)]["paper"]
    pts = sorted((float(r["x"]), float(r["y"])) for r in paper_rows
                 if r["figure"] == pfig and r["panel"] == panels.get(line))
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= d <= x1 and x1 - x0 <= 0.05 * (pts[-1][0] - pts[0][0]):  # no interpolation across a gap
            return y0 + (y1 - y0) * (d - x0) / (x1 - x0) if x1 > x0 else y0
    return None


def stats(a):
    """RMSE, MBE (m/s) and MAPE (%) of wind speed at the masts, per class, for each CSV and the paper's curves."""
    obs = [r for r in csv.DictReader(open(a.masts))]
    paper = [r for r in csv.DictReader(l for l in open(a.paper) if l[0] != "#")] if a.paper else []
    models = {}
    for path in a.csv:
        name = os.path.splitext(os.path.basename(path))[0]
        rows = {(r["site"], r["case"], r["mast"]): r for r in csv.DictReader(open(path)) if r["mast"]}
        models[name] = lambda o, rows=rows: (float(rows[(o["site"], o["case"], o["mast"])]["speedup"])
                                             if (o["site"], o["case"], o["mast"]) in rows else None)
    for series in ("KELU", "COM"):
        if paper:
            models[f"paper {series}"] = lambda o, series=series: paper_at(
                [r for r in paper if r["series"] == series], o["site"], o["case"], o["line"], float(o["d"]))
    out = []
    for name, f in models.items():
        groups = {}
        for o in obs:
            if a.site and o["site"] != a.site:
                continue
            su = f(o)
            if su is None:
                continue
            uref = CASES[(o["site"], o["case"])]["uref"]
            sp, so = uref * (1 + su), float(o["speed"])
            for g in (mast_class(o), "all"):
                groups.setdefault(g, []).append((sp - so, so))
        for g in ("w", "r", "l", "all"):
            e = groups.get(g, [])
            if not e:
                continue
            rmse = math.sqrt(sum(x * x for x, _ in e) / len(e))
            mbe = sum(x for x, _ in e) / len(e)
            mape = 100 * sum(abs(x) / so for x, so in e) / len(e)
            out.append((name, g, len(e), rmse, mbe, mape))
    print("| run | class | n | RMSE (m/s) | MBE (m/s) | MAPE (%) |")
    print("|---|---|---|---|---|---|")
    for name, g, n, rmse, mbe, mape in out:
        print(f"| {name} | {g} | {n} | {rmse:.2f} | {mbe:+.2f} | {mape:.0f} |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("wntools")
    p.add_argument("outdir")
    s = sub.add_parser("sample")
    s.add_argument("plotfile")
    s.add_argument("--site", required=True)
    s.add_argument("--case", required=True)
    s.add_argument("--masts", required=True)
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--terrain", help="terrain.type = xyz file")
    g.add_argument("--erf_terrain", help="ERF's structured terrain file")
    s.add_argument("--z0", type=float, default=0.01, help="roughness of the wall law below the first centre (m)")
    s.add_argument("--out")
    q = sub.add_parser("plot")
    q.add_argument("csv", nargs="+")
    q.add_argument("--masts", required=True)
    q.add_argument("--site", required=True)
    q.add_argument("--case", required=True)
    q.add_argument("--paper", help="digitized curves of the paper (figure,panel,series,x,y)")
    q.add_argument("--out", required=True)
    t = sub.add_parser("stats")
    t.add_argument("csv", nargs="+")
    t.add_argument("--masts", required=True)
    t.add_argument("--paper")
    t.add_argument("--site", help="one site only")
    m = sub.add_parser("smooth")
    m.add_argument("asc")
    m.add_argument("out")
    m.add_argument("--sigma", type=float, required=True, help="Gaussian standard deviation (m)")
    m.add_argument("--cap", type=float, help="soft cap on the elevation (m ASL)")
    m.add_argument("--soft", type=float, default=100.0, help="width of the soft cap (m)")
    fz = sub.add_parser("fit_smoothing")
    fz.add_argument("asc")
    fz.add_argument("--transects", required=True)
    fz.add_argument("--profiles", required=True)
    fz.add_argument("--sigmas", default="0,30,60,90,120,150,200,250,300,400")
    a = ap.parse_args()
    {"prepare": prepare, "sample": sample, "plot": plot, "stats": stats,
     "smooth": lambda a: smooth_asc(a.asc, a.out, a.sigma, a.cap, a.soft), "fit_smoothing": fit_smoothing}[a.cmd](a)


if __name__ == "__main__":
    main()
