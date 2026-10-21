#!/usr/bin/env python3
"""Wind speed 10 m above the ground on the Askervein lines (RT-2.6.1).

    askervein_lines.py sample PLOTFILE --terrain TERRAIN --lines LINES [--height H] [--out CSV]
    askervein_lines.py compare TEST.csv REFERENCE.csv [--rel_tol 0.05] [--lines A AA]
                       [--lee_d D --lee_rel_tol T]
    askervein_lines.py check PLOTFILE REFERENCE.csv --terrain TERRAIN --lines LINES --out CSV [compare options]

sample: the horizontal speed sqrt(u^2 + v^2) at H (default 10 m) above the
ground at each mast of LINES (columns line, d, x, y, zg: Tests/Real_Terrain/
Askervein_210/askervein_lines.csv). PLOTFILE is an ablsteady plotfile (u, v,
z) or an ERF one, or an ERF reference archive's average (x_velocity,
y_velocity, z_phys); a directory with an "average" plotfile inside is taken
as an archive. Both are sampled the same way: in each of the four cell
columns around the mast, linear in ln(z + z0) of the height above the ground z
between the two cell centres that bracket H, or below the first centre by the wall law
S ~ ln((z + z0)/z0) (--z0, default 0.03 m, the Askervein roughness of both
codes); the ground under a column is TERRAIN, ERF's structured terrain file,
bilinear at the column centre (within 1.4 m of the mesh's own ground on the
ERF grid). Then bilinear in x and y between the four columns. Writes "line,d,x,y,speed" rows; prints them.

A PLOTFILE that does not exist is taken as a prefix: the plotfile of that
prefix with the highest step is read (prefix00213, not prefix_level1_00213).

check: sample PLOTFILE into CSV, then compare CSV with REFERENCE.csv.

compare: each mast of the named lines in both files; passes (exit code 0) if
every |speed_test - speed_ref| <= rel_tol * speed_ref, and prints each
mast's relative difference and the largest. Masts are matched by line and d.
With --lee_d D, masts more than D m past the hilltop (d > D) are judged
against --lee_rel_tol instead (RT-2.6.1: the lee, where the advection scheme
sets the wake: ablsteady's second order against ERF's fifth).
"""
import argparse
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

NAMES = {"ablsteady": ("u", "v", "z"), "erf": ("x_velocity", "y_velocity", "z_phys")}


def read_terrain(path):
    """ERF's structured terrain file: nx ny, x values, y values, z (y fastest, as ERF's reader)."""
    with open(path) as fh:
        vals = fh.read().split()
    nx, ny = int(vals[0]), int(vals[1])
    xs = [float(v) for v in vals[2:2 + nx]]
    ys = [float(v) for v in vals[2 + nx:2 + nx + ny]]
    zs = [float(v) for v in vals[2 + nx + ny:2 + nx + ny + nx * ny]]
    if len(zs) != nx * ny:
        sys.exit(f"{path}: expected {nx * ny} heights, got {len(zs)}")
    return xs, ys, zs


def bracket(grid, v):
    """Index i with grid[i] <= v <= grid[i + 1] (clamped at the ends) and the weight of grid[i + 1]."""
    if v <= grid[0]:
        return 0, 0.0
    if v >= grid[-1]:
        return len(grid) - 2, 1.0
    lo, hi = 0, len(grid) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if grid[mid] <= v:
            lo = mid
        else:
            hi = mid
    return lo, (v - grid[lo]) / (grid[lo + 1] - grid[lo])


def ground(terrain, x, y):
    xs, ys, zs = terrain
    i, a = bracket(xs, x)
    j, b = bracket(ys, y)
    ny = len(ys)
    z = lambda ii, jj: zs[ii * ny + jj]
    return ((1 - a) * (1 - b) * z(i, j) + a * (1 - b) * z(i + 1, j)
            + (1 - a) * b * z(i, j + 1) + a * b * z(i + 1, j + 1))


def resolve(plotfile):
    """plotfile itself, or the highest-step plotfile prefix<digits> if it does not exist."""
    if os.path.exists(plotfile):
        return plotfile
    d, prefix = os.path.split(plotfile)
    d = d or "."
    steps = [(int(n[len(prefix):]), n) for n in os.listdir(d)
             if n.startswith(prefix) and n[len(prefix):].isdigit()]
    if not steps:
        sys.exit(f"{plotfile}: no such plotfile, nor any {prefix}<step> in {d}")
    return os.path.join(d, max(steps)[1])


class Columns:
    """Cell columns of a single-level plotfile, read box by box on demand."""

    def __init__(self, plotfile):
        plotfile = resolve(plotfile)
        self.path = plotfile
        if os.path.isdir(os.path.join(plotfile, "average")):
            plotfile = os.path.join(plotfile, "average")
        self.hdr = erf_plotfile._read_header(plotfile)
        names = self.hdr["names"]
        kind = "ablsteady" if all(n in names for n in NAMES["ablsteady"]) else "erf"
        for n in NAMES[kind]:
            if n not in names:
                sys.exit(f"{plotfile}: no field {n} (have {', '.join(names)})")
        self.comps = [names.index(n) for n in NAMES[kind]]
        self.kind = kind
        if self.hdr["finest"] != 0:
            sys.exit(f"{plotfile}: {self.hdr['finest'] + 1} levels; only single-level plotfiles are read")
        self.level_dir = os.path.join(plotfile, "Level_0")
        self.ncomp, _, self.boxes, self.fabs = erf_plotfile._read_cell_h(self.level_dir)
        lo, hi = self.hdr["lo"], self.hdr["hi"]
        for blo, bhi in self.boxes:
            if blo[2] != lo[2] or bhi[2] != hi[2]:
                sys.exit(f"{plotfile}: a box spans z {blo[2]}..{bhi[2]} of {lo[2]}..{hi[2]}; columns split in z are not read")
        # erf_plotfile.py reads 8-byte reals only (it takes the format's field
        # count as the size and "8 ..." as little endian): refuse float data,
        # whose byte order has 4 entries, "FAB ((8, (32 8 23 ...)),(4, (4 3 2 1)))".
        with open(os.path.join(self.level_dir, self.fabs[0][0]), "rb") as fh:
            fh.seek(self.fabs[0][1])
            head = fh.readline().decode(errors="replace")
        order = re.search(r"\),\((\d+), \(([\d ]+)\)\)\)", head)
        if order is None or int(order.group(1)) != 8:
            sys.exit(f"{plotfile}: not double precision; extern_erf/erf_plotfile.py reads only doubles correctly")
        self.cache = {}

    def column(self, i, j):
        """(u, v, z) lists over k of cell (i, j)."""
        for b, ((blo, bhi), (fname, off)) in enumerate(zip(self.boxes, self.fabs)):
            if blo[0] <= i <= bhi[0] and blo[1] <= j <= bhi[1]:
                if b not in self.cache:
                    self.cache[b] = erf_plotfile._read_fab(os.path.join(self.level_dir, fname), off, self.ncomp)
                flo, fhi, data = self.cache[b]
                fnx, fny, fnz = (fhi[d] - flo[d] + 1 for d in range(3))
                npts = fnx * fny * fnz
                out = []
                for c in self.comps:
                    out.append([data[c * npts + ((k - flo[2]) * fny + (j - flo[1])) * fnx + (i - flo[0])]
                                for k in range(blo[2], bhi[2] + 1)])
                return out
        sys.exit(f"cell ({i}, {j}) is in no box of the plotfile")


def speed_at(cols, terrain, x, y, height, z0):
    h = cols.hdr
    lo, dx = h["prob_lo"], h["dx"]
    nx, ny = h["hi"][0] - h["lo"][0] + 1, h["hi"][1] - h["lo"][1] + 1
    fx, fy = (x - lo[0]) / dx[0] - 0.5, (y - lo[1]) / dx[1] - 0.5
    i0, j0 = int(math.floor(fx)), int(math.floor(fy))
    if not (0 <= i0 < nx - 1 and 0 <= j0 < ny - 1):
        sys.exit(f"mast ({x}, {y}) lies outside the cell centres of the plotfile")
    a, b = fx - i0, fy - j0
    total = 0.0
    for di, dj, w in ((0, 0, (1 - a) * (1 - b)), (1, 0, a * (1 - b)), (0, 1, (1 - a) * b), (1, 1, a * b)):
        i, j = i0 + di, j0 + dj
        u, v, z = cols.column(i, j)
        zg = ground(terrain, lo[0] + (i + 0.5) * dx[0], lo[1] + (j + 0.5) * dx[1])
        agl = [zz - zg for zz in z]
        if height > agl[-1] or agl[0] <= 0:
            sys.exit(f"height {height} m is outside the column ({agl[0]:.2f} to {agl[-1]:.2f} m above the ground)")
        if height < agl[0]:
            # Below the first cell centre: the wall law of both codes, S ~ ln((z + z0) / z0).
            s0 = math.hypot(u[0], v[0])
            total += w * s0 * math.log((height + z0) / z0) / math.log((agl[0] + z0) / z0)
            continue
        k, _ = bracket(agl, height)
        # Linear in ln(z + z0) between the two centres: exact for a log profile
        # (linear in z reads up to 1.6 % low at 10 m between centres near 9 and 23 m).
        lz = [math.log(agl[k] + z0), math.log(agl[k + 1] + z0)]
        t = (math.log(height + z0) - lz[0]) / (lz[1] - lz[0])
        s0, s1 = math.hypot(u[k], v[k]), math.hypot(u[k + 1], v[k + 1])
        total += w * ((1 - t) * s0 + t * s1)
    return total


def read_rows(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))]


def cmd_sample(a):
    terrain = read_terrain(a.terrain)
    cols = Columns(a.plotfile)
    rows = []
    for r in read_rows(a.lines):
        x, y = float(r["x"]), float(r["y"])
        rows.append((r["line"], int(r["d"]), x, y, speed_at(cols, terrain, x, y, a.height, a.z0)))
    lines = [f"# {a.height:g} m above the ground; {os.path.basename(cols.path.rstrip('/'))} "
             f"({cols.kind}, t = {cols.hdr['time']:g})", "line,d,x,y,speed"]
    lines += ["%s,%d,%.2f,%.2f,%.6f" % r for r in rows]
    print("\n".join(lines))
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
    return 0


def cmd_compare(a):
    key = lambda r: (r["line"], int(r["d"]))
    test = {key(r): float(r["speed"]) for r in read_rows(a.test)}
    print(f"{a.test} against {a.reference}:")
    ref = {key(r): float(r["speed"]) for r in read_rows(a.reference)}
    worst, worst_lee, n, ok = 0.0, 0.0, 0, True
    for k in sorted(ref, key=lambda k: (k[0], k[1])):
        if k[0] not in a.lines:
            continue
        if k not in test:
            print(f"FAIL: line {k[0]} d = {k[1]} missing from {a.test}")
            return 1
        rel = (test[k] - ref[k]) / ref[k]
        n += 1
        lee = a.lee_d is not None and k[1] > a.lee_d
        if lee:
            worst_lee = max(worst_lee, abs(rel))
        else:
            worst = max(worst, abs(rel))
        bad = abs(rel) > (a.lee_rel_tol if lee else a.rel_tol)
        ok &= not bad
        print(f"  line {k[0]:>2} d = {k[1]:5d} m: {test[k]:7.3f} vs {ref[k]:7.3f} m/s ({100 * rel:+6.2f} %)"
              + (" lee" if lee else "") + ("  > tolerance" if bad else ""))
    if n == 0:
        print("FAIL: no masts compared")
        return 1
    lee = ""
    if a.lee_d is not None:
        lee = f"; lee (d > {a.lee_d:g} m) {100 * worst_lee:.2f} % (tolerance {100 * a.lee_rel_tol:g} %)"
    print(f"{'PASS' if ok else 'FAIL'}: {n} masts on lines {' '.join(a.lines)}, largest difference "
          f"{100 * worst:.2f} % (tolerance {100 * a.rel_tol:g} %){lee}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("plotfile")
    s.add_argument("--terrain", required=True)
    s.add_argument("--lines", required=True)
    s.add_argument("--height", type=float, default=10.0)
    s.add_argument("--z0", type=float, default=0.03, help="roughness of the log law below the first cell centre (m)")
    s.add_argument("--out")
    c = sub.add_parser("compare")
    c.add_argument("test")
    c.add_argument("reference")
    c.add_argument("--rel_tol", type=float, default=0.05)
    c.add_argument("--lines", nargs="+", default=["A", "AA"])
    c.add_argument("--lee_d", type=float, help="masts with d above this (m) take --lee_rel_tol")
    c.add_argument("--lee_rel_tol", type=float)
    k = sub.add_parser("check")
    k.add_argument("plotfile")
    k.add_argument("reference")
    k.add_argument("--terrain", required=True)
    k.add_argument("--lines", dest="lines_file", required=True)
    k.add_argument("--out", required=True)
    k.add_argument("--height", type=float, default=10.0)
    k.add_argument("--z0", type=float, default=0.03)
    k.add_argument("--rel_tol", type=float, default=0.05)
    k.add_argument("--compare_lines", dest="lines_cmp", nargs="+", default=["A", "AA"])
    k.add_argument("--lee_d", type=float)
    k.add_argument("--lee_rel_tol", type=float)
    a = ap.parse_args()
    if a.cmd == "check":
        b = argparse.Namespace(cmd="sample", plotfile=a.plotfile, terrain=a.terrain, lines=a.lines_file, height=a.height,
                               z0=a.z0, out=a.out)
        c = argparse.Namespace(test=a.out, reference=a.reference, rel_tol=a.rel_tol, lines=a.lines_cmp,
                               lee_d=a.lee_d, lee_rel_tol=a.lee_rel_tol)
        a, checks = b, c
    else:
        checks = None
    if a.cmd == "sample" or checks is not None:
        if not (a.height > 0 and math.isfinite(a.height)):
            sys.exit("--height must be finite and > 0")
        if not (a.z0 > 0 and math.isfinite(a.z0)):
            sys.exit("--z0 must be finite and > 0")
        rc = cmd_sample(a)
        if checks is None:
            return rc
        a = checks
    if not (a.rel_tol > 0 and math.isfinite(a.rel_tol)):
        sys.exit("--rel_tol must be finite and > 0")
    if (a.lee_d is None) != (a.lee_rel_tol is None):
        sys.exit("--lee_d and --lee_rel_tol go together")
    if a.lee_d is not None and not (math.isfinite(a.lee_d) and a.lee_rel_tol > 0 and math.isfinite(a.lee_rel_tol)):
        sys.exit("--lee_d must be finite and --lee_rel_tol finite and > 0")
    return cmd_compare(a)


if __name__ == "__main__":
    sys.exit(main())
