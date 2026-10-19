#!/usr/bin/env python3
"""Inflow table from a steady precursor column (bc.inflow_profile = file).

    make_inflow_table.py PLOTFILE TURB_PLOTFILE --ug UG --vg VG
                         [--ustar U] [--zref Z] [--out TABLE]

PLOTFILE is the precursor's steady plotfile (u, v, theta, tke, z; a flat
doubly periodic column such as Tests/Decks/inputs_precursor_askervein_cnbl),
TURB_PLOTFILE its turbulence plotfile (surf_hist.dat gives u*), UG and VG the
geostrophic wind the precursor ran with.

Prints the planar-mean wind direction at --zref (default 10 m, degrees from
+x) and the mean friction velocity, and the geostrophic wind that turns the
--zref wind onto +x and scales u* to --ustar (default 0.6 m/s): the wind
turns with the geostrophic wind, and u* scales with the geostrophic speed up
to the slow change of the geostrophic drag coefficient, so a few reruns
converge. With --out, writes the planar means as "z u v theta k" rows, a
z = 0 row first (u = v = 0, theta and k of the lowest cell), heights above
the ground at the cell centres.

With --check TABLE (regression test of the committed table), the exit code
is the verdict of: the --zref wind within --dir-tol degrees of +x, u* within
--ustar-tol (relative) of --ustar, and every row of TABLE within 1e-3 of the
precursor's (u and v relative to the geostrophic speed, theta in K, k
relative to the largest k), and PLOTFILE's step equal to the step in TABLE's
header. Under a capping inversion theta has no fixed point (the inversion
entrains without end), so the step at the steady verdict, which depends on
rounding, moves theta; the check therefore reads a run stopped at the
table's own step. A PLOTFILE or TURB_PLOTFILE that does not exist is taken
as a prefix, and the last such file (highest step) is read.
"""
import argparse
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def interp(z, f, zq):
    """Linear interpolation of f(z) at zq (z increasing)."""
    if zq <= z[0]:
        return f[0]
    for k in range(1, len(z)):
        if zq <= z[k]:
            t = (zq - z[k - 1]) / (z[k] - z[k - 1])
            return f[k - 1] + t * (f[k] - f[k - 1])
    return f[-1]


def latest(path):
    """PATH itself, or the last PATH* (plotfile prefix) by step number."""
    if os.path.exists(path):
        return path
    hits = [p for p in glob.glob(path + "*") if p[len(path):].isdigit()]
    if not hits:
        raise SystemExit("make_inflow_table: no %s or %s<step>" % (path, path))
    return max(hits, key=lambda p: int(p[len(path):]))


def step_of(name):
    """The step number at the end of a plotfile name, or None."""
    m = re.search(r"(\d+)$", os.path.basename(os.path.normpath(name)))
    return int(m.group(1)) if m else None


def table_step(path):
    """The plotfile step named in a table's header ("planar means of NAME ("), or None."""
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith("#"):
                return None
            m = re.search(r"planar means of (\S+) \(", ln)
            if m:
                return step_of(m.group(1))
    return None


def read_table(path):
    """Rows of a "z u v theta k" table (comments after #)."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.split("#", 1)[0].split()
            if ln:
                rows.append([float(x) for x in ln])
    return rows


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("turb_plotfile")
    ap.add_argument("--ug", type=float, required=True)
    ap.add_argument("--vg", type=float, required=True)
    ap.add_argument("--ustar", type=float, default=0.6)
    ap.add_argument("--zref", type=float, default=10.0)
    ap.add_argument("--out")
    ap.add_argument("--check")
    ap.add_argument("--dir-tol", type=float, default=0.01)
    ap.add_argument("--ustar-tol", type=float, default=1e-3)
    a = ap.parse_args(argv)
    if not (a.ustar > 0 and a.zref > 0 and a.dir_tol > 0 and a.ustar_tol > 0):
        raise SystemExit("make_inflow_table: --ustar, --zref, --dir-tol and --ustar-tol must be > 0")
    a.plotfile, a.turb_plotfile = latest(a.plotfile), latest(a.turb_plotfile)

    _, avg, _ = erf_plotfile.planar_averages(a.plotfile, ["u", "v", "theta", "tke", "z"])
    z, u, v, th, k = avg["z"], avg["u"], avg["v"], avg["theta"], avg["tke"]
    if any(z[i + 1] <= z[i] for i in range(len(z) - 1)) or z[0] <= 0:
        raise SystemExit("make_inflow_table: the plotfile's z is not increasing from above 0 (flat column expected)")
    us = erf_plotfile.read_surf_hist(os.path.join(a.turb_plotfile, "surf_hist.dat"))["u_star"]

    ur, vr = interp(z, u, a.zref), interp(z, v, a.zref)
    alpha = math.degrees(math.atan2(vr, ur))
    g = math.hypot(a.ug, a.vg)
    gdir = math.degrees(math.atan2(a.vg, a.ug))
    gnew = g * a.ustar / us
    dnew = math.radians(gdir - alpha)
    print("wind direction at %g m: %+.4f deg from +x" % (a.zref, alpha))
    print("u* = %.5f m/s (target %.5f, %+.3f %%)" % (us, a.ustar, 100 * (us / a.ustar - 1)))
    print("geostrophic wind %.4f m/s at %+.4f deg; next: forcing.ug = %.6f forcing.vg = %.6f"
          % (g, gdir, gnew * math.cos(dnew), gnew * math.sin(dnew)))

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("# z u v theta k: planar means of %s (ug %.6f, vg %.6f; u* %.5f m/s; "
                     "wind at %g m %+.4f deg)\n" % (os.path.basename(os.path.normpath(a.plotfile)),
                                                    a.ug, a.vg, us, a.zref, alpha))
            fh.write("%.6f %.9e %.9e %.9e %.9e\n" % (0.0, 0.0, 0.0, th[0], k[0]))
            for i in range(len(z)):
                fh.write("%.6f %.9e %.9e %.9e %.9e\n" % (z[i], u[i], v[i], th[i], k[i]))
        print("wrote", a.out)

    if a.check:
        rows = read_table(a.check)
        mine = [[0.0, 0.0, 0.0, th[0], k[0]]] + [[z[i], u[i], v[i], th[i], k[i]] for i in range(len(z))]
        ok = True

        def check(name, value, tol):
            nonlocal ok
            good = abs(value) <= tol
            ok = ok and good
            print("  %-34s %12.4e  (tol %.1e)  %s" % (name, value, tol, "ok" if good else "FAIL"))

        tstep, pstep = table_step(a.check), step_of(a.plotfile)
        same_step = tstep is not None and tstep == pstep
        ok = ok and same_step
        print("  %-34s %12s  (table %s)  %s" % ("plotfile step", pstep, tstep, "ok" if same_step else "FAIL"))
        check("10 m direction (deg)", alpha, a.dir_tol)
        check("u*/target - 1", us / a.ustar - 1, a.ustar_tol)
        if len(rows) != len(mine) or any(len(r) != 5 for r in rows):
            print("  %s has %d rows, the precursor %d  FAIL" % (a.check, len(rows), len(mine)))
            ok = False
        else:
            kmax = max(k)
            check("max |dz| (m)", max(abs(r[0] - m[0]) for r, m in zip(rows, mine)), 1e-3)
            check("max |du|, |dv| / G", max(max(abs(r[1] - m[1]), abs(r[2] - m[2])) for r, m in zip(rows, mine)) / g, 1e-3)
            check("max |dtheta| (K)", max(abs(r[3] - m[3]) for r, m in zip(rows, mine)), 1e-3)
            check("max |dk| / max k", max(abs(r[4] - m[4]) for r, m in zip(rows, mine)) / kmax, 1e-3)
        print("make_inflow_table:", "all checks passed" if ok else "FAILED")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
