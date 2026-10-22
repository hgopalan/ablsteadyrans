#!/usr/bin/env python3
"""Run GUST with a column target wind and check that the 3-D run keeps it.

    check_column_target.py --target U --height Z --z0 Z0 [--mol L] -- COMMAND [ARGS ...]

Runs COMMAND (the application on Tests/Canonical_RANS/Column_Target/inputs_column_target) and checks:
  * the column reached the target ("target wind reached" in its output) and the log lists the restart inputs
    as input lines (forcing.pressure_gradient with its three values, the first the forcing the run took; with
    --mol also bc.surface_heat_flux);
  * the 3-D wind at --height, the planar mean of the last plotfile interpolated linearly in ln(z + z0) between
    the cell centres (as the column does), is the target to 0.5 % (the 3-D discretization differs slightly from
    the column's), and its direction along x to 1e-3 rad;
  * with --mol L: the 3-D surface's Obukhov length, -u*^3 theta0 / (kappa g q) with the mean u* of the turbulence
    plotfile's surf_hist.dat and the flux q the run took from the column, is L to 2 % (theta0 = 300 K,
    g = 9.81 m/s^2, kappa = 0.41).
Prints one row per check; the exit code is the verdict.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def main(argv):
    opts = {"--target": None, "--height": None, "--z0": None, "--mol": 0.0}
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] in opts and i + 1 < len(argv):
            opts[argv[i]] = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if any(v is None for v in opts.values()) or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    print("check_column_target: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(run.stdout)
    if run.returncode != 0:
        sys.stdout.write(run.stderr)
        raise SystemExit("check_column_target: the command failed with code %d" % run.returncode)
    ok = True

    def check(name, value, target, tol, kind):
        nonlocal ok
        err = abs(value - target)
        lim = tol * (abs(target) if kind == "rel" else 1.0)
        good = err <= lim
        ok = ok and good
        print("  %-52s measured %14.8e target %14.8e tol %9.2e %s" % (name, value, target, lim, "pass" if good else "FAIL"))

    out = run.stdout
    check("column reached the target", 1.0 if "target wind reached" in out else 0.0, 1.0, 0.0, "abs")
    # The restart line is an input line: three values, the first the forcing the run reports taking.
    pg = [ln.split("=")[1].split() for ln in out.splitlines() if ln.startswith("  forcing.pressure_gradient = ")]
    check("restart input forcing.pressure_gradient: 3 values", len(pg[0]) if pg else 0, 3, 0.0, "abs")
    took = [ln for ln in out.splitlines() if "the run takes geostrophic wind" in ln]
    if pg and len(pg[0]) == 3 and took:
        fx = float(took[0].split("pressure gradient (")[1].split(",")[0])
        check("restart input forcing.pressure_gradient x vs the forcing taken", float(pg[0][0]), fx, 1e-5, "rel")
    wrote = [ln.split()[-1] for ln in out.splitlines() if ln.startswith("gust: wrote ")]
    plt = [n for n in wrote if n.endswith(tuple("0123456789")) and "_plt" in n]
    turb = [n for n in wrote if "_turb" in n]
    if not plt or not turb:
        raise SystemExit("check_column_target: no plotfile written")
    z, prof, _ = erf_plotfile.planar_averages(plt[-1], ["u", "v"])
    zt, z0 = opts["--height"], opts["--z0"]
    k = max(j for j in range(len(z) - 1) if z[j] <= zt)
    t = (math.log(zt + z0) - math.log(z[k] + z0)) / (math.log(z[k + 1] + z0) - math.log(z[k] + z0))
    u = (1 - t) * prof["u"][k] + t * prof["u"][k + 1]
    v = (1 - t) * prof["v"][k] + t * prof["v"][k + 1]
    check("3-D speed at the target height vs target", math.hypot(u, v), opts["--target"], 5e-3, "rel")
    check("3-D direction at the target height (rad)", math.atan2(v, u), 0.0, 1e-3, "abs")
    if opts["--mol"] != 0:
        q = None
        for ln in out.splitlines():
            if ln.startswith("  bc.surface_heat_flux = "):
                q = float(ln.split("=")[1])
        check("restart input bc.surface_heat_flux listed", 0.0 if q is None else 1.0, 1.0, 0.0, "abs")
        if q is not None:
            us = erf_plotfile.read_surf_hist(os.path.join(turb[-1], "surf_hist.dat"))["u_star"]
            L3 = -us ** 3 * 300.0 / (0.41 * 9.81 * q)
            check("3-D surface Obukhov length vs --mol", L3, opts["--mol"], 2e-2, "rel")
    print("check_column_target:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
