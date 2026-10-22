#!/usr/bin/env python3
"""Run GUST on the 1-D neutral column deck and check the steady state (rt_steady_column).

    check_steady_column.py --ustar U --z0 Z0 --force F --height H --dz DZ --orders N [--floor log_law|one]
                           [--stress-profile] [--kappa K] [--cmu0 C] -- COMMAND [ARGS ...]

Runs COMMAND (the application on Tests/Canonical_RANS/Neutral_Column/inputs_steady_column) and checks:
  * the driver reached the steady verdict ("gust: steady after" in its output);
  * the residual of u fell by at least N orders of magnitude from its reference
    (the <checkpoint>/Driver history: the largest of the first reference
    iterations over the last record);
  * the mean friction velocity of the turbulence plotfile's surf_hist.dat is
    sqrt(F H) to 1e-4 relative: the wall stress balances the force exactly at
    steady (the pinned value; no fitted constant; sqrt(F H) must be --ustar;
    the exact discrete balance is mean(u*^2) = F H, equal to this with identical columns);
  * the layer thickness of the plotfile is --dz;
  * the first cell's speed is the wall law's, (u*/kappa) ln((z + z0)/z0), to
    1e-6 relative, the second cell's within 1 %, the lowest eight within 2 %
    (--floor one, the run's bc.surface_denominator_floor = one: the first
    cell's speed is (u*/kappa) max(ln((z + z0)/z0), 1) to 1e-6, and the log-law
    checks above the first cell are left out: that wall law is not the log
    law's when ln((z + z0)/z0) < 1; --stress-profile: the second cell and the
    lowest eight are measured against the profile of the column's stress
    u*^2 (1 - z/H) with nu_t = kappa u*_l (z + z0), u*_l = u* sqrt(1 - z/H):
    U(z) = U(z_1) + (u*/kappa) [G(z + z0) - G(z_1 + z0)],
    G(x) = (2 r + sqrt(a) ln|(r - sqrt(a))/(r + sqrt(a))|) / sqrt(H),
    r = sqrt(a - x), a = H + z0, from the first cell centre z_1; for a
    shallow column, where 1 - z/H departs from 1 in the lowest cells)
    (the stress of the force-driven column is u*^2 (1 - z/H), not constant,
    and the discrete strain overshoots the log law's in the lowest cells,
    rt_closure_column: both move the profile below the log law by about 1 % above the
    second cell);
  * the wall k is u*^2 / Cmu0^2 to 1e-6 relative;
  * nu_t in the second cell is kappa u* (z + z0) within 6 %;
  * every column is identical (the planar standard deviation of u is below
    1e-10 of u*) and every field finite.
--kappa and --cmu0 (default 0.41 and 0.5562) are the run's turb.kappa and
turb.Cmu0: the wall law, the wall k and nu_t are checked with them, so a wall
model that ignored the closure's constants fails the wall-law checks.
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

CMU0_DEFAULT = 0.5562
KAPPA_DEFAULT = 0.41


def read_history(path):
    """The residual records of <checkpoint>/Driver: list of (iteration, residuals)."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    i = 0
    while i < len(lines) and not lines[i].startswith("records "):
        i += 1
    if i >= len(lines):
        raise SystemExit("check_steady_column: no records in " + path)
    n = int(lines[i].split()[1])
    out = []
    for ln in lines[i + 1:i + 1 + n]:
        parts = ln.split()
        it = int(parts[0])
        nres = int(parts[5])
        res = [float.fromhex(x) for x in parts[6:6 + nres]]
        out.append((it, res))
    return out


def main(argv):
    opts = {"--ustar": None, "--z0": None, "--force": None, "--height": None, "--dz": None, "--orders": None}
    consts = {"--kappa": KAPPA_DEFAULT, "--cmu0": CMU0_DEFAULT}
    floor = "log_law"
    stress_profile = False
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--floor" and i + 1 < len(argv) and argv[i + 1] in ("log_law", "one"):
            floor = argv[i + 1]
            i += 2
        elif argv[i] == "--stress-profile":
            stress_profile = True
            i += 1
        elif argv[i] in consts and i + 1 < len(argv):
            consts[argv[i]] = float(argv[i + 1])
            i += 2
        elif argv[i] in opts:
            opts[argv[i]] = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if any(v is None for v in opts.values()) or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    KAPPA = consts["--kappa"]
    CMU0 = consts["--cmu0"]
    print("check_steady_column: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(run.stdout)
    if run.returncode != 0:
        sys.stdout.write(run.stderr)
        raise SystemExit("check_steady_column: the command failed with code %d" % run.returncode)
    ok = True

    def check(name, value, target, tol, kind):
        nonlocal ok
        err = abs(value - target)
        lim = tol * (abs(target) if kind == "rel" else 1.0)
        good = err <= lim
        ok = ok and good
        print("  %-46s measured %14.8e target %14.8e tol %9.2e %s" % (name, value, target, lim, "pass" if good else "FAIL"))

    check("steady verdict reached", 1.0 if "gust: steady after" in run.stdout else 0.0, 1.0, 0.0, "abs")
    # The files this run wrote (the run directory may hold earlier output).
    wrote = [ln.split()[-1] for ln in run.stdout.splitlines() if ln.startswith("gust: wrote ")]
    turb = [n for n in wrote if "_turb" in n]
    chk = [n for n in wrote if "_chk" in n]
    if not turb or not chk:
        raise SystemExit("check_steady_column: no turbulence plotfile or checkpoint written")
    hist = read_history(os.path.join(chk[-1], "Driver"))
    ref = max(r[0] for _, r in hist[:5])
    last = hist[-1][1][0]
    drop = math.log10(ref / last)
    good = drop >= opts["--orders"]
    ok = ok and good
    print("  %-46s measured %14.8e target %14.8e tol %9.2e %s" % ("orders of residual drop (u), >=", drop, opts["--orders"], 0.0, "pass" if good else "FAIL"))
    fields = ["x_velocity", "y_velocity", "KE", "Kmv", "ustar", "Lturb", "diss"]
    z, prof, _ = erf_plotfile.planar_averages(turb[-1], fields)
    _, full = erf_plotfile.read_fields(turb[-1], ["x_velocity"])
    finite = all(math.isfinite(v) for f in fields for v in prof[f])
    check("all fields finite", 1.0 if finite else 0.0, 1.0, 0.0, "abs")
    check("layer thickness z[1] - z[0] vs dz", z[1] - z[0], opts["--dz"], 1e-10, "rel")
    us = erf_plotfile.read_surf_hist(os.path.join(turb[-1], "surf_hist.dat"))["u_star"]
    us_exact = math.sqrt(opts["--force"] * opts["--height"])
    check("sqrt(F H) vs --ustar (the deck's force)", us_exact, opts["--ustar"], 1e-12, "rel")
    check("mean u* vs sqrt(F H)", us, us_exact, 1e-4, "rel")
    z0 = opts["--z0"]
    def ulog(k):
        return us / KAPPA * math.log((z[k] + z0) / z0)
    if stress_profile:
        a = opts["--height"] + z0
        def g(x):
            r = math.sqrt(a - x)
            return (2 * r + math.sqrt(a) * math.log(abs((r - math.sqrt(a)) / (r + math.sqrt(a))))) / math.sqrt(opts["--height"])
        def uprof(k):
            return ulog(0) + us / KAPPA * (g(z[k] + z0) - g(z[0] + z0))
        prof_name = "stress profile"
    else:
        uprof = ulog
        prof_name = "log law"
    if floor == "one":
        uwall = us / KAPPA * max(math.log((z[0] + z0) / z0), 1.0)
        check("|U|(k=0) vs wall law, floor one", math.hypot(prof["x_velocity"][0], prof["y_velocity"][0]), uwall, 1e-6, "rel")
    else:
        check("|U|(k=0) vs wall law", math.hypot(prof["x_velocity"][0], prof["y_velocity"][0]), ulog(0), 1e-6, "rel")
        check("|U|(k=1) vs " + prof_name, math.hypot(prof["x_velocity"][1], prof["y_velocity"][1]), uprof(1), 1e-2, "rel")
        worst = max(abs(math.hypot(prof["x_velocity"][k], prof["y_velocity"][k]) - uprof(k)) / uprof(k) for k in range(8))
        check("max rel |U - " + prof_name + "| over k < 8", worst, 0.0, 2e-2, "abs")
    check("KE(k=0) vs u*^2/Cmu0^2", prof["KE"][0], us * us / CMU0 ** 2, 1e-6, "rel")
    check("Kmv(k=1) vs kappa u* (z + z0)", prof["Kmv"][1], KAPPA * us * (z[1] + z0), 6e-2, "rel")
    nx, ny, nz = len(full["x_velocity"]), len(full["x_velocity"][0]), len(full["x_velocity"][0][0])
    spread = 0.0
    for k in range(nz):
        vals = [full["x_velocity"][i][j][k] for i in range(nx) for j in range(ny)]
        spread = max(spread, max(vals) - min(vals))
    check("max planar spread of u / u*", spread / us, 0.0, 1e-10, "abs")
    print("check_steady_column:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
