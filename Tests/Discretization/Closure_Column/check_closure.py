#!/usr/bin/env python3
"""Run closure_column and check the full residual of the equilibrium log-law state (rt_closure_column).

    check_closure.py --ustar U --z0 Z0 --direction DEG --dz DZ -- COMMAND [ARGS ...]

Runs COMMAND (Tests/Supporting_Files/Tools/closure_column.cpp on the neutral column deck) and
checks the plotfile it wrote against the closed form of the wall law's
equilibrium state, u = (U/kappa) ln((z + z0)/z0), k = U^2/Cmu0^2, with the
uncapped geometric length and the logarithmic face mean:
  * nu_t = Cmu0 sqrt(k) kappa (z + z0) = kappa U (z + z0) at every cell (1e-9);
  * the surface layer recovers U at every ground cell (1e-10);
  * the momentum flux of the logarithmic profile through every interior face
    is U^2, so R_u and R_v vanish below the top cell (1e-10 of U^2/dz) and the
    top cell, under a slip top that carries no stress, loses U^2/dz (times
    cos, sin of the direction; 1e-10);
  * R_tke = -(P_s - eps) with P_s = nu_t S^2 from the discrete strain, S^2 =
    ((S_{k+1} - S_{k-1})/(2 dz))^2 (the ground ghost from the wall law's
    gradient, the top ghost even), eps = Cmu0^3 k^1.5 / (kappa (z + z0)) =
    U^3/(kappa (z + z0)); every cell above the first to 1e-9 of eps; the first
    cell is the wall-k row, 0 exactly since k = k_w; the discrete strain
    overshoots the log law's in the lowest cells (P_s/eps = 1.45 in the second
    cell), which is the truncation the k equation sees, pinned as a closed form;
  * R_theta, R_w zero (uniform theta, w = 0; 1e-12); every value finite.
Prints one row per check; the exit code is the verdict. The fields are
compared with amrex_fcompare in a separate test.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

CMU0 = 0.5562
KAPPA = 0.41


def main(argv):
    opts = {"--ustar": None, "--z0": None, "--direction": None, "--dz": None}
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] in opts:
            opts[argv[i]] = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if any(v is None for v in opts.values()) or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    print("check_closure: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False)
    if run.returncode != 0:
        raise SystemExit("check_closure: the command failed with code %d" % run.returncode)
    plt = None
    for name in sorted(os.listdir(".")):
        if name.endswith("_plt00000") and os.path.isdir(name):
            plt = name
    if plt is None:
        raise SystemExit("check_closure: no plotfile written")
    names = ["u", "tke", "nut", "length", "R_u", "R_v", "R_w", "R_theta", "R_tke", "Ps", "eps", "ustar"]
    _, fields = erf_plotfile.read_fields(plt, names)
    U, z0, deg, dz = opts["--ustar"], opts["--z0"], opts["--direction"], opts["--dz"]
    ok = True

    def check(name, value, target, tol, kind):
        nonlocal ok
        err = abs(value - target)
        lim = tol * (abs(target) if kind == "rel" else 1.0)
        good = err <= lim
        ok = ok and good
        print("  %-46s measured %14.8e target %14.8e tol %9.2e %s" % (name, value, target, lim, "pass" if good else "FAIL"))

    def all_values(f):
        return [x for plane in fields[f] for col in plane for x in col]

    for f in names:
        if not all(math.isfinite(x) for x in all_values(f)):
            print("  %-46s not finite FAIL" % f)
            ok = False
    nx, ny, nz = len(fields["u"]), len(fields["u"][0]), len(fields["u"][0][0])
    zc = [(k + 0.5) * dz for k in range(nz)]
    S = [U / KAPPA * math.log((z + z0) / z0) for z in zc]
    nut = [KAPPA * U * (z + z0) for z in zc]
    k_eq = U * U / CMU0 ** 2
    # Discrete strain: the ground ghost gives the wall law's gradient U/(kappa (z_0 + z0))
    # at the first centre; the top ghost is even.
    dsdz = []
    for k in range(nz):
        if k == 0:
            dsdz.append(U / (KAPPA * (zc[0] + z0)))
        elif k == nz - 1:
            dsdz.append((S[nz - 1] - S[nz - 2]) / (2 * dz))
        else:
            dsdz.append((S[k + 1] - S[k - 1]) / (2 * dz))
    Ps = [nut[k] * dsdz[k] ** 2 for k in range(nz)]
    eps = [U ** 3 / (KAPPA * (z + z0)) for z in zc]
    rk_pred = [-(Ps[k] - eps[k]) for k in range(nz)]
    rk_pred[0] = 0.0   # the wall-k Dirichlet row: k = k_w exactly
    cosd, sind = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    scale = U * U / dz

    def worst(f, pred, fac, k0=0, k1=None):
        k1 = nz if k1 is None else k1
        return max(abs(fields[f][i][j][k] - fac * pred[k]) for i in range(nx) for j in range(ny) for k in range(k0, k1))

    check("max |ustar - U| / U over the ground", max(abs(fields["ustar"][i][j][0] - U) for i in range(nx) for j in range(ny)) / U, 0.0, 1e-10, "abs")
    check("max |tke - U^2/Cmu0^2| / k", worst("tke", [k_eq] * nz, 1.0) / k_eq, 0.0, 1e-12, "abs")
    check("max |nu_t - kappa U (z + z0)| / nu_t(top)", worst("nut", nut, 1.0) / nut[-1], 0.0, 1e-9, "abs")
    check("max |length - kappa (z + z0)| / L(top)", worst("length", [KAPPA * (z + z0) for z in zc], 1.0) / (KAPPA * (zc[-1] + z0)), 0.0, 1e-9, "abs")
    check("max |R_u| below the top cell / (U^2/dz)", worst("R_u", [0.0] * nz, 1.0, 0, nz - 1) / scale, 0.0, 1e-10, "abs")
    check("max |R_v| below the top cell / (U^2/dz)", worst("R_v", [0.0] * nz, 1.0, 0, nz - 1) / scale, 0.0, 1e-10, "abs")
    # 1e-8: the geometric length's cap at 1e12 m is 2e-11 off, the log mean 1e-14.
    check("top cell R_u / (U^2 cos/dz)", max(fields["R_u"][i][j][nz - 1] for i in range(nx) for j in range(ny)) / (scale * cosd), 1.0, 1e-8, "rel")
    check("top cell R_v / (U^2 sin/dz)", max(fields["R_v"][i][j][nz - 1] for i in range(nx) for j in range(ny)) / (scale * sind), 1.0, 1e-8, "rel")
    check("max |P_s - closed form| / eps(1)", worst("Ps", Ps, 1.0, 1) / eps[1], 0.0, 1e-9, "abs")
    check("max |eps - closed form| / eps(1)", worst("eps", eps, 1.0) / eps[1], 0.0, 1e-9, "abs")
    check("max |R_tke - closed form| / eps(1)", worst("R_tke", rk_pred, 1.0) / eps[1], 0.0, 1e-9, "abs")
    # The discrete strain overshoots the log law's in the lowest cells (the
    # centred difference of ln over 2 dz): the ratio is a closed form, 1.45 here.
    check("P_s / eps in the second cell (strain truncation)", fields["Ps"][0][0][1] / fields["eps"][0][0][1], Ps[1] / eps[1], 1e-9, "rel")
    check("max |R_theta| (K/s)", max(abs(x) for x in all_values("R_theta")), 0.0, 1e-12, "abs")
    check("max |R_w| (m/s^2)", max(abs(x) for x in all_values("R_w")), 0.0, 1e-12, "abs")
    print("check_closure:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
