#!/usr/bin/env python3
"""Run boundary_fluxes and check the boundary fluxes of a log-law state (RT-1.8.1).

    check_bc.py --ustar U --z0 Z0 --direction DEG --dz DZ -- COMMAND [ARGS ...]

Runs COMMAND (Tests/Supporting_Files/Tools/boundary_fluxes.cpp on the Neutral_ABL_Flat deck) and
checks the plotfile it wrote:
  * the surface layer recovers the state's friction velocity at every ground
    cell: |u* - U| <= 1e-10 U (the state is the log law at the first centre
    and the reference height is that centre);
  * the ground stress is u*^2 along the wind: tau = U^2 (cos, sin)(direction)
    to 1e-10 U^2;
  * the wall k is U^2 / Cmu0^2 (Cmu0 = 0.5562);
  * the residual of the state equals its closed form in every cell (R = C - D
    with arithmetic face nu): R_u and R_v to 1e-6 of U^2/dz, whose first cell
    is the face-averaging truncation (-9.4 % of U^2/dz) and whose top cell,
    under a zero-gradient slip top that carries no stress, loses the log law's
    flux (+U^2/dz); R_theta (the sounding's diffusion with the top gradient)
    to 1e-9 K/s; R_tke and R_w zero to 1e-12;
  * every value finite.
Prints one row per check; the exit code is the verdict. The fields are
compared with amrex_fcompare in a separate test.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

CMU0 = 0.5562


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
    print("check_bc: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False)
    if run.returncode != 0:
        raise SystemExit("check_bc: the command failed with code %d" % run.returncode)
    plt = None
    for name in sorted(os.listdir(".")):
        if name.endswith("_plt00000") and os.path.isdir(name):
            plt = name
    if plt is None:
        raise SystemExit("check_bc: no plotfile written")
    names = ["ustar", "tau_x", "tau_y", "k_wall", "R_u", "R_v", "R_w", "R_theta", "R_tke", "u", "theta"]
    _, fields = erf_plotfile.read_fields(plt, names)
    U, z0, deg, dz = opts["--ustar"], opts["--z0"], opts["--direction"], opts["--dz"]
    ok = True

    def check(name, value, target, tol, kind):
        nonlocal ok
        err = abs(value - target)
        lim = tol * (abs(target) if kind == "rel" else 1.0)
        good = err <= lim
        ok = ok and good
        print("  %-44s measured %14.8e target %14.8e tol %9.2e %s" % (name, value, target, lim, "pass" if good else "FAIL"))

    def ground(f):
        return [col[0] for plane in fields[f] for col in plane]   # every (i, j) at k = 0

    def all_values(f):
        return [x for plane in fields[f] for col in plane for x in col]

    for f in names:
        if not all(math.isfinite(x) for x in all_values(f)):
            print("  %-44s not finite FAIL" % f)
            ok = False
    check("max |ustar - U| / U over the ground", max(abs(x - U) for x in ground("ustar")) / U, 0.0, 1e-10, "abs")
    want_x, want_y = U * U * math.cos(math.radians(deg)), U * U * math.sin(math.radians(deg))
    check("max |tau_x - U^2 cos| / U^2", max(abs(x - want_x) for x in ground("tau_x")) / (U * U), 0.0, 1e-10, "abs")
    check("max |tau_y - U^2 sin| / U^2", max(abs(x - want_y) for x in ground("tau_y")) / (U * U), 0.0, 1e-10, "abs")
    kw = ground("k_wall")
    check("k_wall / (U^2 / Cmu0^2)", (sum(kw) / len(kw)) / (U * U / CMU0 ** 2), 1.0, 1e-10, "rel")
    # The residual of the fixed state in closed form (R = C - D, D the face
    # flux difference over dz, face nu the arithmetic mean of the cell values):
    #   nu_k = kappa U (z_k + z0),  S_k = (U/kappa) ln((z_k + z0)/z0),
    #   F_{k+1/2} = (nu_k + nu_{k+1})/2 (S_{k+1} - S_k)/dz,  F_{-1/2} = U^2 (the
    #   wall model), F_{N-1/2} = 0 (the slip top: an even ghost carries no stress),
    #   R_k = -(F_{k+1/2} - F_{k-1/2})/dz, times (cos, sin)(direction) for u, v.
    # The first cell is the face-averaging truncation (-9.4 % of U^2/dz), the
    # top cell loses the log law's flux (+U^2/dz). theta: the sounding's
    # diffusion with nu, zero ground flux (neutral), the top gradient G through
    # the top face (ghost theta = theta_{N-2} + 2 dz G, ghost nu = nu_{N-1});
    # k is uniform (R_tke = 0); w = 0 and theta = theta_ref (R_w = 0).
    kappa, cmu0 = 0.41, CMU0
    nz = len(fields["R_u"][0][0])
    zc = [(k + 0.5) * dz for k in range(nz)]
    nu = [kappa * U * (z + z0) for z in zc]
    S = [U / kappa * math.log((z + z0) / z0) for z in zc]
    flux = [U * U] + [0.5 * (nu[k] + nu[k + 1]) * (S[k + 1] - S[k]) / dz for k in range(nz - 1)] + [0.0]
    ru_pred = [-(flux[k + 1] - flux[k]) / dz for k in range(nz)]
    tz, tt = [0.0, 700.0, 1000.0], [300.0, 300.0, 309.0]

    def theta_of(z):
        if z <= tz[0]:
            return tt[0]
        if z >= tz[-1]:
            return tt[-1]
        m = 1
        while z >= tz[m]:
            m += 1
        w = (z - tz[m - 1]) / (tz[m] - tz[m - 1])
        return tt[m - 1] + w * (tt[m] - tt[m - 1])

    th = [theta_of(z) for z in zc]
    gtop = 0.03
    ftheta = [0.0] + [0.5 * (nu[k] + nu[k + 1]) * (th[k + 1] - th[k]) / dz for k in range(nz - 1)]
    ftheta.append(nu[nz - 1] * ((th[nz - 2] + 2 * dz * gtop) - th[nz - 1]) / dz)
    rt_pred = [-(ftheta[k + 1] - ftheta[k]) / dz for k in range(nz)]
    cosd, sind = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    scale = U * U / dz

    def worst(f, pred, fac):
        return max(abs(fields[f][i][j][k] - fac * pred[k]) for i in range(len(fields[f])) for j in range(len(fields[f][0])) for k in range(nz))

    check("max |R_u - closed form| / (U^2/dz), every cell", worst("R_u", ru_pred, cosd) / scale, 0.0, 1e-6, "abs")
    check("max |R_v - closed form| / (U^2/dz), every cell", worst("R_v", ru_pred, sind) / scale, 0.0, 1e-6, "abs")
    check("first cell R_u / (U^2 cos/dz) (truncation)", ru_pred[0] / scale, -0.09713, 2e-4, "abs")
    check("top cell R_u / (U^2 cos/dz) (lost flux)", ru_pred[nz - 1] / scale, 1.0, 1e-3, "abs")
    check("max |R_theta - closed form| (K/s)", worst("R_theta", rt_pred, 1.0), 0.0, 1e-9, "abs")
    check("max |R_tke| (m^2/s^3)", max(abs(x) for x in all_values("R_tke")), 0.0, 1e-12, "abs")
    check("max |R_w| (m/s^2)", max(abs(x) for x in all_values("R_w")), 0.0, 1e-12, "abs")
    print("check_bc:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
