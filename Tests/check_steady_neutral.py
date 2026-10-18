#!/usr/bin/env python3
"""Run ablsteady on the Neutral_ABL_Flat deck and apply ERF's checks (RT-1.15.2).

    check_steady_neutral.py -- COMMAND [ARGS ...]

Runs COMMAND (the application on Tests/Decks/inputs_steady_neutral_flat),
requires the steady verdict, then runs ERF's own check_neutral.py --physics
(extern_erf/check_neutral.py, copied from
Exec/CanonicalTests/Canonical_RANS/Neutral_ABL_Flat at the pinned ERF commit)
on the turbulence plotfile the stage wrote (io.turb_plot_file: ERF's variable
names and surf_hist.dat) and reports its verdict: every field finite, k >= 0,
the capped kappa (z + z0) length, the AL01 stability functions of Rt, Kmv =
cmu sqrt(k) L, the dissipation Cmu0^3 k^1.5 / L within 5 %, the wall cell's
k balance, and the --physics set: u* within 0.35 +- 0.15 m/s, the wall k
u*^2/Cmu0^2 within 5 %, the log law in the lowest four cells within 10 %,
Kmv(k=1) within 30 % of kappa u* z. The exit code is the verdict of both.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ERF_CHECK = os.path.join(HERE, "..", "extern_erf", "check_neutral.py")


def main(argv):
    if len(argv) < 3 or argv[1] != "--":
        raise SystemExit(__doc__)
    cmd = argv[2:]
    print("check_steady_neutral: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(run.stdout)
    if run.returncode != 0:
        sys.stdout.write(run.stderr)
        raise SystemExit("check_steady_neutral: the command failed with code %d" % run.returncode)
    ok = "ablsteady: steady after" in run.stdout
    print("  steady verdict reached: %s" % ("pass" if ok else "FAIL"))
    wrote = [ln.split()[-1] for ln in run.stdout.splitlines() if ln.startswith("ablsteady: wrote ")]
    turb = [n for n in wrote if "_turb" in n]   # the files this run wrote, not earlier output
    if not turb:
        raise SystemExit("check_steady_neutral: no turbulence plotfile written")
    plt = turb[-1]
    erf = subprocess.run([sys.executable, "-I", ERF_CHECK, "--physics", plt, os.path.join(plt, "surf_hist.dat")],
                         check=False, capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    sys.stdout.write(erf.stdout)
    sys.stdout.write(erf.stderr)
    print("  ERF check_neutral.py --physics: %s" % ("pass" if erf.returncode == 0 else "FAIL"))
    ok = ok and erf.returncode == 0
    print("check_steady_neutral:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
