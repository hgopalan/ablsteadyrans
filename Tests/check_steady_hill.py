#!/usr/bin/env python3
"""Run ablsteady on a hill deck and apply ERF's own checks (RT-1.16.1, RT-1.16.2).

    check_steady_hill.py --erf-check check_hill.py|check_hill3d.py [--orders N] -- COMMAND [ARGS ...]

Runs COMMAND (the application on Tests/Decks/inputs_steady_hill2d or
inputs_steady_hill3d), requires the steady verdict (and, with --orders, a
drop of N orders of the residual of u from its reference, from the Driver
history of the checkpoint), then runs ERF's check script from extern_erf/
with --physics on the turbulence plotfile the stage wrote (io.turb_plot_file:
ERF's variable names): the wall distance against the exact distance to the
hill, the length-scale bound, positivity, the wall cell's k balance, the
crest speed-up against Jackson and Hunt's 2 h/L (1.6 h/L for the hill of
revolution) within a factor of two, the upstream u* within 0.25-0.55 m/s and
the upstream log law in the lowest three cells within 15 %. The exit code is
the verdict of both.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXTERN = os.path.join(HERE, "..", "extern_erf")


def read_history(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    i = 0
    while i < len(lines) and not lines[i].startswith("records "):
        i += 1
    n = int(lines[i].split()[1])
    out = []
    for ln in lines[i + 1:i + 1 + n]:
        parts = ln.split()
        nres = int(parts[5])
        out.append([float.fromhex(x) for x in parts[6:6 + nres]])
    return out


def main(argv):
    script, orders = None, None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--erf-check":
            script = argv[i + 1]
            i += 2
        elif argv[i] == "--orders":
            orders = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if script is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    print("check_steady_hill: running", " ".join(cmd))
    run = subprocess.run(cmd, check=False, capture_output=True, text=True)
    sys.stdout.write(run.stdout)
    if run.returncode != 0:
        sys.stdout.write(run.stderr)
        raise SystemExit("check_steady_hill: the command failed with code %d" % run.returncode)
    ok = "ablsteady: steady after" in run.stdout
    print("  steady verdict reached: %s" % ("pass" if ok else "FAIL"))
    wrote = [ln.split()[-1] for ln in run.stdout.splitlines() if ln.startswith("ablsteady: wrote ")]
    turb = [n for n in wrote if "_turb" in n]
    chk = [n for n in wrote if "_chk" in n]
    if not turb:
        raise SystemExit("check_steady_hill: no turbulence plotfile written")
    if orders is not None:
        if not chk:
            raise SystemExit("check_steady_hill: no checkpoint written (the residual drop needs its Driver history)")
        hist = read_history(os.path.join(chk[-1], "Driver"))
        ref = max(r[0] for r in hist[:100])
        drop = math.log10(ref / hist[-1][0])
        good = drop >= orders
        ok = ok and good
        print("  residual of u dropped %.2f orders from its reference (>= %g): %s" % (drop, orders, "pass" if good else "FAIL"))
    erf = subprocess.run([sys.executable, "-I", os.path.join(EXTERN, script), "--physics", turb[-1]],
                         check=False, capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    sys.stdout.write(erf.stdout)
    sys.stdout.write(erf.stderr)
    print("  ERF %s --physics: %s" % (script, "pass" if erf.returncode == 0 else "FAIL"))
    ok = ok and erf.returncode == 0
    print("check_steady_hill:", "all checks passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
