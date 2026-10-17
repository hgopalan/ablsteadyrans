#!/usr/bin/env python3
"""Run poisson_projection and check it against a gold file (RT-1.10.1).

    check_poisson.py GOLD.json [--write] -- COMMAND [ARGS ...]

Runs COMMAND and checks what it prints:
  * the GMRES iterations, at most 10 % above the gold count;
  * the divergence after the projection, at most 10 x solver.poisson_reltol
    (as the run prints it) times the divergence before;
  * the divergence before, equal to the gold value to a relative 1e-9 (the
    projected wind and the mesh are fixed, so it only changes if they do).
Prints one row per check; the exit code is the verdict. --write stores the
measured values as the new gold (review the diff before committing it).
The potential itself is compared with amrex_fcompare in a separate test.
"""
import json
import re
import subprocess
import sys


def measure(out):
    rt = re.search(r"poisson_projection: reltol (\S+)", out)
    it = re.search(r"poisson_projection: iterations (\d+)", out)
    dv = re.search(r"poisson_projection: divergence before (\S+) after (\S+)", out)
    if not rt or not it or not dv:
        raise SystemExit("check_poisson: the run did not print its tolerance, iterations and divergence")
    return {"reltol": float(rt.group(1)), "iterations": int(it.group(1)), "div_before": float(dv.group(1)),
            "div_after": float(dv.group(2))}


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    write = "--write" in opts
    gold = [o for o in opts if o != "--write"]
    if len(gold) != 1 or not cmd:
        print(__doc__)
        return 2
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    print(proc.stdout)
    if proc.returncode != 0:
        print(f"check_poisson: the run exited {proc.returncode}")
        return 1
    m = measure(proc.stdout)
    if write:
        with open(gold[0], "w") as fh:
            json.dump(m, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"check_poisson: wrote {gold[0]}")
        return 0
    with open(gold[0]) as fh:
        g = json.load(fh)
    rows = [
        ("iterations", m["iterations"], f"<= {1.1 * g['iterations']:.1f}", m["iterations"] <= 1.1 * g["iterations"]),
        ("divergence after / before", m["div_after"] / m["div_before"], f"<= {10 * m['reltol']:g}",
         m["div_after"] <= 10 * m["reltol"] * m["div_before"]),
        ("divergence before", m["div_before"], f"{g['div_before']:.9e} (rel 1e-9)",
         abs(m["div_before"] - g["div_before"]) <= 1e-9 * abs(g["div_before"])),
    ]
    ok = True
    print(f"{'check':30s} {'measured':>18} {'target':>30} pass")
    for name, meas, target, passed in rows:
        print(f"{name:30s} {meas!s:>18} {target:>30} {'pass' if passed else 'FAIL'}")
        ok &= passed
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
