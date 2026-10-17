#!/usr/bin/env python3
"""Run ablsteady and check its mesh against a gold file (RT-1.4.1, RT-1.4.2).

    check_mesh.py GOLD.json [--write] [--rtol R] -- COMMAND [ARGS ...]

Runs COMMAND, then compares
  * the mesh-quality report it prints (every "name = value" line, the STF
    capped count and whether the steep-terrain and top-displacement
    warnings appear), and
  * statistics of the plotfile it wrote (sum, minimum, maximum and sum of
    squares of the cell fields z and detJ)
with GOLD.json. Prints one row per check: measured, target, tolerance, pass;
the exit code is the verdict. Real values are compared to a relative
tolerance of 1e-12 (statistics) or 1e-9 (the report's printed values, with an
absolute floor of 1e-9 for values at rounding level), which
absorbs last-bit differences between compilers that do and do not fuse
multiply-adds; counts and the warnings are compared exactly.

--rtol R replaces both tolerances by R (a single-precision build uses 1e-5:
the golds are double-precision values). --write stores the measured values as
the new gold instead (review the diff before committing it).
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

FIELDS = ("z", "detJ")


def measure(out):
    m = {"report": {}, "stats": {}}
    for line in out.splitlines():
        r = re.match(r"\s{2}(.+?)\s+=\s+([-+0-9.eE]+)$", line)
        if r:
            m["report"][r.group(1)] = float(r.group(2))
        r = re.match(r"\s{2}STF levels capped by terrain\.gamma_m\s+=\s+(\d+) of (\d+)", line)
        if r:
            m["stf_capped"] = int(r.group(1))
    m["warning"] = "WARNING: the STF smoothing stopped early" in out
    m["top_warning"] = "WARNING: the terrain imprint reaches the top" in out
    plt = re.findall(r"ablsteady: wrote (\S+)", out)
    plt = [p for p in plt if not os.path.exists(os.path.join(p, "MANIFEST")) and "chk" not in p]
    if not plt:
        raise SystemExit("check_mesh: the run wrote no plotfile")
    _, data = erf_plotfile.read_fields(plt[0], list(FIELDS))
    for f in FIELDS:
        vals = [v for plane in data[f] for col in plane for v in col]
        m["stats"][f] = {"sum": sum(vals), "min": min(vals), "max": max(vals),
                         "sum_sq": sum(v * v for v in vals)}
    return m


def close(a, b, rtol, atol=0.0):
    return abs(a - b) <= max(rtol * abs(b), atol) or a == b


def compare(m, g, rtol=None):
    rows, ok = [], True
    stats_tol = rtol if rtol is not None else 1e-12
    report_tol = rtol if rtol is not None else 1e-9

    def row(name, meas, target, tol, passed):
        nonlocal ok
        rows.append(f"{name:58s} {meas!s:>24} {target!s:>24} {tol!s:>8} {'pass' if passed else 'FAIL'}")
        ok &= passed

    for f in FIELDS:
        for k in ("sum", "min", "max", "sum_sq"):
            a, b = m["stats"][f][k], g["stats"][f][k]
            row(f"{f}.{k}", f"{a:.15e}", f"{b:.15e}", f"{stats_tol:g}", close(a, b, stats_tol))
    for name, b in g["report"].items():
        a = m["report"].get(name)
        # Absolute floor 1e-9: a value at rounding level (a flat top's
        # displacement, about 1e-13 m) differs between compilers.
        row(name, a, b, f"{report_tol:g}", a is not None and close(a, b, report_tol, 1e-9))
    row("STF levels capped", m.get("stf_capped"), g["stf_capped"], "exact", m.get("stf_capped") == g["stf_capped"])
    row("steep-terrain warning", m["warning"], g["warning"], "exact", m["warning"] == g["warning"])
    row("top-displacement warning", m["top_warning"], g["top_warning"], "exact",
        m["top_warning"] == g["top_warning"])
    print(f"{'check':58s} {'measured':>24} {'target':>24} {'tol':>8} pass")
    print("\n".join(rows))
    return ok


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    write = "--write" in opts
    rtol = None
    if "--rtol" in opts:
        i = opts.index("--rtol")
        rtol = float(opts[i + 1])
        opts = opts[:i] + opts[i + 2:]
    gold = [o for o in opts if o != "--write"]
    if len(gold) != 1 or not cmd:
        print(__doc__)
        return 2
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    print(proc.stdout)
    if proc.returncode != 0:
        print(f"check_mesh: the run exited {proc.returncode}")
        return 1
    m = measure(proc.stdout)
    if write:
        with open(gold[0], "w") as fh:
            json.dump(m, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"check_mesh: wrote {gold[0]}")
        return 0
    with open(gold[0]) as fh:
        g = json.load(fh)
    return 0 if compare(m, g, rtol) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
