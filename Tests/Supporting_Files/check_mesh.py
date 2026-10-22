#!/usr/bin/env python3
"""Run GUST and check its mesh against a gold file (rt_askervein_mesh,
rt_steep_ridge_quality and the grid-hierarchy tests).

    check_mesh.py GOLD.json [--write] [--rtol R] -- COMMAND [ARGS ...]

Runs COMMAND, then compares
  * the mesh-quality report it prints (every "name = value" line, the STF
    capped count and whether the steep-terrain and top-displacement
    warnings appear), and
  * statistics of the plotfile it wrote (sum, minimum, maximum and sum of
    squares of the cell fields z and detJ)
with GOLD.json. With a grid hierarchy (mesh.hierarchy_levels > 1) the run
prints one more report per coarse level ("gust: mesh quality, hierarchy
level l ..."); each is checked the same way against the gold's "levels"
entry, and the number of levels must match. Prints one row per check: measured, target, tolerance, pass;
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
sys.path.insert(0, os.path.join(HERE, "..", "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

FIELDS = ("z", "detJ")


LEVEL_HEADER = re.compile(r"^gust: mesh quality, hierarchy level (\d+)", re.M)
HEADER_FIELDS = re.compile(r"hierarchy level \d+ \((\d+) x (\d+) x (\d+) cells, dx (\S+) m, dy (\S+) m, "
                           r"ground smoothing passes (\d+), boxes (\d+)(, gathered)?\)")
SUMMARY_ROW = re.compile(r"^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)$")
SUMMARY_KEYS = ("level", "nx", "ny", "dx", "dy", "passes", "max_aspect", "min_thick", "max_slope", "capped")


def measure_report(text):
    """The report lines, STF count and warnings of one report block."""
    m = {"report": {}}
    for line in text.splitlines():
        r = re.match(r"\s{2}(.+?)\s+=\s+([-+0-9.eE]+)$", line)
        if r:
            m["report"][r.group(1)] = float(r.group(2))
        r = re.match(r"\s{2}STF levels capped by terrain\.gamma_m\s+=\s+(\d+) of (\d+)", line)
        if r:
            m["stf_capped"] = int(r.group(1))
    m["warning"] = "WARNING: the STF smoothing stopped early" in text
    m["top_warning"] = "WARNING: the terrain imprint reaches the top" in text
    return m


def measure(out):
    # Level 0's report runs up to the first coarse level's header (or the
    # end); each coarse report up to the next header or the summary.
    heads = list(LEVEL_HEADER.finditer(out))
    end0 = heads[0].start() if heads else len(out)
    m = measure_report(out[:end0])
    m["stats"] = {}
    if heads:
        m["levels"] = []
        for n, h in enumerate(heads):
            stop = heads[n + 1].start() if n + 1 < len(heads) else out.find("gust: grid hierarchy summary", h.start())
            block = out[h.start():stop if stop >= 0 else len(out)]
            lev = measure_report(block)
            lev["level"] = int(h.group(1))
            lev["hint"] = f"WARNING: hierarchy level {lev['level']}: the STF smoothing stopped early" in block
            f = HEADER_FIELDS.search(block)
            lev["header"] = None if f is None else {
                "nx": int(f.group(1)), "ny": int(f.group(2)), "nz": int(f.group(3)), "dx": float(f.group(4)),
                "dy": float(f.group(5)), "passes": int(f.group(6)), "boxes": int(f.group(7)),
                "gathered": f.group(8) is not None}
            m["levels"].append(lev)
        m["summary"] = []
        start = out.find("gust: grid hierarchy summary")
        for line in (out[start:].splitlines()[2:] if start >= 0 else []):
            r = SUMMARY_ROW.match(line)
            if not r:
                break
            m["summary"].append({k: (float(v) if "." in v or "e" in v else int(v))
                                 for k, v in zip(SUMMARY_KEYS, r.groups())})
    plt = re.findall(r"gust: wrote (\S+)", out)
    # Level 0's plotfile: not a checkpoint, a reference archive or a hierarchy level's.
    plt = [p for p in plt if not os.path.exists(os.path.join(p, "MANIFEST")) and "chk" not in p
           and "_level" not in p]
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

    def compare_report(prefix, mr, gr):
        for name, b in gr["report"].items():
            a = mr["report"].get(name)
            # Absolute floor 1e-9: a value at rounding level (a flat top's
            # displacement, about 1e-13 m) differs between compilers.
            row(prefix + name, a, b, f"{report_tol:g}", a is not None and close(a, b, report_tol, 1e-9))
        row(prefix + "STF levels capped", mr.get("stf_capped"), gr["stf_capped"], "exact",
            mr.get("stf_capped") == gr["stf_capped"])
        row(prefix + "steep-terrain warning", mr["warning"], gr["warning"], "exact", mr["warning"] == gr["warning"])
        row(prefix + "top-displacement warning", mr["top_warning"], gr["top_warning"], "exact",
            mr["top_warning"] == gr["top_warning"])

    compare_report("", m, g)
    ml, gl = m.get("levels", []), g.get("levels", [])
    row("hierarchy levels below level 0", len(ml), len(gl), "exact", len(ml) == len(gl))
    for a, b in zip(ml, gl):
        p = f"level {b['level']}: "
        row(p + "level", a["level"], b["level"], "exact", a["level"] == b["level"])
        compare_report(p, a, b)
        row(p + "hierarchy hint", a["hint"], b["hint"], "exact", a["hint"] == b["hint"])
        for k, v in (b.get("header") or {}).items():
            got = (a.get("header") or {}).get(k)
            ok_k = got is not None and (close(got, v, report_tol) if isinstance(v, float) else got == v)
            row(p + "header " + k, got, v, "exact" if not isinstance(v, float) else f"{report_tol:g}", ok_k)
    ms, gs = m.get("summary", []), g.get("summary", [])
    row("hierarchy summary rows", len(ms), len(gs), "exact", len(ms) == len(gs))
    # The summary prints 4 to 5 significant digits: reals to 1e-4.
    for a, b in zip(ms, gs):
        for k, v in b.items():
            got = a.get(k)
            ok_k = got is not None and (close(got, v, 1e-4) if isinstance(v, float) else got == v)
            row(f"summary level {b['level']} {k}", got, v, "1e-4" if isinstance(v, float) else "exact", ok_k)
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
