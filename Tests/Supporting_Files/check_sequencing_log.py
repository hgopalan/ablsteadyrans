#!/usr/bin/env python3
"""RT-2.6.2: the stage switches and the grid-sequencing level history of a run log against a gold file.

    check_sequencing_log.py LOG GOLD.json [--iter_rel_tol 0.1] [--write]

LOG is the standard output of an ablsteady run with steady.sequencing = 1.
GOLD.json holds the expected history:

    {"levels": [{"level": 1, "cells": [150, 150], "stop": "every residual ...", "iterations": 100}, ...],
     "switches": [],
     "verdict": "steady"}

Passes (exit code 0) when the levels appear in the gold order with the same
cells and the same stop reason (exact text), every level's iterations within
--iter_rel_tol (relative) of the gold count, the stage switches are the gold list, in order, and the run ends with
the gold verdict. A switch is recorded as "<from> to <to>: <reason>" from the
driver's "ablsteady: iteration N: switching from stage ..." line (the
iteration and the residuals are left out: they move with rounding). Iteration counts move with the compiler's rounding
(gcc contracts multiply-adds), so they get a tolerance; the sequence does not.
With --write, writes GOLD.json from LOG instead (for a reviewed gold change).
"""
import argparse
import json
import re
import sys


def parse(log):
    levels, switches, verdict, inside = [], [], None, False
    for line in log.splitlines():
        if line.startswith("ablsteady: grid sequencing summary"):
            inside, levels = True, []
            continue
        if inside:
            m = re.match(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(.*\S)\s*$", line)
            if m:
                levels.append({"level": int(m.group(1)), "cells": [int(m.group(2)), int(m.group(3))],
                               "iterations": int(m.group(4)), "stop": m.group(5)})
                continue
            if line.startswith("  level") or line.startswith("         prolonged"):
                continue
            inside = False
        m = re.match(r"ablsteady: (?:level \d+: )?iteration \d+: switching from stage (\S+) to (\S+): (.*?)"
                     r"(?: \(scaled residuals:.*)?$", line)
        if m:
            switches.append(f"{m.group(1)} to {m.group(2)}: {m.group(3)}")
        m = re.match(r"ablsteady: (steady|not steady) after \d+ iterations", line)
        if m:
            verdict = m.group(1)
    return {"levels": levels, "switches": switches, "verdict": verdict}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("gold")
    ap.add_argument("--iter_rel_tol", type=float, default=0.1)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    if not (a.iter_rel_tol >= 0):
        sys.exit("--iter_rel_tol must be >= 0")
    with open(a.log) as fh:
        got = parse(fh.read())
    if not got["levels"]:
        print(f"FAIL: no grid sequencing summary in {a.log}")
        return 1
    if a.write:
        with open(a.gold, "w") as fh:
            json.dump(got, fh, indent=1)
            fh.write("\n")
        print(f"wrote {a.gold}")
        return 0
    with open(a.gold) as fh:
        gold = json.load(fh)
    fails = []
    if [lv["level"] for lv in got["levels"]] != [lv["level"] for lv in gold["levels"]]:
        fails.append(f"levels {[lv['level'] for lv in got['levels']]}, gold {[lv['level'] for lv in gold['levels']]}")
    for g, e in zip(got["levels"], gold["levels"]):
        name = f"level {e['level']}"
        if g["cells"] != e["cells"]:
            fails.append(f"{name}: cells {g['cells']}, gold {e['cells']}")
        if g["stop"] != e["stop"]:
            fails.append(f"{name}: stop '{g['stop']}', gold '{e['stop']}'")
        if abs(g["iterations"] - e["iterations"]) > a.iter_rel_tol * e["iterations"]:
            fails.append(f"{name}: {g['iterations']} iterations, gold {e['iterations']} (tolerance {a.iter_rel_tol:g})")
        print(f"  {name}: {g['cells'][0]} x {g['cells'][1]}, {g['iterations']} iterations (gold {e['iterations']}), "
              f"stop: {g['stop']}")
    if got["switches"] != gold["switches"]:
        fails.append(f"switches {got['switches']}, gold {gold['switches']}")
    if got["verdict"] != gold["verdict"]:
        fails.append(f"verdict {got['verdict']}, gold {gold['verdict']}")
    for f in fails:
        print("FAIL: " + f)
    if not fails:
        print(f"PASS: {len(got['levels'])} levels, {len(got['switches'])} switches, verdict {got['verdict']}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
