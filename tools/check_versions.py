#!/usr/bin/env python3
"""Check external/VERSIONS (sub-phase 1.1).

    check_versions.py [--root DIR]

Fails if the file is malformed (a line that is not "name 40-hex-sha", a
missing or repeated entry) or if the amrex entry differs from the commit
checked out in external/amrex.
"""
import argparse
import pathlib
import re
import subprocess
import sys

REQUIRED = ("erf", "amrex", "fastwindterrain", "kynema-sgf")


def parse(text):
    """Return ({name: sha}, [errors])."""
    pins, errors = {}, []
    for n, line in enumerate(text.splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        m = re.fullmatch(r"(\S+)\s+([0-9a-f]{40})", line)
        if not m:
            errors.append(f"line {n}: expected 'name <40-hex sha>', got {line!r}")
            continue
        name, sha = m.groups()
        if name in pins:
            errors.append(f"line {n}: {name} is pinned twice")
        pins[name] = sha
    for name in REQUIRED:
        if name not in pins:
            errors.append(f"{name} is not pinned")
    return pins, errors


def check(root):
    root = pathlib.Path(root)
    pins, errors = parse((root / "external" / "VERSIONS").read_text())
    if "amrex" in pins:
        proc = subprocess.run(["git", "-C", str(root / "external" / "amrex"), "rev-parse", "HEAD"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            errors.append("cannot read the external/amrex submodule commit "
                          "(run git submodule update --init)")
        elif proc.stdout.strip() != pins["amrex"]:
            errors.append(f"amrex is pinned to {pins['amrex']} but external/amrex is at "
                          f"{proc.stdout.strip()}")
    return errors


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(pathlib.Path(__file__).resolve().parents[1]))
    args = ap.parse_args(argv)
    errors = check(args.root)
    if errors:
        print("external/VERSIONS check FAILED:")
        for e in errors:
            print("  " + e)
        return 1
    print("external/VERSIONS check passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
