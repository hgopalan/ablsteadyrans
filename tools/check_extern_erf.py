#!/usr/bin/env python3
"""Drift check for code copied from ERF.

    check_extern_erf.py --erf PATH [--root DIR] [--report-upstream REF]

Reads the manifest table in <root>/extern_erf/EXTRACTED.md and checks, against
the ERF git clone at PATH:

  * an `unchanged` file is byte-identical to its ERF source at its ERF SHA;
  * a `modified` file has <root>/extern_erf/patches/<name>.patch equal to the
    unified diff from the ERF source to the local file;
  * every file under extern_erf/ (except EXTRACTED.md, LICENSE.md and
    patches/) has a manifest row, and every row names a file that exists.

Exit code 1 if any check fails. With --report-upstream REF it also prints the
upstream changes to each ERF source between its pinned SHA and REF; that
report never changes the exit code.

With an empty manifest, the checks pass without reading ERF at all.
"""
import argparse
import difflib
import pathlib
import subprocess
import sys

EXEMPT = {"EXTRACTED.md", "LICENSE.md"}
STATUSES = {"unchanged", "modified"}


def parse_manifest(text):
    """Rows between the manifest markers, as dicts. Raises ValueError."""
    try:
        body = text.split("<!-- manifest:start -->", 1)[1].split("<!-- manifest:end -->", 1)[0]
    except IndexError:
        raise ValueError("EXTRACTED.md has no <!-- manifest:start/end --> markers")
    rows = []
    lines = [ln.strip() for ln in body.strip().splitlines() if ln.strip()]
    if len(lines) < 2 or not lines[0].startswith("| Local file") or not lines[1].startswith("|---"):
        raise ValueError("manifest table must start with its header row and separator")
    for ln in lines[2:]:  # header and separator
        cells = [c.strip() for c in ln.strip("|").split("|")]
        if len(cells) < 4:
            raise ValueError(f"manifest row has fewer than 4 columns: {ln}")
        local, source, sha, status = (c.strip("`") for c in cells[:4])
        if status not in STATUSES:
            raise ValueError(f"{local}: status must be one of {sorted(STATUSES)}, got {status!r}")
        rows.append({"local": local, "source": source, "sha": sha, "status": status})
    return rows


def patch_path(root, local):
    return root / "extern_erf" / "patches" / (pathlib.PurePosixPath(local).name + ".patch")


def make_patch(upstream, local_text, source, local):
    return "".join(difflib.unified_diff(
        upstream.splitlines(keepends=True), local_text.splitlines(keepends=True),
        fromfile=f"a/{source}", tofile=f"b/{local}"))


def erf_show(erf, sha, source):
    proc = subprocess.run(["git", "-C", str(erf), "show", f"{sha}:{source}"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        return None
    return proc.stdout


def check(root, erf, report_ref=None, out=print):
    root = pathlib.Path(root)
    errors = []
    try:
        rows = parse_manifest((root / "extern_erf" / "EXTRACTED.md").read_text())
    except (OSError, ValueError) as e:
        return [str(e)]

    listed = {r["local"] for r in rows}
    ext = root / "extern_erf"
    for f in sorted(p for p in ext.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"):
        rel = f.relative_to(root).as_posix()
        if f.name in EXEMPT and f.parent == ext:
            continue
        if (ext / "patches") in f.parents:
            continue
        if rel not in listed:
            errors.append(f"{rel}: in extern_erf/ but has no manifest row")

    for r in rows:
        local = root / r["local"]
        if not local.is_file():
            errors.append(f"{r['local']}: manifest row names a file that does not exist")
            continue
        upstream = erf_show(erf, r["sha"], r["source"]) if erf else None
        if upstream is None:
            errors.append(f"{r['local']}: cannot read {r['source']} at ERF {r['sha']}")
            continue
        mine = local.read_bytes()
        if r["status"] == "unchanged":
            if mine != upstream:
                errors.append(f"{r['local']}: differs from ERF {r['source']} at {r['sha']} "
                              "but is recorded as unchanged")
        else:
            pp = patch_path(root, r["local"])
            want = make_patch(upstream.decode(), mine.decode(), r["source"], r["local"])
            if not pp.is_file():
                errors.append(f"{r['local']}: recorded as modified but {pp.relative_to(root)} is missing")
            elif pp.read_text() != want:
                errors.append(f"{r['local']}: {pp.relative_to(root)} does not match the current "
                              "local changes; regenerate it and review the diff")
            elif not want:
                errors.append(f"{r['local']}: recorded as modified but identical to ERF")

    if report_ref and erf and rows:
        out(f"Upstream ERF changes since each pinned SHA, up to {report_ref} (report only):")
        for r in rows:
            proc = subprocess.run(["git", "-C", str(erf), "log", "--oneline",
                                   f"{r['sha']}..{report_ref}", "--", r["source"]],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            log = proc.stdout.strip() if proc.returncode == 0 else f"(git log failed: {proc.stderr.strip()})"
            out(f"  {r['source']}: " + ("no changes" if not log else "\n    " + log.replace("\n", "\n    ")))
    return errors


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--erf", help="path to an ERF git clone (needed once the manifest has rows)")
    ap.add_argument("--root", default=str(pathlib.Path(__file__).resolve().parents[1]),
                    help="GUST source tree (default: this script's repository)")
    ap.add_argument("--report-upstream", metavar="REF",
                    help="also list ERF commits touching each source since its pin, up to REF")
    args = ap.parse_args(argv)

    errors = check(args.root, args.erf, args.report_upstream)
    if errors:
        print("extern_erf drift check FAILED:")
        for e in errors:
            print("  " + e)
        return 1
    print("extern_erf drift check passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
