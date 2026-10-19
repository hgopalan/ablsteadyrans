#!/usr/bin/env python3
"""Fixtures for the ablsteady_refavg CTests.

    refavg_fixtures.py snapshots SRC OUT:TIME [OUT:TIME ...]
        Copy plotfile SRC to each OUT with the time in its Header replaced,
        so a single smoke plotfile gives snapshots at distinct times. Any
        earlier OUT, and the archive directory ref_archive*, are removed
        first so the tests can be rerun.
    refavg_fixtures.py damage ARCHIVE OUT
        Copy ARCHIVE to OUT and flip one bit in its first data file.
"""
import pathlib
import shutil
import sys


def set_time(header, time):
    lines = header.read_text().splitlines(keepends=True)
    ncomp = int(lines[1])
    # version, ncomp, ncomp names, ndim, then the time
    lines[2 + ncomp + 1] = f"{time}\n"
    header.write_text("".join(lines))


def snapshots(src, pairs):
    src = pathlib.Path(src)
    for old in src.parent.glob("ref_archive*"):
        shutil.rmtree(old)
    for pair in pairs:
        out, time = pair.rsplit(":", 1)
        out = src.parent / out
        if out.exists():
            shutil.rmtree(out)
        shutil.copytree(src, out)
        set_time(out / "Header", float(time))
    return 0


def damage(archive, out):
    archive, out = pathlib.Path(archive), pathlib.Path(out)
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(archive, out)
    data = sorted((out / "average" / "Level_0").glob("Cell_D_*"))[0]
    raw = bytearray(data.read_bytes())
    raw[-1] ^= 0x01
    data.write_bytes(bytes(raw))
    print(f"flipped one bit in {data.relative_to(out)}")
    return 0


def main(argv):
    if len(argv) >= 3 and argv[0] == "snapshots":
        return snapshots(argv[1], argv[2:])
    if len(argv) == 3 and argv[0] == "damage":
        return damage(argv[1], argv[2])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
