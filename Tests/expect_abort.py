#!/usr/bin/env python3
"""Run a command and check both its exit code and its output.

    expect_abort.py [--success] REGEX [REGEX ...] -- COMMAND [ARGS ...]

By default passes only if COMMAND exits non-zero AND every REGEX matches its
combined stdout/stderr. With --success, COMMAND must exit 0 instead. CTest's
PASS_REGULAR_EXPRESSION alone ignores the exit code, and WILL_FAIL alone
ignores the message; this checks both.
"""
import re
import subprocess
import sys


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    success = argv[0] == "--success"
    if success:
        argv = argv[1:]
    sep = argv.index("--")
    patterns, cmd = argv[:sep], argv[sep + 1:]
    if not patterns or not cmd:
        print(__doc__)
        return 2

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=300)
    out = proc.stdout
    ok = True
    if success and proc.returncode != 0:
        print(f"FAIL: the command exited {proc.returncode}; it was expected to succeed")
        ok = False
    if not success and proc.returncode == 0:
        print("FAIL: the command exited 0; it was expected to abort")
        ok = False
    for pat in patterns:
        if not re.search(pat, out):
            print(f"FAIL: output does not match {pat!r}")
            ok = False
    if not ok:
        print("---- output ----")
        print(out)
        return 1
    print(f"PASS: exit code {proc.returncode}; matched {len(patterns)} pattern(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
