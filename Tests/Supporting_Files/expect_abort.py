#!/usr/bin/env python3
"""Run a command and check both its exit code and its output.

    expect_abort.py [--success] [--log FILE] [--timeout SECONDS] REGEX [REGEX ...] -- COMMAND [ARGS ...]

By default passes only if COMMAND exits non-zero AND every REGEX matches its
combined stdout/stderr. With --success, COMMAND must exit 0 instead. A REGEX
written as !REGEX must NOT match. CTest's
PASS_REGULAR_EXPRESSION alone ignores the exit code, and WILL_FAIL alone
ignores the message; this checks both. --log writes the combined output to
FILE (for a check that reads the run's log afterwards); --timeout (default
300 s) stops COMMAND.
"""
import re
import subprocess
import sys


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    success, log, timeout = False, None, 300.0
    while argv and argv[0] in ("--success", "--log", "--timeout"):
        if argv[0] == "--success":
            success, argv = True, argv[1:]
        elif len(argv) < 2:
            print(__doc__)
            return 2
        elif argv[0] == "--log":
            log, argv = argv[1], argv[2:]
        else:
            timeout, argv = float(argv[1]), argv[2:]
    if "--" not in argv or not timeout > 0:
        print(__doc__)
        return 2
    sep = argv.index("--")
    patterns, cmd = argv[:sep], argv[sep + 1:]
    if not patterns or not cmd:
        print(__doc__)
        return 2

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=timeout)
    out = proc.stdout
    if log:
        with open(log, "w") as fh:
            fh.write(out)
    ok = True
    if success and proc.returncode != 0:
        print(f"FAIL: the command exited {proc.returncode}; it was expected to succeed")
        ok = False
    if not success and proc.returncode == 0:
        print("FAIL: the command exited 0; it was expected to abort")
        ok = False
    for pat in patterns:
        if pat.startswith("!"):
            if re.search(pat[1:], out):
                print(f"FAIL: output matches {pat[1:]!r}, which it must not")
                ok = False
        elif not re.search(pat, out):
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
