"""
The one-line mutant runner shared by tools/transfer_mutants.py and
tools/continuation_mutants.py.

Each mutant replaces one piece of text (it must occur exactly once in its
file), rebuilds unit_tests, runs the given GoogleTest filter on two ranks with
the unit-test flags of Tests/CMakeLists.txt (SNaN fill, FPE trap), restores the
file and prints one row: id | change | outcome | failing tests. The outcome is
killed, killed (timeout), SURVIVED, survived (<reason>) for a mutant listed as
expected to survive, or not compiled. The unmutated code runs first and must
pass. The run exits 1 when a mutant not listed survives or does not compile.
"""
import os
import re
import signal
import subprocess
import sys


def build_and_run(build, gtest_filter, ranks=2):
    """Rebuild unit_tests and run gtest_filter. Returns (result, failing
    tests); result is "build error", "timeout", "crash", "launch error",
    "failed" or "passed"."""
    b = subprocess.run(["cmake", "--build", build, "--target", "unit_tests"], capture_output=True, text=True)
    if b.returncode:
        return "build error", []
    # Its own process group, so a timeout kills mpirun and every rank.
    p = subprocess.Popen(["mpirun", "--oversubscribe", "-np", str(ranks), os.path.join(build, "unit_tests"),
                          f"--gtest_filter={gtest_filter}", "amrex.the_arena_init_size=0", "amrex.init_snan=1",
                          "amrex.fpe_trap_invalid=1"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=build,
                         start_new_session=True)
    try:
        out, _ = p.communicate(timeout=900)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        return "timeout", []
    except BaseException:
        # Ctrl-C: mpirun is in its own session, so SIGINT does not reach it.
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        raise
    failed = sorted(set(re.findall(r"\[  FAILED  \] (\S+\.\S+)", out)))
    if failed:
        return "failed", failed
    if p.returncode:
        # A crash (FPE trap, abort) inside a test: the one that was running.
        started = re.findall(r"\[ RUN      \] (\S+\.\S+)", out)
        done = set(re.findall(r"\[       OK \] (\S+\.\S+)", out))
        running = [n for n in started if n not in done][:1]
        if running:
            return "crash", [f"{n} (crash)" for n in running]
        return "launch error", []
    return "passed", []


def run(mutants, gtest_filter, equivalent, doc):
    """mutants: (id, change, file, old text, new text) rows. equivalent: the
    ids whose survival is expected, as a set (reported "equivalent") or a dict
    id -> reason. Command line: BUILD_DIR [id ...]."""
    if not isinstance(equivalent, dict):
        equivalent = {m: "equivalent" for m in equivalent}
    if len(sys.argv) < 2:
        sys.exit(doc)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    build = os.path.abspath(sys.argv[1])
    only = sys.argv[2:]
    unknown = sorted(set(only) - {m[0] for m in mutants})
    if unknown:
        sys.exit(f"unknown mutant ids: {', '.join(unknown)}")
    # The unmutated code must pass, or every mutant would look killed.
    base, failed = build_and_run(build, gtest_filter)
    if base != "passed":
        sys.exit(f"baseline: {base} {', '.join(failed)}")
    survived = 0
    try:
        for mid, what, f, old, new in mutants:
            if only and mid not in only:
                continue
            path = os.path.join(root, f)
            with open(path) as fh:
                src = fh.read()
            if src.count(old) != 1:
                sys.exit(f"{mid}: the original text occurs {src.count(old)} times in {f}")
            try:
                with open(path, "w") as fh:
                    fh.write(src.replace(old, new))
                res, failed = build_and_run(build, gtest_filter)
            finally:
                with open(path, "w") as fh:
                    fh.write(src)
            if res == "launch error":
                sys.exit(f"{mid}: launch error (no test started)")
            res = {"failed": "killed", "crash": "killed", "timeout": "killed (timeout)", "passed": "SURVIVED",
                   "build error": "not compiled"}[res]
            if res == "SURVIVED" and mid in equivalent:
                res = f"survived ({equivalent[mid]})"
            survived += res in ("SURVIVED", "not compiled")
            print(f"{mid} | {what} | {res} | {', '.join(failed)}", flush=True)
    finally:
        # Leave the binary built from the unmutated files.
        subprocess.run(["cmake", "--build", build, "--target", "unit_tests"], capture_output=True)
    sys.exit(1 if survived else 0)
