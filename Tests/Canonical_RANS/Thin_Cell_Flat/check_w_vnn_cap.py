#!/usr/bin/env python3
"""The w rows' VNN cap (simplec.w_vnn_cap) on the thin-cell flat strip.

    check_w_vnn_cap.py --cap V --min-ratio R --tol T -- COMMAND [ARGS ...]

COMMAND is the application on Tests/Canonical_RANS/Thin_Cell_Flat/inputs_thin_cell_flat
(with its launcher and flags). The script runs it, each run in its own
directory under the current one (emptied first):

  off/       simplec.w_vnn_cap = 0, to the steady verdict;
  on/        simplec.w_vnn_cap = V, to the steady verdict;
  first/     the cap on, 40 iterations; then continued from its checkpoint to 80;
  straight/  the cap on, 80 iterations.

and checks:
  - both runs reach the verdict, and the banner names the cap ("w_vnn_cap off"
    and "w_vnn_cap V"); the capped run prints no terrain warning (flat ground,
    although its aspect ratio of about 134 is above the tested 50.78);
  - the capped run needs at most 1/R of the iterations of the run without it;
  - u, v and w of the two steady states agree to T (m/s, max over the cells).
    The cap changes the pseudo-time path, not the fixed point. T covers the
    convergence error of the run without the cap, whose slow mode leaves
    about 6e-5 m/s at the verdict (the capped run is within 5e-11 m/s of its
    fixed point, a 10000-iteration run);
  - the continued and straight plotfiles are identical in every field (the
    cap keeps the restart bit for bit).
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke"]


def run(cmd, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    full = cmd + extra
    print("check_w_vnn_cap: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_w_vnn_cap: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def last_plotfile(stdout, where):
    wrote = [ln.split()[-1] for ln in stdout.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln]
    if not wrote:
        raise SystemExit("check_w_vnn_cap: no plotfile written in %s/" % where)
    return os.path.join(where, wrote[-1])


def max_diff(a, b, f):
    worst = 0.0
    for pa, pb in zip(a[f], b[f]):
        for ca, cb in zip(pa, pb):
            for x, y in zip(ca, cb):
                worst = max(worst, abs(x - y))
    return worst


def main(argv):
    cap, ratio, tol = None, None, None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] in ("--cap", "--min-ratio", "--tol") and i + 1 < len(argv):
            v = float(argv[i + 1])
            if argv[i] == "--cap":
                cap = v
            elif argv[i] == "--min-ratio":
                ratio = v
            else:
                tol = v
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if cap is None or ratio is None or tol is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    fresh = ["init.first_guess=profile", "steady.max_iterations=20000"]
    ok = True

    off = run(cmd, "off", ["simplec.w_vnn_cap=0"] + fresh)
    on = run(cmd, "on", ["simplec.w_vnn_cap=%g" % cap] + fresh)
    n_off, n_on = verdict(off), verdict(on)
    good = n_off is not None and n_on is not None
    print("  steady verdict without and with the cap: %s, %s: %s" % (n_off, n_on, "pass" if good else "FAIL"))
    if not good:
        raise SystemExit("check_w_vnn_cap: FAIL")
    good = ("w_vnn_cap off," in off) and (("w_vnn_cap %g (slope margin" % cap) in on)
    ok = ok and good
    print("  the stage banner names the cap: %s" % ("pass" if good else "FAIL"))
    good = "WARNING: simplec.w_vnn_cap" not in on
    ok = ok and good
    print("  no terrain warning on flat ground: %s" % ("pass" if good else "FAIL"))
    good = n_on * ratio <= n_off
    ok = ok and good
    print("  iterations %d with the cap, %d without (ratio %.1f >= %g): %s"
          % (n_on, n_off, n_off / n_on, ratio, "pass" if good else "FAIL"))
    _, a = ep.read_fields(last_plotfile(off, "off"), ["u", "v", "w"])
    _, b = ep.read_fields(last_plotfile(on, "on"), ["u", "v", "w"])
    for f in ("u", "v", "w"):
        d = max_diff(a, b, f)
        good = d <= tol
        ok = ok and good
        print("  max |%s(cap) - %s(no cap)| = %.3e m/s (<= %g): %s" % (f, f, d, tol, "pass" if good else "FAIL"))

    capped = ["simplec.w_vnn_cap=%g" % cap]
    first = run(cmd, "first", capped + ["init.first_guess=profile", "steady.max_iterations=40"])
    chk = [ln.split()[-1] for ln in first.splitlines() if ln.startswith("gust: wrote ") and "_chk" in ln]
    if not chk:
        raise SystemExit("check_w_vnn_cap: no checkpoint written in first/")
    cont = subprocess.run(cmd + capped + ["steady.max_iterations=80", "io.restart=" + chk[-1]], cwd="first",
                          check=False, capture_output=True, text=True)
    if cont.returncode != 0:
        sys.stdout.write(cont.stdout[-4000:] + cont.stderr[-4000:])
        raise SystemExit("check_w_vnn_cap: the restart failed with code %d" % cont.returncode)
    straight = run(cmd, "straight", capped + ["init.first_guess=profile", "steady.max_iterations=80"])
    _, c = ep.read_fields(last_plotfile(cont.stdout, "first"), FIELDS)
    _, s = ep.read_fields(last_plotfile(straight, "straight"), FIELDS)
    for f in FIELDS:
        d = max_diff(c, s, f)
        good = d == 0
        ok = ok and good
        print("  restart: max |%s(continued) - %s(straight)| = %.3e (bit for bit): %s"
              % (f, f, d, "pass" if good else "FAIL"))
    if not ok:
        raise SystemExit("check_w_vnn_cap: FAIL")
    print("check_w_vnn_cap: pass")


if __name__ == "__main__":
    main(sys.argv)
