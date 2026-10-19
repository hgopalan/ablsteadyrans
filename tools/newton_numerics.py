#!/usr/bin/env python3
"""Re-derive the two numerical claims of the Newton-Krylov core (sub-phase 3.1).

1. The step of the forward-difference Jacobian-vector product. For
   R(q) = A q + B(q, q) the forward difference
       J_FD v = (R(q + eps v) - R(q)) / eps
   has the truncation error eps B(v, v) and a round-off error of about
   eps_mach |R| / eps. The sum is smallest at eps ~ sqrt(eps_mach), so the
   default h = sqrt(machine epsilon) of the working precision sits near the
   bottom of the V for double and for single precision.

2. The convergence of Newton's method. With exact (or tightly solved) linear
   steps the residual falls quadratically, r_{n+1} <= C r_n^2; with a Krylov
   tolerance eta the rate is linear with factor about eta (inexact Newton).

Usage: python3 tools/newton_numerics.py   (needs numpy only)
"""

import numpy as np


def quadratic_residual(n, dtype, seed=1):
    rng = np.random.default_rng(seed)
    A = (np.eye(n) * 4 + rng.standard_normal((n, n)) / np.sqrt(n)).astype(dtype)
    C = (rng.standard_normal((n, n)) / n).astype(dtype)

    def R(q):
        return A @ q + (C @ q) * q   # B(q, q) = (C q) * q

    def J(q, v):
        return A @ v + (C @ v) * q + (C @ q) * v

    def B(v):
        return (C @ v) * v

    return R, J, B


def eps_sweep(dtype):
    n = 64
    R, J, B = quadratic_residual(n, dtype)
    rng = np.random.default_rng(2)
    q = (1 + rng.random(n)).astype(dtype)
    v = rng.standard_normal(n).astype(dtype)
    jv = J(q.astype(np.float64), v.astype(np.float64))
    mach = np.finfo(dtype).eps
    hs = 10.0 ** np.arange(-1.0, -13.0, -0.5)
    hdef = np.sqrt(mach)
    rows = []
    for h in list(hs) + [hdef]:
        # The step of jvp: eps = h (1 + mean|q|) / rms(v).
        eps = dtype(h * (1 + np.mean(np.abs(q))) / np.sqrt(np.mean(v.astype(np.float64) ** 2)))
        jfd = (R(q + eps * v) - R(q)) / eps
        err = np.linalg.norm(jfd.astype(np.float64) - jv) / np.linalg.norm(jv)
        trunc = float(eps) * np.linalg.norm(B(v.astype(np.float64))) / np.linalg.norm(jv)
        rows.append((h, float(eps), err, trunc))
    return rows, hdef


def newton_column(eta, steps=8):
    """Steady 1-D Burgers-type problem -nu u'' + u u' = f, Dirichlet ends."""
    n, nu = 64, 0.05
    x = np.linspace(0, 1, n + 2)[1:-1]
    dx = 1.0 / (n + 1)
    f = np.sin(np.pi * x) * 2

    def R(u):
        up = np.concatenate(([0.0], u, [0.0]))
        return (-nu * (up[2:] - 2 * up[1:-1] + up[:-2]) / dx**2
                + up[1:-1] * (up[2:] - up[:-2]) / (2 * dx) - f)

    def Jmat(u):
        up = np.concatenate(([0.0], u, [0.0]))
        main = 2 * nu / dx**2 + (up[2:] - up[:-2]) / (2 * dx)
        lower = -nu / dx**2 - u[1:] / (2 * dx)
        upper = -nu / dx**2 + u[:-1] / (2 * dx)
        return np.diag(main) + np.diag(lower, -1) + np.diag(upper, 1)

    u = np.zeros(n)
    res = []
    for _ in range(steps):
        r = R(u)
        res.append(np.linalg.norm(r))
        if res[-1] < 1e-13 * res[0]:
            break
        J = Jmat(u)
        du = np.linalg.solve(J, -r)
        if eta > 0:
            # An inexact step: the linear residual |J du + r| = eta |r|, along
            # a fixed direction (the worst case of GMRES stopping early).
            e = np.ones(n)
            e *= eta * np.linalg.norm(r) / np.linalg.norm(e)
            du = np.linalg.solve(J, -r + e)
        u = u + du
    res.append(np.linalg.norm(R(u)))
    return res


def main():
    for dtype in (np.float64, np.float32):
        rows, hdef = eps_sweep(dtype)
        print(f"FD J v error, {np.dtype(dtype).name}: h, eps, |J_FD v - J v| / |J v|, eps |B| / |J v|")
        for h, eps, err, trunc in rows:
            tag = "  (default sqrt(eps_mach))" if h == hdef else ""
            print(f"  {h:9.3g}  {eps:9.3g}  {err:9.3g}  {trunc:9.3g}{tag}")
        errs = [r[2] for r in rows[:-1]]
        imin = int(np.argmin(errs))
        print(f"  minimum {errs[imin]:.3g} at h = {rows[imin][0]:.3g}; default gives {rows[-1][2]:.3g}"
              f" ({rows[-1][2] / errs[imin]:.2f} x the minimum)\n")
        assert rows[-1][2] < 3 * errs[imin]

    for eta in (0.0, 1e-2, 1e-1):
        res = newton_column(eta)
        print(f"Newton on a 1-D Burgers column, Krylov tolerance {eta:g}: |R_n|, r_(n+1)/r_n, r_(n+1)/r_n^2")
        for k in range(len(res) - 1):
            print(f"  {k}  {res[k]:10.4g}  {res[k + 1] / res[k]:10.3g}  {res[k + 1] / res[k] ** 2:10.3g}")
        print()


if __name__ == "__main__":
    main()
