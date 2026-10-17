#!/usr/bin/env python3
"""Check the boundary conditions of the Poisson tests' manufactured solutions.

    python3 tools/poisson_ms_check.py        (needs sympy; not run by CI)

Tests/Unit/test_Poisson.cpp builds phi on the terrain-following domain over the
ground h(x, y) = H0 X(x) Y(y), with s = (z - h)/(z_top - h):

  general:              phi = f P(s) + s q Q(s),  q = (z_top - h)(h_x f_x + h_y f_y)/(1 + |grad h|^2)
  zero-gradient ground: phi = f G(s),             G(0) = G'(0) = 0

For every lateral BC (periodic, Neumann, Dirichlet) in x and y and both top
BCs (Dirichlet, Neumann), this checks at random points: zero normal derivative
on the sloped ground, the top condition, and each lateral condition. Prints
OK or the first violations.
"""
import random

import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
L, M, H0, ZT = sp.Integer(1000), sp.Integer(800), sp.Integer(60), sp.Integer(500)
X = (1 - sp.cos(2 * sp.pi * x / L)) / 2
Y = (1 - sp.cos(2 * sp.pi * y / M)) / 2
h = H0 * X * Y
s = (z - h) / (ZT - h)
N2 = 1 + sp.diff(h, x) ** 2 + sp.diff(h, y) ** 2


def F(kind, v, Lv):
    return {"periodic": sp.cos(2 * sp.pi * v / Lv) + sp.Rational(3, 10) * sp.sin(2 * sp.pi * v / Lv),
            "neumann": sp.cos(sp.pi * v / Lv), "dirichlet": sp.sin(sp.pi * v / Lv)}[kind]


def phi_of(kx, ky, top, zero_gradient):
    f = F(kx, x, L) * F(ky, y, M)
    dirichlet = top == "dirichlet"
    if zero_gradient:
        return f * (sp.sin(sp.pi * s) ** 2 if dirichlet else 1 - sp.cos(sp.pi * s))
    q = (ZT - h) * (sp.diff(h, x) * sp.diff(f, x) + sp.diff(h, y) * sp.diff(f, y)) / N2
    P = sp.cos(sp.pi * s / 2) if dirichlet else sp.cos(sp.pi * s)
    Q = (1 - s) if dirichlet else (1 - s) ** 2
    return f * P + s * q * Q


def check(kx, ky, top, zero_gradient, rng):
    phi = phi_of(kx, ky, top, zero_gradient)
    g = [sp.diff(phi, v) for v in (x, y, z)]
    n = [-sp.diff(h, x), -sp.diff(h, y), 1]
    bad = []
    for _ in range(6):
        xv, yv = rng.uniform(0, 1000), rng.uniform(0, 800)
        sub = {x: xv, y: yv}
        zg = float(h.subs(sub))
        sub[z] = zg
        v = float(sum(a * b for a, b in zip(n, g)).subs(sub))
        if abs(v) > 1e-12:
            bad.append(("ground", v))
        sub[z] = 500
        v = float(phi.subs(sub)) if top == "dirichlet" else float(g[2].subs(sub))
        if abs(v) > 1e-12:
            bad.append(("top", v))
        zv = rng.uniform(zg, 500)
        for var, val0, kind in [(x, 0, kx), (x, 1000, kx), (y, 0, ky), (y, 800, ky)]:
            sb = {x: xv, y: yv, z: zv}
            sb[var] = val0
            if kind == "neumann":
                v = float((g[0] if var == x else g[1]).subs(sb))
            elif kind == "dirichlet":
                v = float(phi.subs(sb))
            else:
                sb2 = dict(sb)
                sb2[var] = (1000 if var == x else 800) - val0
                v = float(phi.subs(sb) - phi.subs(sb2))
            if abs(v) > 1e-12:
                bad.append((str(var), val0, kind, v))
    return bad


def main():
    rng = random.Random(1)
    ok = True
    for zero_gradient in (False, True):
        for kx in ("periodic", "neumann", "dirichlet"):
            for ky in ("periodic", "neumann", "dirichlet"):
                for top in ("dirichlet", "neumann"):
                    bad = check(kx, ky, top, zero_gradient, rng)
                    ok &= not bad
                    label = "zero-gradient" if zero_gradient else "general"
                    print(f"{label:14s} x {kx:9s} y {ky:9s} top {top:9s}", "OK" if not bad else bad[:3])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
