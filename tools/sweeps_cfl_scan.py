#!/usr/bin/env python3
"""
Sub-phase 1.16 (Pressel lens): the stability of the lagged z-line sweeps.

Von Neumann analysis of one SIMPLEC pseudo-time outer iteration on a 2-D
staggered (x-z) grid, uniform advection U in x, periodic.

State per wavenumber: (u_hat, w_hat, p_hat); u on x-faces, w on z-faces,
p at cell centres (exact staggered phase factors).

Outer iteration:
  predictor : (1/dtau + A1) du = -R(q^n), A1 = 1st-order upwind + diffusion,
              solved with m z-line Jacobi sweeps (z part exact, x neighbours
              lagged, zero at sweep 1)  -> du = S * (-R), S = sum_j (-M^-1 N)^j M^-1
  u*  = u^n + du
  SIMPLEC: 1/d = beta (= dtau unless 'dmode' says otherwise)
           D (beta G phi) = D u*,  u^{n+1} = u* - beta G phi,  p^{n+1} = p^n + phi
Optional extensions (all off in the baseline requested):
  shear  : sigma = dtau*dU/dz, explicit term  sigma*w  in the u residual
           (not in A1)  -- frozen-coefficient local model of a sheared base flow
  dmode  : 'simplec' (beta = dtau) or 'zonly' (d = a_P - sum a_nb(z only),
           i.e. x neighbours not subtracted because they are lagged)
  wscheme: residual scheme for the w equation different from the u equation
"""
import numpy as np
import itertools, sys

def adv_symbol(scheme, th):
    """symbol of  (U/dx) * delta_x  for scalar advection, U>0"""
    e = np.exp(-1j*th)
    if scheme == 'upwind1':
        return 1.0 - e
    if scheme == 'LUD':
        return (3.0 - 4.0*e + e*e)/2.0
    if scheme == 'central':
        return 1j*np.sin(th)
    raise ValueError(scheme)

def lagged_solve_symbol(thx, thz, c, m, dnu, rz):
    """S symbol (times 1/dtau units) for m z-line Jacobi sweeps on
    (1/dtau + A1), A1 = upwind1 advection + diffusion, x neighbours lagged.
    Everything scaled by dtau: M' = dtau*M, N' = dtau*N.  Returns dtau*S."""
    # x neighbour (upwind) : -c e^{-i thx}; diffusion x neighbours: -dnu (e^{i}+e^{-i})
    ex = np.exp(-1j*thx)
    Mx_diag = 1.0 + c + 2.0*dnu                     # 1/dtau + U/dx + 2 nu/dx^2
    dnuz = dnu/(rz*rz)                              # nu dtau/dz^2
    Mz = 2.0*dnuz*(1.0 - np.cos(thz))               # z tridiagonal (exact)
    M = Mx_diag + Mz
    N = -c*ex - 2.0*dnu*np.cos(thx)                 # lagged x neighbours
    r = -N/M
    # sum_{j=0}^{m-1} r^j / M
    if m == 0:
        return 1.0/M  # never used
    s = np.where(np.abs(1.0 - r) > 1e-14, (1.0 - r**m)/(1.0 - r), float(m))
    return s/M

def amplification(thx, thz, c, m, scheme, wscheme=None, dnu=0.0, rz=1.0,
                  sigma=0.0, dmode='simplec'):
    """3x3 amplification matrix for state (u,w,p) (p scaled by dtau/dx so that
    G p enters as c*Gx' p with Gx' = 2i sin(thx/2))."""
    if wscheme is None:
        wscheme = scheme
    # staggered divergence / gradient symbols (dx = 1, dz = rz)
    Gx = 2j*np.sin(thx/2.0)
    Gz = 2j*np.sin(thz/2.0)/rz
    Dx, Dz = Gx, Gz
    L = Dx*Gx + Dz*Gz                                # = -(|Gx|^2+|Gz|^2)
    # residual symbols scaled by dtau: dtau*R
    diff = 2.0*dnu*(1.0 - np.cos(thx)) + 2.0*dnu/(rz*rz)*(1.0 - np.cos(thz))
    au = c*adv_symbol(scheme, thx) + diff
    aw = c*adv_symbol(wscheme, thx) + diff
    # dtau*R = Rm @ (u,w,p)  with p already scaled so that gradient is G p
    Rm = np.array([[au, sigma, Gx],
                   [0.0, aw, Gz]], dtype=complex)
    S = lagged_solve_symbol(thx, thz, c, m, dnu, rz)  # dtau*S
    # predictor:  (u*,w*) = (u,w) - S * dtau R
    Pred = np.zeros((2, 3), dtype=complex)
    Pred[0, 0] = 1.0; Pred[1, 1] = 1.0
    Pred -= S*Rm
    # SIMPLEC beta (scaled by 1/dtau): 1/(d dtau)
    if dmode == 'simplec':
        beta = 1.0
    elif dmode == 'zonly':
        beta = 1.0/(1.0 + c)   # d = a_P - sum a_nb(z) = 1/dtau + U/dx
    else:
        raise ValueError(dmode)
    # phi (scaled) = (D u*)/(beta L); u = u* - beta G phi ; p += phi
    Dop = np.array([Dx, Dz], dtype=complex)
    Gop = np.array([Gx, Gz], dtype=complex)
    div = Dop @ Pred                                 # 1x3
    phi = div/(beta*L)                                # 1x3
    Unew = Pred - np.outer(Gop, phi)                  # 2x3
    Pnew = np.array([0.0, 0.0, 1.0], dtype=complex) + phi
    return np.vstack([Unew, Pnew])

def scan(c, m, scheme, nth=96, **kw):
    thetas = np.linspace(-np.pi, np.pi, nth, endpoint=False)
    worst = (0.0, None, None, None)
    for thx in thetas:
        for thz in thetas:
            if abs(thx) < 1e-12 and abs(thz) < 1e-12:
                continue
            G = amplification(thx, thz, c, m, scheme, **kw)
            ev, V = np.linalg.eig(G)
            k = np.argmax(np.abs(ev))
            rho = abs(ev[k])
            if rho > worst[0]:
                worst = (rho, thx, thz, V[:, k])
    return worst

def fmt(w):
    rho, thx, thz, v = w
    if v is None:
        return f"{rho:7.4f}"
    v = v/np.max(np.abs(v))
    return (f"{rho:7.4f}  (kx dx={thx:+.3f}, kz dz={thz:+.3f})  "
            f"|u|={abs(v[0]):.2f} |w|={abs(v[1]):.2f} |p|={abs(v[2]):.2f}")

def table(title, cfls, ms, schemes, **kw):
    print("\n=== " + title + " ===")
    for scheme in schemes:
        print(f"--- residual scheme: {scheme}")
        print("  CFL   m    rho_max   worst mode")
        for c in cfls:
            for m in ms:
                print(f"  {c:4.1f} {m:3d}  {fmt(scan(c, m, scheme, **kw))}")

def thresholds(schemes=('LUD', 'central'), ms=(2, 4, 8, 16), step=0.05):
    """The first CFL (step 0.05) at which the scalar lagged-sweep factor
    1 - S a exceeds 1 in modulus over theta in (0, pi]: the table of
    docs/predictor.rst (sub-phase 1.16)."""
    th = np.linspace(0.0, np.pi, 721)[1:]
    print("first CFL with an amplification above 1 (scalar factor, pure advection)")
    print("  m    " + "  ".join(f"{s:>8s}" for s in schemes))
    for m in ms:
        row = []
        for scheme in schemes:
            c = step
            while c < 20.0:
                g = 1.0 - lagged_solve_symbol(th, 0.0, c, m, 0.0, 1.0) * c * adv_symbol(scheme, th)
                if np.abs(g).max() > 1.0 + 1e-12:
                    break
                c += step
            row.append(c)
        print(f"  {m:2d}   " + "  ".join(f"{c:8.2f}" for c in row))

if __name__ == '__main__':
    if '--thresholds' in sys.argv:
        thresholds()
        sys.exit(0)
    cfls = [0.5, 1.0, 2.0, 4.0]
    ms = [1, 2, 4, 8, 16, 64]
    schemes = ['upwind1', 'LUD', 'central']
    quick = '--quick' in sys.argv
    nth = 48 if quick else 96

    # sanity: projection exactness -> div-free mode eigenvalue equals 1 - S a
    thx, thz, c, m = 0.7, 1.3, 2.0, 4
    G = amplification(thx, thz, c, m, 'LUD')
    S = lagged_solve_symbol(thx, thz, c, m, 0.0, 1.0)
    g1 = 1.0 - S*c*adv_symbol('LUD', thx)
    ev = np.linalg.eigvals(G)
    print("sanity: eigenvalues of 3x3 =", np.round(ev, 5))
    print("        1 - S a (1-D scalar) =", np.round(g1, 5),
          "  1 - S/dtau (pressure mode) =", np.round(1.0 - S, 5))

    table("baseline: pure advection, dz/dx = 1, SIMPLEC d = 1/dtau",
          cfls, ms, schemes, nth=nth)
    table("baseline, dz/dx = 0.4", cfls, ms, schemes, nth=nth, rz=0.4)
    table("small diffusion nu dtau/dx^2 = 0.05, dz/dx = 0.4",
          cfls, ms, schemes, nth=nth, rz=0.4, dnu=0.05)

    # extensions: what the uniform-flow model leaves out
    table("EXT-1: SIMPLEC d built from z neighbours only (x lagged) ",
          [1.0, 2.0], [4, 16], ['LUD'], nth=nth, rz=0.4, dmode='zonly')
    for sig in [0.25, 0.5, 1.0, 2.0]:
        table(f"EXT-2: explicit shear coupling sigma = dtau dU/dz = {sig}, dz/dx=0.4",
              [0.5, 1.0, 2.0, 4.0], [4, 16, 64], ['LUD'], nth=nth, rz=0.4, sigma=sig)
    table("EXT-3: u residual LUD, w residual central (limiter acting differently)",
          [1.0, 2.0], [4, 16], ['LUD'], nth=nth, rz=0.4, wscheme='central')
