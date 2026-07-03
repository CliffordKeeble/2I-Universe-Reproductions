"""
relax_geometry.py  —  Stage 1 flow, Option 1 (brief v1.1): fixed-volume curvature-
variance minimisation.

    minimise   E(ell) = sum_e w_e (K_e - Kbar)^2 / sum_e w_e       (w_e = V_e^dual)
    subject to total volume fixed.

Round S^3 is the global minimum (constant curvature => E = 0 in the continuum), so the
already-round metric is a fixed point BY CONSTRUCTION — this is the fix for the v1.0
instrument failure (the explicit-Euler Ricci flow neckpinched even a round metric).

Implementation: backtracking line search in the homogenising direction u = -(K - Kbar),
accepting only steps that DECREASE E and keep every tet embeddable, with volume
renormalised each step. No fixed Euler step to overshoot; at a round metric u ~ 0 and the
search makes no move => stays round. Variance is the engine; the SPECTRUM (spectral_s3)
is the verdict.

Standalone:  python relax_geometry.py [N] [iters]
"""

import sys
import numpy as np
from scipy.optimize import minimize
from geometry import generic_s3, total_volume
from fast_geometry import precompute, fast_energy, fast_curvature


def energy(cx, ell, pre=None):
    """Volume-weighted mean-square curvature deviation E, plus (CoV, Kbar, ninv)."""
    if pre is None:
        pre = precompute(cx)
    return fast_energy(ell, pre)


def _fast_volume(ell, pre):
    _, _, Vd, ninv = fast_curvature(ell, pre)
    return float(np.sum(Vd)), ninv     # sum of dual volumes = total volume


def variance_min(cx, ell0, pre=None, iters=120, verbose=True):
    """
    Minimise E over x = log(ell) with L-BFGS-B (numerical gradient, vectorised energy).
    Volume held fixed by renormalising inside the objective (harmless overall-scale flat
    direction). Non-embeddable tets -> large finite penalty so the optimiser steers away.
    """
    if pre is None:
        pre = precompute(cx)
    V_target, _ = _fast_volume(ell0, pre)
    x0 = np.log(ell0)
    hist = []

    def _renorm(x):
        ell = np.exp(x)
        Vn, _ = _fast_volume(ell, pre)
        if Vn > 1e-9:
            ell = ell * (V_target / Vn) ** (1.0 / 3.0)
        return ell

    def f(x):
        E, cov, Kbar, ninv = fast_energy(_renorm(x), pre)
        if ninv > 0 or not np.isfinite(E):
            return 1.0e3 * (1 + ninv)
        return E

    it_box = {"n": 0}

    def cb(xk):
        E, cov, Kbar, ninv = fast_energy(_renorm(xk), pre)
        hist.append(dict(it=it_box["n"], E=float(E), cov=float(cov), Kbar=float(Kbar)))
        if verbose and (it_box["n"] % max(1, iters // 12) == 0):
            print(f"  it={it_box['n']:4d}  E={E:.5e}  CoV(K)={cov:.4f}  Kbar={Kbar:+.4f}")
        it_box["n"] += 1

    E0, cov0, Kbar0, _ = fast_energy(ell0, pre)
    hist.append(dict(it=0, E=float(E0), cov=float(cov0), Kbar=float(Kbar0)))
    res = minimize(f, x0, method="L-BFGS-B",
                   options={"maxiter": iters, "maxfun": iters * (len(x0) + 5)},
                   callback=cb)
    ell = _renorm(res.x)
    Ef, covf, Kbarf, ninvf = fast_energy(ell, pre)
    hist.append(dict(it=it_box["n"], E=float(Ef), cov=float(covf), Kbar=float(Kbarf),
                     note=f"L-BFGS done: {res.message}", success=bool(res.success)))
    return ell, hist


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    cx = generic_s3(N, seed=20260618)
    pre = precompute(cx)
    P = cx.hull_points
    ell_round = np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])
    print(f"Stage 1 Option-1 variance minimisation (N={N})")
    for label, ell0 in (("frustrated ell=1", np.ones(len(cx.edges))),
                        ("round-init (control #3)", ell_round)):
        E0, cov0, Kbar0, _ = fast_energy(ell0, pre)
        print(f"\n[{label}] start E={E0:.5e} CoV={cov0:.4f} Kbar={Kbar0:+.4f}")
        ell, hist = variance_min(cx, ell0, pre=pre, iters=iters)
        print(f"[{label}] end   E={hist[-1]['E']:.5e} CoV={hist[-1]['cov']:.4f} "
              f"note={hist[-1].get('note','max-iters')}")


if __name__ == "__main__":
    main()
