"""
relax_fast.py  —  O(N) local-finite-difference gradient for the Stage-1 FSS ladder.

The Stage-1 energy is  E = A/Vtot - (S/Vtot)^2  with
  A = sum_e (delta_e * ell_e)^2 / Vd_e ,  S = sum_e delta_e ell_e ,  Vtot = sum_e Vd_e
(identical to relax_geometry's variance energy; algebra in RESULTS). Perturbing edge f
changes ONLY the tets containing f, so dE/d ell_f is computed by recomputing those ~5
tets — O(1) per partial, O(N) for the whole gradient (vs O(N^2) for a global numerical
gradient). This does NOT change the scheme; it only makes N=10^3 tractable. Verified
against the global numerical gradient before use.
"""

import numpy as np
from scipy.optimize import minimize
from geometry import tet_embed, tet_dihedrals, tet_volume, TET_EDGES


def precompute_fss(cx):
    T = len(cx.tets)
    tet_edge_idx = np.empty((T, 6), dtype=np.int64)
    for ti, t in enumerate(cx.tets):
        for k, (a, b) in enumerate(TET_EDGES):
            tet_edge_idx[ti, k] = cx.edge_index[tuple(sorted((t[a], t[b])))]
    edge_tets = [[] for _ in range(len(cx.edges))]
    for ti in range(T):
        for e in tet_edge_idx[ti]:
            edge_tets[e].append(ti)
    return dict(tet_edge_idx=tet_edge_idx, edge_tets=edge_tets,
                nE=len(cx.edges), T=T)


def base_geometry(ell, pre):
    """deficits, Vd, per-tet base dihedrals+volumes, and (A,S,Vtot,E)."""
    tei = pre["tet_edge_idx"]
    T, nE = pre["T"], pre["nE"]
    deficits = np.full(nE, 2 * np.pi)
    Vd = np.zeros(nE)
    base_ang = np.empty((T, 6))
    base_vol = np.empty(T)
    for ti in range(T):
        L = ell[tei[ti]]
        ang = tet_dihedrals(L)
        vol = tet_volume(L)
        base_ang[ti] = ang
        base_vol[ti] = vol
        for k in range(6):
            deficits[tei[ti, k]] -= ang[k]
            Vd[tei[ti, k]] += vol / 6.0
    A = np.sum((deficits * ell) ** 2 / Vd)
    S = np.sum(deficits * ell)
    Vtot = np.sum(Vd)
    E = A / Vtot - (S / Vtot) ** 2
    return dict(deficits=deficits, Vd=Vd, base_ang=base_ang, base_vol=base_vol,
                A=A, S=S, Vtot=Vtot, E=E)


def energy_grad(ell, pre, h=1e-6):
    """Return (E, grad) with grad by O(N) local central differences."""
    tei = pre["tet_edge_idx"]
    bg = base_geometry(ell, pre)
    deficits, Vd = bg["deficits"], bg["Vd"]
    base_ang, base_vol = bg["base_ang"], bg["base_vol"]
    A0, S0, V0, E0 = bg["A"], bg["S"], bg["Vtot"], bg["E"]
    grad = np.zeros(pre["nE"])

    def perturbed_E(f, dh):
        dA = dS = dVtot = 0.0
        # affected edges accumulate Δδ, ΔVd
        acc = {}
        for ti in pre["edge_tets"][f]:
            edges_t = tei[ti]
            L = ell[edges_t].copy()
            pos = np.where(edges_t == f)[0][0]
            L[pos] += dh
            ang = tet_dihedrals(L)
            vol = tet_volume(L)
            if np.any(np.isnan(ang)) or np.isnan(vol):
                return np.inf
            dVtot += vol - base_vol[ti]
            for k in range(6):
                e = edges_t[k]
                dd, dvd = acc.get(e, (0.0, 0.0))
                acc[e] = (dd - (ang[k] - base_ang[ti, k]), dvd + (vol - base_vol[ti]) / 6.0)
        A = A0; S = S0
        for e, (dd, dvd) in acc.items():
            le = ell[e] + (dh if e == f else 0.0)
            old = (deficits[e] * ell[e]) ** 2 / Vd[e]
            newVd = Vd[e] + dvd
            new = ((deficits[e] + dd) * le) ** 2 / newVd
            A += new - old
            S += (deficits[e] + dd) * le - deficits[e] * ell[e]
        Vt = V0 + dVtot
        return A / Vt - (S / Vt) ** 2

    for f in range(pre["nE"]):
        Ep = perturbed_E(f, h)
        Em = perturbed_E(f, -h)
        grad[f] = (Ep - Em) / (2 * h)
    return E0, grad


def relax_fss(cx, ell0, pre=None, iters=200, verbose=False):
    """Fixed-volume variance minimisation with the O(N) gradient; volume renormalised."""
    if pre is None:
        pre = precompute_fss(cx)
    V_target = base_geometry(ell0, pre)["Vtot"]
    x0 = np.log(ell0)

    def renorm(x):
        ell = np.exp(x)
        bg = base_geometry(ell, pre)
        if bg["Vtot"] > 1e-12:
            ell = ell * (V_target / bg["Vtot"]) ** (1.0 / 3.0)
        return ell

    def fun(x):
        ell = renorm(x)
        E, g = energy_grad(ell, pre)
        if not np.isfinite(E):
            return 1e6, np.zeros_like(x)
        # chain rule x=log ell, ell=exp(x): dE/dx = dE/dell * ell  (renorm is ~scale, drop)
        return E, g * ell

    res = minimize(fun, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": iters, "maxfun": iters * 3})
    ell = renorm(res.x)
    E_final = base_geometry(ell, pre)["E"]
    return ell, dict(E0=float(base_geometry(ell0, pre)["E"]), Ef=float(E_final),
                     iters=res.nit, msg=str(res.message), success=bool(res.success))


if __name__ == "__main__":
    # verify the local-FD gradient against scipy's global numerical gradient (N=100)
    import time
    from geometry import generic_s3
    from scipy.optimize import approx_fprime
    cx = generic_s3(100, seed=20260618)
    pre = precompute_fss(cx)
    rng = np.random.default_rng(1)
    ell = np.exp(0.1 * rng.standard_normal(len(cx.edges)))  # generic point
    t0 = time.time(); E, g = energy_grad(ell, pre); tg = time.time() - t0
    f = lambda x: base_geometry(x, pre)["E"]
    gn = approx_fprime(ell, f, 1e-6)
    rel = np.linalg.norm(g - gn) / np.linalg.norm(gn)
    print(f"N=100 gradient check: ||g_local - g_numeric|| / ||g_numeric|| = {rel:.2e}")
    print(f"  local-FD gradient time = {tg*1000:.0f} ms  (numeric would be ~{len(cx.edges)*2*1.4:.0f} ms)")
    print(f"  match: {'PASS' if rel < 1e-4 else 'FAIL'}")
