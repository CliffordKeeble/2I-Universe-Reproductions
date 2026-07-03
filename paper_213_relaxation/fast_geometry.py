"""
fast_geometry.py  —  vectorised Regge curvature for Paper 213 Approach B.

Batches the per-tet embedding / dihedral / volume across ALL tets with numpy, so an
energy evaluation is a handful of array ops instead of a Python loop over tets. This is
what makes the Stage-1 L-BFGS variance minimisation practical at N >= 150.

Verified against the scalar geometry.py (edge_deficits / total_volume) — same numbers.
"""

import numpy as np
from geometry import TET_EDGES, TET_EDGE_OPP

# precomputed local (a,b,c,d) for each of the 6 tet edges (b-a is the hinge; c,d opposite)
_ABCD = []
for (a, b) in TET_EDGES:
    c, d = TET_EDGE_OPP[(a, b)]
    _ABCD.append((a, b, c, d))


def precompute(cx):
    """Per-tet array of its 6 global edge indices (TET_EDGES order). Shape (T,6)."""
    idx = np.empty((len(cx.tets), 6), dtype=np.int64)
    for ti, t in enumerate(cx.tets):
        for k, (a, b) in enumerate(TET_EDGES):
            idx[ti, k] = cx.edge_index[tuple(sorted((t[a], t[b])))]
    return dict(tet_edge_idx=idx, n_edges=len(cx.edges))


def _embed_batch(L6):
    """Batched tet embedding from edge lengths L6 (T,6). Returns P (T,4,3) and valid mask."""
    l01, l02, l03, l12, l13, l23 = [L6[:, k] for k in range(6)]
    eps = 1e-12
    valid = l01 > eps
    x2 = np.where(valid, (l01**2 + l02**2 - l12**2) / (2 * np.where(valid, l01, 1)), 0.0)
    y2sq = l02**2 - x2**2
    valid &= y2sq > eps
    y2 = np.sqrt(np.where(valid, y2sq, 1.0))
    x3 = np.where(valid, (l01**2 + l03**2 - l13**2) / (2 * np.where(valid, l01, 1)), 0.0)
    y3 = np.where(valid, (l02**2 + l03**2 - l23**2 - 2 * x2 * x3) / (2 * np.where(valid, y2, 1)), 0.0)
    z3sq = l03**2 - x3**2 - y3**2
    valid &= z3sq > eps
    z3 = np.sqrt(np.where(valid, z3sq, 1.0))
    T = L6.shape[0]
    P = np.zeros((T, 4, 3))
    P[:, 1, 0] = l01
    P[:, 2, 0] = x2; P[:, 2, 1] = y2
    P[:, 3, 0] = x3; P[:, 3, 1] = y3; P[:, 3, 2] = z3
    return P, valid


def fast_curvature(ell, pre):
    """Vectorised (K, deficits, Vd, n_invalid). K_e = delta_e / (Vd_e/ell_e)."""
    idx = pre["tet_edge_idx"]
    L6 = ell[idx]                      # (T,6)
    P, valid = _embed_batch(L6)
    # volumes
    e = P[:, 1:, :] - P[:, 0:1, :]     # (T,3,3)
    V = np.abs(np.linalg.det(e)) / 6.0
    V = np.where(valid, V, 0.0)
    n_invalid = int(np.sum(~valid))
    # dihedral angles, 6 per tet
    ang = np.zeros((P.shape[0], 6))
    for k, (a, b, c, d) in enumerate(_ABCD):
        u = P[:, b] - P[:, a]
        un = u / np.linalg.norm(u, axis=1, keepdims=True).clip(1e-15)
        vv = P[:, c] - P[:, a]
        ww = P[:, d] - P[:, a]
        vp = vv - np.sum(vv * un, axis=1, keepdims=True) * un
        wp = ww - np.sum(ww * un, axis=1, keepdims=True) * un
        nv = np.linalg.norm(vp, axis=1).clip(1e-15)
        nw = np.linalg.norm(wp, axis=1).clip(1e-15)
        cos = np.clip(np.sum(vp * wp, axis=1) / (nv * nw), -1.0, 1.0)
        ang[:, k] = np.where(valid, np.arccos(cos), 0.0)

    nE = pre["n_edges"]
    deficits = np.full(nE, 2 * np.pi)
    Vd = np.zeros(nE)
    for k in range(6):
        np.add.at(deficits, idx[:, k], -ang[:, k])
        np.add.at(Vd, idx[:, k], V / 6.0)
    a_e = np.where(Vd > 1e-12, Vd / ell, np.nan)
    K = deficits / a_e
    return K, deficits, Vd, n_invalid


def fast_energy(ell, pre):
    K, deficits, Vd, ninv = fast_curvature(ell, pre)
    if ninv > 0:
        return np.inf, np.inf, np.nan, ninv
    good = np.isfinite(K) & (Vd > 1e-12)
    w = Vd[good]; Kg = K[good]
    Kbar = np.sum(w * Kg) / np.sum(w)
    E = np.sum(w * (Kg - Kbar) ** 2) / np.sum(w)
    cov = np.sqrt(E) / abs(Kbar) if abs(Kbar) > 1e-12 else np.inf
    return float(E), float(cov), float(Kbar), 0


if __name__ == "__main__":
    # verify against scalar geometry.py
    import time
    from geometry import generic_s3, edge_deficits, total_volume
    from ricci_flow import curvature
    cx = generic_s3(150, seed=20260618)
    pre = precompute(cx)
    ell = np.ones(len(cx.edges))
    Ks, ds, Vds, ninv_s = curvature(cx, ell)
    Kf, df, Vdf, ninv_f = fast_curvature(ell, pre)
    print("VERIFY fast vs scalar (frustrated ell=1):")
    print(f"  max|deficit diff| = {np.nanmax(np.abs(ds - df)):.2e}")
    print(f"  max|Vdual diff|   = {np.nanmax(np.abs(Vds - Vdf)):.2e}")
    print(f"  max|K diff|       = {np.nanmax(np.abs(Ks - Kf)):.2e}")
    print(f"  n_invalid scalar={ninv_s} fast={ninv_f}")
    # round metric too
    P = cx.hull_points
    ellr = np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])
    Ksr, dsr, _, _ = curvature(cx, ellr)
    Kfr, dfr, _, _ = fast_curvature(ellr, pre)
    print(f"  (round) max|deficit diff| = {np.nanmax(np.abs(dsr - dfr)):.2e}")
    # speed
    t0 = time.time()
    for _ in range(200):
        fast_curvature(ell, pre)
    print(f"  fast_curvature: {(time.time()-t0)/200*1000:.2f} ms/call")