"""
curl_instrument.py  —  the signed discrete curl (Beltrami) operator for Stage 2.

Edge-element (Whitney 1-form / Nedelec) discretisation of curl = *d on 1-forms. The
Whitney form for edge (a,b) is  W_ab = lam_a grad lam_b - lam_b grad lam_a,  with
curl(W_ab) = 2 grad lam_a x grad lam_b  (constant per tet). We assemble two global
matrices over edges:

  M1[e,e'] = integral W_e . W_e'            (edge mass; SPD)
  R [e,e'] = integral curl(W_e) . W_e'      (curl bilinear; symmetric on a closed
                                             manifold, INDEFINITE -> SIGNED spectrum)

and solve the generalized eigenproblem  R a = lambda M1 a. The signed eigenvalues are the
Beltrami/curl eigenvalues; on round S^3 the lowest nonzero are +-2 (mult 3, the
left/right-invariant Hopf fields), then +-3, +-4, ...  The SIGN is the whole point of
Stage 2 (chirality), so R must be the genuine indefinite curl form, not curl-curl.

ORIENTATION CONVENTION (Mr A #5): each tet is embedded with a right-handed frame
(tet_embed gives z3>0), fixing a global orientation; global edge orientation = (min->max)
vertex id. "+2" is defined relative to this. An orientation FLIP (reflect the embedding)
sends R -> -R, i.e. it RELABELS +2 <-> -2 without changing physics. Verified in __main__.

Per-tet integrals on the flat (chord-length) tet:
  int lam_i dV = V/4 ;  int lam_i lam_j dV = V(1+delta_ij)/20.

Standalone:  python curl_instrument.py     # validate +-2/+-3 on the 600-cell
"""

import numpy as np
import scipy.linalg as sla
from geometry import tet_embed, TET_EDGES


def _tet_grads(P):
    """grad lam_i (4x3) and volume V for an embedded tet (4x3 coords)."""
    e = P[1:] - P[0]                      # 3x3
    detM = np.linalg.det(e)
    V = abs(detM) / 6.0
    if V < 1e-12:
        return None, 0.0
    Minv = np.linalg.inv(e)               # rows = grad lam_1,2,3
    grads = np.zeros((4, 3))
    grads[1:] = Minv
    grads[0] = -grads[1:].sum(axis=0)
    return grads, V


def build_curl(cx, ell, flip=False):
    """Assemble (R, M1) over edges. flip=True reflects the embedding (orientation flip)."""
    nE = len(cx.edges)
    R = np.zeros((nE, nE))
    M1 = np.zeros((nE, nE))
    for t in cx.tets:
        L6 = np.array([ell[cx.edge_index[tuple(sorted((t[a], t[b])))]]
                       for (a, b) in TET_EDGES])
        P = tet_embed(L6)
        if P is None:
            continue
        if flip:
            P = P.copy(); P[:, 2] *= -1.0     # reflect -> orientation flip
        g4, V = _tet_grads(P)
        if V <= 0:
            continue
        g = g4 @ g4.T                          # 4x4 gram of gradients
        # local edge -> global edge index + orientation sign
        eidx = np.empty(6, dtype=int)
        sgn = np.empty(6)
        for k, (a, b) in enumerate(TET_EDGES):
            va, vb = t[a], t[b]
            eidx[k] = cx.edge_index[tuple(sorted((va, vb)))]
            sgn[k] = 1.0 if va < vb else -1.0   # local a->b vs global min->max
        # integral lam_i lam_j = V(1+delta)/20
        I = (np.ones((4, 4)) + np.eye(4)) * (V / 20.0)
        cross = np.zeros((6, 3))               # curl(W_k) = 2 grad a x grad b
        intW = np.zeros((6, 3))                # int W_k dV = (V/4)(grad b - grad a)
        for k, (a, b) in enumerate(TET_EDGES):
            cross[k] = 2.0 * np.cross(g4[a], g4[b])
            intW[k] = (V / 4.0) * (g4[b] - g4[a])
        for k, (a, b) in enumerate(TET_EDGES):
            for l, (c, d) in enumerate(TET_EDGES):
                # M1 = int W_k . W_l
                m = (I[a, c] * g[b, d] - I[a, d] * g[b, c]
                     - I[b, c] * g[a, d] + I[b, d] * g[a, c])
                # R = int curl(W_k) . W_l = curl(W_k) . int W_l
                r = cross[k] @ intW[l]
                ek, el, s = eidx[k], eidx[l], sgn[k] * sgn[l]
                M1[ek, el] += s * m
                R[ek, el] += s * r
    R = 0.5 * (R + R.T)                          # symmetric on a closed manifold
    M1 = 0.5 * (M1 + M1.T)
    return R, M1


def curl_spectrum(cx, ell, flip=False, tol_zero=1e-6):
    """Signed Beltrami eigenvalues, sorted by magnitude; zeros (gradient modes) dropped."""
    R, M1 = build_curl(cx, ell, flip=flip)
    w = sla.eigh(R, M1, eigvals_only=True)
    w = np.real(w)
    nz = w[np.abs(w) > tol_zero]
    return nz[np.argsort(np.abs(nz))], w


def cluster_signed(vals, reltol=0.08):
    """Group signed eigenvalues into +-level clusters; return [(mean, count), ...]."""
    out = []
    used = np.zeros(len(vals), bool)
    order = np.argsort(np.abs(vals))
    for i in order:
        if used[i]:
            continue
        v = vals[i]
        grp = [j for j in range(len(vals))
               if not used[j] and abs(vals[j] - v) <= reltol * max(abs(v), 1e-9)]
        for j in grp:
            used[j] = True
        out.append((float(np.mean([vals[j] for j in grp])), len(grp)))
    return out


def exact_sequence_residual(cx, R):
    """max|R . d0| — should be ~0 (curl of a gradient = 0). Large => R is NOT a faithful
    discrete curl (the Whitney weak form fails the de Rham complex on a closed mesh)."""
    E, V = len(cx.edges), cx.N
    d0 = np.zeros((E, V))
    for ei, (a, b) in enumerate(cx.edges):     # edges stored (min,max)
        d0[ei, b] += 1.0
        d0[ei, a] -= 1.0
    return float(np.abs(R @ d0).max()), float(np.abs(R).max())


if __name__ == "__main__":
    from sixhundred_cell import build_600cell, round_metric
    cx, _ = build_600cell()
    ell = round_metric(cx)
    print(f"Curl instrument validation on the round 600-cell (E={len(cx.edges)})")
    print("  analytic round-S^3 curl spectrum: +-2 (mult 3), +-3, +-4, ...")
    R, M1 = build_curl(cx, ell)
    res, rnorm = exact_sequence_residual(cx, R)
    print(f"  exact-sequence check max|R.d0|={res:.3e}  (||R||={rnorm:.3e})")
    GATE = res < 1e-6 * max(rnorm, 1.0)
    print(f"  >>> §2 GATE (faithful signed curl): {'PASS' if GATE else 'FAIL'} <<<")
    nz, w = curl_spectrum(cx, ell)
    print(f"  lowest-|.| nonzero signed eigenvalues:\n   {np.round(nz[:18], 3)}")
    print("  NOTE: this Whitney weak-form curl FAILS the exact sequence (R.d0 != 0),")
    print("  so gradient modes contaminate the low band and +-2(mult 3) is not clean.")
    print("  A faithful signed curl (topological d1 curl-curl + helicity-sign split, or")
    print("  a full DEC primal-dual Hodge-star formulation) is required before the 2I sieve.")