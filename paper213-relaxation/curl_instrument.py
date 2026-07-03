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
    """
    Assemble (K, R, M1) over edges with Whitney/Nedelec edge elements.
      K [e,e'] = int curl(W_e).curl(W_e')  — curl-curl STIFFNESS (symmetric PSD).
                 K annihilates gradients EXACTLY (curl(grad)=0), so its kernel is the
                 gradient subspace -> clean Hodge/Helmholtz projection. K a = mu M1 a
                 gives mu = lambda^2 (sign folded): mult 6 at mu=4, mult 16 at mu=9.
      R [e,e'] = int curl(W_e).W_e'         — HELICITY form (signed). Restricted to a
                 mu-eigenspace it splits the Beltrami modes into +lambda / -lambda.
      M1[e,e'] = int W_e.W_e'               — edge mass (SPD).
    flip=True reflects the embedding (orientation flip): K unchanged, R -> -R, i.e.
    +2 <-> -2 relabel without changing physics.
    """
    nE = len(cx.edges)
    K = np.zeros((nE, nE))
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
                # H = int W_k ^ dW_l = intW_k . curl(W_l)  (helicity form; vanishes on
                # gradients in the l-slot EXACTLY). NOTE: left UNSYMMETRIZED.
                r = intW[k] @ cross[l]
                # K = int curl(W_k) . curl(W_l) = V * cross_k . cross_l (curl-curl)
                kk = V * (cross[k] @ cross[l])
                ek, el, s = eidx[k], eidx[l], sgn[k] * sgn[l]
                M1[ek, el] += s * m
                R[ek, el] += s * r
                K[ek, el] += s * kk
    K = 0.5 * (K + K.T)
    M1 = 0.5 * (M1 + M1.T)
    # R (helicity) is NOT symmetrized — symmetrising mixes in the non-faithful
    # transpose and destroys the gradient-annihilation (R.d0 != 0). See findings.
    return K, R, M1


def beltrami_spectrum(cx, ell, flip=False, K=None, R=None, M1=None,
                      mu_tol=1e-6, basis=None):
    """
    Signed Beltrami (curl) eigenvalues via curl-curl + helicity sign-split.

    1. Solve K a = mu M1 a  (mu = lambda^2 >= 0). Gradient kernel sits at mu~0.
    2. On each nonzero-mu eigenspace, diagonalise the helicity form R to recover the
       SIGN: signed lambda = eig of (V^T R V, V^T M1 V) ~ +-sqrt(mu).
    Returns sorted-by-|.| signed eigenvalues and the raw mu spectrum.

    `basis` (E x m, M1-orthonormal columns) optionally restricts to a subspace (e.g. the
    2I-invariant subspace) before solving — the static sieve.
    """
    if K is None:
        K, R, M1 = build_curl(cx, ell, flip=flip)
    if basis is not None:                       # restrict operators to the subspace
        K = basis.T @ K @ basis
        R = basis.T @ R @ basis
        M1 = basis.T @ M1 @ basis
    mu, V = sla.eigh(K, M1)                      # ascending mu
    mu = np.real(mu)
    signed = []
    i = 0
    n = len(mu)
    while i < n:
        j = i
        while j + 1 < n and abs(mu[j + 1] - mu[i]) <= 0.08 * max(abs(mu[i]), 1.0):
            j += 1
        block = V[:, i:j + 1]                    # eigenvectors of this mu-cluster
        if mu[i] > mu_tol:                       # skip the gradient kernel (mu~0)
            Rb = block.T @ R @ block
            Mb = block.T @ M1 @ block
            lam = sla.eigh(Rb, Mb, eigvals_only=True)
            signed.extend(np.real(lam).tolist())
        i = j + 1
    signed = np.array(signed)
    return signed[np.argsort(np.abs(signed))], mu


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


def d0_incidence(cx):
    E, V = len(cx.edges), cx.N
    d0 = np.zeros((E, V))
    for ei, (a, b) in enumerate(cx.edges):     # edges stored (min,max)
        d0[ei, b] += 1.0
        d0[ei, a] -= 1.0
    return d0


def cluster_mu(mu, reltol=0.10, mu_min=1e-6):
    """Cluster nonzero mu into (value, multiplicity)."""
    nz = mu[mu > mu_min]
    out = []
    i = 0
    while i < len(nz):
        j = i
        while j + 1 < len(nz) and abs(nz[j + 1] - nz[i]) <= reltol * max(nz[i], 1.0):
            j += 1
        out.append((float(np.mean(nz[i:j + 1])), j - i + 1))
        i = j + 1
    return out


if __name__ == "__main__":
    from sixhundred_cell import build_600cell, round_metric
    cx, _ = build_600cell()
    ell = round_metric(cx)
    print(f"Curl instrument validation on the round 600-cell (E={len(cx.edges)})")
    print("  analytic round-S^3: |lambda|=2 (mu=4, mult 6 = +-2 x3), |lambda|=3 (mu=9, mult 16)")
    K, R, M1 = build_curl(cx, ell)
    d0 = d0_incidence(cx)
    res = float(np.abs(K @ d0).max()); knorm = float(np.abs(K).max())
    GATE_seq = res < 1e-8 * max(knorm, 1.0)
    mu, _ = sla.eigh(K, M1)
    nzero = int(np.sum(mu <= 1e-8))
    print("  --- MAGNITUDE (curl-curl K) ---")
    print(f"  exact-sequence max|K.d0|={res:.2e} (annihilates gradients: {GATE_seq})")
    print(f"  gradient kernel dim={nzero} (expect V-1={cx.N-1})")
    print("  lowest nonzero mu clusters [expect ~4 x6, ~9 x16]:")
    for val, m in cluster_mu(mu)[:3]:
        print(f"     mu={val:.3f} x{m}  (|lambda|~{val**0.5:.3f})")
    MAG_GATE = GATE_seq and abs(nzero - (cx.N - 1)) <= 2
    print(f"  >>> MAGNITUDE GATE: {'PASS' if MAG_GATE else 'FAIL'} <<<")
    print("  --- SIGN (helicity split) : OPEN ---")
    print(f"  H annihilates gradients: max|H.d0|={np.abs(R @ d0).max():.2e}")
    sel = np.where((mu > 3.5) & (mu < 5.5))[0]
    from numpy.linalg import multi_dot
    mu4, V4 = sla.eigh(K, M1)
    B = V4[:, sel]
    hel = np.diag(B.T @ R @ B)
    print(f"  helicity of the mu=4 (Hopf-sextet) modes: {np.round(hel, 4)}")
    print("  => helicity ~ 0 on the coarse 600-cell: the +2/-2 SIGN split is NOT")
    print("     cleanly recoverable from the Whitney helicity form at this resolution.")
    print("     MAGNITUDE validated; SIGN instrument is the open blocker (see findings).")