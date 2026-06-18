"""
spectral_s3.py  —  the VERDICT instrument for Paper 213 Approach B v1.1.

Per Mr A #2 (via brief v1.1): variance drives the flow, but the SPECTRUM is the
judgment. Convergence to round S^3 = the relaxed scalar-Laplacian spectrum matches the
analytic round-S^3 spectrum

        lambda_k = k(k+2),   multiplicity (k+1)^2,   k = 0,1,2,...
        => 0 (x1), 3 (x4), 8 (x9), 15 (x16), 24 (x25), ...

We build the linear-FEM Laplace-Beltrami operator on the tetrahedral mesh purely from
edge lengths (each tet embedded via geometry.tet_embed; hat-function gradients give the
element stiffness; lumped vertex volumes give the mass matrix), then solve the
generalized eigenproblem L u = lambda M u.

VALIDATION (run before trusting it): on the round-init metric (chord lengths on the
actual S^3 hull points) the low spectrum must approximate 0,3,8,15 with multiplicity
pattern 1,4,9,16. The round-init spectrum on a given mesh is ALSO the empirical "round"
reference for that mesh (removes discretisation bias from the verdict).

Standalone:  python spectral_s3.py [N]
"""

import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from geometry import generic_s3, tet_embed, TET_EDGES


def fem_laplacian(cx, ell):
    """Linear-FEM Laplace-Beltrami (stiffness L, lumped mass M) from edge lengths."""
    N = cx.N
    rows, cols, data = [], [], []
    mass = np.zeros(N)
    n_bad = 0
    for t in cx.tets:
        L6 = np.array([ell[cx.edge_index[tuple(sorted((t[a], t[b])))]]
                       for (a, b) in TET_EDGES])
        P = tet_embed(L6)
        if P is None:
            n_bad += 1
            continue
        e = P[1:] - P[0]                      # 3x3: edge vectors from vertex 0
        detM = np.linalg.det(e)
        V = abs(detM) / 6.0
        if V < 1e-12:
            n_bad += 1
            continue
        Minv = np.linalg.inv(e)               # rows = grad lambda_1,2,3
        grads = np.zeros((4, 3))
        grads[1:] = Minv                       # grad lambda_k, k=1,2,3
        grads[0] = -grads[1:].sum(axis=0)      # grad lambda_0
        Ke = V * (grads @ grads.T)             # 4x4 element stiffness
        for ia in range(4):
            mass[t[ia]] += V / 4.0
            for ib in range(4):
                rows.append(t[ia]); cols.append(t[ib]); data.append(Ke[ia, ib])
    Lmat = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    Lmat = 0.5 * (Lmat + Lmat.T)
    Mmat = sp.diags(mass)
    return Lmat, Mmat, n_bad


def low_spectrum(cx, ell, k=40):
    L, M, n_bad = fem_laplacian(cx, ell)
    k = min(k, cx.N - 2)
    # generalized eigenproblem L u = lam M u, smallest eigenvalues
    for sigma in (-1e-6, -1e-4, 0.0):
        try:
            vals = eigsh(L, k=k, M=M, sigma=sigma, which="LM",
                         return_eigenvectors=False, maxiter=5000)
            return np.sort(vals.real), n_bad
        except Exception:
            continue
    Md = np.array(M.diagonal())
    Dinv = sp.diags(1.0 / np.sqrt(Md))
    vals = np.linalg.eigvalsh((Dinv @ L @ Dinv).toarray())
    return np.sort(vals.real)[:k], n_bad


def normalized_spectrum(eigs):
    """Drop the zero mode; rescale so the first nonzero eigenvalue = 3 (the k=1 level)."""
    e = np.array(eigs)
    e = e[e > 1e-6]
    if len(e) == 0:
        return e
    return e * (3.0 / e[0])


def round_ladder(n_levels=6):
    return np.array([k * (k + 2) for k in range(1, n_levels + 1)], dtype=float)


def spectrum_match_error(eigs, n_levels=4):
    """
    Compare the normalized low spectrum to the round ladder k(k+2). For each of the
    first n_levels round levels, take the mean of the eigenvalues assigned to it by
    nearest-level rounding; report RMS relative error. Lower = closer to round.
    """
    e = normalized_spectrum(eigs)
    ladder = round_ladder(n_levels)
    if len(e) < n_levels + 1:
        return np.nan
    # assign each eigenvalue (up to the n_levels'th level window) to nearest ladder level
    use = e[e <= ladder[-1] * 1.3]
    if len(use) == 0:
        return np.nan
    errs = []
    for lv in ladder:
        near = use[np.abs(use - lv) / lv < 0.5]
        if len(near):
            errs.append((np.mean(near) - lv) / lv)
        else:
            errs.append(1.0)  # a missing level is a full miss
    return float(np.sqrt(np.mean(np.array(errs) ** 2)))


def multiplicity_report(eigs, n_levels=4, tol=0.25):
    """Count eigenvalues within tol of each round level (expect 4,9,16,25)."""
    e = normalized_spectrum(eigs)
    ladder = round_ladder(n_levels)
    expected = [(k + 1) ** 2 for k in range(1, n_levels + 1)]
    counts = [int(np.sum(np.abs(e - lv) / lv < tol)) for lv in ladder]
    return list(zip([int(l) for l in ladder], counts, expected))


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    cx = generic_s3(N, seed=20260618)
    P = cx.hull_points
    ell_round = np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])
    ell_frust = np.ones(len(cx.edges))

    print(f"FEM Laplace-Beltrami spectrum validation (N={N}, tets={len(cx.tets)})")
    print("  analytic round S^3 ladder: 3(x4) 8(x9) 15(x16) 24(x25)")
    for label, ell in (("round-init", ell_round), ("frustrated ell=1", ell_frust)):
        eigs, nbad = low_spectrum(cx, ell, k=40)
        ns = normalized_spectrum(eigs)
        err = spectrum_match_error(eigs, 4)
        mult = multiplicity_report(eigs, 4)
        print(f"\n[{label}] bad_tets={nbad}")
        print(f"  low normalized spectrum: {np.round(ns[:16], 2)}")
        print(f"  spectrum match error (vs k(k+2), 4 levels): {err:.3f}")
        print(f"  multiplicities [level: got/expected]: "
              f"{', '.join(f'{l}:{g}/{e}' for l,g,e in mult)}")
