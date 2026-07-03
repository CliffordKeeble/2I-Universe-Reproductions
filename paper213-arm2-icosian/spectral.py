"""
spectral.py  —  spectral dimension d_s and the <12,20,30> signature score.

============================  PRE-REGISTERED (Pattern 75)  ====================
These constants define the scoring BEFORE the first production run. Do not tune
them to the result.

  M_CLUSTERS      = 8       # number of lowest nonzero L_c clusters scored
  CLUSTER_RELGAP  = 0.03    # new cluster when relative gap to previous > 3%
  LANDING_TOL     = 0.05    # cluster centre within +-5% of a target k(k+2)
  ZERO_THRESH     = 1e-6    # eigenvalues below this are treated as zero modes
  DS_FIT_WINDOW   = (0.05, 1.0)  # intermediate-t window (in units of 1/lambda_med)

Target (Paper 174 / Paper 108): even-sector of S^3/2I scalar Laplacian.
  supported k : numerical semigroup S(12,20,30); lambda_k = k(k+2); k=12 -> 168.
  killed even k (the "missing modes"): {2,4,6,8,10,14,16,18,22,26,28,34,38,46,58}.

SIGNATURE SCORE (pre-registered):
  rescale L_c spectrum so the lowest nonzero cluster centre = 168, then
    frac_supported = (#of lowest M clusters within TOL of some supported k(k+2)) / M
    frac_forbidden = (#of those M clusters within TOL of a killed   k(k+2)) / M
    signature_score = frac_supported - frac_forbidden     in [-1, 1]
  1.0 = perfect reconstruction of the even-sector ladder; ~0 = no signal.
==============================================================================
"""

import numpy as np
from scipy.sparse.linalg import eigsh

# -------- pre-registered constants --------
M_CLUSTERS = 8
CLUSTER_RELGAP = 0.03
LANDING_TOL = 0.05
ZERO_THRESH = 1e-6
DS_FIT_WINDOW = (0.05, 1.0)

# -------- continuum target ladders --------
KILLED_EVEN_K = [2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 28, 34, 38, 46, 58]


def semigroup_supported_k(kmax=120):
    """Even k in S(12,20,30) for 1 <= k <= kmax (k=0 excluded)."""
    gens = (12, 20, 30)
    sup = set()
    for a in range(kmax // 12 + 1):
        for b in range(kmax // 20 + 1):
            for c in range(kmax // 30 + 1):
                k = 12 * a + 20 * b + 30 * c
                if 1 <= k <= kmax:
                    sup.add(k)
    return sorted(sup)


SUPPORTED_K = semigroup_supported_k(120)
SUPPORTED_LADDER = np.array([k * (k + 2) for k in SUPPORTED_K], dtype=float)
KILLED_LADDER = np.array([k * (k + 2) for k in KILLED_EVEN_K], dtype=float)


# ---------------------------------------------------------------------------
# low spectrum
# ---------------------------------------------------------------------------
def low_eigenvalues(L, k=1000):
    """Lowest k eigenvalues via shift-invert eigsh (robust to zero modes)."""
    k = min(k, L.shape[0] - 2)
    for sigma in (0.0, -1e-6, -1e-4):
        try:
            vals = eigsh(L, k=k, sigma=sigma, which="LM",
                         return_eigenvectors=False, maxiter=5000)
            return np.sort(vals.real)
        except Exception:
            continue
    # dense fallback for small matrices
    vals = np.linalg.eigvalsh(L.toarray())
    return np.sort(vals.real)[:k]


# ---------------------------------------------------------------------------
# spectral dimension via heat trace
# ---------------------------------------------------------------------------
def heat_trace_dimension(eigs_Ls, N):
    """
    P(t) = (1/N) sum exp(-t lambda) over the computed low spectrum.
    Fit log P ~ -(d_s/2) log t over the intermediate-t window.
    Returns (d_s, d_s_err, (t_grid, P)).
    """
    eigs = np.asarray(eigs_Ls)
    eigs = eigs[eigs > ZERO_THRESH]
    if len(eigs) < 10:
        return np.nan, np.nan, (np.array([]), np.array([]))
    lam_med = np.median(eigs)
    t0, t1 = DS_FIT_WINDOW
    t_grid = np.geomspace(t0 / lam_med, t1 / lam_med, 40)
    P = np.array([np.exp(-t * eigs).sum() / N for t in t_grid])
    good = P > (1.5 / N)  # stay above the saturation floor of the truncated trace
    lg_t, lg_P = np.log(t_grid[good]), np.log(P[good])
    if len(lg_t) < 5:
        return np.nan, np.nan, (t_grid, P)
    A = np.vstack([lg_t, np.ones_like(lg_t)]).T
    coef, res, *_ = np.linalg.lstsq(A, lg_P, rcond=None)
    slope = coef[0]
    # crude slope uncertainty from residual
    resid = lg_P - A @ coef
    dof = max(1, len(lg_t) - 2)
    serr = np.sqrt((resid @ resid) / dof) / (np.std(lg_t) * np.sqrt(len(lg_t)))
    return -2.0 * slope, 2.0 * serr, (t_grid, P)


# ---------------------------------------------------------------------------
# <12,20,30> signature score
# ---------------------------------------------------------------------------
def cluster_spectrum(eigs):
    """Group sorted eigenvalues into clusters by relative gap; return centres."""
    eigs = np.asarray(eigs)
    eigs = eigs[eigs > ZERO_THRESH]
    if len(eigs) == 0:
        return np.array([])
    centres, cur = [], [eigs[0]]
    for x in eigs[1:]:
        if (x - cur[-1]) / max(x, ZERO_THRESH) > CLUSTER_RELGAP:
            centres.append(np.mean(cur)); cur = [x]
        else:
            cur.append(x)
    centres.append(np.mean(cur))
    return np.array(centres)


def signature_score(eigs_Lc, M=M_CLUSTERS, tol=LANDING_TOL):
    """
    Score the L_c low spectrum against the even-sector k(k+2) ladder.
    Returns a dict (pre-registered metrics).
    """
    centres = cluster_spectrum(eigs_Lc)
    if len(centres) < 2:
        return dict(signature_score=np.nan, frac_supported=np.nan,
                    frac_forbidden=np.nan, lambda1_raw=np.nan,
                    ratio_obs=np.nan, ratio_target=440.0 / 168.0,
                    n_clusters=len(centres))
    lam1 = centres[0]
    scaled = centres * (168.0 / lam1)          # lowest nonzero cluster -> 168
    test = scaled[:M]

    def hits(ladder):
        h = 0
        for x in test:
            rel = np.min(np.abs(ladder - x) / ladder)
            if rel <= tol:
                h += 1
        return h

    n = len(test)
    frac_sup = hits(SUPPORTED_LADDER) / n
    frac_for = hits(KILLED_LADDER) / n
    ratio_obs = (centres[1] / centres[0]) if len(centres) > 1 else np.nan
    return dict(
        signature_score=frac_sup - frac_for,
        frac_supported=frac_sup,
        frac_forbidden=frac_for,
        lambda1_raw=float(lam1),
        ratio_obs=float(ratio_obs),
        ratio_target=440.0 / 168.0,
        n_clusters=int(len(centres)),
        scaled_low=test.tolist(),
    )


if __name__ == "__main__":
    print("supported even k (<=120):", SUPPORTED_K)
    print("killed even k          :", KILLED_EVEN_K)
    print("target ladder (low)    :", SUPPORTED_LADDER[:6].astype(int))
    print("target ratio l2/l1     :", round(440.0 / 168.0, 4))
