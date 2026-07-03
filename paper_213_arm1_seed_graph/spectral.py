"""
spectral.py -- Paper 213 Arm 1: spectral dimension d_s and the <12,20,30>
signature score (brief sec.4).

PRE-REGISTERED CONSTANTS (Pattern 75 -- fixed BEFORE first production run):
  M_CLUSTERS = 10      number of lowest clusters scored
  TAU        = 0.05    relative tolerance for a cluster "landing on" a target
  REL_GAP    = 0.15    new cluster when consecutive eigenvalue gap > 15% of value
  ZERO_REL   = 1e-6    eigenvalues below ZERO_REL * lambda_max treated as zero
  JITTER     = 1e-9    diagonal regulariser for shift-invert on singular L_s

Targets:
  k(k+2) integer test : T_k = k(k+2) for k = 1..M = [3,8,15,24,35,48,63,80,99,120]
  <12,20,30> semigroup: gaps = {2,4,6,8,10,14,16,18,22,26,28,34,38,46,58}
                        (= 2 * gaps of <6,10,15>, Frobenius number 58)
  spectral gap (Paper 108): lambda_1 = 168 on S^3/2I (NOT expected on a growing
                        graph, whose gap -> 0; reported as informative null).
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

M_CLUSTERS = 10
TAU = 0.05
REL_GAP = 0.15
ZERO_REL = 1e-6
JITTER = 1e-9

KK2_TARGETS = [k * (k + 2) for k in range(1, M_CLUSTERS + 1)]
SEMIGROUP_GAPS = [2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 28, 34, 38, 46, 58]


# ---------------------------------------------------------------------------
# Low spectrum
# ---------------------------------------------------------------------------

def low_spectrum(L, k=1000):
    """k eigenvalues of smallest magnitude (the low spectrum), via shift-invert
    sigma=0 with a tiny diagonal jitter so factorisation is stable on a singular
    L_s. Returns eigenvalues ascending. Falls back to dense for small L."""
    n = L.shape[0]
    k = min(k, n - 2)
    if n <= 400:
        w = np.linalg.eigvalsh(L.toarray())
        return np.sort(w)
    Lj = L + JITTER * sp.identity(n, format="csr")
    try:
        w = spla.eigsh(Lj, k=k, sigma=0.0, which="LM",
                       return_eigenvectors=False, maxiter=5000)
        w = np.sort(w) - JITTER
    except Exception:
        # robust fallback: smallest algebraic
        w = spla.eigsh(L, k=k, which="SA",
                       return_eigenvectors=False, maxiter=10000)
        w = np.sort(w)
    return w


# ---------------------------------------------------------------------------
# Spectral dimension via heat-kernel return probability
# ---------------------------------------------------------------------------

def spectral_dimension(eigs, N, n_t=60):
    """d_s from P(t) = (1/N) sum_i exp(-t lambda_i) over the LOW spectrum.

    Fit log P ~ -(d_s/2) log t over the intermediate window
    t in [a/lambda_k, b/lambda_1] where the truncated low spectrum is a faithful
    sample (t large enough that the missing high modes have decayed, small enough
    to sit below the spectral-gap plateau). Returns (d_s, stderr, window)."""
    lam = np.asarray(eigs)
    lam = lam[lam > ZERO_REL * lam.max()]          # drop zero modes
    if lam.size < 10:
        return np.nan, np.nan, (np.nan, np.nan)
    # Robust window from percentiles of the LOW spectrum: stay in the scaling
    # bulk, away from the highest computed mode (small-t truncation error) and
    # away from the vanishing spectral gap (large-t plateau). Using min/max would
    # let lambda_1 -> 0 (as N grows) drag t_hi into the plateau.
    lam_hi = np.percentile(lam, 90)
    lam_lo = np.percentile(lam, 5)
    t_lo = 5.0 / lam_hi             # high modes decayed
    t_hi = 0.5 / lam_lo             # below the gap plateau
    if not (t_hi > t_lo):
        return np.nan, np.nan, (t_lo, t_hi)
    ts = np.logspace(np.log10(t_lo), np.log10(t_hi), n_t)
    P = np.array([np.exp(-t * lam).sum() / N for t in ts])
    x, y = np.log(ts), np.log(P)
    # least squares slope with standard error
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope = coef[0]
    yhat = A @ coef
    dof = max(1, len(x) - 2)
    s2 = np.sum((y - yhat) ** 2) / dof
    se_slope = np.sqrt(s2 / np.sum((x - x.mean()) ** 2))
    d_s = -2.0 * slope
    d_s_se = 2.0 * se_slope
    return d_s, d_s_se, (t_lo, t_hi)


# ---------------------------------------------------------------------------
# <12,20,30> signature score
# ---------------------------------------------------------------------------

def cluster_eigs(eigs, rel_gap=REL_GAP):
    """Cluster positive eigenvalues; split when consecutive gap > rel_gap*value.
    Returns list of cluster centres (means), ascending."""
    lam = np.asarray(eigs)
    lam = np.sort(lam[lam > ZERO_REL * lam.max()])
    if lam.size == 0:
        return []
    clusters, cur = [], [lam[0]]
    for x in lam[1:]:
        if x - cur[-1] > rel_gap * max(cur[-1], 1e-12):
            clusters.append(np.mean(cur)); cur = [x]
        else:
            cur.append(x)
    clusters.append(np.mean(cur))
    return clusters


def _frac_landing(values, targets, tau=TAU):
    """Fraction of `values` within relative tau of SOME target."""
    if not len(values):
        return 0.0
    hits = 0
    for v in values:
        if any(abs(v - t) <= tau * t for t in targets):
            hits += 1
    return hits / len(values)


def signature_score(eigs):
    """Returns a dict of signature diagnostics (brief sec.4 i/ii/iii).

      kk2_score      : fraction of M lowest clusters landing on k(k+2)
                       (lowest nonzero cluster rescaled to 3 = k=1 target)
      semigroup_hit  : fraction of clusters landing on a <12,20,30> element
                       (lowest cluster rescaled to 12)
      gap_avoidance  : fraction of the 15 semigroup gaps left EMPTY
      lambda1        : smallest positive eigenvalue (raw)
      gap_ratio      : lambda_2 / lambda_1  (cluster centres)
      n_negative     : count of eigenvalues < -ZERO_REL*max (non-PSD fingerprint)
    """
    lam = np.asarray(eigs)
    lam_max = lam.max() if lam.size else 1.0
    n_neg = int(np.sum(lam < -ZERO_REL * lam_max))
    clusters = cluster_eigs(eigs)
    out = {"n_clusters": len(clusters), "n_negative": n_neg,
           "lambda1": np.nan, "gap_ratio": np.nan,
           "kk2_score": 0.0, "semigroup_hit": 0.0, "gap_avoidance": 1.0}
    if len(clusters) < 2:
        return out
    c = np.array(clusters)
    out["lambda1"] = float(c[0])
    out["gap_ratio"] = float(c[1] / c[0])

    # (i) k(k+2): rescale lowest cluster -> 3
    c_kk2 = c / c[0] * 3.0
    out["kk2_score"] = _frac_landing(c_kk2[:M_CLUSTERS], KK2_TARGETS)

    # (ii)/(iii) semigroup: rescale lowest cluster -> 12
    c_sg = c / c[0] * 12.0
    # semigroup elements up to max gap+generators
    sg_elems = _numerical_semigroup([12, 20, 30], upto=130)
    out["semigroup_hit"] = _frac_landing(c_sg[:M_CLUSTERS], sg_elems)
    empty = sum(1 for g in SEMIGROUP_GAPS
                if not any(abs(v - g) <= TAU * g for v in c_sg))
    out["gap_avoidance"] = empty / len(SEMIGROUP_GAPS)
    return out


def _numerical_semigroup(gens, upto):
    reach = {0}
    for x in range(1, upto + 1):
        if any((x - g) in reach for g in gens):
            reach.add(x)
    return sorted(reach)


if __name__ == "__main__":
    from growth import grow
    from laplacian import scalar_laplacian, connection_laplacian
    G, meta, em = grow(2000, d_max=4, C="flat-triangles")
    Ls, _ = scalar_laplacian(G)
    Lc, _ = connection_laplacian(G, em)
    es = low_spectrum(Ls, k=600)
    ec = low_spectrum(Lc, k=600)
    ds, se, win = spectral_dimension(es, meta["N"])
    print(f"N={meta['N']}  d_s = {ds:.3f} +/- {se:.3f}  window t in "
          f"[{win[0]:.2e},{win[1]:.2e}]")
    print("L_s signature:", {k: (round(v, 3) if isinstance(v, float) else v)
                             for k, v in signature_score(es).items()})
    print("L_c signature:", {k: (round(v, 3) if isinstance(v, float) else v)
                             for k, v in signature_score(ec).items()})
