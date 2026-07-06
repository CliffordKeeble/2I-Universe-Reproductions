"""
Paper 116 spectral-statistics pipeline (Part 1 + Part 2 shared core).

    Z_Gamma(t) = sum_{l>=1, m(l)>0} m(l) * lambda_l^{-1/2} * cos(t * log lambda_l),
    lambda_l = l(l+2).

l = 0 (lambda = 0, log lambda undefined) is excluded from the sum, matching the
paper's "Active l" convention. Zeros by sign-change bracketing + Brent. Unfolding
by a degree-d polynomial fit to the zero-counting staircase; spacings normalised
to unit mean.

Class discriminant (A1/Q1):
  primary   = CDF-MSD (Cramer-von Mises style, binning-free): mean squared
              difference between the empirical CDF of unit-mean spacings and each
              reference CDF, at the sorted spacings. ratio = MSD_GUE/MSD_GOE.
  secondary = binned-PDF MSD (empirical density vs surmise PDF on s in [0,3]).

Reference moments (A1/F1): Var GOE = 4/pi - 1 = 0.273, GUE = 3pi/8 - 1 = 0.178,
Poisson = 1 -- NOT 116 s2.3's printed 0.405/0.286 (on the erratum list).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# --------------------------------------------------------------------------
# Reference surmises (normalised to unit mean).
# --------------------------------------------------------------------------
def poisson_pdf(s):
    return np.exp(-s)


def goe_pdf(s):
    return (np.pi / 2.0) * s * np.exp(-np.pi * s * s / 4.0)


def gue_pdf(s):
    return (32.0 / np.pi**2) * s * s * np.exp(-4.0 * s * s / np.pi)


REF_PDF = {"Poisson": poisson_pdf, "GOE": goe_pdf, "GUE": gue_pdf}
REF_VAR = {"Poisson": 1.0, "GOE": 4.0 / np.pi - 1.0, "GUE": 3.0 * np.pi / 8.0 - 1.0}

# Reference CDFs by trapezoid integration of the PDF on a fine grid, once.
_SGRID = np.linspace(0.0, 12.0, 48001)


def _cdf_table(pdf):
    y = pdf(_SGRID)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(_SGRID))])
    return c / c[-1]


_CDF = {k: _cdf_table(p) for k, p in REF_PDF.items()}


def ref_cdf(name, s):
    return np.interp(s, _SGRID, _CDF[name])


# --------------------------------------------------------------------------
# Modes and the series Z(t).
# --------------------------------------------------------------------------
def build_modes(m, l_max=None):
    """Return (amp, loglam, ls) for the active modes l>=1 with m(l)>0.

    amp = m(l) * lambda_l^{-1/2},  loglam = log lambda_l,  lambda_l = l(l+2).
    """
    if l_max is None:
        l_max = len(m) - 1
    ls = [l for l in range(1, l_max + 1) if m[l] > 0]
    lam = np.array([l * (l + 2) for l in ls], dtype=float)
    mult = np.array([m[l] for l in ls], dtype=float)
    amp = mult / np.sqrt(lam)
    return amp, np.log(lam), ls


def z_vec(t, amp, loglam):
    """Z(t) over an array t, accumulated mode-by-mode (memory-light)."""
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    for a, w in zip(amp, loglam):
        out += a * np.cos(w * t)
    return out


def rice_density(amp, loglam):
    """Expected zero-crossing rate (Rice): (1/pi) sqrt(<w^2> weighted by amp^2)."""
    w2 = np.sum(amp**2 * loglam**2) / np.sum(amp**2)
    return np.sqrt(w2) / np.pi


# --------------------------------------------------------------------------
# Zeros.
# --------------------------------------------------------------------------
def find_zeros(amp, loglam, t_max, t_min=1e-6, samples_per_period=40):
    """All sign-change zeros of Z on (t_min, t_max), refined by Brent."""
    fmax = float(loglam.max())
    dt = (2.0 * np.pi / fmax) / samples_per_period
    n = int(np.ceil((t_max - t_min) / dt)) + 1
    tg = np.linspace(t_min, t_max, n)
    zg = z_vec(tg, amp, loglam)
    sgn = np.sign(zg)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]

    def f(tt):
        return float(np.sum(amp * np.cos(tt * loglam)))

    zeros = []
    for i in idx:
        a, b = tg[i], tg[i + 1]
        try:
            zeros.append(brentq(f, a, b, xtol=1e-11, rtol=1e-13))
        except ValueError:
            pass
    return np.array(zeros)


def find_zeros_target(amp, loglam, target_N, t_min=1e-6, buffer=1.25,
                      samples_per_period=40, max_grow=8):
    """Grow t_max until at least target_N zeros are bracketed, then take the
    first target_N. Returns (zeros[:target_N], t_max_effective)."""
    rho = rice_density(amp, loglam)
    t_max = buffer * target_N / max(rho, 1e-9)
    for _ in range(max_grow):
        z = find_zeros(amp, loglam, t_max, t_min, samples_per_period)
        if len(z) >= target_N:
            z = z[:target_N]
            return z, float(z[-1])
        t_max *= 1.6
    return z, (float(z[-1]) if len(z) else t_max)  # short group -- flagged upstream


# --------------------------------------------------------------------------
# Unfolding and spacings.
# --------------------------------------------------------------------------
def unfold(zeros, degree=3):
    """Degree-d polynomial fit to the counting staircase; unit-mean spacings."""
    z = np.sort(zeros)
    idx = np.arange(1, len(z) + 1, dtype=float)
    coef = np.polyfit(z, idx, degree)
    x = np.polyval(coef, z)
    s = np.diff(x)
    return s / np.mean(s)


# --------------------------------------------------------------------------
# Statistics and classification.
# --------------------------------------------------------------------------
def cvm_msd(spac, name):
    """Primary discriminant: CDF-MSD against reference `name`."""
    s = np.sort(spac)
    n = len(s)
    emp = (np.arange(1, n + 1) - 0.5) / n
    return float(np.mean((ref_cdf(name, s) - emp) ** 2))


def binned_msd(spac, name, nbins, srange=(0.0, 3.0)):
    """Secondary discriminant: histogram density vs surmise PDF at bin centres."""
    hist, edges = np.histogram(spac, bins=nbins, range=srange, density=True)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    return float(np.mean((hist - REF_PDF[name](ctr)) ** 2))


def classify_cvm(spac):
    """Return (class, ratio=MSD_GUE/MSD_GOE, {msd per ref}) by the primary CDF-MSD."""
    msd = {k: cvm_msd(spac, k) for k in ("Poisson", "GOE", "GUE")}
    cls = min(msd, key=msd.get)
    ratio = msd["GUE"] / msd["GOE"]
    return cls, ratio, msd


def analyse_spacings(spac):
    cls, ratio, msd = classify_cvm(spac)
    return {
        "n_spacings": len(spac),
        "var": float(np.var(spac)),
        "mean": float(np.mean(spac)),
        "class_cvm": cls,
        "ratio_cvm": ratio,
        "cvm_msd": msd,
    }


def run_series(m, target_N=500, l_max=None, degree=3, t_max=None,
               samples_per_period=40):
    """Full pipeline for one multiplicity array m. If t_max is given, use it
    directly (equal-coverage mode); else grow to target_N zeros (equal-N mode)."""
    amp, loglam, ls = build_modes(m, l_max=l_max)
    if t_max is not None:
        zeros = find_zeros(amp, loglam, t_max, samples_per_period=samples_per_period)
        t_eff = t_max
    else:
        zeros, t_eff = find_zeros_target(amp, loglam, target_N,
                                         samples_per_period=samples_per_period)
    spac = unfold(zeros, degree=degree)
    res = analyse_spacings(spac)
    res.update({"n_active": len(ls), "l_max_active": (ls[-1] if ls else 0),
                "n_zeros": len(zeros), "t_max": t_eff})
    return res, zeros, spac


# --------------------------------------------------------------------------
# Smoke test: the switch direction must be right, or the pipeline is wrong.
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import mult

    # sanity: reference surmises have unit mean and the A1/F1 variances.
    for name in ("Poisson", "GOE", "GUE"):
        surv = 1.0 - _CDF[name]
        mean = np.sum(0.5 * (surv[1:] + surv[:-1]) * np.diff(_SGRID))  # int(1-F)ds
        print(f"ref {name:8s}: E[s]={mean:.4f}  Var_ref={REF_VAR[name]:.4f}  "
              f"F(12)={_CDF[name][-1]:.6f}")

    print("\n--- 2I (expect GUE, natural large gap) ---")
    g = mult.group_by_name("2I")
    m2I = mult.m_exact_series(*g[:5], l_max=200)
    res, _, _ = run_series(m2I, target_N=400, l_max=200)
    print({k: res[k] for k in ("n_active", "n_zeros", "var", "class_cvm", "ratio_cvm")})

    print("\n--- Z120 (expect GOE, low modes present) ---")
    mZ = mult.m_exact_series("Z120", "cyclic", 120, None, None, l_max=600)
    res, _, _ = run_series(mZ, target_N=400, l_max=600)
    print({k: res[k] for k in ("n_active", "n_zeros", "var", "class_cvm", "ratio_cvm")})
