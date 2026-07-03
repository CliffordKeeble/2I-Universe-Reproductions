"""
fss_stage1.py  —  the decisive Stage-1 ladder (Paper 213 FSS brief).

For each N in the ladder: relax the frustrated metric to a converged local min (O(N)
local-FD gradient), build the same-mesh round reference, and read the PRIMARY observable

    d_frust(N) = spectral_distance(relaxed frustrated endpoint, same-mesh round)

with the secondary variance gap as a cross-check. Then fit  d_frust(N)=d_inf+c N^-p:
  d_inf ~ 0  -> ARTEFACT (round reachable at fine N; headline falls)
  d_inf > 0  -> FUNDAMENTAL glass (round not a generic attractor; headline survives).

Controls per N: Control #3 (round-init stays round) and reference sanity d_round(N)
(round-init's residual distance to the analytic k(k+2) ladder, should fall toward 0).

Pre-registration committed before this run (RESULTS.md). No tuning of fit/recipe.

Standalone:  python fss_stage1.py "100,200,500,1000"
"""

import sys
import json
import numpy as np

from geometry import generic_s3
from relax_fast import precompute_fss, relax_fss, base_geometry
from spectral_s3 import low_spectrum, normalized_spectrum, spectrum_match_error
from fast_geometry import precompute as precompute_curv, fast_energy


def round_ell(cx):
    P = cx.hull_points
    return np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])


def spectral_distance(eigsA, eigsB, M=12):
    a = normalized_spectrum(eigsA); b = normalized_spectrum(eigsB)
    m = min(M, len(a), len(b))
    if m < 3:
        return np.nan
    a, b = a[:m], b[:m]
    return float(np.sqrt(np.mean(((a - b) / b) ** 2)))


def cov_of(ell, pre_curv):
    _, cov, _, _ = fast_energy(ell, pre_curv)
    return cov


def run_N(N, seed=20260618, k_eig=40, confirm=False):
    cx = generic_s3(N, seed=seed)
    pre = precompute_fss(cx)
    pcv = precompute_curv(cx)
    ell1 = np.ones(len(cx.edges))
    ellr = round_ell(cx)

    relaxed, info = relax_fss(cx, ell1, pre=pre, iters=300)
    cov_frust = cov_of(relaxed, pcv)
    cov_round = cov_of(ellr, pcv)

    # Control #3: round-init must stay round
    _, info_c3 = relax_fss(cx, ellr, pre=pre, iters=300)
    stays_round = bool(info_c3["Ef"] <= info_c3["E0"] * 1.5)

    spec_init, _ = low_spectrum(cx, ell1, k=k_eig)
    spec_relaxed, _ = low_spectrum(cx, relaxed, k=k_eig)
    spec_round, _ = low_spectrum(cx, ellr, k=k_eig)

    rec = dict(
        N=N, edges=len(cx.edges), tets=len(cx.tets),
        E0=info["E0"], Ef=info["Ef"], relax_iters=info["iters"], relax_msg=info["msg"],
        d_frust_init=spectral_distance(spec_init, spec_round),
        d_frust=spectral_distance(spec_relaxed, spec_round),     # PRIMARY observable
        d_round=spectrum_match_error(spec_round, 4),             # reference sanity
        cov_frust=float(cov_frust), cov_round=float(cov_round),
        cov_gap=float(cov_frust - cov_round),
        control3_stays_round=stays_round,
        control3_E0=info_c3["E0"], control3_Ef=info_c3["Ef"],
    )
    if confirm:                                                  # identical-endpoint re-run
        relaxed2, info2 = relax_fss(cx, ell1, pre=pre, iters=300)
        sp2, _ = low_spectrum(cx, relaxed2, k=k_eig)
        rec["confirm_Ef"] = info2["Ef"]
        rec["confirm_d_frust"] = spectral_distance(sp2, spec_round)
    return rec


def fit_dinf(Ns, ds):
    """Fit d_frust = d_inf + c N^-p. Returns (d_inf, c, p, residual)."""
    from scipy.optimize import curve_fit
    Ns = np.array(Ns, float); ds = np.array(ds, float)
    good = np.isfinite(ds)
    Ns, ds = Ns[good], ds[good]
    if len(Ns) < 3:
        return None
    def model(N, dinf, c, p):
        return dinf + c * N ** (-p)
    try:
        popt, pcov = curve_fit(model, Ns, ds, p0=[0.1, 1.0, 0.5],
                               bounds=([-0.5, -50, 0.05], [1.0, 50, 4.0]), maxfev=20000)
        resid = ds - model(Ns, *popt)
        perr = np.sqrt(np.diag(pcov))
        return dict(d_inf=float(popt[0]), d_inf_err=float(perr[0]),
                    c=float(popt[1]), p=float(popt[2]),
                    rms_resid=float(np.sqrt(np.mean(resid ** 2))))
    except Exception as e:
        return dict(error=str(e))


def main():
    ladder = ([int(x) for x in sys.argv[1].split(",")]
              if len(sys.argv) > 1 else [100, 200, 500, 1000])
    print(f"=== Stage 1 FSS ladder: N = {ladder} ===")
    print("    PRIMARY: d_frust(N) = spectral_distance(relaxed frustrated, same-mesh round)")
    results = []
    for i, N in enumerate(ladder):
        confirm = (i == len(ladder) - 1)            # identical-endpoint re-run at largest N
        r = run_N(N, confirm=confirm)
        results.append(r)
        print(f"\nN={N:5d} (E:{r['E0']:.0f}->{r['Ef']:.1f}, iters={r['relax_iters']})")
        print(f"   d_frust_init={r['d_frust_init']:.3f} -> d_frust={r['d_frust']:.3f}"
              f"   d_round(ref sanity)={r['d_round']:.3f}")
        print(f"   CoV gap (frust-round)={r['cov_gap']:.2f}   Control#3 stays_round={r['control3_stays_round']}")
        if confirm:
            print(f"   convergence confirm: Ef={r.get('confirm_Ef'):.1f} "
                  f"d_frust={r.get('confirm_d_frust'):.3f} (vs {r['d_frust']:.3f})")
        json.dump(results, open("results_fss.json", "w"), indent=2)  # incremental save

    fit = fit_dinf([r["N"] for r in results], [r["d_frust"] for r in results])
    print("\n==== FSS VERDICT ====")
    print(f"  d_frust(N): {[round(r['d_frust'],3) for r in results]}")
    print(f"  d_round(N) (ref, should fall): {[round(r['d_round'],3) for r in results]}")
    if fit and "d_inf" in fit:
        print(f"  fit d_frust = d_inf + c N^-p:  d_inf={fit['d_inf']:.3f}±{fit['d_inf_err']:.3f}"
              f"  c={fit['c']:.2f}  p={fit['p']:.2f}  rms_resid={fit['rms_resid']:.4f}")
        if fit["d_inf"] - 2 * fit["d_inf_err"] <= 0:
            print("  => d_inf consistent with 0: ARTEFACT (round reachable at fine N).")
        else:
            print("  => d_inf bounded away from 0: FUNDAMENTAL glass (headline survives).")
    json.dump(dict(ladder=results, fit=fit), open("results_fss.json", "w"), indent=2)
    print("  wrote results_fss.json")


if __name__ == "__main__":
    main()
