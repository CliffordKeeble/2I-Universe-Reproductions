"""
Paper 116 Part 2 -- the 91-group generalisation (pre-registered sweep).

Pre-registration: PRE-REGISTRATION.md, commit 9227629 (Gate G1), committed before
this ran. Substrate: Part 1 merged to main (d08c6f4). Hypotheses H_gap/H_correlate
are A1's frozen text; operational parameters are pre-reg §2. Findings only, no
causal language -- correlations and z-scores.

Run:  python part2_sweep.py   -> writes part2_*.csv, prints the digest.
"""
from __future__ import annotations

import csv
import os

import numpy as np
from scipy.stats import rankdata, spearmanr

import mult
import spectral_stats as ss

DIR = os.path.dirname(os.path.abspath(__file__))
PREREG = "9227629"
SUBSTRATE = "d08c6f4"
L_MAX = 300
TARGET_N = 500
SPP = 64
SEED = 20260706
N_PERM = 10000

EXEMPLARS = ["2T", "2O", "2I", "Z12", "Z60", "Z120", "2D5", "2D15", "2D30"]
FAMILY = {"cyclic": "cyclic", "dihedral": "dihedral", "poly": "polyhedral"}


# ----------------------------------------------------------------------
# Per-group sweep.
# ----------------------------------------------------------------------
def sweep_group(name, kind, param):
    m = mult.m_exact_series(name, kind, param, None, None, L_MAX)
    l_gap = next(l for l in range(1, L_MAX + 1) if m[l] > 0)
    lam1 = l_gap * (l_gap + 2)
    _classes, order = mult.classes_order(name, kind, param)
    res, _z, _s = ss.run_series(m, target_N=TARGET_N, l_max=L_MAX,
                                samples_per_period=SPP)
    return {
        "name": name, "family": FAMILY[kind], "order": order,
        "solvable": name != "2I", "l_gap": l_gap, "lam1": lam1,
        "n_active": res["n_active"], "density": res["n_active"] / L_MAX,
        "n_zeros": res["n_zeros"], "t_max": res["t_max"],
        "var": res["var"], "ratio": res["ratio_cvm"], "class": res["class_cvm"],
        "reached_N": res["n_zeros"] >= TARGET_N,
    }


def run_sweep():
    rows = []
    for (name, kind, param, _cls, _ord, _p1) in mult.population():
        rows.append(sweep_group(name, kind, param))
    return rows


# ----------------------------------------------------------------------
# Statistics: Spearman + permutation z; rank-partial correlation.
# ----------------------------------------------------------------------
def perm_z(x, y, n_perm=N_PERM, seed=SEED):
    """Tie-corrected Spearman (Pearson on average-ranks) + permutation z, p."""
    rx, ry = rankdata(x), rankdata(y)
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return float("nan"), float("nan"), float("nan")
    rho = np.corrcoef(rx, ry)[0, 1]
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1]
                     for _ in range(n_perm)])
    z = (rho - null.mean()) / null.std()
    p = (np.sum(np.abs(null) >= abs(rho)) + 1) / (n_perm + 1)
    return rho, z, p


def rank_partial(a, b, controls):
    """Spearman-partial corr(a, b | controls): correlate rank-residuals."""
    ra, rb = rankdata(a), rankdata(b)
    C = np.column_stack([rankdata(c) for c in controls] + [np.ones(len(a))])
    ba = np.linalg.lstsq(C, ra, rcond=None)[0]
    bb = np.linalg.lstsq(C, rb, rcond=None)[0]
    resa, resb = ra - C @ ba, rb - C @ bb
    if np.std(resa) < 1e-9 or np.std(resb) < 1e-9:
        return float("nan")
    return np.corrcoef(resa, resb)[0, 1]


def outlier_percentile(rows, target="2I"):
    """Where 2I's Var sits vs the lam1->Var rank trend: percentile of its residual."""
    lam1 = np.array([r["lam1"] for r in rows], float)
    var = np.array([r["var"] for r in rows], float)
    C = np.column_stack([rankdata(lam1), np.ones(len(rows))])
    beta = np.linalg.lstsq(C, rankdata(var), rcond=None)[0]
    resid = rankdata(var) - C @ beta
    i = [r["name"] for r in rows].index(target)
    return float(resid[i]), float((np.sum(resid <= resid[i])) / len(resid))


# ----------------------------------------------------------------------
# Sensitivity on the exemplar subset.
# ----------------------------------------------------------------------
def exemplar_spectrum(name):
    if name.startswith("Z"):
        return "cyclic", int(name[1:]), int(name[1:])
    if name.startswith("2D"):
        return "dihedral", int(name[2:]), 4 * int(name[2:])
    return "poly", None, {"2T": 24, "2O": 48, "2I": 120}[name]


def sensitivity(fixed_tmax=150.0, burn=25):
    out = []
    for name in EXEMPLARS:
        kind, param, _order = exemplar_spectrum(name)
        m = mult.m_exact_series(name, kind, param, None, None, L_MAX)
        l_gap = next(l for l in range(1, L_MAX + 1) if m[l] > 0)
        lam1 = l_gap * (l_gap + 2)
        amp, loglam, _ls = ss.build_modes(m, l_max=L_MAX)

        # equal-N baseline (as in the main sweep)
        res_n, _z, _s = ss.run_series(m, target_N=TARGET_N, l_max=L_MAX,
                                      samples_per_period=SPP)
        # equal-coverage: common fixed t_max, zero count varies
        res_c, _zc, _sc = ss.run_series(m, l_max=L_MAX, t_max=fixed_tmax,
                                        samples_per_period=SPP)
        # burn-in: discard the first `burn` zeros, then unfold
        zeros, _tmax = ss.find_zeros_target(amp, loglam, TARGET_N + burn,
                                            samples_per_period=SPP)
        spac_b = ss.unfold(zeros[burn:], degree=3)
        var_b = float(np.var(spac_b))
        _c, ratio_b, _ = ss.classify_cvm(spac_b)

        out.append({
            "name": name, "lam1": lam1,
            "var_equalN": res_n["var"], "var_equalcov": res_c["var"],
            "var_burnin": var_b,
            "nz_equalcov": res_c["n_zeros"],
        })
    return out


def _spear_over(rows, key_x="lam1", key_y="var"):
    x = [r[key_x] for r in rows]
    y = [r[key_y] for r in rows]
    if len(set(x)) < 2:
        return float("nan")
    return spearmanr(x, y).correlation


# ----------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------
def main():
    print(f"pre-reg {PREREG}, substrate {SUBSTRATE}, l_max={L_MAX}, N={TARGET_N}")
    rows = run_sweep()
    flagged = [r["name"] for r in rows if not r["reached_N"]]
    print(f"population: {len(rows)} groups; reached N=500: "
          f"{sum(r['reached_N'] for r in rows)}; flagged short: {flagged or 'none'}")

    # lam1 distribution
    from collections import Counter
    print("lam1 distribution:", dict(sorted(Counter(r["lam1"] for r in rows).items())))
    print("class distribution:", dict(Counter(r["class"] for r in rows)))

    lam1 = [r["lam1"] for r in rows]
    var = [r["var"] for r in rows]
    ratio = [r["ratio"] for r in rows]
    order = [r["order"] for r in rows]
    dens = [r["density"] for r in rows]
    solv = [int(r["solvable"]) for r in rows]

    print("\n=== H_gap (primary: Spearman(lam1, Var)) ===")
    rho, z, p = perm_z(lam1, var)
    print(f"  pooled: rho={rho:+.3f}  perm-z={z:+.2f}  p={p:.4f}")
    for fam in ("cyclic", "dihedral", "polyhedral"):
        sub = [r for r in rows if r["family"] == fam]
        rr = _spear_over(sub)
        print(f"  {fam:11s} (n={len(sub):2d}): Spearman(lam1,Var)="
              f"{'N/A (lam1 constant)' if np.isnan(rr) else f'{rr:+.3f}'}")

    print("\n=== H_gap (secondary: Spearman(lam1, ratio)) ===")
    rho2, z2, p2 = perm_z(lam1, ratio)
    print(f"  pooled: rho={rho2:+.3f}  perm-z={z2:+.2f}  p={p2:.4f}")

    print("\n=== Confound partials (Var vs lam1, controlling ...) ===")
    print(f"  partial(lam1,Var | order)          = {rank_partial(lam1, var, [order]):+.3f}")
    print(f"  partial(lam1,Var | density)        = {rank_partial(lam1, var, [dens]):+.3f}")
    print(f"  partial(lam1,Var | order,density)  = {rank_partial(lam1, var, [order, dens]):+.3f}")

    print("\n=== H_correlate (solvable vs Var | lam1) ===")
    rp = rank_partial(solv, var, [lam1])
    print(f"  rank-partial(solvable, Var | lam1) = {rp:+.3f}  "
          f"[degenerate/low-power: non-solvable class = {{2I}} only]")
    resid_2I, pct_2I = outlier_percentile(rows)
    print(f"  2I residual on the lam1->Var trend: resid={resid_2I:+.2f}, "
          f"percentile={pct_2I:.2f}  (0.5 = sits on the trend)")
    print(f"  solvable GUE-side groups: "
          f"{[r['name'] for r in rows if r['class']=='GUE' and r['solvable']]}")

    print("\n=== Sensitivity (exemplars): H_gap sign under windowing variants ===")
    sens = sensitivity()
    for tag, key in (("equal-N", "var_equalN"), ("equal-cov", "var_equalcov"),
                     ("burn-in", "var_burnin")):
        rr = _spear_over(sens, key_y=key)
        print(f"  Spearman(lam1,Var) [{tag:9s}] = {rr:+.3f}")

    _write_csvs(rows, sens, {
        "H_gap_var": (rho, z, p), "H_gap_ratio": (rho2, z2, p2),
        "partial_order": rank_partial(lam1, var, [order]),
        "partial_density": rank_partial(lam1, var, [dens]),
        "partial_order_density": rank_partial(lam1, var, [order, dens]),
        "H_correlate_partial": rp, "twoI_resid": resid_2I, "twoI_pct": pct_2I,
    })
    print("\nwrote part2_summary.csv, part2_stats.csv, part2_sensitivity.csv")


def _write_csvs(rows, sens, stats):
    with open(os.path.join(DIR, "part2_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gamma", "family", "order", "solvable", "l_gap", "lambda1",
                    "n_active", "density", "n_zeros", "t_max", "Var", "ratio_cvm",
                    "class_cvm", "reached_N"])
        for r in rows:
            w.writerow([r["name"], r["family"], r["order"], r["solvable"],
                        r["l_gap"], r["lam1"], r["n_active"], round(r["density"], 4),
                        r["n_zeros"], round(r["t_max"], 2), round(r["var"], 5),
                        round(r["ratio"], 4), r["class"], r["reached_N"]])

    with open(os.path.join(DIR, "part2_stats.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["statistic", "value", "perm_z", "perm_p"])
        w.writerow(["H_gap_primary_Spearman(lam1,Var)", round(stats["H_gap_var"][0], 4),
                    round(stats["H_gap_var"][1], 3), round(stats["H_gap_var"][2], 5)])
        w.writerow(["H_gap_secondary_Spearman(lam1,ratio)", round(stats["H_gap_ratio"][0], 4),
                    round(stats["H_gap_ratio"][1], 3), round(stats["H_gap_ratio"][2], 5)])
        w.writerow(["partial(lam1,Var|order)", round(stats["partial_order"], 4), "", ""])
        w.writerow(["partial(lam1,Var|density)", round(stats["partial_density"], 4), "", ""])
        w.writerow(["partial(lam1,Var|order,density)",
                    round(stats["partial_order_density"], 4), "", ""])
        w.writerow(["H_correlate_rankpartial(solvable,Var|lam1)",
                    round(stats["H_correlate_partial"], 4), "", ""])
        w.writerow(["2I_residual_on_trend", round(stats["twoI_resid"], 3),
                    "", f"percentile={round(stats['twoI_pct'],3)}"])

    with open(os.path.join(DIR, "part2_sensitivity.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["exemplar", "lambda1", "Var_equalN", "Var_equalcov",
                    "Var_burnin", "nz_equalcov"])
        for s in sens:
            w.writerow([s["name"], s["lam1"], round(s["var_equalN"], 4),
                        round(s["var_equalcov"], 4), round(s["var_burnin"], 4),
                        s["nz_equalcov"]])


if __name__ == "__main__":
    main()
