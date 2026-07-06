"""
Paper 116 Part 1 -- reproduce the DERIVED tier (AT1-AT7).

Reconstructs l_max from the printed Active-l counts (structural inversion) and
t_max from the printed zero counts (reconstruction from a published constraint,
never tuning a tested statistic -- Gate G3). The tested statistics (Var, class)
then fall where they fall. Sensitivity grid: l_max x{1/2,1,2}, t_max x{1/2,1,2},
unfolding degree {2,3,4}; binned discriminant bins {20,40,80}. Verdicts use the
three-word vocabulary.

Run:  python part1_reproduce.py   -> writes part1_*.csv, prints the digest.
"""
from __future__ import annotations

import csv
import functools
import itertools
import os

import numpy as np

import mult
import spectral_stats as ss

DIR = os.path.dirname(os.path.abspath(__file__))
SPP = 64          # samples per shortest period for zero bracketing
L_BIG = 3000      # ceiling for m(l): must exceed 2x the largest reconstructed l_max
BINS = (20, 40, 80)

# name -> (mult-name, kind, param, Active-l target, zero-count target, paper class)
PART1_GROUPS = {
    "2I":   ("2I", "poly", None, 135, 499, "GUE"),
    "2O":   ("2O", "poly", None, 142, 521, "GUE"),
    "Z120": ("Z120", "cyclic", 120, 271, 231, "GOE"),
    "Z60":  ("Z60", "cyclic", 60, 285, 245, "GOE"),
}
# paper's own (Var, ratio, class) -- ratio is the paper's unspecified statistic.
PAPER = {
    "2I": (0.117, 0.67, "GUE"), "2O": (0.122, 0.65, "GUE"),
    "Z120": (0.317, 2.13, "GOE"), "Z60": (0.371, 1.97, "GOE"),
}


@functools.lru_cache(maxsize=None)
def m_array(name, kind, param, l_max):
    """Cached m(l) array. Callers that mutate must copy first (list(m))."""
    return mult.m_exact_series(name, kind, param, None, None, l_max)


# --------------------------------------------------------------------------
# Reconstruction from published constraints.
# --------------------------------------------------------------------------
def reconstruct_lmax(m_big, target_active):
    """Smallest l_max with exactly target_active modes l>=1, m(l)>0."""
    c = 0
    for l in range(1, len(m_big)):
        if m_big[l] > 0:
            c += 1
            if c == target_active:
                return l
    raise ValueError(f"only {c} active modes in l<={len(m_big)-1}; raise L_BIG")


def reconstruct_tmax(amp, loglam, target_zeros):
    z, t_max = ss.find_zeros_target(amp, loglam, target_zeros, samples_per_period=SPP)
    return t_max, len(z)


def group_spectrum(name):
    """Return (m_big, l_max, t_max) for a Part-1 baseline group."""
    nm, kind, param, ta, tz, _cls = PART1_GROUPS[name]
    m_big = m_array(nm, kind, param, L_BIG)
    l_max = reconstruct_lmax(m_big, ta)
    amp, loglam, _ls = ss.build_modes(m_big, l_max=l_max)
    t_max, _nz = reconstruct_tmax(amp, loglam, tz)
    return m_big, l_max, t_max


# --------------------------------------------------------------------------
# Discriminants and verdict.
# --------------------------------------------------------------------------
def binned_class(spac, nbins):
    msd = {k: ss.binned_msd(spac, k, nbins) for k in ("Poisson", "GOE", "GUE")}
    return min(msd, key=msd.get)


def verdict(center_class, paper_class, grid_classes, binned_classes):
    if center_class != paper_class:
        return "DOES-NOT-REPRODUCE"
    all_grid = all(c == paper_class for c in grid_classes)
    all_binned = all(c == paper_class for c in binned_classes)
    if all_grid and all_binned:
        return "ROBUST-REPRODUCES"
    return "REPRODUCES-UNDER-DOCUMENTED-CHOICES"


# --------------------------------------------------------------------------
# AT2 -- baseline four rows + sensitivity grid -> verdict.
# --------------------------------------------------------------------------
def at2_baseline():
    rows, sens_rows = [], []
    for name in ("2I", "2O", "Z120", "Z60"):
        paper_cls = PART1_GROUPS[name][5]
        m_big, l0, t0 = group_spectrum(name)

        res, _z, spac = ss.run_series(m_big, l_max=l0, t_max=t0, degree=3,
                                      samples_per_period=SPP)
        c_cvm, ratio, msd = ss.classify_cvm(spac)
        binned = [binned_class(spac, b) for b in BINS]

        grid_classes = []
        for lf, tf, deg in itertools.product((0.5, 1.0, 2.0), (0.5, 1.0, 2.0),
                                             (2, 3, 4)):
            lmax = max(12, int(round(l0 * lf)))
            tmax = t0 * tf
            r, _zz, sp = ss.run_series(m_big, l_max=lmax, t_max=tmax, degree=deg,
                                       samples_per_period=SPP)
            gc, gratio, _ = ss.classify_cvm(sp)
            grid_classes.append(gc)
            sens_rows.append([name, lf, tf, deg, lmax, round(tmax, 2),
                              r["n_zeros"], round(r["var"], 4), round(gratio, 4), gc])

        v = verdict(c_cvm, paper_cls, grid_classes, binned)
        pv, pr, _pc = PAPER[name]
        rows.append({
            "name": name, "l_max": l0, "t_max": round(t0, 2),
            "n_active": res["n_active"], "n_zeros": res["n_zeros"],
            "var": res["var"], "ratio_cvm": ratio, "class_cvm": c_cvm,
            "class_binned": "/".join(binned), "verdict": v,
            "paper_var": pv, "paper_ratio": pr, "paper_class": paper_cls,
            "grid_flip": sum(1 for gc in grid_classes if gc != paper_cls),
            "spac": spac,
        })
    return rows, sens_rows


# --------------------------------------------------------------------------
# AT3 / AT5 -- impose a gap on a cyclic group (fixed natural window).
# --------------------------------------------------------------------------
def impose_gap(m, l_gap):
    m2 = list(m)
    for l in range(0, l_gap):
        m2[l] = 0
    return m2


def at_gap_experiment(name, gaps):
    m_big, l0, t0 = group_spectrum(name)
    out = []
    # natural row for reference
    r0, _z, sp0 = ss.run_series(m_big, l_max=l0, t_max=t0, samples_per_period=SPP)
    c0, ratio0, _ = ss.classify_cvm(sp0)
    out.append(("natural", 1, _lam1(m_big),
                r0["n_zeros"], r0["var"], ratio0, c0))
    for lg in gaps:
        mg = impose_gap(m_big, lg)
        r, _zz, sp = ss.run_series(mg, l_max=l0, t_max=t0, samples_per_period=SPP)
        c, ratio, _ = ss.classify_cvm(sp)
        out.append((f"gap l={lg}", lg, lg * (lg + 2), r["n_zeros"], r["var"], ratio, c))
    return out, l0, t0


def _lam1(m):
    l = 1
    while m[l] == 0:
        l += 1
    return l * (l + 2)


# --------------------------------------------------------------------------
# AT4 -- fill low modes of 2I (all 2-placements in l<12, plus a 10-mode case).
# --------------------------------------------------------------------------
def at4_fill_2I():
    m_big, l0, t0 = group_spectrum("2I")
    cands = list(range(1, 12))          # sub-gap l = 1..11
    placements = []
    for combo in itertools.combinations(cands, 2):
        m2 = list(m_big)
        for l in combo:
            m2[l] = 1
        r, _z, sp = ss.run_series(m2, l_max=l0, t_max=t0, samples_per_period=SPP)
        c, ratio, _ = ss.classify_cvm(sp)
        placements.append((combo, r["n_zeros"], r["var"], ratio, c))
    # 10-mode case: l = 1..10
    m10 = list(m_big)
    for l in range(1, 11):
        m10[l] = 1
    r, _z, sp = ss.run_series(m10, l_max=l0, t_max=t0, samples_per_period=SPP)
    c10, ratio10, _ = ss.classify_cvm(sp)
    ten = (tuple(range(1, 11)), r["n_zeros"], r["var"], ratio10, c10)
    return placements, ten, l0, t0


# --------------------------------------------------------------------------
# AT6 -- the 2T question: correct spectrum vs counterfactual (mode at l=4).
# --------------------------------------------------------------------------
def at6_2T(l_max=200, target_N=500):
    m2T = m_array("2T", "poly", None, L_BIG)
    out = {}
    for label, mm in (("correct(gap l=6)", m2T),
                      ("counterfactual(+l=4)", _with_mode(m2T, 4))):
        res = {}
        for lmx in (l_max // 2, l_max, l_max * 2):
            r, _z, sp = ss.run_series(mm, target_N=target_N, l_max=lmx,
                                      samples_per_period=SPP)
            c, ratio, _ = ss.classify_cvm(sp)
            res[lmx] = (r["n_zeros"], r["var"], ratio, c, _lam1(mm))
        out[label] = res
    return out


def _with_mode(m, l):
    m2 = list(m)
    m2[l] = max(1, m2[l])
    return m2


# --------------------------------------------------------------------------
# Driver.
# --------------------------------------------------------------------------
def main():
    print("=== AT2 baseline + sensitivity ===")
    rows, sens = at2_baseline()
    for r in rows:
        print(f"{r['name']:5s} l_max={r['l_max']:4d} t_max={r['t_max']:8.2f} "
              f"nz={r['n_zeros']:4d} Var={r['var']:.4f} ratio_cvm={r['ratio_cvm']:.3f} "
              f"class={r['class_cvm']:7s} binned={r['class_binned']:12s} "
              f"flips={r['grid_flip']:2d}/27  -> {r['verdict']}")

    print("\n=== AT3 Z120 gap experiment (fixed natural window) ===")
    z120_gap, _, _ = at_gap_experiment("Z120", (4, 8, 12))
    for lab, lg, lam, nz, var, ratio, c in z120_gap:
        print(f"  {lab:10s} lam1={lam:4d} nz={nz:4d} Var={var:.4f} "
              f"ratio_cvm={ratio:.3f} class={c}")

    print("\n=== AT5 Z60 gap experiment ===")
    z60_gap, _, _ = at_gap_experiment("Z60", (4, 8))
    for lab, lg, lam, nz, var, ratio, c in z60_gap:
        print(f"  {lab:10s} lam1={lam:4d} nz={nz:4d} Var={var:.4f} "
              f"ratio_cvm={ratio:.3f} class={c}")

    print("\n=== AT4 fill 2I low modes ===")
    placements, ten, _, _ = at4_fill_2I()
    n_goe = sum(1 for _c, _nz, _v, _r, c in placements if c == "GOE")
    n_gue = sum(1 for _c, _nz, _v, _r, c in placements if c == "GUE")
    n_poi = sum(1 for _c, _nz, _v, _r, c in placements if c == "Poisson")
    print(f"  55 two-mode placements: GOE={n_goe} GUE={n_gue} Poisson={n_poi}")
    print(f"  10-mode case l=1..10: nz={ten[1]} Var={ten[2]:.4f} "
          f"ratio_cvm={ten[3]:.3f} class={ten[4]}")
    # a few example placements
    for combo, nz, var, ratio, c in placements[:5]:
        print(f"    add {combo}: nz={nz} Var={var:.4f} ratio={ratio:.3f} class={c}")

    print("\n=== AT6 the 2T question (printed row: Var 0.322, ratio 0.85) ===")
    at6 = at6_2T()
    for label, res in at6.items():
        for lmx, (nz, var, ratio, c, lam1) in res.items():
            print(f"  {label:22s} l_max={lmx:4d} lam1={lam1:3d} nz={nz} "
                  f"Var={var:.4f} ratio_cvm={ratio:.3f} class={c}")

    print("\n=== AT7 sub-GUE trend (Var vs gap width) ===")
    trend = [("Z120", z120_gap), ("Z60", z60_gap)]
    for gname, g in trend:
        pts = [(lam, var) for _lab, _lg, lam, _nz, var, _r, _c in g]
        pts_sorted = sorted(pts)
        mono = all(pts_sorted[i][1] >= pts_sorted[i + 1][1]
                   for i in range(len(pts_sorted) - 1))
        print(f"  {gname}: (lam1,Var) {[(l, round(v,4)) for l,v in pts_sorted]} "
              f"monotone_decreasing={mono}")

    _write_csvs(rows, sens, z120_gap, z60_gap, placements, ten, at6)
    print("\nwrote part1_summary.csv, part1_sensitivity.csv, part1_gap_experiments.csv, "
          "part1_at4_placements.csv, part1_at6_2T.csv")


def _write_csvs(rows, sens, z120_gap, z60_gap, placements, ten, at6):
    with open(os.path.join(DIR, "part1_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["group", "l_max", "t_max", "n_active", "n_zeros", "Var",
                    "ratio_cvm", "class_cvm", "class_binned_20_40_80", "verdict",
                    "paper_Var", "paper_ratio", "paper_class", "grid_flips_of_27"])
        for r in rows:
            w.writerow([r["name"], r["l_max"], r["t_max"], r["n_active"],
                        r["n_zeros"], round(r["var"], 5), round(r["ratio_cvm"], 4),
                        r["class_cvm"], r["class_binned"], r["verdict"],
                        r["paper_var"], r["paper_ratio"], r["paper_class"],
                        r["grid_flip"]])

    with open(os.path.join(DIR, "part1_sensitivity.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["group", "l_factor", "t_factor", "degree", "l_max", "t_max",
                    "n_zeros", "Var", "ratio_cvm", "class_cvm"])
        w.writerows(sens)

    with open(os.path.join(DIR, "part1_gap_experiments.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["group", "row", "l_gap_or", "lambda1", "n_zeros", "Var",
                    "ratio_cvm", "class_cvm"])
        for lab, lg, lam, nz, var, ratio, c in z120_gap:
            w.writerow(["Z120", lab, lg, lam, nz, round(var, 5), round(ratio, 4), c])
        for lab, lg, lam, nz, var, ratio, c in z60_gap:
            w.writerow(["Z60", lab, lg, lam, nz, round(var, 5), round(ratio, 4), c])

    with open(os.path.join(DIR, "part1_at4_placements.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["added_modes", "n_zeros", "Var", "ratio_cvm", "class_cvm"])
        for combo, nz, var, ratio, c in placements:
            w.writerow(["+".join(map(str, combo)), nz, round(var, 5),
                        round(ratio, 4), c])
        w.writerow(["+".join(map(str, ten[0])), ten[1], round(ten[2], 5),
                    round(ten[3], 4), ten[4]])

    with open(os.path.join(DIR, "part1_at6_2T.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["spectrum", "l_max", "lambda1", "n_zeros", "Var", "ratio_cvm",
                    "class_cvm"])
        for label, res in at6.items():
            for lmx, (nz, var, ratio, c, lam1) in res.items():
                w.writerow([label, lmx, lam1, nz, round(var, 5), round(ratio, 4), c])


if __name__ == "__main__":
    main()
