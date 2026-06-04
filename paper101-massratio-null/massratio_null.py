"""
massratio_null.py -- reconstruction of the Paper 101 v3.0 section 5.3
mass-ratio null test.

THIS IS A FAITHFUL RE-RUN OF A PUBLISHED, PRE-REGISTERED DESIGN.
The search space, integer-correction set, operations, targets and decision
window are transcribed from Paper 101 section 5.3 (see PRE-REGISTRATION.md).
Nothing here is tuned to improve the match.

Design (locked):
  * search space : every surviving 2I level up to l = 120 (from spectrum.py),
                   eigenvalue lambda_l = l(l+2)
  * corrections  : a single icosahedral integer from the declared set
                   {1,2,6,11,12,13,18,19,20,30,31,41,42,62,168}
  * operations   : + or - (one correction, either sign)
  * targets      : m_p/m_e, m_n/m_e, m_mu/m_e, m_pi/m_e, m_tau/m_e
  * window       : 0.01% relative (sensitivity also reported at 0.05% / 0.005%)
  * decision rule: more than one survivor for a given target => NOT distinguished

Seed 42 is declared in the pre-registration for convention; the enumeration is
exhaustive and deterministic, so no random draw is actually consumed. We seed
anyway so the artifact matches the declared protocol.

Run: python massratio_null.py
"""

import csv
import random
import numpy as np

from spectrum import surviving_levels

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

# ---- declared integer-correction set (icosahedral integers) -----------------
CORRECTIONS = [1, 2, 6, 11, 12, 13, 18, 19, 20, 30, 31, 41, 42, 62, 168]
OPERATIONS = [("+", +1), ("-", -1)]

# ---- targets: dimensionless mass ratios in the >100 (integer-eigenvalue)
#      regime. Values recorded explicitly for reproducibility.
#      m_p/m_e, m_n/m_e, m_mu/m_e : CODATA recommended ratios.
#      m_pi/m_e, m_tau/m_e        : PDG 2024 masses / CODATA m_e.
M_E_MEV = 0.51099895000            # CODATA electron mass (MeV)
TARGETS = {
    "m_p/m_e":   1836.152673426,            # CODATA proton-electron mass ratio
    "m_n/m_e":   1838.68366173,             # CODATA neutron-electron mass ratio
    "m_mu/m_e":  206.7682830,               # CODATA muon-electron mass ratio
    "m_pi/m_e":  139.57039 / M_E_MEV,       # PDG 2024 charged pion 139.57039 MeV
    "m_tau/m_e": 1776.86 / M_E_MEV,         # PDG 2024 tau 1776.86 MeV
}

# declared-then-EXCLUDED before running (documented, not silently dropped):
#   m_n/m_p ~ 1.0014 -- lives near 1, unreachable by integer eigenvalues in the
#   >100 regime. Declared in the pre-registration, excluded as out-of-regime.
EXCLUDED = {"m_n/m_p": 1.0013784}

WINDOWS = [1e-4, 5e-4, 5e-5]       # 0.01%, 0.05%, 0.005%
PRIMARY = 1e-4


def enumerate_combinations():
    levels_raw, _max_dev = surviving_levels(120)
    levels = [(l, lam) for (l, lam, mult) in levels_raw if l > 0]
    combos = []
    for (l, lam) in levels:
        for (opsym, opsign) in OPERATIONS:
            for c in CORRECTIONS:
                value = lam + opsign * c
                combos.append((l, lam, opsym, c, value))
    return levels, combos


def survivors_at(combos, window):
    """Return {target: [ (l,lam,op,c,value,rel) ... ]} within window."""
    out = {t: [] for t in TARGETS}
    for (l, lam, opsym, c, value) in combos:
        for t, tv in TARGETS.items():
            rel = abs(value - tv) / tv
            if rel <= window:
                out[t].append((l, lam, opsym, c, value, rel))
    return out


def main():
    levels, combos = enumerate_combinations()
    n_levels = len(levels)
    n_combos = len(combos)

    print("=" * 70)
    print("Paper 101 section 5.3 mass-ratio null -- faithful reconstruction")
    print("=" * 70)
    print(f"seed                     : {SEED}")
    print(f"surviving levels (l>0)   : {n_levels}")
    print(f"corrections              : {len(CORRECTIONS)}  {CORRECTIONS}")
    print(f"operations               : {len(OPERATIONS)} (+/-)")
    print(f"targets                  : {len(TARGETS)}  {list(TARGETS)}")
    n_full = n_combos * len(TARGETS)
    print(f"(mode x op x correction) : {n_combos}")
    print(f"combination count        : {n_full}   "
          f"(= {n_levels} modes x {len(CORRECTIONS)} integers "
          f"x {len(OPERATIONS)} ops x {len(TARGETS)} targets)")
    print(f"  published ~7500 sanity figure: reconstructed = {n_full} "
          f"(~10% under; anchors l=12->168, l=42->1848 reproduce exactly)")
    print("\ntarget values used:")
    for t, tv in TARGETS.items():
        print(f"  {t:10s} = {tv:.6f}")
    print(f"  EXCLUDED (declared, out-of-regime): "
          f"m_n/m_p = {EXCLUDED['m_n/m_p']:.6f}")

    rows = []
    for window in WINDOWS:
        res = survivors_at(combos, window)
        total = sum(len(v) for v in res.values())
        print("\n" + "-" * 70)
        print(f"WINDOW {window*100:.3f}%  -- total survivors: {total}")
        for t in TARGETS:
            hits = sorted(res[t], key=lambda r: r[5])
            if not hits:
                print(f"  {t:10s}: 0")
            else:
                print(f"  {t:10s}: {len(hits)}")
                for (l, lam, opsym, c, value, rel) in hits:
                    print(f"      l={l:3d} lambda={lam:6d} {opsym}{c:<3d} "
                          f"= {value:9.3f}   rel={rel*100:.4f}%")
                    rows.append((f"{window*100:.3f}%", t, l, lam, opsym, c,
                                 value, f"{rel:.3e}"))

    # decision rule at the primary window
    print("\n" + "=" * 70)
    res = survivors_at(combos, PRIMARY)
    distinguishable = True
    for t in TARGETS:
        n = len(res[t])
        if n > 1:
            distinguishable = False
    n_total_primary = sum(len(v) for v in res.values())
    print(f"DECISION (window {PRIMARY*100:.3f}%): total survivors = {n_total_primary}")
    print(f"  per-target max survivors <= 1 ? {'YES' if distinguishable else 'NO'}")
    if distinguishable and n_total_primary >= 1:
        print("  => DISTINGUISHED: the surviving identification(s) stand alone;")
        print("     no target admits more than one in-window hit.")
    elif n_total_primary == 0:
        print("  => no survivors at the primary window.")
    else:
        print("  => NOT DISTINGUISHED: at least one target has multiple hits.")

    with open("matches.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["window", "target", "l", "lambda", "op", "correction",
                    "value", "rel"])
        w.writerows(rows)
    print("\nwrote matches.csv")


if __name__ == "__main__":
    main()
