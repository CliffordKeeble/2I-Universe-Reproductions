"""
neutrino_mass_null.py -- pre-registered menu null for Paper 121 section 3.1.

Tests Mr Adversary's cycle-1 charge: is m_nu3 = m_e * alpha^3 / 4 a parameter-free
derivation, or a choice from a discrete integer menu (a hidden fit)? The design --
menu, target, windows, decision rule -- is LOCKED in PRE-REGISTRATION.md and is NOT
adjusted in light of results.

Two tests:
  (1) mass form     m_pred(p,q) = m_e * alpha^p / q,  Menu A: p in 1..5, q in 1..6.
  (2) generation    R = Delta m^2_31 / Delta m^2_21 ~= 32.54; forms 2^s and k^(d/2).

Run: python neutrino_mass_null.py   (writes matches.csv)
"""

import csv
import random

random.seed(42)  # programme convention; enumeration is exhaustive/deterministic

# ---- locked constants (see PRE-REGISTRATION.md) -----------------------------
ALPHA = 0.0072973525693          # CODATA fine-structure constant (1/137.035999084)
M_E_EV = 0.51099895000e6         # CODATA electron mass in eV (0.51099895 MeV)

# ---- TEST 1: mass form ------------------------------------------------------
P_MENU = [1, 2, 3, 4, 5]         # exponent menu (Menu A)
Q_MENU = [1, 2, 3, 4, 5, 6]      # divisor menu (Menu A)
CANONICAL = (3, 4)

# target: m_nu3 = sqrt(Delta m^2_31), normal ordering, lightest ~ 0
DM2_31 = 2.45e-3                  # eV^2 (central, normal ordering)
TARGET_EV = DM2_31 ** 0.5        # ~0.0495 eV
TARGET_MEV = TARGET_EV * 1e3
TARGET_ALT_MEV = 50.0            # programme's cited "~50 meV"

WINDOWS = [1e-2, 5e-3, 1e-3]     # 1%, 0.5%, 0.1%

# ---- TEST 2: generation scaling --------------------------------------------
DM2_21 = 7.53e-5                 # eV^2 (central)
R_OBS = DM2_31 / DM2_21          # ~32.54
D = 5                            # icosahedral vertex degree
S_MENU = [1, 2, 3, 4, 5, 6]      # for 2^s
K_MENU = [1, 2, 3, 4, 5, 6]      # for k^(d/2)
GEN_WINDOWS = [1e-2, 5e-3, 1e-3, 2e-2, 5e-2]  # 1%,0.5%,0.1% + 2%,5% locators


def mass_candidates():
    out = []
    for p in P_MENU:
        for q in Q_MENU:
            mpred_ev = M_E_EV * (ALPHA ** p) / q
            rel = abs(mpred_ev - TARGET_EV) / TARGET_EV
            rel_alt = abs(mpred_ev * 1e3 - TARGET_ALT_MEV) / TARGET_ALT_MEV
            out.append((p, q, mpred_ev * 1e3, rel, rel_alt))  # mpred in meV
    return out


def gen_candidates():
    out = []
    for s in S_MENU:
        val = 2.0 ** s
        out.append(("2^s", s, val, abs(val - R_OBS) / R_OBS))
    for k in K_MENU:
        val = float(k) ** (D / 2.0)
        out.append(("k^(d/2)", k, val, abs(val - R_OBS) / R_OBS))
    return out


def main():
    print("=" * 72)
    print("Neutrino-menu null (Paper 121 §3.1) -- pre-registered")
    print("=" * 72)
    print(f"seed 42 | alpha={ALPHA} | m_e={M_E_EV:.5f} eV")
    print(f"\nTARGET (mass): m_nu3 = sqrt({DM2_31}) eV = {TARGET_MEV:.4f} meV "
          f"(alt cited 50 meV)")

    # ---- TEST 1 ----
    cand = mass_candidates()
    print(f"\n[TEST 1] mass form  m_e*alpha^p/q  | Menu A: p in {P_MENU}, "
          f"q in {Q_MENU}  -> {len(cand)} candidates")
    print(f"  canonical (3,4) prediction = "
          f"{[c for c in cand if (c[0],c[1])==CANONICAL][0][2]:.4f} meV")

    rows = []
    for window in WINDOWS:
        surv = sorted([c for c in cand if c[3] <= window], key=lambda c: c[3])
        print(f"\n  WINDOW {window*100:.2f}%  -> {len(surv)} survivor(s) "
              f"(vs primary target {TARGET_MEV:.3f} meV):")
        for (p, q, mpred, rel, rel_alt) in surv:
            tag = "  <-- CANONICAL" if (p, q) == CANONICAL else ""
            print(f"      (p={p}, q={q})  {mpred:10.4f} meV   rel={rel*100:.4f}%{tag}")
            rows.append(["mass", f"{window*100:.2f}%", p, q,
                         f"{mpred:.6f}", f"{rel:.4e}", f"{rel_alt:.4e}"])

    # all candidates near the target, for the distribution picture
    print("\n  full distribution (candidates within a factor ~3 of target):")
    for (p, q, mpred, rel, rel_alt) in sorted(cand, key=lambda c: c[3])[:8]:
        print(f"      (p={p}, q={q})  {mpred:12.4f} meV   rel={rel*100:9.3f}%")

    # ---- TEST 2 ----
    gen = gen_candidates()
    print(f"\n[TEST 2] generation scaling  | R_obs = {R_OBS:.4f} "
          f"(= {DM2_31}/{DM2_21}) | claim 2^d=2^{D}=32")
    for window in GEN_WINDOWS:
        surv = sorted([g for g in gen if g[3] <= window], key=lambda g: g[3])
        listing = ", ".join(f"{fam}({n})={val:.3f}[{rel*100:.2f}%]"
                            for (fam, n, val, rel) in surv) or "none"
        print(f"  WINDOW {window*100:.2f}%  -> {len(surv)}: {listing}")
        for (fam, n, val, rel) in surv:
            rows.append(["gen", f"{window*100:.2f}%", fam, n,
                         f"{val:.6f}", f"{rel:.4e}", ""])

    with open("matches.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["test", "window", "p_or_family", "q_or_n",
                    "value_meV_or_ratio", "rel_primary", "rel_alt50meV"])
        w.writerows(rows)
    print("\nwrote matches.csv")


if __name__ == "__main__":
    main()
