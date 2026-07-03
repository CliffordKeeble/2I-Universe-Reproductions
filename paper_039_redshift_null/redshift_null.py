#!/usr/bin/env python3
"""
Paper 39 redshift-mechanism null test (Pattern 75 used as a sword).

Mechanism under test:  photon loses eps = alpha^k per cell over N = R/l_cell
cells, giving  z = exp(N*eps) - 1.  Paper 39 v1.7 reports z ~ 1096 with k=19
(=5+11+3) at l_cell = 0.5 fm.

This script answers, in exact compute:
  Test 1  what z the mechanism actually produces at the three corpus cell scales.
  Test 2  the null: is z=1091 a special point in the (exponent k, cell l) family,
          or just reachable?  Density, isolation, backsolve transparency.
  Test 3  the 27/5 Omega_c/Omega_b neighbourhood null (independent, small).

Run:  python redshift_null.py
Writes:  z_by_cell.csv, null_grid_summary.csv, backsolve_by_cell.csv to this dir.
All exponential arithmetic in mpmath at 50 digits (z is exp of a tuned exponent;
the round number is where the rot starts -- Paper 66 section 3.2).
"""

import csv
import os
import mpmath as mp

mp.mp.dps = 50  # 50 significant digits; carry ten, throw none away

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# Inputs (constants, full precision -- pinned, never the round figures)
# ----------------------------------------------------------------------
alpha       = mp.mpf("7.2973525693e-3")        # CODATA 2018 fine-structure constant
alpha_inv   = mp.mpf("137.035999084")
hbar_c      = mp.mpf("197.3269804")            # MeV.fm
m_p         = mp.mpf("938.27208816")           # MeV (proton mass)
lambdabar_p = hbar_c / m_p                      # reduced proton Compton wavelength, fm
u           = lambdabar_p / 4                    # corpus ruler (Seam A re-anchor)
R_Gly       = mp.mpf("13.8")                     # v1.7 radius (OBSERVED/contested, fixed input here)
Gly_to_fm   = mp.mpf("9.4607304725808e15") * mp.mpf("1e9") * mp.mpf("1e15")  # 1 Gly in fm
R_fm        = R_Gly * Gly_to_fm

Z_OBS = mp.mpf("1091")          # observed recombination redshift (paper headline 1096; gates use 1091)
E_TARGET = mp.log(1 + Z_OBS)    # exponent N*eps that lands exactly on z=1091


def z_of(k, l_fm):
    """z = exp((R/l) * alpha^k) - 1, full precision."""
    eps = alpha ** k
    N = R_fm / mp.mpf(l_fm)
    E = N * eps
    return E, mp.e ** E - 1


def l_star(k):
    """cell length (fm) that makes z exactly Z_OBS for exponent k."""
    return R_fm * (alpha ** k) / E_TARGET


def fmt(x, n=8):
    return mp.nstr(x, n)


print("=" * 72)
print("DERIVED INPUTS")
print("=" * 72)
print(f"  lambdabar_p = hbar_c/m_p   = {fmt(lambdabar_p,12)} fm")
print(f"  u = lambdabar_p/4          = {fmt(u,12)} fm")
print(f"  R_universe                 = {fmt(R_fm,12)} fm  ({R_Gly} Gly)")
print(f"  eps = alpha^19             = {mp.nstr(alpha**19, 8)}")
print(f"  E_target = ln(1+{int(Z_OBS)})       = {fmt(E_TARGET,10)}")
print()

# ----------------------------------------------------------------------
# Test 1 -- z across the three corpus cell scales (k = 19 fixed)
# ----------------------------------------------------------------------
print("=" * 72)
print("TEST 1 -- z at the three corpus cell scales, exponent k=19 fixed")
print("=" * 72)
cells = [
    ("0.5 fm  (v1.7 as-written)",      mp.mpf("0.5")),
    ("u = lambdabar_p/4 (Seam A)",     u),
    ("lambdabar_p (Compton rung)",     lambdabar_p),
]
t1_rows = []
for label, l in cells:
    E, z = z_of(19, l)
    print(f"  l = {label:32s}  l={fmt(l,10)} fm")
    print(f"      N = R/l        = {mp.nstr(R_fm/l, 8)}")
    print(f"      N*eps (exponent) = {fmt(E,10)}")
    print(f"      z              = {mp.nstr(z, 10)}")
    print()
    t1_rows.append({
        "cell_label": label, "l_fm": fmt(l, 12),
        "N": mp.nstr(R_fm / l, 10), "N_eps_exponent": fmt(E, 12),
        "z": mp.nstr(z, 12),
    })

with open(os.path.join(HERE, "z_by_cell.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(t1_rows[0].keys()))
    w.writeheader()
    w.writerows(t1_rows)

# also: what cell does k=19 NEED to hit z=1091 exactly?
l19 = l_star(19)
print(f"  >> backsolve: k=19 hits z={int(Z_OBS)} exactly at l = {fmt(l19,10)} fm")
print(f"     (paper's 0.5 fm is {fmt(100*(mp.mpf('0.5')-l19)/l19,4)} % away from that)")
print()

# ----------------------------------------------------------------------
# Test 2 -- the null: is z=1091 special, or just reachable?
# ----------------------------------------------------------------------
print("=" * 72)
print("TEST 2 -- null distribution over elegant (k, l) family")
print("=" * 72)

# "elegant" integer exponents: corpus primes 3..23 plus small-prime sums
# (the paper's own 19 = 5+11+3 is in this family).
ELEGANT_K = sorted({3, 5, 7, 11, 13, 17, 19, 23,
                    5 + 3,           # 8
                    11 + 3,          # 14
                    5 + 11,          # 16
                    5 + 11 + 3,      # 19 (dup)
                    7 + 11,          # 18
                    7 + 13,          # 20
                    11 + 11,         # 22
                    })

L_LO, L_HI = mp.mpf("0.01"), mp.mpf("1.0")     # cell grid span (fm)
N_LGRID = 4000                                  # log-uniform points across the span
log_lo, log_hi = mp.log(L_LO), mp.log(L_HI)

# z within +/-5% of Z_OBS  <=>  exponent E within this band
E_lo = mp.log(1 + mp.mpf("0.95") * Z_OBS)
E_hi = mp.log(1 + mp.mpf("1.05") * Z_OBS)
print(f"  +/-5% of z={int(Z_OBS)} -> z in [{mp.nstr(0.95*Z_OBS,6)}, {mp.nstr(1.05*Z_OBS,6)}]")
print(f"  -> exponent band E in [{fmt(E_lo,8)}, {fmt(E_hi,8)}]  (width {fmt(E_hi-E_lo,4)})")
print()

# 2a -- for each elegant k, the unique l*(k) that hits z=1091, and whether
#       it lands inside the plausible cell window [0.01, 1.0] fm.
print("  2a -- l*(k): cell that makes z exactly 1091, per elegant exponent")
print(f"  {'k':>3}  {'l*(k) [fm]':>16}  {'in [0.01,1.0]?':>14}")
in_range_k = []
for k in ELEGANT_K:
    ls = l_star(k)
    inside = (L_LO <= ls <= L_HI)
    if inside:
        in_range_k.append(k)
    print(f"  {k:>3}  {mp.nstr(ls,8):>16}  {('YES' if inside else 'no'):>14}")
print(f"  => elegant exponents whose matching cell is in-range: {in_range_k}")
print()

# 2b -- density over the full (k, l) grid: fraction of pairs within +/-5%.
total_pairs = 0
hit_pairs = 0
for k in ELEGANT_K:
    eps = alpha ** k
    for i in range(N_LGRID):
        t = mp.mpf(i) / (N_LGRID - 1)
        l = mp.e ** (log_lo + t * (log_hi - log_lo))
        E = (R_fm / l) * eps
        total_pairs += 1
        if E_lo <= E <= E_hi:
            hit_pairs += 1
frac = mp.mpf(hit_pairs) / total_pairs
print(f"  2b -- grid density:  {hit_pairs} / {total_pairs} pairs within +/-5% of z=1091")
print(f"        fraction = {mp.nstr(frac,6)}  ({mp.nstr(100*frac,4)} %)")
print()

# 2c -- isolation of the 19/0.5fm point: width of the +/-5% l-band at k=19,
#       and whether any OTHER elegant k also lands near 0.5 fm.
print("  2c -- isolation of the (k=19, l~0.5fm) point")
eps19 = alpha ** 19
# l-interval giving +/-5%-z at k=19:  E in [E_lo,E_hi], E=(R/l)eps -> l=(R eps)/E
l_hi_band = R_fm * eps19 / E_lo   # larger l -> smaller E
l_lo_band = R_fm * eps19 / E_hi
print(f"      +/-5%-z l-band at k=19: [{fmt(l_lo_band,8)}, {fmt(l_hi_band,8)}] fm")
print(f"      band width = {fmt(l_hi_band-l_lo_band,6)} fm "
      f"= {fmt(100*(l_hi_band-l_lo_band)/l19,4)} % of l*(19)")
# nearest other elegant k's matching cell
others = [(k, l_star(k)) for k in ELEGANT_K if k != 19]
near = min(others, key=lambda kv: abs(mp.log(kv[1]) - mp.log(l19)))
print(f"      nearest other elegant k by cell: k={near[0]} at l*={fmt(near[1],8)} fm")
print(f"      ratio l*(that)/l*(19) = {fmt(near[1]/l19,6)}  "
      f"(consecutive elegant k differ by ~1/alpha = {fmt(1/alpha,6)}x in cell)")
print()

# 2d -- backsolve transparency: per cell, the k and Lambda_QCD required for z=1091.
print("  2d -- backsolve: per cell, required exponent k and implied Lambda_QCD")
print(f"  {'cell':>30}  {'l [fm]':>12}  {'req. k':>10}  {'Lambda_QCD=hbar_c/l [MeV]':>26}")
backsolve_cells = [
    ("0.5 fm (v1.7)",            mp.mpf("0.5")),
    ("u = lambdabar_p/4",        u),
    ("lambdabar_p",              lambdabar_p),
    ("l*(19) backsolved",        l19),
]
bs_rows = []
for label, l in backsolve_cells:
    # required continuous k: (R/l) alpha^k = E_target -> k = ln(E_target l / R)/ln(alpha)
    k_req = mp.log(E_TARGET * l / R_fm) / mp.log(alpha)
    Lam = hbar_c / l
    print(f"  {label:>30}  {fmt(l,8):>12}  {fmt(k_req,7):>10}  {fmt(Lam,8):>26}")
    bs_rows.append({"cell_label": label, "l_fm": fmt(l, 12),
                    "k_required_for_z1091": fmt(k_req, 10),
                    "Lambda_QCD_MeV": fmt(Lam, 10)})
print()
with open(os.path.join(HERE, "backsolve_by_cell.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(bs_rows[0].keys()))
    w.writeheader()
    w.writerows(bs_rows)

with open(os.path.join(HERE, "null_grid_summary.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["quantity", "value"])
    w.writerow(["elegant_k_set", ";".join(map(str, ELEGANT_K))])
    w.writerow(["l_grid_span_fm", f"[{L_LO},{L_HI}]"])
    w.writerow(["l_grid_points", N_LGRID])
    w.writerow(["pm5pct_hits", hit_pairs])
    w.writerow(["total_pairs", total_pairs])
    w.writerow(["pm5pct_fraction", mp.nstr(frac, 8)])
    w.writerow(["in_range_elegant_k", ";".join(map(str, in_range_k))])
    w.writerow(["l_star_19_fm", fmt(l19, 12)])
    w.writerow(["pm5pct_l_band_fm", f"[{fmt(l_lo_band,10)},{fmt(l_hi_band,10)}]"])

# ----------------------------------------------------------------------
# Test 3 -- Seam C residual: the 27/5 null for Omega_c/Omega_b
# ----------------------------------------------------------------------
print("=" * 72)
print("TEST 3 -- 27/5 neighbourhood null for Omega_c/Omega_b")
print("=" * 72)
# Planck 2018 base-LCDM (TT,TE,EE+lowE+lensing): Omega_b h^2, Omega_c h^2
ob_h2, ob_err = mp.mpf("0.02237"), mp.mpf("0.00015")
oc_h2, oc_err = mp.mpf("0.1200"),  mp.mpf("0.0012")
ratio = oc_h2 / ob_h2
rel = mp.sqrt((oc_err / oc_h2) ** 2 + (ob_err / ob_h2) ** 2)
ratio_err = ratio * rel
print(f"  Planck 2018: Omega_b h^2 = {ob_h2} +/- {ob_err}")
print(f"               Omega_c h^2 = {oc_h2} +/- {oc_err}")
print(f"  Omega_c/Omega_b = {fmt(ratio,8)} +/- {fmt(ratio_err,4)}")
print(f"  1-sigma band: [{fmt(ratio-ratio_err,6)}, {fmt(ratio+ratio_err,6)}]")
print(f"  2-sigma band: [{fmt(ratio-2*ratio_err,6)}, {fmt(ratio+2*ratio_err,6)}]")
print()
candidates = [
    ("27/5",       mp.mpf(27) / 5),
    ("3*sqrt(3)",  3 * mp.sqrt(3)),
    ("16/3",       mp.mpf(16) / 3),
    ("2*phi^2",    2 * ((1 + mp.sqrt(5)) / 2) ** 2),
    ("21/4",       mp.mpf(21) / 4),
    ("5.4 (=27/5)", mp.mpf("5.4")),
]
# de-dup 5.4 vs 27/5
print(f"  {'form':>12}  {'value':>10}  {'|val-ratio|':>12}  {'within 1s?':>10}  {'within 2s?':>10}")
inside_1s, inside_2s = [], []
for name, val in [("27/5", mp.mpf(27)/5), ("3*sqrt(3)", 3*mp.sqrt(3)),
                  ("16/3", mp.mpf(16)/3), ("2*phi^2", 2*((1+mp.sqrt(5))/2)**2),
                  ("21/4", mp.mpf(21)/4)]:
    d = abs(val - ratio)
    w1 = d <= ratio_err
    w2 = d <= 2 * ratio_err
    if w1: inside_1s.append(name)
    if w2: inside_2s.append(name)
    print(f"  {name:>12}  {fmt(val,7):>10}  {fmt(d,5):>12}  "
          f"{('YES' if w1 else 'no'):>10}  {('YES' if w2 else 'no'):>10}")
closest = min([("27/5", mp.mpf(27)/5), ("3*sqrt(3)", 3*mp.sqrt(3)),
               ("16/3", mp.mpf(16)/3), ("2*phi^2", 2*((1+mp.sqrt(5))/2)**2),
               ("21/4", mp.mpf(21)/4)], key=lambda kv: abs(kv[1] - ratio))
print()
print(f"  closest simple form to {fmt(ratio,6)}: {closest[0]} = {fmt(closest[1],6)}")
print(f"  simple forms inside 1-sigma: {inside_1s}")
print(f"  simple forms inside 2-sigma: {inside_2s}")
print()
print("Done. CSVs written to:", HERE)
