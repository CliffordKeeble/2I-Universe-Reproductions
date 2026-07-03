#!/usr/bin/env python3
"""paper66-ladder-null: exhaustive enumeration under the pre-registered class C.

Class C (locked in preregistration.md):
  - sum of at most 3 terms, each k*m, integer 1 <= |k| <= 7
  - monomials: 1, the five invariants V,E,F,D,chi, and all pairwise
    products including squares (21 monomials total).
Outputs: coverage of [30,200], spelling multiplicities, verification of the
paper's six formulas, near-integer p-values under declared conventions.
"""
from itertools import combinations, product
import math

# ---------------------------------------------------------------- monomials
inv = {"V": 12, "E": 30, "F": 20, "D": 3, "chi": 2}
names = list(inv)
mono = {"1": 1}
mono.update(inv)
for i, a in enumerate(names):
    for b in names[i:]:
        mono[f"{a}*{b}"] = inv[a] * inv[b]
mono_items = sorted(mono.items(), key=lambda kv: kv[1])
M = len(mono_items)
print(f"monomial pool ({M}):")
print("  " + ", ".join(f"{k}={v}" for k, v in mono_items))

COEFFS = [k for k in range(-7, 8) if k != 0]
LO, HI = 30, 200

# ------------------------------------------------- exhaustive enumeration
# spell[t] = list of formal expressions hitting integer t in [LO, HI]
spell = {t: [] for t in range(LO, HI + 1)}

def record(val, expr):
    if LO <= val <= HI:
        spell[val].append(expr)

# 1-term
for (mn, mv) in mono_items:
    for k in COEFFS:
        record(k * mv, ((k, mn),))
# 2-term
for (an, av), (bn, bv) in combinations(mono_items, 2):
    for ka, kb in product(COEFFS, COEFFS):
        record(ka * av + kb * bv, ((ka, an), (kb, bn)))
# 3-term
for (an, av), (bn, bv), (cn, cv) in combinations(mono_items, 3):
    for ka, kb, kc in product(COEFFS, COEFFS, COEFFS):
        record(ka * av + kb * bv + kc * cv, ((ka, an), (kb, bn), (kc, cn)))

covered = [t for t in range(LO, HI + 1) if spell[t]]
n_int = HI - LO + 1
counts = {t: len(spell[t]) for t in range(LO, HI + 1)}
value_multisets = {
    t: len({tuple(sorted((k, mono[m]) for k, m in e)) for e in spell[t]})
    for t in range(LO, HI + 1)
}

print("\n=== R1: promiscuity ===")
print(f"coverage of [{LO},{HI}]: {len(covered)}/{n_int} = {100*len(covered)/n_int:.1f}%")
missing = [t for t in range(LO, HI + 1) if not spell[t]]
print(f"unexpressible integers: {missing if missing else 'none'}")
avg = sum(counts.values()) / n_int
med = sorted(counts.values())[n_int // 2]
print(f"spellings per integer: mean {avg:,.0f}, median {med:,}, "
      f"min {min(counts.values()):,}, max {max(counts.values()):,}")

ladder = [44, 84, 137, 140, 142, 168, 184]
print("\nladder integers — formal spellings (and distinct value-multisets):")
for t in ladder:
    print(f"  {t}: {counts[t]:>7,} formal   ({value_multisets[t]:,} value-distinct)")

print("\npaper-66 formulas verified:")
checks = [
    ("V+E+chi", inv["V"] + inv["E"] + inv["chi"], 44),
    ("2(V+E)", 2 * (inv["V"] + inv["E"]), 84),
    ("7F-D", 7 * inv["F"] - inv["D"], 137),
    ("7F", 7 * inv["F"], 140),
    ("7F+chi", 7 * inv["F"] + inv["chi"], 142),
    ("(F-V)(F+D) = F*F+F*D-V*F-V*D", inv["F"]**2 + inv["F"]*inv["D"] - inv["V"]*inv["F"] - inv["V"]*inv["D"], 184),
]
for s, v, want in checks:
    print(f"  {s} = {v}  {'OK' if v == want else 'MISMATCH'}"
          f"{'' if v != want or s.startswith('(') else '  (in class C)'}")
print("  note: (F-V)(F+D) needs 4 product terms; nearest in-class spellings of 184"
      " exist in the thousands regardless (see count above).")

# ------------------------------------------------- R2: near-integer p-values
print("\n=== R2: near-integer content of the measured exponents ===")
hbar_c = 197.3269804          # MeV fm
mp = 938.27208816             # MeV
mn = 939.56542052             # MeV
GeV = 1e3
MP_GeV = 1.220890e19          # CODATA standard Planck mass, GeV; rel err ~1.1e-5
dMP = 1.1e-5                  # relative

def near_int_report(tag, x, err=None):
    fr = x - math.floor(x)
    d = min(fr, 1 - fr)
    p = 2 * d
    e = f" ± {err:.5f}" if err else ""
    print(f"  {tag}: {x:.5f}{e}   nearest-integer distance {d:.5f}   p = {p:.4f}")
    return p

x44 = math.log(MP_GeV * GeV / mp)
p44 = near_int_report("ln(M_P/m_p)   [standard M_P, proton]", x44, dMP)
x44n = math.log(MP_GeV * GeV / mn)
near_int_report("ln(M_P/m_n)   [standard M_P, neutron]", x44n, dMP)
x44r = x44 - 0.5 * math.log(8 * math.pi)
near_int_report("ln(M̄_P/m_p)  [REDUCED M_P — convention sensitivity]", x44r, dMP)

lP = 1.616255e-35             # m
print("  --- cosmic exponents (definition-dominated; spreads, not errors) ---")
for H0 in (67.4, 73.0):
    R = 2.99792458e5 / H0 * 3.0857e19 * 1e3   # c/H0 in m  (km/s/Mpc -> m)
    near_int_report(f"ln(R_H/l_P)  [Hubble radius, H0={H0}]", math.log(R / lP))
R_ph = 4.4e26
near_int_report("ln(R/l_P)    [particle horizon 4.4e26 m]", math.log(R_ph / lP))
nb = 0.25                      # baryons / m^3
for tag, R in (("Hubble (H0=70)", 2.99792458e5/70*3.0857e22), ("particle horizon", R_ph)):
    Nb = nb * 4 / 3 * math.pi * R**3
    near_int_report(f"ln(N_baryon) [{tag}]", math.log(Nb))

# ------------------------------------------------- R3: the Zoll identity
print("\n=== R3: the system constraint ===")
print("Analytic: with the Zoll closure R = 2GM/(c^2*sqrt(3)) and M = N*m_p,")
print("  ln(R/l_P) = ln N + ln(m_p/M_P) + ln(2/sqrt(3))")
print(f"  i.e. E140 = E184 - E44 + {math.log(2/math.sqrt(3)):.4f}")
print("The identity 184 = 140 + 44 is the closure assumption restated, exact to a")
print(f"residual ln(2/sqrt(3)) = {math.log(2/math.sqrt(3)):.4f}; it can only ever 'bite'")
print("as evidence at a resolution finer than 0.14 in the exponent, which no current")
print("input (R to +/-1, N to +/-2) approaches.")
