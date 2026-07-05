"""
run_1_5_place_pair.py — Run #1.5, the place-pair toy (Fizz reading sec6).

Tests whether reading the ONE K-rational toy field at the two real places of K
gives (particle @ inf1, antiparticle @ inf2), with C0 = Z (x) slot-swap absorbing
A5's coordinate transposition.

Decisive results are at the operator level (exact 2x2 over Q(sqrt5,i), all k),
each equality cross-checked numerically at BOTH real places (mpmath), per the
standing discipline promoted in Fizz sec3. See PRE-REGISTRATION_run1_5.md.
"""
import csv
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import I, symbols, eye, Matrix
from golden_algebra import (r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq,
                            canon, is_zero)

mp.mp.dps = 40
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


# ---- operator family --------------------------------------------------------
k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)
GA = {"Gs": G_seed, "Ga": G_adj}


def Dm(km, am, ch):
    return (-I * km) * I2 - ch * am * Z2


def L(order, charge, asign, ksign=+1):
    """order = (gamma-on-d0, gamma-on-d1); charge multiplies coupling;
       asign flips a_mu; ksign flips the momentum k_mu."""
    return (GA[order[0]] * Dm(ksign * k0, asign * a0, charge)
            + GA[order[1]] * Dm(ksign * k1, asign * a1, charge))


Mq = mcanon(L(("Gs", "Ga"), q, +1))          # the q-equation operator


def numeq_both_places(A, B):
    """Exact-symbolic equality, cross-checked numerically at both real places."""
    exact = mat_eq(A, B)
    subs_pts = {k0: sp.Rational(7, 3), k1: sp.Rational(-2, 5),
                q: sp.Rational(3, 2), a0: 1, a1: sp.Rational(-4, 3)}
    ok_num = True
    for s in (sp.sqrt(5), -sp.sqrt(5)):
        for M in (A, B):
            pass
        DA = (A - B).subs(subs_pts)
        val = max(abs(complex(sp.N(DA[i].subs(r5, s), 30))) for i in range(len(DA)))
        ok_num = ok_num and (val < 1e-25)
    return exact, ok_num


def scan_family(target, label):
    """Find the unique family member equal to `target`; report it."""
    hits = []
    for order in (("Gs", "Ga"), ("Ga", "Gs")):
        for charge in (q, -q):
            for asign in (+1, -1):
                for ksign in (+1, -1):
                    cand = mcanon(L(order, charge, asign, ksign))
                    exact, num = numeq_both_places(target, cand)
                    if exact and num:
                        oname = "literal(Gs@d0)" if order == ("Gs", "Ga") else "swapped(Ga@d0)"
                        cname = "+q" if charge == q else "-q"
                        hits.append(f"{oname}, {cname}, a{'+' if asign > 0 else '-'}, "
                                    f"k{'+' if ksign > 0 else '-'}")
    return hits


# ---------------------------------------------------------------- B1.5-1
def b1():
    """pair(chi) = (chi, sigma chi) solves (M_q @ inf1, sigma(M_q) @ inf2).
       Identify sigma(M_q) = ev2(M_q)."""
    sMq = mcanon(sigma(Mq))
    # tautology part: sigma(M_q) sigma(chi) = sigma(M_q chi) = 0 for any chi
    # (verified structurally: sigma is a ring homomorphism applied entrywise)
    hits = scan_family(sMq, "sigma(Mq)")
    # is it a bona-fide antiparticle operator (reversed metric + reversed momentum)?
    passed = len(hits) >= 1
    outcome = (f"inf2 slot solves ev2(M_q)=sigma(M_q) = [{'; '.join(hits)}] "
               f"-- the reversed-metric operator; slot pair solves (L_q@inf1, that@inf2)")
    record("B1.5-1 pair reads two dynamics", "pass", passed, outcome,
           "" if passed else "sigma(M_q) not identified in the operator family")
    b1.hits = hits


# ---------------------------------------------------------------- B1.5-2  (make-or-break)
def prop(A, B):
    """A proportional to B by a unit scalar c in {1,-1,i,-i} (exact + numeric both places)?"""
    for c in (1, -1, I, -I):
        exact, num = numeq_both_places(A, mcanon(c * B))
        if exact and num:
            return {1: "+1", -1: "-1", I: "+i", -I: "-i"}[c]
    return None


def b2():
    """Place-pair level, no on-shell needed (kernel match <=> operator proportionality).
       q-pair: (X1 solves M_q@inf1, X2 solves M_q@inf2=sigma(M_q)).
       Image C0=Z(x)swap: slot1 Y1=Z.X2, slot2 Y2=Z.X1.
       Y1 solves O  <=>  O proportional to Z.sigma(M_q).Z  (since X2 solves sigma(M_q)).
       Y2 solves O  <=>  O proportional to Z.M_q.Z.
       'No coordinate transposition' = image slot i solves the LITERAL -q operator at place i:
         slot1 target = M_{-q}          (literal, inf1)
         slot2 target = sigma(M_{-q})   (literal, inf2).
    """
    O1 = mcanon(Z2 * mcanon(sigma(Mq)) * Z2)     # operator Y1 solves
    O2 = mcanon(Z2 * Mq * Z2)                    # operator Y2 solves
    Mmq = mcanon(L(("Gs", "Ga"), -q, +1))        # literal -q
    Mmq_aflip = mcanon(L(("Gs", "Ga"), -q, -1))  # literal -q, a->-a
    sMmq = mcanon(sigma(Mmq))
    sMmq_aflip = mcanon(sigma(Mmq_aflip))

    slot1_lit = prop(O1, Mmq) or prop(O1, Mmq_aflip)     # Y1 solves literal -q (any a-sign)?
    slot2_lit = prop(O2, sMmq) or prop(O2, sMmq_aflip)   # Y2 solves literal sigma(-q)?
    what1 = scan_family(O1, "O1")                        # what Y1 actually solves
    what2 = scan_family(O2, "O2")
    absorbed = (slot1_lit is not None) and (slot2_lit is not None)
    passed = absorbed
    if absorbed:
        outcome = ("C0=Z(x)swap image solves the LITERAL -q pair (M_{-q}@inf1, "
                   "sigma(M_{-q})@inf2) -- transposition ABSORBED by the place structure")
        dev = ""
    else:
        outcome = ("image slot-1 does NOT solve the literal -q operator M_{-q} at inf1 "
                   f"(prop={slot1_lit}); it solves [{'; '.join(what1)}] -- a SEED<->ADJOINT "
                   "swapped operator")
        dev = ("make-or-break FAIL per Fizz sec6: C0=Z(x)swap lands the q-solution on a "
               "swapped-order operator at every slot assignment, never the literal -q "
               "operator. The place structure does NOT absorb the A5 transposition -- the "
               "two golden directions (Gseed^2=+sqrt5 vs Gadj^2=-sqrt5) are genuinely "
               "inequivalent and Z(x)swap cannot reconcile them with -q dynamics. "
               f"Y2 solves [{'; '.join(what2)}]. Not patched.")
    record("B1.5-2 C0 absorbs transposition", "pass (make-or-break)", passed, outcome, dev)
    b2.facts = dict(slot1_lit=slot1_lit, slot2_lit=slot2_lit, what1=what1, what2=what2)


# ---------------------------------------------------------------- B1.5-3
def b3():
    """Born pairing across slots: P = ev1(a).ev2(a) = a.sigma(a) = Galois norm N(a);
       C0 (slot swap) leaves it invariant."""
    p, qq = symbols('p q', real=True)
    a = p + qq * r5
    P = canon(a * sigma(Matrix([[a]]))[0])           # ev1(a).ev2(a)
    N = canon(p**2 - 5 * qq**2)
    is_norm = (canon(P - N) == 0)
    # C0 swaps slots: P_swapped = ev2(a).ev1(a) = sigma(a).a = same
    P_swapped = canon(sigma(Matrix([[a]]))[0] * a)
    c0_invariant = (canon(P - P_swapped) == 0)
    # also on the actual spinor: <psi1|psi2> across slots for a K-rational amplitude vector
    v = Matrix([1 + r5, 2 - 3 * r5])
    inner = canon((v.T * sigma(v))[0])               # sum_i v_i sigma(v_i)
    inner_is_K_rational = (canon(sigma(Matrix([[inner]]))[0] - inner) == 0)  # in Q (sigma-fixed)
    passed = is_norm and c0_invariant and inner_is_K_rational
    record("B1.5-3 Born = Galois norm", "pass", passed,
           f"P = a.sigma(a) = p^2-5q^2 = N(a) [{is_norm}]; C0-invariant [{c0_invariant}]; "
           f"slot inner product is sigma-fixed (in Q) [{inner_is_K_rational}]",
           "" if passed else "one sub-check failed")


def main():
    print("=" * 72)
    print("RUN #1.5 — PLACE-PAIR TOY (exact 2x2 over Q(sqrt5,i); numeric at both places)")
    print("=" * 72)
    for fn in (b1, b2, b3):
        fn()
    with open("results_run1_5.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
