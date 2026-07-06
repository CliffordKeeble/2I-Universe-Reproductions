"""
rbench.py — Part C, the reconciliation R-bench (R1-R6).

Demonstrates the "one spine" reconciliation of fizz_reconciliation_seam_paper_226.md:
the internal (cold) construction is the even-sector shadow of the geometric (briefed) one
under the embedding  delta -> Gamma_seed*Gamma_adj = sqrt5*sigma3  (so j -> Z).

Exact over Q(sqrt5,i); operator equalities cross-checked numerically at both real places.
One claim per EXPECT. See PRE-REGISTRATION_rbench.md. Ruling 1: outcome-vs-EXPECT only.
"""
import sys
import csv
import sympy as sp
from sympy import I, symbols, Function, diff, Matrix, cosh, sinh, exp
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq, canon

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROWS = []
k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


def numeq_both_places(A, B):
    """Exact + numeric at both real places (arbiter over the golden field)."""
    ex = mat_eq(A, B)
    pts = {k0: sp.Rational(7, 3), k1: sp.Rational(-2, 5), q: sp.Rational(3, 2),
           a0: 1, a1: sp.Rational(-4, 3)}
    ok = True
    D = (A - B).subs(pts)
    for s in (sp.sqrt(5), -sp.sqrt(5)):
        v = max(abs(complex(sp.N(D[i].subs(r5, s), 30))) for i in range(len(D)))
        ok = ok and (v < 1e-25)
    return ex and ok


def prop(A, B):
    """A proportional to B by a unit scalar in {1,-1,i,-i}? return the scalar-name or None."""
    for c, nm in ((1, "+1"), (-1, "-1"), (I, "+i"), (-I, "-i")):
        if numeq_both_places(A, mcanon(c * B)):
            return nm
    return None


# ------------------------------------------------------------------ embedding
def emb(p, val_q):
    """delta -> sqrt5*Z:  a = p + q*delta  ->  p*I + q*sqrt5*Z."""
    return p * I2 + val_q * r5 * Z2


# ---------------------------------------------------------------- R1
def check_r1():
    p1, q1, p2, q2 = symbols('p1 q1 p2 q2', real=True)
    ab_p, ab_q = p1 * p2 + 5 * q1 * q2, p1 * q2 + q1 * p2      # (p1+q1 d)(p2+q2 d), d^2=5
    homo = mat_eq(emb(ab_p, ab_q), emb(p1, q1) * emb(p2, q2))
    unital = mat_eq(emb(1, 0), I2)
    j_Z = mat_eq(emb(0, 1 / r5), Z2)                          # j = delta/sqrt5 -> Z
    ep, em = (I2 + Z2) / 2, (I2 - Z2) / 2
    ep_ok = mat_eq(emb(sp.Rational(1, 2), 1 / (2 * r5)), ep)  # e+ = (1+j)/2 -> (I+Z)/2
    em_ok = mat_eq(emb(sp.Rational(1, 2), -1 / (2 * r5)), em)
    idem = (mat_eq(ep * ep, ep) and mat_eq(em * em, em) and mat_eq(ep * em, sp.zeros(2))
            and mat_eq(ep + em, I2))
    passed = homo and unital and j_Z and ep_ok and em_ok and idem
    record("R1 embedding delta->sqrt5*sigma3", "unital K-algebra embedding; j->Z; e± idempotents",
           passed,
           "delta->sqrt5*Z is a unital K-algebra homomorphism (a*b preserved); j->Z; "
           "e±->(I±Z)/2 orthogonal complete idempotents")


# ---------------------------------------------------------------- R2
def check_r2():
    th = symbols('theta', real=True)
    phase_emb = cosh(th) * I2 + sinh(th) * Z2                 # e^{th j} -> e^{th Z}
    d00 = csimp((phase_emb[0, 0]).rewrite(exp) - exp(th))
    d11 = csimp((phase_emb[1, 1]).rewrite(exp) - exp(-th))
    diag_ok = (d00 == 0 and d11 == 0)                        # = diag(e^th, e^-th)
    prod_one = csimp(phase_emb[0, 0] * phase_emb[1, 1] - 1) == 0   # u * u^-1 = 1
    # reproduce A4: (u, u^-1) weights preserve the pairing P = s*a; (u,u) does not
    u, s, aa = symbols('u s a', nonzero=True)
    pres = canon((u * s) * (u**-1 * aa) - s * aa) == 0
    notpres = canon((u * s) * (u * aa) - s * aa) != 0
    passed = diag_ok and prod_one and pres and notpres
    record("R2 transported gauge action", "e^{th j} = diag(e^th,e^-th); reproduces A4 weights",
           passed,
           "e^{th j} -> diag(e^th, e^-th) (weights u, u^-1, product 1); (u,u^-1) preserves "
           "the pairing s*a, (u,u) does not — Run #1 A4 reproduced")


# ---------------------------------------------------------------- R3
def check_r3():
    p, val_q = symbols('p q', real=True)                     # sqrt5-free coefficients
    lhs = sigma(emb(p, val_q))                               # entrywise sigma on embedded copy
    rhs = emb(p, -val_q)                                     # 𝕂's Galois sigma: delta -> -delta
    passed = mat_eq(lhs, rhs)
    record("R3 entrywise sigma = Galois sigma", "entrywise sigma|embedded 𝕂 = delta->-delta",
           passed,
           "entrywise sigma(p*I + q*sqrt5*Z) = p*I - q*sqrt5*Z = emb(p - q*delta) — matches "
           "𝕂's Galois sigma on the embedded copy")


# ---------------------------------------------------------------- R4
def check_r4():
    Gp = mcanon(G_seed + G_adj)                              # null-frame gammas
    Gm = mcanon(G_seed - G_adj)

    def C0(X):
        return mcanon(Z2 * sigma(X) * Z2)                    # C0 = Z sigma(.) Z^-1

    s_gp = mat_eq(sigma(Gp), -Gp)                            # sigma(G+) = -G+
    s_gm = mat_eq(sigma(Gm), Gm)                             # sigma(G-) = +G-
    c_gp = mat_eq(C0(Gp), Gp)                                # C0(G+) = +G+
    c_gm = mat_eq(C0(Gm), -Gm)                               # C0(G-) = -G-
    passed = s_gp and s_gm and c_gp and c_gm
    record("R4 null-ray signs", "sigma(G±)=∓? ; C0(G+)=+G+, C0(G-)=-G-", passed,
           f"sigma(G+)=-G+ [{s_gp}], sigma(G-)=+G- [{s_gm}]; C0(G+)=+G+ [{c_gp}], "
           f"C0(G-)=-G- [{c_gm}] — sigma is a chiral sign on one null ray, C0 the compensator")


# ---------------------------------------------------------------- R5  (even-sector control)
def check_r5():
    x, g = symbols('x g', real=True)
    u = Function('u', real=True)(x)
    v = Function('v', real=True)(x)
    Ax = Function('Ax', real=True)(x)

    def reg2(sc, tc):
        return Matrix([[sc, 5 * tc], [tc, sc]])             # regular rep of 𝕂

    JJ = reg2(0, 1 / r5)                                     # j
    af = reg2(u, v)

    def ddx(M):
        return M.applyfunc(lambda e: diff(e, x))

    def sig2(M):
        return Z2 * M * Z2                                   # sigma in the reg rep

    def Dq(M, A):
        return ddx(M) - g * A * (JJ * M)                    # D = d - g A j (charge in A)

    # C = sigma o (a-flip): a -> sigma(a), A -> -A. Solution map iff D_{-A}(sigma a)=sigma(D_A a).
    lhs = Dq(sig2(af), -Ax)
    rhs = sig2(Dq(af, Ax))
    clean = all(csimp((lhs - rhs)[i]) == 0 for i in range(4))
    # confirm NO coordinate factor: the identity uses no x-transformation (x untouched)
    passed = clean
    record("R5 even-sector control", "C=sigma o (a-flip) maps q->-q with no coordinate factor",
           passed,
           "internal scalar EoM: D_{-A}(sigma a) = sigma(D_A a) exactly (uses sigma(j)=-j); "
           "so a q-solution maps to a (-q)-solution with the gauge-field flip alone — NO "
           "coordinate/reflection factor. Even sector is clean")


# ---------------------------------------------------------------- R6  (the dichotomy)
def check_r6():
    GA = {"Gs": G_seed, "Ga": G_adj}

    def Dm(km, am, ch):
        return (-I * km) * I2 - ch * am * Z2

    def L(order, charge, asign, ksign=1):
        return (GA[order[0]] * Dm(ksign * k0, asign * a0, charge)
                + GA[order[1]] * Dm(ksign * k1, asign * a1, charge))

    Mq = mcanon(L(("Gs", "Ga"), q, +1))
    O1 = mcanon(Z2 * sigma(Mq) * Z2)                         # what sigma o (a-flip) image solves

    # half A: WITHOUT the reflection, O1 is NOT the literal -q operator (either a-sign)
    lit = prop(O1, mcanon(L(("Gs", "Ga"), -q, +1))) or prop(O1, mcanon(L(("Gs", "Ga"), -q, -1)))
    fails_literal = (lit is None)
    # ...but it IS a SWAPPED-order (seed<->adjoint reflected) operator -> reflection is needed
    swapped_hit = None
    for asign in (+1, -1):
        for ksign in (+1, -1):
            nm = prop(O1, mcanon(L(("Ga", "Gs"), -q, asign, ksign)))
            if nm:
                swapped_hit = f"swapped-order, a{'+' if asign > 0 else '-'}, k{'+' if ksign > 0 else '-'} ({nm})"
                break
        if swapped_hit:
            break
    half_a = fails_literal and (swapped_hit is not None)

    # half B: the even-projected observable (Born pairing = Galois norm) is swap/C0-invariant
    pp, qq = symbols('pp qq', real=True)
    aK = pp + qq * r5
    Pnorm = canon(aK * sigma(Matrix([[aK]]))[0])            # a*sigma(a)
    norm_is_N = (canon(Pnorm - (pp**2 - 5 * qq**2)) == 0)
    Pswap = canon(sigma(Matrix([[aK]]))[0] * aK)           # swap slots
    norm_invariant = (canon(Pnorm - Pswap) == 0)
    half_b = norm_is_N and norm_invariant

    passed = half_a and half_b
    record("R6 the dichotomy", "sigma o (a-flip) fails on odd sector w/o reflection; ok on even",
           passed,
           f"geometric: O1=Z.sigma(Mq).Z is NOT the literal -q operator [fails_literal="
           f"{fails_literal}] but IS a reflected/swapped one [{swapped_hit}] — reflection "
           f"required (B1.5-2 reconfirmed); even-projected observable a.sigma(a)=N(a) is "
           f"swap-invariant [{half_b}] — succeeds. Dichotomy holds on both halves")


def main():
    print("=" * 72)
    print("PART C — RECONCILIATION R-BENCH (R1-R6), exact over Q(sqrt5,i), both real places")
    print("=" * 72)
    for fn in (check_r1, check_r2, check_r3, check_r4, check_r5, check_r6):
        fn()
    with open("results_rbench.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for row in ROWS if row[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
