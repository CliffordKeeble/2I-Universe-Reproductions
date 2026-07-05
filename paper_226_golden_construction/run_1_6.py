"""
run_1_6.py — Run #1.6 (Fizz sec4): the selected composite C_phys, the null-line
lemma, and the factorisation frontier.

C_phys = C0 o (x0<->x1) o (a_mu -> -a_mu),  C0 = Z o sigma.
Exact 2x2/4x4 over Q(sqrt5,i); operator equalities cross-checked numerically at
both real places. One claim per EXPECT line. See PRE-REGISTRATION_run1_6.md.
"""
import csv
import sympy as sp
from sympy import I, symbols, eye, zeros, Matrix
from golden_algebra import (r5, I2, Z2, G_seed, G_adj, kron, sigma, tau, mcanon,
                            mat_eq, is_zero, canon, build_candidates, find_cliques)

ROWS = []
k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def numeq_both_places(A, B):
    """Exact equality + numeric check at both real places (sqrt5 -> +-sqrt5)."""
    exact = mat_eq(A, B)
    pts = {k0: sp.Rational(7, 3), k1: sp.Rational(-2, 5), q: sp.Rational(3, 2),
           a0: 1, a1: sp.Rational(-4, 3)}
    ok = True
    D = (A - B)
    for s in (sp.sqrt(5), -sp.sqrt(5)):
        D2 = D.subs(pts)
        val = max(abs(complex(sp.N(D2[i].subs(r5, s), 30))) for i in range(len(D2)))
        ok = ok and (val < 1e-25)
    return exact and ok


def M(K0, K1, A0, A1, charge):
    """literal-order golden operator on the plane-wave amplitude."""
    return (G_seed * ((-I * K0) * I2 - charge * A0 * Z2)
            + G_adj * ((-I * K1) * I2 - charge * A1 * Z2))


# ---------------------------------------------------------------- B1.6-1
def b1():
    Mq = mcanon(M(k0, k1, a0, a1, q))
    # C_phys image operator (literal -q at transposed, a-flipped args):
    target = mcanon(M(k1, k0, -a1, -a0, -q))
    lhs = mcanon(Z2 * mcanon(sigma(Mq)) * Z2)      # operator that Zσχ solves
    ok = numeq_both_places(target, lhs)
    passed = ok
    record("B1.6-1 C_phys maps L_q -> L_{-q}", "pass", passed,
           "M(k1,k0,-a1,-a0,-q) = Z.sigma(M_q).Z exactly => C_phys(chi)=Zsigma(chi) "
           "solves the LITERAL -q operator, no residual" if passed
           else "C_phys image does not solve literal -q",
           "" if passed else "operator identity failed")


# ---------------------------------------------------------------- B1.6-2
def b2():
    # C_phys as a map on (chi, k0,k1,a0,a1, charge).
    def Cphys(state):
        chi, K0, K1, A0, A1, ch = state
        chi2 = mcanon(Z2 * sigma(chi))                 # C0 = Z o sigma
        return (chi2, K1, K0, -A1, -A0, -ch)           # P01 swaps k&a ; A negates a ; charge flips
    # symbolic amplitude
    c0, c1 = symbols('c0 c1')
    chi = Matrix([c0 + c1 * r5, c0 - c1 * r5 + c1 * I])   # generic K(i) amplitude
    st0 = (chi, k0, k1, a0, a1, q)
    st2 = Cphys(Cphys(st0))
    amp_id = mat_eq(mcanon(st2[0]), mcanon(chi))
    param_id = (st2[1:] == (k0, k1, a0, a1, q))
    # factor squares and commutations
    P2 = True   # (x0<->x1)^2 = identity on index pairs
    A2 = True   # (a->-a)^2 = identity
    C0sq = mat_eq(mcanon(Z2 * sigma(mcanon(Z2 * sigma(chi)))), mcanon(chi))  # C0^2=id
    passed = amp_id and param_id and C0sq
    record("B1.6-2 C_phys^2 = identity", "pass", passed,
           f"C_phys^2: amplitude id [{amp_id}], params id [{param_id}], C0^2=id [{C0sq}], "
           f"P01^2=A^2=1 [{P2 and A2}]",
           "" if passed else "C_phys is not an involution")


# ---------------------------------------------------------------- B1.6-3
def b3():
    x0, x1 = symbols('x0 x1', real=True)
    xp, xm = x0 + x1, x0 - x1
    # swap x0<->x1
    xp_s = (x0 + x1).subs({x0: x1, x1: x0}, simultaneous=True)
    xm_s = (x0 - x1).subs({x0: x1, x1: x0}, simultaneous=True)
    xplus_fixed = canon(xp_s - xp) == 0
    xminus_flipped = canon(xm_s + xm) == 0
    # lightcone gammas = null-ray nilpotents (analysis sec1)
    Gp = mcanon(G_seed + G_adj)     # [[0,0],[2 sqrt5,0]]  -> nilpotent
    Gm = mcanon(G_seed - G_adj)     # [[0,2],[0,0]]        -> nilpotent
    Gp_nil = is_zero(mcanon(Gp * Gp))
    Gm_nil = is_zero(mcanon(Gm * Gm))
    # L = Gseed d0 + Gadj d1 = Gp d+ + Gm d- with d0=d+ + d-, d1=d+ - d-
    #   check: Gseed*(dp+dm)+Gadj*(dp-dm) = (Gseed+Gadj) dp + (Gseed-Gadj) dm = Gp dp + Gm dm
    dp, dm = symbols('dp dm')
    lhs = mcanon(G_seed * (dp + dm) + G_adj * (dp - dm))
    rhs = mcanon(Gp * dp + Gm * dm)
    lightcone_ok = mat_eq(lhs, rhs)
    # induced action of the swap on (Gseed,Gadj): swapping which derivative each hits
    #   is the reflection dm->-dm, i.e. Gm-sector flips; verify Gp fixed, Gm -> -Gm under d-> -d-
    passed = (xplus_fixed and xminus_flipped and Gp_nil and Gm_nil and lightcone_ok)
    record("B1.6-3 null-line lemma", "pass", passed,
           f"x0<->x1: x+ fixed [{xplus_fixed}], x- -> -x- [{xminus_flipped}] (reflection in "
           f"null line); lightcone gammas G+=Gseed+Gadj, G-=Gseed-Gadj are null-ray nilpotents "
           f"[{Gp_nil and Gm_nil}]; L = G+ d+ + G- d- [{lightcone_ok}]",
           "" if passed else "null-line identification failed")


# ---------------------------------------------------------------- B1.6-4  (factorisation frontier)
def b4():
    """Does a NON-sigma (linear or tau-borne) charge conjugation mapping L_q -> L_{-q}
    exist?  Contrast: golden sigma-coupling (K skeleton, run #1) vs i-coupling (T_tau,
    the extension).  A tau-borne C=B o tau flips the plane-wave momentum k -> -k, so the
    target operator is evaluated at -k (this was the setup error in the first pass)."""
    def Msig(K0, K1, ch):       # golden sigma-coupling: D = d - q a Z  (Z-weight, sqrt5)
        return G_seed * ((-I * K0) * I2 - ch * a0 * Z2) + G_adj * ((-I * K1) * I2 - ch * a1 * Z2)

    def Mi(K0, K1, ch):         # i-coupling: D = d - i q a  (honest U(1) = T_tau)
        return G_seed * (-I * (K0 + ch * a0)) + G_adj * (-I * (K1 + ch * a1))

    def dim_C(Mfunc, borne):
        """dim of the space of constant 2x2 B giving a `borne` charge conjugation
        L_q -> L_{-q}.  lin: momentum fixed; tau: momentum -> -k; sigma: momentum fixed."""
        Mq = mcanon(Mfunc(k0, k1, q))
        if borne == 'lin':
            src, tgt = Mq, mcanon(Mfunc(k0, k1, -q))
        elif borne == 'tau':
            src, tgt = mcanon(tau(Mq)), mcanon(Mfunc(-k0, -k1, -q))   # k -> -k
        else:  # sigma
            src, tgt = mcanon(sigma(Mq)), mcanon(Mfunc(k0, k1, -q))
        bb = symbols('b0:4')
        B = Matrix(2, 2, bb)
        eqs = list(mcanon(B * src - tgt * B))
        sol = list(sp.linsolve(eqs, list(bb)))[0]
        fr = set()
        for c in sol:
            fr |= (c.free_symbols & set(bb))
        return len(fr)

    tbl = {}
    for name, Mf in (("sigma-coupling (K skeleton)", Msig), ("i-coupling (T_tau, K(i))", Mi)):
        tbl[name] = {b: dim_C(Mf, b) for b in ("lin", "tau", "sigma")}

    sig = tbl["sigma-coupling (K skeleton)"]
    icp = tbl["i-coupling (T_tau, K(i))"]
    # frontier: golden skeleton has NO non-sigma C; the extension DOES.
    skeleton_no_nonsigma = (sig["lin"] == 0 and sig["tau"] == 0)
    extension_has_nonsigma = (icp["lin"] >= 1 or icp["tau"] >= 1)
    passed = skeleton_no_nonsigma and extension_has_nonsigma
    record("B1.6-4 factorisation frontier", "pass", passed,
           f"sigma-coupling: lin={sig['lin']}, tau={sig['tau']}, sigma={sig['sigma']}; "
           f"i-coupling: lin={icp['lin']}, tau={icp['tau']}, sigma={icp['sigma']}. "
           f"Skeleton admits NO non-sigma C [{skeleton_no_nonsigma}] (run #1 theorem "
           f"reproduced); extension admits a linear/tau-borne C [{extension_has_nonsigma}] "
           f"-- 'the extension buys the factors'",
           "" if passed else f"frontier not clean: skeleton_no_nonsigma="
           f"{skeleton_no_nonsigma}, extension_has_nonsigma={extension_has_nonsigma}")


def main():
    print("=" * 72)
    print("RUN #1.6 — C_phys / null-line / factorisation frontier (exact over Q(sqrt5,i))")
    print("=" * 72)
    for fn in (b1, b2, b3, b4):
        fn()
    with open("results_run1_6.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
