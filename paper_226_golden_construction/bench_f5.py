"""
bench_f5.py — F5 aligner sigma-type crux (step 1 + structural data) + Mr A's rechecks C1, C6.

Source of record: mr_code_brief_aligner_sigma_type_F5.md (Fizz, 10 Jul 2026), from Mr A's
three-star review. Pre-registration: PRE-REGISTRATION_f5.md (committed 683f9c2, before this verdict).

WHAT THIS SCRIPT DOES (self-contained, exact over K(i)):
  * F5 step 1  : solve the (+,+) aligner family (A Gs, A Ga both symmetric); verify Fizz's
                 antidiagonal family independently.
  * Structural : neutral data on the aligner (nondegeneracy, sigma-action, cross-component
                 bilinear, the member that equals Gamma_s) -- CONJECTURED inputs for the legality
                 test, NOT a legality verdict.
  * C1         : independent re-derivation of the dispersion (hand determinant + explicit kernel
                 vector u(k)), two-place numeric cross-check.
  * C6         : independent construction of the tau-halved symplectic matrix -> rank 4, det 625.

WHAT THIS SCRIPT DOES NOT DO:
  * F5 step 2 (sigma-legality verdict) is BLOCKED. The 191/226 displayed dressings that define the
    legality criterion are uncommitted (Downloads only). The brief forbids inferring them. Requisition
    raised in findings_f5.md; the verdict is not computed here.
"""
import os
import sys
import csv
import sympy as sp
from sympy import sqrt, I, symbols, Matrix, diag, Rational
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

Gs, Ga, Z = G_seed, G_adj, Z2
A0 = diag(r5, 1)                       # the original F2a sigma-pairing (gives (+,-))
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else ("FAIL" if passed is False else "INFO")
    ROWS.append((item, expect, f"{tag} - {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


def is_symmetric(M):
    return mat_eq(M, M.T)


# ================================================================ F5 step 1: aligner family
def f5_step1():
    a, b, c, d = symbols('a b c d')          # entries of A in K(i)
    A = Matrix([[a, b], [c, d]])
    AGs = mcanon(A * Gs)
    AGa = mcanon(A * Ga)
    # symmetry of a 2x2 M  <=>  M[1,0] - M[0,1] = 0  (one scalar eqn each)
    eqs = [csimp(AGs[1, 0] - AGs[0, 1]), csimp(AGa[1, 0] - AGa[0, 1])]
    sol = sp.linsolve(eqs, [a, b, c, d])
    (sa, sb, sc, sd), = sol
    free = (sa.free_symbols | sb.free_symbols | sc.free_symbols | sd.free_symbols) & {a, b, c, d}
    dim = len(free)
    # antidiagonal?  a=d=0, b,c free
    antidiag = (csimp(sa) == 0 and csimp(sd) == 0 and {sb, sc} == {b, c})
    # AGs = diag(b r5, c), AGa = diag(b r5, -c) on the solution
    Asol = Matrix([[sa, sb], [sc, sd]])
    AGs_s = mcanon(Asol * Gs)
    AGa_s = mcanon(Asol * Ga)
    shape_ok = (mat_eq(AGs_s, diag(b * r5, c)) and mat_eq(AGa_s, diag(b * r5, -c)))
    both_sym = is_symmetric(AGs_s) and is_symmetric(AGa_s)
    passed = (dim == 2) and antidiag and shape_ok and both_sym
    record("F5.1 aligner family",
           "antidiagonal A=[[0,b],[c,0]]; A.Gs=diag(b r5,c), A.Ga=diag(b r5,-c) both symmetric",
           passed,
           f"solution family dim={dim}, A={tuple(Asol)} (a=d=0, b,c free) [antidiag={antidiag}]; "
           f"A.Gs=diag(b*r5,c), A.Ga=diag(b*r5,-c) [shape={shape_ok}], both symmetric [{both_sym}] "
           f"-> the (+,+) alignment exists, as Taussky-Zassenhaus guarantees")
    return Asol, (a, b, c, d)


# ================================================================ structural (neutral, CONJECTURED)
def structural(Asol, syms):
    a, b, c, d = syms
    A = Asol                                   # [[0,b],[c,0]]
    # nondegeneracy
    detA = csimp(A.det())                       # = -b c
    # sigma action (entrywise): sigma(A) = [[0,sigma(b)],[sigma(c),0]]; = A iff b,c in Q(i)
    sigA = mcanon(sigma(A))
    sig_fixes_if_Qi = mat_eq(sigA, A.subs({b: b, c: c}))  # structurally equal form; holds if b,c sqrt5-free
    # cross-component bilinear:  (sigma psi)^T A chi  vs the same-component A0 pairing
    p1, p2, q1, q2 = symbols('p1 p2 q1 q2')     # generic spinor components (as opaque symbols)
    psi = Matrix([p1, p2]); chi = Matrix([q1, q2])
    # use opaque sigma via a formal marker: pair structure only (coefficients b,c)
    bil_A = sp.expand((Matrix([[p1, p2]]) * A * chi)[0])     # b p1 q2 + c p2 q1  (cross)
    bil_A0 = sp.expand((Matrix([[p1, p2]]) * A0 * chi)[0])   # r5 p1 q1 + p2 q2   (same-component)
    cross = (bil_A == sp.expand(b * p1 * q2 + c * p2 * q1))
    # the member A(b=1,c=r5) equals Gamma_s exactly
    A_is_Gs = mat_eq(A.subs({b: 1, c: r5}), Gs)
    # the aligner factors as component-swap K2 . diagonal
    K2 = Matrix([[0, 1], [1, 0]])
    factors = mat_eq(A, mcanon(K2 * diag(c, b)))            # [[0,b],[c,0]] = K2 . diag(c,b)
    record("F5.1 structural data (CONJECTURED, not banked)",
           "neutral inputs for the legality test; shape reported, legality NOT adjudicated",
           None,
           f"det A = {detA} (nondegenerate iff b,c != 0); (sigma psi)^T A chi = b*p1*q2 + c*p2*q1 "
           f"CROSS-component [{cross}] vs A0 pairing r5*p1*q1 + p2*q2 same-component; "
           f"member A(b=1,c=r5) = Gamma_s [{A_is_Gs}]; A = K2 . diag(c,b) i.e. component-swap x diagonal "
           f"[{factors}] -> off-diagonal 'place-swap' shape; legality decided by requisitioned axioms, "
           f"NOT by this shape")


# ================================================================ C1: dispersion recheck (independent)
def C1():
    E, k, m = symbols('E k m', real=True, positive=True)
    ok_all, detail = True, []
    for sname, s in (("s=1", sp.Integer(1)), ("s=i", I)):
        M = mcanon(-I * E * Gs + I * k * Ga + m * s * Z)
        # (1) hand determinant, NOT sympy .det()
        det_hand = csimp(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0])
        det_target = csimp(-m**2 * s**2 - r5 * (k**2 - E**2))
        det_ok = (csimp(det_hand - det_target) == 0)
        # (2) explicit kernel vector u(k) = (i(E+k), m s), E on-shell = sqrt(k^2 + m^2 s^2/r5)
        Msq = m**2 * s**2 / r5
        Eshell = sqrt(k**2 + Msq)
        u = Matrix([I * (Eshell + k), m * s])
        Mu = mcanon(M.subs(E, Eshell) * u)
        kernel_ok = (csimp(Mu[0]) == 0 and csimp(Mu[1]) == 0)
        # (3) two-place numeric cross-check (float only here, per pre-reg)
        subs0 = {k: Rational(2), m: Rational(3)}
        num_ok = True
        for placelabel, r5val in (("inf1(+r5)", sqrt(5).evalf()), ("inf2(-r5)", -sqrt(5).evalf())):
            Msq_n = (Rational(3)**2 * s**2 / r5val)
            E_n = sp.sqrt(Rational(2)**2 + Msq_n)
            Mn = (-I * E_n * Gs + I * k * Ga + m * s * Z).subs(subs0)
            Mn = Mn.subs(r5, r5val)
            un = Matrix([I * (E_n + 2), m * s]).subs(subs0)
            resid = complex((Mn * un).norm().evalf())
            num_ok = num_ok and (abs(resid) < 1e-9)
        this_ok = det_ok and kernel_ok and num_ok
        ok_all = ok_all and this_ok
        detail.append(f"{sname}: hand-det=-m^2 s^2 -r5(k^2-E^2) [{det_ok}]; u(k) is exact null vector "
                      f"Mu=0 on-shell [{kernel_ok}]; two-place numeric Mu~0 [{num_ok}]")
    record("C1 dispersion recheck (independent)",
           "E^2=k^2+m^2 s^2/r5; sign from Z^2=+I",
           ok_all, "; ".join(detail))


# ================================================================ C6: symplectic recheck (independent)
def C6():
    # tau-constrained basis vectors (K(i) 2-spinors): f1,g1,f2,g2 directions
    v = [Matrix([1, 0]), Matrix([I * r5, 0]), Matrix([0, 1]), Matrix([0, I * r5])]
    names = ['f1', 'g1', 'f2', 'g2']
    AGs = A0 * Gs                                   # the d0 (Gamma_s) direction pairing, A=diag(r5,1)
    # direct pairing-matrix construction (no Lagrangian diff): P[j,k] = Im[ (sigma v_j)^T (A Gs) v_k ]
    P = sp.zeros(4, 4)
    for j in range(4):
        for kk in range(4):
            val = (sigma(v[j]).T * AGs * v[kk])[0]
            P[j, kk] = sp.simplify(sp.im(sp.expand(val)))
    Omega = sp.simplify((P - P.T) / 2)              # antisymmetric (symplectic) part
    rank = Omega.rank()
    detO = sp.simplify(Omega.det())
    is5_4 = (detO == 625) and (625 == 5**4)
    passed = (rank == 4) and (detO == 625)
    record("C6 symplectic recheck (independent)",
           "tau-halved symplectic matrix: rank 4, det 625 (=5^4)",
           passed,
           f"direct pairing-matrix (basis {names}) -> antisymmetric Omega rank={rank}/4, "
           f"det={detO}=5^4 [{is5_4}]; rank 4 is the basis-independent nondegeneracy claim, "
           f"the value 625 reflects the sqrt5-normalisation of the q-components")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("=" * 72)
    print("BENCH F5 -- aligner sigma-type (step 1 + structural) + rechecks C1, C6")
    print("exact over Q(sqrt5,i); F5 step 2 (legality) BLOCKED on requisition")
    print("=" * 72)
    print("\n[F5 step 1 -- the aligner family]")
    Asol, syms = f5_step1()
    structural(Asol, syms)
    print("\n[C1 -- dispersion recheck, independent route]")
    C1()
    print("\n[C6 -- symplectic recheck, independent construction]")
    C6()
    with open("results_f5.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    ninfo = sum(1 for r in ROWS if r[2].startswith("INFO"))
    print("\n" + "=" * 72)
    print(f"SUMMARY: {npass} PASS, {ninfo} INFO (structural) of {len(ROWS)} items")
    print("F5 step 2 (sigma-legality): NOT RUN -- requisition raised (see findings_f5.md)")
    print("=" * 72)


if __name__ == "__main__":
    main()
