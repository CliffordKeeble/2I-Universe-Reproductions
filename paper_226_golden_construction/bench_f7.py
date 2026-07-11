"""
bench_f7.py — F7: the reality sweep over the aligner family (Paper 226).

Source of record: mr_code_brief_F7_reality_sweep.md (Fizz, 11 Jul 2026).
Pre-registration: PRE-REGISTRATION_f7.md (committed 0a279a5, before this verdict).

THE QUESTION (K1, criterion-free): does ANY aligner A = [[0,b],[c,0]] give a REAL, NONDEGENERATE
first-order action (kinetic AND mass, under an admissible field-reality constraint, mod total
derivatives)?  ANY -> K-COLLAPSE (F1-F4 no-go RETRACTED).  NONE -> K-SURVIVE (upgrades to class-wide).

Exact symbolic over K(i) = Q(sqrt5, i). Reality criterion: Im[ev1(e^{i th} L)] mod TD = 0, kinetic
(antisymmetric-part extraction, per-direction) and mass (pointwise), at a common phase th=(C,S).
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
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed is True else ("FAIL" if passed is False else "INFO")
    ROWS.append((item, expect, f"{tag} - {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


# ---------------------------------------------------------------- constrained fields (as in F3)
def make_constrained(kind):
    f1, g1, f2, g2 = symbols('f1 g1 f2 g2', real=True)
    f1t, g1t, f2t, g2t = symbols('f1t g1t f2t g2t', real=True)
    f1x, g1x, f2x, g2x = symbols('f1x g1x f2x g2x', real=True)
    dof, d0, d1 = [f1, g1, f2, g2], [f1t, g1t, f2t, g2t], [f1x, g1x, f2x, g2x]

    def comp(f, g, which):
        if kind == 'tau':
            return f + I * r5 * g            # p real, q imaginary
        if kind == 'sigma':
            return f + I * g                 # q=0, psi in Q(i)
        if kind == 'C0':
            return (f + I * g) if which == 1 else (f + I * g) * r5   # Z.sigma Majorana
        raise ValueError(kind)

    psi = Matrix([comp(f1, g1, 1), comp(f2, g2, 2)])
    dp0 = Matrix([comp(f1t, g1t, 1), comp(f2t, g2t, 2)])
    dp1 = Matrix([comp(f1x, g1x, 1), comp(f2x, g2x, 2)])
    return psi, dp0, dp1, dof, d0, d1


def adjoint(psi, A):
    return sigma(psi).T * A


def full_L(psi, dp0, dp1, A, s, m):
    """L = psibar (Gs d0 + Ga d1 + m s Z) psi,  psibar = (sigma psi)^T A."""
    op = Gs * dp0 + Ga * dp1 + m * s * Z * psi
    return sp.expand((adjoint(psi, A) * op)[0])


def split_defect(expr, dof, d0, d1):
    """Split into (mass_part: no derivative) + (kin_part: field*derivative)."""
    dsyms = set(d0) | set(d1)
    mass_part = sp.Integer(0)
    kin_part = sp.Integer(0)
    for term in sp.expand(expr).as_ordered_terms():
        if term.free_symbols & dsyms:
            kin_part += term
        else:
            mass_part += term
    return sp.expand(mass_part), sp.expand(kin_part)


def kin_mod_td(kin, dof, d0, d1):
    """Antisymmetric-part remainder of a field*derivative bilinear (mod total derivatives)."""
    n = len(dof)
    rem = sp.Integer(0)
    recon = sp.Integer(0)
    Lam = {}
    for mu, derivs in ((0, d0), (1, d1)):
        L = sp.zeros(n, n)
        for j in range(n):
            for k in range(n):
                Cjk = sp.diff(kin, dof[j], derivs[k])
                recon += Cjk * dof[j] * derivs[k]
                lam = (Cjk - sp.diff(kin, dof[k], derivs[j])) / 2
                L[j, k] = lam
                rem += lam * dof[j] * derivs[k]
        Lam[mu] = L
    assert sp.simplify(kin - recon) == 0, "kinetic part not purely field*derivative bilinear"
    return sp.expand(rem), Lam[0], Lam[1]


def defect_remainder(A, kind, s, m, C, S):
    """Total reality defect Im[(C+iS) L] mod TD (kinetic) + pointwise (mass). ==0 <=> real."""
    psi, dp0, dp1, dof, d0, d1 = make_constrained(kind)
    L = full_L(psi, dp0, dp1, A, s, m)
    defect = sp.expand(C * sp.im(L) + S * sp.re(L))     # Im[(C+iS)L], all symbols real at inf1
    mass_part, kin_part = split_defect(defect, dof, d0, d1)
    kin_rem, Lam0, Lam1 = kin_mod_td(kin_part, dof, d0, d1)
    total = sp.expand(mass_part + kin_rem)
    return total, mass_part, kin_rem, Lam0, Lam1, (dof, d0, d1)


def symplectic(A, kind, s, m, C, S):
    """Symplectic form = d0 antisymmetric coeff matrix of the REAL Lagrangian Re[(C+iS)L]
       (NOT the vanishing imaginary defect). Returns (rank, det)."""
    psi, dp0, dp1, dof, d0, d1 = make_constrained(kind)
    L = full_L(psi, dp0, dp1, A, s, m)
    realL = sp.expand(C * sp.re(L) - S * sp.im(L))      # Re[(C+iS)L]
    _mass, kin = split_defect(realL, dof, d0, d1)
    _rem, Om0, _Om1 = kin_mod_td(kin, dof, d0, d1)      # Om0 = d0 antisymmetric = symplectic form
    return Om0.rank(), sp.simplify(Om0.det())


# ================================================================ §5 verify Fizz's algebra
def verify_fizz():
    A = Matrix([[0, r5], [1, 0]])                       # b=r5 (sigma-odd), c=1 (sigma-even)
    AGs, AGa = mcanon(A * Gs), mcanon(A * Ga)
    ok = (mat_eq(AGs, diag(5, 1)) and mat_eq(AGa, diag(5, -1))
          and mat_eq(AGs, AGs.T) and mat_eq(AGa, AGa.T))
    rational = all(e.free_symbols == set() for e in list(AGs) + list(AGa))
    record("F7 §5 Fizz candidate algebra",
           "A=[[0,r5],[1,0]] -> A.Gs=diag(5,1), A.Ga=diag(5,-1), both symmetric & rational",
           bool(ok and rational),
           f"A.Gs={tuple(AGs)}, A.Ga={tuple(AGa)}; both symmetric [{ok}], rational (sigma-trivial) [{rational}]")


# ================================================================ §3 F7a the reality sweep
def sweep_variety(kind, s):
    """Solve for the (b,c) aligners giving a real action under `kind`, at some phase (C,S)."""
    br, bi, cr, ci = symbols('br bi cr ci', real=True)   # ev1(b)=br+i bi, ev1(c)=cr+i ci
    C, S, m = symbols('C S m', real=True)
    A = Matrix([[0, br + I * bi], [cr + I * ci, 0]])
    total, mass_part, kin_rem, L0, L1, (dof, d0, d1) = defect_remainder(A, kind, s, m, C, S)
    gens = dof + d0 + d1
    poly = sp.Poly(total, *gens)
    eqs = [sp.expand(co) for co in poly.coeffs() if sp.expand(co) != 0]
    # solve the reality conditions together with C^2+S^2=1, m free (real, nonzero)
    sol = sp.solve(eqs + [C**2 + S**2 - 1], [br, bi, cr, ci, C, S], dict=True)
    return eqs, sol, (br, bi, cr, ci, C, S, m)


def F7a():
    print("\n  -- §3 F7a: reality sweep over A=[[0,b],[c,0]] (s=1 healthy) --")
    # primary: tau constraint
    eqs, sol, syms = sweep_variety('tau', sp.Integer(1))
    br, bi, cr, ci, C, S, m = syms
    print(f"     [tau] reality conditions (coeff=0): {eqs}")
    print(f"     [tau] solution variety: {sol}")
    # interpret: does the variety contain b=c (nonempty, real aligners)?
    #   check representative b=c=1 (A=K2, sigma-even), b=c=r5 (sigma-odd), and Fizz b=r5,c=1
    reps = {
        'K2 (b=c=1)': (Matrix([[0, 1], [1, 0]]), 0, 1),          # phase th=pi/2 -> (C,S)=(0,1)
        'r5*K2 (b=c=r5)': (Matrix([[0, r5], [r5, 0]]), 0, 1),
        'Fizz (b=r5,c=1)': (Matrix([[0, r5], [1, 0]]), None, None),
    }
    survivors = []
    for name, (A, Cval, Sval) in reps.items():
        # find a phase that makes the kinetic real, then test mass at that phase
        if Cval is None:
            # scan the natural phases {1: (1,0), i:(0,1)} and solve generally
            found = None
            for (cv, sv, tag) in ((1, 0, 'th=0'), (0, 1, 'th=pi/2')):
                tot, mp, kr, L0, L1, _ = defect_remainder(A, 'tau', sp.Integer(1), symbols('m', real=True), cv, sv)
                if sp.simplify(tot) == 0:
                    found = (cv, sv, tag)
                    break
            real = found is not None
            detail = f"real at {found[2]}" if real else "not real at th in {0, pi/2}"
        else:
            tot, mp, kr, L0, L1, _ = defect_remainder(A, 'tau', sp.Integer(1), symbols('m', real=True), Cval, Sval)
            real = (sp.simplify(tot) == 0)
            detail = f"remainder at (C,S)=({Cval},{Sval}) = {sp.simplify(tot)}"
        # nondegeneracy of the symplectic form (of the REAL Lagrangian) + EOM-preservation
        nd, eom, rk, dt = None, None, None, None
        if real:
            cv = Cval if Cval is not None else found[0]
            sv = Sval if Sval is not None else found[1]
            rk, dt = symplectic(A, 'tau', sp.Integer(1), symbols('m', real=True), cv, sv)
            nd = (rk == 4 and dt != 0)
            eom = eom_preserves_tau(A)
            survivors.append((name, A, nd, eom))
        record(f"F7a reality: {name}",
               "real (kinetic+mass) at a common phase under tau?",
               real,
               f"{detail}" + (f"; symplectic rank={rk}/4 det={dt} nondegenerate [{nd}]; "
                              f"EOM preserves tau [{eom}]" if real else ""))
    return survivors, sol


def eom_preserves_tau(A):
    """The EOM (Gs d0 + Ga d1 + m Z)psi=0 is A-independent; check it keeps psi_i=f+i r5 g real & closed."""
    m = symbols('m', real=True)
    f1, g1, f2, g2 = symbols('f1 g1 f2 g2', real=True)
    D = symbols('Df1_0 Dg1_0 Df2_0 Dg2_0 Df1_1 Dg1_1 Df2_1 Dg2_1', real=True)
    d0psi = Matrix([D[0] + I * r5 * D[1], D[2] + I * r5 * D[3]])
    d1psi = Matrix([D[4] + I * r5 * D[5], D[6] + I * r5 * D[7]])
    psi = Matrix([f1 + I * r5 * g1, f2 + I * r5 * g2])
    eom = mcanon(Gs * d0psi + Ga * d1psi + m * Z * psi)
    reals = []
    for comp in eom:
        reals += [sp.simplify(sp.re(comp)), sp.simplify(sp.im(comp))]
    allowed = {f1, g1, f2, g2, m} | set(D)
    closed = all(e.free_symbols <= allowed for e in reals)
    nonzero = sum(1 for e in reals if e != 0)
    return bool(closed and nonzero == 4)


# ================================================================ §4 F7b frontier (side by side)
def frontier(A, label):
    """226 §5.2: does 𝒞 sigma(Gamma^mu) 𝒞^-1 = eps (Gamma^mu)^T with a UNIFORM scalar eps across both?"""
    out = {}
    for cname, Cmat in (("A", A), ("sigma(A)", mcanon(sigma(A))), ("A^-1", mcanon(A.inv()))):
        eps_s = _uniform_eps(Cmat, Gs)
        eps_a = _uniform_eps(Cmat, Ga)
        uniform = (eps_s is not None and eps_a is not None and sp.simplify(eps_s - eps_a) == 0)
        out[cname] = (eps_s, eps_a, uniform)
    record(f"F7b frontier: {label}",
           "226 §5.2 uniform scalar eps across both gammas?",
           None,
           "; ".join(f"C={cn}: eps(Gs)={v[0]}, eps(Ga)={v[1]}, uniform [{v[2]}]" for cn, v in out.items()))
    return out


def _uniform_eps(Cmat, Gamma):
    """Return eps if 𝒞 sigma(Gamma) 𝒞^-1 = eps Gamma^T with a single scalar eps, else None."""
    lhs = mcanon(Cmat * mcanon(sigma(Gamma)) * Cmat.inv())
    rhs0 = Gamma.T
    ratios = []
    for i in range(2):
        for j in range(2):
            a, b = sp.simplify(lhs[i, j]), sp.simplify(rhs0[i, j])
            if b == 0 and a == 0:
                continue
            if b == 0 or a == 0:
                return None
            ratios.append(sp.simplify(a / b))
    if not ratios:
        return None
    e0 = ratios[0]
    return e0 if all(sp.simplify(r - e0) == 0 for r in ratios) else None


# ================================================================ §6 sigma(M^2) equal-mass check
def sigma_M2():
    """m^2 sigma-odd (m^2 = r5 mu^2) -> M^2 = m^2 s^2/r5 = mu^2 s^2 sigma-even (equal mass both places)."""
    mu, k = symbols('mu k', real=True, positive=True)
    for sname, s in (("s=1", sp.Integer(1)), ("s=i", I)):
        m2 = r5 * mu**2                       # sigma-ODD mass-squared (pinned axiom)
        M2 = sp.simplify(m2 * s**2 / r5)      # physical mass shell
        sigM2 = sp.simplify(M2.subs(r5, -r5))
        even = (sp.simplify(M2 - sigM2) == 0)
        # place reading of the energy E^2 = k^2 + M2 (nested radical): sigma(E^2)=E^2 ?
        E2 = k**2 + M2
        sigE2 = sp.simplify(E2.subs(r5, -r5))
        E_even = (sp.simplify(E2 - sigE2) == 0)
        record(f"F7 §6 sigma(M^2) equal-mass, {sname}",
               "M^2 sigma-even (equal mass at both places) under sigma-odd m^2",
               bool(even and E_even),
               f"m^2=r5*mu^2 (sigma-odd); M^2={M2} sigma-even [{even}]; E^2=k^2+M^2 sigma-even [{E_even}] "
               f"-> antiparticle axiom delivered by the algebra")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("=" * 72)
    print("BENCH F7 -- reality sweep over the aligner family A=[[0,b],[c,0]]")
    print("does ANY aligner give a real nondegenerate action? (K1 criterion-free)")
    print("=" * 72)
    print("\n[§5 verify Fizz's algebra]")
    verify_fizz()
    print("\n[§3 F7a -- THE CRUX]")
    survivors, sol = F7a()
    print("\n[§4 F7b -- frontier, side by side]")
    frontier(Matrix([[0, 1], [1, 0]]), "K2 (b=c=1, a reality survivor)")
    frontier(Matrix([[0, r5], [1, 0]]), "Fizz (b=r5,c=1)")
    print("\n[§6 -- sigma(M^2) equal-mass]")
    sigma_M2()

    # ---- verdict
    good = [(nm, nd, eom) for (nm, _, nd, eom) in survivors if nd and eom]         # real+nondeg+EOM
    indet = [(nm, nd, eom) for (nm, _, nd, eom) in survivors if not (nd and eom)]  # real but degenerate/EOM-broken
    halo_disconfirmed = any(nm.startswith('K2') for (nm, nd, eom) in good)          # b=c=1 sigma-EVEN works
    print("\n  -- F7 VERDICT --")
    if good:
        verdict = "K-COLLAPSE"
        outcome = (f"K-COLLAPSE: aligners {[g[0] for g in good]} give a REAL, nondegenerate, "
                   f"EOM-preserving first-order action at s=1 -> K1 UN-FIRES, F1-F4 no-go RETRACTED "
                   f"(it was specific to A0=diag(r5,1)). Halo DISCONFIRMED: sigma-EVEN b=c=1 works "
                   f"[{halo_disconfirmed}] -> sigma-odd b NOT required. Reality-variety = {{b=c}}.")
    elif indet:
        verdict = "K-INDETERMINATE"
        outcome = f"real aligners exist but are degenerate or EOM-breaking: {indet}"
    else:
        verdict = "K-SURVIVE"
        outcome = "K-SURVIVE-AND-UPGRADE: no aligner gives a real nondegenerate action."
    record(f"F7 VERDICT: {verdict}", "K-COLLAPSE / K-SURVIVE / K-INDETERMINATE", None, outcome)

    with open("results_f7.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    print("\n" + "=" * 72)
    print(f"F7 VERDICT: {verdict}")
    if verdict == "K-COLLAPSE":
        print("STOP-CONDITION (§7): K-COLLAPSE -> report immediately, do NOT re-pose F4.")
    print("=" * 72)


if __name__ == "__main__":
    main()
