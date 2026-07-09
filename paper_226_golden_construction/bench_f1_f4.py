"""
bench_f1_f4.py — Free Golden Dirac Dynamics, benches F1-F4 (Paper 226).

Verification of Fizz's bricks 1-3 and Tor's brick-4 gate in the 1+1 toy.
Source of record: mr_code_brief_free_golden_dynamics_F1-F4.md (Tor, 8 Jul 2026).
Pre-registration: PRE-REGISTRATION_f1_f4.md (committed before this file emits a verdict).

House style: exact symbolic arithmetic over K(i) = Q(sqrt5, i) via sympy sqrt(5), I.
No floating point until the final F4 signature read (gated). Reality criterion (brief
sec.3-4): Im[ev1 L] modulo total derivatives = 0 (mass pointwise; kinetic mod TD).

Objects are transcribed from brief sec.0; arithmetic + Galois maps reused from the
self-tested golden_algebra.py:
    Gamma_s = G_seed, Gamma_a = G_adj, Z = Z2, sigma (sqrt5->-sqrt5), tau (i->-i).
"""
import os
import sys
import csv
import sympy as sp
from sympy import sqrt, I, symbols, Matrix, diag, eye, simplify, expand
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ---------------------------------------------------------------- brief sec.0 objects
Gs = G_seed                       # Gamma_s = [[0,1],[r5,0]],  Gs^2 = +r5 I
Ga = G_adj                        # Gamma_a = [[0,-1],[r5,0]], Ga^2 = -r5 I
Z = Z2                            # Z = diag(1,-1)
A = diag(r5, 1)                   # pairing matrix A = diag(sqrt5, 1)
Ap = diag(-r5, 1)                 # A' = diag(-sqrt5, 1)
eps = Matrix([[0, 1], [-1, 0]])   # antisymmetriser


def adjoint(psi, Amat=A):
    """Dirac adjoint (Fizz): psibar = (sigma psi)^T . A   (a 1x2 row)."""
    return sigma(psi).T * Amat


ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} - {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


# ================================================================ F1
def F1a():
    """Anticommutant of {Gamma_s, Gamma_a}. EXPECT: span(Z), 1-dimensional."""
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    X = Matrix([[x1, x2], [x3, x4]])
    eqs = list(mcanon(X * Gs + Gs * X)) + list(mcanon(X * Ga + Ga * X))
    sol = sp.linsolve(eqs, [x1, x2, x3, x4])
    (s1, s2, s3, s4), = sol
    free = (s1.free_symbols | s2.free_symbols | s3.free_symbols | s4.free_symbols) & {x1, x2, x3, x4}
    dim = len(free)
    Xsol = Matrix([[s1, s2], [s3, s4]])
    # is the solution space exactly span(Z)?  set the single free param -> 1, compare to Z
    is_Z = False
    if dim == 1:
        p = list(free)[0]
        basis = mcanon(Xsol.subs(p, 1))
        is_Z = mat_eq(basis, Z) or mat_eq(basis, -Z)
    # also record the Gamma_s-only anticommutant shape (should be [[a,b],[-r5 b,-a]])
    eqs_s = list(mcanon(X * Gs + Gs * X))
    sol_s = list(sp.linsolve(eqs_s, [x1, x2, x3, x4]))[0]
    passed = (dim == 1) and is_Z
    record("F1a joint anticommutant",
           "span(Z), one-dimensional",
           passed,
           f"dim(joint anticommutant)={dim}, basis = +-Z [{is_Z}]; "
           f"Gamma_s-only shape (x1,x2,x3,x4)={tuple(sol_s)} = [[a,b],[-r5 b,-a]]")


def F1b():
    """Dispersion. M = -iE Gs + ik Ga + m s Z, det M = 0. EXPECT: E^2=k^2+m^2 s^2/r5."""
    E, k, m = symbols('E k m', real=True)
    ok_all = True
    detail = []
    for sname, s in (("s=1", sp.Integer(1)), ("s=i", I)):
        M = mcanon(-I * E * Gs + I * k * Ga + m * s * Z)
        target_M = Matrix([[m * s, -I * (E + k)], [I * r5 * (k - E), -m * s]])
        M_ok = mat_eq(M, target_M)
        detM = csimp(M.det())
        target_det = csimp(-m**2 * s**2 - r5 * (k**2 - E**2))
        det_ok = (csimp(detM - target_det) == 0)
        # solve det=0 for E^2 -> compare to k^2 + m^2 s^2 / r5
        Esol = sp.solve(sp.Eq(detM, 0), E**2)
        # detM linear in E^2 -> one root
        Msq = csimp(m**2 * s**2 / r5)               # claimed M^2
        E2_target = csimp(k**2 + Msq)
        e2_ok = any(csimp(r - E2_target) == 0 for r in Esol) if Esol else False
        # sign of M^2 at this s (m real, nonzero)
        Msq_val = csimp(Msq.subs(m, 1))
        sign = ">0 healthy" if sp.ask(sp.Q.positive(Msq_val)) else ("<0 tachyonic" if sp.ask(sp.Q.negative(Msq_val)) else "?")
        this_ok = M_ok and det_ok and e2_ok
        ok_all = ok_all and this_ok
        detail.append(f"{sname}: M matches [{M_ok}], detM=-m^2 s^2 -r5(k^2-E^2) [{det_ok}], "
                      f"E^2=k^2+m^2 s^2/r5 [{e2_ok}], M^2={Msq_val} ({sign})")
    record("F1b dispersion / M^2 sign",
           "detM=-m^2 s^2 -r5(k^2-E^2); s=1 healthy, s=i tachyonic",
           ok_all, "; ".join(detail))


# ================================================================ F2
def F2a():
    """Transpose intertwining. EXPECT: Gs^T=+A Gs A^-1, Ga^T=-A Ga A^-1; A' flips; never (+,+)."""
    def parity(Gamma, Amat):
        conj = mcanon(Amat * Gamma * Amat.inv())
        if mat_eq(Gamma.T, conj):
            return +1
        if mat_eq(Gamma.T, -conj):
            return -1
        return 0
    ps_A, pa_A = parity(Gs, A), parity(Ga, A)
    ps_Ap, pa_Ap = parity(Gs, Ap), parity(Ga, Ap)
    ok = (ps_A == +1 and pa_A == -1 and ps_Ap == -1 and pa_Ap == +1
          and not (ps_A == 1 and pa_A == 1))
    record("F2a transpose intertwining",
           "A: Gs parity +1, Ga parity -1; A' exchanges; never (+,+)",
           ok,
           f"A=diag(r5,1): (Gs,Ga) parity=({ps_A:+d},{pa_A:+d}); "
           f"A'=diag(-r5,1): ({ps_Ap:+d},{pa_Ap:+d}); diag(+,+) never occurs [{not (ps_A==1 and pa_A==1)}]")


def F2b():
    """No neutral antisymmetriser. EXPECT: N Gs ~ eps forces N=lam diag(1,-1/r5); then N Ga symmetric, not ~eps."""
    n1, n2, n3, n4, lam = symbols('n1 n2 n3 n4 lam')
    N = Matrix([[n1, n2], [n3, n4]])
    # N Gs = lam eps
    eqs = list(mcanon(N * Gs - lam * eps))
    sol = sp.linsolve(eqs, [n1, n2, n3, n4, lam])
    (a1, a2, a3, a4, al), = sol
    Nsol = Matrix([[a1, a2], [a3, a4]])
    # expected N = lam diag(1,-1/r5) (lam=free param); normalise: set free param -> lam symbol
    free = (a1.free_symbols | a2.free_symbols | a3.free_symbols | a4.free_symbols | al.free_symbols)
    param = list(free)[0] if free else sp.Integer(1)
    Nnorm = mcanon(Nsol.subs(param, 1))
    forced_ok = mat_eq(Nnorm, diag(1, -1 / r5)) or mat_eq(Nnorm, -diag(1, -1 / r5))
    # now N Ga : is it ~ eps ?  (i.e. of the form mu*eps = [[0,mu],[-mu,0]])
    NGa = mcanon(Nsol.subs(param, 1) * Ga)
    symmetric = mat_eq(NGa, NGa.T)                    # symmetric => NOT antisymmetric
    prop_eps = mat_eq(NGa, (NGa[0, 1]) * eps)         # exactly proportional to eps?
    passed = forced_ok and symmetric and (not prop_eps)
    record("F2b no neutral antisymmetriser",
           "N Gs~eps forces N=lam diag(1,-1/r5); N Ga symmetric, not ~eps",
           passed,
           f"N (up to scale) = {tuple(Nnorm)} = diag(1,-1/r5) [{forced_ok}]; "
           f"N Ga = {tuple(NGa)} symmetric [{symmetric}], proportional to eps [{prop_eps}] "
           f"-> no constant N antisymmetrises both")


# ================================================================ F3 field builders
# generic (unconstrained) sigma-field: psi_i = (a_i + b_i i) + (c_i + d_i i) r5
def field_generic():
    a1, b1, c1, d1, a2, b2, c2, d2 = symbols('a1 b1 c1 d1 a2 b2 c2 d2', real=True)
    psi = Matrix([(a1 + b1 * I) + (c1 + d1 * I) * r5,
                  (a2 + b2 * I) + (c2 + d2 * I) * r5])
    return psi, (a1, b1, c1, d1, a2, b2, c2, d2)


# A constrained field for the kinetic gate carries independent derivative symbols.
def make_constrained(kind):
    """Return (psi, dpsi0, dpsi1, dof, d0, d1) for a 4-real-dof reality constraint.
       kind='tau'   : psi_i = f + I r5 g   (p real, q imaginary)  [F3b constraint]
       kind='sigma' : psi_i = f + I g      (q=0, psi in Q(i))     [sigma-real]
       kind='C0'    : psi_1 = f+I g in Q(i); psi_2 = (f+I g) r5   [Z sigma Majorana]
    """
    f1, g1, f2, g2 = symbols('f1 g1 f2 g2', real=True)
    f1t, g1t, f2t, g2t = symbols('f1t g1t f2t g2t', real=True)
    f1x, g1x, f2x, g2x = symbols('f1x g1x f2x g2x', real=True)
    dof = [f1, g1, f2, g2]
    d0 = [f1t, g1t, f2t, g2t]
    d1 = [f1x, g1x, f2x, g2x]

    def comp(f, g, kind, which):
        if kind == 'tau':
            return f + I * r5 * g
        if kind == 'sigma':
            return f + I * g
        if kind == 'C0':
            return (f + I * g) if which == 1 else (f + I * g) * r5
        raise ValueError(kind)

    psi = Matrix([comp(f1, g1, kind, 1), comp(f2, g2, kind, 2)])
    dpsi0 = Matrix([comp(f1t, g1t, kind, 1), comp(f2t, g2t, kind, 2)])
    dpsi1 = Matrix([comp(f1x, g1x, kind, 1), comp(f2x, g2x, kind, 2)])
    return psi, dpsi0, dpsi1, dof, d0, d1


def kinetic_density(psi, dpsi0, dpsi1, factor=1):
    """factor * psibar (Gs d0 + Ga d1) psi   (scalar over K(i))."""
    Dpsi = Gs * dpsi0 + Ga * dpsi1
    return sp.expand(factor * (adjoint(psi) * Dpsi)[0])


def imag_mod_td(expr, dof, d0, d1):
    """Im(expr) reduced modulo total derivatives.
       expr is bilinear in (field f in dof)*(derivative in d0 or d1).
       Returns (remainder, Lambda0, Lambda1): the antisymmetric leftover and the
       antisymmetric coefficient matrices per direction. remainder==0 <=> real mod TD."""
    im = sp.expand(sp.im(expr))
    n = len(dof)
    remainder = sp.Integer(0)
    Lam = {}
    recon = sp.Integer(0)
    for mu, derivs in ((0, d0), (1, d1)):
        C = [[sp.diff(im, dof[j], derivs[k]) for k in range(n)] for j in range(n)]
        L = sp.zeros(n, n)
        for j in range(n):
            for k in range(n):
                recon += C[j][k] * dof[j] * derivs[k]
                lam = sp.simplify((C[j][k] - C[k][j]) / 2)
                L[j, k] = lam
                remainder += lam * dof[j] * derivs[k]
        Lam[mu] = L
    # safety: Im(expr) really is field*derivative bilinear (no field^2 / deriv^2 leakage)
    assert sp.simplify(im - recon) == 0, "kinetic density not purely field*derivative bilinear"
    return sp.simplify(sp.expand(remainder)), Lam[0], Lam[1]


def mass_imag(psi, s):
    """Im[ ev1( m s B ) ], B = psibar Z psi.  Non-derivative => pointwise."""
    m = symbols('m', real=True)
    B = sp.expand((adjoint(psi) * Z * psi)[0])
    return sp.simplify(sp.im(sp.expand(m * s * B))), sp.simplify(B)


def F3a():
    """Unconstrained sigma-field: neither s makes Im[ev1 m s B]=0 generically."""
    psi, syms = field_generic()
    res = {}
    for sname, s in (("s=1", sp.Integer(1)), ("s=i", I)):
        im, B = mass_imag(psi, s)
        res[sname] = (im, B)
    # verify B = r5 N(psi1) - N(psi2), N=p^2-5q^2
    a1, b1, c1, d1, a2, b2, c2, d2 = syms
    p1, q1 = a1 + b1 * I, c1 + d1 * I
    p2, q2 = a2 + b2 * I, c2 + d2 * I
    Bformula = sp.expand(r5 * (p1**2 - 5 * q1**2) - (p2**2 - 5 * q2**2))
    B_ok = (sp.simplify(res["s=1"][1] - Bformula) == 0)
    nonzero1 = (res["s=1"][0] != 0)
    nonzeroi = (res["s=i"][0] != 0)
    passed = B_ok and nonzero1 and nonzeroi
    record("F3a unconstrained mass reality",
           "neither s=1 nor s=i gives Im=0 generically (forces a constraint)",
           passed,
           f"B = r5 N(psi1)-N(psi2) [{B_ok}]; Im[ev1 m*1*B] generically nonzero [{nonzero1}]; "
           f"Im[ev1 m*i*B] generically nonzero [{nonzeroi}] -> field-reality constraint mandatory")


def F3b_and_gate():
    """F3b mass reality under tau; sub-checks 1-3; reconciliation; verdict."""
    m = symbols('m', real=True)

    # ---- F3b: mass real & healthy at s=1 under tau
    psi_t, d0_t, d1_t, dof_t, D0_t, D1_t = make_constrained('tau')
    im_mass_t_s1, B_t = mass_imag(psi_t, sp.Integer(1))
    im_mass_t_si, _ = mass_imag(psi_t, I)
    f1, g1, f2, g2 = dof_t
    N1 = sp.simplify((sigma(Matrix([psi_t[0]])) * Matrix([psi_t[0]]))[0])  # (sigma psi1) psi1 = N(psi1)
    N2 = sp.simplify((sigma(Matrix([psi_t[1]])) * Matrix([psi_t[1]]))[0])
    N1_ok = (sp.simplify(N1 - (f1**2 + 5 * g1**2)) == 0)
    N2_ok = (sp.simplify(N2 - (f2**2 + 5 * g2**2)) == 0)
    mass_real_tau_s1 = (im_mass_t_s1 == 0)
    f3b_ok = mass_real_tau_s1 and N1_ok and N2_ok
    record("F3b tau-constrained mass reality",
           "tau -> p real,q imaginary; N=a^2+5d^2>0; s=1 mass real",
           f3b_ok,
           f"Im[ev1 m*1*B] under tau = {im_mass_t_s1} [real={mass_real_tau_s1}]; "
           f"N(psi1)=f^2+5g^2 [{N1_ok}], N(psi2)=f^2+5g^2 [{N2_ok}] (real, positive)")

    # ---- reconciliation: is the tau-constraint exactly "sigma = tau on the physical locus"?
    sig_eq_tau = mat_eq(mcanon(sigma(psi_t)), mcanon(tau(psi_t)))
    record("F3 reconciliation sigma=tau on locus",
           "F3b constraint is exactly sigma=tau on the physical (constrained) locus",
           sig_eq_tau,
           f"sigma(psi)=tau(psi) on the tau-constrained field [{sig_eq_tau}] "
           f"(=> sigma-Majorana and tau-Majorana coincide on-shell; Paper 226 sec.5.3)")

    # ---- sub-check 1 (THE GATE): kinetic reality, swept over constraints x factor
    print("\n  -- F3 sub-check 1: kinetic Im mod TD (the gate) --")
    grid = {}
    for kind in ('tau', 'sigma', 'C0'):
        psi, dp0, dp1, dof, D0, D1 = make_constrained(kind)
        for fname, factor in (("factor=1", sp.Integer(1)), ("factor=i", I)):
            dens = kinetic_density(psi, dp0, dp1, factor)
            rem, L0, L1 = imag_mod_td(dens, dof, D0, D1)
            grid[(kind, fname)] = (rem, L0, L1)
            print(f"     [{kind:5s} | {fname:9s}] Im mod TD = {rem}  -> real={rem == 0}")
    # primary object: brief's literal L is factor=1 under tau
    rem_tau1, L0_tau1, L1_tau1 = grid[('tau', 'factor=1')]
    kinetic_real_tau = (rem_tau1 == 0)
    # which (kind,factor) give a real kinetic term at all?
    real_combos = [key for key, (rem, _, _) in grid.items() if rem == 0]

    # ---- sub-check 2: EOM preserves the tau-constraint (s=1)
    # substitute psi_i = f_i + I r5 g_i into (Gs d0 + Ga d1 + m Z)psi = 0 and check the
    # resulting real system stays real & closed (no term forces f,g complex).
    f1s, g1s, f2s, g2s = symbols('f1 g1 f2 g2', real=True)
    # abstract derivatives as real symbols
    Df1_0, Dg1_0, Df2_0, Dg2_0 = symbols('Df1_0 Dg1_0 Df2_0 Dg2_0', real=True)
    Df1_1, Dg1_1, Df2_1, Dg2_1 = symbols('Df1_1 Dg1_1 Df2_1 Dg2_1', real=True)
    d0psi = Matrix([Df1_0 + I * r5 * Dg1_0, Df2_0 + I * r5 * Dg2_0])
    d1psi = Matrix([Df1_1 + I * r5 * Dg1_1, Df2_1 + I * r5 * Dg2_1])
    psi_c = Matrix([f1s + I * r5 * g1s, f2s + I * r5 * g2s])
    eom = mcanon(Gs * d0psi + Ga * d1psi + m * Z * psi_c)
    # split each of the 2 complex EOM comps into (real, imag) at ev1 -> 4 real equations
    real_eqs = []
    for comp in eom:
        real_eqs.append(sp.simplify(sp.re(comp)))
        real_eqs.append(sp.simplify(sp.im(comp)))
    # closure: every real equation involves ONLY the real dof/derivative symbols (already real=True)
    allowed = {f1s, g1s, f2s, g2s, Df1_0, Dg1_0, Df2_0, Dg2_0,
               Df1_1, Dg1_1, Df2_1, Dg2_1, m}
    closed = all(eq.free_symbols <= allowed for eq in real_eqs)
    # non-triviality: the 4 real equations are independent and nonzero
    nonzero_eqs = sum(1 for eq in real_eqs if eq != 0)
    eom_preserves = closed and (nonzero_eqs == 4)
    record("F3 sub-check 2 EOM preserves tau",
           "golden Dirac evolution keeps the tau-constraint (real, closed system)",
           eom_preserves,
           f"psi_i=f_i+i r5 g_i in EOM -> {nonzero_eqs} real closed equations for (f_i,g_i); "
           f"no term forces complex f,g [{closed}] (a,d sectors evolve by the same golden-Dirac system)")

    # ---- sub-check 3: nondegeneracy of the kinetic/symplectic form on the tau-halved space
    # symplectic form = antisymmetric d0-coefficient matrix of the (real) kinetic Lagrangian.
    # Use the imaginary current L0 from the brief's literal L (factor=1, tau): the d0 antisym part.
    Omega = L0_tau1
    detOmega = sp.simplify(Omega.det())
    rankOmega = Omega.rank()
    nondeg = (detOmega != 0) and (rankOmega == 4)
    record("F3 sub-check 3 nondegeneracy",
           "kinetic/symplectic form nondegenerate on the tau-halved config space",
           nondeg,
           f"d0 antisymmetric (symplectic) matrix rank={rankOmega}/4, det={detOmega} "
           f"-> nondegenerate [{nondeg}] (no K1 from degeneracy)")

    # ---- VERDICT
    # CLEAR requires: exists a single constraint making kinetic (brief's literal factor=1)
    # AND mass real & healthy at s=1, nondegenerate.
    mass_real = {'tau': (im_mass_t_s1 == 0)}
    for kind in ('sigma', 'C0'):
        psi_k, *_ = make_constrained(kind)
        im_k, _ = mass_imag(psi_k, sp.Integer(1))
        mass_real[kind] = (im_k == 0)
    kinetic_real_factor1 = {kind: (grid[(kind, 'factor=1')][0] == 0) for kind in ('tau', 'sigma', 'C0')}
    common = [kind for kind in ('tau', 'sigma', 'C0')
              if mass_real[kind] and kinetic_real_factor1[kind]]
    clear = bool(common) and nondeg

    print("\n  -- F3 VERDICT --")
    print(f"     mass real at s=1 by constraint : {mass_real}")
    print(f"     kinetic real (factor=1) by cons: {kinetic_real_factor1}")
    print(f"     kinetic real (any factor)      : {sorted(real_combos)}")
    verdict = "CLEAR" if clear else "COLLISION"
    if clear:
        outcome = (f"a single constraint {common} makes kinetic (literal L) AND mass real at s=1, "
                   f"nondegenerate -> hand to F4")
    else:
        # is the kinetic term real under a DIFFERENT (sigma-type) constraint than the mass wants?
        mass_c = [k for k in mass_real if mass_real[k]]
        kin_c = [k for k in kinetic_real_factor1 if kinetic_real_factor1[k]]
        outcome = (f"no single constraint reales BOTH kinetic(literal L) and mass at s=1: "
                   f"mass real under {mass_c}, kinetic(factor=1) real under {kin_c or 'none'}; "
                   f"kinetic Im mod TD under tau,s=1 = {rem_tau1} (nonzero antisym d0 current) "
                   f"-> K1/K3 fire, split sigma-theory sick at free level, compact-or-bust")
    record(f"F3 JOINT-GATE VERDICT: {verdict}",
           "CLEAR (hand to F4) or COLLISION (K1/K2/K3 fire at free level)",
           True, outcome)
    return clear, {'rem_tau1': rem_tau1, 'grid': grid, 'kinetic_real_factor1': kinetic_real_factor1,
                   'mass_real': mass_real}


# ================================================================ genericity (sec.1)
def genericity():
    """Re-run F2a parity + sub-check-1 kinetic obstruction with sqrt5 -> sqrt13.
       If the same diag(+,-) parity and the same nonzero kinetic remainder appear, the
       no-go is real-quadratic-generic (not golden)."""
    print("\n  -- Genericity: real-quadratic (sqrt5) vs (sqrt13) --")
    results = {}
    for D in (5, 13):
        r = sqrt(D)
        Gs_r = Matrix([[0, 1], [r, 0]])
        Ga_r = Matrix([[0, -1], [r, 0]])
        Zr = Matrix([[1, 0], [0, -1]])
        Ar = diag(r, 1)

        def sig_r(M):
            return M.applyfunc(lambda e: sp.expand(e.subs(r, -r)))

        # F2a parity
        def parity(Gamma):
            conj = sp.expand(Ar * Gamma * Ar.inv())
            if sp.simplify(Gamma.T - conj) == sp.zeros(2, 2):
                return +1
            if sp.simplify(Gamma.T + conj) == sp.zeros(2, 2):
                return -1
            return 0
        par = (parity(Gs_r), parity(Ga_r))
        # kinetic Im mod TD under tau-analog: psi_i = f + I r g
        f1, g1, f2, g2 = symbols('f1 g1 f2 g2', real=True)
        f1t, g1t, f2t, g2t = symbols('f1t g1t f2t g2t', real=True)
        f1x, g1x, f2x, g2x = symbols('f1x g1x f2x g2x', real=True)
        dof = [f1, g1, f2, g2]
        D0 = [f1t, g1t, f2t, g2t]
        D1 = [f1x, g1x, f2x, g2x]
        psi = Matrix([f1 + I * r * g1, f2 + I * r * g2])
        dp0 = Matrix([f1t + I * r * g1t, f2t + I * r * g2t])
        dp1 = Matrix([f1x + I * r * g1x, f2x + I * r * g2x])
        dens = sp.expand(((sig_r(psi).T * Ar) * (Gs_r * dp0 + Ga_r * dp1))[0])
        im = sp.expand(sp.im(dens))
        rem = sp.Integer(0)
        for derivs in (D0, D1):
            for j in range(4):
                for k in range(4):
                    C = sp.diff(im, dof[j], derivs[k])
                    Ckj = sp.diff(im, dof[k], derivs[j])
                    rem += sp.simplify((C - Ckj) / 2) * dof[j] * derivs[k]
        rem = sp.simplify(rem)
        results[D] = (par, rem)
        print(f"     D={D:2d}: F2a parity(Gs,Ga)={par}; kinetic Im mod TD under tau = {rem} (real={rem==0})")
    same_parity = results[5][0] == results[13][0] == (+1, -1)
    both_nonzero = results[5][1] != 0 and results[13][1] != 0
    generic = same_parity and both_nonzero
    record("Genericity (real-quadratic, not golden)",
           "sqrt5 and sqrt13 give identical diag(+,-) parity and nonzero kinetic obstruction",
           generic,
           f"parity diag(+,-) for both D=5,13 [{same_parity}]; kinetic Im mod TD nonzero for both "
           f"[{both_nonzero}] -> the no-go is a real-quadratic no-go, holds for Q(sqrt13) identically")


# ================================================================ F4 (GATED)
def F4(info):
    """Runs only if F3 CLEARS at healthy s=1. Otherwise gated off (brief sec.6)."""
    print("\n" + "=" * 72)
    print("F4 (conservation + K4 signature) -- GATED on F3 = CLEAR at s=1")
    print("=" * 72)
    E, k, m = symbols('E k m', real=True, positive=True)
    s = sp.Integer(1)
    M2 = m**2 * s**2 / r5
    # u(k) = (i(E+k), m s)^T ; E = sqrt(k^2 + M2)
    def u(sign):
        Ee = sign * sqrt(k**2 + M2)
        return Matrix([I * (Ee + k), m * s]), Ee
    # F4b: pairing <u,u> = (sigma u)^T A Gs u at inf1
    def pairing(ua, ub):
        return sp.simplify(((sigma(ua).T * A) * Gs * ub)[0])
    up, Ep = u(+1)
    um, Em = u(-1)
    norm_p = pairing(up, up)
    norm_m = pairing(um, um)
    print(f"   <u+,u+> = {norm_p}")
    print(f"   <u-,u-> = {norm_m}")
    real_norm_p = sp.simplify(-I * norm_p) if sp.im(norm_p) != 0 else norm_p
    real_norm_m = sp.simplify(-I * norm_m) if sp.im(norm_m) != 0 else norm_m
    record("F4b K4 signature",
           "opposite definite signs on u+/u-; sigma maps +freq -> -freq",
           None if True else True,
           f"<u+,u+>={norm_p}, <u-,u->={norm_m}; extracted real norms {real_norm_p},{real_norm_m}")


# ================================================================ main
def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))  # outputs land next to the script
    print("=" * 72)
    print("BENCH F1-F4 -- Free Golden Dirac Dynamics (Paper 226), exact over Q(sqrt5,i)")
    print("=" * 72)
    print("\n[F1]")
    F1a(); F1b()
    print("\n[F2]")
    F2a(); F2b()
    print("\n[F3]")
    F3a()
    clear, info = F3b_and_gate()
    genericity()

    if clear:
        F4(info)
    else:
        print("\n" + "=" * 72)
        print("F4 GATED OFF: F3 = COLLISION. Stop-condition (brief sec.6) -> do not run F4.")
        print("=" * 72)

    with open("results_f1_f4.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    print("\n" + "=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} items match EXPECT (see results_f1_f4.csv)")
    print(f"F3 JOINT-GATE VERDICT: {'CLEAR' if clear else 'COLLISION'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
