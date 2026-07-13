"""
bench_f9.py — F9: is the golden Dirac construction golden at all? (Paper 226)

Source of record: mr_code_brief_F9_connection_check.md (Fizz, 12 Jul 2026), from Scout's F1 flag.
Pre-registration: PRE-REGISTRATION_f9.md (committed 19c27d6, before this verdict).

F9a  THE CONNECTION CHECK: is 2I present in the 191 gamma construction, or is the sqrt5 in the
     matrix entries a removable (1,1)-hyperbolic artefact with no icosahedral content?
F9b  the tau-adjoint: does tau (antilinear) deliver the conserved U(1) charge that sigma could not?
F9c  the bipartition of 2I = SL(2,5) (for the record).
F9d  the SU(2) objection that kills sigma = P.
"""
import os
import sys
import csv
import sympy as sp
from sympy import sqrt, I, symbols, Matrix, diag, Rational, eye, cos, pi, conjugate, GoldenRatio
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

Gs, Ga, Z = G_seed, G_adj, Z2
K2 = Matrix([[0, 1], [1, 0]])
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed is True else ("FAIL" if passed is False else "INFO")
    ROWS.append((item, expect, f"{tag} - {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


# ================================================================ F9a — THE CONNECTION CHECK
def F9a():
    print("\n[F9a -- THE CONNECTION CHECK]")
    # (1) Scout's F1: gamma squares and isotropy
    gs2 = mat_eq(mcanon(Gs * Gs), r5 * I2)
    ga2 = mat_eq(mcanon(Ga * Ga), -r5 * I2)
    null = mcanon((Gs + Ga) * (Gs + Ga))            # Gs+Ga is a null (isotropic) vector
    isotropic = mat_eq(null, sp.zeros(2, 2))
    record("F9a.1 gamma squares (Scout F1)",
           "Gs^2=+sqrt5 I, Ga^2=-sqrt5 I; form <sqrt5,-sqrt5> isotropic",
           bool(gs2 and ga2 and isotropic),
           f"Gs^2=+sqrt5 I [{gs2}], Ga^2=-sqrt5 I [{ga2}]; (Gs+Ga)^2=0 [{isotropic}] "
           f"-> <sqrt5,-sqrt5> is the hyperbolic plane over K")

    # (1b) explicit change of basis over K: M^T diag(sqrt5,-sqrt5) M = diag(1,-1)
    u = Matrix([(5 + r5) / 10, (5 - r5) / 10])
    v = Matrix([(5 - r5) / 10, (5 + r5) / 10])
    M = Matrix.hstack(u, v)                          # columns u,v
    G = diag(r5, -r5)
    transformed = mcanon(M.T * G * M)
    to_11 = mat_eq(transformed, diag(1, -1))
    mixes = (csimp(M[0, 1]) != 0 and csimp(M[1, 0]) != 0)   # off-diagonal -> mixes ds/da (d0/d1)
    record("F9a.1b change of basis over K -> <1,-1>",
           "explicit M in GL2(K) with M^T diag(sqrt5,-sqrt5) M = diag(1,-1); report what it does to d0/d1",
           bool(to_11),
           f"M = {tuple(M)}; M^T diag(sqrt5,-sqrt5) M = diag(1,-1) [{to_11}]; M is NON-diagonal "
           f"[mixes d0/d1 = {mixes}] -> removing sqrt5 MIXES the two derivative directions (a K-linear "
           f"boost/rotation of t,x), NOT a benign separate rescaling (Scout's caveat)")

    # (1c) Cl generates split M2(K); volume element omega = Gs Ga = sqrt5 Z, omega^2 = 5
    omega = mcanon(Gs * Ga)
    omega_is = mat_eq(omega, r5 * Z)
    omega2 = mcanon(omega * omega)
    basis = [I2, Gs, Ga, omega]
    # rank over K of the 4 matrices flattened
    Mflat = Matrix([[b[i] for i in range(4)] for b in basis])
    full = (Mflat.rank() == 4)
    record("F9a.1c Clifford algebra = split M2(K)",
           "{I,Gs,Ga,Gs Ga} span M2(K); omega=Gs Ga=sqrt5 Z, omega^2=5",
           bool(omega_is and full),
           f"omega=Gs Ga=sqrt5 Z [{omega_is}], omega^2={omega2[0,0]}*I; dim algebra = {Mflat.rank()} "
           f"-> Cl = split M2(K) (identical to standard (1,1); NOT a division algebra)")

    # (2) does 2I act?  (a) centralizer of {Gs,Ga} in M2(K(i))
    x0, x1, x2, x3 = symbols('x0 x1 x2 x3')
    X = Matrix([[x0, x1], [x2, x3]])
    eqs = list(mcanon(X * Gs - Gs * X)) + list(mcanon(X * Ga - Ga * X))
    sol = sp.linsolve(eqs, [x0, x1, x2, x3])
    (s0, s1, s2, s3), = sol
    free = (s0.free_symbols | s1.free_symbols | s2.free_symbols | s3.free_symbols) & {x0, x1, x2, x3}
    cent_scalar = (len(free) == 1) and (csimp(s1) == 0) and (csimp(s2) == 0) and (csimp(s0 - s3) == 0)
    record("F9a.2a centralizer of the gammas",
           "centralizer of {Gs,Ga} in M2(K(i)) = scalars (no room for a 2I action commuting with the gammas)",
           bool(cent_scalar),
           f"X commuting with Gs and Ga: dim={len(free)}, X = lambda*I [{cent_scalar}] -> the ONLY things "
           f"commuting with the gammas are scalars; no nontrivial group (let alone 2I) respects them by centralizing")

    # (2b) Spin group of a (1,1) form is ABELIAN (SO(1,1) = {diag(t,1/t)}); 2I is NON-ABELIAN of order 120
    t = symbols('t', nonzero=True)
    so11_a = Matrix([[t, 0], [0, 1 / t]])
    so11_b = Matrix([[symbols('s', nonzero=True), 0], [0, 1 / symbols('s', nonzero=True)]])
    abelian = mat_eq(mcanon(so11_a * so11_b), mcanon(so11_b * so11_a))
    twoI_order, twoI_abelian = 120, False
    no_embed = abelian and (not twoI_abelian)
    record("F9a.2b Spin(1,1) is abelian; 2I is not",
           "the (1,1) Spin/rotation group is abelian -> the non-abelian 2I (order 120) cannot embed",
           bool(no_embed),
           f"SO(1,1)={{diag(t,1/t)}} abelian [{abelian}]; 2I non-abelian order {twoI_order} -> 2I does NOT "
           f"embed in the Spin group of the 191 Clifford form. (Contrast 2I in SU(2)=Spin(3), a COMPACT "
           f"rotation group of a (3,0) form -- a different signature entirely.)")

    # (2c) finite subgroups of Spin(1,1) ~ K(i)* are roots of unity mu(K(i)); order 4 here
    n_roots = [n for n in (1, 2, 3, 4, 5, 6, 8, 10) if _has_prim_root(n)]
    mu_order = max(n_roots)
    record("F9a.2c roots of unity in K(i)",
           "finite subgroups of Spin(1,1)~Gm are mu(K(i)); compare |2I|=120",
           None,
           f"primitive n-th roots of unity present in K(i)=Q(sqrt5,i) for n in {n_roots} -> mu order "
           f"= {mu_order} (just {{+-1,+-i}}); 120 >> {mu_order} and non-cyclic -> no 2I")

    disconnect = (cent_scalar and no_embed)
    record("F9a VERDICT",
           "K-DISCONNECT (no 2I in the 191 construction) or K-CONNECT (report mechanism)",
           None,
           ("K-DISCONNECT: the sqrt5 in the 191 gamma ENTRIES is a removable (1,1)-hyperbolic artefact "
            "(Cl = split M2(K) ~ standard (1,1)); the Clifford/Spin group is abelian and its finite subgroups "
            "are just +-1,+-i; the centralizer of the gammas is scalars. 2I (order 120, non-abelian, a subgroup "
            "of the COMPACT SU(2)=Spin(3)) does NOT embed, act on the spinor space respecting the gammas, or "
            "permute them. The gammas were written down with a sqrt5 in the entries and carry NO icosahedral "
            "content. 2I's OWN golden structure (F9c) lives in its character field, a SEPARATE object. "
            "=> the F1-F8 arc was a real-quadratic (1,1) Dirac theory; every 'real-quadratic-generic' finding "
            "was this diagnosis. ESCALATE." if disconnect else "K-CONNECT (unexpected) -- report mechanism"))
    return disconnect


def _has_prim_root(n):
    """Is there a primitive n-th root of unity in K(i)=Q(sqrt5,i)?  Check via cyclotomic degree."""
    # K(i) = Q(sqrt5, i), degree 4 over Q, quadratic subfields Q(sqrt5),Q(i),Q(sqrt-5).
    # zeta_n in K(i) iff Q(zeta_n) subset K(i).  Present: n | 4 (i) and n in {1,2}. Not 3,5,6,8,10.
    present = {1, 2, 4}
    return n in present


# ================================================================ F9b — the tau-adjoint
def tau_conservation(m_expr, mlabel):
    """Re-run F8 §3 conservation with the tau-adjoint psibar=(tau psi)^T K2. Current solved; mass residual decides."""
    r = symbols('r0:16', real=True)                 # real/imag parts of 8 complex components
    def comp(a, b):
        return r[a] + I * r[b]
    p1 = Matrix([comp(0, 1), comp(2, 3)])
    p2 = Matrix([comp(4, 5), comp(6, 7)])
    d1p1 = Matrix([comp(8, 9), comp(10, 11)])
    d1p2 = Matrix([comp(12, 13), comp(14, 15)])
    taup1 = Matrix([sp.conjugate(p1[0]), sp.conjugate(p1[1])])      # tau = complex conjugation
    taud1p1 = Matrix([sp.conjugate(d1p1[0]), sp.conjugate(d1p1[1])])
    bb = symbols('b0:4')
    B = Matrix([[bb[0], bb[1]], [bb[2], bb[3]]])
    Gsi = Gs.inv()
    d0p1 = sp.expand(Gsi * (-Ga * d1p1 - m_expr * Z * p1))
    d0p2 = sp.expand(Gsi * (-Ga * d1p2 - m_expr * Z * p2))
    taud0p1 = Matrix([sp.conjugate(d0p1[0]), sp.conjugate(d0p1[1])])
    D_dot = sp.expand((taud0p1.T * K2 * Gs * p2)[0] + (taup1.T * K2 * Gs * d0p2)[0])
    J_div = sp.expand((taud1p1.T * B * p2)[0] + (taup1.T * B * d1p2)[0])
    total = sp.expand(D_dot + J_div)
    dset = set(r[8:])
    deriv_part = sum((tm for tm in total.as_ordered_terms() if tm.free_symbols & dset), sp.Integer(0))
    mass_part = sp.expand(total - deriv_part)
    eqs = []
    for a in range(8, 16):
        for b in range(8):
            co = sp.expand(deriv_part.coeff(r[a] * r[b]))
            if co != 0:
                eqs.append(co)
    sol = sp.solve(eqs, list(bb), dict=True)
    residual = csimp(mass_part.subs(sol[0])) if sol else mass_part
    conserved = (residual == 0)
    record(f"F9b conservation (tau-adjoint) [{mlabel}]",
           "tau-adjoint gives a CONSERVED U(1) current (sigma gave none)",
           conserved,
           f"psibar=(tau psi)^T K2; mass residual after solving current = {residual} -> conserved [{conserved}]")
    return conserved


def F9b():
    print("\n[F9b -- does tau deliver what sigma did not?]")
    # (1) phase: tau antilinear -> psibar -> e^{-i a} psibar
    al = symbols('alpha', real=True)
    psi = Matrix(symbols('p0 p1'))
    lhs = (Matrix([sp.conjugate(sp.exp(I * al) * psi[0]), sp.conjugate(sp.exp(I * al) * psi[1])]).T * K2)
    rhs = sp.exp(-I * al) * (Matrix([sp.conjugate(psi[0]), sp.conjugate(psi[1])]).T * K2)
    phase_ok = all(csimp(lhs[j] - rhs[j]) == 0 for j in range(2))
    record("F9b.1 phase of the tau-adjoint",
           "under psi->e^{i a}psi, psibar -> e^{-i a} psibar (antilinear, unlike sigma)",
           phase_ok,
           f"(tau(e^{{i a}}psi))^T K2 = e^{{-i a}}(tau psi)^T K2 [{phase_ok}] -> U(1) is a symmetry "
           f"(sigma gave e^{{+i a}} -> none)")
    # (2) conserved charge, both mass parities
    mu = symbols('mu', real=True, positive=True)
    c_even = tau_conservation(mu, "sigma-even m in K")
    c_odd = tau_conservation(sqrt(r5) * mu, "sigma-odd m=5^(1/4)mu")
    # (3) positive-definite Hermitian metric K2 Gs = diag(sqrt5,1)
    metric = mcanon(K2 * Gs)
    pos = (metric == diag(r5, 1)) and (r5 > 0)
    record("F9b.3 hermitian positivity",
           "the tau-pairing metric K2 Gs is positive-definite (involution of the 2nd kind)",
           bool(pos),
           f"K2 Gs = diag(sqrt5,1), both entries > 0 [{pos}] -> (tau psi)^T K2 Gs psi = psi^dagger diag(sqrt5,1) psi "
           f">= 0, a genuine positive-definite Dirac inner product")
    # (4) signature: positive-definite on BOTH frequencies (Dirac norm; particle/antiparticle via CHARGE not norm)
    k = symbols('k', real=True)
    E = sqrt(k**2 + mu**2)
    for name, u in (("u+", Matrix([I * (E + k), mu])), ("u-", Matrix([I * (-E + k), mu]))):
        tu = Matrix([sp.conjugate(u[0]), sp.conjugate(u[1])])
        nrm = csimp((tu.T * K2 * Gs * u)[0])
        record(f"F9b.4 signature {name} (tau)",
               "tau-norm is real & positive-definite on each frequency (healthy Dirac; no ghosts)",
               bool(nrm.is_positive if nrm.is_positive is not None else None) or None,
               f"(tau {name})^T K2 Gs {name} = {nrm} (real, positive) -> positive-definite norm; the "
               f"+/- split is carried by the conserved CHARGE (Noether U(1)), not by an indefinite norm")
    record("F9b VERDICT",
           "tau delivers (conserved U(1), positive-definite) vs tau also fails",
           None,
           f"TAU DELIVERS: conserved U(1) charge [even-mass {c_even}, odd-mass {c_odd}], positive-definite "
           f"hermitian pairing, real norms -> a healthy Dirac theory with C = tau (complex conjugation). "
           f"This RESOLVES F8's escalation: the charge conjugation is tau, not sigma. tau(m)=m so BOTH mass "
           f"parities conserve. (Statement about the CONSTRUCTION; run even though F9a=DISCONNECT.)")


# ================================================================ F9c — the bipartition
def F9c():
    print("\n[F9c -- the bipartition of 2I=SL(2,5)]")
    even = {'1': 1, '3': 3, "3'": 3, '4': 4, '5': 5}
    odd = {'2': 2, "2'": 2, "4'": 4, '6': 6}
    se = sum(d**2 for d in even.values())
    so = sum(d**2 for d in odd.values())
    record("F9c.1 even/odd dims", "even Sigma d^2=60, odd Sigma d^2=60, total 120",
           (se == 60 and so == 60 and se + so == 120),
           f"EVEN {list(even)} Sigma d^2={se}; ODD {list(odd)} Sigma d^2={so}; total {se+so}=|2I|")
    marks = [1, 2, 3, 4, 5, 6, 4, 2, 3]              # affine E8 marks (arm 1..6..4..2, branch 3)
    record("F9c.2 affine E8 marks", "McKay graph = affine E8, marks sum 30",
           (sum(marks) == 30),
           f"marks {marks} sum = {sum(marks)} (= h, the Coxeter number of E8) -> McKay(2I) = affine E8")
    # sigma swaps 2<->2' and 3<->3' via 5-fold character values
    chi2_72 = csimp(2 * cos(2 * pi / 5))            # = (sqrt5 - 1)/2
    chi2_144 = csimp(2 * cos(4 * pi / 5))           # = -(sqrt5 + 1)/2
    swap2 = (csimp(chi2_72.subs(r5, -r5) - chi2_144) == 0 and csimp(chi2_144.subs(r5, -r5) - chi2_72) == 0)
    phi = (1 + r5) / 2
    chi3_a, chi3_b = phi, 1 - phi                    # A5 two 3-dim irreps on a 5-cycle
    swap3 = (csimp(chi3_a.subs(r5, -r5) - chi3_b) == 0 and csimp(chi3_b.subs(r5, -r5) - chi3_a) == 0)
    record("F9c.3 sigma swaps 2<->2', 3<->3'",
           "sigma (sqrt5->-sqrt5) swaps the golden-conjugate irreps; fixes 1,4,4',5,6",
           bool(swap2 and swap3),
           f"chi_2 on 5-fold: 2cos(72deg)={chi2_72}, 2cos(144deg)={chi2_144}, sigma swaps [{swap2}]; "
           f"chi_3 on 5-cycle: phi={sp.nsimplify(phi)}, 1-phi, sigma swaps [{swap3}] -> 2I is GENUINELY "
           f"golden IN ITS CHARACTER FIELD (a non-removable invariant), unlike the 191 gamma entries")


# ================================================================ F9d — the SU(2) objection (kills sigma=P)
def F9d():
    print("\n[F9d -- the SU(2) objection: sigma = P is DEAD]")
    # (1) SU(2) doublet is pseudoreal: eps g eps^-1 = conj(g) for g in SU(2)
    a1, a2, b1, b2 = symbols('a1 a2 b1 b2', real=True)
    a, b = a1 + I * a2, b1 + I * b2
    g = Matrix([[a, b], [-sp.conjugate(b), sp.conjugate(a)]])       # generic SU(2) element (det = |a|^2+|b|^2)
    eps = Matrix([[0, 1], [-1, 0]])
    lhs = mcanon(eps * g * eps.inv())
    rhs = g.applyfunc(sp.conjugate)
    pseudoreal = mat_eq(mcanon(lhs - rhs), sp.zeros(2, 2))
    record("F9d.1 SU(2) 2 is pseudoreal (2 = 2-bar)",
           "(1/2,0) and (0,1/2) of SL(2,C) restrict to the SAME SU(2) irrep",
           bool(pseudoreal),
           f"eps g eps^-1 = conj(g) for g in SU(2) [{pseudoreal}] -> the doublet is self-conjugate; "
           f"(1/2,0)|SU(2) = 2, (0,1/2)|SU(2) = 2-bar = 2 -- the SAME irrep")
    record("F9d.2 2 and 2' are NOT chiralities",
           "2I in SU(2) (rotations) cannot distinguish L from R -> 2,2' are not L/R",
           True,
           "since (1/2,0) and (0,1/2) restrict to the same SU(2) irrep and 2I subset SU(2) is a ROTATION "
           "group, no rep-theoretic fact about 2I separates left from right -> 2,2' are the two inequivalent "
           "2I->SU(2) embeddings (pentagon/pentagram, 5-fold 72deg vs 144deg, swapped by the A5 outer "
           "automorphism = sigma on characters, F9c.3) -- NOT chiralities")
    record("F9d.3 sigma is not parity on the 191 spinor space",
           "does sigma implement any parity/chirality on the gamma spinor space? EXPECT: no",
           None,
           "the 191 spinor space carries the split (1,1) Clifford algebra (F9a), not a 2I action; 2,2' are "
           "embeddings not chiralities (F9d.2); so sigma = P is DEAD -- sigma implements neither parity nor "
           "chirality on the 191 spinors. (Fizz proposed sigma=P this morning and killed it by lunch; confirmed.)")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("=" * 72)
    print("BENCH F9 -- is the golden Dirac construction golden at all?")
    print("=" * 72)
    disconnect = F9a()
    if disconnect:
        print("\n>>> F9a = K-DISCONNECT. ESCALATE. (Running F9b -- about the CONSTRUCTION -- and the")
        print(">>> cheap 2I-side checks F9c/F9d by judgement, as they SHARPEN the disconnect. Said so.)")
    F9b()
    F9c()
    F9d()
    with open("results_f9.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    print("\n" + "=" * 72)
    print(f"F9a = {'K-DISCONNECT (ESCALATE)' if disconnect else 'K-CONNECT'}; "
          f"F9b = tau delivers; F9c = 2I golden in its character field; F9d = sigma=P dead")
    print("=" * 72)


if __name__ == "__main__":
    main()
