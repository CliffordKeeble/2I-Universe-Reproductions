"""
bench_f8.py — F8: Brick 4 re-posed on the surviving pairing K2 (Paper 226).

Source of record: mr_code_brief_F8_brick4_on_K2.md (Fizz, 11 Jul 2026).
Pre-registration: PRE-REGISTRATION_f8.md (committed 4902a01, before this verdict).

F7 retracted the no-go (K2 gives a real nondegenerate action). F8 asks the F4-level
existence questions on K2, all previously gated off:
  F8a  conservation of  <psi1,psi2> = int (sigma psi1)^T K2 Gs psi2   (the gate)
  F8b  the K4 frequency signature  sign (sigma u_pm)^T K2 Gs u_pm  + the sigma:+->- map   (crux)
  F8c  ev1-unitarity: does ev1(S) need INDEPENDENT inf2 data?  (K5)
  §6   genericity sqrt2/sqrt3 + the named tension (K2 fails 226 §5.2 frontier), carried RAW.

sigma-odd mass axiom: m^2 = sqrt5 * mu^2 (mu rational).  ENCODED m = sqrt(sqrt5)*mu, so
m is REAL at inf1 (sqrt5>0) and IMAGINARY at inf2 (sqrt5<0).  The pairing uses sigma (Galois),
NOT dagger, so it need not be real -- reported exact, sign extracted per Hack (-i if pure imaginary),
NOT forced.
"""
import os
import sys
import csv
import sympy as sp
from sympy import sqrt, I, symbols, Matrix, diag, Rational, simplify, expand, im, re
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, mcanon, mat_eq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

Gs, Ga, Z = G_seed, G_adj, Z2
K2 = Matrix([[0, 1], [1, 0]])
mu, k = symbols('mu k', real=True, positive=True)
m = sqrt(r5) * mu                       # m^2 = sqrt5 * mu^2  (sigma-ODD);  real at inf1, imaginary at inf2
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed is True else ("FAIL" if passed is False else "INFO")
    ROWS.append((item, expect, f"{tag} - {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


def sig(x):
    """Galois sigma on a sympy expr/Matrix: sqrt5 -> -sqrt5 (handles sqrt(sqrt5) too)."""
    if hasattr(x, 'applyfunc'):
        return x.applyfunc(lambda e: e.subs(r5, -r5))
    return x.subs(r5, -r5)


# ================================================================ §3 F8a — conservation on K2
def F8a_solve(A, m_expr, label, quiet=False):
    """Find the conserved spatial current: solve d0 D + d1 J1 = 0 on shell for J1 = (sigma p1)^T B p2,
       B a general 2x2. Derivative terms fix B; the MASS residual then decides conservation."""
    c = symbols('c0:8')
    p1, p2 = Matrix([c[0], c[1]]), Matrix([c[2], c[3]])
    d1p1, d1p2 = Matrix([c[4], c[5]]), Matrix([c[6], c[7]])
    bb = symbols('b0:4')
    B = Matrix([[bb[0], bb[1]], [bb[2], bb[3]]])
    Gsi = Gs.inv()
    d0p1 = sp.expand(Gsi * (-Ga * d1p1 - m_expr * Z * p1))
    d0p2 = sp.expand(Gsi * (-Ga * d1p2 - m_expr * Z * p2))
    D_dot = sp.expand((sig(d0p1).T * A * Gs * p2)[0] + (sig(p1).T * A * Gs * d0p2)[0])
    J_div = sp.expand((sig(d1p1).T * B * p2)[0] + (sig(p1).T * B * d1p2)[0])
    total = sp.expand(D_dot + J_div)
    deriv = {c[4], c[5], c[6], c[7]}
    deriv_part = sum((t for t in total.as_ordered_terms() if t.free_symbols & deriv), sp.Integer(0))
    mass_part = sp.expand(total - deriv_part)
    # solve deriv_part == 0 (all field*derivative monomials) for B
    eqs = []
    for a in (4, 5, 6, 7):
        for b in (0, 1, 2, 3):
            eqs.append(sp.expand(deriv_part.coeff(c[a] * c[b])))
    eqs = [e for e in eqs if e != 0]
    sol = sp.solve(eqs, list(bb), dict=True)
    if not sol:
        if not quiet:
            record(f"F8a conservation [{label}]", "conserved current exists?", False,
                   f"no spatial current cancels the derivative terms -> not conserved")
        return False, None, mass_part
    Bsol = {k: v for k, v in sol[0].items()}
    residual = csimp(mass_part.subs(Bsol))
    conserved = (residual == 0)
    if not quiet:
        record(f"F8a conservation [{label}]",
               "conserved current exists AND mass residual cancels on shell",
               conserved,
               f"current B solved from derivative terms; mass residual = {residual} "
               f"-> conserved [{conserved}]")
    return conserved, csimp(B.subs(Bsol)) if len(Bsol) == 4 else "B partly free", residual


def structural_root():
    """The root cause behind the non-conservation and the complex signature."""
    from sympy import Symbol
    x = Symbol('x')
    # (1) sigma-odd mass m^2 = sqrt5 mu^2 -> m = 5^(1/4) mu is NOT in K = Q(sqrt5)
    mp = sp.minimal_polynomial(sp.root(5, 4), x)          # x^4 - 5
    deg = sp.degree(mp, x)
    m_notin_K = (deg == 4)
    record("F8 root (1): sigma-odd mass forces m not in K",
           "the first-order mass coefficient leaves the field K",
           None,
           f"m = 5^(1/4) mu; minimal polynomial over Q = {mp} (degree {deg} > 2) -> m NOT in K=Q(sqrt5) "
           f"[{m_notin_K}]. Cliff's axiom M^2 sigma-even + dispersion M^2=m^2 s^2/sqrt5 => m^2=sqrt5 mu^2 "
           f"=> m irrational-of-degree-4: the massive golden Dirac operator is NOT over K.")
    # (2) phase structure: sigma fixes i (no U(1) off locus); tau flips i (U(1) on the locus, sigma=tau)
    al = symbols('alpha', real=True)
    sigma_fixes_phase = (sp.simplify(sigma(Matrix([sp.exp(I * al)]))[0] - sp.exp(I * al)) == 0)
    tau_flips_phase = (sp.simplify(sp.conjugate(sp.exp(I * al)) - sp.exp(-I * al)) == 0)
    record("F8 root (2): the sigma-adjoint carries no phase conjugation",
           "psibar=(sigma psi)^T A -> e^{+i a} psibar, so <psi,psi> -> e^{2i a}: no U(1) charge off locus",
           None,
           f"sigma fixes e^(i a) [{sigma_fixes_phase}] (no phase flip) -> off the physical locus the "
           f"sigma-pairing is a BILINEAR form (scales e^(2i a)), not a conserved charge. ON the tau-locus "
           f"sigma=tau and tau flips e^(i a)->e^(-i a) [{tau_flips_phase}].")
    # (3) but U(1) is ALSO broken ON the tau-locus: e^{i a} psi leaves the locus (Majorana-like)
    f, g = symbols('f g', real=True)
    comp = f + I * r5 * g                                   # a tau-constrained component (p real, q imag)
    rot = sp.expand(sp.exp(I * al) * comp)
    p_part = rot.subs(r5, 0)                                # sqrt5-free part = "p"; tau needs it real
    stays_on_locus = (sp.simplify(sp.im(p_part)) == 0)
    record("F8 root (3): U(1) is broken on the locus too (Majorana-like)",
           "e^{i a} psi leaves the tau-locus -> no U(1) even on-shell -> no conserved Dirac charge",
           None,
           f"e^(i a)*(f + i sqrt5 g) has p-part e^(i a) f with Im = f*sin(a) != 0 -> NOT tau-constrained "
           f"[stays_on_locus={stays_on_locus}]. The tau-reality constraint is a MAJORANA condition (sigma-real); "
           f"such fields carry NO conserved U(1) current -- consistent with the nonzero mass residual above. "
           f"=> the 'conserved pairing' of F8a does not exist for the massive golden field. GATE FAILS.")


def F8a_locus(m_expr, mlabel):
    """Conservation on the PHYSICAL tau-locus (K(i)-valued, real f,g): does the construction conserve
       IN-field?  Separates 'axiom's fault (m not in K)' from 'construction's fault (Majorana)'."""
    r = symbols('a1 b1 a2 b2 a3 b3 a4 b4 e1 h1 e2 h2 e3 h3 e4 h4', real=True)
    p1 = Matrix([r[0] + I * r5 * r[1], r[2] + I * r5 * r[3]])        # tau-constrained psi1
    p2 = Matrix([r[4] + I * r5 * r[5], r[6] + I * r5 * r[7]])        # tau-constrained psi2
    d1p1 = Matrix([r[8] + I * r5 * r[9], r[10] + I * r5 * r[11]])    # tau-constrained d1 psi1
    d1p2 = Matrix([r[12] + I * r5 * r[13], r[14] + I * r5 * r[15]])
    bb = symbols('B0:4')
    B = Matrix([[bb[0], bb[1]], [bb[2], bb[3]]])
    Gsi = Gs.inv()
    d0p1 = sp.expand(Gsi * (-Ga * d1p1 - m_expr * Z * p1))
    d0p2 = sp.expand(Gsi * (-Ga * d1p2 - m_expr * Z * p2))
    D_dot = sp.expand((sig(d0p1).T * K2 * Gs * p2)[0] + (sig(p1).T * K2 * Gs * d0p2)[0])
    J_div = sp.expand((sig(d1p1).T * B * p2)[0] + (sig(p1).T * B * d1p2)[0])
    total = sp.expand(D_dot + J_div)
    dset = set(r[8:])
    deriv_part = sum((t for t in total.as_ordered_terms() if t.free_symbols & dset), sp.Integer(0))
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
    record(f"F8a on tau-locus [{mlabel}]",
           "does the construction conserve IN-field (physical locus)?",
           conserved,
           f"tau-constrained solutions, mass {mlabel}; mass residual after solving current = {residual} "
           f"-> conserved [{conserved}]")
    return conserved


def F8a_diagnostic():
    """2x2: pairing {K2, A0=diag(r5,1)} x mass {sigma-odd, sigma-even}, WITH the current solved for."""
    A0 = diag(r5, 1)
    m_odd, m_even = sqrt(r5) * mu, mu
    print("     conservation 2x2 (pairing x mass parity), current solved each cell:")
    tab = {}
    for An, A in (("K2", K2), ("A0", A0)):
        for mn, me in (("m2-odd", m_odd), ("m2-even", m_even)):
            ok, Bsol, res = F8a_solve(A, me, f"{An}/{mn}", quiet=True)
            tab[(An, mn)] = ok
            print(f"        [{An:3s} | {mn:8s}] conserved = {ok}   (mass residual = {res})")
    record("F8a diagnostic (pairing x mass parity)",
           "with the correct current: does the mass residual cancel?",
           None,
           f"conserved table = {tab}")
    return tab


# ================================================================ §4 F8b — the K4 signature
def onshell_E():
    M2 = sp.simplify(m**2 * 1**2 / r5)      # s=1 -> M^2 = mu^2
    return sqrt(k**2 + M2), M2               # E, M^2


def u_pm(sign, s=sp.Integer(1)):
    E, _ = onshell_E()
    return Matrix([I * (sign * E + k), m * s]), sign * E


def pairing(ua, ub, A):
    """(sigma ua)^T A Gs ub  -- exact K(i,m) element."""
    return sp.expand((sig(ua).T * A * Gs * ub)[0])


def read_sign(val):
    """Report (real, imag) parts and an extracted sign per Hack (-i if pure imaginary)."""
    val = sp.simplify(val)
    rp, ip = sp.simplify(sp.re(val)), sp.simplify(sp.im(val))
    if ip == 0:
        kind = "REAL"
        signexpr = rp
    elif rp == 0:
        kind = "PURE-IMAG (extract -i)"
        signexpr = sp.simplify(-I * val)      # -i * (i*something) = real
    else:
        kind = "COMPLEX (real+imag both nonzero)"
        signexpr = None
    return kind, rp, ip, signexpr


def F8b(A, label):
    print(f"\n  -- F8b K4 signature on {label} --")
    up, Ep = u_pm(+1)
    um, Em = u_pm(-1)
    verdicts = {}
    for name, u in (("u+", up), ("u-", um)):
        P = pairing(u, u, A)
        for place, Pv in (("inf1", P), ("inf2", sig(P))):
            kind, rp, ip, se = read_sign(Pv)
            verdicts[(name, place)] = (kind, sp.simplify(Pv), rp, ip, se)
            record(f"F8b <{name},{name}> at {place} [{label}]",
                   "opposite definite signs on +/- (K4-PASS) vs same/intra-indefinite (K4-FIRE)",
                   None,
                   f"pairing = {sp.simplify(Pv)}  [{kind}]; Re={rp}, Im={ip}"
                   + (f", extracted sign ~ {se}" if se is not None else ""))
    # K4 map: is sigma(u+) a -frequency solution?  compare sigma(u+) to u- (inf2 reading) up to scale
    su = sig(up)                                # Galois image of the +freq inf1 solution
    # -freq inf2 solution:
    um_inf2 = sig(um)
    # proportional?  cross-ratio test:  su[0]*um_inf2[1] - su[1]*um_inf2[0] == 0
    prop_minus = (csimp(su[0] * um_inf2[1] - su[1] * um_inf2[0]) == 0)
    up_inf2 = sig(up)                           # (this is su itself)
    # compare instead sigma(u+) to the +freq vs -freq at the SAME (inf1) frame by re-reading:
    prop_to_uminus = (csimp(su[0] * um[1] - su[1] * um[0]) == 0)
    prop_to_uplus = (csimp(su[0] * up[1] - su[1] * up[0]) == 0)
    record(f"F8b K4 map [{label}]",
           "sigma(+freq inf1 solution) lands in the -freq sector?",
           None,
           f"sigma(u+) proportional to u- [{prop_to_uminus}] / to u+ [{prop_to_uplus}] "
           f"(the frequency-flip claim of C=sigma)")
    return verdicts, prop_to_uminus


# ================================================================ §5 F8c — ev1 unitarity (dof count)
def F8c():
    # tau-constrained field at inf1: psi_i = f_i + i sqrt5 g_i, 4 real dof (f1,g1,f2,g2).
    # sigma = tau on the physical locus (F3b/F5) => sigma(psi) is DETERMINED by (f,g), not independent.
    f1, g1, f2, g2 = symbols('f1 g1 f2 g2', real=True)
    psi = Matrix([f1 + I * r5 * g1, f2 + I * r5 * g2])
    sig_psi = sig(psi)                          # inf2 reading
    tau_psi = psi.applyfunc(lambda e: e.conjugate())   # complex conjugate (tau)
    sig_eq_tau = mat_eq(mcanon(sig_psi), mcanon(tau_psi))
    # independent real dof at inf1: 4 (f1,g1,f2,g2). inf2 (sigma psi) is a function of the SAME 4.
    n_inf1 = 4
    # benign iff sigma psi carries no NEW real dof beyond the 4 (i.e. determined by inf1 data)
    extra = set()
    for comp in sig_psi:
        extra |= (comp.free_symbols - {f1, g1, f2, g2, r5})
    benign = (len(extra) == 0) and sig_eq_tau
    record("F8c ev1-unitarity (dof count)",
           "ev1(S) a functional of inf1 data alone (benign) vs needs independent inf2 data (malign->K5)",
           benign,
           f"tau-constraint: 4 real dof at inf1; sigma(psi) determined by them (no new symbols {extra}); "
           f"sigma=tau on the physical locus [{sig_eq_tau}] -> inf2 is the CONJUGATE of inf1 data, "
           f"benign like ordinary Dirac psi-bar~psi-dagger [{benign}] -> K5 does NOT fire")
    return benign


# ================================================================ §6 genericity + named tension
def genericity_and_tension():
    # signature over sqrt2, sqrt3 (d != 1 mod 4): same complex structure?
    print("\n  -- §6 genericity sqrt2/sqrt3 (signature structure) --")
    res = {}
    for D in (2, 3, 5):
        rD = sqrt(D)
        GsD = Matrix([[0, 1], [rD, 0]])
        K2D = Matrix([[0, 1], [1, 0]])
        mD = sqrt(rD) * mu
        E = sqrt(k**2 + mu**2)

        def sigD(x):
            return x.applyfunc(lambda e: e.subs(rD, -rD)) if hasattr(x, 'applyfunc') else x.subs(rD, -rD)
        up = Matrix([I * (E + k), mD])
        P = sp.expand((sigD(up).T * K2D * GsD * up)[0])
        rp, ip = sp.simplify(sp.re(P)), sp.simplify(sp.im(P))
        res[D] = (sp.simplify(P), rp, ip)
        print(f"     D={D}: <u+,u+> = {sp.simplify(P)}  (Re={rp}, Im={ip})")
    # structure identical (real part ~ -rD(E+k)^2, imag ~ +rD mu^2) across d?
    same_structure = all((res[D][1] != 0 and res[D][2] != 0) for D in (2, 3, 5))
    record("§6 genericity (signature structure over d)",
           "identical signature structure for sqrt2, sqrt3, sqrt5 (no d-dependence -> generic)",
           None,
           f"complex (Re<0 kinetic, Im>0 mass) for all d in {{2,3,5}} [{same_structure}]; "
           f"no d-specific sign flip -> real-quadratic-generic (not golden)")
    # named tension, carried RAW (NOT adjudicated; Fizz's §7.1 resolution off-limits)
    def uniform_eps(Cmat, Gamma):
        lhs = mcanon(Cmat * sig(Gamma) * Cmat.inv())
        ratios = []
        for i in range(2):
            for j in range(2):
                a, b = sp.simplify(lhs[i, j]), sp.simplify(Gamma.T[i, j])
                if a == 0 and b == 0:
                    continue
                if a == 0 or b == 0:
                    return None
                ratios.append(sp.simplify(a / b))
        e0 = ratios[0]
        return e0 if all(sp.simplify(r - e0) == 0 for r in ratios) else None
    frontier_K2 = uniform_eps(K2, Gs) is not None and uniform_eps(K2, Ga) is not None
    record("§6 named tension (RAW, not adjudicated)",
           "K2 fails 226 §5.2 frontier while passing reality -> hand to Mr A / Tor",
           None,
           f"K2 satisfies 226 §5.2 uniform-eps frontier [{frontier_K2}] (False = fails) while F7 reality "
           f"PASSES -> the pairing form and the conjugation matrix disagree about sigma. Reported raw; "
           f"Fizz's operator/amplitude resolution deliberately NOT tested (her §7.1, Mr Code off-limits)")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("=" * 72)
    print("BENCH F8 -- Brick 4 re-posed on K2 (conservation, K4 signature, ev1-unitarity)")
    print("=" * 72)
    print("\n[§3 F8a -- conservation on K2 (the gate)]")
    consK2, Bsol, resK2 = F8a_solve(K2, m, "K2 (sigma-odd mass, the axiom)")
    tab = F8a_diagnostic()
    print("\n     conservation on the physical tau-locus (axiom vs construction):")
    loc_even = F8a_locus(mu, "sigma-EVEN (m in K)")
    loc_odd = F8a_locus(m, "sigma-ODD (m = 5^(1/4) mu, the axiom)")
    print("\n[F8 structural root]")
    structural_root()
    if not consK2:
        print("\nESCALATION (§8): conservation FAILS on K2 (gate) -> do not bank downstream; root cause below.")
    print("\n[§4 F8b -- the K4 signature (crux)]")
    vK2, flipK2 = F8b(K2, "K2")
    v5, flip5 = F8b(mcanon(r5 * K2), "sqrt5*K2")
    # §4.5 do the two b=c representatives agree?
    agree = all(vK2[key][0] == v5[key][0] for key in vK2)   # same kind (REAL/IMAG/COMPLEX) at each slot
    record("F8b.5 representative agreement",
           "K2 and sqrt5*K2 give the SAME signature (else signature not a property of the theory)",
           agree,
           f"signature-kind agrees across K2 and sqrt5*K2 at every (mode,place) [{agree}]"
           + ("" if agree else " -> FLAG: representatives disagree, escalate"))
    print("\n[§5 F8c -- ev1-unitarity (K5)]")
    benign = F8c()
    genericity_and_tension()

    # ---- verdict: the §3 GATE FAILED -> escalation; §4/§5 are downstream (provisional)
    print("\n  -- F8 VERDICT --")
    verdict = "ESCALATE (§3 gate fails): the symmetric Dirac pairing K4 needs does not exist"
    outcome = (
        "GATE FAIL (§3, §8 escalation). The F8a *symmetric* conserved pairing does NOT exist -- TWO "
        "independent structural facts: "
        "(A) CONSTRUCTION (not just the axiom): conservation of (sigma psi1)^T K2 Gs psi2 fails for K2 AND A0, "
        "for sigma-even AND sigma-odd mass, on general fields AND on the physical tau-locus with an IN-field "
        "sigma-even mass -- because the sigma-adjoint fixes the phase and the tau-reality constraint is a "
        "MAJORANA (sigma-real) condition, so U(1) is broken and a sigma-real field carries no conserved "
        "symmetric Dirac charge. (What IS conserved is the ANTISYMMETRIC symplectic form Omega of F7 -- the "
        "'Grassmann-odd shadow' the brief itself names -- not the symmetric pairing K4 asks the sign of.) "
        "(B) AXIOM: Cliff's antiparticle axiom (M^2 sigma-even) + the golden dispersion (M^2=m^2 s^2/sqrt5) "
        "force m=5^(1/4)mu NOT in K -> the massive operator leaves the field (a second, separate obstruction). "
        "DOWNSTREAM (provisional, not banked): F8b signature COMPLEX at both places (symptom of B), no clean "
        "opposite-definite pair -> K4 does NOT cleanly PASS; F8b.5 reps agree in kind; F8c K5 benign on locus. "
        "Named tension (K2 fails 226 §5.2 frontier) carried RAW; Fizz's operator/amplitude resolution NOT tested.")
    record("F8 VERDICT", "K4-PASS / K4-FIRE / K-DEGENERATE / ESCALATE", None, outcome)

    with open("results_f8.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    print("\n" + "=" * 72)
    print(f"F8: conservation gate = FAIL (escalate); K4 signature = COMPLEX (no clean PASS); K5 benign(on locus) = {benign}")
    print("ESCALATION: the sigma-odd mass axiom + the golden dispersion force m NOT in K, and the")
    print("sigma-real (Majorana) structure carries no conserved U(1) charge. Hand to Cliff/Tor/Fizz/Mr A.")
    print("=" * 72)


if __name__ == "__main__":
    main()
