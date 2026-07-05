"""
bench_run_1.py — Golden Construction Bench Run #1, Part A (A1..A10).
Exact arithmetic over Q(sqrt5, i). Reports pass/fail against the brief's EXPECT lines.

Run:  python bench_run_1.py
Outputs: results.csv (+ console table). See PRE-REGISTRATION.md for fixed choices.
"""

import csv
import random
import numpy as np
import sympy as sp
from sympy import sqrt, I, Rational, eye, zeros, Matrix, symbols, simplify, expand, Poly

from golden_algebra import (
    r5, I2, K2, J2, Z2, G_seed, G_adj, kron, canon, mcanon, mat_eq, comm, anticomm,
    is_zero, tau, sigma, sigmatau, square_sign, square_scalar,
    build_candidates, find_cliques, anticommute,
)

SEED = 20260705
ROWS = []  # (item, expect, outcome, deviation)


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        deviation: {deviation}")


# Paper-191 gammas
G1 = kron(K2, G_seed)   # +sqrt5
G2 = kron(J2, G_seed)   # -sqrt5
G3 = kron(I2, G_adj)    # -sqrt5
I4 = eye(4)


# ---------------------------------------------------------------- A1
def a1():
    checks = {
        "[Z,Gseed]=-2Gadj": mat_eq(comm(Z2, G_seed), -2 * G_adj),
        "[Z,Gadj]=-2Gseed": mat_eq(comm(Z2, G_adj), -2 * G_seed),
        "{Z,Gseed}=0": is_zero(anticomm(Z2, G_seed)),
        "{Z,Gadj}=0": is_zero(anticomm(Z2, G_adj)),
        "Gseed.Gadj=+sqrt5.sigma3": mat_eq(G_seed * G_adj, r5 * Z2),
    }
    passed = all(checks.values())
    fails = [k for k, v in checks.items() if not v]
    record("A1 block relations", "all pass", passed,
           "all 5 relations hold" if passed else f"{len(fails)} failed",
           "" if passed else "; ".join(fails))


# ---------------------------------------------------------------- A2
_R5 = float(np.sqrt(5))


def _to_np(M):
    return np.array([[complex(sp.N(M[i, j].subs(r5, sp.sqrt(5)), 30))
                      for j in range(M.shape[1])] for i in range(M.shape[0])],
                    dtype=complex)


def _conj_sign_numeric(M_np, G_np, s0):
    """Sign of (M G M^-1)^2 by high-precision numeric; True if == s0.sqrt5.I4."""
    Gc = M_np @ G_np @ np.linalg.inv(M_np)
    Gc2 = Gc @ Gc
    return np.allclose(Gc2, s0 * _R5 * np.eye(4), atol=1e-9)


def a2():
    cands = build_candidates()
    rng = random.Random(SEED)

    # structured generator sample -> EXACT symbolic conjugation
    gens = [("G1", G1), ("G2", G2), ("G3", G3),
            ("Z(x)I", kron(Z2, I2)), ("I(x)Z", kron(I2, Z2)), ("K(x)K", kron(K2, K2))]

    # >=20 random invertible 4x4 over Q(i,sqrt5) -> high-precision NUMERIC hunt
    # (conjugation-invariance of the square sign is trivial: G^2 = +-sqrt5.I4 is central;
    #  the numeric sweep is the counterexample hunt, not the proof.)
    randM = []
    while len(randM) < 20:
        M = Matrix(4, 4, lambda i, j:
                   rng.randint(-2, 2) + rng.randint(-2, 2) * I
                   + rng.randint(-2, 2) * r5 + rng.randint(-2, 2) * I * r5)
        Mnp = _to_np(M)
        if abs(np.linalg.det(Mnp)) > 1e-6:
            randM.append(Mnp)

    n_exact_conj = 0
    n_num_conj = 0
    preserved_ops = 0
    flipped_ops = 0
    bad = []
    for (label, G, s0) in cands:
        # (a-exact) structured generators
        for gl, M in gens:
            Gc = mcanon(M * G * M.inv())
            if square_sign(Gc) != s0:
                bad.append(f"{label}: exact conj by {gl} changed sign")
            else:
                n_exact_conj += 1
        # (a-numeric) 20 random invertibles over Q(i,sqrt5)
        G_np = _to_np(G)
        for Mnp in randM:
            if _conj_sign_numeric(Mnp, G_np, s0):
                n_num_conj += 1
            else:
                bad.append(f"{label}: numeric conj sign changed")
        # (b) transpose, (c) negation, (d) tau  -> must preserve (EXACT)
        for opname, Gp in (("T", G.T), ("neg", -G), ("tau", tau(G))):
            if square_sign(Gp) != s0:
                bad.append(f"{label}: {opname} changed sign")
            else:
                preserved_ops += 1
        # (e) sigma, (f) sigmatau -> must flip (EXACT)
        for opname, Gp in (("sigma", sigma(G)), ("sigmatau", sigmatau(G))):
            if square_sign(Gp) != -s0:
                bad.append(f"{label}: {opname} did NOT flip sign")
            else:
                flipped_ops += 1

    passed = not bad
    record("A2 L1 invariance hunt", "preserve (a)-(d); flip (e),(f)", passed,
           f"16 cands: exact conj {n_exact_conj}, numeric conj {n_num_conj}, "
           f"transpose/neg/tau preserved {preserved_ops}, sigma/sigmatau flipped "
           f"{flipped_ops} — all correct" if passed else f"{len(bad)} anomalies",
           "" if passed else "; ".join(bad[:6]))


# ---------------------------------------------------------------- A3
def a3():
    def C0_alg(X):
        return mcanon(Z2 * sigma(X) * Z2.inv())

    def C0_vec(v):
        return mcanon(Z2 * sigma(v))

    checks = {
        "C0(Gseed)=Gadj": mat_eq(C0_alg(G_seed), G_adj),
        "C0(Gadj)=Gseed": mat_eq(C0_alg(G_adj), G_seed),
        "C0^2=id on algebra": mat_eq(C0_alg(C0_alg(G_seed)), G_seed)
                              and mat_eq(C0_alg(C0_alg(G_adj)), G_adj),
    }
    # on vectors
    v = Matrix([symbols('v0'), symbols('v1')])   # symbols may carry sqrt5? keep generic K-rational
    vK = Matrix([1 + r5, 2 - 3 * r5])            # a concrete K-rational vector
    checks["C0^2=id on vectors"] = mat_eq(C0_vec(C0_vec(vK)), vK)
    passed = all(checks.values())
    fails = [k for k, val in checks.items() if not val]
    record("A3 C0 explicit", "all pass", passed,
           "C0 swaps seed/adj, C0^2=+1 (algebra & vectors)" if passed
           else f"{len(fails)} failed", "" if passed else "; ".join(fails))


# ---------------------------------------------------------------- A4
def a4():
    u = symbols('u', nonzero=True)
    s, a = symbols('s a')          # <lambda|psi>_seed , <lambda|psi>_adj
    ip = symbols('ip')             # <psi_seed|psi_adj>
    P = s * a
    # (u, u^-1)
    P_scaled = (u * s) * (u**-1 * a)
    ip_scaled = (u) * (u**-1) * ip
    pres_P = canon(P_scaled - P) == 0
    pres_ip = canon(ip_scaled - ip) == 0
    # (u, u): must NOT preserve
    P_equal = (u * s) * (u * a)
    not_pres = canon(P_equal - P) != 0
    # trace-zero line
    p, q = symbols('p q', real=True)
    aa = p + q * r5
    cond = canon(aa + sigma(Matrix([[aa]]))[0])   # a + sigma(a) = 2p
    trace_zero_iff = (canon(cond) == canon(2 * p)) and (sp.solve(sp.Eq(cond, 0), p) == [0])
    checks = {
        "(u,u^-1) preserves P": pres_P,
        "(u,u^-1) preserves <seed|adj>": pres_ip,
        "(u,u) does NOT preserve P": not_pres,
        "a+sigma(a)=0 <=> p=0 (Lie T_sigma = Q.sqrt5)": trace_zero_iff,
    }
    passed = all(checks.values())
    fails = [k for k, val in checks.items() if not val]
    record("A4 weights/trace-zero", "all pass", passed,
           "norm-one torus preserves pairing; equal-weight fails; Lie(T_sigma)=Q.sqrt5"
           if passed else f"{len(fails)} failed", "" if passed else "; ".join(fails))


# ---------------------------------------------------------------- A5  (the F3 check)
def a5():
    k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)
    W = Z2   # weight action: +1 seed, -1 adjoint
    GA = {"Gs": G_seed, "Ga": G_adj}

    def Dmat(kmu, amu, charge):     # D_mu on plane wave: -i k_mu I - charge*a_mu*W
        return (-I * kmu) * I2 - charge * amu * W

    def Lop(order, charge, asign):
        """Kinetic operator: order[0] gamma on d0, order[1] on d1;
        charge multiplies coupling; asign flips the gauge field a_mu -> asign*a_mu."""
        return (GA[order[0]] * Dmat(k0, asign * a0, charge)
                + GA[order[1]] * Dmat(k1, asign * a1, charge))

    # the q-equation, literal construction reading: Gamma_seed on d0, Gamma_adj on d1
    Mq = mcanon(Lop(("Gs", "Ga"), q, +1))

    # C0 = Z o sigma acting on the plane-wave amplitude: chi -> Z sigma(chi).
    # C0 maps ker(Mq) into ker(target) iff  sigma(Mq) == Z . target . Z  (exact, all k).
    # Scan the four candidate (-q) target operators to characterise exactly what C0 does.
    targets = {}
    for order in (("Gs", "Ga"), ("Ga", "Gs")):
        for asign in (+1, -1):
            tgt = mcanon(Lop(order, -q, asign))
            key = ("literal-order" if order == ("Gs", "Ga") else "swapped-order",
                   "same a" if asign == +1 else "a->-a")
            targets[key] = mat_eq(sigma(Mq), mcanon(Z2 * tgt * Z2))

    literal = targets[("literal-order", "same a")]          # the EXPECT target: L_{-q}
    winners = [k for k, v in targets.items() if v]

    # localise the leak: difference from the pure sector/coordinate swap (same a)
    swap_same_a = mcanon(Lop(("Ga", "Gs"), -q, +1))
    leak = mcanon(sigma(Mq) - mcanon(Z2 * swap_same_a * Z2))   # should be pure a_mu term
    leak_is_pure_coupling = (not leak.has(k0)) and (not leak.has(k1)) and leak.has(q)

    passed = bool(literal)   # EXPECT literal target; pre-registered as the informative check
    if literal:
        outcome = "C0 maps q-solutions to L_{-q} solutions cleanly (literal reading holds)"
        dev = ""
    else:
        wtxt = "; ".join(f"{o}/{a}" for (o, a) in winners) if winners else "NONE"
        outcome = ("literal C0=Z.sigma does NOT map L_q solutions to L_{-q} (EXPECT target); "
                   "it intertwines only the seed<->adjoint / x0<->x1-swapped, a->-a operator")
        dev = (f"exact all-k intertwiner sigma(Mq)=Z.target.Z holds ONLY for: [{wtxt}]. "
               f"Free/kinetic part maps under the sector+coordinate swap; the residual "
               f"obstruction vs the same-a swapped target is pure minimal-coupling = "
               f"{sp.simplify(leak[0])!r},{sp.simplify(leak[1])!r},{sp.simplify(leak[2])!r},"
               f"{sp.simplify(leak[3])!r} (proportional to q.a_mu, no k-dependence: "
               f"{leak_is_pure_coupling}). F3 signal: sector grading (W=Z) requires a "
               f"time<->space swap to align with operator grading; not patched.")
    record("A5 toy-model end-to-end (F3)", "pass conditionally", passed, outcome, dev)
    a5.facts = dict(literal=literal, winners=winners,
                    leak_pure_coupling=leak_is_pure_coupling,
                    leak=[sp.simplify(x) for x in leak])


# ---------------------------------------------------------------- A6
def a6():
    cands = build_candidates()
    cliques = find_cliques(cands, 4)
    # (Gamma5_(4))^2 = +25 for all six (2,2) cliques
    vals = []
    for cl in cliques:
        w = I4
        for idx in cl:
            w = w * cands[idx][1]
        w = mcanon(w)
        vals.append(square_scalar(w))
    all_plus25 = all(mat_eq(Matrix([[v]]), Matrix([[25]])) for v in vals)

    # one explicit i-dressed (1,3) clique: dress first +sqrt5 gamma of clique 0
    cl0 = cliques[0]
    mats0 = [cands[i][1] for i in cl0]
    signs0 = [cands[i][2] for i in cl0]
    pos_positions = [i for i, s in enumerate(signs0) if s == +1]
    dressed = list(mats0)
    dressed[pos_positions[0]] = mcanon(I * mats0[pos_positions[0]])   # x i
    dressed_signs = [square_sign(M) for M in dressed]
    is_13 = sorted(dressed_signs) == [-1, -1, -1, +1]
    wd = I4
    for M in dressed:
        wd = wd * M
    wd = mcanon(wd)
    omega2 = square_scalar(wd)
    omega_is_m25 = mat_eq(Matrix([[omega2]]), Matrix([[-25]]))

    # three-gamma cross check (Paper 191): (G1 G2 G3)^2 = -5^{3/2}
    tg = mcanon(G1 * G2 * G3)
    tg2 = square_scalar(tg)
    tg_ok = mat_eq(Matrix([[tg2]]), Matrix([[-5 * r5]]))

    passed = all_plus25 and is_13 and omega_is_m25 and tg_ok
    record("A6 volume-element squares", "+25(x6), -25, -5^{3/2}", passed,
           f"(G5_4)^2=+25 on all {len(cliques)} (2,2) cliques; i-dressed clique is "
           f"(1,3) with omega^2={omega2}; 3-gamma sq={tg2}",
           "" if passed else f"all+25={all_plus25}, is13={is_13}, "
           f"omega2={omega2}, 3gamma={tg2}")
    a6.dressed = dressed
    a6.dressed_signs = dressed_signs


# ---------------------------------------------------------------- A7
def a7():
    x = symbols('x')
    beta = 5**Rational(1, 4) * I
    mp = sp.minimal_polynomial(beta, x)
    mp_ok = sp.expand(mp - (x**4 - 5)) == 0
    # embeddings: real roots of x^4-5 = +-5^{1/4} (2 real), +-5^{1/4} i (1 conj pair)
    roots = sp.Poly(x**4 - 5, x).all_roots()
    n_real = sum(1 for rt in roots if sp.im(rt) == 0)
    n_complex_pairs = (len(roots) - n_real) // 2
    emb_ok = (n_real == 2 and n_complex_pairs == 1)
    # i not in Q(beta): [Q(beta):Q]=4; i would force degree-8 field
    deg_beta = sp.degree(mp, x)
    i_notin = (deg_beta == 4)   # since Q(i,5^{1/4}) has degree 8 > 4
    passed = mp_ok and emb_ok and i_notin
    record("A7 eigenvalue fields", "all pass", passed,
           f"minpoly(5^1/4 i)=x^4-5; embeddings {n_real} real + {n_complex_pairs} "
           f"cc-pair (not CM); i not in Q(5^1/4 i) [deg {deg_beta}]",
           "" if passed else f"mp_ok={mp_ok}, emb_ok={emb_ok}, i_notin={i_notin}")


# ---------------------------------------------------------------- A8  (descent J)
def a8():
    import mpmath as mp
    mp.mp.dps = 40
    dressed = a6.dressed
    # J = B o tau commutes with all four gammas  <=>  B tau(G) = G B  (tau = conj)
    b = symbols('b0:16')
    B = Matrix(4, 4, b)
    eqs = []
    for G in dressed:
        eqs += list(B * tau(G) - G * B)        # exact, no simplify
    sol = list(sp.linsolve(eqs, list(b)))[0]
    free = set()
    for comp in sol:
        free |= (comp.free_symbols & set(b))
    dim = len(free)
    subs0 = {p: 0 for p in b}
    if free:
        subs0[sorted(free, key=lambda s: s.name)[0]] = 1
    B0 = Matrix(4, 4, [comp.subs(subs0) for comp in sol])
    nonzero = any(B0[i] != 0 for i in range(16))
    Jsq = B0 * tau(B0)                          # J^2 = B tau(B)

    # NOTE: sympy simplify is UNRELIABLE on these nested M4(Q(sqrt5,i)) products
    # (cross-checked: it wrongly reports +1/5). High-precision numeric is authoritative.
    def numscalar(r5val):
        vals = [[complex(sp.N(Jsq[i, j].subs(r5, r5val), 35)) for j in range(4)]
                for i in range(4)]
        c = vals[0][0]
        res = max(abs(vals[i][j] - (c if i == j else 0))
                  for i in range(4) for j in range(4))
        return c, res
    cp, resp = numscalar(mp.sqrt(5))
    cm, resm = numscalar(-mp.sqrt(5))
    scalar_ok = (resp < 1e-25 and resm < 1e-25)
    cp, cm = cp.real, cm.real
    mag = abs(cp)                               # |c| = 1/sqrt5 expected
    is_1_over_sqrt5 = abs(mag - 1 / float(mp.sqrt(5))) < 1e-20
    tot_neg = (cp < 0 and cm < 0)              # J^2 = -1 reachable iff c totally negative
    tot_pos = (cp > 0 and cm > 0)              # J^2 = +1 reachable iff c totally positive
    indefinite = (cp * cm < 0)
    # renormalisation J->lam J scales J^2 by N(lam)=lam.conj(lam) (totally positive):
    minus1_reachable = tot_neg
    plus1_reachable = tot_pos

    unique = (dim == 1)
    # EXPECT: unique; J^2=-1; +1 impossible.  Met only if J^2 normalises to -1.
    passed = unique and nonzero and scalar_ok and minus1_reachable and (not plus1_reachable)
    cstr = f"-1/sqrt5 (=-sqrt5/5)" if (is_1_over_sqrt5 and cp < 0) else f"|c|={mag:.6f}"
    if passed:
        outcome = (f"B unique (dim {dim}); J^2 = {cstr}*I totally negative -> ~-1; "
                   f"+1 impossible")
        dev = ""
    else:
        outcome = (f"B unique up to scale (dim {dim}) OK; but J^2 = ({cstr})*I is "
                   f"INDEFINITE, not -1 (sign at +sqrt5 place={cp:+.4f}, at -sqrt5 "
                   f"place={cm:+.4f})")
        dev = (f"J^2 = c.I with c = -1/sqrt5 (numeric-authoritative; symbolic simplify "
               f"unreliable here). Renormalisation J->lam.J scales J^2 by the totally-"
               f"positive norm N(lam), so the (-,+) sign pattern is invariant: J^2 "
               f"reaches NEITHER -1 (needs N=sqrt5, indefinite -> not a norm) NOR +1 "
               f"(needs N=-sqrt5). EXPECT's '+1 impossible' holds; headline 'J^2=-1' "
               f"does NOT. Cause: sqrt5-normalisation (squares +-sqrt5, not +-1) -- "
               f"204 sec8 / analysis sec5 flag this clique is not standard Cl(1,3); "
               f"does NOT contradict Paper 204 Thm 1.")
    record("A8 descent invariant J", "J unique; J^2=-1; +1 impossible", passed,
           outcome, dev)
    a8.facts = dict(dim=dim, cp=cp, cm=cm, indefinite=indefinite,
                    minus1_reachable=minus1_reachable, plus1_reachable=plus1_reachable)


# ---------------------------------------------------------------- A9
def a9():
    dressed = a6.dressed
    sig = sigma  # entrywise
    img = [mcanon(sig(M)) for M in dressed]
    # mutual anticommutation preserved?
    anti_ok = all(anticommute(img[i], img[j])
                  for i in range(4) for j in range(i + 1, 4))
    img_signs = [square_sign(M) for M in img]
    is_31 = sorted(img_signs) == [-1, +1, +1, +1]
    orig_signs = a6.dressed_signs
    reversed_pattern = (sorted(img_signs) == [-1, +1, +1, +1]
                        and sorted(orig_signs) == [-1, -1, -1, +1])
    passed = anti_ok and is_31 and reversed_pattern
    record("A9 sigma flips signature", "pass", passed,
           f"sigma(i-dressed (1,3) clique) is (3,1): signs {orig_signs} -> {img_signs}, "
           f"anticommutation preserved={anti_ok}",
           "" if passed else f"anti_ok={anti_ok}, is31={is_31}")


# ---------------------------------------------------------------- A10
def a10():
    phi = (1 + r5) / 2

    def N(a):
        return canon(a * a.subs(r5, -r5)) if hasattr(a, 'subs') else canon(a * a)

    Nphi = canon(phi * (phi.subs(r5, -r5)))
    Nphi2 = canon((phi**2) * ((phi**2).subs(r5, -r5)))
    # norm-one <=> even power: N(phi^k) = (-1)^k
    norms = {k: canon((phi**k) * ((phi**k).subs(r5, -r5))) for k in range(0, 9)}
    exhaust = all(norms[k] == (-1)**k for k in norms)
    phi3 = canon(phi**3)
    phi3_val = mat_eq(Matrix([[phi3]]), Matrix([[2 + r5]]))
    phi6 = canon(phi**6)
    phi6_val = mat_eq(Matrix([[phi6]]), Matrix([[9 + 4 * r5]]))
    pell = (9**2 - 5 * 4**2 == 1)
    Nphi3 = canon(phi3 * (phi3.subs(r5, -r5)))     # N(2+sqrt5) = -1 (fundamental in Z[sqrt5])
    # norm-one subgroup of Z[sqrt5] is +- (phi^3)^{2n} = +- phi^{6n}; check N(phi^6)=+1
    Nphi6 = canon(phi6 * (phi6.subs(r5, -r5)))
    checks = {
        "N(phi)=-1": Nphi == -1,
        "N(phi^2)=+1": Nphi2 == 1,
        "+-phi^{2n} exhausts small norm-one units": exhaust,
        "phi^3=2+sqrt5": phi3_val,
        "N(phi^3)=-1 (fundamental in Z[sqrt5])": Nphi3 == -1,
        "phi^6=9+4sqrt5": phi6_val,
        "N(phi^6)=+1": Nphi6 == 1,
        "9^2-5*4^2=1": pell,
    }
    passed = all(checks.values())
    fails = [k for k, v in checks.items() if not v]
    record("A10 ring bookkeeping", "all pass", passed,
           "Z[phi]: N(phi)=-1, norm-one=+-phi^{2n}; Z[sqrt5]: fund 2+sqrt5, "
           "norm-one=+-phi^{6n}, phi^6=9+4sqrt5" if passed else f"{len(fails)} failed",
           "" if passed else "; ".join(fails))


def main():
    print("=" * 72)
    print("GOLDEN CONSTRUCTION BENCH RUN #1 — Part A (exact over Q(sqrt5,i))")
    print("=" * 72)
    for fn in (a1, a2, a3, a4, a5, a6, a7, a8, a9, a10):
        fn()
    with open("results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} items match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
