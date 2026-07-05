"""
a8_probe.py — is the A8 J^2=+1 result real or a bug?
Verify the intertwiner machinery on a KNOWN case (all-real (2,2) clique -> must give
J^2=+1, split), then test the i-dressed (1,3) clique under several dressing choices.
Pure diagnostic; does not write repo outputs.
"""
import sympy as sp
from sympy import I, eye, zeros, Matrix, symbols
from golden_algebra import (r5, I2, K2, J2, Z2, G_seed, G_adj, kron, mcanon, mat_eq,
                            is_zero, tau, sigma, square_sign, square_scalar,
                            build_candidates, find_cliques, anticommute)

I4 = eye(4)


def find_B(gammas):
    """Antilinear J=B o tau commutes with all gammas: B tau(G)=G B.
       Returns (dim, B0) for a concrete nonzero B0."""
    b = symbols('b0:16')
    B = Matrix(4, 4, b)
    eqs = []
    for G in gammas:
        tG = tau(G)
        if mat_eq(tG, G):          # tau-fixed -> commute
            eqs += list(mcanon(B * G - G * B))
        elif mat_eq(tG, -G):       # tau-negated -> anticommute
            eqs += list(mcanon(B * G + G * B))
        else:                      # general: B tau(G) = G B
            eqs += list(mcanon(B * tG - G * B))
    sol = list(sp.linsolve(eqs, list(b)))
    free = set()
    for comp in sol[0]:
        free |= (comp.free_symbols & set(b))
    dim = len(free)
    subs0 = {p: 0 for p in b}
    if free:
        subs0[sorted(free, key=lambda s: s.name)[0]] = 1
    B0 = mcanon(Matrix(4, 4, [c.subs(subs0) for c in sol[0]]))
    return dim, B0


def Japply(B, v):
    return mcanon(B * tau(v))


def analyse(name, gammas):
    dim, B0 = find_B(gammas)
    if is_zero(B0):
        print(f"{name}: B space dim={dim}, B0=ZERO -> no J"); return
    # verify B0 solves constraints
    ok_constraint = True
    for G in gammas:
        tG = tau(G)
        if not is_zero(mcanon(B0 * tG - G * B0)):
            ok_constraint = False
    # verify J commutes with each gamma as an operator on the 4 basis vectors
    ok_commute = True
    for G in gammas:
        for e in range(4):
            v = zeros(4, 1); v[e] = 1
            lhs = Japply(B0, mcanon(G * v))      # J(G v)
            rhs = mcanon(G * Japply(B0, v))       # G (J v)
            if not is_zero(mcanon(lhs - rhs)):
                ok_commute = False
    # J^2 two ways
    Jsq_alg = mcanon(B0 * tau(B0))
    c = square_scalar(Jsq_alg)
    # by acting on basis
    ok_Jsq = True
    for e in range(4):
        v = zeros(4, 1); v[e] = 1
        j2v = Japply(B0, Japply(B0, v))
        if c is not None and not is_zero(mcanon(j2v - c * v)):
            ok_Jsq = False
    signs = [square_sign(G) for G in gammas]
    cplus = sp.N(c.subs(r5, sp.sqrt(5))) if c is not None and sp.im(c) == 0 else None
    cminus = sp.N(c.subs(r5, -sp.sqrt(5))) if c is not None and sp.im(c) == 0 else None
    verdict = ("J^2=+1 (split/real)" if (cplus and cplus > 0) else
               "J^2=-1 (quaternionic)" if (cplus and cplus < 0) else "??")
    print(f"{name}: signs={signs} dim(B)={dim} constraint_ok={ok_constraint} "
          f"commute_ok={ok_commute} Jsq_consistent={ok_Jsq}")
    print(f"    J^2 = ({c})*I4 ;  at real places (+sqrt5, -sqrt5) = ({cplus}, {cminus})"
          f"  =>  {verdict}")


def main():
    cands = build_candidates()
    cliques = find_cliques(cands, 4)
    cl0 = cliques[0]
    mats0 = [cands[i][1] for i in cl0]
    signs0 = [cands[i][2] for i in cl0]
    pos = [i for i, s in enumerate(signs0) if s == +1]
    neg = [i for i, s in enumerate(signs0) if s == -1]

    print("KNOWN CASE — undressed (2,2) all-real clique (expect J^2=+1, split):")
    analyse("(2,2) clique0", mats0)

    print("\ni-DRESSED (1,3) cliques (dress one +sqrt5 gamma by i):")
    for p in pos:
        for dsign, dlab in ((I, "+i"), (-I, "-i")):
            dressed = list(mats0)
            dressed[p] = mcanon(dsign * mats0[p])
            analyse(f"(1,3) dress gamma#{p} by {dlab}", dressed)

    print("\nALT — dress a -sqrt5 gamma by i instead (gives (3,1)-type):")
    dressed = list(mats0)
    dressed[neg[0]] = mcanon(I * mats0[neg[0]])
    analyse(f"dress -sqrt5 gamma#{neg[0]} by +i", dressed)

    print("\nsigma-IMAGE of the standard i-dressed (1,3) clique (should be (3,1)):")
    dressed = list(mats0); dressed[pos[0]] = mcanon(I * mats0[pos[0]])
    analyse("sigma((1,3))", [mcanon(sigma(M)) for M in dressed])

    print("\nOTHER CLIQUES (i-dressed first +sqrt5 gamma):")
    for ci in range(1, len(cliques)):
        m = [cands[i][1] for i in cliques[ci]]
        s = [cands[i][2] for i in cliques[ci]]
        p0 = [i for i, x in enumerate(s) if x == +1][0]
        d = list(m); d[p0] = mcanon(I * m[p0])
        analyse(f"(1,3) clique{ci}", d)


if __name__ == "__main__":
    main()
