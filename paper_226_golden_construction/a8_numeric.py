"""
a8_numeric.py — trustworthy descent invariant for the i-dressed (1,3) clique,
computed at high precision (mpmath, 50 digits) to bypass the symbolic-simplify
inconsistency seen in a8_probe. B is solved EXACTLY (linear algebra, reliable),
then J^2 = B.conj(B) is evaluated numerically at BOTH real embeddings of K.
"""
import mpmath as mp
import sympy as sp
from sympy import I, eye, zeros, Matrix, symbols
from golden_algebra import (r5, G_seed, G_adj, kron, mcanon, mat_eq, is_zero, tau,
                            square_scalar, build_candidates, find_cliques)

mp.mp.dps = 50


def solve_B(gammas):
    """Exact B with B.tau(G) = G.B for all gammas (tau = entrywise conj)."""
    b = symbols('b0:16')
    B = Matrix(4, 4, b)
    eqs = []
    for G in gammas:
        eqs += list(B * tau(G) - G * B)      # exact, no simplify
    sol = list(sp.linsolve(eqs, list(b)))[0]
    free = set()
    for comp in sol:
        free |= (comp.free_symbols & set(b))
    subs0 = {p: 0 for p in b}
    if free:
        subs0[sorted(free, key=lambda s: s.name)[0]] = 1
    return Matrix(4, 4, [c.subs(subs0) for c in sol]), len(free)


def to_mp(M, r5val):
    """sympy 4x4 -> mpmath complex matrix with sqrt5 -> r5val (mpf, +/-)."""
    out = mp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            e = M[i, j].subs(r5, r5val)
            re = mp.mpf(str(sp.N(sp.re(e), 40)))
            im = mp.mpf(str(sp.N(sp.im(e), 40)))
            out[i, j] = mp.mpc(re, im)
    return out


def cmatnorm(A):
    return max(abs(A[i, j]) for i in range(A.rows) for j in range(A.cols))


def analyse(name, gammas):
    B, dim = solve_B(gammas)
    r5p, r5m = mp.sqrt(5), -mp.sqrt(5)
    report = []
    for lab, r5val in (("+sqrt5", r5p), ("-sqrt5", r5m)):
        Bn = to_mp(B, r5val)
        Gn = [to_mp(G, r5val) for G in gammas]
        # verify intertwining B conj(G) - G B = 0  (tau(G) = conj(G))
        inter = max(cmatnorm(Bn * Gn[k].conjugate() - Gn[k] * Bn) for k in range(4))
        # J^2 = B conj(B)
        Jsq = Bn * Bn.conjugate()
        cI = Jsq[0, 0]
        scalar_err = cmatnorm(Jsq - cI * mp.eye(4))
        report.append((lab, inter, cI, scalar_err))
    print(f"{name}: dim(B)={dim}")
    for lab, inter, cI, serr in report:
        print(f"    [{lab}] intertwine_resid={mp.nstr(inter,3)} "
              f"J^2 = ({mp.nstr(cI,8)}).I  scalar_resid={mp.nstr(serr,3)}")
    return report


def main():
    cands = build_candidates()
    cliques = find_cliques(cands, 4)
    cl0 = cliques[0]
    mats0 = [cands[i][1] for i in cl0]
    signs0 = [cands[i][2] for i in cl0]
    pos = [i for i, s in enumerate(signs0) if s == +1]

    print("KNOWN — undressed (2,2) all-real clique (expect clean J^2=+1):")
    analyse("(2,2) clique0", mats0)

    print("\ni-DRESSED (1,3) clique (dress first +sqrt5 gamma by i):")
    dressed = list(mats0); dressed[pos[0]] = mcanon(I * mats0[pos[0]])
    analyse("(1,3) dressed", dressed)

    print("\nAll six cliques, i-dressed first +sqrt5 gamma:")
    for ci in range(len(cliques)):
        m = [cands[i][1] for i in cliques[ci]]
        s = [cands[i][2] for i in cliques[ci]]
        p0 = [i for i, x in enumerate(s) if x == +1][0]
        d = list(m); d[p0] = mcanon(I * m[p0])
        analyse(f"clique{ci}", d)


if __name__ == "__main__":
    main()
