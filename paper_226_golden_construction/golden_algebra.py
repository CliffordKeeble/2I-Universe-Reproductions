"""
golden_algebra.py — shared exact-arithmetic substrate for the Golden Construction
Bench Run #1 (Paper 226 target; Fizz brief, 5 Jul 2026).

All definitions are transcribed from the source papers and confirmed against them:
  - Gamma_seed, Gamma_adj, K2, J2, Z2 : Paper 203 v0.3 section 4 (quoting Paper 191 v1.2).
  - The 16 sqrt(5)-normalised candidates and the six (2,2) cliques : Paper 204 v0.4
    section 5 / Appendix A (A tensor B, A,B in {I,K,J,Z,Gamma_seed,Gamma_adj}, squares +-sqrt5.I4).
  - Galois actions tau (i->-i, sqrt5 fixed) and sigma (sqrt5->-sqrt5, i fixed) : brief part A setup.

Exact arithmetic over Q(sqrt5, i) throughout via sympy sqrt(5) and I. No floating point.

Provenance note: the paper .md drafts live in the user's Downloads and are NOT committed
to this repo (they are unreleased). This module vendors only the machine-checkable
definitions, each tagged with its source location above.
"""

import sympy as sp
from sympy import sqrt, I, Rational, eye, zeros, Matrix, simplify, expand

r5 = sqrt(5)

# ---------------------------------------------------------------------------
# 2x2 building blocks  (Paper 203 v0.3 section 4)
# ---------------------------------------------------------------------------
I2 = eye(2)
K2 = Matrix([[0, 1], [1, 0]])        # K^2 = +I
J2 = Matrix([[0, -1], [1, 0]])       # J^2 = -I
Z2 = Matrix([[1, 0], [0, -1]])       # Z^2 = +I,  Z = 1/2 [K,J]
G_seed = Matrix([[0, 1], [r5, 0]])   # Gamma_seed^2 = +sqrt5 . I
G_adj = Matrix([[0, -1], [r5, 0]])   # Gamma_adj^2  = -sqrt5 . I

BLOCKS_4 = {"I": I2, "K": K2, "J": J2, "Z": Z2}          # square +-I
BLOCKS_G = {"Gs": G_seed, "Ga": G_adj}                   # square +-sqrt5.I


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def kron(A, B):
    """Kronecker (tensor) product of two sympy matrices, A (x) B."""
    ra, ca = A.shape
    rb, cb = B.shape
    out = zeros(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            out[i * rb:(i + 1) * rb, j * cb:(j + 1) * cb] = A[i, j] * B
    return out


def canon(expr):
    """Canonicalise a scalar sympy expression over Q(sqrt5,i)."""
    return sp.expand(sp.simplify(expr))


def mcanon(M):
    return M.applyfunc(canon)


def mat_eq(A, B):
    """Exact equality of two matrices over Q(sqrt5,i)."""
    D = mcanon(A - B)
    return all(D[i] == 0 for i in range(len(D)))


def comm(A, B):
    return A * B - B * A


def anticomm(A, B):
    return A * B + B * A


def is_zero(M):
    return mat_eq(M, zeros(*M.shape))


# ---------------------------------------------------------------------------
# Galois actions  (entrywise)
#   tau : i -> -i,  sqrt5 fixed   (complex conjugation of the entries)
#   sigma : sqrt5 -> -sqrt5,  i fixed
#   sigma.tau : both
# ---------------------------------------------------------------------------
def tau(M):
    # i -> -i with sqrt5 real and fixed: elementwise complex conjugate.
    return M.applyfunc(lambda e: canon(sp.conjugate(e)))


def sigma(M):
    return M.applyfunc(lambda e: canon(e.subs(r5, -r5)))


def sigmatau(M):
    return sigma(tau(M))


# ---------------------------------------------------------------------------
# square-sign of a golden gamma:  +1 if G^2 = +sqrt5.I, -1 if -sqrt5.I, else None
# ---------------------------------------------------------------------------
def square_sign(G):
    n = G.shape[0]
    Sq = mcanon(G * G)
    for s in (+1, -1):
        if mat_eq(Sq, s * r5 * eye(n)):
            return s
    # also allow +-25 etc for volume elements -> return the scalar instead
    return None


def square_scalar(G):
    """Return c such that G^2 = c.I  (c a sympy scalar), or None if not scalar."""
    n = G.shape[0]
    Sq = mcanon(G * G)
    c = Sq[0, 0]
    if mat_eq(Sq, c * eye(n)):
        return canon(c)
    return None


# ---------------------------------------------------------------------------
# The 16 sqrt(5)-normalised 4x4 candidates  (Paper 204 v0.4 section 5)
#   A (x) B with exactly one factor a golden block -> square = +- sqrt5 . I4
# ---------------------------------------------------------------------------
def build_candidates():
    cands = []  # list of (label, matrix, sqsign)
    # {I,K,J,Z} (x) {Gs,Ga}
    for an, A in BLOCKS_4.items():
        for gn, G in BLOCKS_G.items():
            M = kron(A, G)
            cands.append((f"{an}(x){gn}", M, square_sign(M)))
    # {Gs,Ga} (x) {I,K,J,Z}
    for gn, G in BLOCKS_G.items():
        for an, A in BLOCKS_4.items():
            M = kron(G, A)
            cands.append((f"{gn}(x){an}", M, square_sign(M)))
    return cands


def anticommute(A, B):
    return is_zero(anticomm(A, B))


def find_cliques(cands, size=4):
    """All size-cliques of pairwise-anticommuting candidates (by index, sorted)."""
    n = len(cands)
    mats = [c[1] for c in cands]
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if anticommute(mats[i], mats[j]):
                adj[i][j] = adj[j][i] = True
    cliques = []

    def extend(current, candidates_idx):
        if len(current) == size:
            cliques.append(tuple(current))
            return
        for k in list(candidates_idx):
            if all(adj[k][m] for m in current):
                nxt = {x for x in candidates_idx if x > k and adj[k][x]}
                extend(current + [k], nxt)

    extend([], set(range(n)))
    # dedup (extend already enforces increasing index -> unique)
    return cliques


# ---------------------------------------------------------------------------
# self-test  (run:  python golden_algebra.py)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ok = True

    def check(name, cond):
        global ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    print("== block squares ==")
    check("K^2 = +I", mat_eq(K2 * K2, I2))
    check("J^2 = -I", mat_eq(J2 * J2, -I2))
    check("Z^2 = +I", mat_eq(Z2 * Z2, I2))
    check("Z = 1/2 [K,J]", mat_eq(Z2, Rational(1, 2) * comm(K2, J2)))
    check("Gseed^2 = +sqrt5 I", mat_eq(G_seed * G_seed, r5 * I2))
    check("Gadj^2  = -sqrt5 I", mat_eq(G_adj * G_adj, -r5 * I2))

    print("== Galois maps validate ==")
    check("sigma(sqrt5) = -sqrt5", canon((r5).subs(r5, -r5) + r5) == 0)
    check("sigma(Gseed) = -Gadj", mat_eq(sigma(G_seed), -G_adj))
    check("sigma(Gadj)  = -Gseed", mat_eq(sigma(G_adj), -G_seed))
    check("sigma involutive on Gseed", mat_eq(sigma(sigma(G_seed)), G_seed))
    check("tau fixes Gseed (real)", mat_eq(tau(G_seed), G_seed))
    check("tau(i.Gseed) = -i.Gseed", mat_eq(tau(I * G_seed), -I * G_seed))
    check("tau involutive", mat_eq(tau(tau(I * G_adj)), I * G_adj))
    check("sigma fixes i", canon(sp.conjugate(sigma(Matrix([[I]]))[0]) + I) == 0)

    print("== 16 candidates + cliques (Paper 204) ==")
    cands = build_candidates()
    check("exactly 16 candidates", len(cands) == 16)
    check("all have square +-sqrt5 I", all(c[2] in (+1, -1) for c in cands))
    cliques = find_cliques(cands, 4)
    check("exactly 6 four-cliques", len(cliques) == 6)
    sigs = []
    for cl in cliques:
        signs = sorted(cands[k][2] for k in cl)
        sigs.append(signs)
    check("all cliques are (2,2) [two +, two -]",
          all(s == [-1, -1, 1, 1] for s in sigs))

    print("\nRESULT:", "ALL SELF-TESTS PASS" if ok else "SELF-TEST FAILURE")
    import sys
    sys.exit(0 if ok else 1)
