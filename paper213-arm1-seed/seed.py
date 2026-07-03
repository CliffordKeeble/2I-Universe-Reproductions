"""
seed.py -- Paper 213 Arm 1: the Paper-191 seed and its connection-admissible
normalisation, with an exact (symbolic) verification of the hinge facts.

This is the local stand-in for bootstrap-universe's golden-dirac/golden_dirac.py
(absent from the reproductions repo). Arm 1 only ever needs the order-2 group
<I, R>, so the heavy ZZ[phi] Dirac machinery is not required -- but we still
verify the seed's algebraic properties symbolically (DERIVED, not just float).

Seed (Paper 191 sec.6):   Gamma_seed = [[0, 1], [sqrt5, 0]],  det = -sqrt5
Normalisation:            R = Gamma_seed / 5^(1/4)

HINGE FACTS (verified below, symbolically):
  (a) det(R) = -1                       -> R is a linear "reflection"
  (b) R^2   = I                         -> R is order 2 (an involution)
  (c) R^T R = diag(sqrt5, 1/sqrt5) != I -> R is NOT orthogonal.

Fact (c) is the documented wrinkle (see RESULTS.md sec."Seed normalisation").
The brief's hinge sentence calls R an "order-2 reflection" implying O(2);
(a)+(b) hold, but (c) shows the 5^(1/4) scaling fixes det and order, NOT
orthogonality. A scalar multiple of Gamma_seed can never be orthogonal because
its columns have unequal norms (sqrt5 and 1). This does NOT move the verdict:
the holonomy group <I, R> = Z/2 is abelian with no 5-torsion *because R^2 = I*,
orthogonal or not -- so it structurally cannot carry 2I representation content.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# Exact (symbolic) seed
# ---------------------------------------------------------------------------
_s5 = sp.sqrt(5)
GAMMA_SEED_SYM = sp.Matrix([[0, 1], [_s5, 0]])
R_SYM = GAMMA_SEED_SYM / sp.root(5, 4)          # 5^(1/4)
I_SYM = sp.eye(2)


def verify_hinge():
    """Symbolically verify the three hinge facts. Returns a dict of results;
    raises AssertionError if (a) or (b) fail (those are load-bearing)."""
    detR = sp.simplify(R_SYM.det())
    R2 = sp.simplify(R_SYM * R_SYM)
    RtR = sp.simplify(R_SYM.T * R_SYM)
    is_orth = (RtR == I_SYM)

    assert detR == -1, f"det(R) expected -1, got {detR}"
    assert R2 == I_SYM, f"R^2 expected I, got {R2}"
    return {
        "det_R": detR,                      # -1            (a)
        "R_squared": R2,                    # I             (b)
        "RT_R": RtR,                        # diag(sqrt5,1/sqrt5)
        "R_is_orthogonal": bool(is_orth),   # False         (c)
        "det_Gamma_seed": sp.simplify(GAMMA_SEED_SYM.det()),   # -sqrt5
    }


# ---------------------------------------------------------------------------
# Float matrices used by the spectral code
# ---------------------------------------------------------------------------
I2 = np.eye(2)


def reflection_from_offdiag(a):
    """The Arm-1 family of order-2 involutions: M = [[0, 1], [a, 0]] / sqrt(a),
    which has det = -1 and M^2 = I for any a > 0. a = 5 gives the golden seed R;
    a in {2, 3} are the non-golden controls (sec.5.1). For a == 1, M is the
    genuine orthogonal reflection [[0,1],[1,0]] (the only orthogonal member)."""
    a = float(a)
    M = np.array([[0.0, 1.0], [a, 0.0]]) / np.sqrt(a)
    return M


# The golden seed normalisation, as floats.
GAMMA_SEED = np.array([[0.0, 1.0], [np.sqrt(5.0), 0.0]])
R = GAMMA_SEED / (5.0 ** 0.25)          # == reflection_from_offdiag(5)

# The two group elements as named O(2)-words. (Stored as 2x2 float blocks.)
ELEMENTS = {"I": I2, "R": R}


def word_product(elems):
    """Product of a sequence of element labels (e.g. ['R','R','I']) -> matrix."""
    M = I2.copy()
    for e in elems:
        M = M @ ELEMENTS[e]
    return M


def is_flat(elems, tol=1e-9):
    """A cycle's holonomy is flat iff the product of its edge elements == I.
    For words in {I, R} with R^2 = I this is exactly 'even number of R's'."""
    return bool(np.allclose(word_product(elems), I2, atol=tol))


if __name__ == "__main__":
    print("Paper 213 Arm 1 -- seed verification (symbolic)\n")
    res = verify_hinge()
    print(f"  Gamma_seed = [[0,1],[sqrt5,0]],  det = {res['det_Gamma_seed']}")
    print(f"  R = Gamma_seed / 5^(1/4)")
    print(f"  (a) det(R)        = {res['det_R']}        [expected -1]")
    print(f"  (b) R^2           = {list(res['R_squared'])}   [expected I]")
    print(f"  (c) R^T R         = {list(res['RT_R'])}")
    print(f"      R orthogonal? = {res['R_is_orthogonal']}   "
          f"[FALSE -> documented wrinkle]")
    print(f"\n  Combinatorial flatness (R^2 = I): "
          f"even #R -> flat, odd #R -> not flat")
    print(f"    word [R,R] flat? {is_flat(['R','R'])}   "
          f"word [R] flat? {is_flat(['R'])}")
    print(f"\n  Group <I,R> = Z/2 (abelian, no 5-torsion) -> "
          f"cannot carry 2I content. Verdict unaffected by (c).")
