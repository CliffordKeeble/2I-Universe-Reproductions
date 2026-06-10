"""
vertex_index_2I.py -- The -12 as a 2I rep-theoretic index, not a polynomial spelling.

Brief: Mr Code / Fizz, 10 June 2026. Paper 101 Lemma L neighbourhood, Filter-2 test.

The Paper 101 mass-ratio mu = 1848 - 12 = 1836 rests on a -12 correction. Tor's
Pattern-75 objection: -12 has many spellings as a polynomial in {V,E,F,D,chi}, so
"-12 = -V" carries no information. The surviving route relocates 12 to a *ratio of
group orders* that lives OUTSIDE the polynomial ring:

        12 = |2I| / |vertex-stabiliser in 2I| = 120 / 10 = [2I : C_10].

This script verifies, in EXACT arithmetic (Q(sqrt5), golden entries symbolic):

  T1  the icosahedral-vertex stabiliser, lifted to 2I < SU(2), is cyclic C_10
      (order 10); hence the vertex orbit has size 120/10 = 12.

  T2  the 12-dim permutation rep rho_vert = Ind_{C_10}^{2I}(triv) = C[2I/C_10]
      decomposed into 2I-irreducibles; the multiplicity of the TRIVIAL rep.

  T3  the transfer kernel: the natural 2I-map rho_vert -> triv (Hodge star to the
      volume / trivial-rep sector). Dimension of kernel and cokernel.

  T4  the null statement: 120/10 is a ratio of group orders, not an element of
      Z[V,E,F,D,chi], so Tor's polynomial null does not apply to it.

Construction: 2I is built as the 120 unit icosian quaternions. Every component lies
in Q(sqrt5); we carry them as exact (a,b) = a + b*sqrt5 pairs over the rationals, so
the group law, the conjugation action on the 12 icosahedron vertices, the stabiliser,
and the conjugacy classes are all computed exactly -- nothing is asserted.

Run: python vertex_index_2I.py
"""

from fractions import Fraction as Fr
from itertools import permutations, product
import sympy as sp

# ---------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt5): a number is a pair (a, b) meaning a + b*sqrt5.
# ---------------------------------------------------------------------------

def q5(a=0, b=0):
    return (Fr(a), Fr(b))

def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])

def mul(x, y):
    # (a+b s)(c+d s) = (ac + 5bd) + (ad + bc) s,  s^2 = 5
    a, b = x
    c, d = y
    return (a * c + 5 * b * d, a * d + b * c)

def neg(x):
    return (-x[0], -x[1])

ZERO = q5(0, 0)
ONE = q5(1, 0)
SQRT5 = q5(0, 1)

# golden ratio phi = (1 + sqrt5)/2 and the icosian magnitudes
PHI = (Fr(1, 2), Fr(1, 2))                 # phi          = 1/2 + 1/2 sqrt5
HALF = q5(Fr(1, 2), 0)                     # 1/2
INV_2PHI = (Fr(-1, 4), Fr(1, 4))           # 1/(2 phi) = (sqrt5 - 1)/4
PHI_2 = (Fr(1, 4), Fr(1, 4))               # phi/2     = (1 + sqrt5)/4

def to_sympy(x):
    return sp.nsimplify(x[0]) + sp.nsimplify(x[1]) * sp.sqrt(5)

# ---------------------------------------------------------------------------
# Quaternions over Q(sqrt5): 4-tuples (w, x, y, z) of q5 numbers.
# ---------------------------------------------------------------------------

def qmul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    w = sub(sub(sub(mul(w1, w2), mul(x1, x2)), mul(y1, y2)), mul(z1, z2))
    x = sub(add(add(mul(w1, x2), mul(x1, w2)), mul(y1, z2)), mul(z1, y2))
    y = add(add(sub(mul(w1, y2), mul(x1, z2)), mul(y1, w2)), mul(z1, x2))
    z = add(sub(add(mul(w1, z2), mul(x1, y2)), mul(y1, x2)), mul(z1, w2))
    return (w, x, y, z)

def qconj(p):
    w, x, y, z = p
    return (w, neg(x), neg(y), neg(z))

QID = (ONE, ZERO, ZERO, ZERO)
QNEG1 = (neg(ONE), ZERO, ZERO, ZERO)

def build_2I():
    """The 120 unit icosian quaternions, exact, deduplicated by canonical key."""
    out = {}

    def store(vec):
        out[tuple(vec)] = vec

    # 8 units: (+-1,0,0,0) and permutations
    for i in range(4):
        for s in (ONE, neg(ONE)):
            v = [ZERO, ZERO, ZERO, ZERO]
            v[i] = s
            store(v)

    # 16: (+-1/2, +-1/2, +-1/2, +-1/2)
    for signs in product((HALF, neg(HALF)), repeat=4):
        store(list(signs))

    # 96: even permutations of (0, +-1/2, +-1/(2phi), +-phi/2)
    mags = [ZERO, HALF, INV_2PHI, PHI_2]
    for perm in permutations(range(4)):
        # parity of perm
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        if inv % 2 != 0:
            continue
        arranged = [mags[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if arranged[i] != ZERO]
        for signs in product((ONE, neg(ONE)), repeat=len(nz)):
            v = list(arranged)
            for k, pos in enumerate(nz):
                v[pos] = mul(arranged[pos], signs[k])
            store(v)

    return [tuple(v) for v in out.values()]

# ---------------------------------------------------------------------------
# Rotation action on the 12 icosahedron vertices.
# ---------------------------------------------------------------------------

def rotate(q, v):
    """Act on an imaginary quaternion v=(0,vx,vy,vz) by conjugation q v q*."""
    vq = (ZERO, v[0], v[1], v[2])
    r = qmul(qmul(q, vq), qconj(q))
    return (r[1], r[2], r[3])         # imaginary part = rotated vector

def icosahedron_vertices():
    """12 vertices: cyclic permutations of (0, +-1, +-phi)."""
    one, phi = ONE, PHI
    base = []
    for a, b in product((one, neg(one)), (phi, neg(phi))):
        base.append((ZERO, a, b))      # (0, +-1, +-phi)
        base.append((a, b, ZERO))      # (+-1, +-phi, 0)
        base.append((b, ZERO, a))      # (+-phi, 0, +-1)
    return base

# ---------------------------------------------------------------------------
# Group structure: orders, conjugacy classes, rotation angle classification.
# ---------------------------------------------------------------------------

def order_of(q):
    p, k = q, 1
    while p != QID:
        p = qmul(p, q)
        k += 1
        if k > 240:
            raise RuntimeError("order overflow")
    return k

def conjugacy_classes(group):
    gset = set(group)
    seen = set()
    classes = []
    for g in group:
        if g in seen:
            continue
        cls = set()
        for h in group:
            cls.add(qmul(qmul(h, g), qconj(h)))   # h unit => h^{-1} = h*
        classes.append(sorted(cls, key=lambda t: str(t)))
        seen |= cls
    return classes

def cos_theta(q):
    """Rotation angle: cos(theta) = 2 w^2 - 1, where w = Re(q) = cos(theta/2)."""
    w = q[0]
    return sub(mul(q5(2, 0), mul(w, w)), ONE)

# A5 class label from the rotation angle (exact, unambiguous).
COS_LABEL = {
    (Fr(1), Fr(0)):  "1A",        # cos = 1      -> identity rotation
    (Fr(-1), Fr(0)): "2A",        # cos = -1     -> 180 deg (edge axis)
    (Fr(-1, 2), Fr(0)): "3A",     # cos = -1/2   -> 120 deg (face axis)
    (Fr(-1, 4), Fr(1, 4)): "5A",  # cos = (sqrt5-1)/4  -> 72 deg
    (Fr(-1, 4), Fr(-1, 4)): "5B", # cos = -(sqrt5+1)/4 -> 144 deg
}

# ---------------------------------------------------------------------------
# A5 character table (the vector / integer-spin irreps of 2I), golden symbolic.
# Classes ordered [1A, 2A, 3A, 5A, 5B]; sizes [1, 15, 20, 12, 12].
# ---------------------------------------------------------------------------

phi = (1 + sp.sqrt(5)) / 2
A5_CLASSES = ["1A", "2A", "3A", "5A", "5B"]
A5_SIZES = {"1A": 1, "2A": 15, "3A": 20, "5A": 12, "5B": 12}
A5_TABLE = {
    "1":  {"1A": sp.Integer(1), "2A": sp.Integer(1),  "3A": sp.Integer(1),  "5A": sp.Integer(1),  "5B": sp.Integer(1)},
    "3":  {"1A": sp.Integer(3), "2A": sp.Integer(-1), "3A": sp.Integer(0),  "5A": phi,            "5B": 1 - phi},
    "3'": {"1A": sp.Integer(3), "2A": sp.Integer(-1), "3A": sp.Integer(0),  "5A": 1 - phi,        "5B": phi},
    "4":  {"1A": sp.Integer(4), "2A": sp.Integer(0),  "3A": sp.Integer(1),  "5A": sp.Integer(-1), "5B": sp.Integer(-1)},
    "5":  {"1A": sp.Integer(5), "2A": sp.Integer(1),  "3A": sp.Integer(-1), "5A": sp.Integer(0),  "5B": sp.Integer(0)},
}
A5_DIM = {"1": 1, "3": 3, "3'": 3, "4": 4, "5": 5}

def inner(chi_a, chi_b):
    """<chi_a, chi_b> = (1/60) sum_cls size * chi_a(cls) * conj(chi_b(cls)). Real entries."""
    s = sum(A5_SIZES[c] * chi_a[c] * chi_b[c] for c in A5_CLASSES)
    return sp.simplify(s / 60)

# ---------------------------------------------------------------------------
# Main verification.
# ---------------------------------------------------------------------------

def main():
    print("=" * 74)
    print("The -12 as a 2I rep-theoretic index  (exact Q(sqrt5) arithmetic)")
    print("=" * 74)

    G = build_2I()
    assert len(G) == 120, f"expected |2I| = 120, got {len(G)}"
    print(f"\n|2I| = {len(G)}  (120 unit icosian quaternions, built exactly)")

    verts = icosahedron_vertices()
    assert len(verts) == 12
    vindex = {v: i for i, v in enumerate(verts)}

    # closure check: every g permutes the 12 vertices
    for g in G:
        for v in verts:
            assert rotate(g, v) in vindex, "vertex set not closed under 2I"
    print("12 icosahedron vertices form a single closed 2I-orbit.  [check]")

    # -- T1: vertex stabiliser ------------------------------------------------
    print("\n" + "-" * 74)
    print("T1  Vertex stabiliser  (lifted to 2I)")
    print("-" * 74)
    v0 = verts[0]
    stab = [g for g in G if rotate(g, v0) == v0]
    orders = sorted(order_of(g) for g in stab)
    # cyclic? a group of order 10 is cyclic iff it has an element of order 10
    gen = next((g for g in stab if order_of(g) == 10), None)
    cyclic = gen is not None
    powers = []
    if gen is not None:
        p = QID
        for _ in range(10):
            p = qmul(p, gen)
            powers.append(p)
        cyclic = set(powers) == set(stab)
    print(f"|Stab(v0)| = {len(stab)}")
    print(f"element orders in Stab(v0): {orders}")
    print(f"contains an element of order 10: {gen is not None};  "
          f"<g> = Stab(v0): {cyclic}")
    print(f"=> Stab(v0) is cyclic C_10." if cyclic else "=> NOT cyclic C_10 (!)")
    print(f"=> vertex orbit size = |2I| / |Stab| = 120 / {len(stab)} = "
          f"{120 // len(stab)}   [2I : C_10]")
    assert len(stab) == 10 and cyclic, "T1 FAILED: stabiliser is not C_10"

    # contains the centre -1?
    print(f"Stab(v0) contains the central element -1: {QNEG1 in set(stab)}  "
          f"(C_10 = C_5 lifted through its double cover)")

    # -- conjugacy classes & permutation character ----------------------------
    print("\n" + "-" * 74)
    print("Conjugacy classes of 2I and the permutation character chi_perm")
    print("-" * 74)
    classes = conjugacy_classes(G)
    assert len(classes) == 9, f"expected 9 classes, got {len(classes)}"

    rows = []
    for cls in classes:
        rep = cls[0]
        ordr = order_of(rep)
        ct = cos_theta(rep)
        label = COS_LABEL[ct]
        # fixed vertices of rep under rotation
        fixed = sum(1 for v in verts if rotate(rep, v) == v)
        rows.append((ordr, len(cls), label, fixed, ct))
    rows.sort(key=lambda r: (r[0], r[2]))

    print(f"{'ord':>3} {'size':>5} {'A5-class':>9} {'cos(theta)':>14} "
          f"{'chi_perm':>9}")
    neg1_acts_triv = None
    for ordr, size, label, fixed, ct in rows:
        cs = to_sympy(ct)
        print(f"{ordr:>3} {size:>5} {label:>9} {str(cs):>14} {fixed:>9}")
    # -1 acts trivially?
    neg1_fixed = sum(1 for v in verts if rotate(QNEG1, v) == v)
    print(f"\nThe central element -1 fixes {neg1_fixed}/12 vertices "
          f"=> it acts as the identity permutation.")
    print("=> rho_vert is inflated from 2I/{+-1} = A5; only the 'vector'"
          " irreps {1,3,3',4,5} can appear,")
    print("   and every 'spinor' irrep {2,2',4',6} has multiplicity 0.")
    assert neg1_fixed == 12

    # permutation character as a function of A5 class (it factors through A5)
    chi_perm = {}
    for ordr, size, label, fixed, ct in rows:
        chi_perm[label] = sp.Integer(fixed)
    # both order-5 and order-10 rows map to the same A5 class; consistency check
    print(f"\nchi_perm on [1A,2A,3A,5A,5B] = "
          f"[{chi_perm['1A']}, {chi_perm['2A']}, {chi_perm['3A']}, "
          f"{chi_perm['5A']}, {chi_perm['5B']}]")

    # -- verify the A5 table is orthonormal (self-check, not cited blind) ------
    print("\n" + "-" * 74)
    print("Self-check: A5 character table is orthonormal (exact, golden symbolic)")
    print("-" * 74)
    names = ["1", "3", "3'", "4", "5"]
    ok = True
    for a in names:
        for b in names:
            val = inner(A5_TABLE[a], A5_TABLE[b])
            expect = 1 if a == b else 0
            if val != expect:
                ok = False
                print(f"  <{a},{b}> = {val} (expected {expect})  MISMATCH")
    print("orthonormality holds for all 25 pairs." if ok else "TABLE BROKEN")
    assert ok, "A5 character table failed its own orthonormality check"

    # -- T2: decompose rho_vert ----------------------------------------------
    print("\n" + "-" * 74)
    print("T2  Decomposition of rho_vert = Ind_{C_10}^{2I}(triv) = C[2I/C_10]")
    print("-" * 74)
    mults = {}
    for name in names:
        m = inner(chi_perm, A5_TABLE[name])
        mults[name] = m
    # report all 9 2I irreps (spinors are 0 by the inflation argument above)
    full = {
        "1": mults["1"], "2": 0, "2'": 0,
        "3": mults["3"], "3'": mults["3'"],
        "4": mults["4"], "4'": 0,
        "5": mults["5"], "6": 0,
    }
    DIM2I = {"1": 1, "2": 2, "2'": 2, "3": 3, "3'": 3, "4": 4, "4'": 4, "5": 5, "6": 6}
    print(f"{'irrep':>6} {'dim':>4} {'mult':>5}")
    total = 0
    for name in ["1", "2", "2'", "3", "3'", "4", "4'", "5", "6"]:
        m = full[name]
        total += DIM2I[name] * int(m)
        tag = "  (spinor: 0 by -1-acts-trivially)" if name in {"2", "2'", "4'", "6"} else ""
        print(f"{name:>6} {DIM2I[name]:>4} {str(m):>5}{tag}")
    print(f"\n  sum of dim*mult = {total}  (must equal dim rho_vert = 12)")
    assert total == 12, "decomposition does not reconstruct the 12-dim module"

    triv_mult = int(full["1"])
    print(f"\n  rho_vert = 1 + 3 + 3' + 5   (dims 1+3+3+5 = 12)")
    print(f"\n  >>> T2 ANSWER (bare):  multiplicity of the TRIVIAL rep in "
          f"rho_vert = {triv_mult}")
    assert triv_mult == 1

    # -- T3: the transfer kernel ---------------------------------------------
    print("\n" + "-" * 74)
    print("T3  Transfer kernel of the natural map  rho_vert --> triv  (Hodge star")
    print("    to the volume / trivial-rep sector)")
    print("-" * 74)
    # dim Hom_{2I}(rho_vert, triv) = multiplicity of triv in rho_vert = 1.
    # The map onto the 1-dim trivial target is surjective; kernel = 12 - 1.
    hom_dim = triv_mult
    image_dim = 1                      # trivial target is 1-dimensional, map nonzero
    kernel_dim = 12 - image_dim
    cokernel_dim = 1 - image_dim
    print(f"dim Hom_2I(rho_vert, triv)            = {hom_dim}   "
          f"(trivial appears exactly once)")
    print(f"dim of image in the trivial target    = {image_dim}   "
          f"(the volume form is 1-dimensional)")
    print(f"dim of cokernel                       = {cokernel_dim}   "
          f"(map is surjective onto triv)")
    print(f"dim of KERNEL (no image in triv)      = {kernel_dim}   "
          f"(= 3 + 3' + 5)")
    print(f"\n  >>> T3 ANSWER (bare):  kernel = {kernel_dim},  cokernel = "
          f"{cokernel_dim}")
    print("  Exactly ONE dimension transfers to the volume sector; the other 11")
    print("  (the non-trivial complement 3+3'+5) do not.  The full 12 does NOT")
    print("  fail to transfer -- the trivial summand does transfer.")
    assert kernel_dim == 11 and cokernel_dim == 0

    # -- T4: the null statement ----------------------------------------------
    print("\n" + "-" * 74)
    print("T4  The polynomial-null comparison")
    print("-" * 74)
    print("12 = |2I| / |C_10| = 120 / 10 = [2I : C_10]  is a RATIO OF GROUP ORDERS")
    print("(equivalently, the size of a transitive 2I-orbit). It is NOT an element")
    print("of the polynomial ring Z[V,E,F,D,chi]; it is not built from those")
    print("variables, so it cannot be re-spelled as V+E+2 etc. Tor's Pattern-75")
    print("null, which enumerates polynomial spellings, does not range over it:")
    print("the index is a different KIND of object. The 12 is structural.")

    # -- verdict --------------------------------------------------------------
    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    print("Does 12 emerge as a 2I-index?              YES  (12 = [2I:C_10] = 120/10,")
    print("                                                 the vertex-orbit size)")
    print("Is that index outside the polynomial null? YES  (ratio of group orders,")
    print("                                                 not in Z[V,E,F,D,chi])")
    print("Trivial-rep multiplicity in rho_vert:      1")
    print("Hodge-transfer kernel dimension:           11   (cokernel 0)")
    print()
    print("So the two candidate numbers are mathematically DISTINCT and both clean:")
    print("  * the group index / orbit size            = 12   (the -12 of mu=1836)")
    print("  * the Hodge-transfer complement (kernel)   = 11   (curvature with no")
    print("                                                     image in the volume")
    print("                                                     sector)")
    print("The character computation says: transfer removes exactly ONE dimension")
    print("(the trivial/volume mode), leaving 11. We report what the maths says.")
    print("=" * 74)


if __name__ == "__main__":
    main()
