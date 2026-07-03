"""
product_divisor.py -- Paper 101 v3.1, Mr A's product-divisor gate (fifth star).

Mr A: verify the divisor of the PRODUCT f_12 * f_30, not just the factor f_12.
Confirm the blocked-in-the-mass-channel locus is the 12-vertex orbit, with the
edge-vanishing of f_30 accounted for.

This script has three parts.

PART A (real R^3, exact Q(sqrt5)) -- and a FLAG.
  Task 2 built the real avatar of f_12 as the degree-6 icosahedral invariant
  f_vert (vanishes on the 12 vertices, nonzero on faces/edges). The brief asks to
  build f_30's avatar "by the same Reynolds construction". We show WHY that does not
  isolate the edges: the unique degree-15 real icosahedral invariant I_15 is (up to
  scale) the product of the 15 two-fold-axis linear forms = the product of the 15
  icosahedral mirror-plane forms. Every special axis -- vertex (5 mirrors), face
  (3 mirrors), edge (2 mirrors) -- lies on mirror planes, so I_15 vanishes on ALL
  THREE orbits. The real degree-15 invariant is therefore NOT a clean edge form.
  >>> The clean edge-form divisor must be stated in the binary picture (Part B). <<<

PART B (binary CP^1, exact over Z) -- the rigorous product divisor.
  The standard Klein generators f_12, H_20, T_30 (integer coefficients) have the 12
  vertices / 20 face-centres / 30 edge-midpoints as their root sets on CP^1. We verify
  the Hessian/Jacobian hierarchy, the icosahedral syzygy, squarefreeness, and pairwise
  COPRIMALITY. Coprimality => the three root sets are pairwise DISJOINT, so:
    div(f_12 . T_30) = (12 vertices) U (30 edges), disjoint;
    f_12 nonzero on the 30 edges, T_30 nonzero on the 12 vertices,
    and the product is nonzero on the 20 faces.

PART C -- the structural mass-channel statement (printed conclusion).

Run: python product_divisor.py
"""

import sympy as sp
from sympy import Matrix, sqrt, simplify, expand, gcd, diff, Poly, symbols
from vertex_index_2I import to_sympy, icosahedron_vertices

PHI = (1 + sqrt(5)) / 2


# ===========================================================================
# PART A -- the real R^3 degree-15 invariant vanishes on ALL special orbits
# ===========================================================================

def sym_vertices():
    return [Matrix([to_sympy(c) for c in v]) for v in icosahedron_vertices()]

def is_zero_vec(u):
    return all(simplify(c) == 0 for c in u)

def parallel(u, w):
    # cross product zero <=> parallel/antiparallel
    cx = Matrix([u[1]*w[2] - u[2]*w[1],
                 u[2]*w[0] - u[0]*w[2],
                 u[0]*w[1] - u[1]*w[0]])
    return is_zero_vec(cx)

def edge_dist2():
    return sp.Integer(4)   # |v_i - v_j|^2 for an icosahedron edge in these coords

def build_orbits():
    V = sym_vertices()
    n = len(V)
    # adjacency: edge iff |v_i - v_j|^2 == 4
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            d2 = simplify((V[i]-V[j]).dot(V[i]-V[j]))
            if d2 == edge_dist2():
                adj[i][j] = adj[j][i] = True

    # edge-midpoints (30): one per adjacent pair
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                edges.append((V[i] + V[j]))            # direction (scale irrelevant)

    # face-centres (20): centroids of mutually-adjacent triples. Each triangular
    # face is exactly one such triple, so no de-duplication is needed (antipodal
    # faces give distinct -- negated -- centroids and must both be kept).
    faces = []
    for i in range(n):
        for j in range(i+1, n):
            if not adj[i][j]:
                continue
            for k in range(j+1, n):
                if adj[i][k] and adj[j][k]:
                    faces.append(V[i] + V[j] + V[k])
    return V, edges, faces

def two_fold_axes(edges):
    """15 distinct lines (one per antipodal pair of the 30 edge-midpoints)."""
    axes = []
    for e in edges:
        if not any(parallel(e, a) for a in axes):
            axes.append(e)
    return axes

def I15_value(axes, p):
    """I_15(p) = product over the 15 two-fold axes of <a_k, p>  (no expansion)."""
    val = sp.Integer(1)
    for a in axes:
        val *= a.dot(p)
    return simplify(val)

def part_A():
    print("=" * 74)
    print("PART A -- real degree-15 icosahedral invariant I_15 (Reynolds-class)")
    print("=" * 74)
    V, edges, faces = build_orbits()
    print(f"orbits built (exact Q(sqrt5)): {len(V)} vertices, "
          f"{len(edges)} edge-midpoints, {len(faces)} face-centres")
    assert len(V) == 12 and len(edges) == 30 and len(faces) == 20

    axes = two_fold_axes(edges)
    print(f"two-fold axes (= mirror normals): {len(axes)}")
    assert len(axes) == 15

    # I_15 is the unique degree-15 invariant (= product of the 15 mirror forms).
    generic = Matrix([sp.Integer(1), sp.Integer(2), sp.Integer(3)])
    print(f"I_15 not identically zero (value at (1,2,3)): "
          f"{I15_value(axes, generic) != 0}")

    v_vals = [I15_value(axes, p) for p in V]
    e_vals = [I15_value(axes, p) for p in edges]
    f_vals = [I15_value(axes, p) for p in faces]
    print(f"\nI_15 = 0 on all 12 VERTICES      : {all(x == 0 for x in v_vals)}")
    print(f"I_15 = 0 on all 30 EDGE-midpoints: {all(x == 0 for x in e_vals)}")
    print(f"I_15 = 0 on all 20 FACE-centres  : {all(x == 0 for x in f_vals)}")
    print("\n>>> FLAG: the real degree-15 edge invariant vanishes on ALL THREE")
    print("    special orbits (every axis lies on icosahedral mirror planes:")
    print("    vertex in 5, face in 3, edge in 2). So the 'same Reynolds")
    print("    construction as Task 2' does NOT yield a clean edge form in R^3.")
    print("    The clean edge-form divisor is a binary (CP^1) statement -> Part B.")
    # the degree-6 avatar f_vert (Task 2) DOES isolate vertices -- recap the contrast
    print("\n    (Contrast, Task 2: the degree-6 avatar f_vert vanishes on the 12")
    print("     vertices only -- f_vert(face)=2/5, f_vert(edge)=1/80, both nonzero.)")
    return V, edges, faces


# ===========================================================================
# PART B -- binary Klein generators: coprime => disjoint divisors
# ===========================================================================

def part_B():
    print("\n" + "=" * 74)
    print("PART B -- binary Klein generators f_12, H_20, T_30 (exact over Z)")
    print("=" * 74)
    z1, z2 = symbols("z1 z2")

    f = z1*z2*(z1**10 + 11*z1**5*z2**5 - z2**10)                       # vertices, deg 12
    H = -(z1**20 + z2**20) + 228*(z1**15*z2**5 - z1**5*z2**15) \
        - 494*z1**10*z2**10                                            # faces,    deg 20
    T = (z1**30 + z2**30) + 522*(z1**25*z2**5 - z1**5*z2**25) \
        - 10005*(z1**20*z2**10 + z1**10*z2**20)                        # edges,    deg 30

    print(f"deg f_12 = {Poly(f, z1, z2).total_degree()}  "
          f"deg H_20 = {Poly(H, z1, z2).total_degree()}  "
          f"deg T_30 = {Poly(T, z1, z2).total_degree()}")

    # Hessian/Jacobian hierarchy: H ~ Hessian(f), T ~ Jacobian(f, H)
    Hess = expand(diff(f, z1, 2)*diff(f, z2, 2) - diff(f, z1, z2)**2)
    cH = simplify(Hess / H)
    Jac = expand(diff(f, z1)*diff(H, z2) - diff(f, z2)*diff(H, z1))
    cT = simplify(Jac / T)
    print(f"\nHessian(f_12) = ({cH}) * H_20        -> H_20 is the face/Hessian form")
    print(f"Jacobian(f_12,H_20) = ({cT}) * T_30  -> T_30 is the edge/Jacobian form")
    assert cH.free_symbols == set() and cT.free_symbols == set(), \
        "Hessian/Jacobian are not scalar multiples of H/T"

    # icosahedral syzygy: find the exact relation among f^5, H^3, T^2
    lhs = expand(H**3 + T**2)
    ratio = simplify(lhs / f**5)
    print(f"\nicosahedral syzygy:  H_20^3 + T_30^2 = ({ratio}) * f_12^5")
    assert ratio.free_symbols == set(), "syzygy is not a scalar multiple of f^5"

    # squarefree (distinct roots): for a homogeneous form F, a multiple root is a
    # common root of BOTH partials (Euler), so F is squarefree iff
    # gcd(dF/dz1, dF/dz2) is constant. A single partial is the wrong test here
    # (e.g. df/dz2 keeps f's z1 factor).
    def squarefree(F):
        return Poly(gcd(diff(F, z1), diff(F, z2)), z1, z2).total_degree() == 0
    print(f"\nf_12 squarefree (12 distinct roots): {squarefree(f)}")
    print(f"H_20 squarefree (20 distinct roots): {squarefree(H)}")
    print(f"T_30 squarefree (30 distinct roots): {squarefree(T)}")
    assert squarefree(f) and squarefree(H) and squarefree(T)

    # pairwise coprimality => disjoint root sets (disjoint divisors)
    g_fT = gcd(f, T)
    g_fH = gcd(f, H)
    g_HT = gcd(H, T)
    def is_const(g):
        return Poly(g, z1, z2).total_degree() == 0
    print(f"\ngcd(f_12, T_30) = {g_fT}   -> coprime: {is_const(g_fT)}  "
          f"(vertices disjoint from edges)")
    print(f"gcd(f_12, H_20) = {g_fH}   -> coprime: {is_const(g_fH)}  "
          f"(vertices disjoint from faces)")
    print(f"gcd(H_20, T_30) = {g_HT}   -> coprime: {is_const(g_HT)}  "
          f"(faces disjoint from edges)")
    assert is_const(g_fT) and is_const(g_fH) and is_const(g_HT)

    print("\n--- product divisor (the bare result Mr A asked for) ---")
    print("div(f_12 . T_30) = (12 vertices) U (30 edge-midpoints), DISJOINT union.")
    print("  on the 12 VERTICES : f_12 = 0  (T_30 != 0)  => product = 0")
    print("  on the 30 EDGES    : T_30 = 0  (f_12 != 0)  => product = 0")
    print("  on the 20 FACES    : f_12 != 0 and T_30 != 0 => product != 0")
    print("The edge-vanishing of the product is carried entirely by the T_30 factor,")
    print("on a locus (30 edges) DISJOINT from the 12-vertex orbit.")


# ===========================================================================
# PART C -- the structural mass-channel statement
# ===========================================================================

def part_C():
    print("\n" + "=" * 74)
    print("PART C -- does the edge-vanishing enter the mass channel?  (structural)")
    print("=" * 74)
    print("The vertex-curvature module rho_vert is supported on the VERTEX ORBIT")
    print("(the 12 vertex delta-functions; the permutation module C[2I/C_10] of")
    print("Task findings_minus12_index). The blocking inner product in the mass")
    print("channel is")
    print("      <f_12 . f_30 , rho_vert element>")
    print("           = sum over the 12 VERTICES of  f_12(v_i).f_30(v_i).c_i")
    print("           = 0   term by term,  because f_12(v_i) = 0 at every vertex")
    print("                                 (Task 2: div f_12 contains the 12 vertices).")
    print()
    print("The product's extra zeros on the 30 edges are a DISJOINT locus (Part B)")
    print("that the vertex-supported module never samples. So the edge-vanishing of")
    print("f_30 does NOT enlarge the blocked module: the blocked content in the mass")
    print("channel is exactly the 12-dimensional vertex orbit.")
    print()
    print("=> Mr A objection answered. Proposed one-sentence addition to section 5.1:")
    print("   'The product f_12.f_30 also vanishes on the 30 edge-midpoints via its")
    print("    f_30 factor; this locus is disjoint from the vertex orbit and is not")
    print("    sampled by the vertex-curvature module, so the blocked content in the")
    print("    mass channel remains exactly the 12-dimensional vertex orbit.'")


def main():
    part_A()
    part_B()
    part_C()
    print("\n" + "=" * 74)
    print("VERDICT: product-divisor gate PASSED. div(f_12.f_30) = vertices U edges")
    print("(disjoint); mass-channel blocking is vertex-supported => exactly 12-dim.")
    print("=" * 74)


if __name__ == "__main__":
    main()
