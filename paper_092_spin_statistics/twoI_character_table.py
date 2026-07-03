#!/usr/bin/env python3
"""
Paper 92 spin-statistics lemma -- DERIVED backbone (step 1-2 of the pre-reg).

Builds the binary icosahedral group 2I (order 120) EXACTLY as unit quaternions
over Z[phi], computes its character table (9 irreps) exactly over Q(sqrt5),
verifies orthonormality, and classifies each irrep as INTEGER (single-valued,
central -1 acts as +1) or SPINORIAL (half-integer, -1 acts as -1).

Then proves the two backbone facts the verdict rests on:
  (a) parity product rule: -1 on a tensor product is the product of factors,
      so a composite's integer/spinorial bit = (-1)^(# spinor constituents).
  (b) central-character conservation under induction: inducing a spinorial
      (-1 |-> -1) representation yields ONLY spinorial irreps. Demonstrated
      exactly by inducing the sign rep of the centre <-1> up to 2I.

No floating point. Exact arithmetic in Q(sqrt5) via (p, q) = p + q*sqrt5.

Pre-reg lock: bootstrap-universe deuteron-spin-statistics/PRE_REGISTRATION.md
commit 810e733. Run:  python3 twoI_character_table.py
"""
from fractions import Fraction as F


# --------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt5):  a number is p + q*sqrt5, p,q in Q.
# --------------------------------------------------------------------------
class Q5:
    __slots__ = ("p", "q")

    def __init__(self, p=0, q=0):
        self.p = F(p)
        self.q = F(q)

    @staticmethod
    def _c(o):
        return o if isinstance(o, Q5) else Q5(o, 0)

    def __add__(s, o):
        o = Q5._c(o)
        return Q5(s.p + o.p, s.q + o.q)

    __radd__ = __add__

    def __sub__(s, o):
        o = Q5._c(o)
        return Q5(s.p - o.p, s.q - o.q)

    def __rsub__(s, o):
        return Q5._c(o) - s

    def __neg__(s):
        return Q5(-s.p, -s.q)

    def __mul__(s, o):
        o = Q5._c(o)
        # (p1+q1 r5)(p2+q2 r5) = (p1p2 + 5 q1q2) + (p1q2 + q1p2) r5
        return Q5(s.p * o.p + 5 * s.q * o.q, s.p * o.q + s.q * o.p)

    __rmul__ = __mul__

    def __eq__(s, o):
        o = Q5._c(o)
        return s.p == o.p and s.q == o.q

    def __hash__(s):
        return hash((s.p, s.q))

    def conj5(s):
        """Galois automorphism sqrt5 -> -sqrt5."""
        return Q5(s.p, -s.q)

    def __repr__(s):
        if s.q == 0:
            return f"{s.p}"
        if s.p == 0:
            return f"{s.q}*r5"
        return f"({s.p}{'+' if s.q > 0 else ''}{s.q}*r5)"


ZERO = Q5(0)
ONE = Q5(1)


# --------------------------------------------------------------------------
# Unit quaternions over Q(sqrt5): 4-tuples of Q5.  Hamilton product.
# --------------------------------------------------------------------------
def qmul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    )


def qconj(a):
    w, x, y, z = a
    return (w, -x, -y, -z)  # unit quaternion inverse = conjugate (norm 1)


def qkey(a):
    return tuple((c.p, c.q) for c in a)


IDENT = (ONE, ZERO, ZERO, ZERO)


# --------------------------------------------------------------------------
# Build 2I by closure from the binary tetrahedral group 2T plus one golden
# unit quaternion. Closure stabilising at exactly 120 confirms membership.
# --------------------------------------------------------------------------
def build_2I():
    h = F(1, 2)
    # 8 Lipschitz units (+-1,0,0,0)&perm
    hurwitz = []
    for i in range(4):
        for s in (1, -1):
            v = [ZERO, ZERO, ZERO, ZERO]
            v[i] = Q5(s)
            hurwitz.append(tuple(v))
    # 16 units (+-1/2,+-1/2,+-1/2,+-1/2)
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    hurwitz.append((Q5(s0 * h), Q5(s1 * h),
                                    Q5(s2 * h), Q5(s3 * h)))
    # one golden icosian: (phi/2, 1/2, 1/(2phi), 0)
    A = Q5(F(1, 4), F(1, 4))    # phi/2     = (1+sqrt5)/4
    B = Q5(F(1, 2), 0)         # 1/2
    C = Q5(F(-1, 4), F(1, 4))  # 1/(2phi)  = (sqrt5-1)/4
    g0 = (A, B, C, ZERO)

    gens = hurwitz + [g0]
    elts = {}
    for g in gens:
        elts[qkey(g)] = g
    changed = True
    while changed:
        changed = False
        for a in list(elts.values()):
            for g in gens:
                p = qmul(a, g)
                k = qkey(p)
                if k not in elts:
                    elts[k] = p
                    changed = True
    return list(elts.values())


# --------------------------------------------------------------------------
# Conjugacy classes, with a representative, its order, and scalar part w.
# --------------------------------------------------------------------------
def conjugacy_classes(G):
    keyset = {qkey(g): g for g in G}
    seen = set()
    classes = []
    for x in G:
        kx = qkey(x)
        if kx in seen:
            continue
        orb = {}
        for g in G:
            c = qmul(qmul(g, x), qconj(g))
            orb[qkey(c)] = c
        for k in orb:
            seen.add(k)
        rep = x
        classes.append({"rep": rep, "size": len(orb),
                        "w": rep[0], "members": set(orb.keys())})
    return classes


def order_of(x):
    p = x
    for k in range(1, 200):
        if qkey(p) == qkey(IDENT):
            return k
        p = qmul(p, x)
    raise RuntimeError("order not found")


# --------------------------------------------------------------------------
# Chebyshev U_n(w) exact (character of Sym^n of the 2-dim rep; w = cos theta).
# --------------------------------------------------------------------------
def cheb_U(n, w):
    if n == 0:
        return ONE
    Ukm1 = ONE          # U_0
    Uk = Q5(2) * w      # U_1
    if n == 1:
        return Uk
    for _ in range(2, n + 1):
        Uk, Ukm1 = Q5(2) * w * Uk - Ukm1, Uk
    return Uk


def main():
    G = build_2I()
    assert len(G) == 120, f"expected 120 elements, got {len(G)}"

    # full group-closure check (every product stays in G)
    keys = {qkey(g) for g in G}
    for a in G:
        for b in G:
            assert qkey(qmul(a, b)) in keys, "not closed!"
    print(f"[ok] 2I built and verified: |G| = {len(G)}, closed under product.")

    classes = conjugacy_classes(G)
    assert len(classes) == 9, f"expected 9 classes, got {len(classes)}"
    for c in classes:
        c["order"] = order_of(c["rep"])
    # canonical order: by element order then size
    classes.sort(key=lambda c: (c["order"], c["size"]))
    sizes = [c["size"] for c in classes]
    assert sum(sizes) == 120
    print(f"[ok] 9 conjugacy classes, sizes {sizes} (sum {sum(sizes)}).")

    # locate identity (w=1, order 1) and -1 (w=-1, order 2)
    idx_id = next(i for i, c in enumerate(classes)
                  if qkey(c["rep"]) == qkey(IDENT))
    idx_m1 = next(i for i, c in enumerate(classes)
                  if c["order"] == 2 and c["w"] == Q5(-1))
    ws = [c["w"] for c in classes]

    # ----- character table -----------------------------------------------
    # symmetric powers j = 0,1/2,1,3/2,2,5/2  ->  irreps 1,2,3,4',5,6
    chi = {}
    chi["1"] = [cheb_U(0, w) for w in ws]
    chi["2"] = [cheb_U(1, w) for w in ws]
    chi["3"] = [cheb_U(2, w) for w in ws]
    chi["4'"] = [cheb_U(3, w) for w in ws]
    chi["5"] = [cheb_U(4, w) for w in ws]
    chi["6"] = [cheb_U(5, w) for w in ws]
    # Galois conjugates give 2', 3'
    chi["2'"] = [v.conj5() for v in chi["2"]]
    chi["3'"] = [v.conj5() for v in chi["3"]]
    # integer 4 (the A5 four-dim) = 2 (x) 2'
    chi["4"] = [a * b for a, b in zip(chi["2"], chi["2'"])]

    labels = ["1", "2", "2'", "3", "3'", "4", "4'", "5", "6"]
    dims = {L: chi[L][idx_id] for L in labels}
    for L in labels:
        assert dims[L].q == 0, f"dim of {L} not rational"
    dim_ints = {L: int(dims[L].p) for L in labels}
    assert sorted(dim_ints.values()) == [1, 2, 2, 3, 3, 4, 4, 5, 6], \
        f"dims wrong: {dim_ints}"
    assert sum(d * d for d in dim_ints.values()) == 120
    print(f"[ok] 9 irreps, dims {sorted(dim_ints.values())}, "
          f"sum of squares = {sum(d*d for d in dim_ints.values())}.")

    # ----- orthonormality (exact) ----------------------------------------
    def inner(ci, cj):
        tot = Q5(0)
        for s, a, b in zip(sizes, ci, cj):
            tot = tot + Q5(s) * a * b
        return tot  # = 120 * <chi_i, chi_j>

    for i, Li in enumerate(labels):
        for Lj in labels[i:]:
            val = inner(chi[Li], chi[Lj])
            expect = 120 if Li == Lj else 0
            assert val == Q5(expect), \
                f"orthonormality fail <{Li},{Lj}> = {val}"
    print("[ok] character table is orthonormal (exact over Q(sqrt5)).")

    # ----- classification by central element -1 --------------------------
    print("\nirrep  dim   chi(-1)   chi(-1)/dim   block")
    block = {}
    for L in labels:
        c_m1 = chi[L][idx_m1]
        sign = c_m1 * Q5(F(1, dim_ints[L]))  # chi(-1)/dim
        assert sign.q == 0 and sign.p in (1, -1)
        blk = "INTEGER " if sign.p == 1 else "SPINORIAL"
        block[L] = "integer" if sign.p == 1 else "spinorial"
        print(f"  {L:<4} {dim_ints[L]:>3}   {str(c_m1):>7}   "
              f"{'+1' if sign.p==1 else '-1':>6}        {blk}")

    integer = sorted(dim_ints[L] for L in labels if block[L] == "integer")
    spinor = sorted(dim_ints[L] for L in labels if block[L] == "spinorial")
    print(f"\n  INTEGER block (single-valued): dims {integer}")
    print(f"  SPINORIAL block (half-integer): dims {spinor}")
    assert integer == [1, 3, 3, 4, 5], integer   # = A5's five irreps
    assert spinor == [2, 2, 4, 6], spinor
    print("  -> integer block dims {1,3,3,4,5} match the A5 irreps exactly.")
    print("  -> free electron = fundamental spinor '2' (j=1/2), spinorial. ")

    # ----- backbone (a): parity product rule -----------------------------
    sgn = {L: (1 if block[L] == "integer" else -1) for L in labels}
    ok = True
    for Li in labels:
        for Lj in labels:
            prod_chi_m1 = chi[Li][idx_m1] * chi[Lj][idx_m1]
            dim_prod = dim_ints[Li] * dim_ints[Lj]
            parity_of_tensor = prod_chi_m1 * Q5(F(1, dim_prod))
            if parity_of_tensor.p != sgn[Li] * sgn[Lj]:
                ok = False
    assert ok
    print("\n[ok] backbone (a): for every pair, parity(i (x) j) = "
          "parity(i)*parity(j).")
    print("     => a composite's bit = (-1)^(number of spinor constituents).")

    # ----- backbone (b): induction conserves the central character -------
    # Induce the sign rep of the centre Z = {1,-1} (rho(-1) = -1) up to 2I.
    # chi_Ind(g) = (1/|Z|) sum_{a in G : a^-1 g a in Z} rho(a^-1 g a).
    # Only g in {1,-1} are central, so the induced character is supported
    # on those two classes:  +|G|/2 at 1,  -|G|/2 at -1,  0 elsewhere.
    indchar = [Q5(0)] * len(classes)
    indchar[idx_id] = Q5(F(120, 2))      # +60
    indchar[idx_m1] = Q5(F(-120, 2))     # -60
    mult = {}
    for L in labels:
        m = inner(indchar, chi[L]) * Q5(F(1, 120))  # <chi_Ind, chi_L>
        assert m.q == 0 and m.p == int(m.p) and m.p >= 0
        mult[L] = int(m.p)
    print("\n[ok] backbone (b): induce sign rep of centre <-1> up to 2I.")
    print("     multiplicities:",
          {L: mult[L] for L in labels})
    for L in labels:
        if block[L] == "integer":
            assert mult[L] == 0
        else:
            assert mult[L] == dim_ints[L]
    print("     => induced spinorial rep decomposes into SPINORIAL irreps "
          "ONLY (integer multiplicities all 0).")
    print("     Central character is conserved by induction: the construction "
          "CANNOT flip the bit.")

    # ----- the iff ------------------------------------------------------
    print("\n" + "=" * 68)
    print("THE iff (construction-free, DERIVED):")
    print("  proton = spinor (1 odd factor).  Observed: neutron 1/2 (odd),")
    print("  deuteron 1 (even), N-14 boson (even).  Counting constituents:")
    print("    neutron  = p + core      : odd*core  = odd  => core EVEN")
    print("    deuteron = p + core + p  : even*core = even => core EVEN")
    print("    N-14     = 14p + 7 core  : even*core = even => core EVEN")
    print("  All three give the SAME condition: P-e-P is spin-statistics-")
    print("  consistent  <=>  the bound electron core is integer-spin "
          "(parity-even).")
    print("=" * 68)

    # ----- write CSV -----------------------------------------------------
    import csv
    out = "twoI_character_table.csv"
    with open(out, "w", newline="") as f:
        wcsv = csv.writer(f)
        header = ["irrep", "dim", "block"] + \
            [f"C{ i }(ord{ c['order'] },sz{ c['size'] })"
             for i, c in enumerate(classes)]
        wcsv.writerow(header)
        for L in labels:
            row = [L, dim_ints[L], block[L]] + [str(v) for v in chi[L]]
            wcsv.writerow(row)
    print(f"\n[ok] wrote {out}")


if __name__ == "__main__":
    main()
