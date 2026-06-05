#!/usr/bin/env python3
"""
Paper 92 -- the Hodge-dual reframe of the bridge node. Compute on the LOCKED
inputs (PRE-REGISTRATION.md, this dir; lock commit recorded there).

Locked from the corpus, before this ran:
  M = S^3/2I, n=3, Riemannian; proton = 0-form (Paper 121 Table 2); and the
  corpus's own Hodge star *(proton) = electron (Paper 101 s4.3, Lemma L).

This script:
  (1) base-case assert  **=+1 on n=3 Riemannian for all k (W-002);
  (2) reads *(proton) = 3-form = electron -> moment scale (Gate 1(c), via D1);
  (3) tests the form-character claim (Gate 4, via D2): builds 2I (reused from
      ../paper92-spin-statistics), the vector rep V3, the exterior powers
      Lambda^k(V3) and the full form bundle Lambda*(V3); proves the form bundle
      carries ZERO spinorial content -- so "node is a form" => integer-spin is
      true for GENUINE forms, but the spin-1/2 proton/electron (the spinorial
      irrep '2', ABSENT from Lambda*) are not genuine forms;
  (4) Dirac-Kahler check;
  (5) states the landing.

Exact arithmetic over Q(sqrt5); every identity an assert. Reuses the 2I machine.
Run:  python3 hodge_dual.py        (needs the sibling paper92-spin-statistics dir)
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
SIB = os.path.join(HERE, "..", "paper92-spin-statistics")
sys.path.insert(0, SIB)
from twoI_character_table import (   # noqa: E402  (reuse the locked backbone)
    Q5, build_2I, conjugacy_classes, cheb_U, qkey, IDENT, order_of,
)

HALF = Q5(F(1, 2))
THIRD = Q5(F(1, 3))


# --------------------------------------------------------------------------
# (1) Hodge involution sign on the locked manifold.
# --------------------------------------------------------------------------
def involution_sign(k, n, t):
    """**  = (-1)^{k(n-k)+t} * id  on k-forms."""
    return (-1) ** (k * (n - k) + t)


def section_involution():
    print("=" * 70)
    print("(1) BASE CASE (W-002): Hodge involution on M = S^3/2I, n=3, Riemannian")
    n, t = 3, 0   # Riemannian => t = 0 minus signs
    for k in range(n + 1):
        s = involution_sign(k, n, t)
        print(f"    **  on {k}-forms = {s:+d} * id   "
              f"[(-1)^(k(n-k)+t), k(n-k)={k*(n-k)}]")
        assert s == 1, "involution sign != +1 -- (M,signature,degree) malformed"
    print("    => ** = +1 for all k: * is a real involution on n=3. [ok]")
    return n


# --------------------------------------------------------------------------
# 2I backbone (reused) + the 9-irrep character table, rebuilt inline.
# --------------------------------------------------------------------------
def build_table():
    G = build_2I()
    assert len(G) == 120
    classes = conjugacy_classes(G)
    for c in classes:
        c["order"] = order_of(c["rep"])
    classes.sort(key=lambda c: (c["order"], c["size"]))
    sizes = [c["size"] for c in classes]
    ws = [c["w"] for c in classes]
    idx_id = next(i for i, c in enumerate(classes)
                  if qkey(c["rep"]) == qkey(IDENT))
    idx_m1 = next(i for i, c in enumerate(classes)
                  if c["order"] == 2 and c["w"] == Q5(-1))
    chi = {}
    chi["1"] = [cheb_U(0, w) for w in ws]
    chi["2"] = [cheb_U(1, w) for w in ws]
    chi["3"] = [cheb_U(2, w) for w in ws]
    chi["4'"] = [cheb_U(3, w) for w in ws]
    chi["5"] = [cheb_U(4, w) for w in ws]
    chi["6"] = [cheb_U(5, w) for w in ws]
    chi["2'"] = [v.conj5() for v in chi["2"]]
    chi["3'"] = [v.conj5() for v in chi["3"]]
    chi["4"] = [a * b for a, b in zip(chi["2"], chi["2'"])]
    labels = ["1", "2", "2'", "3", "3'", "4", "4'", "5", "6"]
    dims = {L: int(chi[L][idx_id].p) for L in labels}
    block = {}
    for L in labels:
        sign = chi[L][idx_m1] * Q5(F(1, dims[L]))
        block[L] = "integer" if sign.p == 1 else "spinorial"
    return dict(classes=classes, sizes=sizes, ws=ws, idx_id=idx_id,
                idx_m1=idx_m1, chi=chi, labels=labels, dims=dims, block=block)


def inner(T, ci, cj):
    """<ci, cj> = (1/120) sum_classes size * ci * cj  (real characters)."""
    tot = Q5(0)
    for s, a, b in zip(T["sizes"], ci, cj):
        tot = tot + Q5(s) * a * b
    return tot * Q5(F(1, 120))


def decompose(T, char):
    mult = {}
    for L in T["labels"]:
        m = inner(T, char, T["chi"][L])
        assert m.q == 0 and m.p == int(m.p) and m.p >= 0, \
            f"non-integer multiplicity for {L}: {m}"
        mult[L] = int(m.p)
    return mult


# --------------------------------------------------------------------------
# (2) *(proton) and Gate 1(c).
# --------------------------------------------------------------------------
def section_gate1c(n):
    print("\n" + "=" * 70)
    print("(2) *(proton) AND GATE 1(c) -- magnetic moment (dictionary D1)")
    k_proton = 0                       # locked: Paper 121 Table 2
    k_dual = n - k_proton              # Hodge star: 0 -> 3
    print(f"    proton form-degree k = {k_proton} (0-form/H^0, locked).")
    print(f"    *(proton): degree {k_proton} -> n-k = {k_dual} (3-form/H^3).")
    print("    corpus identity of the 3-form: ELECTRON "
          "(Paper 101 s4.3; Lemma L: *:Omega^0->Omega^3, proton<->electron).")
    print("    => the Hodge dual of the proton IS the electron, not a new "
          "proton-scale object.")
    # D1: moment scale = that of the dual object (the electron) = Bohr.
    mp_over_me = 1836.15267343
    mu_B_in_muN = mp_over_me           # mu_B = (m_p/m_e) mu_N
    mu_d_measured = 0.8574382          # nuclear magnetons (deuteron)
    mu_p, mu_n = 2.792847, -1.913043
    print(f"\n    D1 moment scale of *(proton)=electron: mu_B = (m_p/m_e) mu_N "
          f"= {mu_B_in_muN:.1f} mu_N")
    print(f"    deuteron measured moment           : {mu_d_measured:.4f} mu_N")
    print(f"    conventional mu_p + mu_n           : {mu_p+mu_n:.4f} mu_N")
    print(f"    miss factor (electron node / target): "
          f"{mu_B_in_muN/mu_d_measured:,.0f}x")
    print("    => GATE 1(c) NOT CURED. The Hodge dual reproduces the electron")
    print("       core and its ~1836 mu_N (Bohr) moment -- the 3-orders-of-")
    print("       magnitude problem the reframe hoped to remove.")
    return k_proton, k_dual


# --------------------------------------------------------------------------
# (3) Gate 4 -- the form-character claim (dictionary D2).
# --------------------------------------------------------------------------
def section_gate4(T):
    print("\n" + "=" * 70)
    print("(3) GATE 4 -- does 'node is a form' confer integer spin? (D2)")

    chi = T["chi"]
    ws = T["ws"]
    idx_id, idx_m1 = T["idx_id"], T["idx_m1"]

    # V3 = the vector (rotation) rep = the integer irrep '3' = U_2(w).
    chiV3 = chi["3"]
    assert (chiV3[idx_m1] * Q5(F(1, 3))).p == 1, "V3 must be integer-spin"

    # Newton's identities for exterior powers of the 3-dim rep.
    # p_m = chi_V3(g^m); scalar part of g^m has cos(m*half) = Cheb-T_m(w).
    def scalar_pow(w, m):
        # cos(m*theta) via Chebyshev T: T0=1, T1=w, T_{k+1}=2w T_k - T_{k-1}
        Tkm1, Tk = Q5(1), w
        if m == 0:
            return Q5(1)
        if m == 1:
            return w
        for _ in range(2, m + 1):
            Tk, Tkm1 = Q5(2) * w * Tk - Tkm1, Tk
        return Tk

    lam = {0: [], 1: [], 2: [], 3: []}  # exterior-power characters
    for j, w in enumerate(ws):
        p1 = cheb_U(2, scalar_pow(w, 1))
        p2 = cheb_U(2, scalar_pow(w, 2))
        p3 = cheb_U(2, scalar_pow(w, 3))
        e0 = Q5(1)
        e1 = p1
        e2 = (e1 * p1 - p2) * HALF
        e3 = (e2 * p1 - e1 * p2 + p3) * THIRD
        lam[0].append(e0)
        lam[1].append(e1)
        lam[2].append(e2)
        lam[3].append(e3)

    # exact identities: Lambda^1 = V3, Lambda^2 = V3, Lambda^3 = trivial(det=+1)
    assert lam[1] == chiV3
    assert lam[2] == chiV3, "Lambda^2(V3) != V3"
    assert lam[3] == chi["1"], "Lambda^3(V3) != trivial (det != +1)"
    print("    exterior powers of V3 (exact): "
          "Lambda^0=1, Lambda^1=3, Lambda^2=3, Lambda^3=1 (det=+1). [ok]")

    # full form bundle Lambda*(V3) = sum_k Lambda^k
    lam_star = [lam[0][j] + lam[1][j] + lam[2][j] + lam[3][j]
                for j in range(len(ws))]
    dim_star = int(lam_star[idx_id].p)
    assert dim_star == 8
    par = int((lam_star[idx_m1] * Q5(F(1, dim_star))).p)
    print(f"    form bundle Lambda*(V3): dim {dim_star}, "
          f"central -1 acts as {par:+d} -> "
          f"{'INTEGER (parity-even)' if par == 1 else 'spinorial'}.")

    # decompose Lambda* and each Lambda^k; show ZERO spinorial content
    dec = decompose(T, lam_star)
    print(f"    Lambda*(V3) decomposition: "
          f"{ {L: dec[L] for L in T['labels'] if dec[L]} }")
    spin_content = sum(dec[L] for L in T["labels"]
                       if T["block"][L] == "spinorial")
    assert spin_content == 0
    print("    => the de Rham form bundle carries ZERO spinorial content.")
    print("       A GENUINE k-form is parity-even: 'node is a form => integer-")
    print("       spin' is TRUE -- for genuine forms.")

    # the physical spin-1/2 particle is the spinorial '2', ABSENT from Lambda*
    spin2 = int((chi['2'][idx_m1] * Q5(F(1, 2))).p)
    print(f"\n    BUT the proton/electron are spin-1/2 = the spinorial irrep "
          f"'2' (central -1 acts as {spin2:+d}),")
    print(f"       which appears in Lambda*(V3) with multiplicity "
          f"{dec[chr(50)] if '2' in dec else 0}.")
    print("    => representing a spin-1/2 particle 'as a form' places it in a")
    print("       bundle with no spinorial content -- stripping its spin-1/2.")
    print("       The form-degree label (cohomology sector) and the spin label")
    print("       (2I central character) are ORTHOGONAL classifications.")
    return spin_content


# --------------------------------------------------------------------------
# (4) Dirac-Kahler check.
# --------------------------------------------------------------------------
def section_dirac_kahler():
    print("\n" + "=" * 70)
    print("(4) DIRAC-KAHLER CHECK (before any Landing-3 call)")
    print("    Dirac-Kahler identifies the Clifford module (Dirac spinors) with")
    print("    the inhomogeneous form space Lambda* under the LEFT Clifford")
    print("    action -- NOT the wedge/de-Rham form-degree action computed in")
    print("    (3). The physical fermion's spin lives in the left action")
    print("    (spinorial, parity-odd); the form-degree k is a separate grading")
    print("    (parity-even, shown above). DK therefore relates forms and")
    print("    spinors but does NOT make a homogeneous k-form an integer-spin")
    print("    physical fermion: it confirms the two gradings are distinct,")
    print("    exactly the boundary. The clash is real, not an artefact of")
    print("    picking the wrong construction.")


def main():
    n = section_involution()
    T = build_table()
    section_gate1c(n)
    spin_content = section_gate4(T)
    section_dirac_kahler()

    print("\n" + "=" * 70)
    print("LANDING")
    print("=" * 70)
    print("  LANDING 3 -- CLASH BITES.")
    print()
    print("  Gate 1(c): NOT cured. *(proton) = electron (corpus, robust under")
    print("    the matter-block swap) -> Bohr-scale moment ~1836 mu_N, three")
    print("    orders off the deuteron's 0.857 mu_N. The reframe reproduces the")
    print("    electron-core reading; it yields no proton-scale dual. (And k is")
    print("    not uniquely pinned -- swap 'chosen not forced' -- an independent")
    print("    Landing-4 reason 1(c) cannot be claimed cured.)")
    print()
    print("  Gate 4: clash bites. The form bundle on S^3/2I has zero spinorial")
    print("    content, so 'node is a form => integer-spin' holds only for")
    print("    GENUINE forms. The spin-1/2 proton/electron are the spinorial")
    print("    '2', absent from Lambda*; calling them 'forms' strips their")
    print("    spin-1/2. Dirac-Kahler confirms form-degree and spin are distinct")
    print("    gradings -- it does not rescue the reframe.")
    print()
    print("  Net: the Hodge reframe does NOT improve on the spin lemma's")
    print("  neutrino-ejection result (W-102 stays conditional) and does NOT")
    print("  cure 1(c). The electron-core reading and its 1836x magnetic-moment")
    print("  problem stand as the likely-fatal gate -- the paper's ceiling,")
    print("  exactly as the brief anticipated. Landing 3, not 1.")
    print("=" * 70)
    assert spin_content == 0


if __name__ == "__main__":
    main()
