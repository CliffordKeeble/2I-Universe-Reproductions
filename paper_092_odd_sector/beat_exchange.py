#!/usr/bin/env python3
"""
Paper 92 odd sector -- Stage 3: the beat (spin-0 dressing) and the exchange
statistics (deuteron vs dineutron). Exact over Q(sqrt5); reuses the backbone.

These computations are CONDITIONAL on Stage 2 (proton = spinor 2). They are
reported as rep-theory facts that hold IF the proton is the spinor 2; the
grounding of that "if" is Stage 2's open crux (see findings.md).

  (3a) the beat: decompose 2 (x) 2* and 2 (x) 2'. The pre-registered open
       question -- which conjugate yields the scalar [1] -- is answered by
       computation, not assumed.
  (3b) deuteron vs dineutron: decompose 2 (x) 2 into symmetric (Sym^2) and
       antisymmetric (Lambda^2) parts; test whether p-n-binds / n-n-doesn't
       falls out of Fermi exchange or must be fitted.

Run:  python3 beat_exchange.py     (needs sibling paper_092_spin_statistics)
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "paper_092_spin_statistics"))
from twoI_character_table import (   # noqa: E402
    Q5, build_2I, conjugacy_classes, cheb_U, qkey, IDENT, order_of,
)

HALF = Q5(F(1, 2))


def build_table():
    G = build_2I()
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
    block = {L: ("integer" if (chi[L][idx_m1] * Q5(F(1, dims[L]))).p == 1
                 else "spinorial") for L in labels}
    return dict(sizes=sizes, ws=ws, idx_id=idx_id, idx_m1=idx_m1,
                chi=chi, labels=labels, dims=dims, block=block)


def decompose(T, ch):
    out = {}
    for L in T["labels"]:
        tot = Q5(0)
        for s, a, b in zip(T["sizes"], ch, T["chi"][L]):
            tot = tot + Q5(s) * a * b
        v = tot * Q5(F(1, 120))
        assert v.q == 0 and v.p >= 0 and v.p == int(v.p), f"{L}:{v}"
        if v.p:
            out[L] = int(v.p)
    return out


def fmt(dec):
    return " + ".join(f"{c}*[{L}]" if c > 1 else f"[{L}]"
                      for L, c in dec.items())


def main():
    T = build_table()
    chi, ws = T["chi"], T["ws"]
    # squared-element character of the 2-dim rep: chi_2(g^2) = 2 cos 2a = 2(2w^2-1)
    chi2 = chi["2"]
    chi2_sq = [Q5(2) * (Q5(2) * w * w - Q5(1)) for w in ws]   # at g^2

    print("=" * 70)
    print("(3a) THE BEAT -- which conjugate yields the spin-0 scalar [1]?")
    print("  The proton is self-dual as a rep: 2* = 2 (its character 2cos(a) is")
    print("  real), so the representation-dual beat is 2 (x) 2.")
    two_dual = [a * a for a in chi2]            # 2 (x) 2*  =  2 (x) 2
    two_galois = [a * b for a, b in zip(chi2, chi["2'"])]   # 2 (x) 2'
    dec_dual = decompose(T, two_dual)
    dec_gal = decompose(T, two_galois)
    print(f"    2 (x) 2*  = 2 (x) 2  = {fmt(dec_dual)}")
    print(f"    2 (x) 2'  (Galois)   = {fmt(dec_gal)}")
    has_dual = dec_dual.get("1", 0)
    has_gal = dec_gal.get("1", 0)
    assert has_dual == 1 and has_gal == 0
    print()
    print("    => the spin-0 scalar [1] comes from the REPRESENTATION DUAL")
    print("       (2 (x) 2* = 2 (x) 2 contains [1] once), NOT from the GALOIS")
    print("       conjugate (2 (x) 2' = [4], no trivial). The brief's named")
    print("       '2 (x) 2'' does not carry the scalar; mode x its own dual does.")
    print("    => spin-0 beat-dressing EXISTS and is DERIVED (the |amplitude|^2")
    print("       scalar), but it is the dual beat, not the Galois beat. The")
    print("       Galois symmetry (Stage 1) and the scalar beat are DIFFERENT")
    print("       true facts; the picture conflated them. (hard-rule new clause)")

    print("\n" + "=" * 70)
    print("(3b) DEUTERON vs DINEUTRON -- predict or relabel?")
    # Sym^2 and Lambda^2 of the spinor 2
    sym2 = [(a * a + b) * HALF for a, b in zip(chi2, chi2_sq)]   # Sym^2
    alt2 = [(a * a - b) * HALF for a, b in zip(chi2, chi2_sq)]   # Lambda^2
    dec_sym = decompose(T, sym2)
    dec_alt = decompose(T, alt2)
    print(f"    Sym^2(2)    (symmetric, spatial-even pair) = {fmt(dec_sym)}")
    print(f"    Lambda^2(2) (antisymmetric)                = {fmt(dec_alt)}")
    # SU(2): Sym^2(V_1/2)=V_1=[3] triplet (spin 1); Lambda^2=V_0=[1] singlet (0)
    assert dec_sym == {"3": 1} and dec_alt == {"1": 1}
    print("    => Sym^2(2) = [3]  : the spin-1 TRIPLET (symmetric).")
    print("       Lambda^2(2) = [1]: the spin-0 SINGLET (antisymmetric).")
    print()
    print("    Two identical spin-1/2 fermions (Fermi, from the spin lemma:")
    print("    proton/neutron = spinor 2, parity-odd) need a TOTALLY")
    print("    ANTISYMMETRIC state. In s-wave (symmetric spatial) the spin part")
    print("    must be antisymmetric = the singlet [1] (spin-0). The bound")
    print("    nuclear channel is the spin-1 triplet [3] (the deuteron's 3S1).")
    print()
    print("      p-n (distinguishable): can occupy the [3] triplet -> BINDS (deuteron)")
    print("      n-n (identical)      : s-wave forced to [1] singlet -> NO 3S1 -> unbound")
    print()
    print("    VERDICT (3b): the contrast FALLS OUT of Fermi exchange + the")
    print("    Sym^2/Lambda^2 split [3]/[1] -- PREDICTION from proton=spinor-2,")
    print("    NOT a relabel. No 'beat commensurability' condition is needed;")
    print("    imposing one would BE the relabel. (This is the textbook Pauli")
    print("    reason, now grounded in the 2I spinor assignment -- conditional")
    print("    on Stage 2 forcing proton=spinor-2, which it does not, so the")
    print("    grounding rides on that open input.)")
    print()
    print("    3-body caveat (honest): 3H (pnn) is BOUND (beta-unstable, not")
    print("    unbound); the brief's 'p-n-n unstable' is beta-stability, a")
    print("    different question. The clean, decided test is the 2-body case.")

    print("\n" + "=" * 70)
    print("STAGE 3 SUMMARY")
    print("  3a: spin-0 beat is real and DERIVED -- but via the dual (2(x)2),")
    print("      not the Galois conjugate (2(x)2'=[4]). Correction logged.")
    print("  3b: deuteron-binds/dineutron-doesn't falls out of Fermi exchange")
    print("      ([3] triplet sym / [1] singlet antisym) -- prediction, not")
    print("      relabel -- conditional on proton=spinor-2 (Stage 2, unforced).")
    print("=" * 70)


if __name__ == "__main__":
    main()
