#!/usr/bin/env python3
"""
Paper 92 odd sector -- Stage 1: the Galois structure and the SPINORIAL spectrum
that Paper 174 skipped. Exact over Q(sqrt5); reuses the 2I backbone.

Paper 174 computed m_k = <chi_trivial, chi_{Sym^k V2}> = dim Sym^k(C^2)^{2I},
the BOSONIC (invariant) spectrum -- Molien (1-t^60)/((1-t^12)(1-t^20)(1-t^30)),
all even k (odd harmonics killed by the central -1). This script computes the
COMPLEMENTARY twisted Molien series for every irrep, in particular the four
SPINORIAL irreps {2, 2', 4', 6} that live at ODD k -- "where the actors live."

  (1a) Galois structure: 2<->2' over Q(sqrt5); Galois fixes the central
       character on every irrep, so it preserves the integer/spinorial split.
  (1b) m_{rho,k} = <chi_rho, chi_{Sym^k V2}> for all irreps; validate the
       trivial row against 174; closed-form generating functions; the
       integer-even / spinorial-odd selection rule.

Run:  python3 spinorial_spectrum.py     (needs sibling paper_092_spin_statistics)
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "paper_092_spin_statistics"))
from twoI_character_table import (   # noqa: E402  (reuse the locked backbone)
    Q5, build_2I, conjugacy_classes, cheb_U, qkey, IDENT, order_of,
)

KMAX = 200   # enough to confirm the covariant numerators truncate


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
        block[L] = "integer" if (chi[L][idx_m1] * Q5(F(1, dims[L]))).p == 1 \
            else "spinorial"
    return dict(classes=classes, sizes=sizes, ws=ws, idx_id=idx_id,
                idx_m1=idx_m1, chi=chi, labels=labels, dims=dims, block=block)


def stage1a_galois(T):
    print("=" * 70)
    print("(1a) GALOIS STRUCTURE over Q(sqrt5)")
    chi, labels, idx_m1 = T["chi"], T["labels"], T["idx_m1"]
    # Galois conjugate of the whole table = swap classes that carry sqrt5? No:
    # Galois acts on character VALUES; a row is Galois-fixed iff all its values
    # are rational. Compute each row's Galois image as a multiset-compatible row.
    galois_pairs = []
    fixed = []
    for L in labels:
        row_conj = [v.conj5() for v in chi[L]]
        # find which label has this conjugated row
        match = next(M for M in labels if chi[M] == row_conj)
        if match == L:
            fixed.append(L)
        else:
            galois_pairs.append((L, match))
    pairset = {frozenset(p) for p in galois_pairs}
    print(f"    Galois-swapped pairs : "
          f"{sorted(tuple(sorted(p)) for p in pairset)}")
    print(f"    Galois-fixed (rational char): {fixed}")
    assert frozenset({"2", "2'"}) in pairset, "2 and 2' must be Galois conjugates"
    assert frozenset({"3", "3'"}) in pairset
    # central character is rational (Galois-fixed) on every irrep
    for L in labels:
        assert chi[L][idx_m1].q == 0, f"central char of {L} not rational"
    print("    central character chi(-1) is rational on EVERY irrep "
          "(Galois-fixed).")
    # therefore Galois preserves the integer/spinorial split
    spinor = {L for L in labels if T["block"][L] == "spinorial"}
    for (A, B) in galois_pairs:
        assert (A in spinor) == (B in spinor), "Galois crossed the split!"
    img = set()
    for L in spinor:
        img.add(next(M for M in labels if chi[M] == [v.conj5() for v in chi[L]]))
    assert img == spinor
    print(f"    => Galois preserves the integer/spinorial split; the spinorial")
    print(f"       set {sorted(spinor)} maps to itself (2<->2', 4' & 6 fixed).")
    print("    [ok] a spinor's Galois conjugate is always a spinor.\n")


def stage1b_spectrum(T):
    print("=" * 70)
    print("(1b) THE SPINORIAL SPECTRUM (twisted Molien series)")
    chi, labels, sizes = T["chi"], T["labels"], T["sizes"]
    ws, dims, block = T["ws"], T["dims"], T["block"]

    # chi_{Sym^k V2}(class) for k=0..KMAX, per class, via the U-recurrence.
    nC = len(ws)
    symchar = [[None] * (KMAX + 1) for _ in range(nC)]
    for j, w in enumerate(ws):
        symchar[j][0] = Q5(1)
        if KMAX >= 1:
            symchar[j][1] = Q5(2) * w
        for k in range(2, KMAX + 1):
            symchar[j][k] = Q5(2) * w * symchar[j][k - 1] - symchar[j][k - 2]

    # m_{rho,k} = (1/120) sum_classes size * chi_rho(class) * chi_{Sym^k}(class)
    m = {L: [] for L in labels}
    for L in labels:
        for k in range(KMAX + 1):
            tot = Q5(0)
            for j in range(nC):
                tot = tot + Q5(sizes[j]) * chi[L][j] * symchar[j][k]
            val = tot * Q5(F(1, 120))
            assert val.q == 0 and val.p >= 0 and val.p == int(val.p), \
                f"bad multiplicity m[{L}][{k}] = {val}"
            m[L].append(int(val.p))

    # global check: sum_rho dim(rho) * m_{rho,k} = dim Sym^k = k+1
    for k in range(KMAX + 1):
        assert sum(dims[L] * m[L][k] for L in labels) == k + 1
    print("    [ok] global check: sum_rho dim(rho)*m_{rho,k} = k+1 for all k.")

    # selection rule: integer irreps only at EVEN k, spinorial only at ODD k
    for L in labels:
        for k in range(KMAX + 1):
            if m[L][k] != 0:
                parity_ok = (k % 2 == 0) if block[L] == "integer" else (k % 2 == 1)
                assert parity_ok, f"{L} appears at wrong-parity k={k}"
    print("    [ok] selection: INTEGER irreps live at EVEN k, SPINORIAL at ODD k")
    print("         (central-character parity = the spin-lemma bit).")

    # validate trivial row against Paper 174
    killed = {2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 28, 34, 38, 46, 58}
    for k in killed:
        assert m["1"][k] == 0
    for k in (0, 12, 20, 24, 30, 32, 42):
        assert m["1"][k] == 1
    assert m["1"][60] == 2
    for k in range(1, KMAX + 1, 2):
        assert m["1"][k] == 0
    print("    [ok] trivial row reproduces Paper 174 exactly: m12=m20=m30=1,")
    print("         m24=m42=1, m60=2, the 15 killed evens, all odd k = 0.")
    print("         (174 Molien (1-t^60)/((1-t^12)(1-t^20)(1-t^30)) confirmed.)")

    # closed form: M_rho(t) = N_rho(t)/((1-t^12)(1-t^20)); N truncates (free
    # module over C[f12,f20]).  N coeff: n_k = m_k - m_{k-12} - m_{k-20} + m_{k-32}
    def numerator(seq):
        n = []
        for k in range(KMAX + 1):
            v = seq[k]
            if k >= 12:
                v -= seq[k - 12]
            if k >= 20:
                v -= seq[k - 20]
            if k >= 32:
                v += seq[k - 32]
            n.append(v)
        return n

    def polystr(n):
        terms = []
        for k, c in enumerate(n):
            if c == 0:
                continue
            mono = "1" if k == 0 else (f"t^{k}")
            terms.append(("" if c == 1 and k == 0 else
                          (f"{c}*" if c != 1 else "")) + (mono if k else f"{c}"))
        return " + ".join(terms) if terms else "0"

    print("\n    closed-form generating functions  M_rho(t) = N_rho(t) / "
          "((1-t^12)(1-t^20)):")
    cutoff = 80  # numerators must vanish well before KMAX if denominator right
    for L in labels:
        n = numerator(m[L])
        assert all(c == 0 for c in n[cutoff:]), \
            f"numerator for {L} did not truncate -- denominator wrong"
        deg = max((k for k, c in enumerate(n) if c != 0), default=0)
        tag = block[L].upper()
        print(f"      {L:<3} [{tag:<9}] N = {polystr(n)}   (deg {deg})")

    # the "where the actors live" map: first appearance level of each irrep
    print("\n    WHERE THE ACTORS LIVE -- first level each irrep appears:")
    for L in labels:
        first = next(k for k in range(KMAX + 1) if m[L][k] != 0)
        print(f"      {L:<3} [{block[L]:<9}] first at k = {first:<3} "
              f"(spin j = {first}/2);  levels<=29: "
              f"{[k for k in range(30) if m[L][k]]}")

    print("\n    THE SPINORIAL SECTOR (the actors), first levels:")
    for L in ("2", "2'", "4'", "6"):
        lv = [(k, m[L][k]) for k in range(KMAX + 1) if m[L][k]][:8]
        print(f"      {L:<3}: " + ", ".join(f"k={k}(x{c})" for k, c in lv))
    return m


def main():
    T = build_table()
    stage1a_galois(T)
    stage1b_spectrum(T)
    print("\n" + "=" * 70)
    print("STAGE 1 COMPLETE: Galois structure clean; the odd/spinorial spectrum")
    print("174 skipped is computed and validated against 174's bosonic series.")
    print("The actors live at ODD k in {2, 2', 4', 6}. Downstream: Stages 2-3.")
    print("=" * 70)


if __name__ == "__main__":
    main()
