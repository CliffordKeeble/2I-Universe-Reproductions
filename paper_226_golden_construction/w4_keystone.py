"""
w4_keystone.py — Part D / W4 (the keystone cancellation), verification + charge sectors.

Machine-verifies the cold derivation committed in W4_derivation.md (W4-i), then diffs
against the adjudication memo's stated values (W4-ii), then checks the J-charge sectors
under C0 and tau (W4b). Exact over Q(sqrt5, i) with numeric confirmation at both real
places (sqrt5 -> +-sqrt5); the both-places check is the arbiter where sympy simplify is
unreliable over the golden field.

Conventions (from golden_algebra, transcribed from Paper 203 sec4):
  Gamma_seed=[[0,1],[r5,0]], Gamma_adj=[[0,-1],[r5,0]], Z=diag(1,-1)=sigma_3.
  sigma: r5->-r5 (i fixed);  tau: i->-i (r5 fixed), entrywise.
  Compact generator J := i*Z = diag(i,-i), J^2=-I  (NOT golden_algebra.J2).
  On operators X:  C0(X) = Z sigma(X) Z^-1  (Z^-1=Z).
  On states   v:  C0(v) = Z sigma(v)        (the induced state map; run #1.6 C0=Z o sigma).

Memo diff targets are hard-coded ONLY in the W4-ii diff step, never used to derive.
"""
import csv
import sympy as sp
from sympy import I, Matrix, eye
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq, canon

ROWS = []
J = I * Z2                       # compact generator i*Z = diag(i,-i)


def both_places_eq(A, B):
    """Exact equality AND numeric agreement at both real places sqrt5 -> +-sqrt5."""
    exact = mat_eq(A, B)
    D = A - B
    ok = True
    for s in (sp.sqrt(5), -sp.sqrt(5)):
        for e in D:
            if abs(complex(sp.N(e.subs(r5, s), 30))) > 1e-25:
                ok = False
    return exact, ok


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def C0_op(X):
    """C0 on an operator: Z sigma(X) Z^-1, Z^-1 = Z."""
    return mcanon(Z2 * sigma(X) * Z2)


def C0_state(v):
    """C0 on a state vector: Z sigma(v)."""
    return mcanon(Z2 * sigma(v))


def proportional(u, w):
    """Return the scalar c with u = c*w if u is a scalar multiple of w (w != 0), else None."""
    idx = next((i for i in range(len(w)) if canon(w[i]) != 0), None)
    if idx is None:
        return None
    c = canon(u[idx] / w[idx])
    return c if mat_eq(mcanon(u - c * w), sp.zeros(*u.shape)) else None


# ============================ W4-ii : verify + diff ============================
def w4_ii():
    # 1. sigma on the golden gammas (ground-truth signs)
    s_seed, s_adj = sigma(G_seed), sigma(G_adj)
    ex1, np1 = both_places_eq(s_seed, -G_adj)
    ex2, np2 = both_places_eq(s_adj, -G_seed)
    record("W4-ii-1 sigma(Gamma) signs", "sigma(Gs)=-Ga, sigma(Ga)=-Gs",
           ex1 and np1 and ex2 and np2,
           "sigma(Gamma_seed)=-Gamma_adj and sigma(Gamma_adj)=-Gamma_seed "
           "(both signs minus, exact+both places)")

    # 2. bivector B = Gseed*Gadj = sqrt5 Z
    B = mcanon(G_seed * G_adj)
    exB, npB = both_places_eq(B, r5 * Z2)
    record("W4-ii-2 bivector B=Gs*Ga", "B = sqrt5 Z",
           exB and npB, f"B = Gamma_seed*Gamma_adj = sqrt5*Z (= {B.tolist()})")

    # 3. C0 on the bivector: expect flip
    C0B = C0_op(B)
    exf, npf = both_places_eq(C0B, -r5 * Z2)
    record("W4-ii-3 C0(sqrt5 Z) flips", "C0(sqrt5 Z) = -sqrt5 Z",
           exf and npf, f"C0(B) = -sqrt5*Z (bivector flips); C0(B)={C0B.tolist()}")

    # 4. C0 on the compact generator J=iZ: expect fixed
    C0J = C0_op(J)
    exj, npj = both_places_eq(C0J, J)
    record("W4-ii-4 C0(J) fixed", "C0(J) = +J",
           exj and npj, f"C0(iZ) = +iZ = +J (compact generator fixed); C0(J)={C0J.tolist()}")

    # 5. DIFF vs adjudication memo sec2 (values quoted here ONLY for the diff, not the derivation)
    memo_bivector = mat_eq(C0B, -r5 * Z2)      # memo: C0(sqrt5 Z) = -sqrt5 Z
    memo_generator = mat_eq(C0J, J)            # memo: C0(J) = +J
    agree = memo_bivector and memo_generator
    record("W4-ii-5 diff vs memo sec2", "cold derivation agrees with memo (flip + fixed)",
           agree,
           "cold derivation MATCHES memo: C0(sqrt5 Z)=-sqrt5 Z (bivector flips) AND "
           "C0(J)=+J (generator fixed) — the keystone holds, fork stays closed" if agree
           else "SIGN MISMATCH vs memo — highest priority; the fork reopens",
           "" if agree else "cold vs memo disagree: bivector_match=%s generator_match=%s"
           % (memo_bivector, memo_generator))


# ============================ W4b : charge sectors ============================
def w4b():
    e1, e2 = Matrix([1, 0]), Matrix([0, 1])    # J=diag(i,-i): eigenvalues +i, -i
    # sanity: e1,e2 are the J-eigenvectors with eigenvalues +i,-i
    ev1 = mat_eq(mcanon(J * e1), mcanon(I * e1))
    ev2 = mat_eq(mcanon(J * e2), mcanon(-I * e2))

    # (a) C0 preserves each J-charge sector
    c1 = proportional(C0_state(e1), e1)
    c2 = proportional(C0_state(e2), e2)
    preserves = (c1 is not None) and (c2 is not None) and ev1 and ev2
    record("W4b-a C0 preserves J-sectors", "C0(e1) prop e1, C0(e2) prop e2",
           preserves,
           f"C0(e1)={c1}*e1, C0(e2)={c2}*e2 — C0 keeps each compact-charge sector "
           f"(does NOT flip the compact charge), consistent with C0(J)=+J")

    # (b) tau flips the compact charge at the OPERATOR level: tau(J) = -J
    tJ = tau(J)
    ext, npt = both_places_eq(tJ, -J)
    record("W4b-b tau(J) = -J", "tau(J) = -J (=> compact coupling q -> -q)",
           ext and npt,
           "tau(J) = -J: entrywise conjugation sends q A.J -> -q A.J, i.e. the compact "
           "charge q -> -q. This is 'tau flips the compact charge' (sharpened at solution "
           "level in W3)")

    # (c) HONEST check of the eigenvector-gloss: does tau EXCHANGE e1<->e2 as vectors?
    te1, te2 = mcanon(tau(e1)), mcanon(tau(e2))
    exchanges = (proportional(te1, e2) is not None) and (proportional(te2, e1) is not None)
    fixes = mat_eq(te1, e1) and mat_eq(te2, e2)
    # Report the ground truth: for J=iZ the eigenvectors are REAL, tau fixes them; the
    # "exchange" is the charge-level statement (b), not a permutation of these eigenvectors.
    record("W4b-c eigenvector-gloss, precisely", "clarify sense of 'tau exchanges sectors'",
           fixes and not exchanges,
           "tau FIXES the J=iZ eigenvectors as vectors (they are real): tau(e1)=e1, "
           "tau(e2)=e2 — it does NOT permute them. The 'exchange' in the brief is the "
           "charge-level fact tau(J)=-J (b): +i-eigenspace of J = -i-eigenspace of tau(J). "
           "Verified, gloss clarified — not a permutation of eigenvectors.",
           "clarification, not a failure: substance of W4b (C0 keeps compact charge, tau "
           "flips it) holds; only the parenthetical eigenspace wording is imprecise for J=iZ")


def main():
    print("=" * 72)
    print("W4 — keystone cancellation (verify + diff + charge sectors), exact over Q(sqrt5,i)")
    print("=" * 72)
    w4_ii()
    w4b()
    with open("results_w4.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for row in ROWS if row[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
