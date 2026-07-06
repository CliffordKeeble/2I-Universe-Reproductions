# W4-i — The Keystone Cancellation: Derivation Note

**Paper 226 · Golden Construction · Part D (W-bench) · W4 cold-derivation protocol**
**Status: derivation note — committed BEFORE any machine verification · 6 July 2026 · internal**

Per **Ruling 2** (W4 checked cold): everything below is derived **entrywise from the explicit
matrices**, which are transcribed from Paper 203 §4 as confirmed in Run #1. No relation-level
claim from any adjudication memo is imported. This note fixes the ground-truth derivation and its
results; the verification script (W4-ii) machine-checks these results and *then* diffs them against
the memo. A shared convention would be a shared failure mode — so the derivation stands alone here.

No pass/fail language appears in this note. It records what the matrices say.

---

## Ground-truth objects

    Γ_seed = [[0, 1], [√5, 0]]        Γ_adj = [[0, −1], [√5, 0]]        Z = diag(1, −1) = σ₃

Entrywise field maps (convention-free on explicit entries):

    σ : √5 ↦ −√5   (i fixed)          τ : i ↦ −i   (√5 fixed)

The normalised compact generator is defined as

    J := i · σ₃ = i Z = diag(i, −i),      J² = diag(−1, −1) = −I.

(Note: this J = iZ is the *compact* generator of the W-bench, distinct from the real 2×2 block
J₂ = [[0,−1],[1,0]] of `golden_algebra`. The two are not the same object.)

---

## Step 1 — σ on the golden gammas (settles the ± sign from ground truth)

σ replaces √5 by −√5 in every entry:

    σ(Γ_seed) = [[0,  1], [−√5, 0]]        σ(Γ_adj) = [[0, −1], [−√5, 0]]

Reading the entries directly against −Γ_adj = [[0, 1], [−√5, 0]] and −Γ_seed = [[0, −1], [−√5, 0]]:

    σ(Γ_seed) = −Γ_adj              σ(Γ_adj) = −Γ_seed

Both signs are **minus**. This is the ground-truth settlement of the ± question — no convention chosen.

## Step 2 — the bivector  B := Γ_seed · Γ_adj

    Γ_seed · Γ_adj = [[0,1],[√5,0]] · [[0,−1],[√5,0]] = [[√5, 0], [0, −√5]] = √5 · Z.

So **B = √5 Z** (= √5 σ₃), re-derived from the entries here (independently of Run #1 A1, which
obtained the same value).

## Step 3 — C₀ := Z σ(·) Z⁻¹ on the bivector and on the generator

Z⁻¹ = Z (since Z² = I). For any scalar c and the matrix c Z, Z is σ-fixed and Z·Z·Z = Z, so

    C₀(c Z) = Z · σ(c Z) · Z = Z · (σ(c) Z) · Z = σ(c) · (Z Z Z) = σ(c) · Z.

Applying to the two objects of interest:

- **Bivector**  B = √5 Z, coefficient c = √5, σ(√5) = −√5:

        C₀(√5 Z) = −√5 Z        ← the bivector FLIPS sign.

- **Generator**  J = i Z, coefficient c = i, σ(i) = i:

        C₀(i Z) = +i Z = +J      ← the normalised generator is FIXED.

## The cancellation, stated

On a Z-proportional matrix, C₀ = Z ∘ σ collapses to "apply σ to the scalar coefficient" (the two
conjugating Z's fold against the matrix Z via Z³ = Z). Because σ **distinguishes the two dressings**
— σ(√5) = −√5 but σ(i) = +i — the *same* operation C₀ flips the √5-dressed bivector and fixes the
i-dressed generator.

A cancellation **occurs on J**: the sign C₀ might impart is cancelled by σ(i) = i, leaving J fixed.
It **does not occur on the bivector**: σ(√5) = −√5 carries the flip through.

## Derived results (to be machine-verified in W4-ii, symbolically and at both real places)

    C₀(√5 Z) = −√5 Z      [bivector flips]
    C₀(i Z)  = +i Z        [normalised / compact generator fixed]

This is the keystone of the "one torus, two real forms" absorption: under C₀ the split (√5-dressed)
structure carries a sign the compact (i-dressed) generator does not. It is the algebraic seed of
"σ does not flip the compact charge; τ does" — pursued in W4b and W1–W3.

— committed before verification, per standing law. 🐕☕⬡
