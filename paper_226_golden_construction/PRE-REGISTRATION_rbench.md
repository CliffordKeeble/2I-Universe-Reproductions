# PRE-REGISTRATION — Part C, the reconciliation R-bench (R1–R6)

**Paper 226 · Golden Construction · Part C**
**Status: pre-registration · committed before rbench.py runs · 6 July 2026 · internal**

**Source.** These six items are transcribed verbatim from Fizz's working note
`fizz_reconciliation_seam_paper_226.md` §4 (in the user's Downloads, uncommitted per the standing
paper-draft rule). That note works the "one spine" reconciliation: the internal (cold) construction
is the even-sector shadow of the geometric (briefed) one under the embedding δ ↦ Γ_seed·Γ_adj =
√5·σ₃, so j ↦ Z. One σ, dressed with the a-flip, is charge conjugation on both sectors — clean on
the even/internal sector, carrying the null-line reflection compensator on the odd/geometric one.
The R-bench converts that reconciliation from warm to demonstrated. One claim per EXPECT.

- **R1:** δ ↦ √5σ₃ is a unital K-algebra embedding; j ↦ Z; e± ↦ (1 ± Z)/2 orthogonal idempotents.
  **EXPECT: pass.**
- **R2:** the transported gauge action e^{θj}·a equals diag(e^θ, e^{−θ}) on spinor components and
  reproduces Run #1 A4's weight verification. **EXPECT: pass.**
- **R3:** entrywise σ restricted to the embedded 𝕂 equals 𝕂's Galois σ (δ ↦ −δ). **EXPECT: pass.**
- **R4:** σ(Γ₊) = −Γ₊, σ(Γ₋) = +Γ₋; C₀(Γ₊) = +Γ₊, C₀(Γ₋) = −Γ₋. **EXPECT: pass.**
- **R5 (even-sector control):** scalar toy on the internal Lagrangian's equation of motion:
  C = σ ∘ (a-flip) maps q-solutions to (−q)-solutions **with no coordinate factor**.
  **EXPECT: pass — a fail breaks the one-spine dichotomy.**
- **R6 (the dichotomy):** geometric toy: σ ∘ (a-flip) **without** the reflection fails as a solution
  map (re-confirming Run #1.5's B1.5-2 in this framing), and succeeds on the even-projected
  observables. **EXPECT: pass on both halves.**

Per Ruling 1 (badge discipline): outcomes are reported against EXPECT only; no DERIVED/STRUCTURAL
language in the findings. Results fold into v0.2's reconciliation section after merge + Adversary.

🐕☕⬡
