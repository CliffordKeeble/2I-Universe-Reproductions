# Findings — F7: the reality sweep over the aligner family. **The no-go is RETRACTED.**

**Paper 226 · Golden Construction · F7**
**11 July 2026 · exact symbolic over K(i) = ℚ(√5, i) · `bench_f7.py` · `results_f7.csv` · `run_f7.log`**
**Pre-registration:** `PRE-REGISTRATION_f7.md` (committed `0a279a5`, before this verdict).
**Brief of record:** `mr_code_brief_F7_reality_sweep.md` (Fizz, 11 Jul 2026).

Status tags: **DERIVED** (exact proof), **STRUCTURAL** (mechanism), **OBSERVED** (numeric). Nothing banks until the cold W-108 pass.

---

## Headline (Pattern 39 — reported flat; this retracts Mr Code's own F1–F4 result)

> **F7 = K-COLLAPSE. The F1–F4 COLLISION no-go is RETRACTED. K1 un-fires.**
>
> Sweeping the aligner family A = [[0, b], [c, 0]] symbolically, the reality-variety is exactly **b = c**. The representatives **A = K₂ = [[0,1],[1,0]]** (b=c=1) and **A = √5·K₂** yield a **real, nondegenerate (symplectic rank 4, det 125 / 3125), EOM-preserving** first-order action at the healthy s = 1 — kinetic **and** mass, under the τ constraint, at phase θ = π/2 (factor i). The F1–F4 COLLISION was **presentation-specific to A₀ = diag(√5, 1)** all along; with an aligning pairing the free split σ-theory has a healthy first-order action. Per §7 I **stop here and do not re-pose F4** — the retraction is the result; Cliff/Tor decide the next move.
>
> **Two load-bearing caveats, both reported straight (below): (1)** the surviving pairing K₂ **fails** 226 §5.2's frontier criterion — reality and the frontier **disagree**, and that disagreement is a genuine open structural question, not adjudicated here. **(2)** the "real action exists" criterion is the K1 question the brief posed; the *physical* sensibility of the antidiagonal (place-swap) pairing — positivity, the ± frequency split — is the re-posed F4 question, deferred per §7.

---

## §5 — Fizz's warm lean: algebra correct, but the lean is WRONG (reported straight)

- **Algebra verified (DERIVED):** A = [[0,√5],[1,0]] gives AΓ_s = diag(5,1), AΓ_a = diag(5,−1), both symmetric and rational. Fizz's pencil holds.
- **But the candidate does NOT give a real action (DERIVED).** With b = √5 ≠ c = 1, the mass term is not real at any phase — Fizz's candidate **fails F7a**. Her preferred aligner passes her preferred *kinetic-parity* test and fails the *actual reality* test. She called this exactly ("watch me on the σ-odd b; I want it too much").
- **⚠️ HALO FLAG — DISCONFIRMED, logged, not banked.** The working variety is **b = c**, which contains **σ-EVEN** members: **b = c = 1 (K₂, rational)** gives a real nondegenerate action. **A σ-odd b is NOT required.** The "the irrational balances / odd pays" rhyme with the σ-odd mass parity is **not a fact of this pairing** — the data actively refutes it. It must not appear in any paper, summary, or handoff as a finding or a unity. (Recorded here solely to close the flag.)

## §3 F7a — the reality sweep (THE CRUX, DERIVED)

For the general aligner A = [[0, b], [c, 0]], ψ̄ = (σψ)ᵀA, L = ψ̄(Γ_s∂₀ + Γ_a∂₁ + msZ)ψ, reality criterion **Im[ev₁(e^{iθ}L)] mod TD = 0** (kinetic mod TD + mass pointwise), swept over {τ, σ-real, C₀} × the phase circle:

- **Solved symbolically for (b, c):** the reality conditions force **br = cr and bi = ci, i.e. b = c** (the sympy solve returns exactly this variety; the phase θ is then fixed by arg(b), θ = π/2 for real b). The aligner must be a **scalar multiple of K₂**: A = b·K₂.
- **Mechanism (STRUCTURAL).** For any aligner, AΓ_s = diag(b√5, c) and AΓ_a = diag(b√5, −c) are **diagonal** → the kinetic term is **same-component** (each σψᵢ pairs with ψᵢ), unlike A₀'s cross/antisymmetric structure. Same-component ⟹ Re[σψᵢ∂ψᵢ] = ½∂|ψᵢ|² is a total derivative and Im is a current — the ordinary-Dirac situation, made real by the **factor i** (θ = π/2). The mass term AZ = [[0,−b],[c,0]] is cross-component B = −b·σψ₁ψ₂ + c·σψ₂ψ₁; under τ it is **purely imaginary iff b = c** (then B = −2i·Im(ψ₁*ψ₂)), so i·L_mass = 2m·Im(ψ₁*ψ₂) is real and non-trivial. **b = c reales kinetic and mass at the same phase.** This is precisely the parity divorce of F3 *healed* by the aligner.
- **Survivors (DERIVED):**

| aligner | b, c | real? (τ, s=1) | symplectic | EOM keeps τ | class |
|---|---|---|---|---|---|
| **K₂ = [[0,1],[1,0]]** | b=c=1 (σ-even) | ✅ at θ=π/2 | rank 4, det 125 | ✅ | **K-COLLAPSE** |
| **√5·K₂** | b=c=√5 (σ-odd) | ✅ at θ=π/2 | rank 4, det 3125 | ✅ | **K-COLLAPSE** |
| Fizz [[0,√5],[1,0]] | b=√5, c=1 | ❌ | — | — | fails |

- **Nondegeneracy (DERIVED):** the symplectic form (∂₀-antisymmetric part of the *real* Lagrangian Re[iL]) has **rank 4** for the survivors — genuine first-order dynamics, no K1 degeneracy. (det values 125 = 5³, 3125 = 5⁵ are basis-dependent; rank 4 is the invariant.)
- **EOM-preservation (DERIVED):** the golden Dirac EOM is A-independent (δS/δψ̄ = (Γ∂+msZ)ψ = 0) and, as in F3 sub-check 2, keeps ψᵢ = fᵢ + i√5 gᵢ real and closed. The reality constraint is dynamically consistent.

⟹ **K-COLLAPSE, not K-INDETERMINATE:** the survivors are real **and** nondegenerate **and** EOM-consistent.

## §4 F7b — the frontier, side by side (THE NAMED TENSION — reported, not adjudicated)

226 §5.2's displayed criterion **𝒞·σ(Γ^μ)·𝒞⁻¹ = ε(Γ^μ)ᵀ** with a **uniform scalar ε**, tested with 𝒞 ∈ {A, σ(A), A⁻¹}:

- **K₂ FAILS the frontier** (all three conventions): even within Γ_s the entry-ratios are {−1, +1} — no uniform ε. σ(K₂) = K₂, K₂⁻¹ = K₂ (rational), same verdict.
- Fizz's candidate also fails (ratios {√5, −√5/5}, as she reported).

> **The disagreement is the finding.** The **reality** criterion (physically mandated, the one that fires K1) says K₂ gives a healthy action → **no-go retracted**. The **frontier** criterion (a condition on the *conjugation matrix* 𝒞, a different object) says K₂ is not an admissible σ-conjugation → an obstruction persists, but a **different** one. **These two criteria disagree about σ for the golden pair.** This goes to Mr A and Tor as an open structural question — *is the K1-relevant object the action's pairing form A (K-COLLAPSE) or the conjugation matrix 𝒞 (frontier still obstructs)?* **Not adjudicated here.**

## §6 — σ(M²) equal-mass under the σ-odd mass axiom (DERIVED, EXPECT met)

With m² σ-odd (m² = √5·μ², Cliff's antiparticle axiom), **M² = m²s²/√5 = μ²s² is σ-even** — identical at ∞₁ and ∞₂ — and **E² = k² + M² is σ-even** (the nested-radical place-reading Tor declined to read warm, read here exact). **The equal-mass antiparticle axiom is delivered by the algebra, not merely asserted.** (s = i recorded: M² = −μ², still σ-even.)

## The B lemma (owed regardless — written down)

**Claim.** Let X₀ := Im[ev₁ L] mod TD and X₁ := Re[ev₁ L] mod TD. Then the reality defect of the phase-rotated Lagrangian is **defect(e^{iθ}L) = cos θ · X₀ + sin θ · X₁** (ℝ-linearity of Im in the real/imaginary split, since e^{iθ}L = (cos θ + i sin θ)L). **Proof of the circle from two points:** X₀ and X₁ are bilinear forms; if they are **direction-separated** (X₀ supported on ∂₀-currents, X₁ on ∂₁-currents — verified for the 191 pair: `findings_f1_f4.md` gives A₀'s X₀ pure-∂₀, X₁ pure-∂₁) and independent mod TD (§3b D-lemma: distinct derivative directions cannot cancel, as first-order bilinear total derivatives are exactly the per-direction symmetric parts), then cos θ·X₀ + sin θ·X₁ = 0 as a form requires **cos θ = 0 and sin θ = 0 separately — no θ.** Two evaluations (θ = 0 → X₀; θ = π/2 → X₁) therefore decide the entire circle. **Dressing rigidity (Mr A):** the σ-type dressing is rigid up to sign (σ fixes i; σ(D)D = 1 pins the phase to ±1), so the continuous dressing family collapses to this θ-circle. — For the *aligner* K₂, X₁ = 0 identically (the diagonal AΓ makes Re a total derivative), so the circle is real at θ = π/2; the lemma's "no θ" branch was the A₀ story, now superseded.

## Genericity (STRUCTURAL)

The reality-variety **b = c** and the working representative **K₂ (rational)** carry **no √5** — √5 is load-bearing nowhere in the collapse. The result is real-quadratic-generic (the whole construction is field-level; K₂ works over any ℚ(√d, i)). **No d-dependence, no cloud.** (The only √5 in sight is the σ-odd member √5·K₂, which is an *alternative* representative, not a requirement — see the disconfirmed halo flag.)

---

## What this means / what's owed

- **`findings_f1_f4.md` gets its second correction banner:** the COLLISION no-go is **RETRACTED** — it was specific to the pairing A₀ = diag(√5,1); an aligning pairing (K₂) gives a real nondegenerate EOM-consistent action, so **K1 does not fire** for the theory.
- **Stopped at F7 (§7 stop-condition).** F4 (the ± frequency signature / K4 / unitarity on the surviving pairing) is **not** re-posed — that is Cliff/Tor's call, and it is where the antidiagonal pairing's physical sensibility (positivity, place-swap structure) must be vetted.
- **Open question handed up (the §4 tension):** reality and the 226 §5.2 frontier disagree about K₂. Whether the free split σ-theory is *fully* healthy or merely *reality-healthy-but-frontier-obstructed* is unresolved and belongs with Mr A / Tor.
- **Halo flag closed:** σ-odd b is **not** required (K₂ is σ-even). No "odd pays" unity enters anywhere.

## Reproduce

```
cd paper_226_golden_construction
python bench_f7.py   # -> results_f7.csv; §5 PASS, F7a K2/√5K2 real+nondeg+EOM, Fizz fails, §6 PASS; VERDICT K-COLLAPSE
```
Environment: `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic; no RNG; floats only in §6's two-place read.

🐕☕⬡
