# Findings — Run #1.5, the place-pair toy

**Spec:** Fizz reading `fizz_reading_bench_run_1_paper_226.md` §6 (pre-registered). Executor: Mr Code. Pre-reg: [PRE-REGISTRATION_run1_5.md](PRE-REGISTRATION_run1_5.md).
**Method:** exact 2×2 over ℚ(√5,i), every operator equality cross-checked numerically at **both** real places (mpmath) per the standing discipline promoted in Fizz §3. Script: [run_1_5_place_pair.py](run_1_5_place_pair.py); output: [results_run1_5.csv](results_run1_5.csv), [run_1_5.log](run_1_5.log).
**Result: 2 / 3 match EXPECT.** The two-place *picture* is partly validated (B1.5-1, B1.5-3 pass); the specific claim that **C₀ = Z⊗swap implements charge conjugation between the slots fails** (B1.5-2, the make-or-break).

## Summary

| Item | EXPECT | Outcome |
|---|---|---|
| B1.5-1 pair reads two dynamics | pass | **PASS** |
| B1.5-2 C₀ absorbs transposition | pass (make-or-break) | **FAIL** |
| B1.5-3 Born = Galois norm | pass | **PASS** |

## B1.5-1 — PASS (with a bonus)

Reading the one K-rational field at the two places gives a slot-pair solving (L_q at ∞₁, ev₂(M_q)=σ(M_q) at ∞₂). Identifying σ(M_q) in the operator family: it is the **swapped-order (Γ_adj on ∂₀), reversed-momentum (k→−k)** operator — i.e. the ∞₂ slot natively carries **reversed metric *and* reversed momentum**, which is exactly the antiparticle kinematics. Fizz's §1 synthesis in propagating form holds: the ∞₂ evaluation of the same data is a genuine conjugate-dynamics reading.

## B1.5-2 — FAIL (the make-or-break)

Tested at the place-pair level (kernel match ⟺ operator proportionality; no on-shell needed):
- The C₀=Z⊗swap image has slot-1 = Z·(∞₂ reading), slot-2 = Z·(∞₁ reading).
- **Image slot-1 does NOT solve the literal −q operator M_{−q}@∞₁** (`prop = None`, at either a-sign). It solves a **seed↔adjoint swapped** operator [swapped(Γ_adj@∂₀), ±q, ...].
- Image slot-2 solves a *literal*-order operator but momentum-reversed (k−) — so the image pair is **mixed-order** (swapped at slot 1, literal at slot 2), never the clean literal −q pair.

**Verdict (per Fizz §6's registered criterion):** the place structure does **not** absorb the A5 coordinate transposition. The two golden directions (Γ_seed²=+√5 timelike-like vs Γ_adj²=−√5 spacelike-like) are genuinely inequivalent, and C₀=Z⊗swap cannot reconcile them with −q dynamics at any slot assignment. **Not patched.**

**Constructive reading (for Fizz — data selects among your own §2 candidates):** this is *not* a refutation of the two-place synthesis (B1.5-1/B1.5-3 stand). It is a selection: **C₀=Z⊗swap alone is not the antiparticle map.** The working map is the composite A5 already identified — **C₀ ∘ (coordinate transposition x⁰↔x¹) ∘ (a_μ→−a_μ)** — which is precisely your §2 candidate "C₀ ∘ (coordinate transposition), the toy's closing symmetry." The transposition is physical (it swaps the timelike/spacelike golden directions), not bookkeeping, and the place structure does not remove it. Run #1.5 says: build the antiparticle map as that composite, not as Z⊗swap.

## B1.5-3 — PASS

Born pairing across slots P = ψ⁽¹⁾·ψ⁽²⁾ = a·σ(a) = p²−5q² = **N(a)**, the Galois norm (192 §5's DERIVED reduction), and it is C₀(slot-swap)-invariant; the slot inner product ⟨ψ⁽¹⁾|ψ⁽²⁾⟩ is σ-fixed (lands in ℚ). The place-level Born rule reproduces the golden norm exactly.

## Net

The place-pair reformulation **succeeds as an interpretation** (two slots = two signatures = particle/antiparticle-with-reversed-metric-and-momentum; Born rule = Galois norm) and **fails as a simplification of C** (Z⊗swap does not remove the A5 transposition). The v0.2 "place-pair formulation" spine is viable for states and Born structure, but the charge-conjugation operator must still carry the coordinate transposition — the metric-activity of σ is a real feature of the dynamics, not an artifact removable by relabelling slots.

## Reproducibility

`python run_1_5_place_pair.py` → `results_run1_5.csv`, `run_1_5.log`. Exact symbolic + numeric-at-both-places. Environment: sympy 1.14.0, mpmath 1.3.0, numpy 2.3.4, Python 3.13.
