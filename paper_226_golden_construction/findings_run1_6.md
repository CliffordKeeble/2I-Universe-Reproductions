# Findings — Run #1.6 (C_phys, null-line lemma, factorisation frontier)

**Spec:** Fizz reading `fizz_reading_run_1_5_paper_226.md` §4 (pre-registered). Executor: Mr Code. Pre-reg: [PRE-REGISTRATION_run1_6.md](PRE-REGISTRATION_run1_6.md).
**Method:** exact 2×2 / 4×4 over ℚ(√5,i), operator equalities cross-checked numerically at both real places (standing discipline). Script: [run_1_6.py](run_1_6.py); output: [results_run1_6.csv](results_run1_6.csv), [run_1_6.log](run_1_6.log).
**Result: 4 / 4 match EXPECT.** Fizz's §3 reframe survives its own falsifier (B1.6-4).

## Summary

| Item | EXPECT | Outcome |
|---|---|---|
| B1.6-1 C_phys maps L_q → L_{−q} | pass | **PASS** |
| B1.6-2 C_phys² = identity | pass | **PASS** |
| B1.6-3 null-line lemma | pass | **PASS** |
| B1.6-4 factorisation frontier | pass | **PASS** |

## B1.6-1 — the composite is a clean charge conjugation

`C_phys = C₀ ∘ (x⁰↔x¹) ∘ (a→−a)` sends a plane-wave solution `(χ; k₀,k₁,a₀,a₁; q)` of L_q to `(Zσχ; k₁,k₀,−a₁,−a₀; −q)`, and
> **M(k₁,k₀,−a₁,−a₀,−q) = Z·σ(M_q)·Z exactly** (verified symbolically and numerically at both real places).

So `C_phys(χ)=Zσ(χ)` solves the **literal** −q operator with **no residual** — the composite that A5 identified as the unique all-k intertwiner is confirmed end-to-end. This is the working antiparticle map that Z⊗swap (run #1.5) failed to be.

## B1.6-2 — C_phys is an involution

Applying C_phys twice returns the amplitude (Zσ(Zσχ)=χ), the parameters, and the charge exactly. The three factors each square to 1 (transposition²=1, (a→−a)²=1, C₀²=1) and compose to the identity. C_phys² = 𝟙 on solutions.

## B1.6-3 — the transposition is a null-line reflection

In lightcone coordinates x± = x⁰±x¹, the swap x⁰↔x¹ **fixes x⁺ and negates x⁻** — reflection in the null line x⁻=0. The lightcone gammas Γ± = Γ_seed ± Γ_adj are exactly the **null-ray nilpotents** ([[0,2],[0,0]] and [[0,0],[2√5,0]], as in analysis §1: Γ±²=0), and the operator is L = Γ₊∂₊ + Γ₋∂₋, so the reflection ∂₋→−∂₋ is chirality-selective (one null mover reflected, the other fixed). This is the geometric content of Fizz §2's "charge-flip ∘ single reflection = CRT-shaped."

## B1.6-4 — the factorisation frontier (the falsifier that could have killed the reframe)

Searching for a charge conjugation C mapping L_q → L_{−q} as a constant-matrix `lin`ear map, a `tau`-borne map B∘τ (which flips the plane-wave momentum k→−k), or a `sigma`-borne map B∘σ:

| coupling | linear | τ-borne | σ-borne |
|---|---|---|---|
| **σ-coupling** (golden, K skeleton; D=∂−q·a·Z) | 0 | 0 | 0 |
| **i-coupling** (T_τ / honest U(1), K(i); D=∂−i·q·a) | 0 | **2** | 0 |

- **The K skeleton admits no factorized charge conjugation at all** — linear, τ-borne, and even σ-borne alone all give dim 0. This reproduces run #1's theorem (the golden coupling flips only through the full composite C_phys, never a single-factor C).
- **The extension revives it:** the i-coupling admits a 2-parameter family of **τ-borne** charge conjugations (standard complex-conjugation charge conjugation). This is Fizz §3's "the extension buys the factors," now DERIVED in the toy.

**EXPECT met** (Fizz's falsifier: "a fail would mean the factorisation claim is wrong and §3 overreaches"). §3 does not overreach.

**Deviation (reported, transparent):** the pre-registered B1.6-4 sketch tested the τ-borne intertwiner at fixed momentum. That is a setup error — a τ-borne C sends e^{−ikx}↦e^{+ikx}, so the target operator must be evaluated at **−k**. The first run FAILed on that ill-posed sub-test; the corrected test (target at −k, and framed as the σ-vs-i contrast the spec intends) is the faithful one and passes. The correction could have left the frontier dirty (it did not).

## Net — the reframe, benched

Both bench theorems of Fizz §3 are now machine-confirmed in the toy:
1. **No standalone C at the golden skeleton** (B1.6-4, σ-coupling column all-zero; run #1 theorem reproduced).
2. **The only conjugation is the σ-borne composite C_phys**, which is a genuine involutive charge conjugation (B1.6-1/2) whose coordinate factor is a **null-line reflection** (B1.6-3) — the CRT shape.
3. **The extension buys the factors** (B1.6-4, i-coupling admits a τ-borne C).

So at the golden skeleton, C/P/T do not exist separately — only the CRT-shaped composite is ℚ(√5)-rational; the factorisation is purchased by the closure K→K(i). The programme's old "CPT-broken" claim has no K-rational subject to break; the well-posed post-closure map carries a reflection, and reflections mirror — consistent with the BESIII / 1958-positron pins. Routed to the v0.2 input pile (still gated on the cold construction run).

## Reproducibility

`python run_1_6.py` → `results_run1_6.csv`, `run_1_6.log`. Exact symbolic + numeric-at-both-places. Environment: sympy 1.14.0, mpmath 1.3.0, numpy 2.3.4, Python 3.13.
