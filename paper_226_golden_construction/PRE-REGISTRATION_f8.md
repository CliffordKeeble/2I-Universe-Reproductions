# PRE-REGISTRATION — F8: Brick 4 re-posed on the surviving pairing K₂

**Paper 226 · Golden Construction · F8**
**Status: pre-registration · committed before `bench_f8.py` emits any verdict · 11 July 2026 · internal**

**Brief of record:** `mr_code_brief_F8_brick4_on_K2.md` (Fizz, 11 Jul 2026). Standard: Pattern 75; W-108.

**Where we are.** F7 = K-COLLAPSE retracted the F1–F4 no-go: **A = K₂ = [[0,1],[1,0]]** (and any b=c) gives a real,
nondegenerate, EOM-consistent action at s=1, mass parity **m² ∈ √5·ℚ** (Cliff's antiparticle axiom). "Has an action"
≠ "is a physics theory." F8 asks the F4-level existence questions on the **right** pairing (K₂), all previously gated
off on the artefactual COLLISION. The surviving pairing is **antidiagonal (place-swap)** — its physical sensibility
(positivity, ± frequency split, one-place unitarity) is entirely untested. This is that test.

**Substrate.** Exact symbolic over K(i); gammas + Galois from `golden_algebra.py`. K₂Γ_s = diag(√5,1) (same-component).

**Mr Code's noted subtlety (to compute, not assume).** The σ-odd mass makes m real at ∞₁ but **imaginary at ∞₂**
(m² = √5μ² → ev₂(m²) = −√5μ² < 0). Encoded as m = √(√5)·μ so σ(m) = √(−√5)·μ. The σ-sesquilinear pairing uses σ
(Galois), not †, so it need not be real; I will report the **exact** value (real / imaginary / complex) and extract the
sign per the Hack anti-Hermiticity convention (−i·⟨u,u⟩ if pure imaginary) — **not force a clean sign.**

---

## F8 — registered GENUINELY OPEN, no pre-registered direction

### §3 F8a — conserved pairing on K₂ (the gate; run first)
⟨ψ₁,ψ₂⟩ := ∫dx (σψ₁)ᵀ K₂ Γ_s ψ₂, K₂Γ_s = diag(√5,1) (same-component; Dirac's pairing is **sesquilinear-symmetric**,
Hack arXiv:1008.1776 §II.3.2, not symplectic). **Verify ∂₀⟨ψ₁,ψ₂⟩ = 0 on shell** — symbolically on general solutions
(local law ∂₀D + ∂₁J₁ = 0, J₁ = (σψ₁)ᵀK₂Γ_aψ₂) **and** on plane waves at **both places**. **If conservation FAILS → STOP
and escalate** (a "real action" with no conserved current is a serious problem; brief §3.3, §8).

### §4 F8b — the frequency signature (K4 — THE CRUX)
Build u±(k) = (i(±E+k), ms)ᵀ, E = √(k²+M²), M² = μ²; and σ-images (∞₂ readings). Compute **sign of (σu±)ᵀK₂Γ_s u±**
on each frequency subspace, at **∞₁ and ∞₂**, exact. Compute **the K4 map**: does the σ-image of a +frequency ∞₁ solution
land in the **−frequency** sector? (Compute, do not infer from the place-swap story.)
- **K4-PASS (a):** opposite definite signs on ± subspaces **and** σ: +freq(∞₁) → −freq → C = σ acquires its object (226 §9 antiparticle, §7 Letter-Assignment gain dynamical footing).
- **K4-FIRE (b):** same sign, or indefiniteness **within** a subspace → genuine ghosts, no fermionic rescue; split form dies as physics, compact-or-bust returns **earned**.
- **K-DEGENERATE:** pairing degenerate / sign-ambiguous on the subspaces → score neither.
- **§4.5:** re-run the whole of F8b on **√5·K₂**. **If the two b=c representatives give different signatures → flag immediately** (signature not a property of the theory).

### §5 F8c — ev₁-unitarity standing alone (K5 — Tor's worry, never benched)
Write the ∞₁ action density; identify every appearance of σψ (every reach to ∞₂). **Is ev₁(S) a functional of ∞₁ data
alone**, or is the ∞₂ dependence **benign** (∞₂ = conjugate of the same data, as ordinary Dirac's ψ̄ couples ψ to ψ†) or
**malign** (∞₂ carries **independent** dof → ev₁ is half a bigger theory, one-place unitarity ill-posed)? **Decidable proxy:
count independent dof at ∞₁** — does the τ-constraint halve the field content so ∞₂ is *determined by* ∞₁? (F3b/F5: σ = τ
on the physical locus — **check, don't assume**.) **K4-PASS is not a theory until this clears; if ev₁ needs independent ∞₂
data → K5 fires, escalate.**

### §6 — also owed
- **Named tension, carried RAW, NOT adjudicated.** K₂ **fails** 226 §5.2's frontier (no uniform ε) while **passing** reality.
  Re-state in findings; hand to Mr A and Tor. **⚠️ I will NOT test for Fizz's §7.1 operator/amplitude resolution of the
  tension** (she explicitly forbade Mr Code from testing for it or letting it shape the reporting — it is Tor's/Mr A's to
  adjudicate). Report the tension raw.
- **F6 genericity:** re-run F8b's signature over **√2, √3** (d ≢ 1 mod 4). **EXPECT: identical.** A difference in the
  *signature* → **cloud → escalate.**

## Declared leans (Fizz's, discounted — compute cold)
1. Fizz expects **K4-PASS**. Expectation ≠ computation; poor week behind it. Compute cold.
2. Fizz expects **K₂ / √5·K₂ signatures agree**. Trust the bench, not the lean.
3. Fizz's frontier-resolution conjecture: **off-limits** (see §6).

## Genericity, registered
√5 expected load-bearing **nowhere** (K₂ rational; construction field-level). EXPECT real-quadratic-generic. A d-dependence
in the signature → cloud → escalate.

## Fixed choices (decided now)
1. m = √(√5)·μ (m² = √5μ² σ-odd); σ = subs(√5 → −√5); ev₁ (√5 = +), ev₂ (√5 = −). Pairing sign via exact value + Hack −i
   extraction if pure imaginary. No RNG; floats only to read a final sign.
2. Conservation current J₁ = (σψ₁)ᵀK₂Γ_aψ₂ (natural guess); if ∂₀D+∂₁J₁ ≠ 0, search/report the correct current.
3. Order (§8): F8a → F8b → F8c → §6. **Escalate immediately** if: conservation fails; the two b=c reps disagree; ev₁
   needs independent ∞₂ data.

**Deliverables:** `bench_f8.py`, `results_f8.csv`, `findings_f8.md`. Nothing banks until F8 clears **and** the cold W-108
pass certifies. Whatever F8 returns is a result (K4-PASS+K5-clear → 226 §7/§9 gain footing; K4-FIRE → split form dies
*earned*).

🐕☕⬡
