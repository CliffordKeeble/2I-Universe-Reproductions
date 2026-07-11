# Findings — F8: Brick 4 on K₂. **§3 gate FAILS → ESCALATION.**

**Paper 226 · Golden Construction · F8**
**11 July 2026 · exact symbolic over K(i) = ℚ(√5, i) · `bench_f8.py` · `results_f8.csv` · `run_f8.log`**
**Pre-registration:** `PRE-REGISTRATION_f8.md` (committed `4902a01`, before this verdict).
**Brief of record:** `mr_code_brief_F8_brick4_on_K2.md` (Fizz, 11 Jul 2026).

Status tags: **DERIVED** (exact proof), **STRUCTURAL** (mechanism). Nothing banks until the cold W-108 pass.

---

## Headline (Pattern 39 — reported flat; escalated per §8, F4 not banked)

> **F8 = ESCALATE. The §3 gate fails: the *symmetric* conserved pairing that K4 asks the sign of does not exist for the golden field.** Two **independent** structural facts, both DERIVED:
>
> **(A) The construction, not just the axiom.** ⟨ψ₁,ψ₂⟩ = ∫(σψ₁)ᵀK₂Γ_sψ₂ is **not conserved** — for K₂ *and* A₀, for σ-even *and* σ-odd mass, on general fields *and* **on the physical τ-locus with an in-field σ-even mass** (residual `2μ(−a₂a₃ − √5·i·a₂b₃ + √5·i·a₃b₂ − 5b₂b₃) ≠ 0`). Reason: the σ-adjoint **fixes the phase** (σ fixes i, so ψ̄ = (σψ)ᵀA → e^{+iα}ψ̄, not e^{−iα}), and the τ-reality constraint is a **Majorana (σ-real) condition** (e^{iα}ψ leaves the τ-locus). A σ-real field carries **no conserved symmetric U(1) Dirac charge.**
>
> **(B) The axiom, separately.** Cliff's antiparticle axiom (M² σ-even) + the golden dispersion (M² = m²s²/√5) force the first-order mass coefficient **m = 5^{1/4}μ ∉ K** (minimal polynomial x⁴−5, degree 4). The massive golden Dirac operator **leaves the field K.**
>
> **What IS conserved** is the *antisymmetric* symplectic form **Ω** of F7 (rank 4, det 125 — the "Grassmann-odd shadow" the brief itself names) — **not** the symmetric pairing K4 needs. So the theory is a **symplectic / Majorana-like** object, not a charged Dirac theory. **K4 cannot be posed as stated; hand back to Cliff/Tor/Fizz/Mr A.**

Fizz pre-registered a lean toward **K4-PASS**. The cold bench finds the *premise* of K4 — a conserved symmetric charge whose sign splits by frequency — **does not exist**. That is the "compute it cold, do not let me tell you what they'll be" outcome.

---

## §3 F8a — conservation (the gate): FAILS (DERIVED)

Density D = (σψ₁)ᵀK₂Γ_sψ₂ (K₂Γ_s = diag(√5,1), same-component; sesquilinear-**symmetric**, Hack §II.3.2). I solved for the spatial current B (from the derivative terms), leaving a **B-independent mass residual** — the mass term in ∂₀ψ that no current can absorb. It is nonzero in every case:

| pairing | mass | conserved? | mass residual |
|---|---|---|---|
| K₂ | σ-odd (axiom) | ❌ | `5^{1/4}μ(c₀c₃(1−i) − c₁c₂(1+i))` |
| K₂ | σ-even | ❌ | `−2μ c₁c₂` |
| A₀ = diag(√5,1) | σ-odd | ❌ | `−5^{3/4}μ c₀c₂(1+i) + 5^{1/4}μ c₁c₃(1−i)` |
| A₀ | σ-even | ❌ | `−2√5 μ c₀c₂` |
| **K₂, τ-locus** | **σ-even (in-field)** | ❌ | `2μ(−a₂a₃ − √5 i a₂b₃ + √5 i a₃b₂ − 5b₂b₃)` |

The last row is the decisive one: **even on the physical locus, with an honest in-field (σ-even) mass, there is no conserved symmetric charge.** So this is the **construction**, not merely Cliff's axiom.

**Root cause (three DERIVED facts).**
1. **m ∉ K.** m = 5^{1/4}μ has minimal polynomial x⁴−5 over ℚ (degree 4 > 2). The σ-odd-mass axiom takes the massive operator out of K. *(Fact B.)*
2. **σ fixes the phase.** σ(e^{iα}) = e^{iα}, so ψ̄ = (σψ)ᵀA → e^{+iα}ψ̄, and ⟨ψ,ψ⟩ → e^{2iα}⟨ψ,ψ⟩ — a **bilinear** form, not a phase-invariant charge (off the locus).
3. **The τ-constraint is Majorana.** e^{iα}·(f + i√5 g) has p-part e^{iα}f with Im = f sinα ≠ 0 → **leaves the τ-locus.** U(1) is broken *on* the locus too. σ-real fields have no conserved symmetric Dirac charge. *(Facts 2+3 = A.)*

**The conserved object.** F7 already found the *antisymmetric* symplectic form Ω (rank 4, det 125), and it **is** conserved / nondegenerate — the theory has a good symplectic structure. What it lacks is the **symmetric** Dirac inner product. This is the Majorana signature: a real symplectic form, no U(1) charge.

## §4 F8b — the K4 signature (provisional; downstream of the failed gate, NOT banked)

Computed exactly per the brief, and reported because it is diagnostic — but it is **not a K4 verdict**, because K4's conserved-charge premise (§3) fails.

- **(σu±)ᵀK₂Γ_s u± is COMPLEX** at both places, e.g. `⟨u₊,u₊⟩|∞₁ = √5(−2k² − 2k√(k²+μ²) − μ² + iμ²)`. Real part (kinetic) negative and **same sign** on u₊ and u₋; imaginary part (mass) from σ(m)·m = i√5μ² — a direct symptom of **m ∉ K** (Fact B). **Not** a clean opposite-definite pair → **K4 does NOT PASS**; but with no conserved symmetric charge the "signature" has no invariant meaning.
- **K4 map:** σ(u₊) is proportional to **neither** u₊ nor u₋ — the clean "+freq → −freq" flip is **not** exhibited (again because m ∉ K muddies the place reading).
- **§4.5 representative agreement:** K₂ and √5·K₂ give the **same signature kind** at every (mode, place) — consistent, so no §4.5 alarm.

## §5 F8c — ev₁-unitarity (K5): benign on the locus (DERIVED)

On the τ-locus there are **4 real dof** at ∞₁, and σ(ψ) is **determined by them** (σ = τ on the locus, F3b/F5 — verified). So ∞₂ is the *conjugate* of ∞₁ data, not independent — the **benign** case (like ordinary Dirac's ψ̄ ~ ψ†). **K5 does not fire.** (This is the one piece of good news, and it is exactly what would matter *if* a conserved charge existed.)

## §6 — genericity + the named tension

- **Genericity (STRUCTURAL):** the signature structure (Re < 0 kinetic, Im > 0 mass) is **identical for √2, √3, √5** — no d-dependence, real-quadratic-generic, not golden. **No cloud.**
- **Named tension, carried RAW (NOT adjudicated).** K₂ **fails** 226 §5.2's uniform-ε frontier while **passing** F7 reality. The pairing form and the conjugation matrix disagree about σ. Handed to Mr A / Tor. **Fizz's §7.1 operator/amplitude resolution was deliberately NOT tested** (her instruction; Mr Code off-limits). Reported raw.

---

## Escalation (§8) — what this means, and what's owed

Per §8: "escalate immediately, do not proceed, if conservation fails on K₂." It does. **F4/K4 are not banked.** The decisions belong to Cliff/Tor/Fizz/Mr A:

1. **The golden σ-construction is Majorana-like** (σ-real → no conserved symmetric U(1) charge; a conserved *symplectic* form instead). If 226 wants a **charged** Dirac field with a particle/antiparticle **charge**, this construction does not provide one at the free level. If 226 is content with a **Majorana / real** field (particle = antiparticle, symplectic structure), F7's real nondegenerate action **stands** and the "antiparticle" of §9 must be read as the Majorana (self-conjugate) kind, not a charged partner.
2. **The antiparticle axiom (m² σ-odd) is separately problematic:** it forces m ∉ K, taking the first-order operator off the field. Either the mass enters differently (not as m·s·Z with m ∉ K), or the equal-mass-antiparticle requirement is met by the Majorana structure (1) rather than by a σ-odd mass. **This is Cliff's axiom meeting the algebra — a decision, not a bench call.**
3. The 226 §5.2 frontier vs reality tension (K₂) remains open (Fizz's §7.1 resolution un-adjudicated).

**Positive residue:** F7 is **not** retracted — the real, nondegenerate action and its conserved symplectic Ω stand. F8's finding is about the *nature* of the surviving theory (Majorana/symplectic, not charged Dirac), not a defect in F7.

## Reproduce

```
cd paper_226_golden_construction
python bench_f8.py   # -> results_f8.csv; §3 conservation FAILS (all cells incl tau-locus/in-field mass);
                     #    roots m∉K + Majorana; §4 signature COMPLEX; §5 K5 benign; VERDICT ESCALATE
```
Environment: `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic; no RNG; floats only to read a sign.

🐕☕⬡
