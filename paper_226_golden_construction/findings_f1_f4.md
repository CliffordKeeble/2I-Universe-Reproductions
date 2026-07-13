# Findings — Benches F1–F4 (Free Golden Dirac Dynamics)

**Paper 226 · Golden Construction · free-dynamics benches**
**9 July 2026 · exact symbolic over K(i) = ℚ(√5, i) · `bench_f1_f4.py` · `results_f1_f4.csv` · `run_f1_f4.log`**
**Pre-registration:** `PRE-REGISTRATION_f1_f4.md` (committed `211833b`, before this verdict).
**Brief of record:** `mr_code_brief_free_golden_dynamics_F1-F4.md` (Tor, 8 Jul 2026).

Status tags per programme convention: **DERIVED** (exact symbolic proof), **STRUCTURAL** (mechanism-level, argued + computationally exhibited), **OBSERVED** (numeric). Nothing here is banked until merged and cold-checked by Mr Adversary (W-108).

> ## 🔭 SCOPE — 12 July 2026 (F9a connection check; see `findings_f9.md`)
>
> **The object of this entire arc carried NO icosahedral (2I) content.** F9a (DERIVED): the Paper 191 Clifford form ⟨√5,−√5⟩ is isotropic ≅ ⟨1,−1⟩ over K → Cl = split M₂(K), identical to standard (1,1); 2I (order 120, a subgroup of the *compact* SU(2)) does not embed in the abelian Spin(1,1) and does not respect the gammas (their centralizer is scalars). **The √5 in the gamma entries was placed there by hand and is a removable artefact** — every "real-quadratic-generic" result below (and it was *every* result) was this diagnosis. The theory is a real-quadratic (1,1) Dirac theory, **not golden**. Read all "golden" language as "over ℚ(√5)", not "icosahedral".

---

> ## ⚠️ CORRECTION — 10 July 2026 (Mr A three-star review; see `findings_f5.md`)
>
> Two claims below are demoted by Mr A's review of the no-go and must be read with these corrections:
>
> 1. **The "never diag(+,+)" mechanism is NOT a theorem.** By Taussky–Zassenhaus a symmetric (+,+) aligner exists for *any* anticommuting invertible pair; F5 step 1 confirmed it concretely for the 191 pair (the antidiagonal family **A = [[0,b],[c,0]]**). So the F1–F4 COLLISION is a **bench result about the presentation A = diag(√5,1), not a proof.** Whether the no-go survives now hinges on whether *any* (+,+) aligner is a **legal σ-sesquilinear pairing** (F5 step 2 — REQUISITIONED, not yet decided). Everything below tagged "STRUCTURAL / never diag(+,+)" is downgraded to **presentation-specific, pending F5**.
> 2. **"Real-quadratic-generic" is an overclaim on {5, 13}.** Both 5 and 13 are ≡ 1 (mod 4); genericity across d ≢ 1 (mod 4) (√2, √3) is **F6, not yet run**. Read the C8/genericity headline as **"generic on the fields tested {5, 13}"** until F6 lands.
>
> The exact-symbolic F1/F2/F3a/F3b/reconciliation/sub-check results and the C1/C6 numbers below are **unaffected** and were re-confirmed independently in `findings_f5.md`.

---

> ## ⛔ RETRACTION — 11 July 2026 (F7 reality sweep; see `findings_f7.md`)
>
> **The COLLISION headline below is RETRACTED. K1 does not fire for the theory.**
>
> F7 ran the reality test — *does the pairing give a real, nondegenerate first-order action?*, the actual K1 question — over the **whole aligner family** A = [[0,b],[c,0]], not just A₀. The reality-variety is **b = c**: the pairing **A = K₂ = [[0,1],[1,0]]** (and any b=c) gives a **real, nondegenerate (rank 4), EOM-preserving** action at healthy s = 1. **The F1–F4 COLLISION was an artefact of the specific pairing A₀ = diag(√5,1)** — exactly the presentation Mr A demoted. With an aligning pairing the free split σ-theory has a healthy first-order action, so **the no-go is withdrawn.**
>
> Everything computed below remains **arithmetically correct for A₀** — A₀ genuinely gives no real action. What is wrong is the *scope*: A₀ was never privileged, and the COLLISION does not hold class-wide. **Do not cite the COLLISION as a result.**
>
> **Two caveats (see `findings_f7.md`):** (1) K₂ **fails** 226 §5.2's frontier criterion — reality and the frontier **disagree** about σ (open question for Mr A/Tor, not adjudicated); (2) the physical sensibility of the antidiagonal/place-swap pairing (positivity, ± frequency split) is the re-posed F4 question, **deferred** — the retraction is "a real action exists", not "the QFT is fully vetted".

---

## Headline

> **F3 = COLLISION at the free level, at healthy s = 1.**
> There is **no field-reality condition** — τ, σ-real, or the C₀ Majorana, under **either** overall factor (1 or i) — that makes the free golden Dirac **kinetic** term real modulo total derivatives. The **mass** term is real (only) under the τ-constraint at s = 1, but the kinetic term is real under **nothing**. So no nondegenerate first-order *real* action exists: **K1 fires** (with K2, K3 in train). The free split σ-theory is sick at the free level; 226 §3.3's fork is resolved **from below** by a forced no-go. **The no-go is real-quadratic-generic** (identical for ℚ(√13)) — a rescope, **not** a golden discovery.

F4 was **gated off** per the §6 stop-condition (do not run F4 on a COLLISION). All 11 recorded items match their pre-registered EXPECT; the two "genuinely open" gate checks resolved as reported below.

---

## F1 — dispersion and anticommutant *(Tor's build — confirmed)*

- **F1a (DERIVED).** The joint anticommutant of {Γ_s, Γ_a} is exactly **span(Z)**, one-dimensional. Anticommuting with Γ_s alone gives X = [[a, b], [−√5 b, −a]]; imposing Γ_a forces b = 0. Confirms EXPECT.
- **F1b (DERIVED, convention-free).** M = −iEΓ_s + ikΓ_a + msZ has
  det M = **−m²s² − √5(k² − E²)**, giving **E² = k² + m²s²/√5**, i.e. **M² = m²s²/√5**.
  - **s = 1: M² = +1/√5 > 0 — healthy.**
  - **s = i: M² = −1/√5 < 0 — tachyonic.**
  This inverts the retracted §2.1 of Fizz's construction, exactly as Tor computed. The determinant is bare (no pairing convention enters), so the PASS also certifies the §0 matrices are transcribed correctly.

## F2 — pairing obstruction *(Fizz's build — confirmed)*

- **F2a (DERIVED).** With A = diag(√5, 1): **Γ_sᵀ = +A Γ_s A⁻¹** (parity **+1**) and **Γ_aᵀ = −A Γ_a A⁻¹** (parity **−1**). With A′ = diag(−√5, 1) the two signs exchange. The intertwined transpose action is diag(+,−) or diag(−,+) — **never diag(+,+)**. *This single fact is the engine of the F3 collision (see below).*
- **F2b (DERIVED).** N·Γ_s ∝ ε forces **N = λ·diag(1, −1/√5)**; then N·Γ_a = λ·[[0,−1],[−1,0]] is **symmetric**, not ∝ ε. No constant N antisymmetrises both ⟹ no neutral (Galois-trivial) bilinear supports a first-order action; **the σ-dressing is forced**.

## F3 — the σ/τ decider and the joint gate *(the crux)*

**F3a (DERIVED).** For the unconstrained σ-field, B(ψ) = ψ̄Zψ = √5·N(ψ₁) − N(ψ₂) with N(ψᵢ) = pᵢ² − 5qᵢ². Neither s = 1 nor s = i makes Im[ev₁(m s B)] vanish for generic configurations ⟹ a field-reality constraint is **mandatory**. Confirms EXPECT (not itself a kill).

**F3b (DERIVED).** The τ-constraint ev₂(ψ) = τ(ev₁(ψ)) forces **pᵢ ∈ ℝ, qᵢ ∈ iℝ**; then N(ψᵢ) = aᵢ² + 5dᵢ² (**real, positive**) and ev₁(B) is real. **With s = 1 the mass term is real and healthy.**

**Reconciliation (DERIVED).** On the τ-constrained locus, **σ(ψ) = τ(ψ) exactly** (each ψᵢ = aᵢ + i dᵢ√5 has σψᵢ = τψᵢ = aᵢ − i dᵢ√5). So the F3b constraint **is** precisely "σ = τ on the physical locus": the σ-Majorana and τ-Majorana conditions coincide on-shell (Paper 226 §5.3). F2's σ-sesquilinear *bilinear* and F3b's τ-type *field reality* are compatible on the locus. Had the kinetic gate cleared, this would have been the clean Fork B.

### The gate — sub-check 1 (kinetic reality): **FAILS under everything**

The decisive computation. Im[ev₁(ψ̄(Γ_s∂₀ + Γ_a∂₁)ψ)] modulo total derivatives, swept over three field-reality constraints × two overall factors:

| constraint | factor = 1 (leftover, mod TD) | factor = i (leftover, mod TD) |
|---|---|---|
| **τ** (F3b) | `5(f₁∂₀g₂ − ∂₀f₁·g₂ + f₂∂₀g₁ − ∂₀f₂·g₁)` ≠ 0 | `√5(−f₁∂₁f₂ + ∂₁f₁·f₂ − 5g₁∂₁g₂ + 5∂₁g₁·g₂)` ≠ 0 |
| **σ-real** | `√5(−f₁∂₁g₂ + …)` ≠ 0 | `√5(−f₁∂₁f₂ + …)` ≠ 0 |
| **C₀** (Zσ) | `5(f₁∂₀g₂ − … )` ≠ 0 | `5(f₁∂₀f₂ − …)` ≠ 0 |

**Every cell is non-zero.** The pattern is the whole story:

- **factor = 1** leaves a **∂₀ current** (the Γ_s, parity-**+1** direction); the Γ_a (∂₁) piece is a total derivative.
- **factor = i** leaves a **∂₁ current** (the Γ_a, parity-**−1** direction); the Γ_s (∂₀) piece is a total derivative.

Whatever overall phase reales one derivative direction makes the **other** imaginary. This is a direct consequence of **F2a**: AΓ_s is symmetric (parity +1) but AΓ_a is antisymmetric (parity −1), and no single reality convention can accommodate both parities at once. **STRUCTURAL.** (For calibration: ordinary Dirac/Schrödinger pass this test at factor = i, because there *both* kinetic directions carry the same transpose parity. The golden σ-adjoint's diag(+,−) — F2a's "never diag(+,+)" — is exactly what breaks it. The result is also robust to Hermitian symmetrisation: adding the h.c. is a total derivative iff the parities match, so the symmetrised action inherits the same obstruction.)

### Sub-check 2 (DERIVED) — EOM *does* preserve τ

Substituting ψᵢ = fᵢ + i√5 gᵢ into (Γ_s∂₀ + Γ_a∂₁ + mZ)ψ = 0 yields **four real, closed** equations for (fᵢ, gᵢ); no term forces the fields complex. The (a-sector) and (d-sector) evolve by the same golden-Dirac system. **The τ-reality road is dynamically available** — the constraint is not what breaks.

### Sub-check 3 (DERIVED) — the form is nondegenerate

The ∂₀ antisymmetric (symplectic) coefficient matrix on the τ-halved space has **rank 4, det = 625 ≠ 0** — **nondegenerate**. **The kinetic form is not degenerate.**

### What kind of no-go this is

Sub-checks 2 and 3 passing is the sharp part: the constraint evolves fine and the form is nondegenerate. **The failure is reality, and reality alone.** So K1 fires *specifically via its "no *real* action" clause* — a nondegenerate first-order symplectic structure exists, but no field-reality dressing makes it real. Precisely:

- **K1 (no dynamics):** no nondegenerate first-order **real** action — CONFIRMED (nondegenerate ✓, real ✗ under every dressing).
- **K2 (wrong home):** there is no healthy real action to home C = σ in ⟹ σ loses its dynamical home; the Fourth Face rescopes σ to **kinematics**. CONFIRMED in train.
- **K3 (no split):** the split σ-theory dies as physics before the ±-frequency structure can be posed; F4b (which would have tested the ± split) is gated off. Seagull-consistent; not independently tested here.

**Sharper than the brief's anticipated shape (FLAG).** The brief framed COLLISION as "the s making the action real is incompatible with healthy s = 1." What the bench actually finds is stronger and cleaner: the kinetic term carries **no** s (s lives only in msZ), so kinetic reality is s-independent and fails for **every** constraint and factor. It is not "mass wants s = 1, kinetic wants s = i"; it is "mass real under τ, **kinetic real under nothing**." The obstruction is intrinsic to the σ-adjoint's opposite transpose-parities (F2a), not to a clash of preferred s-values.

### Genericity (DERIVED) — real-quadratic, not golden

Re-running F2a's transpose parity and the sub-check-1 kinetic obstruction with √5 → √13 gives the **identical** diag(+,−) parity and the **same** non-zero ∂₀ current (`13(f₁∂₀g₂ − …)`). **The no-go holds for ℚ(√13) identically.** Per §1 this is a **real-quadratic no-go, not a golden one** — a rescope, and it must not be reported as a golden discovery. √5 is load-bearing nowhere in F1–F4, exactly as registered.

## F4 — gated off

Per §6: on a COLLISION, **do not run F4**. The pairing conservation (F4a), the K4 frequency signature (F4b), and ev₁(S)-unitarity (F4c) are not evaluated — the free split σ-theory they would characterise is sick at the free level. `bench_f1_f4.py` contains the F4a/F4b scaffolding, gated behind the computed CLEAR flag, ready if CinC re-poses the reality question (e.g. on the compact / s = i route, or with a different adjoint).

---

## Flags for CinC / Tor

1. **The verdict is COLLISION, against the reconciliation section's hopeful "Fork B lives" framing.** Sub-check 1 was pre-registered as genuinely open ("no pre-registered direction"), so this is within the registered space — but it lands on the *opposite* side of Tor's un-banked lean. Reporting it straight, per Pattern 39.
2. **The collision is s-independent and reality-only** (sub-checks 2, 3 pass). The mechanism is F2a's diag(+,−) transpose parity — AΓ_s symmetric vs AΓ_a antisymmetric. This is a **cleaner and more general kill** than an s-clash: it says the σ-sesquilinear *split* adjoint cannot support a real first-order action for *any* real-quadratic field. Fizz's F2 obstruction and Tor's F1∧F3 collision are the *same* fact seen from two sides.
3. **Real-quadratic-generic, confirmed on √13.** This is a **rescope**, not a result about golden. It supports 226 rescoping the Fourth Face / C = σ to **kinematics** (the split σ-theory has no free dynamical home), and points the surviving dynamical content to the **compact / s = i route** ("compact-or-bust").
4. **Convention caveat, stated for the cold check.** I applied the brief's literal reality criterion (Im[ev₁ L] mod TD = 0) and *additionally* swept the overall factor ∈ {1, i} so the collision cannot be dismissed as a missing-i artifact — ordinary Dirac passes the factor-i test; the golden split adjoint fails both. If Mr Adversary prefers the Hermitian-symmetrised action as the primary object, the verdict is unchanged (the symmetrisation is a total derivative iff F2a's parities match, which they don't).

## Reproduce

```
cd paper_226_golden_construction
python golden_algebra.py     # substrate self-test — all PASS
python bench_f1_f4.py        # -> results_f1_f4.csv, verdict COLLISION, 11/11 vs EXPECT
```
Environment: see `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic throughout; no RNG; no floating point (F4's float-to-report-sign step was gated off).

🐕☕⬡
