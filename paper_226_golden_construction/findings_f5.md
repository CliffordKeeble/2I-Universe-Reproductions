# Findings — F5 aligner σ-type crux (step 1 + rechecks); the crux itself is REQUISITIONED

**Paper 226 · Golden Construction · F5 battery**
**10 July 2026 · exact symbolic over K(i) = ℚ(√5, i) · `bench_f5.py` · `results_f5.csv` · `run_f5.log`**
**Pre-registration:** `PRE-REGISTRATION_f5.md` (committed `683f9c2`, before this verdict).
**Brief of record:** `mr_code_brief_aligner_sigma_type_F5.md` (Fizz, 10 Jul 2026), from Mr A's three-star review.

Status tags: **DERIVED** (exact proof), **STRUCTURAL/CONJECTURED** (shape-level, for interpretation, not banked).

---

## Headline

> **Fizz's aligner family is correct — verified independently (DERIVED).** The (+,+) aligner exists and is exactly the antidiagonal family **A = [[0, b], [c, 0]]**, b, c ∈ K(i). So Mr A is right and my F1–F4 linchpin ("never (+,+)") is correctly demoted: **the parity divorce was a fact about the presentation A = diag(√5, 1), not a theorem.**
>
> **The crux — F5 step 2, the σ-legality verdict — is NOT computed. It is BLOCKED on a requisition the brief itself names:** the σ-sesquilinear legality axioms must be taken "from the displayed dressings, do not infer them" (§2.2), and the 191/226 drafts are uncommitted (Downloads only; `cold_boot/PROVENANCE.md`). **I will not guess the crux warm.** Requisition below.
>
> The §7 rechecks (run regardless) both **confirm**: C1 dispersion, C6 symplectic rank 4 / det 625.

---

## F5 step 1 — the aligner family (DERIVED)

Solving **(AΓ_s)ᵀ = AΓ_s and (AΓ_a)ᵀ = AΓ_a** for constant A ∈ M₂(K(i)), by an independent linear solve:

- The solution family is **2-dimensional, antidiagonal**: **A = [[0, b], [c, 0]]** (a = d = 0 forced; b, c free).
- **AΓ_s = diag(b√5, c)** and **AΓ_a = diag(b√5, −c)** — both diagonal, hence symmetric: the **(+,+) alignment**.
- Contrast the original σ-pairing **A₀ = diag(√5, 1)**, which gives AΓ_s symmetric / AΓ_a antisymmetric — F2a's **(+,−)**.

**Fizz's pencil-work reproduces exactly.** Taussky–Zassenhaus is confirmed in the concrete: a symmetric (+,+) aligner exists for the 191 pair. Whether it is a *legal* C = σ pairing is F5 step 2 (below).

## Structural data on the aligner (CONJECTURED — neutral inputs for the legality test, NOT a verdict)

Reported for Tor/CinC to interpret; the legality criterion decides, **not** this shape:

- **Nondegenerate**: det A = −bc ≠ 0 whenever b, c ≠ 0. Every non-trivial aligner is invertible.
- **Cross-component (place-swap-shaped)**: the bilinear (σψ)ᵀAχ = **b·(σψ)₁·χ₂ + c·(σψ)₂·χ₁** pairs component 1 with component 2 — versus A₀'s **same-component** norms √5·(σψ)₁χ₁ + (σψ)₂χ₂. If the two spinor components track the two real places, this off-diagonal pairing has the shape of a **place-swap**.
- **A = K₂·diag(c, b)**: literally component-swap × diagonal.
- **The member A(b=1, c=√5) is Γ_s itself.** The seed gamma is a (+,+) aligner. (Also A(b=1, c=−√5) = the "reversed" seed.)

This is exactly the fork Fizz flagged: an antidiagonal/place-swap aligner **either** collapses straight into the C = σ structure (no-go dies) **or** is a σ-odd / degenerate-for-C=σ object (no-go upgrades to arithmetic). **The shape alone does not decide it** — the requisitioned axioms do.

## §7 rechecks (run regardless — both DERIVED, EXPECT met)

- **C1 (dispersion).** Independent route — **hand 2×2 determinant** (not sympy `.det()`) **plus** an explicit **kernel-vector** check that u(k) = (i(E+k), ms)ᵀ satisfies **Mu = 0 on-shell**, cross-checked numerically at **both real places** (√5 → ±√5). Confirms **E² = k² + m²s²/√5** for s ∈ {1, i}; the sign is carried by **Z² = +I** (msZ squares to +m²s²). ✓
- **C6 (symplectic).** Independent construction — the τ-halved symplectic matrix built as a **direct pairing-matrix** Ω_{jk} = Im[(σvⱼ)ᵀ(AΓ_s)v_k] on the constrained basis {f₁, g₁, f₂, g₂}, **not** the double-`diff` extraction of `bench_f1_f4.py`. Result: **rank 4, det = 625 = 5⁴.** ✓ (Note: **rank 4 / nondegeneracy is the basis-independent claim**; the value 625 reflects the √5-normalisation of the q-components — a basis artefact, not a physical invariant. Checked, as Mr A asked.)

---

## ⛔ Requisitions (blocking — the crux and §4 cannot proceed without these)

Per the brief's own instruction ("take these from the displayed dressings, do not infer them"; "226/Tor's to pin, requisition it, don't invent it"), and CLAUDE.md (no reconstructing physics from papers not in hand):

1. **F5 step 2 — the σ-sesquilinear legality criterion.** I need the **displayed dressings from Papers 191/226** that define when ψ̄ = (σψ)ᵀA is a *legal* Dirac pairing for C = σ — i.e. the pairing axioms (the "three dressings as displayed equations" §8 also asks Cliff/CinC to hand Mr A). With those on the page I test each aligner A = [[0, b], [c, 0]] and return **σ-legal (→ K-COLLAPSE, no-go retracted)** or **σ-odd (→ survive-and-upgrade to arithmetic)**. Without them I do not compute the verdict.
2. **§4 — the σ-parity convention for m² and s.** Which convention do 191/226 adopt for σ's action on the mass-squared parameter m² and on s? This one ledger entry decides C, G, and much of C9. **226/Tor to pin on record.** I have both branches ready to compute the instant it lands:
   - m² σ-**even** (m² ∈ ℚ) ⟹ mass shell m²s²/√5 is σ-**odd** ⟹ σ swaps healthy (s=1) ↔ tachyonic (s=i) ⟹ C = σ was **never** a symmetry of the massive free dynamics (the 30-second dispersion reading; the first-order no-go merely confirms it by heavy labour).
   - m² σ-**odd** (m² ∈ √5·ℚ) ⟹ mass shell σ-**even** ⟹ reopens the symmetry question and the second-order question (§5).

## Corrections owed regardless of F5 (applied)

- **C2** restated. "Never (+,+)" is **false** as an unrestricted claim (Taussky–Zassenhaus; verified here). The live claim narrows to the **admissible (σ-legal) pairing class**, pending F5 step 2. A dated correction banner is added to the top of **`findings_f1_f4.md`**.
- **C8** headline softened. It may **not** read "real-quadratic-generic" on {5, 13} alone — **both are ≡ 1 (mod 4)**. Corrected to **"generic on the fields tested {5, 13}"**; d ≢ 1 (mod 4) confirmation (√2, √3) is **F6, gated on SURVIVE**, not run this session. Banner records this too.

## Gated / not run this session (per §8 ordering)

- **§3 lemmas** (B: all-phases-from-two-points; D: mixed total derivatives; C7 repair/withdrawal) — gated on **F5 = SURVIVE**. *Note, ready:* §3a's premise is already in hand — the F3 defect X₀ (factor 1) is supported **purely on ∂₀**, X₁ (factor i) **purely on ∂₁** (see `findings_f1_f4.md` / `run_f1_f4.log`), so the two-bench-points-plus-lemma proof-over-the-circle is writable the moment the no-go survives.
- **§4 → §5** (σ/m² ledger → second-order σ-invariance) — gated on **requisition 2**.
- **§6 F6** (√2, √3 genericity) — gated on **SURVIVE** (per §8; also substantiates the C8 correction).

**Nothing banks** until the battery clears **and** the cold W-108 Adversary pass certifies (booted without Mr A's report in context). Separately, per §8, Cliff/CinC owe Mr A's desk: `findings_f1_f4.md`, Paper 191's gamma basis, and the three dressings as displayed equations — the same artefacts I've requisitioned for step 2.

## Reproduce

```
cd paper_226_golden_construction
python bench_f5.py     # -> results_f5.csv; F5.1 PASS, C1 PASS, C6 PASS; step 2 not run (requisition)
```
Environment: `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic; no RNG; floats only in C1's two-place cross-check.

🐕☕⬡
