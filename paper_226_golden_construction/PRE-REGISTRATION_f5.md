# PRE-REGISTRATION — F5 aligner σ-type crux + Mr A's requisition battery

**Paper 226 · Golden Construction · F5 + lemmas + rechecks**
**Status: pre-registration · committed before `bench_f5.py` emits any verdict · 10 July 2026 · internal**

**Brief of record:** `mr_code_brief_aligner_sigma_type_F5.md` (Fizz, 10 Jul 2026), itself born of
Mr A's three-star review of `adversary_brief_free_golden_dirac_nogo.md`. Standard: Pattern 75
(pre-register before verdict); W-108 (cold Adversary boots without Mr A's report in context).

**What changed (brief §0).** Mr A demoted my linchpin: **C2's unrestricted "never (+,+)" is FALSE.**
By Taussky–Zassenhaus (Pacific J. Math 9, 1959) every square matrix is similar to its transpose
through a nonsingular **symmetric** intertwiner, so a (+,+) aligner **always exists** for any
anticommuting invertible pair. The F1–F4 no-go is therefore a **bench result about one presentation
(A = diag(√5,1)), not a theorem.** The live question is now F5.

**Substrate.** Exact symbolic over K(i) = ℚ(√5, i) via sympy; Galois maps + gamma basis reused from
the self-tested `golden_algebra.py`. Objects per brief §2/§0.

---

## Registered items

### F5 — the aligner σ-type (THE CRUX) · registered GENUINELY OPEN, no direction

- **F5 step 1 (run this session).** Solve (AΓ_s)ᵀ = AΓ_s **and** (AΓ_a)ᵀ = AΓ_a for constant
  A ∈ M₂(K(i)); report the full solution family. **EXPECT to VERIFY (Fizz pencil-work, not trusted):**
  the family is **antidiagonal** A = [[0, b], [c, 0]] (b, c ∈ K(i)), giving AΓ_s = diag(b√5, c),
  AΓ_a = diag(b√5, −c) — both diagonal, hence the (+,+) alignment. Contrast the original σ-pairing
  A = diag(√5, 1) → (+,−). *If this family is wrong, all of §2 restarts.*
- **F5 step 2 (σ-legality verdict) — BLOCKED, REQUISITION.** For each family member, apply 191/226's
  σ-sesquilinear **legality criterion** — "the axioms that make ψ̄ = (σψ)ᵀA a legal Dirac pairing
  for C = σ; **take these from the displayed dressings, do not infer them**" (brief §2.2). The 191/226
  drafts are **uncommitted** (Downloads only; see `cold_boot/PROVENANCE.md`), so the authoritative
  dressings are **not on disk**. I will **not infer** the criterion. **Requisition raised**; step 2 is
  not run this session. Kills gated on it:
  - **K-COLLAPSE:** any aligner is a legal σ-sesquilinear pairing → no-go **retracted**, K1 un-fires,
    re-pose F3/F4 on the legal aligner (Pattern 39 — report straight if it fires).
  - **Survive-and-upgrade:** every aligner is σ-odd → no-go survives, upgrades linear-algebraic →
    arithmetic; restate C2 as a class-statement; C7 becomes provable (co-descent).
- **Place-swap structural tell — CONJECTURED, not banked.** The aligner is off-diagonal (pairs
  component 1 with component 2). I report **neutral structural data only** (family, nondegeneracy,
  σ-action, the cross-component bilinear (σψ)ᵀAχ, the special member that equals Γ_s) as **inputs**
  to the legality test. I do **not** adjudicate legality from the shape.

### §7 rechecks — run regardless · registered EXPECT

- **C1.** Re-derive det(−iEΓ_s + ikΓ_a + msZ) by an **independent route** (hand-formula + explicit
  kernel vector u(k), not sympy `.det()`). **EXPECT:** E² = k² + m²s²/√5; the sign is carried by
  Z² = +I (msZ squares to +m²s²).
- **C6.** Re-verify the τ-halved symplectic coefficient matrix by an **independent construction**
  (direct unit-basis plug, not the diff-based extraction of `bench_f1_f4.py`). **EXPECT:** rank 4,
  det 625 (= 5⁴). *Note registered now:* the determinant **value** is basis-dependent (it reflects
  the √5-normalisation of the q-components); **rank 4 (nondegeneracy) is the basis-independent claim.**

### Corrections owed regardless of F5 outcome (brief §6, §8)

- **C2 restatement.** "Never (+,+)" is false as an unrestricted claim. Restate as a statement about
  the **admissible pairing class** (σ-legal pairings), pending F5. Annotate `findings_f1_f4.md`.
- **C8 headline.** May **not** say "real-quadratic-generic" on {5, 13} alone (both ≡ 1 mod 4). Soften
  to "generic on the fields tested"; full field-genericity (√2, √3; d ≢ 1 mod 4) is F6, **gated on
  SURVIVE** (brief §8), not run this session.

### Downstream — GATED, not run this session

- **§3 lemmas (B / D / C7-repair)**, **§4 σ/m² ledger** (needs the σ-action-on-m²-and-s convention —
  **226/Tor's to pin, requisition, do not invent**), **§5 second-order σ-invariance**, **§6 F6
  genericity**. All gated on **F5 = SURVIVE** and/or the two requisitions. Per §8 ordering.

## Fixed choices (decided now, reported)

1. Aligner family solved with 4 free K(i) entries; solution space dimension + basis reported exactly.
2. C1 independence = 2×2 hand determinant **and** verification that u(k) = (i(E+k), ms)ᵀ is an exact
   null vector on-shell (Mu = 0), cross-checked numerically at both real places √5 → ±√5.
3. C6 independence = build the 4×4 ∂₀ antisymmetric (symplectic) matrix by evaluating the reality
   defect on unit basis fields, then rank + det — a different mechanism from `bench_f1_f4.py`'s
   double-`diff`. Same (fᵢ, gᵢ) basis so the value 625 is comparable.
4. No RNG. No floating point except C1's two-place numeric cross-check (a sign/consistency read).

**Deliverables:** `bench_f5.py`, `results_f5.csv`, `findings_f5.md`, correction banner on
`findings_f1_f4.md`. Nothing banks until the battery clears **and** the cold Adversary pass certifies.

🐕☕⬡
