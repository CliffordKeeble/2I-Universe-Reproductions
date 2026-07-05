# Findings — Golden Construction Bench Run #1 (Part A)

**Target ("Paper 226"):** the Golden Interacting Theory construction (`golden_interacting_theory_construction_v0_1.md`, Fizz, 5 Jul 2026) and its analysis predecessor (`fizz_analysis_cpt_golden_clifford_v0_1.md`). Part A is that document's §7 Mr-Code bench list, and it blocks the construction's v0.2.
**Method:** exact arithmetic over ℚ(√5, i) via sympy (`sqrt(5)`, `I`); high-precision mpmath used where symbolic simplification proved unreliable (A8, and the A2 conjugation hunt). Pre-registration: [PRE-REGISTRATION.md](PRE-REGISTRATION.md), committed before any verdict.
**Result:** **8 / 10 items match EXPECT.** The two that do not (A5, A8) are *results*, not errors — both were pre-registered as the most informative outcomes, and both were verified free of computational artifact (see notes).

## Summary table

| Item | EXPECT | Outcome |
|---|---|---|
| A1 block relations | all pass | **PASS** |
| A2 L1 invariance hunt | preserve (a)–(d); flip (e),(f) | **PASS** |
| A3 C₀ explicit | all pass | **PASS** |
| A4 weights / trace-zero | all pass | **PASS** |
| A5 toy-model end-to-end (F3) | pass conditionally | **FAIL of EXPECT — F3 fired (see below)** |
| A6 volume-element squares | +25(×6), −25, −5^{3/2} | **PASS** |
| A7 eigenvalue fields | all pass | **PASS** |
| A8 descent invariant J | J unique; J²=−1; +1 impossible | **FAIL of EXPECT — J²=−1/√5 indefinite (see below)** |
| A9 σ flips signature | pass | **PASS** |
| A10 ring bookkeeping | all pass | **PASS** |

Full machine output: [results.csv](results.csv), [run.log](run.log).

---

## The eight passes (DERIVED, exact over ℚ(√5, i))

- **A1** — [Z,Γ_seed]=−2Γ_adj, [Z,Γ_adj]=−2Γ_seed, {Z,Γ_seed}={Z,Γ_adj}=0, Γ_seed·Γ_adj=+√5·σ₃. All five hold exactly. (Re-confirms Paper 203's machine result and its companion.)
- **A2** — For each of the 16 √5-normalised candidates, the square-sign is preserved under conjugation (96 exact structured + 320 high-precision-numeric random invertibles over ℚ(i,√5)), transpose, negation, and τ (48 checks); and flipped under σ and στ (32 checks). Zero anomalies. This is the invariance-lemma **L1** the construction's §3 Step 2 rests on — no counterexample found. (The proof route via Skolem–Noether is separate, per W-107; this is the counterexample hunt, and it is clean.)
- **A3** — C₀ = Z∘σ: C₀(Γ_seed)=Γ_adj, C₀(Γ_adj)=Γ_seed, C₀²=id on both the algebra and on K-rational vectors. The construction's §3 Step 3 block form is exactly reproduced.
- **A4** — The norm-one torus scaling (u, u⁻¹) preserves both P(λ)=⟨λ|ψ⟩_seed·⟨λ|ψ⟩_adj and ⟨ψ_seed|ψ_adj⟩; equal-weight (u,u) does not; and a=p+q√5 with a+σ(a)=0 ⟺ p=0, i.e. **Lie(T_σ)=ℚ·√5**. The construction's §1 centrepiece stands.
- **A6** — (Γ⁵₍₄₎)² = +25·I₄ on all **six** (2,2) cliques (re-confirms Paper 204); one explicit i-dressed (1,3) clique over K(i) has ω² = −25·I₄; three-gamma (Γ¹Γ²Γ³)² = −5^{3/2}·I₄ (Paper 191/203). All three targets hit.
- **A7** — minpoly of 5^{1/4}·i over ℚ is x⁴−5 (Eisenstein at 5); ℚ(5^{1/4}·i) has 2 real + 1 complex-conjugate-pair embeddings (**not CM**); i∉ℚ(5^{1/4}·i) by degree. Grounds the analysis's PR-4 resolution (§5).
- **A9** — σ (entrywise) sends the i-dressed (1,3) clique to a (3,1) clique: square-sign pattern (−,−,−,+) → (+,+,+,−), anticommutation preserved. Adjacent to Q204-2; recorded as data.
- **A10** — ℤ[φ]: N(φ)=−1, N(φ²)=+1, and N(±φ^k)=(−1)^k so ±φ^{2n} are exactly the (small) norm-one units. ℤ[√5]: φ³=2+√5 (N=−1, fundamental), φ⁶=9+4√5, 9²−5·4²=1, norm-one subgroup ±φ^{6n}. The two rings kept distinct.

---

## A5 — the F3 leak (the informative one)

**EXPECT:** pass conditionally; *"a failure or an ambiguity here is the most informative outcome on the whole sheet — it is exactly falsifier F3; report the precise point of breakage verbatim; do not patch."*

**Result: F3 fired.** The literal charge-conjugation operator **C₀ = Z∘σ does not map solutions of the q-equation L_q to solutions of the (−q)-equation L_{−q}** in the 1+1 golden Dirac toy model
> (Γ_seed D₀ + Γ_adj D₁)ψ = 0, D_μ = ∂_μ − q·a_μ·Z, W = Z = diag(+1,−1).

Testing C₀ as an exact all-k intertwiner (`σ(M_q) = Z·target·Z`) against the four candidate (−q) targets, the **only** target it intertwines is the operator with **both** (i) Γ_seed↔Γ_adj swapped relative to the two derivatives — equivalently the coordinate relabel x⁰↔x¹ — **and** (ii) the gauge field flipped a_μ→−a_μ.

- The **free/kinetic part** (a=0) maps correctly *only under the sector+coordinate swap*: because Γ_seed²=+√5 (timelike-like) and Γ_adj²=−√5 (spacelike-like) are **inequivalent** and sit on different derivatives, a sector-swapping C is a symmetry of the kinetic operator **only if time and space are also swapped**.
- The gauge-field flip a_μ→−a_μ is benign (standard QED C-oddness of the photon). The residual obstruction against the same-a swapped target is **purely the minimal-coupling term**:
  > residual = ( 0, 2q(a₀−a₁), 2√5·q(a₀+a₁), 0 ) — ∝ q·a_μ, no k-dependence.

**Reading (STRUCTURAL, for CinC/Fizz — not a Mr Code claim to bank):** A1–A4 all hold, so the construction's conditional antecedent is satisfied; the failure therefore localises to the construction's named assumption **A2** (state grading ↔ operator grading alignment). In this toy model the alignment requires an *extra time↔space swap* that the construction's §3 does not license. This is precisely "a sign that decouples the state grading from the operator grading." **Not patched.** Whether the intended 1+1 equation, or the intended action of C on coordinates/gauge field, differs from the literal brief wording is a construction-level decision (a re-brief), not a bench fix.

## A8 — the descent invariant is indefinite, not −1

**EXPECT:** J exists uniquely up to scale; J²=−1; renormalisation to +1 impossible. *"A J²=+1 outcome would contradict Paper 204's Theorem 1 via descent — if found, stop everything and report."*

**Result (numeric-authoritative — see methodological note):** For the i-dressed (1,3) clique, the antilinear J = B∘τ commuting with all four gammas exists and is **unique up to scale (dim B = 1)** ✓. But
> **J² = (−1/√5)·I** — an **indefinite** element of K: J² < 0 at the +√5 real place, J² > 0 at the −√5 place.

Renormalisation J→λJ (λ∈K(i)^×) scales J² by the norm N(λ)=λ·τ(λ), which is **totally positive** (a sum of two squares in K). The (−,+) sign pattern of −1/√5 is therefore invariant, and J² can be normalised to **neither −1** (would need N(λ)=√5, indefinite → not a norm) **nor +1** (would need N(λ)=−√5). So:

- EXPECT's "renormalisation to +1 impossible" — **holds** ✓ (and it is the same obstruction, refined).
- EXPECT's headline "J²=−1" — **does not hold**; the invariant is indefinite (−1/√5).
- **No J²=+1 was found**, so the "stop everything, contradicts Paper 204" alarm does **not** trigger. Paper 204 Thm 1 concerns *standard* Cl(1,3) (squares ±1); this clique is **√5-normalised** (squares ±√5), which Paper 204 §8 and the analysis §5 both explicitly flag as *not* standard Cl(1,3) ("rescaling to standard relations would leave K" — 5^{1/4}∉K(i)). The √5 in the squares is what makes the descent invariant indefinite. The finding **qualifies** the analysis §4.2 [STRUCTURAL] expectation "J²=−1": that value is the standard-normalisation value; the golden √5-normalised clique gives −1/√5.

**Methodological note (honest flag):** sympy's `simplify` is **unreliable** on these nested M₄(ℚ(√5,i)) products — it reported J²=+1/5 (totally positive) where the truth is −1/√5 (indefinite). The discrepancy was caught by an internal consistency check (`B·τ(B)` vs applying J twice disagreed) and resolved with a high-precision mpmath computation (scalar residual 0, intertwining residual ~1e-40, identical across all six cliques and both dressing signs, and validated against the known all-real (2,2) case which gives a clean J²=+1). The A8 verdict rests on the numeric computation, not the symbolic one. Diagnostics: [a8_probe.py](a8_probe.py), [a8_numeric.py](a8_numeric.py).

---

## Resolved brief ambiguity (reported, not silently fixed)

The brief sources the 16 √5-normalised candidates as `{I,K,J,Z}⊗{Γ_seed,Γ_adj}` — that enumerates only **8**. Paper 204 §5's authoritative 16 also includes the mirror ordering `{Γ_seed,Γ_adj}⊗{I,K,J,Z}`. Per the brief's own instruction ("the paper is the authority") I used all 16; `golden_algebra.py` independently reproduces Paper 204's count (exactly 16 with square ±√5·I₄, exactly 6 four-cliques, all (2,2)).

## Recommendations for the construction v0.2 (routing to CinC / Fizz)

1. **A5 / F3 is the headline.** The σ-semilinear C₀ derived in §3 does not, as literally written, act as charge conjugation on the 1+1 toy dynamics without an accompanying time↔space swap. Before v0.2 banks §3–§4's CPT verdict, the construction needs to state explicitly how C acts on coordinates and on a_μ, and confirm the resulting operator is an acceptable symmetry (a time↔space swap is *not*). This is the assumption-A2 gap the bench was built to find.
2. **A8 qualifies the descent claim.** The analysis's "J²=−1" is the standard-Cl(1,3) value; the golden √5-normalised (1,3) clique has J²=−1/√5 (indefinite). Any Θ²/CPT arithmetic that wants to source a −1 from the descent (analysis §3.5, wound #3) should not assume −1 from this clique. The clean −1 lives in standard Cl(1,3), which is **not representable over K(i)** at √5-normalisation.
3. Everything the construction leans on that is *not* the C-derivation — the gauge-torus centrepiece (A4, Lie(T_σ)=ℚ√5), the L1 invariance lemma (A2), the C₀ block algebra (A3), the volume-element and eigenvalue-field structure (A6/A7), the ring bookkeeping (A10) — **passed exactly**.

## Reproducibility

- `python golden_algebra.py` — substrate self-test (definitions, Galois maps, 16 candidates, 6 cliques).
- `python bench_run_1.py` — A1–A10, writes `results.csv`, `run.log`. Random seed 20260705.
- `python a8_numeric.py` — high-precision descent cross-check (authoritative for A8).
- Environment: sympy 1.14.0, mpmath 1.3.0, numpy 2.3.4, Python 3.13.
- Paper .md drafts are **not** committed (unreleased, per Cliff); definitions are vendored into `golden_algebra.py` with source citations.
