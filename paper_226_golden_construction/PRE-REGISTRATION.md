# Pre-registration — Golden Construction Bench Run #1 (Part A)

**Target:** the Golden Interacting Theory construction (`golden_interacting_theory_construction_v0_1.md`, Fizz, 5 Jul 2026) — the "Paper 226" the brief scores against — and its analysis predecessor (`fizz_analysis_cpt_golden_clifford_v0_1.md`).
**Brief:** `mr_code_brief_bench_run_1_and_hecke_locks.md` (Fizz, 5 Jul 2026), Part A.
**Discipline:** each item carries an EXPECT line authored by Fizz; I report pass/fail against EXPECT verbatim. A failed EXPECT is a result, not an error. Committed **before** any verdict output (`bench_run_1.py`) is produced.

**Pass criterion:** exact symbolic equality over ℚ(√5, i) via sympy (`sqrt(5)`, `I`). No floating point, no tolerance, except A7/A8 embedding-sign checks which evaluate a fixed exact element at the two real places of K and read its sign.

**Sources of definitions (confirmed against text before running; paper is authority):**
- Γ_seed, Γ_adj, K, J, Z: Paper 203 v0.3 §4 (quoting Paper 191 v1.2). Confirmed: Z = ½[K,J], Γ_seed·Γ_adj = √5·Z, [Z,Γ_seed] = −2Γ_adj, σ(Γ_seed) = −Γ_adj (all re-checked in `golden_algebra.py` self-test — PASS).
- The 16 √5-normalised candidates and the six (2,2) cliques: Paper 204 v0.4 §5 / App. A.
- Galois: τ (i↦−i, √5 fixed) = entrywise complex conjugation; σ (√5↦−√5, i fixed) = subs √5→−√5; στ = both.

**Resolved brief ambiguity (reported, not silently fixed):** the brief sources the 16 candidates as `{I,K,J,Z}⊗{Γ_seed,Γ_adj}`, which is only **8**. Paper 204 §5's authoritative 16 also includes the mirror ordering `{Γ_seed,Γ_adj}⊗{I,K,J,Z}`. I use all 16 per Paper 204 ("the paper is the authority"). `golden_algebra.py` confirms exactly 16 with square ±√5·I₄ and exactly 6 four-cliques, all (2,2).

---

## EXPECT sheet (verbatim from the brief) + fixed construction choices

- **A1** Block relations: [Z,Γ_seed]=−2Γ_adj; [Z,Γ_adj]=−2Γ_seed; {Z,Γ_seed}=0; {Z,Γ_adj}=0; Γ_seed·Γ_adj=+√5·σ₃. **EXPECT: all pass.**
- **A2** Invariance-lemma hunt (L1): for each of 16 candidates Γ (Γ²=±√5·I₄), sign of Γ² before/after (a) MΓM⁻¹, (b) transpose, (c) negation, (d) τ, (e) σ, (f) στ. **EXPECT: sign preserved under (a)–(d) every case; flipped under (e),(f) every case.**
  - *Choice:* generator sample = the three Paper-191 gammas {Γ¹=K⊗Γ_seed, Γ²=J⊗Γ_seed, Γ³=I⊗Γ_adj} plus building-block tensors {Z⊗I, I⊗Z, K⊗K}. Random invertibles: **seed 20260705**, 20 matrices with entries a+b·i+c·√5+d·i√5, coefficients drawn uniformly from {−2,−1,0,1,2}, singular matrices rejected and redrawn.
- **A3** C₀: X↦Z·σ(X)·Z⁻¹ on algebra, v↦Z·σ(v) on vectors. Verify C₀(Γ_seed)=Γ_adj, C₀(Γ_adj)=Γ_seed, C₀²=id (algebra and vectors). **EXPECT: all pass.**
- **A4** Weights/trace-zero: (u,u⁻¹) preserves P(λ)=⟨λ|ψ⟩_seed·⟨λ|ψ⟩_adj and ⟨ψ_seed|ψ_adj⟩; (u,u) does not; a=p+q√5 with a+σ(a)=0 ⟺ p=0. **EXPECT: all pass.**
  - *Choice:* amplitudes ⟨λ|ψ⟩_seed, ⟨λ|ψ⟩_adj and ⟨ψ_seed|ψ_adj⟩ carried as symbols; sector scaling acts as scalar multiplication. Symbolic u≠0.
- **A5** Toy model end-to-end (F3): (Γ_seed D₀ + Γ_adj D₁)ψ=0, D_μ=∂_μ − q a_μ W, a_μ constant, **W = Z = diag(+1,−1)** (weight +1 seed / −1 adjoint, per construction §2 and 192 §5 state pair). Plane-wave ansatz ψ=exp(−i(k₀x⁰+k₁x¹))·χ. Construct explicit q-solution, apply C₀=Z∘σ (acting (C₀ψ)(x)=Z·σ(ψ(x)); σ on K(i)-coefficients only, coordinates and i untouched), verify image solves the **(−q) equation** L_{−q}ψ=0, L_{−q}=Γ_seed(∂₀+q a₀ Z)+Γ_adj(∂₁+q a₁ Z). **EXPECT: pass — conditionally; a failure or ambiguity is the most informative outcome (falsifier F3); report the precise point of breakage verbatim; do not patch.**
  - *Choice:* the literal reading above is fixed now. On-shell residual = L_{−q}(C₀ψ_q) evaluated on the dispersion locus det M_q(k)=0, simplified exactly. Free case (a_μ=0) computed separately to localise any leak to the kinetic vs coupling sector. Any interpretive ambiguity in how coordinates / sectors / a_μ transform under C is itself recorded as a finding about grading alignment (per brief).
- **A6** Volume-element squares: (Γ⁵₍₄₎)²=+25·I₄ for all six (2,2) cliques; one explicit i-dressed (1,3) clique over K(i) with ω²=−25·I₄; three-gamma (Γ¹Γ²Γ³)²=−5^{3/2}·I₄. **EXPECT: +25 (×6), −25, −5^{3/2}.**
  - *Choice:* i-dressed clique = clique index 0 (deterministic ordering from `find_cliques`); dress the **first** of its two +√5 gammas by ×i.
- **A7** Eigenvalue fields: minpoly of 5^{1/4}·i over ℚ is x⁴−5; ℚ(5^{1/4}·i) has 2 real + 1 complex-conjugate-pair embeddings (not CM); i∉ℚ(5^{1/4}·i). **EXPECT: all pass.**
- **A8** Descent invariant: for the A6 i-dressed clique, solve antilinear J=B∘τ commuting with all four gammas; verify J²=−1; attempt renormalisation J→λJ (λ∈K(i)^×) to J²=+1 and verify failure (norm negative at both real places). **EXPECT: J unique up to scale; J²=−1; renormalisation to +1 impossible.** *(J²=+1 would contradict Paper 204 Thm 1 via descent — if found, stop and report.)*
  - *Choice:* B via homogeneous linear solve (16 unknowns): B commutes with the 3 τ-fixed gammas, anticommutes with the τ-negated (i-dressed) one. J² = B·τ(B); report B·τ(B)=c·I₄ and evaluate c at both real embeddings (√5→±√5).
- **A9** σ flips signature: σ (entrywise) of the A6 i-dressed (1,3) clique is a (3,1) clique (square-sign pattern reversed). **EXPECT: pass.**
- **A10** Ring bookkeeping: N(φ)=−1, N(φ²)=+1, ±φ^{2n} exhausts small norm-one units of ℤ[φ]; ℤ[√5]: φ³=2+√5 fundamental, norm-one subgroup ±φ^{6n}, 9²−5·4²=1 with 9+4√5=φ⁶. **EXPECT: all pass as stated.**
  - *Choice:* norm N(x)=x·σ(x); "exhausts small units" checked by N(±φ^k)=(−1)^k for k=0..8 (norm-one ⟺ k even).

**Deliverable:** one table (item · EXPECT · outcome · deviation verbatim), `results.csv`, and `findings.md`. Scripts + outputs committed; paper .md drafts are NOT committed (unreleased, per Cliff).
