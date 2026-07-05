# Pre-registration — Run #1.6 (C_phys, null-line lemma, factorisation frontier)

**Spec author:** Fizz, `fizz_reading_run_1_5_paper_226.md` §4 (pre-registered there, 5 Jul 2026). **Executor:** Mr Code. EXPECT lines are Fizz's, verbatim; this file pins construction choices. Committed before any verdict output.

**The composite under test (Fizz §2):** `C_phys = C₀ ∘ (x⁰↔x¹) ∘ (a_μ → −a_μ)`, C₀ = Z∘σ. On a plane-wave solution `(χ; k₀,k₁,a₀,a₁)` of L_q it acts as amplitude χ ↦ Zσ(χ), coordinate swap k₀↔k₁ **and** a₀↔a₁ (both are μ-indexed), and a→−a; net image `(Zσχ; k₁,k₀,−a₁,−a₀)` at charge −q.

**Operator family:** `M(k₀,k₁,a₀,a₁,charge) = Γ_seed(−i k₀ − charge·a₀·Z) + Γ_adj(−i k₁ − charge·a₁·Z)` (literal order). Exact 2×2/4×4 over ℚ(√5,i); every operator equality cross-checked numerically at **both** real places (standing discipline, Fizz §3). One claim per EXPECT line (per Fizz §0 note).

**EXPECT sheet (verbatim, Fizz §4):**
- **B1.6-1** — C_phys, end-to-end as the composite, maps every L_q solution to an L_{−q} solution (same operator family, literal order, no residual). **EXPECT: pass.**
  - *Test:* verify `M(k₁,k₀,−a₁,−a₀,−q) = Z·σ(M(k₀,k₁,a₀,a₁,q))·Z` exactly ⇒ Zσχ solves the literal −q operator.
- **B1.6-2** — C_phys² = identity on solutions (incl. the coordinate factor: transposition² = 1; check the factors' commutation exactly). **EXPECT: pass.**
  - *Test:* implement C_phys as a map on `(χ,k₀,k₁,a₀,a₁,charge)`, apply twice, verify identity; check P₀₁²=A²=C₀²=id and the pairwise commutation of the three factors.
- **B1.6-3** — the transposition equals the reflection x⁻→−x⁻ with x⁺ fixed (x± = x⁰±x¹); verify on coordinates and on the induced action on (Γ_seed, Γ_adj). **EXPECT: pass.**
  - *Test:* x⁺,x⁻ under x⁰↔x¹; lightcone gammas Γ± = Γ_seed±Γ_adj (the null-ray nilpotents, analysis §1); operator L = Γ₊∂₊ + Γ₋∂₋; reflection ∂₋→−∂₋.
- **B1.6-4 (factorisation frontier)** — over K(i) (the extension), a linear / τ-borne charge conjugation *does* exist (the standard mechanism revives), in contrast to the K-rational skeleton where run #1 proved only the σ-borne composite exists. **EXPECT: pass** — a fail would mean §3's "extension buys the factors" overreaches.
  - *Choices:* (a) i-coupled 1+1 toy `M_q^i = Γ_seed(−i(k₀+q a₀)) + Γ_adj(−i(k₁+q a₁))` (the T_τ / honest-U(1) coupling, i.e. the i-dressing of the *coupling*): solve for a τ-borne C = B∘τ (B a constant matrix, τ = complex conjugation) intertwining M_q^i with M_{−q}^i for all k (B τ(M_q^i) = M_{−q}^i B); report existence (nonzero B). (b) 4d i-dressed (1,3) clique from A6: solve for a linear and a τ-borne intertwiner C with C γ_μ C⁻¹ = −γ_μ; report existence and whether it needs i (lives in M₄(K(i))\M₄(K)). Contrast: the same non-σ conjugation is unavailable at the all-real K skeleton (run #1 theorem).

**Deliverable:** `run_1_6.py`, `results_run1_6.csv`, `findings_run1_6.md`. Verbatim execution, deviations loud, no patching.
