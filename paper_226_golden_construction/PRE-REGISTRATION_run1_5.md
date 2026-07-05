# Pre-registration — Run #1.5, the place-pair toy

**Spec author:** Fizz, `fizz_reading_bench_run_1_paper_226.md` §6 (pre-registered there, 5 Jul 2026).
**Executor:** Mr Code. This file pins the construction choices; the EXPECT lines are Fizz's, transcribed verbatim. Committed before any verdict output.

**Idea under test.** A single ℚ(√5)-rational (K-rational) toy spinor field, read at the two real places of K (√5↦+√5 in slot 1, √5↦−√5 in slot 2), is claimed to carry the particle at ∞₁ and the antiparticle at ∞₂; σ = the place swap; C₀ = Z⊗(slot swap). This relocates the A5 "coordinate transposition" into "which slot is read."

**Setup (from Run #1, fixed):** 1+1 operator `L_q = Γ_seed D₀ + Γ_adj D₁`, `D_μ = ∂_μ − q a_μ Z`, plane wave ⇒ `M_q(k) = Γ_seed(−i k₀ I − q a₀ Z) + Γ_adj(−i k₁ I − q a₁ Z)` on the amplitude χ. Place evaluation ev₁ = identity (√5→+√5), ev₂ = σ (√5→−√5). C₀ = Z∘σ ⇒ pair action (ψ⁽¹⁾,ψ⁽²⁾) ↦ (Z ψ⁽²⁾, Z ψ⁽¹⁾) = Z⊗swap.

**Method / choices.**
- Decisive results at the **operator level** (exact symbolic 2×2 over ℚ(√5,i), all k): match σ(M_q) and Z·σ(M_q)·Z against the operator family `L(order, charge, asign, ksign)` where order ∈ {literal (Γ_seed on ∂₀), swapped (Γ_adj on ∂₀)}, charge ∈ {+q,−q}, asign ∈ {+,−} (a_μ→asign·a_μ), ksign ∈ {+,−} (k→ksign·k). The unique family member matched names exactly what each slot solves.
- Cross-check every operator equality **numerically at BOTH real places** (mpmath, √5→±√5) per the standing discipline promoted in Fizz §3.
- A concrete on-shell K(i)-rational solution χ (null vector of M_q on det M_q = 0) is exhibited as a witness where a solution is needed.
- Born pairing demonstrated on a symbolic K-rational amplitude a = p + q√5 (192 §5's reduction).

**EXPECT sheet (verbatim, Fizz §6):**
- **B1.5-1** — the K-rational solutions of Run #1's L_q, evaluated at both places, give a slot-pair solving (L_q at ∞₁, L_{−q}-with-reversed-metric at ∞₂) — i.e. the ∞₂ slot natively carries the conjugate dynamics. **EXPECT: pass** (failure kills the place-pair reading).
- **B1.5-2** — C₀ = Z⊗swap maps the solution pair to a solution pair with q→−q **without any coordinate transposition in the pair formulation** (the transposition absorbed into which slot is read). **EXPECT: pass — the make-or-break.** Failure ⇒ metric-activity not absorbable by the place structure ⇒ the fork's relocation was wrong.
- **B1.5-3** — the Born pairing at place level P = ψ⁽¹⁾·ψ⁽²⁾ reproduces the Galois norm N (192 §5's DERIVED reduction) and is C₀-invariant. **EXPECT: pass.**

**Discipline:** pre-registered, verbatim execution, deviations reported loud, no patching. Deliverable: `run_1_5_place_pair.py`, `results_run1_5.csv`, `findings_run1_5.md`.
