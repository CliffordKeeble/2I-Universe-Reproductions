# PRE-REGISTRATION — F9: is the golden Dirac construction golden at all?

**Paper 226 · Golden Construction · F9**
**Status: pre-registration · committed before `bench_f9.py` emits any verdict · 12 July 2026 · internal**

**Brief of record:** `mr_code_brief_F9_connection_check.md` (Fizz, 12 Jul 2026), from Scout's F1 flag.
Standard: Pattern 75; W-108.

**Why this exists.** Scout's cold F1 flag: Γ_s² = +√5·I, Γ_a² = −√5·I → the Clifford quadratic form is
⟨√5, −√5⟩ over K = ℚ(√5). But ⟨a, −a⟩ is **isotropic** (hyperbolic plane) for any a, so ⟨√5, −√5⟩ ≅ ⟨1, −1⟩
over K and the Clifford algebra is **split M₂(K)** — identical to standard (1,1). Every F1–F8 result came back
real-quadratic-generic; that may be a **diagnosis**, not a rescope: the 191 gamma construction may carry **no
icosahedral (2I) content at all**, while 2I's representation theory (genuinely golden) is a *separate* object we
have been conflating with it.

**Substrate.** Exact symbolic over K = ℚ(√5) and K(i); gammas from `golden_algebra.py`.

---

## F9a — THE CONNECTION CHECK (run first) · registered GENUINELY OPEN, no direction

- **K-DISCONNECT:** no 2I structure in the 191 construction → F1–F8 was a real-quadratic Dirac theory with no
  icosahedral content (explains every genericity finding). **Major negative result; escalate immediately.**
- **K-CONNECT:** 2I embeds / acts / is forced → report **exactly how** (the load-bearing thing never checked).

Compute: (1) confirm Γ² signs; ⟨√5,−√5⟩ isotropic; **exhibit an explicit change of basis over K** to ⟨1,−1⟩;
confirm Cl ≅ M₂(K) split; **report what the isometry does to the two derivative directions** (Scout's caveat: it
mixes them). (2) Does 2I (order 120) embed in the unit/Clifford/**Spin** group over K(i)? Any 2I action on the
2-spinor space the gammas respect (centralizer/normalizer)? (3) Report the mechanism **or its absence** plainly —
"the √5 was placed in the entries by hand and carries no icosahedral content" is a valid, important finding; do not soften.
**On DISCONNECT: escalate to Cliff/Tor/Fizz.**

## F9b — does τ deliver what σ did not? · registered EXPECT: τ delivers (declared lean, discount)

Theorem (Fizz): σ fixes i ⟹ σ is ℚ(i)-linear ⟹ σ cannot implement C (which needs i→−i, antilinear). F8's
symptom (ψ̄→e^{+iα}ψ̄, no U(1)) is *why*. τ (i→−i, √5→√5) IS antilinear → the C-candidate. Build the **τ-adjoint**
ψ̄ := (τψ)ᵀK₂ and re-run F8 §3: (1) ψ→e^{iα}ψ ⟹ ψ̄→e^{−iα}ψ̄? (2) conserved U(1) charge? (σ gave none — the
decider) (3) Hermitian/positive-definite (involution of the second kind, KMRT)? (4) ± frequency signature splits
cleanly (K4)? (5) s=1, **both** mass parities (σ-even, σ-odd) — which does τ prefer? **If τ ALSO fails → no C at
all → Majorana by necessity → escalate.** *(F9b is about the construction, worth running even on DISCONNECT — I will run it and say so.)*

## F9c — the bipartition · registered EXPECT: pass (textbook, for the record)

2I ≅ SL(2,5): EVEN {1,3,3′,4,5} Σdim²=60; ODD {2,2′,4′,6} Σdim²=60; total 120. McKay graph (⊗ defining 2)
bipartite = affine E₈, marks (1,2,3,4,5,6,4,2 arm; 3 branch) Σ=30. σ (√5→−√5) swaps 2↔2′, 3↔3′; fixes 1,4,4′,5,6.

## F9d — the SU(2) objection (kills Fizz's σ=P) · registered EXPECT: σ=P is DEAD, 2/2′ are the two embeddings

(1) (½,0) and (0,½) of SL(2,ℂ) restrict to the **same** SU(2) irrep (the 2 is pseudoreal, 2≅2̄). (2) ⟹ 2,2′ are
**not** L/R; test that they are the two inequivalent 2I↪SU(2) embeddings (pentagon/pentagram, 72°↔144°, exchanged by
the A₅ outer automorphism). (3) Report whether σ implements any parity/chirality on the 191 spinor space.

## Genericity, registered
√5 expected load-bearing nowhere in the *construction* (that is the point); genuinely golden content lives in **2I's
rep theory** (F9c). If F9a finds 2I *does* act on the gammas, that is K-CONNECT — escalate the mechanism.

## Fixed choices (decided now)
1. Change of basis exhibited as an explicit M ∈ GL₂(K) with Mᵀ diag(√5,−√5) M = diag(1,−1); report whether M is
   diagonal (benign rescale) or mixes the ∂₀/∂₁ directions.
2. 2I-embedding decided via: (a) centralizer of {Γ_s,Γ_a} in M₂(K(i)); (b) finite subgroups of the Spin group of a
   (1,1) form = roots of unity μ(K(i)) — compare |2I|=120.
3. τ = entrywise complex conjugation (`golden_algebra.tau`); τ-adjoint conservation via the F8 current-solve method,
   σ-even and σ-odd mass.
4. F9c σ-swap shown via 5-fold character values 2cos(2π/5), 2cos(4π/5) and the φ-valued 3-dim characters.
5. No RNG; exact symbolic.

**Deliverables:** `bench_f9.py`, `results_f9.csv`, `findings_f9.md`, a **scope banner** on the F1–F8 findings
recording whether the arc's object had icosahedral content. Nothing banks until the cold W-108 pass.

🐕☕⬡
