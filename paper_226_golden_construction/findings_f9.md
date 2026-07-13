# Findings — F9: is the golden Dirac construction golden at all? **K-DISCONNECT.**

**Paper 226 · Golden Construction · F9**
**12 July 2026 · exact symbolic over K = ℚ(√5) and K(i) · `bench_f9.py` · `results_f9.csv` · `run_f9.log`**
**Pre-registration:** `PRE-REGISTRATION_f9.md` (committed `19c27d6`, before this verdict).
**Brief of record:** `mr_code_brief_F9_connection_check.md` (Fizz, 12 Jul 2026), from Scout's F1 flag.

Status tags: **DERIVED** (exact proof), **STRUCTURAL** (mechanism). Nothing banks until the cold W-108 pass.

---

## Headline (say it loudly — Fizz's instruction)

> **F9a = K-DISCONNECT. The Paper 191 gamma construction carries NO icosahedral (2I) content. The √5 in the gamma entries is a removable (1,1)-hyperbolic artefact.**
>
> Γ_s² = +√5·I, Γ_a² = −√5·I ⟹ the Clifford quadratic form is ⟨√5, −√5⟩, which is **isotropic** (Γ_s+Γ_a is null) — the hyperbolic plane — **≅ ⟨1, −1⟩ over K** (explicit change of basis exhibited). So **Cl ≅ split M₂(K)**, identical to standard (1,1) Minkowski. The Spin group is **abelian** (SO(1,1) = {diag(t,1/t)}), its finite subgroups are just the roots of unity μ(K(i)) = {±1, ±i} (order 4), and the **centralizer of {Γ_s, Γ_a} is scalars**. The non-abelian **2I (order 120)**, a subgroup of the *compact* SU(2) = Spin(3) of a **(3,0)** form, **cannot embed, cannot act on the spinor space respecting the gammas, and does not permute them.** The gammas were **written down with a √5 in the entries**, and that √5 carries no icosahedral content.
>
> **This is a diagnosis, not a rescope.** Every "real-quadratic-generic" finding across F1–F8 was this fact showing through. The fortnight studied a **real-quadratic (1,1) Dirac theory**; 2I's genuinely-golden structure (F9c) is a **separate object** we conflated with it. **ESCALATE to Cliff / Tor / Fizz.**
>
> **But the fortnight was not wasted, and F8's escalation is resolved (F9b):** with the **τ-adjoint** (complex conjugation, the antilinear C), the construction is a **healthy Dirac theory** — conserved U(1) charge, positive-definite Hermitian norm, C = τ. F8's "no conserved charge" was because it used the wrong (σ) adjoint. **σ was never C; τ is.**

---

## F9a — the connection check (DERIVED): DISCONNECT

1. **Scout's F1 confirmed.** Γ_s² = +√5·I, Γ_a² = −√5·I; **(Γ_s+Γ_a)² = 0** — a null vector — so ⟨√5,−√5⟩ is the hyperbolic plane.
2. **Explicit change of basis over K.** M = [[½+√5/10, ½−√5/10],[½−√5/10, ½+√5/10]] ∈ GL₂(K) satisfies **Mᵀ diag(√5,−√5) M = diag(1,−1)**. M is **non-diagonal** ⟹ removing the √5 **mixes the two derivative directions** (a K-linear boost/rotation of t and x), **not** a benign separate rescaling — Scout's caveat, confirmed and reported.
3. **Cl = split M₂(K)** (DERIVED). ω = Γ_sΓ_a = √5·Z, ω² = 5·I; {I, Γ_s, Γ_a, ω} span M₂(K) (dim 4). Identical to the standard (1,1) Clifford algebra — a matrix algebra, not a division algebra; the √5 is invisible to the isomorphism class.
4. **2I does not act** (DERIVED, three ways):
   - **Centralizer** of {Γ_s, Γ_a} in M₂(K(i)) = **scalars λI** (dim 1). Nothing but scalars commutes with the gammas.
   - **Spin(1,1) is abelian** (SO(1,1) = {diag(t,1/t)}); 2I is **non-abelian** of order 120 → **cannot embed**. (Contrast: 2I ⊂ SU(2) = Spin(3) is a compact rotation group of a **(3,0)** form — a different signature.)
   - **Roots of unity** in K(i) = ℚ(√5,i) are {±1, ±i} (order 4); 120 ≫ 4 and non-cyclic → no 2I even set-theoretically.
5. **Mechanism, stated plainly (not softened):** *the √5 was placed in the gamma matrix entries by hand; it is a removable feature of the presentation and carries no icosahedral content.* The 191 construction and 2I's representation theory are **two unrelated programmes** that share the symbol "√5."

## F9b — does τ deliver what σ did not? (DERIVED): τ DELIVERS — resolves F8

The σ-linearity theorem (Fizz): σ fixes i ⟹ σ is ℚ(i)-**linear** ⟹ σ **cannot** implement charge conjugation (which needs i → −i, antilinear). F8's "ψ̄ → e^{+iα}ψ̄, no U(1)" was the *symptom*; this is the *reason*. **τ (i → −i, √5 → √5) is antilinear** — the correct C. With the **τ-adjoint** ψ̄ := (τψ)ᵀK₂:

- **Phase (DERIVED):** ψ → e^{iα}ψ ⟹ ψ̄ → **e^{−iα}ψ̄** — U(1) *is* a symmetry (σ gave e^{+iα}).
- **Conserved U(1) charge (DERIVED):** the mass residual after solving the current is **0**, for **both** mass parities (σ-even m ∈ K *and* σ-odd m = 5^{1/4}μ — because **τ(m) = m**, unlike σ). σ gave **none**; τ gives a genuine conserved current.
- **Positive-definite Hermitian (DERIVED):** the metric K₂Γ_s = **diag(√5, 1)** is positive-definite; (τψ)ᵀK₂Γ_sψ = ψ†diag(√5,1)ψ ≥ 0 — a real **positive-definite Dirac inner product** (involution of the second kind, KMRT sense).
- **Signature (DERIVED):** the τ-norm is **real and positive on both frequencies** — a healthy Dirac theory with **no ghosts**. (The particle/antiparticle ± split is carried by the conserved **charge**, not by an indefinite norm; so K4's "opposite definite signs" premise dissolves — the τ theory is *better* than K4-PASS, it is positive-definite.)

**Verdict: τ delivers. C = τ (complex conjugation), not σ.** This is a **construction** result (valid regardless of F9a), and it **resolves F8's escalation**: F8 correctly found σ gives no charge; the fix is that C was never σ. I ran F9b despite F9a = DISCONNECT, per the brief's §3 allowance — said so here.

## F9c — the bipartition of 2I ≅ SL(2,5) (DERIVED, textbook, for the record)

- **EVEN** {1,3,3′,4,5} Σdim² = **60**; **ODD** {2,2′,4′,6} Σdim² = **60**; total **120 = |2I|**.
- **McKay marks** (1,2,3,4,5,6,4,2 arm; 3 branch) **sum = 30** — the affine **E₈** Dynkin diagram.
- **σ (√5 → −√5) swaps 2 ↔ 2′ and 3 ↔ 3′** (fixes 1,4,4′,5,6): the 2-dim character on the 5-fold classes is 2cos72° = (√5−1)/2 and 2cos144° = −(√5+1)/2, which σ interchanges; the two A₅ 3-dim characters φ and 1−φ likewise. **2I is genuinely golden in its character field** — a *non-removable* Galois invariant — in sharp contrast to the *removable* √5 of the 191 gamma entries. **That contrast is the whole disconnect.**

## F9d — the SU(2) objection: σ = P is DEAD (DERIVED)

- **The SU(2) doublet is pseudoreal:** ε g ε⁻¹ = conj(g) for g ∈ SU(2) (DERIVED) → **2 ≅ 2̄**. So **(½,0) and (0,½) of SL(2,ℂ) restrict to the same SU(2) irrep** — they differ only by boosts.
- Since **2I ⊂ SU(2) is a rotation group**, no rep-theoretic fact about 2I separates left from right → **2 and 2′ are NOT chiralities**; they are the **two inequivalent 2I ↪ SU(2) embeddings** (pentagon/pentagram, 5-fold 72° vs 144°, swapped by the A₅ outer automorphism = σ on characters, F9c).
- **σ implements neither parity nor chirality on the 191 spinor space** (which carries the split (1,1) Clifford algebra, not a 2I action). **σ = P is dead** — Fizz proposed it and killed it the same morning; confirmed flat.

---

## Escalation (§6) and what's owed

**F9a = DISCONNECT is the escalation.** Per the brief this outranks everything: *the object the F1–F8 arc studied (a real-quadratic (1,1) Dirac theory) was not the programme's assumed object (an icosahedral / 2I construction).* This is handed to **Cliff / Tor / Fizz**. The scope banner is added to the F1–F8 findings (`findings_f1_f4.md`, `findings_f5.md`, `findings_f7.md`, `findings_f8.md`) recording that the arc's object carried **no** icosahedral content.

**What survives, stated fairly:**
- The **construction** is a consistent **real-quadratic (1,1) Dirac theory** with **C = τ** (F9b) — F7's real nondegenerate action and F9b's conserved positive-definite charge stand. It is simply **not golden/icosahedral**.
- **2I's golden structure is real** (F9c) but lives in its **representation theory** (character field ℚ(√5), σ-conjugate irreps, affine E₈), which the 191 gamma construction **does not touch**.
- The programme's open task is now explicit: **if** an icosahedral Dirac theory is wanted, 2I must be put into the construction *on purpose* (e.g. gammas built as 2I-equivariant intertwiners, or a (3,0)/SU(2)-signature Clifford algebra where 2I ⊂ Spin), not inherited from a √5 in the entries. **That is a design decision for Cliff/Tor, not a bench call.**

## Reproduce

```
cd paper_226_golden_construction
python bench_f9.py   # -> results_f9.csv; F9a DISCONNECT (split M2(K), abelian Spin(1,1), centralizer scalars);
                     #    F9b tau delivers (conserved U(1), positive-definite); F9c 2I golden in char field; F9d sigma=P dead
```
Environment: `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic; no RNG.

🐕☕⬡
