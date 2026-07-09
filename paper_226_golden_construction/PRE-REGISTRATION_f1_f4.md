# PRE-REGISTRATION — Benches F1–F4 (Free Golden Dirac Dynamics)

**Paper 226 · Golden Construction · Free-dynamics benches F1–F4**
**Status: pre-registration · committed before `bench_f1_f4.py` produces any verdict output · 9 July 2026 · internal**

**Source of record.** `mr_code_brief_free_golden_dynamics_F1-F4.md` (Tor, 8 Jul 2026, `bench brief · pre-registered · v1.0`),
handed to Mr Code by Cliff. The brief is self-contained; §0's objects are used exactly as written. Provenance
per brief §0/§7: bricks 1–3 are Fizz's build; the F1 dispersion sign, the F1∧F3 collision framing, and the
reduction of the gate to K4 are Tor's builds. Nothing here is banked until merged and cold-checked by Mr Adversary.

**Discipline.** Each item carries an EXPECT line authored by Tor/Fizz; I report outcome-vs-EXPECT. A failed EXPECT
is a result, not an error. F3 sub-checks 1–3 are registered as *genuinely open* — the brief states "no pre-registered
expectation on the direction" for the kinetic-reality gate; I report what it does.

**Pass criterion.** Exact symbolic equality over K(i) = ℚ(√5, i) via sympy (`sqrt(5)`, `I`). No floating point until
the final F4 signature read (gated). Reality criterion, per brief §3–§4: **Im[ev₁ L] modulo total derivatives = 0**
(the mass term, carrying no derivative, must vanish *pointwise*; the kinetic term may vanish *modulo total derivatives*).

**Substrate.** Arithmetic + Galois maps reused from `golden_algebra.py` (self-tested: `python golden_algebra.py` → all PASS).
Γ_s = `G_seed`, Γ_a = `G_adj`, Z = `Z2`; σ (√5→−√5, i fixed) = `sigma`, τ (i→−i, √5 fixed) = `tau`, both entrywise.
Objects specific to this brief — A = diag(√5, 1), A′ = diag(−√5, 1), ψ̄ = (σψ)ᵀA, the plane wave u(k), and the
pairing ⟨·,·⟩ — are defined fresh in `bench_f1_f4.py` per §0.

---

## Registered kills (brief §1, verbatim)

- **K1 (no dynamics):** no Galois dressing yields a nondegenerate first-order real action → 226's dynamical layer dies at the free level.
- **K2 (wrong home):** the only healthy nondegenerate real action forces a non-σ dressing → C = σ loses its dynamical home; Fourth Face rescopes to kinematics.
- **K3 (no split):** the conserved pairing admits no clean ± frequency decomposition → seagull pathology at the free level; compact-or-bust.
- **K4 (wrong map):** the split exists but σ maps +frequency → +frequency → the Feynman–Stückelberg reading dies dynamically; 226 §9 kinematic sentence dies.
- **K5 (dissolution):** the K-valued action resists a single-place energy → "bounded below" must be re-posed per place.
- **Genericity:** √5 is expected load-bearing *nowhere* in F1–F4. Any no-go is a **real-quadratic** no-go (holds for ℚ(√13) identically), a rescope not a golden discovery.

## EXPECT sheet (brief §2–§5, transcribed)

- **F1a** anticommutant of {Γ_s, Γ_a}: **EXPECT** joint anticommutant = span(Z), one-dimensional. (Γ_s alone ⟹ X=[[a,b],[−√5 b,−a]]; adding Γ_a forces b=0.)
- **F1b** M = −iEΓ_s + ikΓ_a + msZ, det M = 0: **EXPECT** M = [[ms, −i(E+k)],[i√5(k−E), −ms]], det M = −m²s² − √5(k²−E²), ⟹ E² = k² + m²s²/√5, M² = m²s²/√5. s=1 healthy (M²>0); s=i tachyonic (M²<0). Convention-free; a FAIL ⟹ §0 transcription error.
- **F2a** transpose intertwining with A=diag(√5,1): **EXPECT** Γ_sᵀ = A Γ_s A⁻¹ (parity +1), Γ_aᵀ = −A Γ_a A⁻¹ (parity −1); with A′ the signs exchange; the pair is diag(+,−) or diag(−,+), **never diag(+,+)**.
- **F2b** neutral antisymmetriser N (N·Γ_s ∝ ε, N·Γ_a ∝ ε, ε=[[0,1],[−1,0]]): **EXPECT** N·Γ_s ∝ ε forces N = λ·diag(1,−1/√5), for which N·Γ_a is symmetric not ∝ ε ⟹ **no constant N antisymmetrises both**; σ-dressing forced; σ-adjoint solves the place-swapped operator.
- **F3a** unconstrained σ-field, mass density m·s·B, B = ψ̄Zψ = √5 N(ψ₁) − N(ψ₂), N(ψᵢ)=pᵢ²−5qᵢ²: **EXPECT** *neither* s∈{1,i} gives Im[ev₁(m s B)]=0 for generic configs. Not a kill — forces F3b.
- **F3b** τ-constraint ev₂(ψ)=τ(ev₁(ψ)): **EXPECT** forces pᵢ∈ℝ (bᵢ=0), qᵢ∈iℝ (cᵢ=0); then N(ψᵢ)=aᵢ²+5dᵢ² real & positive; ev₁(B) real; with s=1 mass term real and healthy ("no collision, τ, Fork B lives").
- **F3 sub-check 1 (the gate, OPEN):** Im[ev₁(ψ̄(Γ_s∂₀+Γ_a∂₁)ψ)] mod total derivatives, under the τ-constraint, s=1. **No pre-registered direction.** If real only under a different (σ-type) constraint → clash → COLLISION.
- **F3 sub-check 2 (OPEN):** golden Dirac evolution (Γ_s∂₀+Γ_a∂₁+msZ)ψ=0 preserves the τ-constraint. If evolution breaks it → τ road not dynamically available.
- **F3 sub-check 3 (OPEN):** kinetic form nondegenerate on the τ-halved configuration space. If it degenerates → K1.
- **F3 reconciliation:** state whether the F3b constraint is exactly "σ = τ on the physical locus" (σ- and τ-Majorana coincide on-shell; Paper 226 §5.3).
- **F3 verdict line:** CLEAR (field-reality condition makes kinetic **and** mass real, nondegenerate, EOM-preserved, at s=1) → hand to F4; or COLLISION (the s/constraint making the action real is incompatible with healthy s=1) → K1/K2/K3 fire, stop, report whether real-quadratic-generic.
- **F4 (GATED — run only if F3 CLEARS at s=1):** F4a ⟨ψ₁,ψ₂⟩=∫(σψ₁)ᵀAΓ_sψ₂ conserved on shell (**EXPECT** pass); F4b sign of pairing on ±frequency u(k) at ∞₁ and its σ-image at ∞₂ (**EXPECT/lean** opposite definite signs; σ maps +freq→−freq, K4 passes); F4c ev₁(S)-unitarity standing alone (K5).

## Fixed construction choices (decided now, reported)

1. **Field d.o.f. count.** §0 lists per component aᵢ,bᵢ,cᵢ,dᵢ (4 real) via pᵢ=aᵢ+bᵢi, qᵢ=cᵢ+dᵢi. The parenthetical "eight real symbols per component" is read as **eight real symbols for the two-spinor** (4 per component). Reported ambiguity, not silently fixed.
2. **Modulo-total-derivatives reduction.** A kinetic density is bilinear in (undifferentiated field f)·(∂_μ g). Its coefficient matrix C^μ splits C = S + Λ (symmetric S, antisymmetric Λ). The symmetric part is ½∂_μ(Σ S_{fg} f g) = a total derivative; the antisymmetric part Σ Λ_{fg} f∂_μg is **not** removable. "Im mod TD" = the antisymmetric-part remainder, reconstructed and simplified per direction μ∈{0,1}. Real mod TD ⟺ remainder = 0.
3. **σ-type constraints tested at the gate** (besides τ, F3b): σ-real (qᵢ=0, ψᵢ∈ℚ(i)) and the C₀ Majorana ψ = Z·σ(ψ) (Paper 226 A3). Reported: which constraint (if any) makes the kinetic term real, and whether it clashes with the mass term's constraint.
4. **Genericity confirmation.** F2a transpose-parity and the sub-check-1 kinetic obstruction are re-run with √5 → √13 (same construction, A=diag(√13,1)) to test whether any no-go is real-quadratic-generic per §1.
5. **s ∈ {1, i}; healthy = s = 1** (F1). ev₁: √5 kept as `sqrt(5)` (positive real). ev₂ = ev₁∘σ. No RNG. No floats except the final F4 sign read (gated).
6. **Deliverables:** `bench_f1_f4.py`, `results_f1_f4.csv` (item · expect · outcome · deviation), `findings_f1_f4.md`. Every numerical claim tagged DERIVED / OBSERVED / STRUCTURAL.

🐕☕⬡
