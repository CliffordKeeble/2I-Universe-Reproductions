# PRE-REGISTRATION — F7: the reality sweep over the aligner family

**Paper 226 · Golden Construction · F7**
**Status: pre-registration · committed before `bench_f7.py` emits any verdict · 11 July 2026 · internal**

**Brief of record:** `mr_code_brief_F7_reality_sweep.md` (Fizz, 11 Jul 2026). Standard: Pattern 75; W-108.

**What changed (brief §0).** The F5-step-2 blocker is **dissolved, not authored.** "Is the aligner σ-legal?"
was only ever a **proxy** for the question that actually fires K1 — **does the pairing give a real,
nondegenerate first-order action?** That is criterion-free and directly benchable. And the gap: F3's
COLLISION was computed **only with A₀ = diag(√5,1)** — the demoted presentation. **The reality sweep has
never been run over the aligner family** A = [[0,b],[c,0]] (F5.1, DERIVED). That is F7.

**Substrate.** Exact symbolic over K(i); gammas + Galois maps from `golden_algebra.py`. Objects per §2.

---

## Fixed inputs (previously benched — do not re-derive)

- Γ_s = [[0,1],[√5,0]], Γ_a = [[0,−1],[√5,0]], Z = diag(1,−1), Z² = +I.
- Dispersion (F1b, DERIVED): M² = m²s²/√5 (A-independent — it comes from the EOM operator, not the pairing).
- Aligner family (F5.1, DERIVED): A = [[0,b],[c,0]], b,c ∈ K(i), det A = −bc; AΓ_s = diag(b√5,c), AΓ_a = diag(b√5,−c).
- **MASS PARITY — AXIOM, pinned by Cliff 11 Jul 2026 (226 §9): ∞₂ is a healthy equal-mass antiparticle**
  ⟹ M²_phys must be **σ-even** ⟹ (since √5 is σ-odd) the Lagrangian **m² is σ-ODD: m² ∈ √5·ℚ.** Not a free
  choice; owned by the antiparticle. All reality tests use it. *Note (mine, to verify in §6):* at a single
  place ∞₁ the mass appears as a real positive coefficient regardless of parity; the parity is a two-place
  statement, tested in §6.

## F7 — registered GENUINELY OPEN, no pre-registered direction

- **K-COLLAPSE:** if **any** aligner A = [[0,b],[c,0]] yields a **real, nondegenerate** first-order action
  (kinetic **and** mass, under an admissible field-reality constraint, mod total derivatives) → **K1 UN-FIRES,
  the F1–F4 no-go is RETRACTED** (it was presentation-specific to A₀). **STOP and report immediately; do not
  re-pose F4** (§7). *This retracts Mr Code's own headline COLLISION. Report it flat (Pattern 39).*
- **K-SURVIVE-AND-UPGRADE:** if **no** A in the family yields a real nondegenerate action → the no-go
  **survives and strengthens** to class-wide over the aligner family; C2 becomes a genuine class-statement.
- **K-INDETERMINATE:** if reality holds only under a constraint the EOM does **not** preserve (re-run F3
  sub-check 2) or only with a **degenerate** symplectic form (F3 sub-check 3) → report as such, score neither.

## §5 — Fizz's warm lean, DECLARED so it can be discounted

Candidate **A = [[0,√5],[1,0]]** (b = √5 σ-odd, c = 1 σ-even) → AΓ_s = diag(5,1), AΓ_a = diag(5,−1) (both
plain-symmetric, rational). **Registered as prime suspect, NOT trusted.** Fizz's three declarations, adopted:
1. Verify her algebra first; if wrong, say so, brief still stands on §3's symbolic sweep.
2. It **fails the §4 frontier** (she found ratios {√5, −√5/5}, no uniform ε). Confirm.
3. **⚠️ HALO FLAG.** A working aligner needing **b σ-odd** rhymes with the pinned σ-odd mass parity. **Rhyme,
   not evidence.** Must NOT enter any paper/summary/handoff as a finding or unity. **Test explicitly whether
   σ-odd b is *required*** (e.g. does a σ-even representative like b=c=1 work?). Log the flag; do not let it
   become a sentence.

## §4 — F7b frontier, run alongside (the named tension)

For the same family, test 226 §5.2's displayed criterion **𝒞·σ(Γ^μ)·𝒞⁻¹ = ε(Γ^μ)ᵀ** with a **uniform scalar
ε** across both gammas, using 𝒞 ∈ {A, σ(A), A⁻¹} (report which). Report frontier verdict and reality verdict
**side by side.** If they **disagree**, **that disagreement is the finding** — the action's pairing form and the
conjugation matrix are different objects and may not agree about σ. **Report, do not adjudicate.**

## §6 — also owed (cheap)

- **σ(M²) confirmation.** With m² σ-odd, compute σ(M²) on u±(k) at both places; **EXPECT: M² σ-even (equal mass
  ∞₁ = ∞₂)** — the pinned axiom delivered by the algebra, not merely asserted. If it fails, the axiom is
  unachievable → major escalation.
- **B lemma written down:** defect(e^{iθ}L) = cos θ·X₀ + sin θ·X₁ (ℝ-linearity); if X₀, X₁ are direction-separated
  (∂₀ vs ∂₁) and independent mod TD, two phase points prove the whole circle. Owed to the paper regardless.

## Genericity watch

√5 expected load-bearing nowhere. **Exception to watch:** σ-parity conditions on b, c may want a σ-odd field
element (every ℚ(√d) has one → still generic). **If a d-dependence appears, escalate — cloud, not footnote.**

## Fixed choices (decided now)

1. Reality criterion = **Im[ev₁(e^{iθ}L)] mod total derivatives = 0**, kinetic (mod TD, antisymmetric-part
   extraction) **and** mass (pointwise), at a common phase θ. Full phase circle via the B-lemma (two points).
2. Constraints swept: {τ, σ-real (q=0), C₀ = Zσ} — same three as F3, apples-to-apples.
3. Mass coefficient m taken **real positive at ∞₁** (place value); σ-odd m² is the two-place statement (§6).
4. s = 1 (healthy) primary; s = i recorded for completeness. No RNG; floats only in §6's two-place read.
5. **Order (§7):** §5 verify → §3 sweep → §4 frontier → §6, in parallel. If K-COLLAPSE, STOP before F4.
6. **Second correction banner owed on `findings_f1_f4.md` whichever way F7 lands.**

**Deliverables:** `bench_f7.py`, `results_f7.csv`, `findings_f7.md`, banner on `findings_f1_f4.md`.
Nothing banks until F7 clears **and** the cold W-108 pass certifies.

🐕☕⬡
