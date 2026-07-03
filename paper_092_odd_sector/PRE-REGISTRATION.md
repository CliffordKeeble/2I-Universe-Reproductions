# Pre-registration — the odd sector: where the actors live (Papers 92 + 174)

**Status:** LOCKED before compute. This commit (own SHA) is the witness that the
manifold, the load-bearing proton-sector *input*, and the dictionary were fixed
before the spectrum and the beats decided anything.

**Provenance:** CinC brief `BRIEF_odd_sector_20260605.md`, 5 June 2026 (evening).
Builds on the spin-statistics lemma (lock `810e733`) and the Hodge-dual run
(Landing 3, lock `71bf669` / compute `1e93f63`). Reuses the 2I backbone +
character table from `paper92-spin-statistics`.

## The frame (what the walk found)

Landing 3 proved the **even/bosonic** sector (de Rham forms; integer irreps
{1,3,3′,4,5}=A₅) carries **zero spinorial content**. Read forward: the spin-½
actors were never in the evens. Paper 174 computed the bosonic spectrum —
Klein's invariants of degrees **12, 20, 30 (all even)**, the trivial-rep Molien
series, and explicitly *kills all odd harmonics by the central −1* (174 §5). The
actors (proton, neutron, the beat) live in the **odd/spinorial** sector
{2, 2′, 4′, 6}, which 174 never computed. **Job: compute the odds.**

Cliff's "reverse Hodge conjugate" = the **Galois conjugate over ℚ(√5)**: in the
character table, **2 ↔ 2′** are √5↦−√5 images. Galois fixes the central character
(−1 is rational), so a spinor's Galois conjugate is a spinor — it never leaves
the odd sector.

## LOCKED INPUTS (before compute)

1. **Manifold / signature:** S³/2I, n=3, Riemannian (as locked in the Hodge run
   `71bf669`; the round metric S³=SU(2)). Nothing in this brief forces otherwise.

2. **The proton's sector — the load-bearing input, RECOVERED not chosen.** The
   Hodge Quartet (Paper 121, Table 2) assigns proton = 0-form/H⁰ = the
   *trivial/integer* rep. The Galois/beat picture needs proton = **spinor 2**.
   This is the same form-vs-spinor split that bit in Landing 3. **Pre-committed:
   the proton is not assigned to the spinor because the picture wants it there.**
   Stage 2 recovers whether Dirac–Kähler *forces* a consistent
   proton=spinor-2-with-0-form-label reading, only *permits* it, or *forbids* it
   (Landings A / B / C). An unforced answer is Landing B and is reported as such.

3. **The dictionary (flagged, not asserted; not bent afterwards):**
   - spinorial irrep ⟹ half-integer spin (central −1 = −1); integer irrep ⟹
     integer spin (−1 = +1). [the spin-lemma backbone]
   - a **beat** (mode × conjugate-mode) is a scalar amplitude modulation ⟹
     integer-spin / no intrinsic dipole — **iff** the beat actually contains the
     trivial [1]. **Pre-registered open question:** *which* conjugate yields the
     scalar — the **Galois conjugate** 2′ (the brief's named operation) or the
     **representation dual** 2* ? These are a-priori different maps; Stage 3a
     computes the decompositions of both 2⊗2′ and 2⊗2* and reports which (if
     either) contains [1]. The scalar is read off the computation, not assumed.

## THE COMPUTATIONS — staged

**Stage 1 (mechanical; `spinorial_spectrum.py`).**
- (1a) Galois structure: 2↔2′ are ℚ(√5)-Galois conjugates; verify Galois fixes
  the central character on every irrep, hence preserves the integer/spinorial
  split; verify the spinorial set {2,2′,4′,6} maps to itself. Clean asserts.
- (1b) The odd-sector spectrum 174 skipped: for every irrep ρ compute
  m_{ρ,k} = ⟨χ_ρ, χ_{Sym^k V₂}⟩ (the multiplicity of ρ in the degree-k SU(2)
  harmonics), k=0..N. Validate the **trivial-rep** series reproduces 174 exactly
  (m₁₂=m₂₀=m₃₀=1, m₆₀=2, odd-k=0; closed form (1−t⁶⁰)/((1−t¹²)(1−t²⁰)(1−t³⁰))).
  Then give, for each **spinorial** ρ∈{2,2′,4′,6}, the first levels and the
  closed-form generating function M_ρ(t)=N_ρ(t)/((1−t¹²)(1−t²⁰)). Verify the
  selection: integer irreps appear only at even k, spinorial only at odd k
  (central-character parity = the spin-lemma bit). Global check
  Σ_ρ dim(ρ)·M_ρ(t)=1/(1−t)².

**Stage 2 (the crux; reasoning anchored on Stage 1 + the Hodge Λ\* result).**
Does Dirac–Kähler **force** proton = spinor 2 while carrying the 0-form de-Rham
label, only **permit** it, or **forbid** it? Anchor: the de-Rham form bundle
(integer) and the spinor bundle (spinorial) are disjoint in 2I-content (Landing
3; reconfirmed by the even/odd split in 1b). DK relates Λ\* to spinors via the
*left Clifford action* — a different grading from the de-Rham degree. Recover
whether that bridge forces the identification. Do not pick the convenient answer.

**Stage 3 (downstream; conditional on Stage 2; `beat_exchange.py`).**
- (3a) The beat as spin-0 dressing: decompose 2⊗2* and 2⊗2′ exactly; report
  which contains [1]; the spin-0 dressing is the trivial component **derived**,
  not named. State plainly which conjugate delivers it.
- (3b) Deuteron vs dineutron (predict-vs-relabel): the clean two-body case is
  p–n binds, n–n does not. Decompose 2⊗2 into symmetric (Sym²) and
  antisymmetric (Λ²) parts and read the spin content; test whether the
  deuteron/dineutron contrast **falls out** of exchange antisymmetry of two
  identical spinor-2 fermions, or must be **fitted** via a separate "beat
  commensurability" condition. State which. (The 3-body extension is noted with
  its subtleties — ³H is β-unstable but *bound*, not unbound.)

## HONEST LANDINGS (state which, with the SHA)

- **A — frame holds:** Stage 2 forces proton=spinor-2; the beat gives a derived
  spin-0 dressing; the contrast predicts deuteron≠dineutron. The big one.
- **B — partial:** spectrum + Galois clean, but Stage 2 leaves the proton sector
  *unforced* (permitted, not required). Galois move structurally apt but
  ungrounded; claim nothing about Gate 1(c).
- **C — crux fails:** DK cannot reconcile 0-form-proton with spinor-2-proton.
- **Sub-verdict (3b), independent:** "fit predicts deuteron≠dineutron" or "only
  relabels it."

**Honest expectation (so it can't bias):** Stage 1 lands cleanly (mechanical).
Stage 2 is genuinely open — not betting on a clean force (likely B). Stage 3 is
the prize and most likely to slip on derived-vs-asserted.

## THE HARD RULE (carried, plus the walk's new clause)

- The proton's sector is **recovered, never chosen** to make the Galois move work.
- The beat's spin-0 character must be **derived** (the dressing *is* the trivial
  component), not named — and the *correct conjugate* identified by computation.
- The deuteron/dineutron contrast must **fall out**, not be fitted.
- **New clause:** many threads land in one place (Landing 3, Galois, beat,
  even/odd, spin-0 core). **Convergence is not derivation.** The more beautifully
  it assembles, the more carefully each piece earns its assert. A thread that
  "obviously fits" still gets its assert.

🐕☕⬡
