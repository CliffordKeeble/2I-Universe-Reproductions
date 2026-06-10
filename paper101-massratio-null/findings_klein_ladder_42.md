# Findings — f₁₂ on the order-5 fibres, and f₃₀·u_vert → level 42

**Date:** 10 June 2026
**Brief:** Tor ⚛️ via CinC (relayed by Cliff). Paper 101 Lemma L neighbourhood;
companion to `vertex_index_2I.py` / `findings_minus12_index.md`.
**Script:** `klein_ladder_42.py` (exact: Klein forms over ℤ via SymPy; order-5
fibre bridge in exact ℚ(√5) from the 120 icosians).
**Verdict (one line):** f₁₂ (the ℓ=12 invariant) vanishes precisely on the
order-5 exceptional fibres; **f₃₀·u_vert = f₃₀·f₁₂ is the unique level-42
invariant** (a₄₂ = 1), landing exactly on the mass anchor λ = 1848. The
trivial/vertex mode's bare Hodge transfer is the non-mass volume sector; only
the vertex-blind edge form f₃₀ carries it to the mass mode. Ledger 12+30=42 closed.

## Declared interpretation (flagged — this is Mr Code's reading, not Tor's words)

The brief's objects are read against the repo's established vocabulary:

- **"the ℓ=12 invariant" = u_vert = f₁₂**, Klein's degree-12 icosahedral form —
  the unique level-12 2I-invariant (a₁₂ = 1), the H⁰/0-form "vertex mode" that
  Lemma L (`paper92-hodge-dual`: "∗:Ω⁰→Ω³ … places the proton at l=12") sits at
  level 12.
- **"f₃₀·u_vert" = f₃₀·f₁₂**, a product in the 2I-invariant ring ℂ[f₁₂,f₂₀,f₃₀].
  This is the only non-degenerate reading: f₃₀ is *itself* a 2I-invariant, so the
  abstract-rep reading f₃₀ ⊗ triv = triv is empty.
- **"the invariant spectrum"** = the graded invariant ring / Molien series
  `(1−t⁶⁰)/((1−t¹²)(1−t²⁰)(1−t³⁰))`, multiplicities a_ℓ, eigenvalue λ_ℓ = ℓ(ℓ+2)
  (the §4.3 / `spectrum.py` spectrum).

**Seam status — u_vert confirmed.** The `Fizz → Tor` −12/−11 seam brief
(10 June 2026) fixes u_vert independently: "proton-as-0-form and electron-as-
volume-form are both trivial-rep objects", and ⋆:Ω⁰→Ω³ "sends the vertex module
into the volume (trivial) sector … the one transferring dimension is the
proton↔electron channel itself, not a stray curvature unit." So u_vert = the
trivial/0-form vertex mode (the proton), and "the trivial mode's transfer is
non-mass" = that proton↔electron ⋆ channel — exactly the ledger (ii) closes.
The remaining reconstruction is only that "f₃₀·u_vert" means the *ring product*
f₃₀·f₁₂ (today's Tor brief, downstream of the seam brief); if Tor intended a
different operation the computations below still hold verbatim, only the labelling
moves.

## (i) — f₁₂ vanishes on the order-5 exceptional fibres (DERIVED)

Exact ℚ(√5) bridge, all 12 vertices:

- Every vertex has 2I-stabiliser **C₁₀** (order 10) — the order-5 rotation lifted
  through the double cover. All 12 checked.
- Each order-5/order-10 element fixes **exactly one antipodal vertex pair**; the 12
  vertices form a **single 2I-orbit** of size 120/10 = 12 = [2I:C₁₀].

Standard Klein coordinates certify f₁₂ as the genuine invariant and exhibit the
vanishing:

- **Syzygy f₂₀³ + f₃₀² − 1728·f₁₂⁵ = 0** (exact) — f₁₂, f₂₀, f₃₀ are the Klein
  icosahedral invariants, not lookalikes.
- **f₁₂(z=0) = f₁₂(z=∞) = 0**, and {0,∞} = Fix(S), S: z↦e^{2πi/5}z the order-5
  rotation. The remaining 10 roots solve z⁵ = (−11±5√5)/2 = {φ⁻⁵, −φ⁵} (two
  5-cycles under S). So f₁₂'s 12-root divisor **is** the order-5 vertex orbit.

> **(i) verified:** the ℓ=12 invariant vanishes precisely on the order-5 fibres.

## (ii) — f₃₀·u_vert lands on level 42 = the mass anchor (DERIVED)

- **deg(f₃₀·f₁₂) = 30 + 12 = 42.**
- The only generator-monomial f₁₂^a f₂₀^b f₃₀^c of total degree 42 is
  **(a,b,c) = (1,0,1)** — i.e. f₁₂·f₃₀. The syzygy lives at degree 60 and does not
  touch degree 42. Hence the **level-42 invariant space is 1-dimensional, spanned
  by f₁₂·f₃₀**.
- **a₄₂ = 1** confirmed two exact ways: the closed Molien series, and the character
  sum a₄₂ = (1/120)Σ_g U₄₂(Re g) computed exactly over ℚ(√5) from our 120 icosians.
  (Same for a₁₂ = a₃₀ = 1.)
- **λ₄₂ = 42·44 = 1848** — the published mass anchor. μ = 1848 − 12 = **1836**.

**The ledger (closed):**

> **12 (vertex mode f₁₂) + 30 (edge form f₃₀) = 42 (mass mode, λ=1848).**

The trivial/vertex mode sits at ℓ=12. By itself its Hodge transfer ∗:Ω⁰→Ω³ is the
volume/trivial sector — **non-mass** (cf. `findings_minus12_index.md`: the trivial
summand transfers, the 11-dim complement does not). It reaches the mass mode **only**
dressed by f₃₀, and lands on the unique level-42 invariant.

**Mechanism (the box's "vertex-blindness"):** f₃₀ is **vertex-blind** — f₃₀(z=0) =
f₃₀(z=∞) = 1 ≠ 0; its 30 zeros are the edge midpoints, not the vertices. Being
nonzero on the order-5 vertex divisor is exactly what lets it carry the vertex mode
up two levels without re-imposing the vertex zeros as the obstruction.

## (iii) — what this furnishes for the §5.1 box (REPORTED, not edited)

The computation supplies the two box fields the brief names:

- **correction = orbit index** — the −12 is [2I:C₁₀] = 120/10 = 12 (the vertex orbit
  size; outside the polynomial null, per `findings_minus12_index.md` T4).
- **mechanism = vertex-blindness of the target mode** — f₃₀ is nonzero on the order-5
  fibres, so f₃₀·f₁₂ reaches λ=1848 as the unique level-42 invariant.

**Out of Mr Code scope:** the actual §5.1 status change (STRUCTURAL-posit →
CONDITIONAL-on-L′) is a paper edit, and the precise statement of **L′** is CinC's
to write — the present run furnishes its content, it does not define it. Paper 101
v3.1 source is not in this repo; nothing here edits it.

## Status flags

- (i) all 12 vertices C₁₀, single orbit, f₁₂ vanishes on Fix(S): **DERIVED** (exact).
- Syzygy f₂₀³+f₃₀²=1728 f₁₂⁵: **DERIVED** (exact SymPy identity).
- (ii) level-42 invariant 1-dim = ⟨f₁₂·f₃₀⟩, a₄₂=1, λ=1848: **DERIVED** (exact, two
  routes).
- f₃₀ vertex-blind: **DERIVED** (exact).
- Ledger 12+30=42 ↔ μ=1848−12: **STRUCTURAL** (grading identity in the invariant ring).
- Identification of the level-42 invariant *as the physical mass mode*, and the
  status of L′: **out of scope** (Tor's seam; not decided here).

🐕☕⬡
