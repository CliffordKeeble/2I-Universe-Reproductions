# Pre-registration — Paper 92, the Hodge-dual reframe of the bridge node

**Status:** LOCKED before compute. This commit is the timestamped witness that
the manifold, signature, proton form-degree, and the form↔observable dictionary
were fixed *from the corpus* before the dual was computed.

**Provenance:** CinC brief `BRIEF_paper92_hodge_dual_20260605.md`, 5 June 2026,
relayed by Cliff. Follows the spin-statistics lemma (lock `810e733`,
2I-Universe-Reproductions `c9afdac`), whose result — *P-e-P consistent ⟺ bound
core parity-even* — is the spine this run either improves on or does not.

## The reframe under test (Cliff's, CONJECTURED — make it pass or fail honestly)

The bridge node is not a confined electron but the **Hodge-star dual of the
proton, p\*** = ∗(proton). The hope is two gates cured by one construction:
- **Gate 1(c) — magnetic moment.** A free-electron core would carry ≈ 1836 μ_N
  (Bohr); the deuteron's measured moment is 0.857 μ_N. A *proton-dual* node might
  carry a *nuclear*-magneton moment — **if** the dual inherits proton-scale
  inertia. That "if" is the whole job.
- **Gate 4 — the W-102 conditional.** If the bound node is a **form** (integer-
  spin tensor object), it is parity-even by construction — the iff's condition —
  possibly deriving the parity-even core with no neutrino-count needed.

## LOCKED INPUTS (recovered from the corpus, firsthand-cited, before compute)

**1. Manifold, dimension, signature.**
- M = **S³ / 2I** (Poincaré dodecahedral space), dimension **n = 3**,
  **Riemannian** (the round metric the Hodge star requires; the t=0 spatial
  slice). Corpus: Paper 121 §1, RESULT 1.1 (**DERIVED**) — "S³ has exactly four
  de Rham cohomology groups… (H⁰,H³) b=1, (H¹,H²) b=0"; Betti (1,0,0,1).
- Consequence locked: on n=3 Riemannian, ∗∗ on a k-form = (−1)^{k(n−k)+t},
  t=0, and k(3−k) ∈ {0,2,2,0} is always even ⟹ **∗² = +1 for all k** (∗ is a
  real involution). This is the base-case assert (W-002).

**2. Proton form-degree k.**
- k = **0** (proton = 0-form / H⁰). Corpus: Paper 121 §2, Table 2
  (status **CONJECTURED**): "0-form (point) | H⁰ | Vertex mode | Proton".
- **Recorded caveat (decisive, not cosmetic):** Paper 121 §2.1 states the
  matter-block *which-is-which* is **"chosen, not forced"** — "nothing in the
  cohomology forbids the swap" of proton↔electron between H⁰ and H³. So k is
  **not uniquely pinned**: k ∈ {0, 3}. Per the brief's Q2 rule, an unforced k
  means **Gate 1(c) cannot be *claimed* cured** (Landing-4 condition).
- **Robust regardless of the swap:** the Hodge star pairs each matter particle
  with its Poincaré partner. ∗ : Ω⁰ ↔ Ω³ is the *proton↔electron* pairing.
  Corpus: Paper 101 §4.3 — "the Hodge star connects only the k=0 and k=3
  cohomology sectors, so only two modes are particle-bearing (proton,
  electron)"; Lemma L (open) — "∗ : Ω⁰ → Ω³ … places the proton at l=12 and the
  electron at l=42." So **∗(proton) = electron**, swap or no swap.

**3. The form↔observable dictionary (stated before compute, FLAGGED not asserted).**
- **D1 (degree → moment scale).** A node's magnetic-moment scale is that of the
  *particle the corpus assigns to its form-degree*, because the corpus ties each
  degree to one specific particle with one specific inertia (Paper 121 Table 2):
  0-form ↔ proton ⟹ nuclear scale μ_N; 3-form ↔ electron ⟹ Bohr scale
  μ_B = (m_p/m_e)·μ_N ≈ 1836 μ_N. (Honest reading of the brief's "if it inherits
  proton-scale inertia": moment tracks the *inertia of the dual object*, not the
  word "p\*".)
- **D2 (form-character → spin-statistics).** A genuine homogeneous k-form
  transforms in Λᵏ(V₃) (V₃ = the 3-dim vector rep of 2I), an **integer-spin**
  rep (2I central element −1 acts as +1) ⟹ parity-even. A spinor / Dirac–Kähler
  object transforms in the **spinorial** sector (−1 acts as −1) ⟹ parity-odd.
  Discriminator = the 2I central character, reusing the spin-lemma backbone.

## HARD RULE (permission ≠ requirement, transposed)

- The proton's degree is **recovered, never assigned** to reach the answer. The
  convenient k is never chosen to make μ come out nuclear-scale.
- **"Hodge dual → nuclear-scale moment" must be DERIVED, not read off the name
  p\*.** Naming the node "proton dual" does not give it proton inertia. Earned
  only if the construction *forces* the dual to carry proton-scale inertia.
- The dictionary is stated here and **not bent afterwards**.

## COMPUTE (on the locked inputs)

1. **∗² base-case (W-002):** assert ∗² = +1 for k=0..3 on n=3 Riemannian; if the
   involution sign is wrong the (M, signature, degree) don't fit — fail loudly.
2. **∗(proton):** degree 0 → n−k = 3; identity = electron (corpus). Read the
   **moment scale** via D1 → **Gate 1(c)**.
3. **Form-character (D2):** build 2I (reuse `paper92-spin-statistics`), the V₃
   rep, and Λ⁰…Λ³(V₃); verify each is integer-spin via the central character.
   Decompose the full form bundle Λ\*(V₃) into 2I irreps; check its **spinorial
   content** → whether "node is a form" can legitimately confer integer spin on
   the spin-½ proton/electron → **Gate 4**.
4. **Dirac–Kähler check (before any Landing-3 call):** does representing the
   proton as an inhomogeneous Dirac–Kähler form make proton-as-form legitimate
   and integer-spin, or does the form-grading stay orthogonal to the spin-
   grading? Report the precise boundary.

## FOUR HONEST LANDINGS (state which, with the SHA)

1. **Cure** — p\* nuclear-scale moment AND integer-spin by form character → 1(c)
   cured, Gate 4 derived (W-102 → survived). Only then is the deuteron
   magnetic-moment derivation worth starting.
2. **Half** — p\* integer-spin (Gate 4 derived) but no nuclear-scale moment →
   1(c) still owed.
3. **Clash bites** — proton cannot be consistently a form (spinor/form
   incompatibility); reframe bounded; fall back to electron-core reading, 1(c)
   likely-fatal. Check Dirac–Kähler before calling this.
4. **Underdetermined** — corpus doesn't pin (M, signature, k) uniquely; claim
   nothing, report the missing assignment.

**Honest expectation (stated so it can't bias):** the likely landing is 2 or 3,
not 1. The recovered input ∗(proton) = electron already points away from a
nuclear-scale moment; the hard rule forbids bending it. Cliff's conviction rides
along; the construction decides.

🐕☕⬡
