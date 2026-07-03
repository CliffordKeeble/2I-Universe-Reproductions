# Findings — Paper 92 spin-statistics lemma

**Claim under test (the gate on Paper 92).** Paper 92 models the deuteron as a
Proton–electron core–Proton (P-e-P) standing wave, with the "neutron" being a
proton whose 108-bridge is occupied by a shared electron core — i.e. *the
neutron is not a separate particle.* That picture is the pre-1932
Rutherford/Pauli proton-electron nuclear model, which was abandoned over the
**nitrogen-14 spin-statistics catastrophe**. This investigation asks whether
the icosahedral / 2I standing-wave reading evades that objection or repeats it.

**Verdict: OUTCOME 1 (conditional).** The P-e-P neutron **survives** the
spin-statistics gate, by a single mechanism that cures the neutron, the
deuteron, and nitrogen-14 together — *conditional on one independently-motivated
model commitment*, which is named below for adversary review.

Pre-registration lock: `bootstrap-universe/deuteron-spin-statistics/PRE_REGISTRATION.md`,
commit **810e733** (frozen before any compute). A copy is in this directory.

---

## What was computed

Two standalone scripts, exact arithmetic, no floating point, no fitting.

### 1. `twoI_character_table.py` — the DERIVED backbone

Builds the binary icosahedral group **2I** (order 120, the double cover of the
icosahedral group I ≅ A₅) explicitly as 120 unit quaternions over ℤ[φ], by
closure from the binary tetrahedral group plus one golden icosian; verifies
|G| = 120 and full closure. Computes the **complete character table** (9 irreps)
exactly over ℚ(√5) via symmetric powers of the defining 2-dim representation
(Chebyshev Uₙ) plus Galois conjugation, and **verifies orthonormality exactly**.

Results (all machine-verified, exact):

| quantity | value | status |
|---|---|---|
| group order / closure | 120, closed | **DERIVED** |
| conjugacy classes (sizes) | 9: [1, 1, 20, 30, 12, 12, 20, 12, 12] | **DERIVED** |
| irrep dimensions | 1, 2, 2′, 3, 3′, 4, 4′, 5, 6 (Σ² = 120) | **DERIVED** |
| orthonormality | exact, all 45 inner products | **DERIVED** |
| **INTEGER block** (−1 ↦ +1) | dims {1, 3, 3′, 4, 5} = the A₅ irreps | **DERIVED** |
| **SPINORIAL block** (−1 ↦ −1) | dims {2, 2′, 4′, 6} | **DERIVED** |

2I is exactly the object that distinguishes integer from half-integer spin: the
central element −1 (the 2π rotation) acts as +1 on integer-spin irreps and −1 on
spinorial ones. The free electron is the fundamental spinor **2** (j = ½),
parity-odd.

Two backbone facts, both proved in-script:

- **(a) Parity product rule.** −1 acts on a tensor product as the product of its
  factors, so a composite's integer/spinorial bit = **(−1)^(number of spin-½
  constituents)**. Verified for all 81 irrep pairs.
- **(b) Central-character conservation under induction.** Inducing a spinorial
  (−1 ↦ −1) representation up to 2I yields **only** spinorial irreps —
  demonstrated exactly by inducing the sign rep of the centre ⟨−1⟩, whose
  decomposition has multiplicity 0 on every integer irrep. **A construction
  cannot flip the bit.**

### The *iff* (construction-free, DERIVED)

Because the bit is just spinor-count parity, and the proton is a spinor (one odd
factor each), the observed spins force one and the same condition:

| composite | constituents | observed spin | ⇒ bound core must be |
|---|---|---|---|
| bare neutron | p + core | ½ (odd) | **parity-even (integer-spin)** |
| deuteron | p + core + p | 1 (even) | **parity-even (integer-spin)** |
| nitrogen-14 | 14 p + 7 core | boson (even) | **parity-even (integer-spin)** |

> **The P-e-P model is spin-statistics-consistent if and only if the bound
> electron core is integer-spin (parity-even).** This is rigorous and
> independent of how the bridge mode is constructed (backbone b). All three
> stand or fall together on one bit.

This also resolved the construction question raised in review: the two
candidate constructions (induced-from-stabiliser vs Clebsch–Gordan) **provably
agree on the bit** by central-character conservation, so the verdict cannot be
an artifact of a construction choice. The cross-check is an audit of the
bookkeeping, not a degree of freedom.

### 2. `composite_parity.py` — the hinge and the verdict

Tests the **named forcing mechanism**: neutrino ejection, recovered from Paper
92 §3.1 ("the extended electron field — the neutrino — is ejected at
formation"), which is in the model for *independent energetic reasons*
(the core contracting from Bohr to Compton scale), not invented to dodge N-14.

- **Reading A — no mechanism** (bound core = the electron, a spinor): all three
  composites are **CONTRADICTED** (neutron predicted boson, deuteron/N-14
  predicted fermion). The naive P-e-P count reproduces the 1932 catastrophe
  exactly.
- **Reading B — neutrino ejection.** By fermion-number conservation,
  electron(odd) = neutrino(odd) ⊗ remnant ⟹ **remnant is parity-even**. With a
  parity-even bound core, all three composites match observation:
  neutron ½, deuteron 1, N-14 boson — **all consistent**.
- **Conservation audit.** Both `p + p + e → d + ν` (formation) and
  `n → p + e + ν̄` (β-decay) conserve total fermion parity with the
  parity-even bound core. The even-core reading is the conservation-consistent
  one.

---

## The single load-bearing assumption (flagged for Mr Adversary)

Per the pre-registered **hard rule (permission ≠ requirement):** a single-valued
*spatial* standing-wave mode only *permits* even parity — a bound electron keeps
its intrinsic spin-½ however its spatial mode is confined. That alone would be
**outcome 3**, not 1. What lifts it to outcome 1 is that neutrino ejection
**removes one real spin-½ factor** from the bound count, so conservation
*forces* the remnant parity-even.

That forcing is **conditional on exactly one spin-½ neutrino being ejected per
electron core.** This is the model's standing commitment (Paper 92 §3.1 + the
β-decay phenomenology of one antineutrino per neutron), independent of the
spin-statistics question — which is the bar the pre-reg set for a "requirement-
type" mechanism. **But it is the load-bearing assumption**, and the verdict is
honestly conditional on it:

- if formation could eject the extended field as a **spin-0** object, or as a
  number of neutrinos other than one per core, the forcing weakens to mere
  permission and the honest verdict **drops to OUTCOME 3 (conditional)**.

**Status flags:**
- the 2I representation theory, the integer/spinorial classification, the parity
  rule, central-character conservation, and the *iff*: **DERIVED**.
- "neutrino ejection forces the bound core parity-even": **DERIVED, conditional**
  on the §3.1 one-neutrino-per-core ejection commitment.
- "the §3.1 ejection necessarily carries exactly one spin-½ factor out":
  **OBSERVED/STRUCTURAL** in the model, **not independently derived here** — the
  isolated claim a referee should press.

## Why this is the strong result (and where it could still fail)

The same neutrino ejection that Paper 92 introduced for energetic/contraction
reasons is exactly what rescues the spin-statistics count — without being tuned
for it, and curing the neutron, the deuteron, and nitrogen-14 with one stroke.
A mechanism doing double duty it was not fitted for is the signature of a real
one. The objection that buried the proton-electron neutron for 90 years becomes,
under this reading, the headline: *the geometry evades nitrogen-14 because the
electron's spinor factor leaves with the neutrino.*

It is not unconditional. The whole result rests on the one bit above, and the
honest failure mode is explicit: show that formation need not eject exactly one
spin-½ neutrino per core, and outcome 1 reverts to outcome 3. The *iff* itself
(all three coupled to one bit) is unconditional either way.

## Reproduce

```
python3 twoI_character_table.py    # backbone + iff; writes twoI_character_table.csv
python3 composite_parity.py        # hinge test + verdict
```

Pure-stdlib Python (`fractions`, `csv`); no numpy/sympy required. Every claim is
an `assert` in the scripts — they fail loudly if any exact identity breaks.
