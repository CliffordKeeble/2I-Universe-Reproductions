# Findings — the odd sector: where the actors live (Papers 92 + 174)

**Frame (from the walk).** Landing 3 proved the even/bosonic sector (de Rham
forms; integer irreps {1,3,3′,4,5}=A₅) carries zero spinorial content — the
spin-½ actors were never in the evens. Paper 174 computed the bosonic spectrum
(Klein invariants of degrees 12, 20, 30; all even k; odd harmonics killed by the
central −1). This run computes the **complementary spinorial spectrum 174
skipped**, tests the **Galois/beat** mechanism, and asks whether the proton can
be the **spinor 2** that the picture needs.

**Verdict: LANDING B (partial).** The spinorial spectrum and Galois structure
are clean and new (Stage 1). But the load-bearing input — the proton as the
spinor 2 — is **not forced** by Dirac–Kähler (Stage 2): it is *permitted and
consistent*, not *required*. So the Galois/beat picture is structurally apt but
ungrounded, and no claim is made on Gate 1(c). Stage 3 added one correction and
one conditional prediction. Pre-reg lock (own SHA): **554a6fa**.

---

## Stage 1 — Galois structure + the spinorial spectrum (`spinorial_spectrum.py`) — CLEAN

**(1a) Galois over ℚ(√5).** The √5↦−√5 automorphism swaps **2↔2′** and **3↔3′**,
fixes {1, 4, 4′, 5, 6} (rational characters). The central character χ(−1) is
rational on *every* irrep, so Galois **preserves the integer/spinorial split**;
the spinorial set {2,2′,4′,6} maps to itself (2↔2′; 4′, 6 fixed). **A spinor's
Galois conjugate is always a spinor** — the odd sector is closed under the
arithmetic symmetry. DERIVED (exact asserts).

**(1b) The spinorial spectrum.** For every irrep ρ, m_{ρ,k} = ⟨χ_ρ, χ_{Sym^k V₂}⟩
(multiplicity of ρ in the degree-k SU(2) harmonics), exact for k=0..200.

- **Validation against Paper 174:** the trivial row reproduces 174 exactly —
  m₁₂=m₂₀=m₃₀=1, m₂₄=m₄₂=1, m₆₀=2, the 15 killed even modes, **all odd k = 0**;
  closed form (1+t³⁰)/((1−t¹²)(1−t²⁰)) = (1−t⁶⁰)/((1−t¹²)(1−t²⁰)(1−t³⁰)).
- **Selection rule (DERIVED):** integer irreps appear **only at even k**,
  spinorial **only at odd k** — the central-character parity, i.e. the spin-lemma
  bit. The "odd sector" is literally the odd harmonic degrees.
- **Closed-form generating functions** M_ρ(t)=N_ρ(t)/((1−t¹²)(1−t²⁰)); the
  spinorial numerators are the **E8 exponents** split between the two 2-dim reps:

  | ρ | block | N_ρ(t) | first level (spin j) |
  |---|---|---|---|
  | 2 | spinorial | t+t¹¹+t¹⁹+t²⁹ | k=1 (j=½) |
  | 2′ | spinorial | t⁷+t¹³+t¹⁷+t²³ | k=7 (j=7/2) |
  | 4′ | spinorial | t³+t⁹+t¹¹+t¹³+t¹⁷+t¹⁹+t²¹+t²⁷ | k=3 (j=3/2) |
  | 6 | spinorial | t⁵+t⁷+t⁹+…+t²⁵ | k=5 (j=5/2) |

  Global check Σ_ρ dim(ρ)·m_{ρ,k} = k+1 (all k). **This is the "where the actors
  live" map:** the proton-mode spinor **2 first appears at k=1** (the
  fundamental harmonic), its Galois partner 2′ at k=7.

---

## Stage 2 — the crux: is the proton forced to be the spinor 2? — UNFORCED (Landing B)

The candidate bridge is **Dirac–Kähler** (inhomogeneous forms Λ\* ↔ spinors via
the left Clifford action). The decisive structural facts:

- The de-Rham/integer sector (even k) and the spinor/spinorial sector (odd k)
  are **disjoint** — proven in Stage 1b (and = Landing 3's Λ\*=2·[1]⊕2·[3], zero
  spinorial content). **No object is simultaneously a pure de-Rham form and a
  spinor.**
- Dirac–Kähler *relates* the two bundles, but via a grading (the left Clifford
  action) **distinct** from the de-Rham form-degree. It therefore **permits** a
  consistent reading — the proton modelled as a spinor field living in the odd
  sector, carrying a cohomology class — but does **not force** that the proton's
  spinor sit in the H⁰ (0-form) slot specifically.
- Independently, Paper 121 §2.1 already states the proton↔H⁰ assignment (vs the
  electron) is **"chosen, not forced."** The input is unforced from both sides.

**There is no genuine conflict (not Landing C):** DK shows forms and spinors are
bridged, so "spinor-proton with a cohomology label" is consistent. **But there
is no force (not Landing A):** nothing requires the proton's spinor to be the
H⁰ 0-form. **Landing B** — the Galois/beat picture is structurally apt but its
foundation (proton = spinor 2) is permitted, not derived. Per the hard rule, the
convenient assignment is **not** adopted to ground the picture. STRUCTURAL.

---

## Stage 3 — the beat and the exchange (`beat_exchange.py`) — one correction, one conditional prediction

**(3a) The beat — which conjugate gives the spin-0 scalar? (a correction.)**
The pre-registered open question, answered by computation:

- 2 ⊗ 2\* = **2 ⊗ 2 = [1] ⊕ [3]** — contains the trivial **[1]** once.
- 2 ⊗ 2′ (Galois) = **[4]** — contains **no** trivial.

So the spin-0 scalar (the |amplitude|² beat-dressing) comes from the
**representation dual** (and 2 is self-dual, 2\*=2), **not** the **Galois
conjugate** the brief named. The spin-0 beat-dressing **exists and is DERIVED**
— but it is the dual beat 2⊗2, while the Galois symmetry (Stage 1) is a separate
arithmetic fact. **The picture conflated two different true operations.** This is
the new hard-rule clause biting: convergence is not derivation. (The spin-0-core
lead for Gate 1(c) is still grounded — by 2⊗2\*⊃[1] — *if* proton=spinor-2; the
"Galois beat" wording was imprecise.)

**(3b) Deuteron vs dineutron — predict or relabel?**
Exact: **Sym²(2) = [3]** (the spin-1 **triplet**, symmetric) and **Λ²(2) = [1]**
(the spin-0 **singlet**, antisymmetric). Two identical spin-½ fermions
(proton/neutron = spinor 2, parity-odd, from the spin lemma) need a totally
antisymmetric state, so in s-wave the spin must be the antisymmetric singlet [1]:

- **p–n** (distinguishable): can occupy the [3] triplet → **binds** (deuteron ³S₁).
- **n–n** (identical): s-wave forced to the [1] singlet → no ³S₁ → **unbound**.

**The contrast FALLS OUT of Fermi exchange + the [3]/[1] split — a PREDICTION
from proton=spinor-2, not a relabel.** No "beat commensurability" condition is
needed; **imposing one would be the relabel.** This is the textbook Pauli reason,
now grounded in the 2I spinor assignment — but **conditional on Stage 2** (which
only permits proton=spinor-2), so the grounding rides on that open input.

*Honest 3-body caveat:* ³H (pnn) is **bound** (β-unstable, not unbound); the
brief's "p-n-n unstable" is β-stability, a different question. The clean, decided
test is the 2-body case.

---

## Status flags

- Galois structure, the spinorial spectrum, the 174 validation, the
  even/odd selection, 2⊗2 / 2⊗2′ / Sym²/Λ² decompositions: **DERIVED** (exact).
- Proton = spinor 2 (the load-bearing input): **UNFORCED** — permitted by
  Dirac–Kähler, not required (and "chosen not forced" already in Paper 121).
- The spin-0 beat-dressing: **DERIVED via the dual** (2⊗2⊃[1]); the Galois-beat
  framing: **corrected** (2⊗2′=[4], no scalar).
- Deuteron≠dineutron from Fermi exchange: **DERIVED conditional on proton=spinor-2**.

## Bottom line

The odd sector is now mapped — a clean, new, 174-complementary result, and the
even/odd "clues vs actors" split is vindicated at the spectral level. But the
mechanism's keystone, proton = spinor 2, is **permitted, not forced** (Landing
B), so the Galois/beat picture and its Gate-1(c) payoff remain ungrounded pending
a real force at Stage 2. Two honest course-corrections fell out: the spin-0 beat
is the **dual** beat, not the Galois beat; and deuteron-vs-dineutron is **Fermi
exchange**, which predicts the contrast without any added "commensurability"
condition. The walk found a beautiful shape; the computation kept the parts that
earned their assert and corrected the two that didn't.

## Reproduce
```
python3 spinorial_spectrum.py    # Stage 1: Galois + spinorial spectrum (vs 174)
python3 beat_exchange.py         # Stage 3: beat decompositions + exchange
```
Pure-stdlib Python; exact over ℚ(√5); every identity an `assert`. Reuses the 2I
backbone from `../paper92-spin-statistics`.
