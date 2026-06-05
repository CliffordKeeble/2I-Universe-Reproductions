# Findings — Paper 92, the Hodge-dual reframe of the bridge node

**Reframe under test (Cliff's, CONJECTURED).** The deuteron's bridge node is not
a confined electron but the **Hodge-star dual of the proton, p\* = ∗(proton)**.
The hope was one construction curing two gates: the magnetic-moment gate 1(c)
(a proton-dual should carry a nuclear, not Bohr, moment) and the spin-statistics
W-102 conditional (a *form* is parity-even by construction, deriving the bound
core's even parity with no neutrino-count).

**Verdict: LANDING 3 — the clash bites.** The reframe does not cure 1(c) and
does not derive Gate 4. It reproduces the electron-core reading, whose
three-orders-of-magnitude magnetic-moment problem remains the paper's ceiling —
exactly the brief's honest expectation (landing 2 or 3, not 1).

Pre-registration lock (own SHA, inputs fixed from the corpus before compute):
commit **71bf669**, `paper92-hodge-dual/PRE-REGISTRATION.md`.

---

## Locked inputs (recovered from the corpus, firsthand-cited)

| input | value | source | status |
|---|---|---|---|
| manifold | S³/2I, n=3, **Riemannian** | Paper 121 §1, RESULT 1.1 | DERIVED |
| proton form-degree k | **0** (0-form/H⁰) | Paper 121 §2, Table 2 | CONJECTURED |
| k uniquely pinned? | **No** — proton↔electron swap "chosen, not forced" | Paper 121 §2.1 | — |
| ∗(proton) | **= electron** (3-form/H³), robust under swap | Paper 101 §4.3, Lemma L | DERIVED-given-assignment |

The decisive recovered fact: the corpus's *own* Hodge star is the
proton↔electron pairing — Paper 101 §4.3, *"the Hodge star connects only the
k=0 and k=3 cohomology sectors, so only two modes are particle-bearing (proton,
electron)"*; Lemma L, *"∗ : Ω⁰ → Ω³ … places the proton at l=12 and the electron
at l=42."* So **∗(proton) is the electron**, not a new proton-scale object.

---

## What was computed (`hodge_dual.py`, exact over ℚ(√5), reuses the 2I backbone)

### (1) Base case — the Hodge involution (W-002)
On n=3 Riemannian, ∗∗ on a k-form = (−1)^{k(n−k)+t}, t=0, and k(3−k) ∈ {0,2,2,0}
is always even ⟹ **∗² = +1 for all k** (asserted). The (manifold, signature,
degree) fit together; ∗ is a real involution.

### (2) Gate 1(c) — magnetic moment (dictionary D1) → NOT cured
∗(proton): degree 0 → 3 = the electron. By D1 the node's moment scale is that of
its dual object's inertia:

| quantity | value |
|---|---|
| ∗(proton) = electron → μ_B = (m_p/m_e)·μ_N | ≈ **1836 μ_N** |
| deuteron measured moment | 0.857 μ_N |
| conventional μ_p + μ_n | 0.880 μ_N |
| miss factor | **≈ 2141×** |

The Hodge dual of the proton carries the electron's Bohr-scale moment. The
1836× problem the reframe hoped to dissolve **returns unchanged**. Independently,
because k is not uniquely pinned (the swap is "chosen, not forced"), the brief's
own Q2 rule says 1(c) **cannot be claimed cured** even before the moment number —
a Landing-4 condition sitting underneath Landing 3.

### (3) Gate 4 — does "node is a form" confer integer spin? (dictionary D2) → no
Built 2I (reused from `paper92-spin-statistics`), the vector rep V₃ = the
integer irrep **3**, and its exterior powers (exact):

- Λ⁰=**1**, Λ¹=**3**, Λ²=**3**, Λ³=**1** (det=+1) — verified by Newton's
  identities; the SO(3) identities Λ²(V₃)≅V₃ and Λ³(V₃)≅trivial hold exactly.
- the full form bundle **Λ\*(V₃) = 2·[1] ⊕ 2·[3]** (dim 8), central −1 acts
  as **+1** → entirely integer-spin, **zero spinorial content**.

So **"node is a form ⟹ integer-spin" is TRUE — for genuine forms.** But the
physical proton/electron are spin-½ = the **spinorial irrep 2** (central −1 acts
as −1), which appears in Λ\*(V₃) with **multiplicity 0**. Representing a spin-½
particle "as a form" places it in a bundle with no spinorial content — i.e. it
silently strips the spin-½. The form-degree label (de Rham / cohomology sector)
and the spin label (2I central character) are **orthogonal classifications**;
one cannot be read off the other.

### (4) Dirac–Kähler check (done before calling the clash)
Dirac–Kähler identifies Dirac spinors with the inhomogeneous form space Λ\*
under the **left Clifford action**, *not* the wedge / de-Rham form-degree action
of (3). The physical fermion's spin lives in the left action (spinorial,
parity-odd); the form-degree k is a separate grading (parity-even, as computed).
Dirac–Kähler therefore *relates* forms and spinors but does **not** make a
homogeneous k-form an integer-spin physical fermion — it confirms the two
gradings are distinct. The clash is a real structural boundary, not an artefact
of picking the wrong construction.

---

## Why Landing 3, stated for Mr Adversary

- **The reframe's premise is contradicted by the corpus it invokes.** It posits
  p\* as a distinct proton-scale dual; the corpus's Hodge star makes ∗(proton)
  the electron. There is no third object — on n=3, ∗ maps the 0-form to the
  3-form, and the 3-form is the electron (swap-robust).
- **Gate 4 needs an illegitimate conflation.** "The bound node is a form, hence
  integer-spin" is only true if the node is a *genuine* (bosonic) form. The
  proton/electron are spinors that the Hodge Quartet *labels* by a cohomology
  sector; the label does not change their parity. The form bundle provably
  carries no spinorial content, so the spin-½ cannot come from the form-degree.
- **The corpus describes the node as "an electron core," not p\*.** Paper 121
  §9 (PREDICTION 1) states the deuteron is "a P-e-P standing-wave configuration
  of two protons sharing **an electron core**" — the Hodge Quartet paper itself
  reads the bridge node as the electron, not a proton-dual. And §2.2(ii) calls
  proton and electron "the same handle viewed from opposite ends" (explicit
  Poincaré/Hodge duals); §6 ties the 0-form/3-form pair to μ = m_p/m_e = 1836.
  Every corroboration runs the same way: ∗(proton) is the electron.
- **The Hodge-star facts used are textbook.** ∗ acts on differential forms
  (sections of Λᵏ T\*M — bosonic/tensor objects), *not* spinors; ∗∗ =
  (−1)^{k(n−k)}·s with s=+1 Riemannian, giving +1 for all k on n=3; ∗ induces
  the isomorphism Hᵏ → H^{n−k} (Poincaré duality) on harmonic forms. (Standard
  Hodge theory, e.g. Hodge 1941; the ∗∗ sign and the form-vs-spinor domain are
  the two pillars of the clash, and both are bedrock, not framework-specific.)
- **Status flags.** ∗²=+1, the exterior-power characters, Λ\*'s zero spinorial
  content, and ∗(proton)=electron's moment scale: **DERIVED** (exact). The
  Dirac–Kähler boundary: **STRUCTURAL** (reasoned from the representation facts,
  not a single exact identity). The reframe's failure to cure 1(c) or Gate 4:
  **DERIVED given the locked corpus inputs.**

## Bottom line

The spin lemma's result is **unchanged**: P-e-P consistent ⟺ bound core
parity-even, with neutrino ejection forcing it **conditionally** (W-102). The
Hodge reframe does not lift that conditional and does not cure the magnetic
moment. A clean negative — the construction decided, the feeling did not.

**Lead, not a claim — the live route to Gate 1(c) is the spin-0 core, not p\*.**
Gate 1(c)'s "≈1836 μ_N" problem *presupposes the core carries an unpaired
spin-½* (a Bohr magneton). But the spin lemma's Outcome 1 says the bound core is
**parity-even (spin-0)** — and a spin-0 object has **no intrinsic magnetic
moment at all**. So the Bohr disaster may not arise: not because the moment is
nuclear-scale (the failed Hodge hope), but because a spin-0 core contributes
**zero**, leaving the deuteron's 0.857 μ_N to the two protons. Whether two
protons + a spin-0 core actually reproduce 0.857 μ_N is the **deuteron
magnetic-moment derivation** — still owed, now motivated by Outcome 1, and
distinct from the dead p\* route. Flagged for CinC's call; not started here.

## Reproduce
```
python3 hodge_dual.py     # needs the sibling paper92-spin-statistics dir (2I backbone)
```
Pure-stdlib Python (`fractions`); exact over ℚ(√5); every identity an `assert`.
