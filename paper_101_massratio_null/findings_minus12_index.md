# Findings — The −12 as a 2I rep-theoretic index (Filter-2 test)

**Date:** 10 June 2026
**Brief:** Fizz 🌀 via CinC, "The −12 as a 2I rep-theoretic index (NOT a polynomial
spelling)." Paper 101 Lemma L neighbourhood.
**Script:** `vertex_index_2I.py` (exact ℚ(√5) arithmetic, golden entries symbolic).
**Verdict (one line):** −12 emerges as the 2I-index **[2I : C₁₀] = 120/10 = 12**
(orbit size, **YES**), and that index is **outside** the polynomial null ℤ[V,E,F,D,χ]
(**YES**). The Hodge-transfer kernel is **11**, not 12 — the trivial/volume mode
transfers, the non-trivial complement (3 ⊕ 3′ ⊕ 5) does not.

## What was computed (all exact, nothing asserted)

2I is built as the 120 unit icosian quaternions, components carried as exact
`a + b√5` pairs over ℚ. The conjugation action on the 12 icosahedron vertices
(cyclic permutations of `(0, ±1, ±φ)`), the stabiliser, the 9 conjugacy classes,
and the permutation character are all computed from the explicit group.

## T1 — Vertex stabiliser → C₁₀ (DERIVED)

`|Stab(v₀)| = 10`. Element orders in the stabiliser: `[1, 2, 5, 5, 5, 5, 10, 10,
10, 10]`. It contains an element of order 10 whose powers exhaust the stabiliser, so

> **Stab(v₀) = C₁₀** (cyclic, order 10) — the 5-fold rotation C₅ lifted through 2I's
> double cover. It contains the central element −1, as required (C₅ has odd order, so
> its preimage in 2I is the cyclic C₁₀, not C₅ × C₂).

Hence the vertex orbit has size `|2I| / |Stab| = 120 / 10 = 12 = [2I : C₁₀]`.

## T2 — Decomposition of ρ_vert (DERIVED)

The central element −1 fixes all 12 vertices ⇒ it acts as the **identity**
permutation ⇒ ρ_vert is inflated from `2I/{±1} = A₅`. Only the *vector* irreps
{1, 3, 3′, 4, 5} can appear; every *spinor* irrep {2, 2′, 4′, 6} has multiplicity 0.

Permutation character on `[1A, 2A, 3A, 5A, 5B] = [12, 0, 0, 2, 2]` (fixed vertices:
identity 12; edge-axis 180° 0; face-axis 120° 0; vertex-axis 72° and 144° each 2).
Inner products against the (self-checked, orthonormal) A₅ table give:

| 2I irrep | dim | mult in ρ_vert |
|---|---|---|
| **1** | 1 | **1** |
| 2 | 2 | 0 (spinor) |
| 2′ | 2 | 0 (spinor) |
| **3** | 3 | **1** |
| **3′** | 3 | **1** |
| 4 | 4 | 0 |
| 4′ | 4 | 0 (spinor) |
| **5** | 5 | **1** |
| 6 | 6 | 0 (spinor) |

`Σ dim·mult = 1+3+3+5 = 12` ✓ (reconstructs the full module; dimension exhaustion
independently forces every other multiplicity to 0).

> **ρ_vert = 1 ⊕ 3 ⊕ 3′ ⊕ 5.**
> **T2 bare number: multiplicity of the TRIVIAL rep in ρ_vert = 1.**
> (Confirms Frobenius reciprocity ⟨Ind triv, triv⟩ = ⟨triv, triv⟩_{C₁₀} = 1 by direct
> character inner product, not by citation.)

## T3 — The transfer kernel (DERIVED)

The natural 2I-equivariant map ρ_vert → triv (Hodge ⋆: Ω⁰ → Ω³, the volume/trivial
sector). Since the trivial rep appears exactly once, `dim Hom_{2I}(ρ_vert, triv) = 1`;
the map onto the 1-dimensional volume target is surjective.

> **T3 bare numbers: kernel = 11, cokernel = 0.**

Exactly **one** dimension (the trivial summand) transfers to the volume sector; the
other **11** — the non-trivial complement 3 ⊕ 3′ ⊕ 5 — do not. So of the brief's two
options, it is the former: the part with no image in the trivial target is exactly
`12 − 1 = 11`-dimensional. The full 12 does **not** fail to transfer.

## T4 — The null comparison (STRUCTURAL)

`12 = |2I|/|C₁₀| = 120/10 = [2I : C₁₀]` is a **ratio of group orders** (equivalently,
the size of a transitive 2I-orbit). It is **not** an element of the polynomial ring
ℤ[V,E,F,D,χ] — it is not built from those variables and cannot be re-spelled as
`V + E + 2` etc. Tor's Pattern-75 null enumerates polynomial spellings; it does not
range over a group index. The index is a **different kind of object**, so the −12,
relocated from spelling to group index, is promoted from "suggestive" to "structural."

## Bottom line — what the maths says vs. what 1836 wants

Two clean, mathematically **distinct** numbers fall out:

- **the group index / orbit size = 12** — the unambiguous [2I : C₁₀]. This is the −12
  of μ = 1848 − 12 = 1836, and it is outside the polynomial null.
- **the Hodge-transfer complement (kernel) = 11** — the curvature with no image in the
  volume sector, because the trivial/volume mode transfers.

The character computation does not pick 12 *over* 11 or vice versa for the physics;
it states cleanly that **the orbit index is 12** and **the transfer removes exactly
one dimension, leaving 11**. If the corpus's −12 is "the vertex-orbit / group index,"
the answer is 12 and it beats the null. If the corpus's −12 is meant to be "curvature
that fails to transfer to the volume form," the honest answer is 11. This is a Filter-2
test and both branches are reportable; the seam between them (which number is the mass
correction) is Tor's, not the character table's.

## Status flags

- T1 stabiliser = C₁₀, orbit = 12: **DERIVED** (explicit group, exact).
- T2 trivial multiplicity = 1, ρ_vert = 1⊕3⊕3′⊕5: **DERIVED** (exact inner products).
- T3 kernel = 11, cokernel = 0: **DERIVED**.
- T4 index ∉ ℤ[V,E,F,D,χ]: **STRUCTURAL** (it is a group-order ratio by construction).
- Whether transferred/untransferred curvature *is* the mass correction: **out of scope**
  (Tor's seam; not decided here).
