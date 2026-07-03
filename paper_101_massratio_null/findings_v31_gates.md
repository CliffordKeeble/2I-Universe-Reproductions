# Findings — Paper 101 v3.1, the two computational gates

**Date:** 11 June 2026
**Brief:** Fizz 🌀 via CinC, "Two computational gates for Paper 101 v3.1" — gates Mr A
objections (1) and (3). Both pre-registered, bounded, exact; commit regardless of outcome.
**Scripts:** `orbit_index_null.py` (Task 1), `vertex_form_divisor.py` (Task 2).
**Verdict (one line each):**
- **Task 1 — PASS (positive).** The §5.3 null survives the change of generative process:
  with corrections drawn from genuine 2I orbit indices, **1848 − 12 = 1836 on m_p/m_e is
  the unique survivor at 0.01%.** Fizz's on-record prediction holds.
- **Task 2 — PASS (positive).** **div(f₁₂) is exactly the 12 vertices:** the vertex form
  vanishes on all 12 vertices (exactly) and is nonzero on the face-centre (2/5) and
  edge-midpoint (1/80) representatives. "One divisor, two jobs" confirmed.

---

## Task 1 — Orbit-index null (Mr A objection 1)

### Step 1 — realised subgroups of 2I (DERIVED, exact, complete)

All subgroups enumerated from the explicit 120-element group via an integer Cayley
table (every subgroup of 2I is 2-generated, so closures of all singletons and pairs are
exhaustive — the "not realised" claims are computed, not asserted).

- **76 distinct subgroups.** Realised orders **{1, 2, 3, 4, 5, 6, 8, 10, 12, 20, 24, 120}**
  (per-order counts 1,1,10,15,6,10,5,6,10,6,5,1).
- Divisors of 120 that are **not** a subgroup order: **{15, 30, 40, 60}** (computed). 2I is
  perfect, so no index-2/3 subgroup (rules out 60, 40); no element of order 15 (rules out
  15, 30).
- **Realised orbit-index set [2I:H] = 120/|H| = {1, 5, 6, 10, 12, 15, 20, 24, 30, 40, 60, 120}**
  (12 distinct indices). All 12 fall inside the old swept range (≤168).

### Step 2/3 — the null with corrections = orbit indices (DERIVED)

Same surviving 2I spectrum (45 levels, l ≤ 120), same 5 mass-ratio targets, same ±, same
windows. Combination count 5400 = 45 × 12 × 2 × 5.

| window | total survivors | detail |
|---|---|---|
| **0.01%** | **1** | `1848 − 12 = 1836` on m_p/m_e (rel 0.0083%) |
| 0.05% | 2 | adds `1848 − 10 = 1838` on **m_n/m_e** (rel 0.0372%) |
| 0.005% | 0 | 1836 exact vs measured 1836.153 falls outside |

> **Bare result: at the primary 0.01% window, exactly ONE survivor overall, and it is
> 1848 − 12 = 1836 on m_p/m_e. 12 is the unique orbit index that hits.**

Note the 0.05% second survivor differs from the flat-set run: the flat set had `1848 − 11
= 1837` (11 ∉ the orbit-index set); the orbit-index family instead produces `1848 − 10 =
1838` on m_n/m_e. Either way the primary-window verdict is unchanged — exactly one.

**Reading against Mr A's charge.** The charge was that "exactly one survivor" might be
measuring spectral sparsity rather than the orbit-index *mechanism*. The test discriminates:
even though all 12 orbit indices lie in the swept range, only one (12) produces a hit, and
it is the same identification. So the verdict does **not** collapse to "any dense correction
set would do this" — the orbit-index family, a genuinely different generative process,
reproduces the unique survivor. **Clause (c) keeps its null support.**

Caveat for the record (does not change the verdict): the survival is jointly driven by the
sparse spectrum *and* by 12 ∈ orbit indices; the test shows the orbit-index restriction does
not *introduce* spurious competitors, which is exactly what objection (1) asked.

---

## Task 2 — div(f₁₂) = exactly the 12 vertices (Mr A objection 3)

### Naming note (flagged for CinC / Mr A)

Klein's f₁₂ is canonically the degree-12 *binary* form on ℂP¹ (a product of 12 linear
forms, one per vertex point). The brief asks to evaluate it at the real vertex
*coordinates* (0,±1,±φ) **exact over ℚ(√5)** — which is the real ℝ³ avatar of f₁₂: the
unique-up-to-scale degree-**6** icosahedral invariant `f_vert` vanishing on the 12 vertices
(Klein deg-12 in z₁,z₂ ↔ deg-6 in x,y,z, via the 2:1 spinor map). I built `f_vert` by the
**Reynolds operator** — averaging a seed monomial over the 60 icosahedral rotations taken
exactly from the 2I quaternions (entries in ℚ(√5)), landing in the 2-dim invariant space
{r⁶, I₆}, then subtracting the r⁶ multiple that makes it vanish at one vertex. Crucially it
is **not** built as a product of vertex factors (that would make vertex-vanishing
tautological) and **not** numerically fitted. Invariance + vanishing at one vertex ⇒
vanishing on the whole orbit; the invariant is unique up to scale, so the seed choice does
not affect the zero/nonzero verdict.

### The checks (DERIVED, exact ℚ(√5))

`f_vert` (seed x⁶) is a genuine degree-6 icosahedral invariant. Evaluations:

- **All 12 vertices `(0,±1,±φ)` & cyclic: `f_vert = 0`, exactly, term by term.** ✓
- **Face-centre (3-fold axis) `(1,1,1)`: `f_vert = 2/5` ≠ 0.** ✓
- **Edge-midpoint (2-fold axis) `(0,0,1)`: `f_vert = 1/80` ≠ 0.** ✓

(The √5 terms cancel at the face/edge points, leaving clean nonzero rationals.)

> **Falsifier NOT triggered.** f₁₂ vanishes on the vertex orbit and on neither the 20
> face-centres nor the 30 edge-midpoints. **div(f₁₂) is exactly the 12 vertices.** Clause
> (c)'s "one divisor, two jobs" correction mechanism is confirmed.

Honest caveat on "divisor": on the sphere a single degree-6 invariant has a *curve* as its
zero locus, not 12 isolated points; "divisor = exactly the 12 vertices" is precise in the
binary-form (ℂP¹) sense, where the degree-12 f₁₂ has exactly 12 roots. The computational
gate Mr A asked for — vanish on vertices, nonzero on the face/edge orbits — is what is
verified here, and it is the content that the correction mechanism rests on.

---

## Status flags

- Realised subgroup orders / orbit indices: **DERIVED** (exhaustive Cayley-table enumeration).
- 1848 − 12 unique survivor under orbit-index corrections at 0.01%: **DERIVED** (exact group
  theory; final mass comparison floating-point, identical to §5.3).
- f₁₂ vanishes on all 12 vertices, nonzero on face/edge reps: **DERIVED** (exact ℚ(√5)).
- f_vert is the deg-6 real avatar of Klein's deg-12 binary f₁₂: **STRUCTURAL** (naming;
  flagged above).
- Whether the transferred curvature *is* the mass correction: **out of scope** (Tor's seam).

## Outputs
- `orbit_index_matches.csv` — survivors at each window for the orbit-index family.

---

# Addendum — the product-divisor gate (Mr A, fifth star)

**Date:** 12 June 2026
**Brief:** Fizz 🌀 via CinC, "Divisor of the PRODUCT f₁₂·f₃₀ (not the factor)."
**Script:** `product_divisor.py`.
**Verdict (one line):** **PASS.** div(f₁₂·f₃₀) = (12 vertices) ⊔ (30 edge-midpoints),
a **disjoint** union; f₁₂ is nonzero on the edges and f₃₀ nonzero on the vertices; the
product is nonzero on the 20 faces; and the mass-channel blocking, being sampled on the
vertex orbit only, remains exactly the 12-dimensional vertex orbit. **One flag** (below)
on the real-avatar construction.

## FLAG — the real ℝ³ "edge form" does not isolate the edges

The brief asked to build f₃₀'s avatar "by the same Reynolds construction Task 2 used for
f₁₂." That construction does **not** transfer to the edges, and this is worth the paper
knowing:

- Task 2's vertex avatar `f_vert` is degree **6**, where the invariant space {r⁶, I₆} is
  2-dimensional — enough freedom to vanish on the vertices and stay nonzero on faces/edges.
- The edge avatar would be the unique degree-**15** real invariant `I₁₅`, which is (up to
  scale) the **product of the 15 two-fold-axis = mirror-plane linear forms**. Every special
  axis lies on icosahedral mirror planes (a vertex on 5, a face on 3, an edge on 2), so
  **I₁₅ vanishes on all three orbits — 12 vertices, 30 edges, AND 20 faces** (verified
  exactly: all zero). It is therefore **not** a clean edge form in ℝ³.

> The clean edge-form divisor is intrinsically a **binary (ℂP¹) statement**, where the
> vertex/edge/face point-sets are distinct and the Klein generators are coprime. The
> mirror-coincidence that spoils the real ℝ³ picture does not occur on ℂP¹.

## The rigorous product divisor (binary Klein generators, exact over ℤ)

The standard Klein generators (integer coefficients) `f₁₂` (vertices), `H₂₀` (faces),
`T₃₀` (edges) were verified to be the genuine icosahedral forms and pairwise coprime:

- **Hierarchy:** `Hessian(f₁₂) = 121·H₂₀`, `Jacobian(f₁₂, H₂₀) = 20·T₃₀`.
- **Icosahedral syzygy:** `H₂₀³ + T₃₀² = 1728·f₁₂⁵` (the classical relation, exact).
- **Squarefree:** f₁₂, H₂₀, T₃₀ have 12, 20, 30 **distinct** roots respectively.
- **Pairwise coprime:** `gcd(f₁₂,T₃₀) = gcd(f₁₂,H₂₀) = gcd(H₂₀,T₃₀) = 1`.

> **Bare result.** Coprimality ⇒ the three root-sets are pairwise **disjoint**. Hence
> `div(f₁₂·T₃₀) = (12 vertices) ⊔ (30 edges)`, with:
> - **12 vertices:** f₁₂ = 0, T₃₀ ≠ 0 → product **0**.
> - **30 edges:** T₃₀ = 0, f₁₂ ≠ 0 → product **0** (this is the edge-vanishing Mr A flagged).
> - **20 faces:** f₁₂ ≠ 0 and T₃₀ ≠ 0 → product **≠ 0**.

## The structural resolution (step 3 — does the edge locus enter the mass channel?)

The vertex-curvature module ρ_vert is supported on the **vertex orbit** (the 12-point
permutation module ℂ[2I/C₁₀] of `findings_minus12_index`). The blocking inner product is

  ⟨f₁₂·f₃₀, ρ_vert element⟩ = Σ over the **12 vertices** of f₁₂(vᵢ)·f₃₀(vᵢ)·cᵢ = **0**,

term by term, because **f₁₂(vᵢ) = 0** at every vertex (Task 2). The product's extra zeros
on the 30 edges are a **disjoint** locus (coprimality, above) that a vertex-supported module
never samples. So the edge-vanishing of f₃₀ does **not** enlarge the blocked module: the
blocked content in the mass channel is exactly the **12-dimensional vertex orbit**.

**Confirmed §5.1 addition (as proposed in the brief):**
> *The product f₁₂·f₃₀ also vanishes on the 30 edge-midpoints via its f₃₀ factor; this
> locus is disjoint from the vertex orbit and is not sampled by the vertex-curvature
> module, so the blocked content in the mass channel remains exactly the 12-dimensional
> vertex orbit.*

## Status flags (addendum)

- I₁₅ vanishes on all three orbits (real avatar fails to isolate edges): **DERIVED** (exact ℚ(√5)) — **flag**.
- Klein hierarchy, syzygy, squarefreeness, pairwise coprimality: **DERIVED** (exact over ℤ).
- div(f₁₂·f₃₀) = vertices ⊔ edges, nonzero on faces: **DERIVED** (from coprimality).
- Mass-channel blocking is vertex-supported ⇒ 12-dim: **STRUCTURAL** (rests on Task 2's
  f₁₂(vᵢ)=0 and the orbit disjointness).

---

# Addendum — the live-cell recount (Mr A's null-arithmetic fix, FIX D)

**Date:** 12 June 2026
**Brief:** Fizz 🌀 via CinC, "The live-cell recount" — gates the fifth star (FIX D).
**Script:** `live_cell_recount.py`. **Output:** `live_cell_recount.csv`.
**Verdict (one line):** **PASS.** Mr A is right that "5,400 cells" flatters the null:
of the five targets only **m_p/m_e** is reachable by an integer-minus-integer value at
0.01%, so **N_live = 1,080** (orbit-index family) of the 5,400-cell grid — the other
4,320 cells sat on four unreachable targets and were dead before the run. Exactly one
live cell fired (1848 − 12 = 1836). One refinement-flag below on the tau.

## Why the cells die: candidates are integers, the window is sub-integer

Every candidate value λ_l ± c is an **integer** (integer eigenvalue ∓ integer correction).
The 0.01% window around a target t is [t(1−w), t(1+w)], a band of half-width t·w. At the
proton scale t·w ≈ 0.18. An integer can land in that band **only if t has an integer
within t·w of it.** Mr A's necessary condition.

## Per-target nearest-integer distances (window 0.01%)

| target | value | nearest int | dist to int | window t·w | integer in window? |
|---|---|---|---|---|---|
| m_p/m_e | 1836.152673 | 1836 | **0.15267** | 0.18362 | **YES** |
| m_n/m_e | 1838.683662 | 1839 | 0.31634 | 0.18387 | no |
| m_μ/m_e | 206.768283 | 207 | 0.23172 | 0.02068 | no |
| m_π/m_e | 273.132440 | 273 | 0.13244 | 0.02731 | no |
| m_τ/m_e | 3477.228280 | 3477 | 0.22828 | **0.34772** | **YES** |

Mr A's envelope reproduced exactly: m_p/m_e ≈ 1836.153 sits 0.153 from 1836, inside its
0.184 window. The neutron sits 0.316 from 1839 — **outside** its window (it only enters at
0.05%, where 1848 − 10 = 1838 hits; that is the 0.05% second survivor of Task 1, dead at
0.01%). Muon and pion sit ≫ window from their nearest integers (their windows are tiny,
~0.02–0.03, because they live near ~200–270). **Dead before the run, exactly as charged.**

## The tau refinement (FLAG — Mr A's envelope is slightly too generous to itself)

The brief's envelope says "only m_p/m_e sits ~0.15 from an integer." That is true for the
*distance*, but liveness is distance-relative-to-window, and **m_τ/m_e's window is wide**
(t·w = 0.348 because t ≈ 3477). So m_τ/m_e **does** have an integer (3477) inside its
0.01% window — it passes Mr A's stated near-integer condition. It is nonetheless **dead**,
for a stronger reason: 3477 is **not a producible candidate.** The surrounding surviving
eigenvalues are λ = 3363, 3480, 3599; reaching 3477 needs a correction of 3, 114 or 122,
none of which is in either family. So the full liveness test (definition (2): *an actual
λ_l ± c lands in window*) kills the tau where the near-integer heuristic would have spared
it. Net effect on the verdict: **none** — the tau was dead either way; the honest reason is
"the in-window integer is unreachable," not "no integer is in window."

This is why the script reports **two** liveness columns:
- **(A) near-integer** (Mr A's necessary condition): live targets = {m_p/m_e, m_τ/m_e}.
- **(B) reachable** (the full condition): live targets = {m_p/m_e} only.
(B) ⊂ (A); the tau is the gap. The recount uses **(B)** — a cell is live iff its target is
actually reachable by a candidate.

## The recounted denominators

| family | grid (full) | cells/target | reachable targets | **N_live** | dead-before-run | fired |
|---|---|---|---|---|---|---|
| flat declared-integer (15) | 6,750 = 45×15×2×5 | 1,350 | {m_p/m_e} | **1,350** | 5,400 | 1 |
| **orbit-index (12)** | **5,400** = 45×12×2×5 | 1,080 | {m_p/m_e} | **1,080** | 4,320 | 1 |

Mr A's "~4/5 of the 5,400 were dead" is confirmed exactly for the orbit-index family:
4,320 / 5,400 = 4/5 dead (the four unreachable targets), 1,080 / 5,400 = 1/5 live (the one
reachable target). One of those 1,080 live cells fired.

## The one-line headline to replace "5,400 cells, one survivor" in §5.3

> Of **1,080 live cells** (those on a target an integer-minus-integer value can reach at
> 0.01% — only m_p/m_e among the five), exactly one fired: 1848 − 12 = 1836 on m_p/m_e.
> The full grid is 5,400 cells, but the other four targets sit no closer to a *reachable*
> integer than 0.23 (neutron, muon, tau) at sub-integer windows, so 4,320 cells were dead
> before the run. "5,400" is the full grid (most cells unreachable), not the live
> denominator.

## What this does NOT change

The decision rule (one expected hit per target) and the verdict (one survivor) stand
unchanged; this corrects the *reported denominator*, not the result. The −12 correction's
structural distinction (Task 1) and the divisor mechanism (Task 2, addendum) are untouched.
The ex-ante frozen-registry run is separate and not affected.

## Status flags (FIX D)

- Per-target nearest-integer distances: **DERIVED** (exact arithmetic).
- m_p/m_e the unique reachable target at 0.01%, both families: **DERIVED** (exhaustive
  enumeration, identical to §5.3 / Task 1).
- N_live = 1,080 (orbit) / 1,350 (flat); 4,320 / 5,400 dead-before-run: **DERIVED**.
- Tau near-integer-but-unreachable (3477 not a producible candidate): **DERIVED** — **flag**
  (refines Mr A's envelope; verdict unchanged).
