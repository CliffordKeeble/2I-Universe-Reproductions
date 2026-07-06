# Paper 116 — GUE Statistics from the Spectral Gap: Part 1 findings

**Brief:** 116-CODE-1 Part 1 (Fizz, via CinC, 5 Jul 2026; Amendment A1). Verify the
DERIVED tier of Paper 116 v2.1. Findings only — no interpretation of FRAMEWORK/
OBSERVED claims, no erratum wording.

**Substrate:** multiplicities imported from the verified Paper 118 modules
(`eleven_barrier.py`, `solvable_control_null.py`; Gate G5). Pipeline, discriminants
and reconstruction are new and live in this folder.

**Status:** all seven acceptance tests executed; pytest 18/18 green. Three rows
reproduce the paper's *class* robustly, the abelian rows reproduce under documented
choices, and three of the paper's specific numeric claims do **not** reproduce and
are flagged below (AT3 switch point, AT4 placement universality, AT6 2T variance).

---

## Method (what was built)

- **Series.** `Z_Γ(t) = Σ_{l≥1, m(l)>0} m(l)·λ_l^{−½}·cos(t·log λ_l)`, `λ_l = l(l+2)`.
  l = 0 (λ = 0, log λ undefined) is excluded — this is the paper's "Active l" set.
- **Zeros.** Sign-change bracketing on a grid of 64 samples per shortest period,
  refined by Brent. Determinism verified (no RNG in Part 1).
- **Unfolding.** Degree-3 polynomial fit to the zero-counting staircase; spacings
  normalised to unit mean. Degree ∈ {2,3,4} in the sensitivity grid.
- **Discriminant (A1/Q1).** *Primary* = CDF-MSD (Cramér–von Mises style, binning-
  free): mean squared gap between the empirical spacing CDF and each reference CDF
  at the sorted spacings; `ratio = MSD_GUE / MSD_GOE`, `<1 ⇒ GUE-side`. *Secondary*
  = binned-PDF MSD, bins ∈ {20,40,80}, s ∈ [0,3]. A verdict is ROBUST only when
  primary and secondary agree across the grid. Three-way class (Poisson/GOE/GUE) by
  minimum CDF-MSD.
- **Reference moments (A1/F1).** GOE Var = 4/π − 1 = 0.273, GUE = 3π/8 − 1 = 0.178,
  Poisson = 1 — the corrected values, **not** 116 §2.3's printed 0.405/0.286.
- **Reconstruction (Gate G3, no tuning).** `l_max` is the structural inverse of the
  printed Active-l count; `t_max` grows until the printed zero count is reached
  (reconstruction from a published constraint). The tested statistics (Var, class)
  then fall where they fall. Every inferred parameter carries a ×{½,1,2} sensitivity
  band.
- **Character reconciliation (assertion, not comment).** 116 §2.1 writes χ_l at the
  *rotation angle* α: `sin((l+1)α/2)/sin(α/2)`. The repo stores θ = pπ/q with the
  rotation angle α = 2pπ/q, so **θ = α/2** (half-angle). The pytest
  `test_rotation_angle_half_angle_reconciliation` asserts the two forms agree to
  1e-40 over the 2I classes and l ≤ 24. **[DERIVED]**

---

## AT1 — Multiplicity identity **[DERIVED, PASS]**

Route A (mpmath character sum, ≥50 dps, integrality residual < 1e-30) equals route B
(exact integer: Molien for polyhedral, weight-count for cyclic, half-formula for
dihedral) for all tested groups across the three families (2I, 2O, 2T, Z12, Z120,
2D5, 2D15) to l = 60. The two independent *modules* agree for 2I to l = 120
(`eleven_barrier` Molien vs `solvable_control_null` Molien), and the 11 barrier is
intact: m(10) = 0, m(12) = 1. The m(l) feeding Z_Γ are the repo's, exactly.

## AT2 — Table 1 baseline **[OBSERVED]**

Reconstructed l_max/t_max from the printed (Active-l, zeros); class by primary CvM,
robustness across the 27-cell grid (l_max×{½,1,2} · t×{½,1,2} · degree{2,3,4}) and
the three bin counts.

| Γ | l_max | t_max | zeros | Var(s) | ratio_cvm | class | grid flips | verdict | paper |
|---|---|---|---|---|---|---|---|---|---|
| 2I | 300 | 200.0 | 498 | **0.118** | 0.49 | GUE | 0/27 | **ROBUST-REPRODUCES** | GUE, 0.117 |
| 2O | 296 | 199.9 | 520 | 0.110 | 0.54 | GUE | 0/27 | **ROBUST-REPRODUCES** | GUE, 0.122 |
| Z₁₂₀ | 542 | 185.1 | 231 | 0.369 | 2.93 | GOE | 5/27 | **REPRODUCES-UNDER-DOCUMENTED-CHOICES** | GOE, 0.317 |
| Z₆₀ | 570 | 148.7 | 245 | 0.472 | 2.34 | GOE | 7/27 | **REPRODUCES-UNDER-DOCUMENTED-CHOICES** | GOE, 0.371 |

All four classes reproduce. The non-abelian GUE assignment is robust (0 grid flips,
binned agrees on all bin counts); the abelian GOE assignment is window-sensitive
(5–7 of 27 grid cells flip), consistent with the abelian groups sitting nearer the
GOE/Poisson boundary. 2I's variance reproduces the paper to three places (0.118 vs
0.117); the abelian variances run higher than the paper's (0.369 vs 0.317) but on
the same side of the class boundary — reported against, not gated on.

## AT3 — Impose a gap on Z₁₂₀ (the load-bearing experiment) **[OBSERVED — FLAG]**

Fixed natural Z₁₂₀ window (t_max = 185.1); modes below l_gap zeroed.

| row | λ₁ | zeros | Var(s) | ratio_cvm | class (CvM / binned) | paper |
|---|---|---|---|---|---|---|
| natural | 8 | 231 | 0.369 | 2.93 | GOE / GOE | GOE, 0.317, 2.13 |
| gap l=4 | 24 | 322 | 0.203 | 1.20 | **GOE / GOE** | **GUE**, 0.156, 0.60 |
| gap l=8 | 80 | 384 | 0.133 | 0.22 | GUE / GUE | GUE, 0.096, 0.60 |
| gap l=12 | 168 | 427 | 0.098 | 0.46 | GUE / GUE | GUE, 0.074, 0.71 |

The zero counts track the paper (322≈309, 384≈377, 427≈435): the window is faithful.
The switch GOE→GUE is real, but **the switch point is not l=4**. At gap l=4 both
independent discriminants (CvM and binned, all three bin counts) return **GOE**,
though the variance (0.203) has already dropped between the surmises. The class
switch is robust only at **gap l=8**. The paper's Result 1 — "removing three modes
(gap at l=4) switches Z₁₂₀ to GUE" — is not reproduced under either discriminant;
the reproduced statement is "the gap switches Z₁₂₀ to GUE by l=8, with the variance
falling monotonically from l=4."

## AT4 — Fill low modes of 2I **[OBSERVED — FLAG: placement-dependent]**

All C(11,2)=55 placements of two added modes in l ∈ {1..11}, fixed 2I window:

- **GUE: 15 placements** — both modes in the upper gap (l ≥ 6).
- **GOE: 29 placements** — mid-gap (l ∈ {2..5}).
- **Poisson: 11 placements** — any placement including l = 1 (plus {2,3}).

The outcome is a clean monotone gradient in the *added-mode height*: the lower the
mode, the louder it is (amplitude λ^{−½}) and the more it degrades the statistics
GUE → GOE → Poisson. The paper's Result 2 — "adding **any** two low modes switches
2I to GOE" — is **placement-dependent**, not universal.

The paper's Table 3 "2I+2" row (365 zeros, Var 0.220, GOE) is reproduced **exactly**
by the mid-gap pair **{5, 8}** (365 zeros, Var 0.220) — identifying the likely
unstated placement. The "2I+10" row (324 zeros, Var 0.208, GOE) reproduces in
*class* for upper placements (l=2..11 → GOE, 258 zeros, 0.280) but not in exact
zero-count; l=1..10 collapses to Poisson (176 zeros). The 10-mode placement is
unstated and not uniquely recoverable.

## AT5 — Z₆₀ replication **[OBSERVED, REPRODUCES]**

Fixed natural Z₆₀ window (t_max = 148.7):

| row | λ₁ | zeros | Var(s) | ratio_cvm | class | paper |
|---|---|---|---|---|---|---|
| natural | 8 | 245 | 0.472 | 2.34 | GOE | GOE, 0.371 |
| gap l=4 | 24 | 331 | 0.234 | 0.45 | GUE | GUE, 0.192 |
| gap l=8 | 80 | 397 | 0.139 | 0.39 | GUE | GUE (tightens) |

Reproduces §3.3: natural GOE, gap l=4 → GUE, gap l=8 tightens. Note the asymmetry
with AT3 — **Z₆₀ switches at l=4, Z₁₂₀ only by l=8** — a group-dependence the paper's
"the transition is sharp" framing does not capture.

## AT6 — The 2T question **[the finding]**

Printed Table 4 row: 2T, λ₁ = 24, Var = 0.322, ratio = 0.85, GUE. The correct 2T
gap is at l = 6 (λ₁ = 48); λ₁ = 24 corresponds to l = 4.

| spectrum | λ₁ | Var(s) | ratio_cvm | class |
|---|---|---|---|---|
| correct (gap l=6) | 48 | **0.10** (0.095–0.104 over l_max 100/200/400; 0.10–0.12 over N 80–500) | 0.58 | GUE |
| counterfactual (+mode at l=4) | 24 | 0.15–0.17 | 0.46–0.51 | GUE |
| **printed** | 24 | **0.322** | 0.85 | GUE |

**Outcome: A on the label, plus a separate variance error.** The correct spectrum
reproduces the printed *class* (GUE) and ratio-sign (<1), so the λ₁ = 24 is a
mislabel of the correct λ₁ = 48 (Outcome A). But the printed **Var 0.322 is
reproduced by neither spectrum at any sample size** — the correct spectrum is stable
at ≈0.10, ruling out a small-sample artefact. So the 2T row carries two independent
errors: the λ₁ label (24 → 48) and the variance (0.322 → ≈0.10).

Consequence for §5: with the correct 2T variance (~0.10), the printed polyhedral
hierarchy (Table 4: 2T 0.322 > 2O 0.122 > 2I 0.117) **flattens** — all three binary
polyhedral groups have Var ≈ 0.10–0.12, and 2T is no longer the loose outlier. The
"larger gap → tighter" ordering across 2T/2O/2I is not supported once 2T is
corrected. (Interpretation is Mr A / CinC / Cliff; numbers only here.)

## AT7 — Sub-GUE trend **[OBSERVED, confirmed]**

Within the imposed-gap runs, Var(s) decreases monotonically with gap width:

- Z₁₂₀: (λ₁, Var) = (8, 0.369) → (24, 0.203) → (80, 0.133) → (168, 0.098). Monotone ✓
- Z₆₀: (8, 0.472) → (24, 0.234) → (80, 0.139). Monotone ✓

The within-group gap→tightness trend is robust and monotone. (Distinct from the
*across-polyhedral-group* hierarchy of §5, which AT6 shows does not hold once 2T is
corrected.)

---

## §3 ruling — Paper 118 §12 "power spectrum crossover at l ≈ 1112" — **DECLINE**

Read 116 §2.2 + §4 against Paper 118 §12 (now in hand). 118 §12's only relevant
sentence is the "1100 coincidence" (§12, third open question): *"The Hopf-to-generic
crossover in the power spectrum of Z_2I(t) occurs at l ≈ 1112."* It asserts the
number with **no formula and no definition** of "power spectrum of Z_2I(t)". Paper
116 defines Z_Γ(t) (§2.2) but **no power spectrum**, and 1112 appears nowhere in it;
§8's "sinc² connection" is an explicitly open gesture, not an operational definition.

The crossover statistic is therefore **not derivable from what either text writes**
without inventing a power-spectrum definition — which the no-tuning discipline
forbids. Ruling: **decline again** (unchanged from 118-CODE-1). Path B — a Flint text
edit to 118 §12 — activates. Both texts are cited so the record is complete.

---

## Flags for CinC / Fizz (before the 130 erratum leans on these tables)

1. **AT3 — Z₁₂₀ gap-4 does not switch to GUE.** Under both discriminants it is GOE;
   the switch is robust only at gap-8. The paper's "gap at l=4 switches" (Result 1,
   Table 2 row 2) is a boundary claim not reproduced here. Zero counts confirm the
   window is faithful, so this is a genuine class discrepancy, not a windowing one.
2. **AT4 — the 2I fill switch is placement-dependent, not universal.** 15/55 stay
   GUE, 11/55 go Poisson. The paper's Result 2 ("adding two low modes → GOE") holds
   only for mid-gap placements; {5,8} reproduces Table 3 exactly.
3. **AT6 — the 2T Table 4 row has two errors, not one:** λ₁ (24 → 48) *and* Var
   (0.322 → ≈0.10). Correcting the variance flattens the §5 polyhedral hierarchy.
   This is a third 116 numeric error beyond the λ₁ label and the §2.3 reference
   variances (F1).
4. The core mechanism survives: the gap *does* switch the class and *does* tighten
   the variance monotonically (AT7); 2I/2O are robustly GUE and sub-GUE. The
   "correlate, not a cause" sentence the 130 erratum wants is exactly what Part 2
   tests on the full 91-group population.

## Parameters and decisions (documented, Gate G3)

- Samples per shortest period for bracketing: 64 (zero-count convergence checked).
- l_max = structural inverse of printed Active-l; t_max = grow-to-printed-zero-count.
- Gap surgery: `m(l)=0` for l < l_gap (AT3/AT5), fixed natural window. Fills:
  `m(l)=1` at the added l (AT4), fixed natural window.
- 2T (AT6): equal-N (N=500 nominal) plus l_max ∈ {100,200,400} and N ∈ {80..500}
  sensitivity, since the paper prints no Active-l/zero count for 2T.
- Discriminant reference CDFs built by trapezoid integration of the surmise PDFs on
  a 48001-point grid to s=12 (binning-free; avoids closed-form erf edge cases).
- No random seeds in Part 1 (deterministic). The classifier self-consistency test
  seeds np.random with 20260706.

*Bootstrap Universe Programme — Paper 116 Part 1 companion. 🐕☕⬡*
