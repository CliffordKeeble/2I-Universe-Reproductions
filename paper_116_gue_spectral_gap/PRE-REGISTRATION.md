# Paper 116 Part 2 — PRE-REGISTRATION (91-group generalisation)

**Committed before the sweep runs (Gate G1).** Substrate: Part 1, merged to main
(commit d08c6f4; pytest 18/18). Population and multiplicity machinery imported from
the verified Paper 118 modules (Gate G5).

**Discipline (Fizz, post-Part-1):** the stake below is the A1 text, frozen 5 Jul
2026 — it *predates* Part 1's findings. Part 1's threshold result (the Z₁₂₀ switch
is robust from gap-8, not gap-4) does **not** reword the hypotheses it will be
tested against. The frozen hypotheses are transcribed verbatim in §1; my mechanical
operationalisation is §2; nothing in §2 relaxes or redirects §1.

---

## 1. Frozen hypotheses (A1, verbatim — not editable by Part-1 knowledge)

**Population:** the 118-CODE-1 Part 2 population (90 solvable + 2I), natural spectra.

**H_gap (A1/Q3):** across all 91 natural spectra, spacing rigidity increases with
gap width — **primary**: Spearman(λ₁, Var(s)) negative; **secondary**: Spearman(λ₁,
GUE/GOE ratio) negative — each with a permutation-test z-score (2I-075: the claim is
separation from null). Var(s) is primary because the GUE/GOE ratio degrades as a
continuous variable on Poisson-like rows (large MSD to both references), and 59 of 91
groups are cyclic with λ₁ = 8, where Poisson-like statistics are a live possibility
(noted as possibility, not prediction).

**H_correlate:** after conditioning on λ₁ (rank-partial correlation),
solvable-vs-non-solvable adds no predictive power for Var(s) (primary; ratio
secondary). This is 116's "correlate, not a cause" made testable on the whole
population — and it is precisely the sentence the 130 erratum wants to lean on.

**Stratified reporting (A1):** all correlations reported pooled *and* within families
({cyclic}, {dihedral}, {polyhedral}) so family-composition effects are visible; 59/91
rows are cyclic and must not silently dominate the pooled statistic.

**Confound note (A1):** |Γ|, λ₁, and spectral density co-vary; report partials on all
three. **No causal language in the findings** — correlations and z-scores only;
interpretation is Mr A / CinC / Cliff.

**Windowing (A1/Q2):** equalise expected zero-count — a common N target per group
(nominal N = 500), per-Γ t_max recorded in the CSV; groups exceeding a compute bound
are flagged, never silently truncated. Sensitivity: the equal-coverage variant and a
25-zero burn-in run on the exemplar subset (2T, 2O, 2I, Z₁₂, Z₆₀, Z₁₂₀, 2D₅, 2D₁₅,
2D₃₀).

## 2. Operationalisation (Mr Code's mechanical choices; frozen here, pre-sweep)

- **Discriminant:** the Part-1 primary — CDF-MSD (Cramér–von Mises style, binning-
  free). ratio = MSD_GUE / MSD_GOE. 3-way class = argmin CDF-MSD over {Poisson, GOE,
  GUE}. Reference moments GOE 0.273 / GUE 0.178 / Poisson 1 (A1/F1).
- **λ₁ / gap:** l_gap = first l ≥ 1 with m(l) > 0; λ₁ = l_gap(l_gap + 2).
- **N (target zeros):** **500** per group (A1's nominal; set with compute in view).
- **l_max (mode cut-off):** **300** for every group — a common angular-momentum reach
  so the different mode *densities* remain a genuine per-group property (the series
  amplitudes m(l)·λ^{−½} do not decay, so the cut-off is a real regulariser and must
  be uniform). Density confound = (active modes, l ≥ 1, m > 0, l ≤ 300) / 300.
- **Unfolding:** degree-3 polynomial to the zero-counting staircase; unit-mean
  spacings. Zeros by sign-change bracketing (64 samples/shortest period) + Brent.
- **Spearman ties:** 59 cyclic share λ₁ = 8, dihedral share λ₁ = 24 — tie-corrected
  Spearman (scipy) and a permutation null that preserves the tie structure.
- **Permutation test:** 10000 permutations of the response, **seed 20260706**;
  report z = (ρ_obs − mean_null)/std_null and two-sided permutation p.
- **Rank-partial correlation:** ranks of each variable, residualise the response and
  the predictor on the control(s) by least squares, correlate residuals.
- **Sensitivity (exemplars 2T, 2O, 2I, Z₁₂, Z₆₀, Z₁₂₀, 2D₅, 2D₁₅, 2D₃₀):**
  (a) equal-coverage variant — common fixed t_max, zero-count allowed to vary;
  (b) 25-zero burn-in — discard the first 25 zeros before forming spacings.
  Both must leave the sign of H_gap unchanged. Z₁₂₀ (order 120) is outside the
  Z₂–Z₆₀ population and is constructed directly from the same cyclic machinery.

## 3. Pre-registered decision rules

- H_gap **supported** iff primary Spearman(λ₁, Var) < 0 with permutation z ≤ −2 (2I-075
  separation from null), sign stable across both sensitivity variants.
- H_correlate **supported** (the "correlate, not a cause" reading) iff the rank-partial
  correlation of the solvable indicator with Var(s) given λ₁ is not distinguishable
  from 0 (|z| < 2). *Power caveat stated in advance:* the non-solvable class is a
  single group (2I), so this test can show 2I sits on the λ₁→Var trend but cannot
  strongly separate; reported as such, no over-claim.
- Any result is reported as it falls. Contradiction of the working direction is a
  finding, not a failure (Pattern 39).

*Bootstrap Universe Programme — Paper 116 Part 2 pre-registration. 🐕☕⬡*
