# Paper 116 Part 2 — the 91-group generalisation: findings

**Pre-registration:** `PRE-REGISTRATION.md`, commit **9227629** (Gate G1), committed
before this sweep ran. **Substrate:** Part 1 merged to main (d08c6f4). Hypotheses are
A1's frozen text; operational parameters are pre-reg §2 (l_max = 300, N = 500, CvM
discriminant, permutation seed 20260706). **Findings only — no causal language;
interpretation is Mr A / CinC / Cliff.**

**Population:** 91 natural spectra (59 cyclic Z₂–Z₆₀, 29 dihedral 2D₂–2D₃₀, 2T/2O/2I).
All 91 reached N = 500 zeros; none flagged short. λ₁ takes five values:
**{8 ×59 (cyclic), 24 ×29 (dihedral), 48 (2T), 80 (2O), 168 (2I)}**. Class split:
**33 GUE, 58 GOE, 0 Poisson**.

---

## H_gap — SUPPORTED (primary), with a structural decomposition

**Primary — Spearman(λ₁, Var(s)), pooled:** ρ = **−0.53**, permutation z = **−5.0**,
p = **0.0001**. Sign stable across all three windowing variants on the exemplars
(equal-N −0.67, equal-coverage −0.81, burn-in −0.64). By the pre-registered rule
(z ≤ −2, sign stable), **H_gap is supported**: across the population, spacing rigidity
tends to increase with gap width.

**Secondary — Spearman(λ₁, ratio):** ρ = −0.02, z = −0.22, p = 0.83 — **null**. The
CvM ratio carries no monotone relation to λ₁ pooled. This **vindicates A1's
pre-registered choice of Var(s) as primary**: the ratio degrades as a continuous
variable on the many borderline/Poisson-like rows, exactly as anticipated.

**Stratified reporting (A1-mandated) shows the pooled effect is between-family, and
the within-family picture is more complex:**

| family | n | λ₁ | Spearman(λ₁, Var) within family | Var(s) range |
|---|---|---|---|---|
| cyclic | 59 | 8 (constant) | N/A — λ₁ constant | 0.09 → 0.46 |
| dihedral | 29 | 24 (constant) | N/A — λ₁ constant | 0.09 → 0.28 |
| polyhedral | 3 | 48, 80, 168 | **+1.00** (reversed) | 0.099 → 0.118 |

- **Between families**, small-gap cyclic/dihedral are looser on average and the
  large-gap polyhedral are tightest — this drives the pooled −0.53.
- **Within the polyhedral family the trend reverses**: 2T (λ₁=48, Var 0.099) <
  2O (80, 0.11) < 2I (168, 0.118). Over this flat 0.10–0.12 band, 2I is marginally
  the *loosest* of the three — the AT6 flattening, now confirmed on the sweep.
- **Within the cyclic family (constant λ₁=8), Var spans 0.09–0.46, driven by group
  order, not gap** — see the observed second axis below.

**Confound partials (rank-partial of λ₁ with Var, controlling):**

| control | partial(λ₁, Var) |
|---|---|
| group order \|Γ\| | **−0.82** |
| spectral density | −0.46 |
| \|Γ\| and density | −0.58 |

The λ₁–Var association survives all three controls (stays negative) — it is not merely
a proxy for group order or spectral density. (Controlling for order it *strengthens*,
to −0.82.)

## Observed second rigidity axis (NOT pre-registered; robustness-checked)

At constant gap (cyclic, all λ₁ = 8), Var(s) increases monotonically with group order:
Z₂–Z₈ ≈ 0.10 (GUE-tight), rising to Z₆₀ ≈ 0.46 (GOE/Poisson-loose); odd orders run
looser than adjacent even orders. Dihedral (all λ₁ = 24) shows the same shape,
0.09 → 0.28. **So at fixed gap, a different variable — group order / spectral density —
drives the full GUE-to-GOE range.**

Robustness (l_max ∈ {150, 300, 600}): the tight small-order result is **stable**
(Z₂: 0.11/0.10/0.09; Z₆: 0.11/0.10/0.09 — GUE throughout), so it is not a truncation
artefact. High-order cyclic is l_max-sensitive (Z₆₀: 0.42/0.46/0.49, GOE→Poisson at
L=600) — reported, not smoothed. The paper's four-group sample used only high-order
cyclic (Z₁₂₀, Z₆₀, both GOE); this second axis is invisible at that sampling and
emerges only across the full population.

## H_correlate — supported in spirit (with the pre-declared caveat)

**Rank-partial(solvable, Var | λ₁) = +0.03** — solvable-vs-non-solvable adds no
predictive power for Var beyond λ₁. 2I sits at the **33rd percentile** of trend
residuals (residual −5.3 on the λ₁→Var rank fit): roughly on the trend, marginally
tighter than its gap predicts, not an outlier.

**32 solvable groups classify GUE-side**, including **2T (λ₁=48) and 2O (λ₁=80)** — so
GUE membership is not gated by non-solvability. This is 116's "correlate, not a cause"
borne out on the population: the algebraic character does not determine the class.
(Fizz's expected coherence dividend — solvable 2O, GUE-side — is confirmed.)

**Pre-declared power caveat (from the pre-registration) holds:** the non-solvable class
is the single group 2I, and "solvable" is rank-collinear with λ₁ (the unique
non-solvable group is the unique maximum-gap group), so the formal partial is
low-power. The substantive evidence for "character is not the cause" is the 32 solvable
GUE-side groups here, plus Part 1's gap surgery (abelian Z₆₀ + imposed gap → GUE;
2I with gap filled → GOE), which breaks the collinearity that the natural-spectrum
sweep cannot.

## Net (numbers only — interpretation deferred)

- H_gap (Var vs λ₁) supported pooled and robust to windowing and to the |Γ|/density
  confounds; but the effect is between-family, reverses within the polyhedral triple
  (flat 0.10–0.12), and coexists with a second, order/density axis that spans the
  whole class range at constant gap.
- H_gap on the ratio is null — Var is the right primary variable (A1 confirmed).
- H_correlate: solvable adds nothing beyond λ₁; many solvable groups are GUE-side.
  Character is a correlate, not the class-determining variable — with the single-group
  power caveat stated in advance.

## Parameters (pre-reg §2, unchanged)

l_max = 300 (uniform, since amplitudes m(l)·λ^{−½} do not decay — the cut-off is a
real regulariser); N = 500 equal-zero-count windowing, t_max per Γ in the summary CSV;
degree-3 unfolding; CvM primary discriminant; permutation N = 10000, seed 20260706;
tie-corrected Spearman; exemplar sensitivity (equal-coverage + 25-zero burn-in).

**CSVs:** `part2_summary.csv` (91 rows), `part2_stats.csv` (the correlations + z),
`part2_sensitivity.csv` (exemplar windowing variants).

*Bootstrap Universe Programme — Paper 116 Part 2 companion. 🐕☕⬡*
