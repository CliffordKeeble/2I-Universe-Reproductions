# Candidate B Results — 2I ⊕ ⋆2I golden-content quantities

Paper 208 (the E8/icosahedron ghost) · Candidate B reproduction · 10 June 2026

Construction per Dechant 2016, Sec 6. E8 roots = 2I ∪ scaling·(I·2I), the two H4
copies. All inner products computed exactly as `a + τb` over ℤ[τ]; no floating point
until the final comparison line. Grounding check (`candidate_b.py`): 120 icosians,
unit norm, full multiplicative closure verified over 14,400 products — the construction
is the binary icosahedral group, not a 120-point lookalike.

The reported quantities measure golden content via `b` = the τ-coefficient of each
pairwise inner product.

- **Q1**  = Σ|b| over all pairs (total golden content)
- **Q1×** = Σ|b| over cross-sector pairs only (one root from each copy)
- **within** = Q1 − Q1× (same-copy pairs)
- **Q2** = Q1×/within
- **Q3** = histogram of the `b` spectrum

| Quantity | Convention A (τ-scaled 2nd copy) | Convention B (bare) |
|---|---|---|
| pairs | 28680 | 28680 |
| Q1 (Σ\|b\|, raw) | 8940 | 5760 |
| Q1× (cross-sector, raw) | 4080 | 2880 |
| within (Q1−Q1×) | 4860 | 2880 |
| Q2 = Q1×/within | 68/81 ≈ 0.8395 | exactly 1 |
| Q1 per-pair | 0.311715 | 0.200837 |
| Q1× per-pair | 0.142259 | 0.100418 |
| Q1 per-root | 37.25 | 24 |
| Q1× per-root | 17 | 12 |
| Q3 b-spectrum | {−1, −½, 0, +½, +1} | {−½, 0, +½} |

## Verdict

> **NULL against both targets (μ = 6π⁵ = 1836.118…, α⁻¹ = π+π²+4π³ = 137.036…),
> under all three normalisations, in both conventions. Predicted null (Fizz,
> 10 Jun 2026) confirmed. Every quantity is rational (integer or simple fraction);
> the targets are transcendental; a rational structure cannot match a π-target
> except by rounding coincidence. The null is structural, not marginal.**

## Two notes that ride with this result (honesty load-bearing)

**(a) Q2 is convention-contaminated — exploratory only.** The within-sector
denominator is altered by the τ-scaling of the second copy itself: within goes
2880 (B) → 4860 (A) purely because one copy was multiplied by τ, injecting golden
content into within-pairs that isn't structural. So Q2 measures partly the
convention, not the geometry. It came back null so the point is moot here, but Q2
may not be cited as a clean structural quantity. Treat as Q3-grade (exploratory).

**(b) The Q3 half-integer spectrum is a labelled curiosity, NOT a result.**
Convention B's golden content is quantised at {−½, 0, +½} — half-integer or nothing.
This is structurally tidy and CONJECTURED-interesting, unbanked, pending Mr Scout:
*is "golden content of bare 2I⊕2I quantised at ±½" a new observation or a known fact
about the 600-cell?* The histogram is committed as data; the interpretation is
committed nowhere.

## Reproduce

```
python candidate_b.py             # grounding: 120 icosians, closure over 14,400 products
python candidate_b_quantities.py  # both conventions, exact ℤ[τ] quantities + targets
```
