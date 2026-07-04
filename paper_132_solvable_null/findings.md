# Findings — Solvable-Case Null (Withholding Test)

**Pattern-75 falsifier. Design locked by Flint's brief (3 Jul 2026); pre-registration at
`PRE-REGISTRATION.md`. Reproduction: Mr Code, 4 Jul 2026.** `python solvable_null.py` → **14/14 PASS**.

## Verdict (one line)

The framework's weight-2 selection machinery **withholds on genus-0 conics** (β₁ = 0, only the
Hasse–Minkowski yes/no) and **selects on the genus-1 grammar curves** (Nagao sum, ranks separated
`r0 < r1 < r2`). The selection is **not vacuous** — it does not manufacture a rank where there is
no β₁ topology. That is the one thing this test establishes.

## Premise gate (verified before results — brief's hard requirement)

All three premises hold; two bugs were caught **at the gate, before the acceptance code**:

1. **genus 0 ⇒ β₁ = 0** — P¹(ℂ)=S², Betti (1,0,1); Jacobian a point, no positive-rank MW. **[DERIVED]**
2. **Hasse–Minkowski decides conic solubility** — Serre IV; only ∞, 2, odd p|2abc obstruct. **[DERIVED]**
   - *Bug caught:* a first Hilbert-symbol implementation stripped the sign of the unit part,
     giving reciprocity-violating odd-count obstruction sets. Fixed (signed valuation).
   - *Brief correction (Flint's, owned):* conic B `x²+y²=3` was described as failing "ℚ₃ / ℝ." It
     is **soluble over ℝ**; the obstruction is 2-adic + 3-adic, `{2,3}`.
3. **Mestre–Nagao tracks rank.** **[OBSERVED]** Nagao's conjecture is *pending Scout ticket S1* —
   stated as the statistic, not cited as settled.

## Selection — genus-1 grammar curves (β₁ = 2)

`S(X) = (1/log X)·Σ_{good p≤X}(−a_p log p)/p`, a_p by point counting from Weierstrass coefficients
only. C# (Icosian acceptance, since retired) and this Python agree to 6 digits at X = 10⁴.

| curve | rank | **S(10⁴)** (frozen oracle) | S(10⁵) | dist. to rank: @10⁴ → @10⁵ |
|---|---|---|---|---|
| 11a1  | 0 | **−0.192472** | −0.269805 | 0.19 → 0.27 |
| 37a1  | 1 | **+0.823544** | +0.747900 | 0.18 → 0.25 |
| 389a1 | 2 | **+2.056365** | +1.809933 | 0.06 → 0.19 |

- **Separation `r0 < r1 < r2` holds at both X** — the asserted claim.
- X=10⁴ values are **Pattern-75 regression oracles** (reproduce-this-number), **not targets**.
  They **retire the manifest's from-memory quotes** (−0.26 / +0.75 / +2.06). Provenance:
  *reproduced by Mr Code, 4 Jul 2026, X=10⁴, from Weierstrass coefficients.*

### ⚠ Flag — the X=10⁵ trend did NOT confirm convergence-to-rank

The working expectation was a trend *toward* the integer rank (11a1→0, 37a1→1, 389a1→near 2). The
data shows the **opposite direction over this decade**: all three **decreased**, moving *away* from
their ranks — 389a1 from a near-exact 2.056 (0.06 off) down to 1.810 (0.19 off). The values stay
within ~0.27 of the ranks and stay correctly ordered — consistent with the known **slow, noisy,
non-monotone** convergence of first-moment Nagao sums — but a clean trend-to-integer is **not
demonstrated at X ≤ 10⁵** and is **not asserted**.

**Recommendation (feeds G10 / Paper 136's Nagao paragraph):** rest the claim on **separation** (a
rank-ordered selection signal), not on convergence-to-integer-rank, unless a much-larger-X or
smoothed/averaged computation earns the latter. Reproduce the trend with `python solvable_null.py
--x 100000` (~30 min).

## Withholding — genus-0 conics (β₁ = 0)

The conic carries no loop register; its only output is the HM binary. In the instrument a genus-0
controller has no rank channel at all — `rank_signal` returns `WITHHELD`, never a number. The
withholding is **structural, not stipulated**; a number here would be the falsification, and none
can arise.

| conic | verdict | obstruction places |
|---|---|---|
| x²+y²=1  | soluble   | — |
| x²+y²=3  | insoluble | {2, 3}  (2-adic + 3-adic; ℝ-soluble) |
| x²+y²=5  | soluble   | — |
| x²+y²=21 | insoluble | {3, 7}  (odd-prime pair) |
| x²+y²=6  | insoluble | {2, 3} |
| x²+y²=7  | insoluble | {2, 7} |

All obstruction sets are even-count (Hilbert reciprocity) and match the sum-of-two-squares
criterion. Both local flavours (p = 2 and odd p) are exercised. **[DERIVED]**

## Honest scope (so the green tick is not oversold)

Passing establishes **one** thing: the selection machinery is **not vacuous**. It does **not** earn
the independent identification of the survivor, and it is **not** the divergent-prediction
falsifier — both remain open and honestly conceded in the papers. This answers *"is the framework
empty?"* (no); it does **not** answer *"does the framework predict anything BSD doesn't?"* (still
no). A real, bankable partial win — one of Mr A's four path-to-four items moved from *conceded in
text* to *demonstrated on a bench* — not the whole fourth star.

## Reproduce

```
python solvable_null.py            # X = 10^4, ~25 s → nagao_sums.csv, conic_solubility.csv
python solvable_null.py --x 100000 # the X=10^5 trend (~30 min in Python)
```
