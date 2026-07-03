# Paper 39 — Redshift-mechanism Pattern-75 null (job iii)

**Mr Code, 13 June 2026.** Brief: Tor ⚛️ via CinC via Cliff.
Script: [`redshift_null.py`](redshift_null.py) · outputs: `z_by_cell.csv`,
`null_grid_summary.csv`, `backsolve_by_cell.csv`. All exponential arithmetic in
mpmath at 50 digits. Inputs pinned at full precision (CODATA 2018 α, hbar_c, m_p;
R = 13.8 Gly held as a fixed input per brief, itself OBSERVED/contested).

Mechanism under test: `z = exp(N·ε) − 1`, `ε = α^k`, `N = R / l_cell`.
Paper 39 v1.7: k = 19 (= 5+11+3), l = 0.5 fm, claims z ≈ 1096.
Observed recombination redshift used for the gate: **z_obs = 1091**.

---

## Test 1 — z at the three corpus cell scales (k = 19 fixed) — *DERIVED*

| cell | l (fm) | N·ε (exponent) | **z** |
|---|---|---|---|
| 0.5 fm (v1.7 as-written) | 0.5 | 6.561075174 | **706.0314679** |
| u = ƛ_p/4 (Seam A re-anchor) | 0.0525772276 | 62.39464762 | **1.252135019 × 10²⁷** |
| ƛ_p (Compton rung, 108's ruler) | 0.2103089103 | 15.59866191 | **5 948 571.93** |

Backsolve check: k = 19 hits z = 1091 **exactly** at l = **0.468931853 fm** — the
paper's 0.5 fm sits **6.6 % away** from that cell.

**Reading.** Confirms the whiteboard on the headline point and sharpens it on one:
- The paper's own honest cell (0.5 fm) gives **z = 706, not 1096/1091**. The "N·ε ≈ 7"
  headline is a rounding *in the exponent* (true value 6.561), and z is exp of it — so at
  the stated cell the mechanism is already ~35 % low in z.
- The Seam-A re-anchored ruler u overshoots to **z ≈ 1.25 × 10²⁷** — not the
  "10⁶-ish" the whiteboard guessed, but **twenty-one orders of magnitude worse**. The
  exponent is unforgiving (the brief's own warning): u is ~9.5× smaller than 0.5 fm, so
  the exponent is ~9.5× larger, and z is exponential in that.
- The Compton rung gives z ≈ 5.9 × 10⁶.

**Outcome: (B).** No corpus-honest cell lands near 1091. The mechanism does not close at
any of the three real corpus length scales. **Verdict: does not survive.**

---

## Test 2 — the null: is 1091 special, or just reachable? — *DERIVED / OBSERVED*

Elegant exponent family (corpus primes 3–23 plus small-prime sums, incl. the paper's
19 = 5+11+3): `{3,5,7,8,11,13,14,16,17,18,19,20,22,23}`. Cell grid: log-uniform,
4000 points across [0.01, 1.0] fm. ±5 % of z = 1091 ⇔ exponent band E ∈ [6.9445, 7.0445].

**2a — which elegant exponents can reach 1091 inside a plausible cell window?**
For each k there is exactly one cell `l*(k) = R·α^k / ln(1+1091)` that hits z = 1091.
Consecutive elegant exponents shift l* by a factor of **1/α ≈ 137**, so the window
[0.01, 1.0] fm admits **exactly one** exponent:

| k | l*(k) (fm) | in [0.01,1.0]? |
|---|---|---|
| 18 | 64.26 | no |
| **19** | **0.46893185** | **YES** |
| 20 | 0.003422 | no |

**2b — grid density.** **12 / 56000 = 0.021 %** of elegant (k, l) pairs land within
±5 % of z = 1091. The neighbourhood is **sparse, not crowded** — hitting 1091 is *not*
reaching into a dense net.

**2c — isolation of the (19, ~0.5 fm) point.** The ±5 %-z band at k = 19 is the cell
interval [0.4657, 0.4724] fm — width **1.43 % of l\***. So the point is **sharp in l**,
not a plateau: the cell must be tuned to ~1 % to hold z near 1091. The nearest other
elegant exponent's matching cell is k = 18 at 64.3 fm — a factor 137 away.

**2d — backsolve transparency (the decisive panel).**

| cell | l (fm) | exponent k required for z = 1091 | implied Λ_QCD = ħc/l (MeV) |
|---|---|---|---|
| 0.5 fm (v1.7) | 0.5 | **18.987** | 394.65 |
| u = ƛ_p/4 | 0.0526 | 19.445 | 3753.1 |
| ƛ_p | 0.2103 | 19.163 | 938.27 |
| l\*(19) backsolved | 0.46893 | 19.000 | **420.80** |

**Reading — why isolation here is not content.** The sparse, sharp, single-exponent
picture (2a–2c) looks like a PASS on the brief's first two clauses. It is not, because
the brief's PASS requires the isolated pair's cell to *match a corpus ruler derived
elsewhere, not backsolved*. It fails that clause three ways, all visible in 2d:

1. **The matching cell is not a corpus ruler.** k = 19 hits 1091 only at l = 0.4689 fm.
   The corpus rulers are u = 0.0526 fm and ƛ_p = 0.2103 fm. 0.4689 fm is neither — it is
   backsolved to land on 1091.
2. **Holding the cell at the paper's 0.5 fm forces a non-integer exponent, k = 18.987** —
   the "elegant integer 19" story breaks; 0.5 fm + integer 19 gives z = 706 (Test 1).
3. **The "Λ_QCD ≈ 420 MeV prediction" is this backsolve read forwards.** The cell that
   makes k = 19 hit 1091 (0.4689 fm) implies Λ_QCD = **420.8 MeV** — exactly v1.7's
   quoted figure. It is not a prediction; it is the redshift target re-expressed as an
   energy scale. The circularity is now on the page.

**Verdict: FAIL — the redshift number is selection, not separation from null.** z = 1091
is reachable only by a backsolved cell (0.4689 fm) with no independent derivation; the
v1.7 Λ_QCD "prediction" is that backsolve. §3 ships as **OPEN**: mechanism-class +
epistemics only, redshift number demoted to a stated open problem.

---

## Test 3 — Seam C residual: the 27/5 null for Ω_c/Ω_b — *OBSERVED*

Planck 2018 base-ΛCDM (TT,TE,EE+lowE+lensing): Ω_b h² = 0.02237 ± 0.00015,
Ω_c h² = 0.1200 ± 0.0012 (full-table values; abstract rounds to 0.0224 / 0.120).
Propagated: **Ω_c/Ω_b = 5.3643 ± 0.0646.** 1σ band [5.2997, 5.4289]; 2σ [5.2352, 5.4935].

| simple form | value | \|Δ to 5.3643\| | within 1σ? | within 2σ? |
|---|---|---|---|---|
| **16/3** | 5.33333 | **0.03099** | **YES** | YES |
| 27/5 | 5.40000 | 0.03567 | **YES** | YES |
| 2·φ² | 5.23607 | 0.12826 | no | YES |
| 21/4 | 5.25000 | 0.11433 | no | YES |
| 3√3 | 5.19615 | 0.16817 | no | no |

**Reading.** **Two** simple closed forms sit inside the 1σ bar — and **16/3 is closer to
the data than 27/5** (0.031 vs 0.036). At 2σ, four forms qualify. 27/5 is therefore *not*
the uniquely closest simple form, and the neighbourhood is not sparse.

**Verdict: FAIL (flag).** The numerical match "27/5 ≈ Ω_c/Ω_b" does not separate 27/5
from rival simple forms — 16/3 fits better. Paper 59 v2.0's Kuramoto mode-counting
derivation must carry the 27/5 claim **on its own**; the number alone does not earn it.
Hand the count up to CinC against 59's standing.

---

## Verdicts (one line each)

- **Test 1:** redshift mechanism **does not survive** re-anchoring — outcome (B). No
  corpus-honest cell lands near 1091 (0.5 fm → 706; u → 1.3×10²⁷; ƛ_p → 5.9×10⁶).
- **Test 2:** **FAIL** — 1091 is selection, not separation. Reachable only via a
  backsolved 0.4689 fm cell that matches no independent ruler; the Λ_QCD ≈ 420 MeV
  "prediction" is that backsolve. §3 → OPEN.
- **Test 3:** **FAIL (flag)** — 16/3 fits Ω_c/Ω_b better than 27/5, both inside 1σ. The
  mode-counting derivation must carry 27/5 alone.

Tor's warm picture — (B) and FAIL — holds under cold compute, and is sharpened on two
points the whiteboard could not see: the u-cell overshoot is 10²⁷ not 10⁶, and 16/3
beats 27/5 on the data. The programme tells the true thing.

🐕☕⬡
