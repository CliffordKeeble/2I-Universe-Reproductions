# Findings — Paper 66 Scale-Ladder Null

**2I Universe Programme · paper66-ladder-null · 10 June 2026**
*Enumeration exhaustive (no sampling). Verdicts per the locked rules of preregistration.md.*

## R1 — Promiscuity: the formula column is DECORATIVE

Coverage of [30, 200] under class C: **171/171 = 100.0%**. Mean spellings per integer:
**≈ 2,000** (median 1,779; min 621; max 5,713). The ladder integers specifically:

| integer | formal spellings | value-distinct |
|---|---|---|
| 44 | 2,762 | 2,466 |
| 84 | 4,601 | 3,779 |
| 137 | 828 | 749 |
| 140 | 2,824 | 2,303 |
| 142 | 1,668 | 1,497 |
| 168 | 4,130 | 3,320 |
| 184 | 2,104 | 1,770 |

Every formula Paper 66 uses is in class (the 184 product form expands to four in-class
product terms; 184 has thousands of in-class spellings regardless). **Rule R1 fires:**
expressibility in icosahedral invariants carries no evidential weight at this scale.
In Paper 66 v2.0 the formula column is re-badged **DECORATIVE** pending a forced
derivation (the b₀ = 11 − 2n_f/3 route).

## R2 — Near-integer content: one survivor, weak, convention-anchored

| exponent | value | dist. to integer | p |
|---|---|---|---|
| ln(M_P/m_p), standard M_P | 44.01241(1) | 0.0124 | **0.025** |
| ln(M_P/m_n), standard M_P | 44.01103(1) | 0.0110 | 0.022 |
| ln(M̄_P/m_p), reduced M_P | 42.40033(1) | 0.4003 | 0.80 |
| ln(R_H/ℓ_P), H₀ = 67.4–73 | 140.21–140.29 | 0.21–0.29 | 0.43–0.59 |
| ln(R/ℓ_P), particle horizon | 141.46 | 0.46 | 0.92 |
| ln(N_b), Hubble volume | 180.48 | 0.48 | 0.97 |
| ln(N_b), particle horizon | 184.09 | 0.093 | 0.19 |

**Rule R2 verdicts:** the Planck–proton exponent lands in the SUGGESTIVE-weak band
(0.001 < p ≤ 0.05) — and only under the standard (non-reduced) Planck convention; under
the reduced convention the near-integerness vanishes entirely. Folding the declared
convention freedom into the look-elsewhere puts the effective p near 0.05, the edge of
the band; the prior peek (disclosed in the pre-registration) means SUGGESTIVE-weak is a
**ceiling**, not a floor. Every cosmic exponent returns **no claim**: 140 is a quarter
off its integer however H₀ is chosen; the baryon exponent moves by 3.6 with the horizon
definition and never enters the band.

## R3 — The system constraint is the assumption restated

Analytically, with the Zoll closure R = 2GM/(c²√3) and M = N·m_p:

> ln(R/ℓ_P) = ln N + ln(m_p/M_P) + ln(2/√3), i.e. **E₁₄₀ = E₁₈₄ − E₄₄ + 0.1438.**

The identity 184 = 140 + 44 is therefore not independent evidence; it is the closure
assumption rewritten in logarithms, exact to a residual 0.144 that no current input
(R known to ±1 in the exponent, N to ±2) can resolve. **Rule R3 verdict:** the
"system test" cannot bite at any achievable resolution; it is consistency by
construction.

## Net verdict for Paper 66 v2.0

The ladder's entire evidential content reduces to **one weak fact**: ln(M_P/m_p) sits
0.012 from an integer, p ≈ 0.025 under one declared convention, ≈ 0.05 with the
convention freedom priced in, void under the alternative convention. The icosahedral
spellings add nothing (R1). The cosmic rows add nothing (R2). The Zoll identity adds
nothing (R3). The 44-row of v2.0 should read, in substance:

> The measured exponent is 44.0124(1). Its proximity to an integer is a 1-in-40
> curiosity under the standard Planck convention and absent under the reduced one.
> The spelling V + E + 2 is one of ~2,500 in-class spellings and is DECORATIVE.
> The incumbent mechanism for this exponent is dimensional transmutation,
> ln(M_P/m_p) ≈ 2π/(b₀α_s) + O(1); the programme's only route to a claim here is to
> derive b₀'s 11 from the structure, at which point the integer part inherits a
> mechanism and the fractional part becomes threshold physics.

This is the same posture Paper 101 earned with its null: the paper carries its own
adversarial computation and is stronger for it.

## Mr Code check-in brief

Directory `paper66-ladder-null/` alongside `paper101-massratio-null/` — three files:
`preregistration.md` (locked first), `enumerate.py` (exhaustive, deterministic,
runtime ≈ minutes, stdlib only), `findings.md` (this file). Hard gates: the
pre-registration must be committed with content unaltered from the run; the script
output must reproduce the tables above verbatim. No other repo writes.

🐕⚛️
