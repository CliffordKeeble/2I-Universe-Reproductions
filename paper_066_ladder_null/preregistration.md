# Pre-registration — Paper 66 Scale-Ladder Null

**2I Universe Programme · paper66-ladder-null**
*Drafted by Tor ⚛️, 10 June 2026, locked before any enumeration was run.*
*Pattern 75 discipline. Companion to `paper101-massratio-null/`.*

## Purpose

Paper 66 (The Bootstrap Spacetime Lattice, v1.1) presents a scale ladder of six exponents
{44, 84, 137, 140, 142, 184}, each written as a small polynomial in the icosahedral
invariants V = 12, E = 30, F = 20, D = 3, χ = 2, and identifies 44 = V + E + 2 with the
measured ln(M_Planck/m_p). This null quantifies, before Paper 66 v2.0 is drafted:

1. **Promiscuity of the formula class** — what fraction of integers in the ladder's range
   is expressible at all, and with what multiplicity of spellings;
2. **The near-integer content of the measured exponents** — how unusual it is that the
   measured values sit close to integers, under declared conventions;
3. **The system constraint** — whether the Zoll identity 184 = 140 + 44 is evidence or
   assumption.

## Declared conventions (locked)

- **Planck mass:** the standard (non-reduced) M_P = √(ħc/G), CODATA-2022
  G = 6.67430(15)×10⁻¹¹ m³ kg⁻¹ s⁻². Sensitivity to the reduced convention
  M̄_P = M_P/√(8π) is **reported**, not selected.
- **Hadron anchor:** the proton, m_p c² = 938.27208816 MeV. Sensitivity to the neutron
  anchor is reported.
- **Cosmic radius:** reported across the definitional spread (Hubble radius c/H₀ for
  H₀ ∈ [67.4, 73.0] km/s/Mpc; particle horizon ≈ 4.4×10²⁶ m), not at a single value.
- **Baryon count:** n_b = 0.25 m⁻³ over the same two volumes; reported as a spread.

## Declared disclosure

The value ln(M_P/m_p) = 44.0125(1) was already known to the bench before this
pre-registration (it motivated the test). It cannot be un-peeked. This pre-registration
therefore locks the **class definition, the conventions, and the decision thresholds**;
the near-integer p-value for the 44-row is reported as a quantity whose target was
known in advance, with the convention freedom declared as additional look-elsewhere.

## Declared formula class C (locked)

An expression is a sum of **at most 3 terms**, each term k·m with integer coefficient
1 ≤ |k| ≤ 7 and monomial m drawn from:

- the unit **1**;
- the five invariants **V, E, F, D, χ**;
- all pairwise products including squares: **V², E², F², D², χ², VE, VF, VD, Vχ, EF,
  ED, Eχ, FD, Fχ, Dχ**.

(21 monomials total.) This class is chosen to be the smallest that contains every
formula Paper 66 actually uses: V+E+2 (writing 2 as χ or as 2·1), 2(V+E), 7F, 7F−D,
7F+2, and (F−V)(F+D) = F² + FD − VF − VD. Spellings are counted as distinct formal
expressions over the monomial list (so χ and 2·1 count separately; numeric collisions
such as FD = Eχ = 60 count separately). A secondary count collapses to distinct
(coefficient, value) multisets.

## Locked decision rules

- **R1 (promiscuity).** If ≥ 95% of integers in [30, 200] are expressible in C, the
  formula column of the ladder is re-badged **DECORATIVE** in Paper 66 v2.0 —
  expressibility carries no evidential weight — pending a forced derivation
  (the b₀ route named in the v2.0 build sheet).
- **R2 (near-integer content).** For each independently measured exponent, report
  p = 2·min(frac, 1−frac), the two-sided probability that a random real lands as close
  to an integer. Banking thresholds: p ≤ 0.001 → OBSERVED-strong; 0.001 < p ≤ 0.05 →
  SUGGESTIVE-weak; p > 0.05 → no claim. Convention sensitivity is reported beside
  every p.
- **R3 (system constraint).** Derive analytically whether 184 = 140 + 44 is independent
  evidence or the Zoll closure assumption restated; state the exact residual term and
  the resolution at which the integer identity could ever bite.

## Method

Exhaustive enumeration (no sampling, no seed needed): all 1-, 2-, and 3-term
expressions over the 21 monomials with coefficients in ±{1..7}; coverage and spelling
multiplicity tabulated for every integer in [30, 200]; the six ladder formulas verified;
p-values computed from CODATA inputs. Script: `enumerate.py`. Findings: `findings.md`.
Whatever the rules return is reported as-is.

🐕⚛️
