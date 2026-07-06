# Paper 174 — The Icosahedral Spectrum: Appendix A reproduction

**Status:** all acceptance targets reproduced (15/15 PASS). See
`icosahedral_spectrum_results.md` for the raw run log.

## What this is

Paper 174 predates the Reproductions repo; its Appendix A computation existed
only as "source available from the author" and was never merged. This folder
is a **first-add**: an independent, auditable re-derivation of the paper's
numerical claims, committed as an acceptance test. It is *not* a paper-text
edit and does not touch any adjacent finding (Weeks-2006 citation scope,
ring-series-vs-harmonic-series phrasing) — pure reproduction + repo hygiene.

Language note: the brief requested C#. The repo is Python-first (CLAUDE.md;
71 existing Python files, 0 C#). Cliff confirmed Python for consistency with
the existing reproductions, so this is Python.

## The quantity

    m_k = dim( Sym^k(C^2) )^{2I}

the multiplicity of the trivial representation of the binary icosahedral
group 2I (order 120) in the degree-k SU(2) harmonics, for k = 0..120.
Geometrically: the dimension of the space of degree-k harmonics on S^3 that
are invariant under 2I — i.e. descend to the Poincaré homology sphere
S^3 / 2I. A survivor at degree k sits at Laplacian eigenvalue λ = k(k+2).

## Four independent methods (all cross-checked at every k)

| # | Method | Route | Tier |
|---|--------|-------|------|
| 1 | Character sum, paper eq. 2: `m_k = (1/120) Σ_classes size·χ_k(α)`, `χ_k(α)=sin((k+1)α)/sin(α)` | mpmath sin/cos, rounded | **OBSERVED** |
| 2 | Reduced monomial basis, paper eq. 6: `#{(a,b,c): 12a+20b+30c=k, a,b≥0, c∈{0,1}}` | exact integer | **DERIVED** |
| 3 | Klein complete-intersection Hilbert series, paper eq. 5: `(1−t⁶⁰)/((1−t¹²)(1−t²⁰)(1−t³⁰))` | exact integer | **DERIVED** |
| 4 | Kostant form: `(1+t³⁰)/((1−t¹²)(1−t²⁰))` | exact integer | **DERIVED** |

The three exact routes (2, 3, 4) agreeing across k=0..120 is a **DERIVED**
statement; the character route (1) reproducing them after rounding (max
deviation from an integer **1.68e-49**, at 50-digit precision) corroborates
it independently.

## Results (all PASS)

- **Four-method agreement** for k = 0..120: **zero mismatches**.
- **λ₁ = 168**: first non-trivial survivor at **k = 12**, λ = 12·14 = 168.
- **All odd k killed**: `m_k = 0` for every odd k (all generators 12,20,30
  are even, so odd degrees are unreachable). — STRUCTURAL, confirmed.
- **Exactly 15 killed even modes**:
  `{2,4,6,8,10,14,16,18,22,26,28,34,38,46,58}`, largest = 58. These are
  `2 ×` the 15 gaps of the numerical semigroup ⟨6,10,15⟩.
- **Frobenius(⟨6,10,15⟩) = 29** (computed independently by dynamic
  programming), so 2·29 = **58** is the largest even non-member; every even
  k > 58 survives (verified through k = 120).
- **Series identity**: methods 3 and 4 produce identical coefficients for all
  k ≤ 120 — the concrete form of `(1−t⁶⁰)/(1−t³⁰) = 1 + t³⁰`, the
  Kostant/Klein equivalence. — DERIVED.
- **Spot values**: m₁₂ = 1, m₂₀ = 1, m₂₄ = 1, m₆₀ = 2. All match.

46 survivors in 0..120; the first with multiplicity ≥ 2 is m₆₀ = 2, and
m₁₂₀ = 3.

## How to reproduce

```
python icosahedral_spectrum.py
```

Standalone, no arguments. Prints the PASS/FAIL table and the survivor
spectrum, writes `icosahedral_spectrum_results.md`, and exits non-zero on any
acceptance failure (a reproduction that disagrees with the paper must not
pass silently).

## Environment

- Python 3.13.14, mpmath 1.3.0 (character sum at 50 decimal digits).
- Methods 2–4 use only Python integer arithmetic — exact, no rounding.
