# Findings — Neutrino-Menu Null (Paper 121 §3.1)

**Date:** 4 June 2026
**Brief:** CinC, "Neutrino-Menu Null (Paper 121 §3.1)" — the null Mr Adversary requested
(cycle-1 CIRCULAR/NULL finding on the neutrino section).
**Pre-registration:** locked and committed *before* this run at `e8f758d` (menu, target,
windows, decision rule frozen first; see `PRE-REGISTRATION.md`).

## Verdict (one line)

- **Mass form `(3,4)`: DISTINGUISHED** — the *unique* survivor in the primary 1% window;
  the next candidate is 19.8% away. The "combinatorial fishing" charge is **refuted** for
  the mass form. "Zero free parameters" should be recast as "uniquely selected from the
  integer menu" (defensible). Status stays **OBSERVED** (a clean null clears promiscuity; it
  does not derive the mechanism — same as the 101 null).
- **Generation scaling `2^d = 32`: near-unique but a ~1.7% match, not sub-percent.** Not a
  fishing problem (only `32` lies within 2% of the observed ratio), but it *fails the 1%
  window* against the central observed ratio. The caveat here is **match quality**, not menu
  promiscuity — §3.1 should disclose the scaling match is ~1.7%, not sub-percent.

## Test 1 — mass form `m_e · α^p / q` (Menu A)

- Menu A (integer-only, no privilege to the canonical pair): `p ∈ {1..5}`, `q ∈ {1..6}`
  → **30 candidates** (the sanity figure).
- Constants: α = 0.0072973525693 (CODATA), m_e = 0.51099895 MeV.
- Target: m_ν3 = √(Δm²₃₁) = √(2.45e-3) = **49.50 meV** (normal ordering, lightest ≈ 0).
  Secondary cited value 50 meV also reported.
- Canonical prediction `(3,4)` = **49.6428 meV** → **0.294%** vs 49.50 meV (0.714% vs 50 meV).

| window | survivors | detail |
|--------|-----------|--------|
| 1.0%   | **1**     | `(3,4)` 49.643 meV, rel 0.294% — CANONICAL, unique |
| 0.5%   | **1**     | `(3,4)` only |
| 0.1%   | 0         | `(3,4)` at 0.294% falls outside |

Nearest non-canonical candidates: `(3,5)` = 39.71 meV (**19.8%**), `(3,6)` = 33.10 meV
(33.1%), `(3,3)` = 66.19 meV (33.7%). Every other (p,q) is orders of magnitude off, because
each unit of exponent p multiplies by α ≈ 1/137 — so only `p=3` places `m_e·α^p` in the
~10²-meV band at all, and within that band only `q=4` reaches the target window.

**Why this is a strong pass:** unlike the 101 null's *additive* ±integer corrections (a
dense set), the *multiplicative* α^p here spans ~9 orders of magnitude across the menu, so
the candidate set is intrinsically sparse near any fixed target. `(3,4)` is not one of
several near-misses — it is alone in the window by a factor of ~70× in relative error over
the runner-up. Menu A, which gives `(3,4)` no privilege, still selects it uniquely.

## Test 2 — generation scaling (R = Δm²₃₁/Δm²₂₁)

- R_obs = 2.45e-3 / 7.53e-5 = **32.54** (central values; see sensitivity below).
- Claim: `2^d = 2^5 = 32`. Families: `2^s` (s∈{1..6}), `k^(d/2)=k^2.5` (k∈{1..6}), d=5.

| window | survivors | detail |
|--------|-----------|--------|
| 1.0%   | 0         | nothing within 1% of 32.54 |
| 0.5%   | 0         | — |
| 0.1%   | 0         | — |
| 2.0%   | 2 (degenerate value 32) | `2^5 = 32` and `4^2.5 = 32`, both rel **1.65%** |
| 5.0%   | 2 (same)  | still only the value 32 |

The only small form anywhere near R_obs is **32** (reachable two ways, `2^5` and `4^2.5` —
the same number). So the scaling form is *near-uniquely selected* (no competing forms → not
fishing), **but** `32` sits **1.65%** below the observed `32.54` — it does not reach the
sub-percent window the mass form clears.

**Sensitivity (Δm²₂₁ ambiguity):** with Δm²₂₁ = 7.42e-5 (newer global fits) R_obs = 33.0 and
the gap to 32 widens to ~3.1%; with 7.53e-5 it is 1.65%. Either way the scaling claim is a
**~2–3% match**, not sub-percent. This is honest and should be stated: the `2^d` reading is
the only natural small form (good against the fishing charge) but is a coarser match than
§3.1's mass-form match.

## Bottom line for Paper 121 §3.1

- **Mass form:** keep the match; recast "no free parameters" → "the icosahedral form
  `(v/f, F/d) = (3,4)` is **uniquely** selected from the integer menu (null:
  `neutrino-mass-null/`)." Status OBSERVED. Mr A gets his distribution; the promiscuity
  charge does not hold.
- **Generation scaling:** disclose the match quality. `2^d = 32` is the only natural small
  form near the observed ratio (not fishing), **but the match is ~1.7% (2–3% under
  Δm²₂₁ uncertainty), not sub-percent** — word §3.1 so it does not imply the scaling ratio
  matches at the same precision as the mass.

Outputs: `matches.csv` (survivors per test/window with relative errors).
