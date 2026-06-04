# PRE-REGISTRATION — Neutrino-Menu Null (Paper 121 §3.1)

**Locked before running. Seed 42. Committed as its own step before `findings.md` exists.**

Brief: CinC, "Neutrino-Menu Null (Paper 121 §3.1)", 4 June 2026. The null Mr Adversary
requested (cycle-1 CIRCULAR/NULL finding on the neutrino section). Direct analogue of
`paper101-massratio-null/`.

## What is under test

Paper 121 §3.1 claims the heaviest neutrino mass is

    m_nu3 = m_e * alpha^(v/f) / (F/d) = m_e * alpha^3 / 4

with exponent 3 read as v/f and divisor 4 as F/d, and asserts "no free parameters."
Mr Adversary's charge: choosing the pair (p,q)=(3,4) that lands on the observed mass is a
**choice from a discrete menu**, and a choice from a menu is a fitted parameter even when
every item is an integer. This null asks: **how many candidate forms `m_e alpha^p / q`,
with p,q from a declared small-integer menu, land within sub-percent of the observed
neutrino mass?** If several, the match is combinatorial fishing; if (3,4) is unique or
near-unique, the form is near-uniquely selected and the claim is defensible (recast).

## Candidate form (LOCKED)

    m_pred(p, q) = m_e * alpha^p / q

- `alpha` = 0.0072973525693  (CODATA 2018 fine-structure constant, 1/137.035999084) — the
  same source used elsewhere in the programme.
- `m_e`  = 0.51099895000 MeV = 510998.95000 eV  (CODATA electron mass; same value as the
  101 mass-ratio null used).
- Claimed canonical pair: `(p, q) = (3, 4)`.

## Integer menu (LOCKED — Menu A, integer-only)

Per CinC's recommendation: **Menu A**, the most hostile-reader-proof choice. It gives the
canonical pair no privilege — the icosahedral readings v/f=3, F/d=4, d=5 all live inside
the menu on equal footing with every other small integer.

- exponent `p in {1, 2, 3, 4, 5}`           (5 values)
- divisor  `q in {1, 2, 3, 4, 5, 6}`         (6 values)

Total candidate count: **5 x 6 = 30** (the sanity figure, analogous to the 101 null's 6750).

Menu A chosen over Menu B (icosahedral-ratio-only) because the selection of *which* ratios
count in B is itself a choice Mr A could attack. Menu A is declared on the flat principle
"all small integers in range, none privileged" and cannot be accused of post-hoc favouring.
**No rational/non-integer entries** (e.g. V/d = 12/5) are admitted — integer-only is the
declared rule.

## Target (LOCKED)

Observed heaviest neutrino mass, normal ordering, lightest mass ~ 0:

    m_nu3 ~= sqrt(Delta m^2_31) = sqrt(2.45e-3 eV^2) = 0.049497 eV = 49.50 meV

- Source: Delta m^2_31 = 2.45e-3 eV^2 (Paper 121 §3.1 / global-fit central value, normal
  ordering). Primary target = **49.50 meV** (the sqrt of the central Delta m^2_31).
- Programme also cites "~50 meV"; recorded as a secondary check (the prediction's relative
  error is reported against both 49.50 and 50.00 meV).
- Experimental/ordering caveat: the lightest neutrino mass is not exactly zero and ordering
  is not fully settled; we test against the central value and report window sensitivity.

## Decision window and rule (LOCKED)

- **Primary window: 1%** relative (Mr A's word was "sub-percent"; 1% is the honest reading
  of §3.1's own quoted match quality, 0.4%–0.9%).
- **Sensitivity:** report **0.5%** and **0.1%**, mirroring the 101 null's three-window report.
- **Decision rule:** the canonical (3,4) is **DISTINGUISHED** if it is the *unique* survivor
  in the primary 1% window, or near-unique (at most one other, on a stated tie-break).
  **More than ~2 survivors in the primary window => NOT DISTINGUISHED**: the menu selection
  is a hidden fit and §3.1's "no free parameters" is false and must be reworded.

## Generation-scaling sub-test (LOCKED — §2.3)

§3.1 also claims `r = 2^(d/2) = 2^2.5` for generation scaling, giving
`Delta m^2_31 / Delta m^2_21 = 2^d = 32`. Same menu logic applied to the observed ratio:

    R_obs = Delta m^2_31 / Delta m^2_21 = 2.45e-3 / 7.53e-5 = 32.54
    (source: Paper 121 §3.1 / global-fit central values; sensitivity to Delta m^2_21 noted)

Candidate form families (LOCKED), d = 5 fixed (icosahedral vertex degree):

- family F1: `2^s`        for `s in {1,2,3,4,5,6}`   (claimed: s = d = 5 -> 32)
- family F2: `k^(d/2)`    for `k in {1,2,3,4,5,6}`   (claimed: k = 2 -> 2^2.5 = r; ratio k^d)

Note `2^5 = 4^2.5 = 32` are degenerate. Total: 12 candidates (10 distinct values).
Windows 1%/0.5%/0.1% as above; additionally report 2% and 5% to locate the nearest form,
since the claim's own match quality may exceed 1%. Same DISTINGUISHED/NOT rule.

## Anti-tuning commitment

The menu and rule above are frozen. Survivor counts will not be changed by adjusting the
menu after seeing results. Any defensible menu revision would be a NEW pre-registration and
a NEW run, documented as such. This is the one null where post-hoc menu-tuning would be
self-falsifying — guarded hardest.
