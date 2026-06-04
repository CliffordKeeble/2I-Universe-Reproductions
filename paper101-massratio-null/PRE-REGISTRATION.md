# PRE-REGISTRATION — Paper 101 §5.3 Mass-Ratio Null

**Status: faithful transcription of an already-published, pre-registered design.**

This file is **not** a new pre-registration. The design below is locked by what
Paper 101 v3.0 (concept DOI 10.5281/zenodo.18817872; v3 DOI 10.5281/zenodo.20546831),
§5.3 and reference [10], *already published*. It is transcribed here verbatim so the
published null can be re-run reproducibly. Per the reconstruction brief (CinC, 4 June
2026): the search space, target set, integer set, operations and decision window are
**not to be adjusted to improve the match** — doing so would defeat the pre-registration.

See the honesty note in `findings.md`: this artifact is a *reconstruction* of the
published pre-registered design, committed after the fact because the original null was
run in-session and never committed (Pattern 74). The design itself was pre-registered in
the published paper; the code here is its faithful re-implementation.

## Locked design (transcribed from §5.3)

- **Seed:** 42 (convention; the enumeration is exhaustive and deterministic, so no random
  draw is consumed — seeded anyway to match the declared protocol).

- **Search space:** every surviving 2I-filtered mode up to `l = 120`. Surviving levels are
  computed independently from binary-icosahedral (2I) character theory in `spectrum.py`
  (Molien series of 2I acting on C^2), and must reproduce §4.3's spectrum — in particular
  the first mode at `l = 12`, eigenvalue `λ = l(l+2) = 168`, and `l = 42`, `λ = 1848`.

- **Corrections:** ± a single icosahedral integer drawn from the **declared set**
  `{1, 2, 6, 11, 12, 13, 18, 19, 20, 30, 31, 41, 42, 62, 168}` (15 integers).

- **Operations:** `+` and `−` (a single correction, either sign).

- **Targets (five):** `m_p/m_e`, `m_n/m_e`, `m_μ/m_e`, `m_π/m_e`, `m_τ/m_e` — dimensionless
  mass ratios in the integer-eigenvalue (`>100`) regime. Values used (recorded for
  reproducibility):
  - `m_p/m_e   = 1836.152673426`  (CODATA proton–electron mass ratio)
  - `m_n/m_e   = 1838.68366173`   (CODATA neutron–electron mass ratio)
  - `m_μ/m_e   = 206.7682830`     (CODATA muon–electron mass ratio)
  - `m_π/m_e   = 139.57039 / 0.51099895 = 273.132440`  (PDG 2024 charged pion / CODATA m_e)
  - `m_τ/m_e   = 1776.86 / 0.51099895 = 3477.228280`   (PDG 2024 τ / CODATA m_e)

- **Excluded before running (declared, not silently dropped):** `m_n/m_p ≈ 1.0014` —
  declared in the pre-registration but excluded as unreachable by integer eigenvalues
  (it lives near 1, not in the `>100` regime). Documented as declared-then-excluded.

- **Decision window:** `0.01%` relative. Sensitivity also reported at `0.05%` and `0.005%`.

- **Decision rule (locked):** **more than one survivor for a given target ⇒ NOT
  distinguished.** A clean unique identification within the window is the criterion for the
  match being distinguished from what a dense-enough correction set would produce by chance.

## Combination-count sanity figure

§5.3 reports "roughly 7,500" combinations:
`(surviving modes ≤ l=120) × 15 integers × 2 operations × 5 targets`.
This is a checkable sanity figure, not part of the locked decision. If the reconstructed
count is *wildly* off, the spectrum computation is suspect and the run must stop. (See
`findings.md` for the reconstructed value and the small, benign discrepancy.)
