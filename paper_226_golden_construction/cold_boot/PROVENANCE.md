# cold_boot/ — provenance

**Paper 226 · Golden Construction · Part A (cold-script check-in)**

## What is here (committed, verbatim)

- `verify_c_from_sigma.py` — the cold-boot verification script, byte-identical to the received artifact.
- `verification_run.log` — its run log, byte-identical.

Both are the machine-checkable outputs of a cold-boot construction and are committed unmodified
as the candidate reproduction the draft's Appendix B asks to be re-derived before merge.

## What is NOT here (by standing rule)

The two cold-boot **paper `.md` drafts** remain in the user's Downloads, **uncommitted**, per the
programme's standing rule that unreleased paper drafts stay out of the repo until Zenodo (see the
provenance note atop `../golden_algebra.py`). Their claim statements are instead *vendored with
source citations* into `../verify_c_from_sigma_repro.py`, exactly as Run #1 vendored its
definitions.

## Model provenance (two independent cold boots)

- **Fable 5 max** produced the scripted construction: the draft
  `paper_2XX_fourth_face_charge_conjugation_v0.1.md` (formal Lemma 1 / Theorem 1 / M1–M3 / Appendix
  A+B) **and** the two artifacts committed here.
- **Opus 4.8 max** produced a prose construction, `cold_boot_charge_conjugation_C_equals_sigma.md`
  — a **second independent cold witness** to C = σ, reached by a different route.

The two witnesses **converge** on the headline C = σ but **diverge** on the gauge structure: Fable
gauges the non-compact split torus SO(1,1)⁺ over K (charge quantised by the unit lattice); Opus
gauges a compact U(1) over K(i) (i required for the gauging). This fork is not adjudicated in this
reproduction — it is benched in `../W4_derivation.md` / `../w4_keystone.py` (the "one torus, two
real forms" keystone) and in Part E's eight-cell 𝒞 frontier.

## Assertion-count resolution (A-3)

| Source | Count | What it is |
|---|---|---|
| `verify_c_from_sigma.py` | **49** | `check()` calls |
| `verification_run.log` | **49** | `[PASS]` lines (51 physical lines: 49 + blank + banner) |
| Fable draft (abstract, App A, App B) | 50 | claimed "50 exact-symbolic checks" ×3 |
| relay | 52 | line-oriented miscount |

Script and log **agree exactly at 49** — no assertion present in one and absent in the other. The
draft's "50" is a documentation off-by-one (repeated three times, a write-time slip); the relay's
"52" counts lines, not checks. **Resolves to bookkeeping.** True verified count = 49, all PASS,
sympy 1.14.0.

## Independence caveat (A-1)

`verify_c_from_sigma_repro.py` was written **after** this cold script was read and run (the cold
package arrived piecemeal during the check-in). Independence is therefore preserved at the
**derivation level** — the reproduction uses a different representation (the regular representation
of 𝕂 on basis (1, δ), where N = det and σ = Z-conjugation), and each assertion is re-derived from
the draft's Appendix-A claim statements with citations — **not** at the script-blind level. Stated
plainly rather than claiming a blindness we no longer have.

🐕☕⬡
