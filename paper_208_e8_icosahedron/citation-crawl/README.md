# Paper 208 §8 item 2 — Inspire-HEP / ADS cross-citation crawl

Converts Paper 208's four negative-existential cross-citation claims (§§3.4, 4.4,
4.6, 5.2) and the §4.8 Dechant-forward claim from "consistent across a broad
sample" to "exhaustive (to the depth of the citation corpus)". A database job, not
a reading job: this directory produces the citation network and the cross-reference
tests. Reading bodies and adjudicating seam status are **not** done here (Scout/Fizz).

**Database coverage: Inspire-exhaustive, ADS-pending.** ADS walled at HTTP 401
(needs a Bearer token); the ADS leg was not run. See `crosscite_findings.md` →
"Database coverage" before citing anything as exhaustive.

## Headline

- §3.4 P↔A — **EMPTY** (holds) · §4.4 S↔F — **EMPTY** (holds)
- §4.6 E↔U — **1 hit** (recid 2876584, Truini–Marrani–Rios–de Graaf) → wording call
- §4.8 Dechant-forward — **holds** (no physical-parameter attacher)
- §5.2 — deferred (not citation-testable from metadata)

## Files

- `inspire_crawl.py` — the crawl: resolves anchors, pulls Dechant forward-citers,
  community-tags, and runs the four claim intersection tests. Inspire-HEP only.
- `inspire_ads_citers_dechant.csv` — deduped, community-tagged citer spine
  (Inspire leg; ADS union pending).
- `crosscite_findings.md` — the four-claim table, filled, with purpose notes. No
  seam adjudication (hard gate 1).
- `claim_findings.json` — machine-readable claim outputs.
- `raw/` — cached Inspire API responses (anchors, citer sets, side-anchor citers).
