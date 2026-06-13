# Paper 208 §8 item 2 — Inspire-HEP cross-citation crawl: findings

Issued from the Fizz/CinC brief, 12 Jun 2026 · run by Mr Code, 13 Jun 2026.
Plain findings, no seam adjudication (hard gate 1). Purpose-qualifier calls are
Fizz/CinC's; this file reports only what the citation graph contains.

## Database coverage — READ FIRST (hard gate 2)

**Inspire-exhaustive, ADS-pending.** The NASA ADS API (`api.adsabs.harvard.edu`)
returned **HTTP 401** without a Bearer token; no ADS leg was run. Per the brief's
gate 2 and network note, this deliverable is **not yet "exhaustive" in the
two-database sense** — it is exhaustive to the depth of the Inspire-HEP corpus
only. ADS is the decorrelating second resource (W-107) and remains outstanding.
To close it: supply an ADS API token and re-run the ADS leg, or run this leg from
a surface with ADS access. **Do not cite these as two-database-exhaustive until
ADS is run.**

## Anchors (step zero, reconciled)

| Ref | Paper | Inspire recid | arXiv | DOI / venue |
|---|---|---|---|---|
| [4] | Dechant 2016, *Birth of E8 out of spinors of the icosahedron* | **1422619** | 1602.05985 | 10.1098/rspa.2015.0504 (RSPA 472:**20150504** — same object) |
| [31] | Dechant 2017, *E8 geometry from a Clifford perspective* | **1427769** | 1603.04805 | 10.1007/s00006-016-0675-9 (AACA 27) |
| [22] | Chaudhuri–Polchinski 1995 (CHL) | 396016 | hep-th/9506048 | PRD 52:7168 |
| [36] | Kronheimer–Nakajima 1990 (ALE/ADHM) | 415625 | — (1990) | Math. Ann. 288 |

arXiv:1602.05985 ↔ DOI 10.1098/rspa.2015.0504 reconciled: the "20150504" in the
paper's own reference is the RSPA article number, not a separate eprint.

## The spine — forward citers of the two Dechant anchors

**14 forward citations (7 + 7), 10 distinct records** after dedup across the two
anchors. Full table in `inspire_ads_citers_dechant.csv`. Community spread: 8 of 10
tag **P** (polytope / Clifford-spinor); 2 untagged are a Riabchenko 2025
preprint/proceedings near-duplicate pair (AdS/CFT integrability + ML).

The population is, without exception, **applied-Clifford-algebra / root-system /
ML-of-Coxeter / AdS-CFT-integrability** work — largely Dechant self-citation
(5 of 10) plus Sen, Breuils, Chen, Moxness, Riabchenko. No string-theory,
no field-theory-phenomenology, no octonion-unification citer appears.

## The four claims

Method: each side of a pair is anchored to canonical Inspire recid(s); the full
citing set of each side is pulled and intersected. A paper in the intersection
cites **both** sides. Whether it does so **for the qualified purpose** is a Scout
read, not ours. **Empty intersection = claim holds at exhaustive Inspire depth.**

| Claim | Pair | Sides (recids) | |left| ∩ |right| | Result |
|---|---|---|---|---|
| **§3.4** | P ↔ A | Dechant {1422619,1427769} ∩ ALE/KN {415625} | 10 ∩ 61 → **0** | **EMPTY — holds** |
| **§4.4** | S ↔ F | heterotic-E8×E8 {213013,207031} ∩ twin-Higgs {685922,699848} | 1529 ∩ 908 → **0** | **EMPTY — holds** |
| **§4.6** | E ↔ U | emergent-E8/Coldea {842444} ∩ EJA/octonion {2095962,2667356,1658866,555166} | 294 ∩ 108 → **1** | **1 HIT — wording call** |
| **§5.2** | internal | — | — | **DEFERRED (not citation-testable)** |

### §3.4 P↔A — EMPTY (holds)
No citer of either Dechant anchor also cites Kronheimer–Nakajima. The known
instanton bridges (Allen–Sutcliffe 1302.4664, Choi–Lee Symmetry 10(8):326) do not
appear in the Dechant forward-citation set at all, so they did not even need
purpose-exclusion. The S³/2I orbit-space / Poincaré-homology-sphere bridge is
**not surfaced** in Inspire.

### §4.4 S↔F — EMPTY (holds)
Of 1529 distinct citers of the heterotic-E8×E8 anchors and 908 citers of the
twin-Higgs anchors, **zero papers cite both**. The two-E8-factor-as-particle-pairing
bridge is **not surfaced** in Inspire — for *any* purpose, let alone the qualified
one. The string-phenomenology and twin/mirror-matter literatures are citation-disjoint
at this depth.

### §4.6 E↔U — 1 HIT → flag for wording decision
One paper cites both the emergent-E8 (Coldea CoNb₂O₆) side and the EJA/octonion side:

> **recid 2876584** — Truini, Marrani, Rios, de Graaf, *Exceptional Periodicity and
> Magic Star algebras* (no arXiv; citation_count 0).

This is a genuine candidate bridge: the authors are the Magic-Star / exceptional-
periodicity programme (E8-as-fundamental-algebra, U-side) and the record cites
Coldea's emergent-E8 experiment (E-side). **Per gate 1, no adjudication here** —
whether this constitutes the qualified "substantive interaction" or a passing
"E8 also appears in nature" citation is Fizz/CinC's call. **Flagged for a §4.6
wording decision.**

### §5.2 internal — DEFERRED
The five involutions of §5.2 are internal to Paper 208; deciding whether a single
external source unifies them requires the involution list and a body read
(Scout/Dig), not citation metadata. Out of scope for this crawl, declared.

## §4.8 — Dechant-forward physical-parameter exhaustiveness

**Holds (Inspire depth).** All 10 forward citers are applied-Clifford / root-system /
ML / integrability math. Six carry a *secondary* hep-th cross-list, but **none is
primary hep-ph or hep-th**, and none attaches a physical parameter (mass ratio,
coupling) to the Dechant construction by title/venue. The 6 hep-th-crosslisted:

| recid | primary | title |
|---|---|---|
| 1427769 | math.RT | The E8 Geometry from a Clifford Perspective |
| 1514749 | physics.gen-ph | Emergence of an aperiodic Dirichlet space … icosahedral internal space |
| 1707613 | math-ph | From the Trinity (A3,B3,H3) to an ADE correspondence |
| 1851687 | math.GR | Clifford Spinors and Root System Induction: H4 and the Grand Antiprism |
| 2705182 | cs.LG | Machine Learning Clifford Invariants of ADE Coxeter Elements |
| 2718815 | math.GR | The Isomorphism of H4 and E8 |

Confirming **none** attaches a physical parameter needs a Scout body-read of these
six, but the categories and venues (all Adv. Appl. Clifford Algebras / Mathematics /
RSPA-math) leave no hep-ph/hep-th phenomenology paper in the set. Scout's "broad
sample found none" → **exhaustive-at-Inspire-depth: still none.**

## Verdict summary

- §3.4 — holds (EMPTY), Inspire-exhaustive, ADS-pending.
- §4.4 — holds (EMPTY), Inspire-exhaustive, ADS-pending.
- §4.6 — **1 hit (recid 2876584)** needing a purpose/wording call; Inspire-exhaustive, ADS-pending.
- §4.8 — holds (no physical-parameter attacher), pending Scout read of 6 hep-th-crosslisted; Inspire-exhaustive, ADS-pending.
- §5.2 — deferred (not citation-testable).

## Reproduce

```
python inspire_crawl.py    # Inspire-HEP only; writes CSV + claim_findings.json + raw/
```

ADS leg (when a token is available): set the token and add a `citations(bibcode:…)`
query against `api.adsabs.harvard.edu/v1/search` for recids 1422619 / 1427769,
dedupe into the CSV, re-test the four intersections. Until then: **ADS-pending.**
