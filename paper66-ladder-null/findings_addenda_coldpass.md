# Findings — Cold-Pass Addenda (robustness + the 19)

**2I Universe Programme · paper66-ladder-null · 10–11 June 2026**
*Two addendum computations run in response to the two cold passes on Paper 66 v2.0.
Provenance disclosed below; verdicts per rules declared in each script's header
before its run.*

## Provenance — disclosed plainly

Both addenda were run at **Tor's bench** (not Mr Code's) to feed a live revision
between cold passes: the pre-registration of each lives in its script docstring,
written before the run; the outputs were captured verbatim and are already cited in
the v2.0-w3 draft (§5.1 and §3.2). Per W-074 the work is not real until merged: the
repo's run of record is **Mr Code's verification re-run**, which must reproduce both
output files byte-for-byte. Any mismatch stops the check-in and comes back as a flag,
not a silent reconciliation.

## Addendum 1 — class-boundary robustness (`robustness.py`)

First cold pass, fifth-star condition 1: show DECORATIVE is not an artefact of a
loose class. Settings declared before running: T1 = (≤2 terms, |k| ≤ 3);
T2 = (≤3 terms, |k| ≤ 3); T0 = original reference.

| Setting | Coverage [30,200] | Mean spellings | 44 | 137 |
|---|---|---|---|---|
| T1 (≤2 terms, \|k\|≤3) | 76.6% | 8 | 12 | **0** |
| T2 (≤3 terms, \|k\|≤3) | 100.0% | 242 | 339 | 63 |
| T0 (≤3 terms, \|k\|≤7) | 100.0% | 2,021 | 2,762 | 828 |

**Verdict:** DECORATIVE survives every tightening that still contains the paper under
audit. Honest wrinkle, in the paper: at T1, 137 itself becomes unexpressible — the
class was widened to contain v1.1's own formulas (7F − D needs |k| = 7), not to
manufacture promiscuity. Cited at draft §5.1.

## Addendum 2 — the 19 through the null (`null_19.py`)

Second cold pass, fifth-star condition 1: turn the sword on the paper's favourite
integer. Decision rule R4, declared before running: a count is bankable only if
anchored to a **physically defined** length; a match produced only by a round decimal
figure is a fact about the metric system and is carried at NO CLAIM
(definition-dominated — the cosmic-row rule).

| Candidate for "the range" | r (fm) | N = r/u | Nearest int (dist) |
|---|---|---|---|
| Proton rms charge radius | 0.8409 | 15.99 | 16 (0.01) — the §3.2 lock, already booked |
| Charged-pion Yukawa ħ/m_π c | 1.4138 | 26.89 | 27 (0.11) |
| Neutral-pion Yukawa | 1.4619 | 27.81 | 28 (0.19) |
| **Decimal convention 1.000 fm** | 1.0000 | **19.02** | **19 (0.02)** |

In-class spellings of the integer 19 under class C: **1,324** — the same promiscuity
that convicted the ladder (neighbours 12–29 range from 1,297 to 5,048).

**Verdict (R4):** only the round decimal figure lands on 19; the physically defined
candidates land on the existing lock (15.99) or on nothing (26.9, 27.8). The 19 is a
fact about the metric system's round number, not about confinement → **NO CLAIM,
definition-dominated.** Cited at draft §3.2; the §3.4 table row amended to match.

## Mr Code check-in brief

Directory `paper66-ladder-null/`, five files, two commits:

**Commit 1 (declaration):** `robustness.py`, `null_19.py`,
`findings_addenda_coldpass.md` (this file) — the scripts carry their pre-registrations
in their docstrings; this file discloses the bench-first provenance.

**Verification (between commits):** re-run both scripts; stdout must reproduce
`robustness_output.txt` and `null_19_output.txt` **byte-for-byte** (line endings
aside). Both scripts are believed ASCII-safe on output per the v1.1 console ruling
(executables print ASCII; prose may dress) — if the console disagrees, flag it, do
not patch the committed script.

**Commit 2 (verified outputs):** `robustness_output.txt`, `null_19_output.txt`, plus
a one-line verification note in the commit message naming the re-run as the run of
record.

**Hard gates:** (1) scripts committed unaltered from the bench run; (2) outputs
committed only if the verification re-run reproduces them; mismatch = stop and flag.
No other repo writes.

🐕⚛️
