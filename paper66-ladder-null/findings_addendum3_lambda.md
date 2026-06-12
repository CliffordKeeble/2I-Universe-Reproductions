# Findings — Addendum 3: the Λ_QCD rival sweep (`lambda_sweep.py`)

**2I Universe Programme · paper66-ladder-null · 11 June 2026**

## Provenance — disclosed plainly

Numbers originated in the briefed Adversary pass on Paper 66 (Mr A's own steelman
run, J2 of the joint list); reproduced independently at **Tor's bench** before
absorption (this script, candidates pre-registered in the docstring before the run);
cited in the paper at §3.2 with the phrase "commit pending" — which this check-in
makes true, after which the canonical Mr Library build strikes the phrase. Per W-074
the repo's run of record is **Mr Code's verification re-run**.

## Result (verified twice: briefed pass + Tor's bench)

| Λ_QCD (MeV) | length ħc/Λ (fm) | N = length/u | nearest int (dist) |
|---|---|---|---|
| 210 | 0.9397 | 17.872 | 18 (0.128) |
| 250 | 0.7893 | **15.012** | **15 (0.012)** |
| 263 | 0.7503 | 14.270 | 14 (0.270) |
| 290 | 0.6804 | 12.942 | 13 (0.058) |
| 332 | 0.5944 | 11.304 | 11 (0.304) |
| 340 | 0.5804 | 11.038 | 11 (0.038) |

**Verdict (extends rule R4):** the nearest physical rival length sweeps ≈11–18 units
across standard scheme choices, landing a spuriously sharp 15.012 at one of them —
sharper than the 19-row's 0.02. Definitional promiscuity demonstrated on the rival;
the §3.2 no-claim demotion of the 19 is strengthened, not weakened.

## Mr Code check-in brief — addendum 3

Directory `paper66-ladder-null/`, three files, two commits, same gates as the
five-file check-in:

**Commit 1 (declaration):** `lambda_sweep.py` + `findings_addendum3_lambda.md`
(this file). Script unaltered from the bench run, hash-matched.

**Verification:** re-run the script; stdout must reproduce
`lambda_sweep_output.txt` byte-for-byte (line endings aside). Output is ASCII.

**Commit 2 (verified output):** `lambda_sweep_output.txt`, with the one-line
run-of-record note in the commit message.

**Hard gates:** (1) script unaltered + hash-matched; (2) output committed only on
exact reproduction; mismatch = stop and flag. No other repo writes — the §3.2
"commit pending" strike belongs to the Mr Library build, not to this check-in.

🐕⚛️
