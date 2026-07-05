# Findings — Paper 118 Part 2: the solvable-control null (results note)

**Pre-registration:** `4d74af8` (Gate G1, frozen before this ran). **Substrate:** Part 1
machinery `96a290f` (AT1–AT9 green). `python solvable_control_null.py` → **landed on H_unique**.
Findings only — no interpretation (Mr A / CinC / Cliff).

## Verdict (one line)

Across the whole population of finite Γ < SU(2) (**91 groups: 90 solvable + 2I**), the barrier
predicate is **true for Γ = 2I and false for every other Γ**. The sweep lands on the
pre-registered **H_unique**. **H_arith is not realised** — no solvable Γ produces a barrier.

## Machinery-validation anchors (classical — passed)

First l > 0 with m_Γ(l) ≥ 1: **2T → 6, 2O → 8, 2I → 12**, all as required. Machinery validated;
run proceeded.

## Two-route exactness (Gate G3)

Route A (mpmath 60 dps character sum) and route B (exact integer: Molien for {2T,2O,2I},
weight-count for Z_n, `½[m_{Z_2n}+χ(π/2)]` for 2D_n) **agree for every l on every group**;
max integrality residual < 1e-30. Every reported multiplicity is confirmed both ways.

## The barrier predicate, per Γ (key rows; full table in `barrier_summary.csv`)

| Γ | \|Γ\| | solvable | l_gap | λ₁ | λ₁ vs \|Γ\| | p₁ | m(p₁) | Barrier | gap vs p₁ |
|---|---|---|---|---|---|---|---|---|---|
| Z₅ | 5 | yes | 2 | 8 | > | 10 | 3 | **False** | gap<p₁ |
| Z₁₀ | 10 | yes | 2 | 8 | < | 20 | 5 | **False** | gap<p₁ |
| 2D₅ | 20 | yes | 4 | 24 | > | 10 | 1 | **False** | gap<p₁ |
| 2D₁₅ | 60 | yes | 4 | 24 | < | 30 | 1 | **False** | gap<p₁ |
| 2T | 24 | yes | 6 | 48 | > | 6 | 1 | **False** | gap==p₁ |
| 2O | 48 | yes | 8 | 80 | > | 8 | 1 | **False** | gap==p₁ |
| **2I** | **120** | **no** | **12** | **168** | **>** | **10** | **0** | **True** | **gap>p₁** |

Only 2I: at its top-fibre perfect return (l = p₁ = 10, all vertex characters back to 1) the
mode is killed, m(10) = 0, and the gap is pushed up to l = 12, λ₁ = 168.

**Nearest solvable approach:** the odd-n binary dihedral family — H_arith's flagged candidate
— gets closest, m(p₁) = 1 for every 2D_{2j+1} (2D₃, 2D₅, 2D₇, …), but never crosses to 0.
The minimum m(p₁) over all 90 solvable groups is 1; only 2I reaches 0.

## Pre-registered confound (applies; stated at pre-reg)

Within this population non-solvability ⟺ Γ = 2I ⟺ stabiliser triple (2, 3, 5) — the three
properties coincide on a population of one. Uniqueness therefore **cannot separate** "A₅
simplicity" from "(2, 3, 5) coprime arithmetic" as cause. The dihedral family partially
decouples (many solvable groups, varying period arithmetic); full decoupling is impossible
here and none is claimed. (Per pre-reg; not a new interpretation.)

## EXPLORATORY (badged — NOT part of the H_unique/H_arith stake)

Turning §5's residue argument into a general computed statement across the family:

> **Barrier(Γ) ⟺ p₁ < l_gap** — verified for all 91 groups — i.e. the top-fibre perfect
> return undershoots the first invariant. Equivalently `2·(largest stabiliser order) <
> (lowest invariant degree)`: 2T gives 6 = 6, 2O gives 8 = 8, **2I gives 10 < 12** — the
> unique undershoot. For Z_n and 2D_n the gap sits far *below* the top-fibre return
> (`gap<p₁`), so by the time the top fibre returns the mode has long since survived.

The barrier is the phenomenon of a top-fibre return that lands strictly below the gap. Among
finite SU(2) subgroups only 2I is built that way.

## Files

- `PRE-REGISTRATION.md` — frozen stake (`4d74af8`).
- `solvable_control_null.py` — population, routes A/B, barrier predicate, sweep.
- `barrier_summary.csv` — all 91 rows (the per-Γ predicate table).
- `m_tables.csv` — full m(l), l = 0..max(120, 2·p₁), per Γ.
- `fibre_slots_at_p1.csv` — per-class character values at p₁ (the perfect-return record).

## Reproduce

```
python solvable_control_null.py   # anchors, verdict, CSVs
```

## Scope

Findings only. The meaning of H_unique for Paper 118 §6.2's simplicity reading, and any
§9/§8 flip in Papers 132/135, are **not** Mr Code's to decide (Mr A / CinC / Cliff / Flint).
The conic leg (Paper 132 solvable-null, `paper_132_solvable_null/`, 14/14) and this spectral
A₄/S₄-geometric leg are complementary legs of §9's own control sentence, confirmed
non-overlapping (W-107 grep clean).
