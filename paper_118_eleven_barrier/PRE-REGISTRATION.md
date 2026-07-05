# PRE-REGISTRATION — Paper 118 Part 2: the solvable-control null

**Brief 118-CODE-1, Part 2 (Amendment A1: Python).** Gate G1: this file is committed
as its own commit **before** the sweep executes; the results note cites its commit hash
and never absorbs it. Frozen once staked.

**Substrate.** The machinery is the Part 1 reproduction (commit `96a290f`, AT1–AT9 green):
`m_Γ(l) = (1/|Γ|) Σ_C |C| χ_l(θ_C)`, `χ_l(θ) = sin((l+1)θ)/sin(θ) = U_l(cos θ)`, verified on
the 2I instance three independent ways. Part 2 runs the same machinery across the whole
population of finite Γ < SU(2).

## Population (locked)

- Cyclic `Z_n`, 2 ≤ n ≤ 60 (|Γ| = n).
- Binary dihedral `2D_n`, 2 ≤ n ≤ 30 (|Γ| = 4n).
- Binary tetrahedral `2T` (24), binary octahedral `2O` (48), binary icosahedral `2I` (120).

Every Γ except 2I has solvable image in SO(3); 2I is the unique non-solvable case
(image A₅). For each Γ: assert Σ class sizes = |Γ| and m_Γ(0) = 1 before use.

## Machinery-validation anchors (locked, classical — not hypotheses)

First l > 0 with m_Γ(l) ≥ 1 must be **6 for 2T, 8 for 2O, 12 for 2I** (invariant degrees
6/8/12, 8/12/18, 12/20/30). If any anchor fails, the machinery is wrong — **stop and report**,
run no further.

## The barrier predicate (locked)

For each Γ, `p₁` is the period of the bounded fibre of longest period:
- polyhedral Γ: `p₁ = 2k`, k = largest stabiliser order — 2T: 6, 2O: 8, 2I: 10;
- `Z_n`: `p₁ = 2n`;
- `2D_n`: `p₁ = 2n` (top fibre; against two period-4 fibres).

**Barrier(Γ)** := at `l = p₁` — the first perfect return of the top fibre, all its class
characters returning to their identity value (χ = 1) — `m_Γ(p₁) = 0`.

For every Γ report: |Γ|; the gap `l_gap` (first l > 0 with m ≥ 1) and `λ₁ = l_gap(l_gap + 2)`;
λ₁ vs |Γ|; the top-fibre return `p₁`; the fibre-slot (per-class character) values at p₁;
Barrier(Γ) true/false; and where the gap sits relative to p₁.

## Pre-registered hypotheses (the stake)

- **H_unique:** Barrier(Γ) holds for Γ = 2I and for no other Γ in the population.
- **H_arith:** Barrier(Γ) holds for at least one solvable Γ (the odd-n binary dihedral family
  is the natural place a coprime-clash could occur — but **no expected outcome is registered
  for it; compute, don't assume**).

## Pre-registered interpretations (staked before the run)

- If **H_arith**: the barrier phenomenon occurs without simplicity, and Paper 118 §6.2's
  simplicity reading is falsified *as stated* — an honest result its own §6.1 concession
  anticipates.
- If **H_unique**: uniqueness is consistent with §6.2 **but attributes nothing**, because of
  the confound below.

## Confound (stated at pre-reg, not discovered after)

Within spherical space forms, non-solvability ⟺ Γ = 2I ⟺ stabiliser triple (2, 3, 5): the
properties coincide on a population of one, so uniqueness **cannot separate** "A₅ simplicity"
from "(2, 3, 5) coprime arithmetic" as cause. The binary dihedral family is the partial
decoupler (many solvable groups, varying period arithmetic); full decoupling is impossible
inside this population and **no claim of it may be made**.

## Decision rule (locked)

The sweep lands on **H_unique** iff Barrier is true for 2I and false for every other Γ; on
**H_arith** iff Barrier is true for at least one solvable Γ. Report the verdict as findings
only — no interpretation (that is Mr A / CinC / Cliff). Any anchor failure aborts the run.

## Exactness gates (locked, G3)

Two independent routes per load-bearing number (the anchors and every barrier verdict):
- route A — mpmath ≥ 50 dps character sum, integrality assertion `|m − round(m)| < 1e-30`;
- route B — exact integer: Molien series for {2T, 2O, 2I}
  ((1+t¹²)/((1−t⁶)(1−t⁸)), (1+t¹⁸)/((1−t⁸)(1−t¹²)), (1+t³⁰)/((1−t¹²)(1−t²⁰))),
  weight-counting for `Z_n`, and `m_{2D_n}(l) = ½[m_{Z_{2n}}(l) + χ_l(π/2)]` for `2D_n`.
Routes A and B must agree for every l on every Γ.

## Exploratory (flagged, NOT pre-registered)

Characterise computationally which period-multisets produce a killed first return — turning
§5's residue argument into a general computed statement across the family. Reported separately,
badged **EXPLORATORY**, outside the H_unique/H_arith stake.
