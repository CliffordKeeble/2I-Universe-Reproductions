# PRE-REGISTRATION — Solvable-Case Null (Withholding Test)

**Pattern-75 falsifier. Design locked by Flint's brief (3 July 2026) before any computation;
this document records the locked design. Reproduction run: Mr Code, 4 July 2026.**

**Origin:** Mr A's cold review of the 132/135 pair (pair-level NULL / Pattern-75 finding). The
programme waves the solvable cases (conic, A₄, S₄) away as "too simple" *by stipulation*. Mr A:
"A framework that cannot be shown to fail where it should fail is not being tested." This test
converts the stipulation into a run.

## What is under test

Run the framework's weight-2 selection machinery on **solvable / genus-0 controllers** and
**show** — not stipulate — that it produces **no rank-like selection**, only a Hasse–Minkowski
yes/no; and that the same machinery **does** produce a rank-tracking signal on the **genus-1**
grammar curves. Withholding where it should withhold; selecting where it should select.

**Scope (locked, so the green tick is not oversold):** passing shows the selection machinery is
**not vacuous** — it does not manufacture a rank where there is no β₁ topology. It does **not**
earn the independent identification of the survivor, and it is **not** the divergent-prediction
falsifier. Both remain open. This answers exactly one question: *does the framework withhold on
solvable controllers when run honestly?*

## The mathematics it rests on (verified on the bench before results — see findings)

1. A smooth conic over ℚ (genus 0) has trivial Jacobian ⇒ no Mordell–Weil group of positive rank.
   Over ℂ it is P¹(ℂ) = S², Betti (1, 0, 1), so **β₁ = 0** — no loop register for the quintic
   filter to read. **[DERIVED]**
2. The only arithmetic output of a conic is binary: does it have a ℚ-point? By **Hasse–Minkowski**,
   a ternary form is isotropic over ℚ iff over ℝ and every ℚ_p; only ∞, 2, and odd p | 2abc can
   obstruct — decidable and finite. **[DERIVED]**
3. The genus-1 grammar curves have β₁ = 2, Jacobian = the curve, Mordell–Weil rank 0/1/2. The
   **Mestre–Nagao sum** S(X) = (1/log X) · Σ_{good p≤X} (−a_p·log p)/p tracks the rank (Nagao's
   conjecture — **pending Scout ticket S1**; stated as the statistic, not cited as settled).
   **[OBSERVED]**

## Inputs (locked)

**Positive controls (genus 1, β₁ = 2) — Weierstrass coefficients [a₁,a₂,a₃,a₄,a₆]:**

| curve | Weierstrass | expected rank |
|---|---|---|
| 11a1  | [0,−1,1,−10,−20] | 0 |
| 37a1  | [0,0,1,−1,0]     | 1 |
| 389a1 | [0,1,1,−2,0]     | 2 |

**Null cases (genus 0, β₁ = 0) — conics x²+y²=n, i.e. the form [1,1,−n]:**
`n = 1` (soluble), `n = 3` (insoluble), `n = 5` (soluble), `n = 21` (insoluble), plus `n = 6, 7`
to exercise both local flavours (an obstruction at 2 and at an odd prime).

## Procedure (locked)

- **Nulls:** assert genus 0 and β₁ = 0; run Hasse–Minkowski, emit the **binary** result and the
  obstruction places; assert the pipeline returns **no rank register** — the output is *only* the
  binary. **Falsification:** if a non-trivial rank-like value is returned on a genus-0 controller,
  the withholding test **FAILS**.
- **Positive controls:** compute a_p for good p ≤ 10⁴ by point counting on the minimal model;
  compute S(10⁴); assert it **separates the ranks in order (r0 < r1 < r2)**; assert the β₁
  register is filled (a rank-tracking signal exists where there *is* topology).

## Decision rule (locked)

PASS iff: (a) all conics give the correct HM binary with **no** rank selection; (b) the three
grammar curves' Nagao sums separate `r0 < r1 < r2`. The deliverable is the **contrast** —
withholding on genus 0, selecting on genus 1.

## Two hard gates (locked)

1. **LMFDB is a test oracle, never an input.** Ranks, Nagao values, solubility — computed from
   Weierstrass coefficients / the conic equation alone.
2. **Minimality guards.** The three grammar-curve models are minimal (LMFDB-minimal); bad primes
   are those dividing Δ, computed here, not looked up.

## On the targets

The manifest's Nagao figures (−0.26 / +0.75 / +2.06) were **from-memory quotes** this run
replaces with a real computation. They are **retired**. The reproduced S(10⁴) values are frozen as
**Pattern-75 regression oracles** (reproduce-this-number), **not mathematical targets**: the
mathematical claim is the **separation**, and — as an open question, not asserted — trend toward
the integer rank. Nothing hangs on any single value.
