# B1.6-4b — the evaluation point, derived from antilinearity

**Paper 226 · Golden Construction · Part B · derivation note (committed BEFORE re-execution)**
**Status: derivation note · 6 July 2026 · internal**

This note fixes the evaluation point of the factorisation-frontier target operator **as a
consequence of the map's definition**, before any re-run. It contains no reference to whether a
charge-conjugation operator exists, to any dimension, or to any pass/fail outcome. The point is
that the evaluation momentum is determined by *which field automorphism a candidate carries*, not
by what result that choice produces.

## Setup

A charge-conjugation candidate on the plane-wave sector has the shape **C = B ∘ g**, where B is a
constant matrix part and g ∈ {1 (linear), τ, σ} is the field automorphism it carries. It is asked
to intertwine the momentum-k, charge-q operator L_q with the charge-(−q) operator: C · L_q = L_{−q} · C
on solutions. The solutions are plane waves with factor e^{i(k·x)} (here k·x = k₀x₀ + k₁x₁; k, x real).

## The derivation

A field automorphism acts on the plane-wave factor through its action on the imaginary unit i:

- **Linear (g = 1):** no field map. e^{i k x} is unchanged. Momentum stays **+k**.
- **σ-borne (g = σ):** σ sends √5 ↦ −√5 and **fixes i**. The plane-wave factor carries i, not √5,
  so e^{i k x} is unchanged. Momentum stays **+k**.
- **τ-borne (g = τ):** τ is antilinear — it conjugates the imaginary unit, i ↦ −i. Therefore

      τ(e^{i k x}) = e^{−i k x},

  which is the plane wave of momentum **−k**. A τ-borne candidate carries the k-momentum sector to
  the −k sector.

Consequently the target operator against which a candidate of type g is tested is evaluated at

      g = 1   → +k        g = σ   → +k        g = τ   → −k,

each read off from g's action on e^{i k x}. For the τ-borne candidate the target sits at −k **as a
consequence of antilinearity alone** — it is fixed the moment one says "the candidate carries τ",
and does not depend on the outcome of the subsequent search for B.

## Consequence for the re-run

In the re-run (B-ii) the evaluation momentum is not supplied by hand. It is obtained by applying
the candidate's field automorphism to e^{i k x} and reading the resulting momentum — so the −k of
the τ-borne row is produced by the code from τ(e^{ikx}) = e^{−ikx}, not chosen. The linear and
σ-borne rows are evaluated at +k by the same procedure.

— committed before re-execution, per standing law. 🐕☕⬡
