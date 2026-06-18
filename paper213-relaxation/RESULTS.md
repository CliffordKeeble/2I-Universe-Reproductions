# Paper 213 — Approach B (Relaxation): does a frustrated field relax to round S³?

**Approach B**, the complement to the failed Approach A (GR-1 growth, `d_s ≈ 0.5` — a
local rewrite cannot grow the closure). Here the closed-3-sphere **topology is GRANTED**
("first closure" axiom); we test whether the **geometry relaxes to round S³** and a
chirality-pumped field relaxes to the **+2 Hopf Beltrami** steady state.

- **Main arm (hard):** a *generic* combinatorial 3-sphere, **2I un-imported**.
- **Control arm:** the 600-cell ({3,3,5}, 120 vertices = icosians), **2I imported, labelled**.
- **Anti-circularity gate (§4):** any 2I structure on the generic arm → STOP (import leak;
  the `z=+13.9 → +0.64` Arm-2 episode is the cautionary tale).

Staged, debugger-style (one variable at a time): Stage 1 geometry alone → Stage 2 field
on round S³ → Stage 3 coupled. **Stage N runs only after Stage N−1 passes.**

---

## STAGE 1 — PRE-REGISTRATION (committed BEFORE the run; hard gate §5)

**Question:** given the closed-S³ topology and a frustrated initial metric (all edge
lengths `ℓ_e = 1` on a generic triangulation), does **normalized (volume-fixed) discrete
Ricci flow** relax the geometry to **round S³** (constant curvature)?

**Scheme (frozen):** curvature-homogenizing combinatorial Ricci flow (Regge + Chow–Luo).
Per-edge discrete sectional curvature `K_e = δ_e / a_e`, with dual area
`a_e = V_e^dual / ℓ_e` and barycentric dual `V_e^dual = (1/6) Σ_{t⊃e} V_t`. Flow
`d(log ℓ_e)/dτ = −(K_e − K̄)` (K̄ = volume-weighted mean), total volume rescaled to a
constant each step. Fixed point `K_e = K̄ ∀e` = round S³. *(Literal gradient-descent on
the Regge action `S = Σ ℓ_e δ_e` is a saddle at the round metric — conformal-mode problem
— so the brief-permitted combinatorial Ricci flow is used. Sign = standard homogenizing
direction, fixed before the run; no metric/sign tuning.)*

**Convergence metric (frozen):** curvature coefficient of variation
`CoV(K) = std(K_e) / |mean(K_e)|` (volume-weighted). `CoV → 0` ⇒ round.

**PRE-REGISTERED PREDICTION:** `CoV(K)` drops substantially and plateaus low
(target: ≥ 5× reduction from the frustrated initial value, and `CoV_final ≲ 0.3`),
i.e. the flow reaches *approximately* round S³ (a fixed generic triangulation may be
resolution-limited from reaching `CoV = 0` exactly — finer mesh, **not** a tuned fix).

**PRE-REGISTERED NULL / FAILURE:** curvature localizes into a defect tangle (discrete
neckpinch: a tet becomes non-embeddable, or `CoV` stays high / rises). A real result —
it would call for a Perelman-surgery analogue or finer mesh, **not** a tuned patch.

**Controls for Stage 1:** (3) already-round init stays round (sanity — deferred to the
600-cell/analytic round check). **Anti-circularity:** confirm the substrate is 2I-free
(`looks_2I_regular = False`) before flowing.

**Status tags** on every reported quantity: DERIVED / OBSERVED / CONJECTURED / STRUCTURAL.

---

## STAGE 1 — RESULTS

*Run of record: `seed=20260618`, `python ricci_flow.py 150 2000`. Artefact:
`results_stage1_N150.json`.*

### Verdict: NO VALID VERDICT — the flow instrument FAILED its own sanity control

> **The pre-registered normalized-Ricci-flow scheme is unfaithful: it neckpinches an
> already-round metric. Control #3 invalidated the instrument before any Stage 1 verdict
> could be banked. The frustrated-arm "neckpinch" is therefore an artefact of the flow
> scheme, NOT a geometric result. Stage 1 is HELD pending a faithful flow.** `[OBSERVED]`

This is the discipline working — a mirror of the Arm-2 `z=+13.9 → +0.64` episode the
brief names: a result that *looked* like the pre-registered null (neckpinch/glass) is
suspect until the control validates the instrument. It did not.

### What IS validated (substrate + curvature diagnostic)

| metric (N=150) | initial CoV(K) | initial K̄ | reading |
|---|---|---|---|
| **round-init** (chord lengths on actual S³ points) | **0.73** | **+3.56** | low scatter, positive curvature — *correctly reads as round* |
| **frustrated** (all ℓ_e = 1) | **55.7** | **−0.32** | high scatter, net-negative — *correctly reads as far-from-round* |

So `geometry.py` (substrate, Regge deficits, the curvature measure `K_e = δ_e/a_e`) is
sound: a round metric registers round, a frustrated one registers frustrated. The
**floor** achievable on this coarse mesh is the round-init value `CoV ≈ 0.73`, not 0 —
my pre-registered `CoV_final ≲ 0.3` target was too strict for N=150 (a discretisation
resolution effect; higher N lowers the round-init floor). `[OBSERVED]`

### What FAILED (the flow)

| flow from | CoV(K): start → end | outcome |
|---|---|---|
| frustrated ℓ=1 | 55.7 → 21.9 | neckpinch (min tet vol 0.118 → 0.010 over 5 steps, accelerating) |
| **round-init** | **0.73 → 2.71** | **neckpinch — CoV got WORSE; an already-round metric should stay round** |

The flow `d(log ℓ_e)/dτ = −(K_e − K̄)` with `K_e = δ_e/a_e` (barycentric dual) does
**not** have the round metric as a stable fixed point — it drives even round S³ into a
degenerate tet, and `K̄` runs away during the frustrated flow (−0.32 → −3.31). The crude
barycentric dual area and the explicit-Euler integration make the scheme unfaithful as a
discrete Ricci flow. This is a scheme problem, not a step-size problem (smaller `cap` did
not save it; the round-init failure is the proof). `[OBSERVED]`

### Stop point and proposed fix (CinC's call — non-trivial method change)

Per discipline, I am **not** swapping the pre-registered scheme on my own initiative.
The flow needs to be made faithful first; candidate fixes, in order of preference:

1. **Fixed-volume curvature-variance minimisation** — minimise `E = Σ_e w_e (K_e − K̄)²`
   at fixed total volume with a robust optimiser (L-BFGS). Round S³ is the global min, so
   round-init stays round by construction; the question becomes whether the frustrated
   metric reaches the same low-CoV basin or gets stuck (glass) — a clean, stable test.
2. **Glickenstein / circumcentric-dual discrete Ricci flow** — replace the crude
   barycentric `a_e` with a geometrically rigorous dual; round becomes a stable fixed
   point. Heavier to implement, closer to the literal "Ricci flow."
3. **Implicit (backward-Euler) integration** of the current flow — may damp the
   instability, but won't fix an incorrect fixed point if the dual is the root cause.

I recommend (1) as the crude-and-honest first move, re-running Control #3 as the gate.
**Stage 2 (field) and Stage 3 (coupled) remain gated and untouched** — Stage 1 has not
passed.
