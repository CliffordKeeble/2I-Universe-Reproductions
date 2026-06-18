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

*(filled after this pre-registration is committed.)*
