# Paper 213 — ARM 2 (icosian / 2I imported): Reconstruction Test

**2I IS IMPORTED.** This arm couples by unit icosians (the 120 elements of 2I, the
imported object built and verified in `icosians.py`). It asks one question: *does a
discrete relational graph carrying genuine 2I holonomy reconstruct the continuum
S³/2I scalar even-sector* — the eigenvalues `λ_k = k(k+2)` supported on the numerical
semigroup `S(12,20,30)` (Paper 174), with spectral gap `λ₁ = 168` at `k=12` (Paper
108)? A success is a **faithful discretisation of a known 20-year-old spectrum**
(Luminet 2003 / Lachièze-Rey & Caillerie 2005 / Weeks 2005 / Kramer 2004), **not an
ignition.** Every figure and result carries the "2I imported" label.

---

## PART 1 — PRE-REGISTRATION (Pattern 75)

*Written and committed BEFORE the first production run. The commit timestamp is the
audit trail. None of the constants below may be tuned to the result.*

### What is fixed across the experiment

- **Topology engine (`growth.py`):** GR-1 self-reference growth, coupling-agnostic.
  Under `closure="no_constraint"` the topology is generated from the seeded RNG alone,
  so it is **identical for Arm 1 and Arm 2 at the same seed**. This is what makes the
  comparison clean: `L_s` (hence `d_s`) is a property of topology, not coupling
  (Arm 1 brief §4). `compare_arms.py` asserts `L_s` byte-equality before scoring.
- **Coupling = the only variable.** Arm 1: holonomy alphabet `{I, R}`, fibre ℝ²
  (`R = Γ_seed/5^¼`, the Paper-191 seed; **involution, det −1, R²=I, but NOT
  orthogonal** — `RᵀR = diag(√5, 1/√5)`; flagged, matches the parallel Arm 1 build).
  Arm 2: holonomy alphabet `{s, t} ⊂ 2I`, fibre ℝ⁴ (4×4 SO(4) icosian matrices,
  **2I imported**).
- **Read-out (`spectral.py`):** scalar `d_s` from the heat-trace of `L_s`; the
  `⟨12,20,30⟩` signature score applied to the connection Laplacian `L_c`.

### Pre-registered scoring constants

| constant | value | meaning |
|---|---|---|
| `M_CLUSTERS` | 8 | lowest nonzero `L_c` clusters scored |
| `CLUSTER_RELGAP` | 0.03 | new cluster when relative gap > 3% |
| `LANDING_TOL` | 0.05 | cluster centre within ±5% of a target `k(k+2)` |
| `ZERO_THRESH` | 1e-6 | eigenvalues below this are zero modes |
| `DS_FIT_WINDOW` | (0.05, 1.0)/λ_med | intermediate-t window for the `d_s` fit |

*(The parallel Arm 1 window pre-registered `M=10`; the cross-arm comparison here runs
BOTH arms through this directory's `spectral.py` with `M=8`, so the head-to-head is
internally consistent. The constant difference is cosmetic to the comparison.)*

### Continuum target (Paper 174 / Paper 108) — DERIVED, classical

- Supported even `k ∈ S(12,20,30)`: `{12, 20, 24, 30, 32, 36, 40, 42, 44, 48, 50, 52,
  54, 56, 60, …}` and every even `k ≥ 60`.
- Killed even `k` (the "missing modes"): `{2,4,6,8,10,14,16,18,22,26,28,34,38,46,58}`
  — exactly 15, largest 58.
- Eigenvalue ladder `λ_k = k(k+2)`: `168, 440, 624, 960, 1088, 1368, …`;
  `λ₁ = 168` (k=12); first ratio `λ₂/λ₁ = 440/168 = 2.619`.

### Signature score (pre-registered)

Rescale the `L_c` spectrum so the lowest nonzero cluster = 168, then
`signature_score = frac_supported − frac_forbidden` over the lowest `M=8` clusters
(fraction landing within ±5% of a supported `k(k+2)` minus fraction landing on a
killed `k(k+2)`). Range `[−1, 1]`; **1.0 = perfect even-sector reconstruction,
~0 = no signal.**

### PRE-REGISTERED VERDICT CRITERIA

Decided now, before the run:

- **RECONSTRUCTION SUCCESS** — `signature_score(Arm 2) → 1` as `N` grows (monotone
  improvement across the finite-size sequence) **AND** `signature_score(Arm 1) ≈ 0`
  (no signal) at every `N`. Interpretation: *faithful discretisation of the known
  S³/2I even-sector; the signature is **imported with 2I**, not emergent from golden
  coupling.* This is the predicted outcome and the experiment's real result. **It is
  NOT an ignition and does NOT trip the W-108 stop** (Arm 2 brief §6/§7).
- **PARTIAL** — Arm 2 score materially exceeds Arm 1 but does not converge to 1
  (discretisation captures some but not all of the ladder). A legitimate result.
- **NULL** — Arm 2 score indistinguishable from Arm 1 (≈0). The discretisation simply
  fails to converge; also a result. Report as such.
- **ANOMALY (escalate)** — if Arm 1 (`{I,R}`, ℤ/2 holonomy) shows a `⟨12,20,30⟩`
  signature, something carried 2I into the topology. STOP and hand to cold Adversary
  (W-108) — this would be the Arm 1 brief §7 trip, surfaced from the Arm 2 side.

### Controls (pre-registered)

1. **Arm-1 head-to-head** (the decisive control): same engine, same seed, same
   topology — only the coupling differs.
2. **Scrambled 2I:** random icosian per edge (not the `{s,t}` rule). Confirms the
   signature tracks 2I-import, not the specific growth rule.
3. **Bare vs connection:** `L_s` signature score vs `L_c` — `L_s` must show no
   signature (it is identical across arms).
4. **`d_max` sweep** `{3,4,5,6}`: confirm `d_s` and the score are not artefacts of the
   branching cap. Logged to `logs/rule_trials.jsonl`.

---

## PART 2 — RESULTS

*Run of record: `seed=20260618`, `python run_arm2.py` / `compare_arms.py` / `null_check.py`.
Artefacts: `results.json`, `logs/rule_trials.jsonl`, `figures/arm2_convergence.png`.*

### Headline — NULL (no reconstruction), and it is the substrate, not the coupling

> **Under GR-1, the 2I-imported connection Laplacian does NOT reconstruct the
> `⟨12,20,30⟩` even-sector. The grown graph is sub-1D (`d_s ≈ 0.5`), not a
> 3-manifold — so the reconstruction premise is structurally unmet by this growth
> rule. Neither arm shows a real signature once each is tested against its OWN
> matched null. No W-108 trip.** `[OBSERVED]`

This is the brief's pre-registered **NULL / PARTIAL** branch (§6): *"the
discretisation may simply fail to converge — also a result."* It is **not** an
ignition and does **not** trip W-108.

### Spectral dimension (shared `L_s`, identical across arms)

| N | `d_s` ± err | d_max | `d_s` ± err |
|---|---|---|---|
| 1000 | 0.545 ± 0.031 | 3 | 0.512 ± 0.030 |
| 2000 | 0.535 ± 0.032 | 4 | 0.535 ± 0.032 |
| 4000 | 0.541 ± 0.033 | 5 | 0.539 ± 0.032 |
| | | 6 | 0.539 ± 0.033 |

`d_s` is **flat across N and across the branching cap d_max** — far below 3, with no
trend toward it. `[OBSERVED]` *(Caveat: the absolute value is biased low by spectrum
truncation — only the lowest ~250 eigenvalues enter the heat-trace fit; the robust,
truncation-insensitive content is "stable, nowhere near 3, insensitive to d_max",
not the precise figure 0.5.)*

### Signature score vs each arm's OWN matched null (the decisive evidence)

`L_c` low spectrum rescaled so lowest cluster = 168, scored on the `k(k+2)` ladder
(pre-registered, `M=8`, ±5%). Matched null = scrambled coupling on the same topology.

| arm (N=2000) | structured score | matched null (mean ± σ, max) | separation | verdict |
|---|---|---|---|---|
| **Arm 2 — 2I imported `{s,t}`** | **−0.167** | +0.017 ± 0.073, max +0.333 | **z = −2.5, 0th pct** | **below its own scramble null → NO reconstruction** |
| Arm 1 — seed `{I,R}` | +0.750 | **+0.500 ± 0.389, max +1.000** | z = +0.6, 80th pct | **inside its own null → artefact, NOT a signature** |

The raw `compare_arms` table looked like the pre-registered ANOMALY (Arm 1 high, Arm 2
low — the *opposite* of the prediction). The matched-null test dissolves it: a random
`{I,R}` coupling scores +0.5 on average (the 2-D abelian structure lands on the dense
upper `k(k+2)` ladder regardless of content), so Arm 1's +0.75 carries **no** `2I`
content. **Nothing carried the icosahedron into the topology; no W-108 escalation.**
`[OBSERVED]` Arm 2's low clusters `[168, 209, 216, 230, 243, 321]` show no ladder at
all (target `168, 440, 624, 960, …`).

### Rule-space map (`logs/rule_trials.jsonl`)

- **d_max sweep {3,4,5,6}:** `d_s` ≈ 0.51–0.54 throughout; Arm 2 score ∈ {−0.20, −0.17,
  0.00, +0.33}, no trend — noise around the null.
- **closure {no_constraint, flat_triangles, flat_shortest}:** `d_s` ≈ 0.48–0.54; Arm 2
  score ∈ {−0.17, +0.33, −0.25}. `flat_shortest` is holonomy-gated so it breaks the
  cross-arm topology match (expected; excluded from the head-to-head).
- **ρ {parity, hash}:** Arm 2 score −0.167 / 0.000 — the context-hash into the full 120
  does not help. No variant produces a signature.
- **N=4000:** Arm 2 score `nan` — near-degenerate low spectrum at the 16000-dim `L_c`
  with only 250 eigenvalues requested. A pipeline limitation, flagged; not load-bearing
  (the null is already clear at N=1000/2000). `[FLAG]`

### Methodological finding — the pre-registered score is asymmetric in power

The score is **informative for Arm 2** (random 2I scrambles score ≈ 0, so a real
reconstruction would stand out) but **vacuous for Arm 1** (any `{I,R}` pattern scores
~0.5–1.0). Root cause: the supported semigroup ladder `k(k+2)` is co-dense with the
reals for `k ≳ 40` (consecutive allowed `λ` differ by `< 5%` tolerance), so the upper
ladder cannot discriminate; only the well-separated low rungs (168, 440, 624, 960) can.
A sharper future score should weight only `k ∈ {12,20,24,30}` and penalise the low
forbidden rungs `{8,24,48,80}`. `[OBSERVED — recommend score v2 before any rescale-up]`

### Honest framing — situated against the literature (Scout brief)

Weeks (2005), Lachièze-Rey & Caillerie (2005) and Kramer (2004) constructed the exact
eigenmodes of the **actual 3-manifold** S³/2I; the `⟨12,20,30⟩` even-sector is theirs,
classical (Paper 174). Our GR-1 graph is **not a discretisation of that manifold** —
its spectral dimension never leaves ≈0.5. So the honest statement of *what golden /
seed-style relational growth can and cannot do* is:

- **It does not build the geometric substrate.** GR-1 self-reference growth yields a
  thin, sub-1D graph, not a 3-manifold — independent of coupling, d_max, or closure.
- **Importing 2I as holonomy on a non-3D graph does not conjure the 3-manifold
  spectrum.** The signature is absent in Arm 2 exactly as it is in Arm 1; the icosahedron
  lives in the *group* and in the *manifold*, and neither is reconstructed by hanging 2I
  holonomies on a stringy graph. A discrete echo of Paper 174's point that the spectrum
  comes from the **quotient of a 3-sphere**, not from aperiodic/relational growth.

### Verdict (pre-registered mapping)

- **NULL** for the reconstruction, with the stronger diagnosis that the GR-1 substrate
  is not 3-dimensional. `[OBSERVED]`
- **No W-108 trip** — Arm 1's apparent signature is a matched-null artefact. `[OBSERVED]`
- **2I remains IMPORTED throughout** — nothing emerged; the test simply did not
  reconstruct.

### Open / handed up to CinC (NOT decided here — borders W-107)

The natural next move — search for a growth rule with `d_s ≈ 3` and re-test — is **not
taken unilaterally**: tuning the build until it both reaches 3D *and* shows the signature
is exactly the W-107 "manufacture the spectrum from inside the build chain" trap. Whether
to (a) accept this NULL for GR-1 as the result, (b) commission a *pre-registered*
3D-substrate growth rule as a distinct experiment, or (c) redesign the score (v2 above)
and rescale to N=10⁵–10⁶, is CinC's call.
