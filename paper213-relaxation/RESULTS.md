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

---

## STAGE 1 v1.1 — PRE-REGISTRATION (committed BEFORE the forcing run; hard gate §5)

*Per brief v1.1 (Mr Adversary cycle-1 + the v1.0 instrument finding). Supersedes the
v1.0 Stage 1 above; the v1.0 record is kept as honest history.*

**Flow (decided):** Option 1 — fixed-volume **curvature-variance minimisation**,
`min_ℓ E = Σ_e w_e (K_e − K̄)²` (`w_e = V_e^dual`), volume held fixed by renormalisation,
solved by L-BFGS-B on `log ℓ` over the vectorised energy. Round S³ is the global minimum,
so **Control #3 (round-init stays round) holds by construction** — the fix for the v1.0
neckpinch-of-a-round-metric instrument failure.

**VERDICT = the SPECTRUM, not the variance (Mr A #2):** the relaxed scalar-Laplacian
spectrum (FEM Laplace–Beltrami, `spectral_s3.py`, validated to reproduce `k(k+2)` on
round-init) must approach round S³. Operationalised as the **spectral distance to the
round-init reference on the same mesh**: `dist_init → dist_relaxed`. Variance/CoV is the
engine and a diagnostic only (CoV is unreliable near `K̄=0` and a glass can be locally
uniform — hence the spectrum judges).

**Stage 0 gate (committed, Mr A #4):** generic complex must be 2I-free —
WL-colour-refinement bound `|Aut| = 1` (all colour classes singletons). Abort Stage 1 if
it fails.

**Controls:** (#3) round-init stays round — by construction; **(NULL) defective mesh**
(clustered points, sliver tets, known-can't-reach-round) — its relaxed spectrum must
*not* approach round; this gives "the variance/spectrum improved" its discriminating
teeth.

**PREDICTION — CONJECTURED (Mr A #3; the discrete flow cannot borrow Hamilton's smooth
guarantee):** on the generic arm, `E` drops substantially and the relaxed spectral
distance to round-init falls toward the N-dependent round-init floor; on the defective
null it does not.

**NULL:** generic arm traps with a spectrum far from round (glass) — a real result,
calling for finer mesh / surgery, **not** a tuned fix.

**Anti-circularity (§4):** the relaxed generic spectrum must show NO 2I-sieve structure
(no ⟨12,20,30⟩ missing-even-degree pattern). If it does → STOP (leak hunt).

**FSS observable (Mr A #7):** the round-init spectral floor and the relaxed
spectrum-match error → 0 as N grows (reported where tractable).

---

## STAGE 1 v1.1 — RESULTS

*Run of record: `seed=20260618`, `python run_stage1.py 100 50` (artefact
`results_stage1_v11_N100.json`). A full-budget re-run (`iters=250`) reached the
IDENTICAL generic endpoint (E→61.5, dist 0.249) and L-BFGS reported a local minimum
(line-search exhausted) — so the result is CONVERGED, not budget-limited.*

### Verdict: NULL — GLASS (variance minimises, spectrum does NOT reach round) `[OBSERVED]`

> **The instrument is now sound (Control #3 holds: round-init stays round). On the
> generic frustrated metric, fixed-volume curvature-variance minimisation drives the
> energy down ×5.6 to a local minimum — but the scalar-Laplacian SPECTRUM does not
> converge to round S³ (dist-to-round 0.280 → 0.249, ~11% of the way). Variance
> collapses; the spectrum stays glassy. This is the pre-registered NULL, and it is
> exactly the failure mode Mr A #2 built the spectrum-verdict to catch.** A fixed
> generic triangulation does not relax to round by edge-length flow alone.

### The numbers

| arm (N=100) | E: start → end | CoV | spectral dist-to-round: init → relaxed | verdict |
|---|---|---|---|---|
| **GENERIC** (test) | 342.7 → 61.5 (**×5.6**) | 275 → 7 | **0.280 → 0.249** | variance down, **spectrum NOT round** → glass |
| **round-init** (Control #3) | 8.42 → 3.05 | 0.73 → low | 0 (reference) | **stays round ✓ — instrument valid** |
| **DEFECTIVE** (null) | 326.3 → 72.5 (×4.5) | 156 → 8 | 0.150 → 0.123 | also glasses |

- **Stage 0 gate:** PASS — `|Aut| = 1` (all 100 WL colour classes singleton), 2I-free.
- **Anti-circularity (§4): CLEAN** — generic relaxed spectrum `[3.0, 3.6, 3.8, 4.1, 6.5,
  6.8, 7.2, …]` shows the plain `k(k+2)`-ish ladder, **no 2I sieve / missing-even-degree
  pattern**. No leak. `[OBSERVED]`
- The relaxed spectrum diverges from the round reference in the **higher modes** (relaxed
  jumps 4.1→6.5; round rises smoothly 3.9→4.6→5.1→…) — the glass keeps a different mode
  structure even after variance is minimised.

### Honest limitations (flagged)

- **Defective null did not give a clean contrast** — both generic and defective glass
  (neither reaches round). So the discriminative content is "low variance ⇏ round"
  (both achieved low variance, neither round), realised as *both-glass* rather than
  *generic-separates-from-defective*. The negative control still does its job (it proves
  small variance is not sufficient) but does not isolate the generic arm as special.
- **Single N (=100), coarse.** No finite-size scaling yet. The spectral-distance metric
  has compressed dynamic range here (the round reference itself match-errs 0.275 to the
  analytic `k(k+2)` at N=100), so 0.280→0.249 should be read qualitatively (a small move),
  not over-quantified. FSS across 10²–10⁴ is owed (Mr A #7).
- **L-BFGS terminated `ABNORMAL`** (numerical-gradient line-search failure near the flat
  minimum) — benign (it marks the local min), but an analytic gradient would let it probe
  deeper.

### What this licenses (and what it does not)

- **Licenses:** the v1.0 instrument failure is fixed (round is a stable fixed point by
  construction); and a clean, pre-registered observation — *edge-length curvature-variance
  minimisation on a FIXED generic triangulation does not reach round S³; it glasses*.
  `[OBSERVED]` Per the brief's own pre-registered response, reaching round calls for the
  **connectivity degrees of freedom** (finer mesh / remeshing / a Perelman-surgery
  analogue), **not** a tuned fix.
- **Does NOT license:** any claim that round S³ is unreachable in principle (this is one
  fixed mesh at one N), nor any chirality/2I claim (Stage 2 territory).

### Handed up to CinC — staging decision (NOT taken here)

Stage 1 did **not** pass (generic arm glasses, no relaxed-round metric produced).
Per debugger discipline ("Stage N after N−1 passes"), I am **not** proceeding to Stage 2
unilaterally. The options, for CinC:
1. **Proceed to Stage 2 on an analytic round S³ / the 2I arm** (the brief explicitly
   permits "an analytic round S³ at this stage, for cleanliness") — the chirality test
   does not strictly need the generic arm to have reached round.
2. **Address the substrate first** — finer mesh / remeshing moves (edge flips) so the
   generic arm *can* reach round, then re-run Stage 1. This is the brief's prescribed
   response to the glass null, but it is a materially larger build.
3. **Accept the glass null as the Stage 1 result** for the generic arm and record it.

Recommendation: (1) — run Stage 2 (chirality on the 2I arm, handedness out of the drive)
on an analytic round S³, since that is the next *independent* question and does not
depend on the generic arm reaching round; carry the Stage 1 glass-null as an honest,
separately-reported finding. **Stage 3 (back-reaction) remains gated.**

---

## STAGE 2 (chirality, 2I arm) — INSTRUMENT STATUS (curl gate)

*Per the Stage-2 brief §6: instrument-first. Substrate + curl-MAGNITUDE validated; the
SIGNED curl (the chirality crux) is an open instrument blocker — held before the 2I
sieve/dynamics.*

- **Substrate (600-cell):** built as the convex hull of the 120 unit icosians; verified
  V=120, E=720, F=1200, C=600, χ=0; fully 2I-symmetric (5 tets/edge, edge=1/φ). Scalar
  `k(k+2)` ladder clean (mult 4, ~9). `[VALIDATED]`
- **Curl MAGNITUDE (Nédélec curl-curl `K=∫curl(W)·curl(W)`):** annihilates gradients
  exactly (`max|K·d₀|=4e-15`), clean kernel dim=119=V−1, and lowest nonzero
  `μ=4.43 (mult 6)`, `μ≈9 (mult ~16, smeared)` — the analytic `|λ|=2,3` with `n(n+2)`
  multiplicities. **MAGNITUDE GATE PASS.** `[VALIDATED]`
- **Curl SIGN (helicity split): OPEN — blocker.** The Whitney helicity form
  `H[e,e']=∫W_e∧dW_e'` (now correctly gradient-annihilating, `max|H·d₀|=1e-16`) gives
  helicity ≈ ±0.01 — i.e. **~0** — on the μ=4 Hopf-sextet modes, instead of ±2. The
  symmetric part is tiny and the non-symmetric solve gives complex eigenvalues. So the
  +2/−2 chirality is **not cleanly recoverable** from this form on the coarse 600-cell.
  This gates Stage 2 §3–§5 (2I sieve, dynamics, amphichiral control) — not run on an
  unvalidated signed curl (W-108). `[HELD]`

**Proposed fixes (CinC's call):** (a) a proper DEC **primal–dual signed curl**
(`⋆₁,⋆₂` Hodge stars + topological `d₁`), the literal `⋆d`, rather than the Whitney
helicity quadratic form; (b) a **2I-symmetric refinement** of the 600-cell for
resolution (the brief allows it) — the coarse 120-vertex mesh may simply not resolve the
helicity sign; (c) a known closed-form signed-curl construction if Fizz has one. The
magnitude result (the `n(n+2)` fingerprint) stands regardless.

---

## STAGE 2 CHIRALITY — PRE-REGISTRATION (committed BEFORE the run; Fizz reroute)

The signed curl is OFF the critical path. curl-curl is sign-blind (its μ=4 eigenspace is
the full sextet +2⊕−2, so a real solver returns helicity-cancelling mixtures — that is
why the helicity form read ~0). The chirality is read instead from the **exact** left/right
2I action: **+2 modes are left-invariant, −2 modes are right-invariant**. Project the
validated μ=4 Hopf sextet onto the left- and right-2I-invariant subspaces (averaging the
exact signed edge-permutation reps over the 120 group elements). The one-sided quotient
S³/2I keeps one half, kills the mirror. No signed operator, no "+2/−2" label needed
(respects Mr A #5 by construction); the 2I action is an exact permutation so dimensions
are clean integers at N=120.

**PRE-REGISTERED PREDICTION:**
`dim left-2I-invariant(μ=4) = 3`, `dim right-2I-invariant(μ=4) = 3`, `intersection = 0`.
→ the sextet splits into two chiral triples; the S³/2I quotient retains one (+2 Hopf),
kills the mirror = the Beltrami sieve, confirmed discretely by exact symmetry.

**NULL / surprise:** any other split (e.g. 6/0, 0/0, nonzero intersection) → STOP and
hunt the cause (group-action bug, wrong eigenspace) before interpretation.

## STAGE 2 CHIRALITY — RESULT: CLEAN SIEVE (prediction confirmed) `[OBSERVED-as-consistency]`

*Run of record: `python chirality_sieve.py` on the round 600-cell, seed-free (exact).*

> **The μ=4 Hopf sextet splits exactly into two disjoint chiral triples by the exact
> left/right 2I action. The one-sided S³/2I quotient keeps one (+2 Hopf, left-invariant)
> and kills the mirror (−2, right-invariant). The Beltrami sieve is confirmed discretely
> — by exact symmetry, with no signed curl.**

| quantity (μ=4 sextet) | predicted | observed |
|---|---|---|
| dim **left**-2I-invariant | 3 | **rank 3** (proj. trace 2.968) |
| dim **right**-2I-invariant | 3 | **rank 3** (proj. trace 2.968) |
| intersection (both-invariant) | 0 | **0** (principal cos ≈ 0.003) |

- The **amphichiral control is intrinsic:** on the *full* round S³ both triples are present
  (sextet = 3⊕3, mirror-symmetric → no chirality preference). The chirality appears *only*
  under the one-sided quotient — exactly the contrast that proves the selection is
  **geometric**, not put in by any drive (there is no drive here).
- **Orientation convention (Mr A #5):** "left = +2" is the convention; an orientation flip
  relabels left↔right. The convention-independent content is *one chirality survives the
  one-sided quotient, its mirror does not*.
- Trace 2.968 (not 3.000) is the coarse-mesh smearing of the μ=4 block (μ∈[4.27,4.53]);
  the **rank is the clean integer 3** and the intersection is exactly 0 because the 2I
  action is an *exact* edge permutation, not a discretised operator.

**Status:** CONJECTURED-from-rep-theory (the +2 Hopf survives S³/2I, −2 projected out),
**confirmed discretely** here. Consistency on imported 2I — **not** an emergence, not a
discovery. Licenses: *"on S³/2I the geometry selects one Beltrami chirality; the handedness
is in the one-sided manifold, not in any pump."* Stage 3 (back-reaction) remains gated.

---

## STAGE 1 FSS — PRE-REGISTRATION (committed BEFORE the ladder; hard gate §2/§6)

**The hinge.** Stage 1 gave a glass null at N=100 (d_frust ≈ 0.249). Is the generic-arm
glass FUNDAMENTAL (round S³ not a generic attractor — headline survives, unhedged) or a
COARSE-MESH ARTEFACT (round reachable at fine N — headline falls, honest negative on
generality)? The FSS ladder decides.

**Primary observable:** `d_frust(N) = spectral_distance(relaxed frustrated endpoint,
same-mesh round reference)` — same-mesh round-init, NOT the analytic k(k+2) ladder
(topology-dominated, blunt). **Secondary (cross-check, not a vote):** the variance gap
`CoV_frust(N) − CoV_round(N)`.

**Pre-registered fit (frozen):**  `d_frust(N) = d_∞ + c · N^(−p)`.
- **d_∞ consistent with 0** → **ARTEFACT.** Headline falls; paper carries the negative on
  generality (round *is* reachable on a generic substrate at fine N) + stands on the
  methodology and the exact-symmetry chirality (Stage 2).
- **d_∞ a positive constant, bounded away from 0** → **FUNDAMENTAL glass.** Headline
  survives, unhedged: round S³ is not a generic attractor.

No tuning of the fit, the fit form, or the frustration recipe to favour either outcome.

**Ladder:** N = 100, 200, 500, 1000 (floor for a verdict); extend to 3000, 10⁴ as
tractable. Report where compute bites.

**Hygiene:** (i) consistent frustration recipe across N — generic_s3(N) (4D-hull of N
random S³ points, coords discarded) with frustrated ℓ_e=1 — same *character*, only
resolution differs; (ii) same-mesh round reference at each N; (iii) converge the
relaxation at each N (L-BFGS to a reported local min; confirm by identical-endpoint
re-run at the largest N). Efficiency: an O(N) local finite-difference gradient
(perturbing edge f touches only the tets containing f), verified against the numerical
gradient before use — this is what makes N=10³ tractable; it does not change the scheme.

**Controls:** Control #3 at each N (round-init stays round — instrument validity
precondition; if it fails at some N, STOP and report); reference sanity `d_round(N)`
(round-init's residual distance to analytic round) must fall toward 0 as N grows, else
the reference is suspect — report before reading d_frust.

**Status:** the verdict (fundamental vs artefact) is OBSERVED once the ladder is in; the
extrapolation is the evidence, fit residuals reported.

---

## STAGE 1 FSS — RESULTS

*(filled after this pre-registration is committed.)*

## STAGE 1 FSS — RESULTS (3 of 4 rungs; N=1000 owed) `[OBSERVED, verdict leaning ARTEFACT]`

*Run of record: `python fss_stage1.py`, seed=20260618, O(N) local-FD gradient
(verified vs numerical to 3e-6). Each rung's relaxation converged to a reported local
min; Control #3 holds at every rung.*

> **The generic-arm glass is a COARSE-MESH ARTEFACT, not a fundamental trap (strong, on
> 3 converged rungs; N=1000 confirmation owed).** d_frust collapses toward 0 as the mesh
> resolves round S³ — so the frustrated metric *does* reach round at fine resolution, and
> the headline *"round S³ is not a generic attractor"* does NOT survive. Caught before
> print — exactly what the ladder was for.

| N | d_frust | d_round (ref, →0?) | Control #3 |
|---|---|---|---|
| 100 | 0.250 | 0.275 | round ✓ |
| 200 | 0.164 | 0.272 | round ✓ |
| 500 | **0.037** | 0.251 | round ✓ |
| 1000 | *owed* | — | — |

- **Trend:** `d_frust ∝ N^(−1.20)` (log-log), dropping ×1.5 then ×4.4 — a steep power-law
  decline toward 0; power-law extrapolation gives `d_frust(1000) ≈ 0.018`. At N=500 the
  relaxed frustrated metric is essentially *at* the same-mesh round (d_frust=0.037).
- **Pre-registered fit `d_frust=d_∞+c N^(−p)`:** with 3 points (3 params) the fit is
  exactly determined / degenerate, pushing `d_∞` to ≤0 — consistent with **ARTEFACT
  (d_∞≈0)** but not yet a committed value. The 4th rung (N=1000) is needed for the
  non-degenerate fit and the formally-committed verdict.
- **Reference sanity:** `d_round(N)` falls slowly (0.275→0.251) — the same-mesh round
  reference improves with N, as required (not suspect).
- **Variance cross-check:** the CoV gap is *not* informative here (223 at N=500 — K̄
  crosses 0); the spectral distance is the verdict, exactly per Mr A #2.

**N=1000 is owed and BLOCKED by infrastructure, not physics:** the N=1000 relaxation
(~8 min) is auto-backgrounded by the harness, and every background job has been torn down
at the session boundary before completing (5 attempts). To land it: a stable ~10-min
window, or an O(1)-per-partial vectorised/analytic gradient (≈30 s total) — flagged for
CinC. The 3-rung trend already makes the verdict clear; N=1000 firms the fit.

**Verdict (pending N=1000):** **ARTEFACT** — round S³ *is* reachable on a generic
substrate at fine N; the N=100 glass was resolution. The paper carries the honest negative
on generality and stands on the methodology + the exact-symmetry chirality (Stage 2). `[OBSERVED]`
