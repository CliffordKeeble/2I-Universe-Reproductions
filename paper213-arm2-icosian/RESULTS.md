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

*(filled by `run_arm2.py` and `compare_arms.py` after the pre-registration commit.)*
