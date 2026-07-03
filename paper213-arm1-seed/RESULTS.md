# Paper 213 — ARM 1: Seed/dilation coupling, import nothing

**Non-circular falsification of v1.0's premise.** Grow a coupling graph by the
Paper-191 seed (the √5 dilation), import *nothing* about 2I, and measure
(a) the spectral dimension `d_s` of the grown graph and (b) whether any
S³/2I ⟨12,20,30⟩ signature appears in either Laplacian.

Status labels per claim: **DERIVED / OBSERVED / CONJECTURED**.

---

## PART 1 — PRE-REGISTRATION (Pattern 75)

*Written and committed before the first production run. The commit timestamp is
the audit trail. Nothing below this section's criteria was tuned after seeing
production numbers.*

### 1.1 The seed and its normalisation (with a documented wrinkle)

Seed (Paper 191 §6): `Γ_seed = [[0,1],[√5,0]]`, `det = −√5`.
Normalisation: `R = Γ_seed / 5^¼`.

Symbolically verified (`seed.py`, DERIVED):

| fact | value | brief's wording | status |
|---|---|---|---|
| (a) `det R` | `−1` | "det −1" | ✅ holds |
| (b) `R²` | `I` | "R² = I, order-2" | ✅ holds |
| (c) `Rᵀ R` | `diag(√5, 1/√5) ≠ I` | "order-2 **reflection**" (implies O(2)) | ⚠️ **R is NOT orthogonal** |

**The wrinkle.** The 5^¼ scaling fixes the determinant and the order, **not**
orthogonality. No scalar multiple of `Γ_seed` can be orthogonal because its
columns have unequal norms (√5 and 1). So `R` is a *linear* reflection (det −1
involution), not a *Euclidean/O(2)* one — the brief's hinge sentence ("the
holonomy must be orthogonal, so the seed must be normalised") is imprecise.

**Why the verdict is unaffected (in fact more robust).** The pre-registered
"no signature" expectation rests on `⟨I,R⟩ ≅ ℤ/2` being abelian with no
5-torsion. That holds because `R² = I` — orthogonal or not. And the connection
Laplacian `L_c` is still **real-symmetric** (we set the `(j,i)` block to
`U[i,j]ᵀ`), so its spectrum is real regardless of (c). We build exactly the
operator the brief defines and carry (c) forward as a flag for the cold
Adversary (W-108). One *measurable* consequence of (c): `L_c` need not be
positive-semidefinite, so negative eigenvalues, if they appear, are the
fingerprint of the non-orthogonal seed — reported, not suppressed.

### 1.2 What is fixed before running

- **Growth:** GR-1 self-reference rewrite (`growth.py`). The brief sketches GR-1
  but leaves the *cycle-formation* mechanism open; that mechanism is made
  explicit here (triangle closure gated by the closure predicate `C`). This is
  **Mr Code's instantiation** and it moves `d_s` — which the brief pre-registers
  as OPEN — but **not** the signature verdict, which rests on the ℤ/2 holonomy
  and is instantiation-independent.
- **Coupling:** every edge holonomy ∈ `{I, R}` (Arm 1 = abelian by construction).
- **Read-outs:** `L_s = D − A` (scalar) and `L_c` (connection), both built and
  scored at every `N`.
- **Pre-registered signature constants** (`spectral.py`, frozen):
  `M_CLUSTERS = 10`, `TAU = 0.05`, `REL_GAP = 0.15`, `ZERO_REL = 1e-6`,
  `JITTER = 1e-9`.
- **Targets:** k(k+2) for k=1..10 = `{3,8,15,24,35,48,63,80,99,120}`;
  ⟨12,20,30⟩ semigroup gaps `{2,4,6,8,10,14,16,18,22,26,28,34,38,46,58}`;
  Paper-108 spectral gap `λ₁ = 168` (not expected on a growing graph, whose
  gap → 0 — reported as informative null).
- **Controls** (`controls.py`): non-golden seeds (√2, √3, generic a=7),
  scrambled `{I,R}`, bare-vs-connection, `d_max` sweep `{3,4,5,6,8}`.
- **Rule-space sweep:** ρ ∈ {GR-1, all-R} × C ∈ {flat-triangles, flat-shortest,
  none} × d_max, logged to `logs/rule_trials.jsonl`.

### 1.3 Pre-registered verdict criteria

The signature score band of the **controls** (scrambled + non-golden) is the
null distribution for the signature. The verdict is read against it:

- **CLEAN NULL (expected).** Signal (golden seed) signature scores on `L_c` and
  `L_s` are statistically indistinguishable from the control band — i.e. the
  golden seed produces no ⟨12,20,30⟩ structure that scrambled/non-golden growth
  does not also produce by chance. `d_s` may take any value (OPEN). **This
  closes the v1.0 hinge honestly.**
- **SIGNATURE PRESENT → STOP + ESCALATE (W-108).** If the golden seed's `kk2`
  *or* `semigroup_hit` on either Laplacian sits materially **above** the control
  band (and `λ₁` approaches 168 in the L_c normalisation), something carried 2I
  in. Do **not** bank an emergent icosahedron from inside the build chain
  (W-107). Hand to a cold Adversary to find the door (graph automorphism? a C/ρ
  variant that built symmetry in?).

**"Materially above"** is pinned, pre-run, as: golden score exceeds
`mean(control) + 3·sd(control)` on the same N **and** exceeds the best control
by ≥ 0.25 absolute. Anything less is within-null.

---

## PART 2 — RESULTS

*(populated by `run_arm1.py`; written after Part 1 was committed)*

<!-- RESULTS_PLACEHOLDER -->
