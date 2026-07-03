# Paper 214 — The Proton Mass on S³/2I (Skyrmion Energetics)

**Phase 1 — Foundation / validation.** Status: **COMPLETE.**
Brief: `mr-code-brief-proton-mass-S3-2I.md` (CinC, 22 June 2026).

> One-line summary: a numerical B=1 Skyrme soliton solver, validated against the
> textbook numbers, plus an explicit numerical confirmation of the brief's §0 —
> the classical proton mass is **quotient-independent**. This is validation, not
> discovery: §0 proves the classical energy cannot shift under the 2I quotient,
> and the computation confirms it to machine precision. The genuine 2I test lives
> in the **quantisation** (Phase 2), not the classical mass.

---

## What was built

| File | Purpose |
|---|---|
| `skyrme_hedgehog.py` | B=1 hedgehog soliton: profile ODE (BVP), dimensionless energy, virial check, baryon number, Bogomolny bound, ANW MeV mass. |
| `icosian_2I.py` | The 120 unit icosians (binary icosahedral group 2I) as SU(2) target rotations; verified group of order 120. |
| `quotient_invariance.py` | §0 numerical confirmation: independent 3D FD energy cross-check + 2I-invariance of the energy to machine precision. |
| `run_phase1.py` | Orchestrator: runs all three, writes `results.json` + `figures/profile.png`. |

Every script is standalone (`python <file>.py`) with no cross-directory dependency.

---

## 1. The validated B=1 Skyrmion

**Model (brief §1).** Static energy
`E = ∫d³x [ −(F_π²/16)Tr(LᵢLᵢ) − (1/32e²)Tr([Lᵢ,Lⱼ]²) ]`, `Lᵢ = U†∂ᵢU`,
hedgehog `U = exp(i τ·n̂ F(r))`, `F(0)=π, F(∞)=0`. Non-dimensionalising with
`r = (2/F_π e)·x` reduces the energy to

```
E = (π F_π/e) · I[F],
I[F] = ∫₀^∞ [ x²F'² + 2sin²F(1+F'²) + sin⁴F/x² ] dx.
```

The Euler–Lagrange profile equation is solved as a boundary-value problem
(`scipy.integrate.solve_bvp`) on `x ∈ [10⁻³, 30]`.

### Validation anchors (Pattern 75 — validate before trust)

| Quantity | Computed | Reference | Status |
|---|---|---|---|
| Dimensionless mass `M e/F_π = πI` | **36.46** | ANW 1983: **36.5** | **OBSERVED** ✓ |
| Bogomolny ratio `I/(3π)` | **1.2315** | textbook **1.232** | **OBSERVED** ✓ |
| Virial theorem `I2 = I4` | gap `1.3×10⁻⁴` | exact (Derrick) | **DERIVED** ✓ |
| Baryon number `B` | `0.99999993` | `1` (topological) | **DERIVED** ✓ |
| Bogomolny completion identity | `1.5×10⁻¹⁶` | exact | **DERIVED** ✓ |

- `I = 11.6066` (`I2 = 5.8025`, `I4 = 5.8041`). The equality of the quadratic
  (σ-model) and quartic (Skyrme) contributions is the Derrick virial theorem and
  is a convention-independent internal check; it holds to `1.3×10⁻⁴`.
- The **Bogomolny ratio is the gold-standard, convention-independent number**: the
  identity `I = ∫[(xF'+sin²F/x)² + 2sin²F(1+F')²] dx + 3π` gives the bound `I ≥ 3π`,
  and the true profile sits at `1.2315×` the bound — the textbook 1.232.
  The completion identity reproduces `I` to machine precision (`1.5×10⁻¹⁶`),
  confirming both the algebra and the quadrature.

### The MeV value — and what is INPUT vs OUTPUT (brief discipline)

`F_π` and `e` are **MEASURED/ASSUMED** scales, not derived. With the ANW best-fit
inputs `F_π = 129 MeV`, `e = 5.45`:

```
M_classical = (π F_π/e)·I = 863 MeV     [scale F_π/e is INPUT; I = 36.46/π is OUTPUT]
```

This is the classical soliton mass (ANW quote ≈865 MeV); the physical nucleon
(938 MeV) requires the rotational quantisation energy added on top — that is
Phase 2. The model's accuracy on baryon observables is the well-known ~30% level.

**The geometry's genuine output is the dimensionless number `I` (→ `πI = 36.46`).**
The dimensionful mass is `INPUT-scale × OUTPUT-number`. No MeV value is claimed
to come from pure geometry.

---

## 2. The binary icosahedral group 2I

The 120 unit icosians (= 600-cell vertices) are built by the verified construction
(8 + 16 + 96 quaternions) and confirmed to form a group: **order 120, closed under
quaternion multiplication, contains the identity, closed under inverse.** As SU(2)
matrices they are unitary with `det = 1` to machine precision (`2×10⁻¹⁶`).
Status: **STRUCTURAL** (the 2I group is exact, not numerical).

---

## 3. §0 confirmed — the classical mass is quotient-independent

**The argument (brief §0, provable).** The Skyrme energy density depends on `U`
only through the left current `Lᵢ = U†∂ᵢU`, which is *invariant* under a constant
left target rotation `U → qU` (because `(qU)†∂ᵢ(qU) = U†q†q∂ᵢU = U†∂ᵢU`). The
quotient `S³/2I` identifies field values `U ~ qU` for `q ∈ 2I`, so the energy
density is **literally identical on every 2I-orbit** — the functional descends to
the quotient unchanged. The classical energy *cannot* shift. This is provable, not
discovered.

**The numerical confirmation (joint check on code + argument).**

**(A) Independent 3D cross-check.** The full 3D Skyrme energy, computed by finite
differences on the hedgehog field as SU(2) matrices — an entirely separate code
path from the 1D ODE — Richardson-extrapolates (convergence order `p ≈ 1.9`) to:

```
E_tilde(3D) → 36.42  vs  1D target πI = 36.46     (rel. error 0.11%)
|B|(3D)     → 1.0008  vs  target 1
```

(The 3D baryon converges to `−1`; the sign is the orientation convention of the
`ε^{ijk}` contraction relative to the SU(2) embedding handedness — `|B| → 1` is the
physics.) Two independent computations agree → the machinery is sound.

**(B) 2I-invariance to machine precision.** Acting with **every one of the 120**
elements of 2I (`U → qU`) and recomputing the energy:

```
max |ΔE| over all 120 elements = 1.4×10⁻¹⁴     (§0 predicts exactly 0)
max |ΔB| over all 120 elements = 1.1×10⁻¹⁶
```

The energy is identical on every 2I-orbit to machine precision. **§0 confirmed.**
Status: **DERIVED** (the invariance is exact; the residual is numerical noise).

---

## Phase 1 deliverable — checklist against the brief

- [x] Validated B=1 Skyrme soliton solver.
- [x] Dimensionless soliton energy (`πI = 36.46`), matching ANW 36.5 and the
      Bogomolny ratio 1.232.
- [x] Textbook nucleon-scale mass to the model's known accuracy (863 MeV classical,
      ANW inputs labelled INPUT).
- [x] Explicit confirmation that the classical mass is **quotient-independent**
      (§0), both by argument and to machine precision numerically.
- [x] **Not** a 2I-shifted classical mass (excluded by §0).
- [x] **Not** a MeV value from pure geometry (scale `F_π/e` is INPUT).

---

## What Phase 1 does NOT claim, and what comes next

- No 2I effect on the **classical** mass is claimed or found — §0 forbids it, and
  the numerics confirm the prohibition. A shift here would have meant a bug.
- **Phase 2 (not started, per brief):** the genuine 2I test is the **quantisation**
  of the soliton's collective coordinates on the quotient. The quotient reduces the
  target's chiral symmetry to a 2I-related subgroup, changing the **spectrum of
  baryon states** without moving the classical mass — the aether-safe site for a
  real 2I effect. Required Phase-2 check (flagged by Mr A): confirm imposing 2I on
  the zero-mode quantisation does not covertly select a preferred (Type-1) frame.
- The scale-free target is the mass **ratio** `m_p/m_e`, which cancels the input
  scale and meets the existing Paper 101 (1836) claim — the over-determination test.
  This is a Phase-2 deliverable; Phase 1 does not touch it.

**Discipline note (brief).** Built *forwards*: solve, then read the number off. No
2I-combination was hunted to match 938 or 1836. The classical computation was free
to fail (a quotient-dependent classical energy would have been a bug or a
refutation of §0); it did not, exactly as §0 predicts.
