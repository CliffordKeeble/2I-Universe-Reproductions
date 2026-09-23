# Paper 213 Stage-1 — FSS ladder (arm 1, .NET `Icosian.Regge`)

*instrument: Icosian@1790cc2 (DIRTY) · seed 20260618 · 2026-07-01T07:56:24Z*  
*Observable: `d_frust` (normalised-spectrum distance, relaxed vs round reference). Committed fit `d_frust = d_∞ + c·N^(−p)`; `d_∞ → 0 ⇒ artefact`, `d_∞ > 0 ⇒ fundamental glass`. Criterion frozen (arm 2) — not re-tuned. Seed 20260618.*

| N | d_frust | Control #3 | Ef | ‖∇E‖ | conv | exit | min tet vol |
|---|---|---|---|---|---|---|---|
| 100 | 0.0745 | 0.0278 | 59.32 | 2.43E+002 | False | linesearch-fail | 1.43E-009 |
| 200 | 0.1494 | 0.0075 | 53.08 | 5.13E+001 | False | maxIters | 1.10E-007 |
| 500 | 0.0620 | 0.0021 | 74.53 | 4.32E+001 | False | linesearch-fail | 2.47E-009 |

**Convergence gate:** 0/3 rungs reached `‖∇E‖ < tol` and enter the fit.


**Verdict:** NO VERDICT — 0/3 rungs converged (‖∇E‖<tol). The generic-glass relaxation reaches constrained boundary minima (a tet driven to the embeddability boundary; see minVol), so ‖∇E‖<tol is structurally unreachable and no fit is earned. See STATUS.

A near-zero **min tet vol** means the relaxation drove a tet onto the embeddability boundary — a constraint-active minimum where `‖∇E‖` is nonzero by construction, so the gradient-norm criterion is structurally unreachable. Control #3 (round reference, same instrument) stays small ⇒ the reader is sound; this is a relaxation/observable issue.

---

**23 September 2026 — superseding note, appended; nothing above is edited** `[Fizz 🌀; Paper 213 v0.2 §3]`. *"The generic-glass relaxation reaches constrained boundary minima … so ‖∇E‖<tol is structurally unreachable"* is **not supported.** CinC's ruling of 26 August found a fully feasible descent direction at the N = 100 endpoint (`59.32 → 54.67`, no cell made non-embeddable). An exact-derivative barrier–Newton minimiser also releases the near-flat cells and descends far lower. **The rungs are unconverged wall-strike values, not constraint-active minima.** The observable awaits a certified minimum.
