# Target 1 (Paper 213 Stage 3) — Component Inventory & Re-scope Input

**For CinC, via Cliff. From Mr Code.** 29 June 2026.
**Purpose:** the kickoff brief tells Mr Code to *reuse* the "paper213 Stage-1 metric
relaxer." It is not present to reuse in the target instrument. This document is the
accurate component state so Phase A can be re-issued against real repo reality rather
than the "already-landed" assumption. **No build code was written; this is the pause.**

---

## 0. Baseline (verified this session)

- **.NET `Icosian` repo**, HEAD `c171b43`, branch `main`, clean tree.
- `dotnet build` → 0 warnings, 0 errors. `dotnet test` → **Core 84 · TwoI.Core 11 ·
  Lattice 21 · SixHundredCell 14 = 130 passed, 0 failed.** (SixHundredCell eigensolve
  ~2m20s, the one slow test.) The instrument the brief points at is green.

## 1. The two-repo reality

The brief's component list mixes pieces that live in **two different repos and two
different languages**:

| Brief says "reuse…" | Actually lives in | Language | State |
|---|---|---|---|
| `Icosian.Core` (the 120 icosians) | `C:\Dev\Icosian` | .NET | present, validated |
| `Icosian.SixHundredCell` (read-out) | `C:\Dev\Icosian` | .NET | present, validated — **but fixed-metric, curl-only (see §3)** |
| "paper213 Stage-1 **metric relaxer**" | `2I-Universe-Reproductions\paper213-relaxation` | **Python** | **source lost from VC; bytecode only (see §4)** |
| Beltrami / +2 Hopf field (Stage 2) | both | .NET (curl op) + Python (curl_instrument) | present |

The brief's hard gates — **BCL-only, the `SixHundredCell` reader, hand-rolled
Cholesky+Jacobi, the xUnit suite, `check_sixhundredcell.py` as *witness*** — all name
the **.NET** instrument. So the intended build home is .NET/`Icosian`. The metric
relaxer it tells us to reuse there was **never ported to .NET**.

## 2. What the .NET port actually covers

The "213 port" into `Icosian.SixHundredCell` ported the **Stage-2 spectral / curl**
side on the **fixed round metric** only:

- `SixHundredCell.cs` — the 600-cell as a **combinatorial** complex (V120/E720/F1200/
  C600 + exact left/right 2I vertex permutations). Edges are index pairs with one
  uniform inner product `cos36°`. **No per-edge length DOFs, no deficit angles.** Its
  own doc-comment: *"the space layer … it carries no field."* [VERIFIED by read]
- `WhitneyCurlOperator.cs`, `BeltramiSieve.cs`, `SpectralReader.cs`, `Linalg.cs` —
  signed curl-curl operator + 2I sieve + hand-rolled generalized eigensolver.
  `SpectralReader.Read(SixHundredCell cell)` takes **only the combinatorial cell** and
  reads the **curl/Beltrami** spectrum (the +2 Hopf sextet, 3⊕3 sieve). [VERIFIED]

**Two consequences for Target 1:**
1. The read-out is **fixed-metric**: it has no way to accept a relaxed (back-reacted)
   edge-length vector. Reading a *moved* metric requires re-parameterising the operators
   on per-edge lengths.
2. The read-out currently exposes the **curl** spectrum, not the **scalar `n(n+2)`
   Laplace** ladder that P1's verdict (the m-ladder split) needs. (A linear-FEM scalar
   Laplacian *exists in the Python `spectral_s3.py`*; in .NET it must be confirmed
   present/parameterised on edge lengths.)

## 3. What Target 1 actually needs that does **not** exist in .NET

Everything that makes the metric *move* — the whole reason Stage 3 is "the first arm
not pre-forced" (pre-reg §0):

- **Regge geometry over edge-length DOFs** on the 600-cell: tet embedding from 6 edge
  lengths, dihedral angles, deficit angles `δ_e = 2π − Σθ`, dual volumes, Regge action.
  *(Python `geometry.py` / `fast_geometry.py` — not ported.)*
- **The metric relaxer**: fixed-volume curvature-variance minimiser
  `E(ℓ)=Σ w_e (K_e−K̄)²/Σ w_e`, `K_e = δ_e/a_e`, `a_e = V_e^dual/ℓ_e`. *(Python
  `relax_geometry.py` + `relax_fast.py` — not ported. This is the "do not resurrect the
  neckpinch Euler scheme" replacement the brief refers to.)*
- **The scalar-Laplace read-out parameterised on relaxed edge lengths** *(Python
  `spectral_s3.py` — the verdict instrument; .NET has the curl analogue, not this.)*
- **The generic (2I-free) substrate + WL automorphism gate** for the P2 arm *(Python
  `generic_s3` in `geometry.py` + `stage0_gate.py` — not ported.)*
- **NEW regardless (already flagged in the brief, §6/B1):** the matter-coupled source
  `∂S_matter/∂ℓ` and the co-relaxation loop. This was always "new"; it now sits on top
  of an *also-unported* metric layer.

**Net:** Phase A is **build-from-scratch of the Regge metric layer in .NET**, not a
reuse. That roughly doubles the "new" surface relative to the brief's effort model.

## 4. The Stage-1 relaxer is recovered (design intact), but its source was lost

`paper213-relaxation/` is **entirely untracked** in `2I-Universe-Reproductions`
(`git status` → `?? paper213-relaxation/`), and `.gitignore` excludes `*.pyc` /
`__pycache__`. The `.py` sources have been **deleted from the working tree** — only
`__pycache__/*.cpython-313.pyc` survive. So the relaxer was **never committed** and the
text source is gone.

**Recovered this session** (interpreter Python 3.13.14 matches the `cpython-313`
bytecode): module docstrings live intact in the bytecode, so the **complete design
narrative, full API (signatures), exact formulas, the seed `20260618`, and convergence
criteria** of all nine modules are recovered, plus a faithful `dis` disassembly of each.
Saved durably here:
- `RECOVERED_module_map.txt` — docstrings + signatures of all 9 modules.
- `RECOVERED_disassembly/*.dis.txt` — full bytecode disassembly per module.

The nine modules (the actual Stage-1 + Stage-2 Python instrument):
`geometry.py` (Regge machinery + generic_s3), `fast_geometry.py` (vectorised curvature),
`ricci_flow.py` (normalized discrete Ricci flow, CoV metric, Control #3 round_metric),
`relax_geometry.py` (fixed-volume variance min, L-BFGS-B over log ℓ — the neckpinch
fix), `relax_fast.py` (O(N) local-gradient FSS version), `spectral_s3.py` (FEM
Laplace verdict vs `k(k+2)`), `sixhundred_cell.py` (2I substrate), `stage0_gate.py`
(WL automorphism gate), `curl_instrument.py` (signed curl — the .NET `WhitneyCurlOperator`
is its port).

This is enough to **port faithfully** without resurrecting anything wrong; a byte-exact
`.py` reconstruction (from the saved `dis`) is available as a follow-up if wanted.

## 5. Version-control hygiene flag (independent of Target 1)

The 213-relaxation reproduction — the source behind Paper 213's published Stage-1 glass
null and FSS — exists **only as untracked, gitignored bytecode**. One `__pycache__`
clear and it is gone for good. **Recommend:** commit the recovered artifacts (and, if a
`.py` reconstruction is commissioned, the sources) into `2I-Universe-Reproductions` so
the published result stays reproducible. *(Mr Code did not commit the untracked tree
unilaterally — adding a whole tree to VC in a second repo is an ask-first action.)*

## 6. Re-scope options (for CinC's call)

1. **Build Target 1 in .NET, Phase A = build the Regge layer + relaxer (recommended).**
   Honours every hard gate. Sequence: (A0) edge-length-DOF Regge geometry on the
   600-cell; (A1/C-A) curvature-variance relaxer, validate round-stays-round; (A3) field
   stress-energy T; then Phase B coupling. Port the *method* from the recovered Python
   (formulas in §4) — not the neckpinch Euler scheme. Larger than briefed, but lands the
   metric arm inside the one validated, gated instrument.
2. **Build in Python (`paper213-relaxation`), extending the existing relaxer.** Reuses
   the lineage directly (after a `.py` reconstruction), but contradicts the brief's
   .NET/BCL-only/`SixHundredCell` instructions and uses the ungated Python read-out.
3. **Re-scope the brief's Phase A** to explicitly commission the .NET Regge-layer port
   as its own gated deliverable *before* the coupled solver — i.e. option 1, but with the
   effort acknowledged up front and pre-registered as a build, not a reuse.

**Mr Code's recommendation:** option 1 / 3 (they are the same path, honestly scoped).
The pre-registration (P1–P3, controls, kill conditions) is unaffected and still binds;
only the *engineering route* and the *effort estimate* change.

🐕☕⬡

*Component inventory & re-scope input, v0.1 — the pause before the build.*
