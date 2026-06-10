# Paper 208 — Candidate B

The E8/icosahedron ghost. **Candidate B** tests whether the golden content discarded
in the E8 = 2I ⊕ ⋆2I split (Dechant 2016, Sec 6) encodes either programme target
(μ = 6π⁵, α⁻¹ = π+π²+4π³).

**Result: predicted NULL, confirmed.** Every golden-content quantity is rational; the
targets are transcendental. The null is the publishable artefact — see
[`candidate_b_results.md`](candidate_b_results.md) for the full table and verdict.

Pre-registration of record: Candidate B freeze, `Session_Handoff_E8_Hunt_20260610.md` §5
(commit-regardless-of-outcome was a frozen output requirement; the outcome is null;
the commit stands).

## Files

- `candidate_b.py` — builds the 120 icosians (2I = H4 roots), verifies count, unit
  norm, and full multiplicative closure over 14,400 products. The grounding warrant.
- `candidate_b_quantities.py` — both conventions, all pairwise inner products as exact
  `a + τb` over ℤ[τ], reports Q1/Q1×/Q2/Q3 under all three frozen normalisations.
- `candidate_b_results.md` — results table, verdict, and the two honesty notes.
