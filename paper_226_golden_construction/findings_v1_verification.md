# V1 — Verification of the Arc Report against the benches

**Paper 226 · Golden Construction · verification pass (NOT banked)**
**13 July 2026 · `bench_v1.py` · `results_v1.csv` · `run_v1.log`**
**Object verified:** `arc_report_golden_dirac_construction.md` v1.0 (Fizz; in Downloads, uncommitted).
**Brief of record:** `mr_code_brief_V1_verify_arc_report.md` (Fizz, 13 Jul 2026).
**Against:** committed `findings_f1_f4/f5/f7/f8/f9.md` + `scout_report_galois_dirac.md` (Downloads).

**Mode:** mechanical check — does the report say what the benches say? Where they differ, **the findings win.**
No physics adjudication, no framing improvements. **Nothing here is banked.**

---

## Verdict

**The Arc Report is faithful to the benches on the physics — every F1–F9 number checks out.** Four drifts,
all in the *narrative/attribution* layer, none in a bench result:

1. **The survival claim overstates novelty.** The report says the σ-linearity theorem is "NOT FOUND in the
   literature (Q3, Q4)"; **Scout marked Q3 and Q4 FOUND** (textbook corollary). *(§1-D2 — the important one.)*
2. **The retraction ledger count is wrong: seven, not six.** *(§1-D3 — the item the brief asked about.)*
3. **Scout attributions are imprecise / not verifiable** from the provided Scout report (the icosian
   `(−1,−1)_ℚ(√5)` division-algebra detail and the Wilson `𝔽₉`/CPT/fiat claims; "Scout Q7" mislabelled). *(§1-D4.)*
4. **`det = 125` is not in the report** (it cites rank 4 only); the findings carry it. *(§1-D1 — minor.)*

**The one new computation confirms the frontier and resolves the co-descent flag:** σ and στ cells are empty
(dim 0); **the two obstructions are INDEPENDENT** (against Fizz's suspicion) and may be counted as two
confirmations; M6 and M7 both stand.

---

## §1 — line-check (findings win)

### Matches (re-verified; report = findings)
All of the following in the report agree with the committed findings, and the cheap ones recompute exactly
(`bench_v1.py` §1 spot-recompute: **det = −m²s² − √5(k²−E²), E² = k² + m²s²/√5, (Γ_s+Γ_a)² = 0, Z² = +I** all PASS):

- **F1–F4:** the gamma/Z matrices; det and E²; A₀ = diag(√5,1) parity **diag(+,−)**. ✅
- **F5:** aligner family **A = [[0,b],[c,0]]** antidiagonal, 2-dim. ✅ **Taussky–Zassenhaus citation confirmed
  correct:** *On the similarity transformation between a matrix and its transpose*, **Pacific J. Math. 9 (1959),
  893–896.** ✅
- **F7:** reality-variety **b = c**; **K₂** real/nondegenerate/EOM-preserving at s=1 (θ=π/2, τ); **symplectic
  rank 4** (recomputed, matches); Fizz's candidate [[0,√5],[1,0]] **fails** (b≠c); **halo refuted** (K₂ σ-even,
  rational). ✅
- **F8:** σ fixes i ⟹ ψ̄→**e^{+iα}ψ̄** ⟹ no U(1); conservation fails **even on the τ-locus with an in-field
  σ-even mass**; **m = 5^{1/4}μ, min poly x⁴−5, m ∉ K**. ✅
- **F9:** (Γ_s+Γ_a)²=0; **M** and Mᵀdiag(√5,−√5)M = diag(1,−1), M non-diagonal; Cl ≅ split M₂(K); Spin(1,1)
  abelian; centraliser = scalars; roots of unity in K(i) = {±1,±i} order 4; |2I|=120; τ-adjoint ψ̄→**e^{−iα}ψ̄**,
  conserved U(1) both mass parities (τ(m)=m), metric diag(√5,1) positive-definite, no ghosts; 2I bipartition
  **60/60/120**, E₈ marks **Σ=30**, σ swaps 2↔2′/3↔3′, **2cos72°=(√5−1)/2, 2cos144°=−(√5+1)/2**; (½,0)/(0,½)
  same SU(2) irrep, pseudoreal. ✅

### Drifts (report is wrong / imprecise where it differs)

- **D1 (minor). `det = 125` is not in the report.** The brief §1 says "the report cites rank 4 as the invariant
  and 125 as basis-dependent — confirm." The report cites **rank 4** (lines 55, 120) but **never states 125**.
  The findings state both (rank 4 invariant, **det 125 = 5³** basis-dependent) and 125 recomputes. So the report
  is *consistent* (rank 4) but **omits** the 125 detail — not an error, an omission. No fix needed unless the
  report wants the basis-dependence caveat in it.

- **D2 (important). The survival claim overstates novelty.** Report §5 (line 135) and exec summary (line 25):
  the σ-linearity theorem is *"NOT FOUND in the literature (Q3, Q4)."* **Scout's report marks Q3 = FOUND and
  Q4 = FOUND**, both *as corollaries of textbook material*: Q3 is the **first-kind/second-kind involution
  dichotomy** (KMRT, *Book of Involutions*) — Scout writes he "did not find it stated in the specific words …
  but **it is a corollary of textbook material and any algebraist would call it such**"; Q4 is the standard
  **Majorana-existence obstruction** (De Andrade–Toppan, hep-th/9904134). So the *exact sentence* is not in the
  literature verbatim, but the **result is textbook**. Calling Q3/Q4 "NOT FOUND" and the theorem "one nobody has
  written down" (§9) **drifts from Scout's own FOUND verdict.** *Findings/Scout win: the theorem's honest status
  is "a first/second-kind corollary, not stated in these words," not "novel/absent."*

- **D3 (the flagged one). The retraction ledger says six of Fizz's; it lists SEVEN.** §4 headline: "Six of the
  ten are mine." The "claimed by = Fizz" (solely) rows are **R1, R2, R4, R5, R6, R7, R8 = seven** (R3 is shared
  Mr Code+Tor+Fizz; R9 Paper 226; R10 the programme). §4's own breakdown is also internally inconsistent: "four
  flagged as wants/leans + R1, R6, R8 genuine errors" = **seven**, not six. *The brief asked directly ("if there
  are seven, say so"): there are seven.*

- **D4 (attribution / not fully verifiable). Scout attributions are imprecise.** In `scout_report_galois_dirac.md`
  as provided:
  - The report (line 65, §5) attributes to **"Scout (second sweep), Q7"** the icosian ring being *a maximal order
    in `(−1,−1)_ℚ(√5)`, totally definite, a division algebra*, and *"the split M₂(K) is a different quaternion
    algebra … 2I cannot live in it."* **Scout's Q7 is about φ in quantum theory** (2cos π/5), not the icosian
    ring. Scout's icosian content is in **Q1** ("a quaternion algebra over ℚ(√5), **definite at both real
    places**") and **Q6**, stated **more weakly** — **no** verbatim "(−1,−1)_ℚ(√5)", "maximal order", or "division
    algebra", and **no** "2I cannot live in the split M₂(K)" comparison. The **Wilson** claims (L/R **by fiat**,
    conjugation from an **𝔽₉ automorphism negating i** that **does not lift**, **no CPT**) are **not in this Scout
    file** at all (Scout's Q2 discusses Wilson only as chirality-not-C, and notes he read Wilson only partially).
  - **These may live in a second Scout sweep not on disk** (the brief lists "both sweeps"; only one file is
    present). Per the brief ("report faithfully; do not re-derive his literature") I do **not** adjudicate the
    physics — but I flag them **"not checkable from the provided Scout report."** Do not wave them through.
  - The **Q-number labels are loose**: R7's kill is credited to "Scout Q7", but the σ-preserves-positivity /
    total-positivity content is Scout's **Q5/Q6** (and Surprise 4), not Q7.

---

## §2 — the 226 §5.2 C-frontier (the new computation)

### §2.1 — eight cell dimensions (DERIVED)
`𝒞·g(Γ^μ)·𝒞⁻¹ = ε(Γ^μ)ᵀ`, g ∈ {1, τ, σ, στ}, ε ∈ {±1}, dim of the 𝒞-space:

| g \ ε | +1 | −1 |
|---|---|---|
| **1** | 1 | 1 |
| **τ** | 1 | 1 |
| **σ** | **0** | **0** |
| **στ** | **0** | **0** |

**The σ and στ cells are empty (dim 0); the 1 and τ cells are non-empty (dim 1).** The report's frontier claim
is **confirmed**.

### §2.2 — M6 reading: STANDS (DERIVED)
`𝒞 Γ 𝒞⁻¹ = ε Γᵀ` is exactly the standard **transpose charge-conjugation-matrix condition** (the g-dressing
supplies the antilinear/field part, C = 𝒞∘g). An empty σ-cell ⟹ **no C-matrix for the σ-dressing** ⟹ "σ cannot
be C." **Fizz has not misread its role; M6 stands.**

### §2.3 — the co-descent (W-108) verdict: **INDEPENDENT** (PRIORITY, DERIVED)
On the **i-free** gammas, dress(1) = dress(τ) and dress(σ) = dress(στ). The empty cells are exactly **{σ, στ}**,
which tracks the **√5-SWAP** {σ, στ} — **not** the **i-MOTION** {τ, στ}. The decisive tell:

> **στ MOVES i yet its cell is EMPTY (like σ); τ MOVES i yet its cell is NON-EMPTY (like 1).**

So the frontier emptiness is driven by the **√5 place-swap of the gamma pair** (σ acting on the √5 in the
entries), **not** by σ fixing i. The **σ-linearity theorem** is driven by σ **fixing i** (its action on the
field). **Different aspects of σ ⟹ the two results are INDEPENDENT.** They **may** be reported as **two
confirmations.**

*This is against Fizz's registered suspicion (she suspected co-descent and would rather lose a confirmation).
Reported flat: they do not co-descend. The one thing she asked to guard against — double-counting — is not
triggered, because the double-count she feared isn't there.*

### §2.4 — M7 (dissolved tension): STANDS (DERIVED)
K₂ **as a conjugation matrix** fails the σ-frontier (no ε works — this *is* the "named tension"). K₂ **as a
pairing** gives the metric K₂Γ_s = diag(√5,1) used in ψ̄=(σψ)ᵀK₂. A **sesquilinear pairing form** and a
**conjugation intertwiner** are **different structures**; K₂ passing reality (as a pairing) and failing the
frontier (as a conjugation) is **no contradiction. M7 stands** — the tension was a category confusion.

---

## §3 — what the brief asked for, collected

- **Every drift:** D1 (det 125 omitted, minor), **D2 (survival claim overstates novelty — Scout Q3/Q4 = FOUND),
  D3 (ledger is seven of Fizz's, not six)**, D4 (Scout attributions imprecise / from an absent second sweep).
  The physics numbers all match.
- **Eight frontier cell dimensions:** 1/±1 = 1, τ/±1 = 1, **σ/±1 = 0, στ/±1 = 0.**
- **Co-descent verdict (priority):** **INDEPENDENT** — the frontier emptiness tracks the √5-swap, not σ-fixes-i;
  στ moves i yet is empty. The frontier and the σ-linearity theorem are two results, not one.
- **Not checkable from the benches:** the Scout icosian-`(−1,−1)_ℚ(√5)`/division-algebra detail and the Wilson
  `𝔽₉`/CPT/fiat claims (not in the provided Scout file; possibly a second sweep). External literature (Varlamov,
  Berg–DeWitt-Morette et al., Conway–Sloane, Baez, De Andrade–Toppan, KMRT) is reported faithfully but not
  benched. The narrative/interpretive lines (§0, §6 M1–M5, §7, §9) are prose, not bench claims — not verified,
  not disputed.
- **Fizz-right-that-she-wasn't:** the retraction ledger credits **seven** retractions to Fizz solely
  (R1,R2,R4,R5,R6,R7,R8), while §4 says six. **It is seven.**

## Reproduce
```
cd paper_226_golden_construction
python bench_v1.py   # -> results_v1.csv; frontier 8 cells, co-descent INDEPENDENT, M6/M7 stand, §1 spot-recompute
```
Environment: `environment.txt` (Python 3.13.14, sympy 1.14.0). Exact symbolic. **Nothing here is banked** — the
report goes to Tor; Mr A gets the findings cold, never the narrative.

🐕☕⬡
