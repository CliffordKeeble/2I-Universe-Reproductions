# Findings — Paper 118, "The 11 Barrier" (Part 1: DERIVED-tier reproduction)

**Brief 118-CODE-1, Part 1** (Amendment A1: Python, reproductions repo, pytest).
`python -m pytest test_eleven_barrier.py` → **12/12 PASS**. `python eleven_barrier.py`
writes `m_2I.csv` (l = 0..120) and `fibre_cycles.csv` (l = 0..59).

## Verdict (one line)

Every DERIVED-tier claim of Paper 118 reproduces from first principles: the spectral
multiplicity `m_2I(l)` computed three independent ways agrees exactly for `0 ≤ l ≤ 120`;
the barrier `m(10) = 0` and the gap `λ₁ = 12·14 = 168` fall out with no tuning.

## Machinery

For finite Γ < SU(2), the multiplicity of the trivial rep of Γ in the (l+1)-dim irrep
of SU(2) is `m_Γ(l) = (1/|Γ|) Σ_C |C| χ_l(θ_C)`, with `χ_l(θ) = sin((l+1)θ)/sin(θ) = U_l(cos θ)`
(Chebyshev U), `χ_l(0) = l+1`, `χ_l(π) = (−1)^l (l+1)`.

2I class data (Σ sizes = 120, asserted before anything else): θ = 0 (1), π (1),
π/5, 2π/5, 3π/5, 4π/5 (12 each — Vtx), π/3, 2π/3 (20 each — Face), π/2 (30 — Edge).

**Three routes, agreement is the verification gate:**

- **(a) character sum** — mpmath at 60 dps; integrality asserted `|m − round(m)| < 1e-30`. **[OBSERVED→DERIVED]**
- **(b) Molien** — exact integer power series `Σ m(l) t^l = (1 + t^30)/((1 − t^12)(1 − t^20))`,
  i.e. `m(l) = p(l) + p(l−30)`, `p(n) = #{(a,b)≥0 : 12a + 20b = n}` (Klein degrees 12, 20, 30). **[DERIVED]**
- **(c) exact Q(√5)** — Chebyshev-U recurrence in the ring Q(√5); each fibre sum is Galois-stable,
  so its contribution is an exact rational, and the l = 12 vertex values land on {φ, ψ} exactly. **[DERIVED]**

Routes (a), (b), (c) agree for all `0 ≤ l ≤ 120` (`test_three_routes_agree_and_integral`, AT6).

## Acceptance tests — all green

| Test | Claim | Result |
|---|---|---|
| **AT1** | m(0)=1; m(l)=0 for 1≤l≤11; m(12)=1; λ₁ = 168 | ✅ gap at l=12, λ = 12·14 = 168 |
| **AT2** | m(10)=0 with (Gen,Vtx,Face,Edge) = (11/60, +2/5, −1/3, −1/4); 22+48−40−30 = 0 (÷120) | ✅ exact 0 |
| **AT3** | all four vertex characters = 1 at l=10 | ✅ [1,1,1,1] in Q(√5) |
| **AT4** | l=12 = (13/60, +1/5, +1/3, +1/4), total 1; vertex chars {ψ,φ,φ,ψ}, Σ = 2 in Q(√5) | ✅ multiset {φ,φ,ψ,ψ}, Σ = 2 exact |
| **AT5** | Table 5: m at l = 0,10,20,30,40,50 = 1,0,1,1,1,1 with stated face/edge slots+values | ✅ all six rows |
| **AT6** | Molien identity, coefficient-by-coefficient to l=120 | ✅ (a)=(b)=(c) ∀ l≤120 |
| **AT7** | fibre cycles: Vtx period 10, Face period 6, Edge period 4; Gen(l)=(l+1)/60 even, 0 odd | ✅ l = 0..59 exact |
| **AT8** | among l≡0 (mod 10), l≤120, only l=10 has m=0; m(70)=1 via 71+24−20−15 = 60 (÷60) | ✅ zeros = [10]; m(70)=1 |
| **AT9** | l=12 is smallest l>0 with l≡0 (mod 4), l≡0 (mod 6), Vtx(l mod 10)>0 | ✅ smallest = 12 |

## The barrier, exactly

At **l = 10** the vertex fibre returns perfectly (all four χ = 1, Vtx = +2/5) while the
face and edge fibres sit at their maximally destructive slots — 10 mod 6 = 4 → −1/3,
10 mod 4 = 2 → −1/4 — and the generic baseline 11/60 is absorbed:

```
m(10) = 22/120 + 48/120 − 40/120 − 30/120 = 0/120 = 0        [DERIVED]
```

At **l = 12** the vertex fibre carries golden values {ψ, φ, φ, ψ} summing to 2(φ+ψ) = 2
(not the flat 4 of l = 10), face and edge return constructive, and all four fibres reach
consensus: `13/60 + 12/60 + 20/60 + 15/60 = 60/60 = 1`. The gap is `λ₁ = 168`.

**Uniqueness (AT8):** l = 10 is the *only* vertex return in l ≤ 120 that vanishes. At l = 70
the destructive face/edge combination recurs identically, but the generic baseline has grown
to 71/60 and compensates: `71/60 + 24/60 − 20/60 − 15/60 = 1`.

## Optional §12 item — deferred, not run

The brief flags the §12 "Hopf-to-generic power-spectrum crossover" (paper states l ≈ 1112)
as **optional, OBSERVED-tier**. It is **not reproduced here**: the defining statistic — the
"power spectrum of Z_2I(t)" and what "crossover" means on it — is not specified in the supplied
Paper 118 text (§12 asserts l ≈ 1112 with no formula) or in the brief. Reproducing it would
require guessing the definition, and reverse-engineering a definition until it lands on 1112
is exactly the move the two-route discipline forbids. Flagged for CinC: supply the exact
statistic and it computes in one pass. **[not run — definition unavailable]**

## Files

- `eleven_barrier.py` — class data, three routes, exact fibre decomposition, table writer.
- `test_eleven_barrier.py` — AT1–AT9 + substrate (12 tests).
- `m_2I.csv` — l, λ, m (all three routes), fibre fractions, l = 0..120.
- `fibre_cycles.csv` — fibre contributions and slots, l = 0..59.

## Reproduce

```
python -m pytest test_eleven_barrier.py -v   # 12/12
python eleven_barrier.py                      # summary + CSVs
```

## Scope

This is Part 1 only — the machinery's validation instance (Gate G2). Part 2 (the
solvable-control null across all finite SU(2) subgroups) is **held** pending the W-107
diff clearance and Cliff's explicit go; no Part 2 artefact (including its pre-registration)
is written yet.
