# Findings — Cold-Script Check-In, R-Bench, W-Bench, and the 𝒞 Frontier

**Paper 226 · Golden Construction · one sitting: W4 → A → B → C → W1–W3 → E**
**6 July 2026 · exact over ℚ(√5, i), numeric cross-check at both real places**

**Disciplines held:** every pre-registration / derivation note committed *before* its verification;
cold artifacts committed verbatim; printed scripts run as written; deviations reported, not patched.
Per **CinC Ruling 1 (badge discipline)** this findings file reports **outcome-vs-EXPECT only** — no
DERIVED/STRUCTURAL language; promotion happens after merge and Mr Adversary, not here. Per **Ruling 2**
W4 was derived cold from the explicit matrices before any diff against the adjudication memo.

Run order followed the addendum: **W4 (keystone) → A → B → C → W1–W3 → E.**

---

## Part A — cold-script check-in

Cold artifacts (Fable-5-max `verify_c_from_sigma.py` + `verification_run.log`) committed verbatim
under `cold_boot/`; the two paper `.md` drafts stay uncommitted in Downloads per the standing rule,
their Appendix-A claims vendored with citations into the reproduction.

| Item | EXPECT | Outcome |
|---|---|---|
| A-1 independent reproduction (`verify_c_from_sigma_repro.py`) | every claim reproduces exactly | **53/53 pass**, exact + both real places, via an independent mechanism (regular representation of 𝕂: N=det, places=eigenvalues, σ=Z-conjugation) |
| A-2 coverage diff vs cold script | coverage superset; zero disagreements | **superset**: all 49 cold claim-families reproduced + 4 extras (N=det, place-coords=eigenvalues, coords(Φ²), δ²-faithfulness); **zero assertion disagreements** |
| A-3 count resolution (50 vs 52) | resolves to bookkeeping | **resolves to 49**: script `check()`=49 and log `[PASS]`=49 agree exactly; draft "50" = documentation off-by-one (×3); relay "52" = line count (log is 51 physical lines). Bookkeeping, not a substantive gap |
| A-4 both green + merge | both scripts green | cold **49/49**, repro **53/53**; merged |

**Independence caveat (recorded, not hidden):** the reproduction was written after the cold script
was read (the package arrived mid-check-in). Independence is at the derivation level (different
representation, each claim re-derived from the draft's stated maths), **not** script-blind.

**Two cold boots (provenance flag):** Fable-5-max produced the scripted construction; Opus-4.8-max
produced a prose second witness (`cold_boot_charge_conjugation_C_equals_sigma.md`). They **converge**
on C = σ but **diverge** on the gauge structure — Fable gauges the non-compact split torus over K,
Opus a compact U(1) over K(i). This fork is benched in W4 and Part E, not adjudicated here.

## Part B — B1.6-4b, the −k re-run (asterisk-remover)

B-i derivation note committed first (no pass/fail language): τ is antilinear, τ(e^{ikx}) = e^{−ikx},
so a τ-borne candidate carries +k → −k; linear and σ-borne candidates fix i and stay at +k.

| Item | EXPECT | Outcome |
|---|---|---|
| B-i evaluation point | τ → −k (antilinear); lin, σ → +k | **derived**: momentum signs read off from applying each map to e^{ikx} (lin +1, σ +1, τ −1) |
| B-ii frontier re-run | dims 0/0/0 over K; τ-borne dim 2 over K(i) | **reproduces Run #1.6 exactly**: σ-coupling 0/0/0; i-coupling τ-borne dim 2. The −k was produced by antilinearity, not chosen |

The corrected-in-flight asterisk on B1.6-4 is resolved: the evaluation point follows from the map's
definition independent of outcome.

## Part C — reconciliation R-bench (one spine)

Registered against `fizz_reconciliation_seam_paper_226.md` §4 (`PRE-REGISTRATION_rbench.md`).

| Item | EXPECT | Outcome |
|---|---|---|
| R1 embedding δ→√5σ₃ | unital K-algebra embedding; j→Z; e± idempotents | **pass** |
| R2 transported gauge action | e^{θj}=diag(e^θ,e^{−θ}); reproduces A4 weights | **pass** |
| R3 entrywise σ = Galois σ | equal on embedded 𝕂 | **pass** |
| R4 null-ray signs | σ(Γ₊)=−Γ₊, σ(Γ₋)=+Γ₋; C₀(Γ₊)=+Γ₊, C₀(Γ₋)=−Γ₋ | **pass** |
| R5 even-sector control | C=σ∘(a-flip) maps q→−q, no coordinate factor | **pass** (D_{−A}(σa)=σ(D_A a) exactly) |
| R6 the dichotomy | fails on odd sector w/o reflection; ok on even | **pass on both halves**: geometric image is a reflected operator (B1.5-2 reconfirmed), even-projected norm N(a) swap-invariant |

The internal construction is the even-sector shadow of the geometric one; the reflection is the toll
the odd part pays. j = Z is the recognition.

## Part D — W-bench

W4 derived cold first (`W4_derivation.md`), then verified and diffed.

| Item | EXPECT | Outcome |
|---|---|---|
| W4-ii keystone | cold derivation agrees with memo (bivector flips, generator fixed) | **8/8**: σ(Γ_seed)=−Γ_adj, σ(Γ_adj)=−Γ_seed; B=Γ_seed·Γ_adj=√5Z; C₀(√5Z)=−√5Z; C₀(iZ)=+iZ; **matches memo — fork stays closed on independent evidence** |
| W4b charge sectors | C₀ keeps compact charge, τ flips it | **pass**; the brief's eigenspace gloss clarified — τ *fixes* the real J=iZ eigenvectors, the "exchange" is the charge-level τ(J)=−J |
| W1 J=iZ properties | J²=−I; e^{θJ}=diag(e^{iθ},e^{−iθ}); τ(J)=−J; σ(J)=+J | **pass** |
| W2 two real forms | exp(ℝZ) split; exp(ℝJ) compact; both in exp(ℂZ) | **pass** |
| W3 compact reverser | τ maps q→−q with no matrix part | **pass**: entrywise τ with B=I; B=I lies in B1.6-4's τ-borne dim-2 family |

W4 is the keystone of the "one torus, two real forms" absorption of the Fable/Opus fork, and it held
cold — reproduced from ground truth, not from Fizz's path.

## Part E — the 𝒞 frontier (eight cells)

On the i-dressed (1,3) clique over K(i) (Run #1 A6/A8). Every reported representative **verified**
(numeric, both places) to satisfy 𝒞·g(Γ^μ)·𝒞⁻¹ = ε·(Γ^μ)ᵀ.

| Cell (g, ε) | dim | (𝒞∘g)² invariant |
|---|---|---|
| 1, + | 1 | 𝒞²=−I (scalar) |
| 1, − | 1 | 𝒞²=−I (scalar) |
| **τ, +** | **1** | (𝒞τ)²=diag(5,1,5,1), NON-scalar (√5-graded 5:1), place-definite |
| **τ, −** | **1** | (𝒞τ)²=−diag(5,1,5,1), NON-scalar, place-definite |
| σ, + | **0** | — (empty) |
| σ, − | **0** | — (empty) |
| στ, + | **0** | — (empty) |
| στ, − | **0** | — (empty) |

| Item | EXPECT | Outcome |
|---|---|---|
| E-1 τ-type nonempty | nonempty in ≥1 sign | **pass** (dim 1 in both signs) |
| E-1 σ-type / στ-type | OPEN (no expectation) | **empty (dim 0)** — the 4d form of "σ needs the compensator"; C_phys carries the reflection factor absent from the pure matrix frontier |
| E-2 squares | report value + sign at both places | linear 𝒞²=−I (scalar); τ (𝒞τ)²=±diag(5,1,5,1), √5-graded, place-definite (4d echo of A8's √5-normalisation; scale is convention, 5:1 ratio and sign invariant) |
| E-3 interlock | τ-cell ↔ B1.6-4 τ-family; σ-cell ↔ C_phys | τ-cell is the 4d generalisation of B1.6-4's τ-borne family; σ-cell empty = the 4d "σ needs the compensator", matching R6 |

**The Letter-Assignment data:** on this clique the standard C-matrix relation is satisfiable for the
**τ-type** and **not** for σ / στ. This is data for the Letter-Assignment theorem (Fizz's warm work,
next), not an adjudication of it. It coheres with the two-faced picture: σ is charge conjugation at
the field-automorphism level (needs the reflection compensator, R6 / empty σ-cell), while τ carries
the matrix-level reversal (nonempty τ-cell, the compact-coupling reverser of W3).

---

## Flags for CinC

1. **The Fable/Opus fork was absorbed, not papered over.** W4 (cold) and Part E give it independent
   footing: the "one torus, two real forms" cancellation holds from ground truth, and the eight-cell
   frontier separates τ (matrix-bearing) from σ (compensator-needing) cleanly.
2. **The τ-square is non-scalar** (√5-graded diag(5,1,5,1)) — the Letter-Assignment "C²=±1" question
   needs the A8-style renormalisation discussion on this √5-normalised clique before a clean
   letter-square can be quoted. Reported at the strength shown, no more.
3. **Two brief wordings clarified, not failed:** W4b's τ-eigenspace gloss (τ fixes the real J=iZ
   eigenvectors; the exchange is the charge-level τ(J)=−J) and A-3's "50 vs 52" (true count 49).
4. **Not on disk:** the adjudication memo (W4/W1–W3 cite §2/§5) and the two paper drafts — the
   addendum and seam note served as the pre-registrations of record; drop the memo in if you want it
   in the repo lineage.

**Routes (per the briefs):** A → the cold construction's April 9 standing + v0.2 provenance; B → the
asterisk resolution in v0.2's bench appendix; C → v0.2's reconciliation section; D/E → the
Letter-Assignment theorem input.

🐕☕⬡
