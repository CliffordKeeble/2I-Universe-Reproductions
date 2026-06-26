# Paper 63 retirement — independent verification of the three G estimates

**Author:** Mr Code (independent compute), June 2026
**Brief:** `Brief_MrCode_paper_63_calc_retire_v0_1.md` (Tor ⚛️, carried by Cliff)
**Paper retired:** *G and the Hierarchy Dilemma* v2.2 (Zenodo concept — Bootstrap/2I series)
**Status of this record:** OBSERVED / physics order-of-magnitude estimates — **NOT** 2I derivations.

---

## Why this exists

Paper 63 argued that the gravitational constant *G* is a **local statistical
quantity** — the local amplitude of a collective baryonic standing wave on S³ —
and read four published anomalies (precision ceiling, >5σ inter-lab spread, the
Anderson 5.9-yr oscillation, BBN consistency) as evidence for it.

The retirement rests on three order-of-magnitude estimates, all from standard
galactic/cosmological inputs, that **self-falsify the local thesis** and leave a
global (Machian/Sciama) picture standing. This file is the independent check of
those three estimates, run on inputs sourced by Mr Code rather than lifted from
the brief (W-107 decorrelation: the independent resource must pull its own
inputs, not re-run the architect's arithmetic on the architect's numbers).

Reproduce with `python verify_G.py`.

## Independently sourced inputs (Mr Code, web, 2026-06)

| Input | Value used | Independent source | Brief's value |
|---|---|---|---|
| Oort total local density ρ₀ | 0.100 ± 0.010 M⊙/pc³ | Holmberg & Flynn / Hipparcos (Oort limit) | ~0.1 ✓ |
| Baryonic local ρ_b | ~0.09 M⊙/pc³ | known matter 0.108 − DM ~0.01 | ~0.09 ✓ |
| Local dark matter ρ_DM | 0.3 GeV/cm³ = 0.008 M⊙/pc³, **smooth, subdominant** | dynamical, 0.3 ± 0.1 GeV/cm³ | ~0.01 (~10%) ✓ |
| Disk scale height | ~300 pc | thin-disk canonical | ~300 ✓ |
| LLR bound Ġ/G | (4 ± 9)×10⁻¹³ /yr → ~10⁻¹²/yr cap | lunar laser ranging (JPL) | <~10⁻¹²/yr ✓ |

Every independent input corroborated the brief within the order-of-magnitude
tolerance. **The inputs were decorrelated and they agreed — corroboration, not a
chorus.**

---

## Check A — terrestrial envelope (kills lab-variation)

If *G* tracked a wave whose smallest source scale is the galactic disk
(L ~ 300 pc ≈ 9×10¹⁸ m), two labs separated by d ≈ 900 km (Paris ↔ Florence)
sample it at:

- **δG/G = d/L = 9.7×10⁻¹⁴** (≈10⁻¹³; ≈9×10⁻¹³ if L = 10¹⁸ m).
- To produce the observed inter-lab spread ~5×10⁻⁴, the wave would need a source
  scale **L_eff = d/(spread) = 1.8×10⁹ m** — a scale with no physical source.
- Predicted spread is **9.7 orders of magnitude below** the observed spread.

**PASS.** Terrestrial labs are predicted to agree to ~12 figures; the observed
>5σ inter-lab spread is systematics, not the wave. The paper's Consequence 2
(structured lab disagreement) is excluded.

## Check B — galactic orbit (the discriminator)

**Local model** (G tracks local baryon density). Over a ±100 pc vertical bob
through a disk of scale height 300 pc, an exponential profile drops by 28%; the
~8% smooth DM component damps the baryonic swing only slightly:

- **local δG/G ≈ 0.26** (≈0.3; up to ~2 with gas/arms). DM does **not** rescue
  this — it is ~10% and smooth, damping by ~10%, not by orders.

**Global model** (Machian, G⁻¹ ∝ Σ mᵢ/rᵢc²). The swinging local slice's
potential-weighted share of the cosmic sum:

- bob-slice term (M_slice ≈ 4.2×10⁵ M⊙ at r ≈ 100 pc): **2.1×10⁻¹⁰**
- eccentric-orbit MW-depth term [GM_MW/(R_orbit c²)]·ecc = 6.0×10⁻⁷ × 0.06: **3.6×10⁻⁸**
- **global total δG/G ≈ 3.6×10⁻⁸ per orbit.**

**Bounds:** LLR ~10⁻¹²/yr × 2.5×10⁸ yr (galactic year) = **2.5×10⁻⁴ per orbit**.

- Local (0.26) exceeds the LLR bound by **~3 orders → EXCLUDED.**
- Global (3.6×10⁻⁸) sits **~4 orders below** the LLR bound → **CONSISTENT.**

**PASS.** *G is global, not local.* This is the discriminator and the number the
whole verdict leans on; the independently-computed slice fractions sit right
next to the brief's (2.1×10⁻¹⁰ vs ~10⁻¹⁰; 3.6×10⁻⁸ vs ~10⁻⁸). The verdict
survives factors of a few.

## Check C — Mach / Sciama estimate (the surviving result)

- **G ≈ R c² / (N m_p) = 6.99×10⁻¹¹** vs measured 6.674×10⁻¹¹ — **4.7% high**
  (baryon budget; within a factor ~1.05).
- With total matter (× Ω_b/Ω_m ≈ 1/6): **1.16×10⁻¹¹** — **5.7× low.**

**PASS.** The baryon budget is the one that lands, mild support for
baryon-sourcing — **no more.** The Mach coefficient is fuzzy to a factor of a
few; the 4.7% does not wear a badge.

---

## Verdict

| Check | Result | Means |
|---|---|---|
| A — terrestrial envelope | **PASS** | predicted ~10⁻¹³ vs observed ~10⁻⁴, ~9 orders short → lab spread is systematics |
| B — galactic orbit | **PASS** | local model excluded by LLR/BBN at ~3 orders; global ~10⁻⁸ consistent → **G is global** |
| C — Mach estimate | **PASS** | G ≈ Rc²/(N m_p) ≈ 7×10⁻¹¹ → G is a global constant; baryon budget lands within ~2× |

**All three pass on independently sourced inputs.** The local statistical thesis
of Paper 63 is self-falsified three independent ways; the surviving content is
the Mach order-of-magnitude estimate.

## The honesty line

These are physics estimates from standard numbers, **not** 2I derivations. This
record verifies that the local model is excluded and the global/Mach picture
survives. It does **not** verify a derivation of N or R — that is the open
question (see `summary.md`). The estimate still **imports the universe (N, R) to
estimate G**; until N and R are geometric, `G ≈ Rc²/(N m_p)` is the Sciama–Mach
relation with imported cosmic inputs, not a 2I retrodiction.

## Flags for CinC / Cliff

1. **`paper_63_v2_2_findings.md` (the v0.2 record the brief's §2a points to)
   never reached my surface** — not in Downloads, not in the repo. This file is
   Mr Code's *independent verification* record, authored from the paper, the
   brief's calculation spec, and my own compute. It deliberately does **not**
   reproduce CinC's per-section history (I never saw it). If you want CinC's full
   v2.2 record archived alongside this, hand me the file and I'll add it.
2. **This repo has no central catalogue/index.** The Retired status + OPEN
   question live in `summary.md` here (per the per-folder-summary convention that
   replaced the catalogue). The master `paper_index_v0_16.md` is in Downloads /
   the private papers repo — outside this repo — so updating *that* index is a
   papers-repo job, not this commit.
3. **Novelty of the Mach line** (Sciama / Dirac / Brans–Dicke prior art) is
   Mr Scout's, cold, before the estimate becomes any paper. Flagged, not
   certified — as the brief instructs.
