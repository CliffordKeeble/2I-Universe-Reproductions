# 2I-Universe-Reproductions

Reproduction artifacts for the numerical claims in the **2I Universe Programme**
Zenodo papers. One directory per claim: each holds a pre-registration, the
enumeration/derivation scripts, the null distribution or survivor set, and a findings
write-up. Everything here is meant to be re-run from a clean checkout with nothing more
than Python and numpy.

This repo is the public, code-side companion to the (private) **2I-Universe-Papers**
repo. The papers state the claims; this repo lets you check them.

## The discipline

Every OBSERVED numerical match in a programme paper — "signal X lands on target Y at rate
R" — is a claim about **separation from chance**, not pattern-spotting. It is not a real
claim until its null has been tested and reported (the programme calls this Pattern 75).
This repo is where those nulls live. If a paper quotes a sub-percent coincidence, you
should be able to land in the matching directory here and reproduce the number, the search
space it was found in, and how many other things that same search space would have hit.

If you are a sceptical reader (the kind who wants to know whether a "0.01% match" is one
hit out of a handful or one hit out of thousands of tries): that is exactly the question
each `findings.md` is written to answer. Read the `PRE-REGISTRATION.md` first — it fixes
the search space *before* the result — then run the script.

## Layout

```
2I-Universe-Reproductions/
├── README.md
├── LICENSE
├── paper_101_massratio_null/     # Paper 101 §5.3: the 2I mode l=42 (λ=1848)
│   ├── PRE-REGISTRATION.md        #   minus icosahedral 12 lands on m_p/m_e = 1836,
│   ├── spectrum.py                #   uniquely, at the 0.01% window.
│   ├── massratio_null.py
│   ├── findings.md
│   └── environment.txt
├── paper_092_spin_statistics/    # Paper 92 gate: does the 2I P-e-P standing
│   ├── PRE-REGISTRATION.md        #   wave give the right spin & statistics?
│   ├── twoI_character_table.py    #   DERIVED: consistent iff bound core is
│   ├── composite_parity.py        #   integer-spin; neutrino ejection forces it
│   ├── findings.md                #   -> OUTCOME 1 (conditional). Cures neutron,
│   ├── twoI_character_table.csv   #   deuteron, and nitrogen-14 together.
│   └── environment.txt
├── paper_092_hodge_dual/         # Paper 92 follow-up: is the bridge node the
│   ├── PRE-REGISTRATION.md        #   Hodge dual of the proton, p* = *(proton)?
│   ├── hodge_dual.py              #   Corpus: *(proton) = electron. LANDING 3
│   ├── findings.md                #   (clash bites): magnetic moment not cured
│   └── environment.txt            #   (~1836x), form-degree != spin grading.
└── paper_092_odd_sector/         # Papers 92+174: the spinorial spectrum 174
    ├── PRE-REGISTRATION.md        #   skipped -- "where the actors live."
    ├── spinorial_spectrum.py      #   Spinorial Molien (E8 exponents), Galois
    ├── beat_exchange.py           #   2<->2'. LANDING B: proton=spinor-2 is
    ├── findings.md                #   permitted not forced. Spin-0 beat is the
    └── environment.txt            #   DUAL (2x2), not Galois; d!=nn from Fermi.
```

Each directory follows the programme's `paper_NNN_shortname` convention (matching the
`2I-Universe-Papers` repo), so a paper and its reproduction share the same number prefix.

## Running

```
cd paper_101_massratio_null
python spectrum.py        # reproduces the §4.3 2I spectrum (DERIVED)
python massratio_null.py  # reproduces the §5.3 null (writes matches.csv)
```

## Links

- 2I Universe Programme on Zenodo: https://zenodo.org/communities/ (programme community)
- Paper 101 — *Spectral Geometry of S³/2I*: concept DOI 10.5281/zenodo.18817872

## Licence

See `LICENSE`.
