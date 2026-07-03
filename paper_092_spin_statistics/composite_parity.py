#!/usr/bin/env python3
"""
Paper 92 spin-statistics lemma -- the hinge and the verdict (step 4-5).

The DERIVED backbone (twoI_character_table.py) established, construction-free:
  * 2I splits into INTEGER {1,3,3',4,5} and SPINORIAL {2,2',4',6} irreps by
    how the central element -1 acts (+1 vs -1).
  * a composite's statistics bit = (-1)^(number of spin-1/2 constituents).
  * P-e-P is spin-statistics-consistent  <=>  the bound electron core is
    integer-spin (parity-even).   [the iff]

This script tests the NAMED forcing mechanism from the pre-reg (neutrino
ejection, recovered from Paper 92 s3.1) against the HARD RULE
(permission != requirement), and reports the verdict.

No empirical fitting -- this is exact integer (Z/2) parity bookkeeping.
Run:  python3 composite_parity.py
"""

# spin-1/2 fermion = parity-odd (-1). boson/integer-spin = parity-even (+1).
ODD, EVEN = -1, +1


def parity(*factors):
    """Statistics bit of a composite = product of constituent parities."""
    out = 1
    for f in factors:
        out *= f
    return out


def label(bit):
    return "fermion (half-integer)" if bit == ODD else "boson (integer)"


# Intrinsic parities (observed / standard-model facts, not assumptions
# about the bound core):
PROTON = ODD       # spin 1/2
ELECTRON = ODD     # spin 1/2
NEUTRINO = ODD     # spin 1/2

OBSERVED = {
    "neutron":  ODD,   # spin 1/2
    "deuteron": EVEN,  # spin 1
    "N-14":     EVEN,  # spin 1 (the historical nitrogen anomaly target)
}


def composites(core_parity):
    """P-e-P statistics bit for each system given the bound-core parity."""
    return {
        # bare neutron = proton + bound core
        "neutron":  parity(PROTON, core_parity),
        # deuteron = proton + shared bound core + proton
        "deuteron": parity(PROTON, core_parity, PROTON),
        # N-14 = 7 protons + 7 "neutrons" = 14 protons + 7 bound cores
        "N-14":     parity(*([PROTON] * 14 + [core_parity] * 7)),
    }


def report(title, core_parity):
    print(f"\n{title}")
    print(f"  bound-core parity assumed: "
          f"{'EVEN (integer-spin)' if core_parity == EVEN else 'ODD (spin-1/2)'}")
    res = composites(core_parity)
    allok = True
    for name in ("neutron", "deuteron", "N-14"):
        ok = res[name] == OBSERVED[name]
        allok &= ok
        print(f"    {name:<9}: predicted {label(res[name]):<22} "
              f"observed {label(OBSERVED[name]):<22} "
              f"{'OK' if ok else 'CONTRADICTION'}")
    print(f"  -> {'ALL THREE CONSISTENT' if allok else 'MODEL INCONSISTENT'}")
    return allok


def main():
    print("=" * 70)
    print("PAPER 92 SPIN-STATISTICS LEMMA -- HINGE TEST AND VERDICT")
    print("=" * 70)

    # ---- Reading A: no mechanism. Bound core IS the electron (a spinor). ----
    a_ok = report(
        "Reading A  --  core retains the electron's spinor identity "
        "(no ejection):",
        ELECTRON)

    # ---- Reading B: neutrino-ejection mechanism (Paper 92 s3.1). ----
    print("\n" + "-" * 70)
    print("THE NAMED MECHANISM -- neutrino ejection (Paper 92 s3.1, "
          "independently motivated):")
    print("  Formation ejects the electron's extended field as the neutrino")
    print("  (energetic reason: core contracts Bohr -> Compton). The neutrino")
    print("  is spin-1/2 (observed). Fermion-number conservation:")
    print("        electron(ODD) = neutrino(ODD)  x  remnant")
    print("    =>  remnant parity = ODD / ODD = EVEN.")
    remnant = ELECTRON * NEUTRINO          # (-1)*(-1) = +1  (exact Z/2)
    assert remnant == EVEN
    print(f"  Computed remnant (bound core) parity = "
          f"{'EVEN' if remnant == EVEN else 'ODD'}  "
          f"[electron x neutrino = {remnant:+d}]")

    b_ok = report(
        "Reading B  --  core = parity-even remnant after neutrino ejection:",
        remnant)

    # ---- conservation audit at formation -------------------------------
    print("\n" + "-" * 70)
    print("CONSERVATION AUDIT (total fermion parity at formation/decay):")
    # deuteron formation: p + p + e  ->  deuteron + neutrino
    lhs = parity(PROTON, PROTON, ELECTRON)
    rhs = parity(OBSERVED["deuteron"], NEUTRINO)
    print(f"  deuteron formation  p+p+e -> d + nu : "
          f"LHS {lhs:+d}  RHS {rhs:+d}  "
          f"{'CONSERVED' if lhs == rhs else 'VIOLATED'}")
    assert lhs == rhs
    # neutron beta decay: n -> p + e + antinu
    lhs = OBSERVED["neutron"]
    rhs = parity(PROTON, ELECTRON, NEUTRINO)
    print(f"  neutron beta decay  n -> p+e+nu     : "
          f"LHS {lhs:+d}  RHS {rhs:+d}  "
          f"{'CONSERVED' if lhs == rhs else 'VIOLATED'}")
    assert lhs == rhs
    print("  Both consistent: the parity-even bound core is the conservation-")
    print("  consistent reading; one spin-1/2 neutrino per core carries the")
    print("  electron's spinor factor out.")

    # ---- hard rule: permission != requirement --------------------------
    print("\n" + "=" * 70)
    print("HARD RULE CHECK (permission != requirement):")
    print("  * 'single-valued spatial standing-wave mode' only PERMITS even")
    print("    parity (a bound electron keeps spin-1/2 regardless) -> not enough.")
    print("  * neutrino ejection REMOVES one real spin-1/2 factor from the")
    print("    bound count; by conservation the remnant parity is FORCED even.")
    print("    This is forcing, not permission -- IF the model ejects exactly")
    print("    one spin-1/2 neutrino per core.")
    print("  * That ejection is in Paper 92 for INDEPENDENT energetic reasons")
    print("    (Bohr->Compton contraction) and matches beta-decay (one anti-")
    print("    neutrino per neutron): requirement-type, not tuned for N-14.")

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    assert not a_ok and b_ok
    print("  Reading A (no mechanism): all three CONTRADICTED -> the naive")
    print("    P-e-P count reproduces the 1932 catastrophe.")
    print("  Reading B (neutrino ejection): all three CONSISTENT.")
    print()
    print("  OUTCOME 1 (conditional) -- the icosahedral P-e-P neutron SURVIVES")
    print("  the spin-statistics gate, by ONE mechanism that cures the neutron")
    print("  (1/2), the deuteron (1), and nitrogen-14 (boson) together.")
    print()
    print("  Status: DERIVED that consistency <=> core parity-even (the iff,")
    print("  construction-free). DERIVED that neutrino ejection forces parity-")
    print("  even, CONDITIONAL on the model's independent commitment of exactly")
    print("  one spin-1/2 ejecta per core. That single load-bearing assumption")
    print("  is named for Mr Adversary: if formation could eject the extended")
    print("  field as spin-0 (or as != 1 neutrino), the forcing weakens to")
    print("  permission and the honest verdict drops to OUTCOME 3.")
    print("=" * 70)


if __name__ == "__main__":
    main()
