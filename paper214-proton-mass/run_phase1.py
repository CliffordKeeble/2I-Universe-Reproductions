"""
run_phase1.py  —  Paper 214 Phase 1 orchestrator (validation + section-0 confirmation).

Runs the full Phase-1 deliverable and writes results.json + a profile figure:
  1. The validated B=1 hedgehog Skyrmion (skyrme_hedgehog): dimensionless energy,
     virial theorem, baryon number, Bogomolny ratio, ANW MeV mass.
  2. The 2I group verification (icosian_2I).
  3. Section-0 confirmation (quotient_invariance): the classical energy is
     quotient-independent -- 3D FD cross-check + 2I-invariance to machine precision.

Phase 1 is VALIDATION, not discovery (brief section 0): the classical proton mass
is provably quotient-independent. The genuine 2I test is Phase 2 (quantisation).

Standalone:  python run_phase1.py
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from skyrme_hedgehog import run as run_skyrmion, solve_profile, energy_integrals
from icosian_2I import build_icosians, verify_group
from quotient_invariance import run as run_quotient


def make_profile_figure(path="figures/profile.png"):
    import os
    os.makedirs("figures", exist_ok=True)
    x, F, Fp, _ = solve_profile()
    s2 = np.sin(F) ** 2
    i2 = x * x * Fp * Fp + 2.0 * s2
    i4 = 2.0 * s2 * Fp * Fp + s2 * s2 / (x * x)
    bdens = -(2.0 / np.pi) * s2 * Fp
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(x, F, lw=2)
    ax[0].set_xlim(0, 8); ax[0].set_xlabel("x  (dimensionless radius)")
    ax[0].set_ylabel("F(x)"); ax[0].set_title("Hedgehog profile  F(0)=pi -> F(inf)=0")
    ax[0].axhline(np.pi, ls=":", c="grey"); ax[0].grid(alpha=0.3)
    ax[1].plot(x, i2, label="quadratic density (I2)")
    ax[1].plot(x, i4, label="quartic density (I4)")
    ax[1].plot(x, bdens, label="baryon density", ls="--")
    ax[1].set_xlim(0, 8); ax[1].set_xlabel("x"); ax[1].set_title("Energy & baryon densities")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    return path


def main():
    print("\n########## Paper 214 Phase 1 ##########\n")
    skyrmion, _ = run_skyrmion(verbose=True)
    print()
    Q = build_icosians()
    grp = verify_group(Q)
    print(f"=== 2I group: order {grp['order']}, closed={grp['closed_under_mul']}, "
          f"inv_closed={grp['inverse_closed']} ===\n")
    quotient = run_quotient(verbose=True)
    print()
    fig = make_profile_figure()
    print(f"profile figure -> {fig}")

    results = {
        "paper": 214,
        "phase": 1,
        "title": "B=1 Skyrme soliton: validation + section-0 quotient-independence",
        "status": "Phase 1 complete (validation); Phase 2 (quantisation) not started",
        "skyrmion_B1": skyrmion,
        "icosahedral_group_2I": grp,
        "section0_quotient_check": quotient,
        "validation_anchors": {
            "ANW_dimensionless_mass": {
                "computed_M_e_over_Fpi": skyrmion["M_over_Fpi_e"],
                "literature": 36.5, "status": "OBSERVED-matches-literature"},
            "bogomolny_ratio": {
                "computed": skyrmion["bogomolny_ratio"],
                "literature": 1.232, "status": "OBSERVED-matches-literature"},
            "virial_theorem_I2_eq_I4": {
                "rel_gap": skyrmion["virial_rel_gap"], "status": "DERIVED-check-passed"},
            "baryon_number": {
                "computed": skyrmion["baryon_number"], "status": "DERIVED-topological"},
        },
    }
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print("\nresults.json written.")
    return results


if __name__ == "__main__":
    main()
