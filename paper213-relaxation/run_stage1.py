"""
run_stage1.py  —  Paper 213 Approach B v1.1, Stage 1 forcing run.

Pipeline (one run per pre-registration; RESULTS.md committed first — hard gate §5):
  Stage 0 gate  : generic complex must be 2I-free (|Aut| bounded) -> PASS or abort.
  Flow          : fixed-volume curvature-variance minimisation (Option 1) on
                    (a) frustrated ell=1  [the TEST]
                    (b) round-init        [control #3: must stay round, by construction]
                    (c) defective mesh    [the NULL: known-can't-reach-round]
  VERDICT       : SPECTRUM (Mr A #2) — spectral distance of each relaxed metric to the
                  round-init reference on its own mesh, and match to the analytic k(k+2)
                  ladder. Variance is the engine; the spectrum judges.
  Anti-circularity (§4): scan the relaxed generic spectrum for 2I-sieve structure
                  (the 15 missing even degrees). Expected NONE.

Prediction badged CONJECTURED. Target = approach the N-dependent round-init floor.

Standalone:  python run_stage1.py [N] [iters]
"""

import sys
import json
import numpy as np
from geometry import generic_s3, defective_s3, verify_sphere, symmetry_probe
from fast_geometry import precompute, fast_energy
from relax_geometry import variance_min
from stage0_gate import gate
from spectral_s3 import (low_spectrum, normalized_spectrum, spectrum_match_error,
                         multiplicity_report)


def round_ell(cx):
    P = cx.hull_points
    return np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])


def spectral_distance(eigsA, eigsB, M=12):
    """RMS relative difference of the first M normalized eigenvalues."""
    a = normalized_spectrum(eigsA); b = normalized_spectrum(eigsB)
    m = min(M, len(a), len(b))
    if m < 3:
        return np.nan
    a, b = a[:m], b[:m]
    return float(np.sqrt(np.mean(((a - b) / b) ** 2)))


def run_one(cx, label, iters, k_eig=40):
    pre = precompute(cx)
    ell0 = np.ones(len(cx.edges))
    ell_r = round_ell(cx)
    E0, cov0, Kbar0, _ = fast_energy(ell0, pre)
    print(f"\n[{label}] frustrated start: E={E0:.4e} CoV={cov0:.3f} Kbar={Kbar0:+.3f}")
    ell_relaxed, hist = variance_min(cx, ell0, pre=pre, iters=iters, verbose=True)
    Ef = hist[-1]["E"]
    print(f"[{label}] frustrated end:   E={Ef:.4e}  (E reduction x{E0/Ef:.1f})  "
          f"note={hist[-1].get('note','')}")

    spec_init, _ = low_spectrum(cx, ell0, k=k_eig)
    spec_relaxed, _ = low_spectrum(cx, ell_relaxed, k=k_eig)
    spec_round, _ = low_spectrum(cx, ell_r, k=k_eig)

    d_init = spectral_distance(spec_init, spec_round)
    d_relaxed = spectral_distance(spec_relaxed, spec_round)
    err_relaxed = spectrum_match_error(spec_relaxed, 4)
    err_round = spectrum_match_error(spec_round, 4)
    print(f"[{label}] SPECTRUM verdict: dist-to-round  init={d_init:.3f} -> "
          f"relaxed={d_relaxed:.3f}  (round-ref match-err={err_round:.3f}, "
          f"relaxed match-err={err_relaxed:.3f})")
    return dict(label=label, E0=float(E0), Ef=float(Ef),
                E_reduction=float(E0 / Ef) if Ef > 0 else np.inf,
                cov0=float(cov0), cov_f=float(hist[-1]["cov"]),
                dist_init=float(d_init), dist_relaxed=float(d_relaxed),
                match_err_round=float(err_round), match_err_relaxed=float(err_relaxed),
                spec_relaxed=normalized_spectrum(spec_relaxed)[:16].tolist(),
                spec_round=normalized_spectrum(spec_round)[:16].tolist(),
                note=hist[-1].get("note", ""))


def control_round_stays_round(cx, iters):
    """Control #3: flow from round-init must NOT blow up (E stays low)."""
    pre = precompute(cx)
    ell_r = round_ell(cx)
    E0, cov0, _, _ = fast_energy(ell_r, pre)
    ell, hist = variance_min(cx, ell_r, pre=pre, iters=iters, verbose=False)
    Ef = hist[-1]["E"]
    print(f"\n[control #3 round-init] E {E0:.4e} -> {Ef:.4e}  "
          f"(stays round: {Ef <= E0 * 1.5})  note={hist[-1].get('note','')}")
    return dict(E0=float(E0), Ef=float(Ef), stays_round=bool(Ef <= E0 * 1.5))


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    print(f"=== Paper 213 Approach B v1.1 — Stage 1 forcing run (N={N}, iters={iters}) ===")
    print("    topology GRANTED (closed S^3); spectrum is the verdict; CONJECTURED.")

    cx = generic_s3(N, seed=20260618)
    chk = verify_sphere(cx, max_link_checks=120)
    g = gate(cx, verbose=True)
    if not g["PASS"]:
        print("STAGE 0 GATE FAILED — generic complex not 2I-free. ABORT Stage 1.")
        json.dump(dict(stage0=g, aborted=True), open(f"results_stage1_v11_N{N}.json", "w"), indent=2)
        return
    print(f"  substrate ok: chi={chk['euler_char']} closed={chk['closed_no_boundary']} "
          f"links_ok={chk['links_ok']}")

    res_generic = run_one(cx, "GENERIC arm", iters)
    res_ctrl = control_round_stays_round(cx, iters)

    cxd = defective_s3(N, seed=20260618)
    res_defect = run_one(cxd, "DEFECTIVE null", iters)

    # anti-circularity: 2I-sieve missing-degree structure on the relaxed generic spectrum
    print("\n[anti-circularity §4] generic relaxed spectrum should show NO 2I sieve "
          "(no <12,20,30> missing-even-degree pattern).")
    print(f"   generic relaxed normalized low spectrum: "
          f"{np.round(res_generic['spec_relaxed'], 2)}")

    out = dict(N=N, iters=iters, stage0=g, substrate=chk,
               generic=res_generic, round_control=res_ctrl, defective=res_defect)
    json.dump(out, open(f"results_stage1_v11_N{N}.json", "w"), indent=2)
    print(f"\nwrote results_stage1_v11_N{N}.json")
    # headline
    print("\n==== STAGE 1 READ ====")
    print(f"  GENERIC : dist-to-round {res_generic['dist_init']:.3f} -> "
          f"{res_generic['dist_relaxed']:.3f}  (E x{res_generic['E_reduction']:.1f})")
    print(f"  DEFECTIVE: dist-to-round {res_defect['dist_init']:.3f} -> "
          f"{res_defect['dist_relaxed']:.3f}  (E x{res_defect['E_reduction']:.1f})")
    print(f"  CONTROL #3 round-init stays round: {res_ctrl['stays_round']}")


if __name__ == "__main__":
    main()
