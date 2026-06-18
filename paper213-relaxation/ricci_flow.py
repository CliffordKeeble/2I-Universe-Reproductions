"""
ricci_flow.py  —  Stage 1 of Paper 213 Approach B.

Normalized (volume-fixed) discrete Ricci flow on a frustrated metric: does the
GEOMETRY relax to round S^3?

SCHEME (documented, per brief §2 Stage 1 — Glickenstein/combinatorial discrete
Ricci flow is the chosen route; literal gradient-descent on the Regge action is a
saddle at the round metric — the conformal-mode problem — so we use the textbook
curvature-homogenizing flow instead):

  - discrete sectional curvature on edge e:   K_e = delta_e / a_e
       a_e = dual area to edge e = V_e^dual / ell_e,  V_e^dual = (1/6) sum_{t⊃e} V_t
       (barycentric dual; crude-and-honest). K_e is constant <=> round S^3.
  - normalized Ricci flow (multiplicative, keeps ell>0):
       d(log ell_e)/dtau = -(K_e - Kbar),   Kbar = volume-weighted mean curvature
       => high-curvature edges shrink, low-curvature grow; mean curvature preserved.
  - fix total volume each step (rescale ell). Fixed point: K_e = Kbar for all e = round.

CONVERGENCE METRIC (pre-registered): curvature coefficient of variation
       CoV = std(K_e) / |mean(K_e)|   ->  0  means round.
We do NOT tune the metric or the flow sign to the prediction; the chosen sign is the
standard homogenizing direction (shrink positive curvature), fixed before the run.

A tet becoming non-embeddable mid-flow = discrete neckpinch/glass signal (reported,
not patched).

Standalone:  python ricci_flow.py [N] [iters]
"""

import sys
import json
import numpy as np
from geometry import (generic_s3, verify_sphere, symmetry_probe,
                      edge_deficits, total_volume, tet_volume, TET_EDGES)


def edge_dual_volume(cx, ell):
    """V_e^dual = (1/6) sum_{t ⊃ e} V_t  (barycentric split of tet volume to its edges)."""
    Vd = np.zeros(len(cx.edges))
    for t in cx.tets:
        L = np.array([ell[cx.edge_index[tuple(sorted((t[a], t[b])))]]
                      for (a, b) in TET_EDGES])
        v = tet_volume(L)
        if np.isnan(v):
            continue
        for (a, b) in TET_EDGES:
            Vd[cx.edge_index[tuple(sorted((t[a], t[b])))]] += v / 6.0
    return Vd


def curvature(cx, ell):
    """K_e = delta_e / a_e, with a_e = V_e^dual / ell_e. Returns (K, deficits, Vdual, n_invalid)."""
    deficits, ninv = edge_deficits(cx, ell)
    Vd = edge_dual_volume(cx, ell)
    a_e = np.where(Vd > 1e-12, Vd / ell, np.nan)
    K = deficits / a_e
    return K, deficits, Vd, ninv


def cov_curvature(K, Vd):
    """Volume-weighted mean and coefficient of variation of K."""
    good = np.isfinite(K) & (Vd > 1e-12)
    w = Vd[good]
    Kg = K[good]
    Kbar = np.sum(w * Kg) / np.sum(w)
    var = np.sum(w * (Kg - Kbar) ** 2) / np.sum(w)
    cov = np.sqrt(var) / abs(Kbar) if abs(Kbar) > 1e-12 else np.inf
    return Kbar, cov


def run_flow(cx, iters=300, eta=0.05, tol=1e-5, seed_ell=None, verbose=True):
    ell = np.ones(len(cx.edges)) if seed_ell is None else seed_ell.copy()
    V_target = total_volume(cx, ell)
    hist = []
    prev_cov = None
    for it in range(iters):
        K, deficits, Vd, ninv = curvature(cx, ell)
        Kbar, cov = cov_curvature(K, Vd)
        hist.append(dict(it=it, cov=float(cov), Kbar=float(Kbar),
                         volume=float(total_volume(cx, ell)),
                         invalid_tets=int(ninv),
                         ell_min=float(ell.min()), ell_max=float(ell.max())))
        if verbose and (it % max(1, iters // 15) == 0 or it == iters - 1):
            print(f"  it={it:4d}  CoV(K)={cov:.4f}  Kbar={Kbar:+.4f}  "
                  f"V={hist[-1]['volume']:.3f}  invalid={ninv}  "
                  f"ell∈[{ell.min():.3f},{ell.max():.3f}]")
        if ninv > 0:
            hist[-1]["note"] = "non-embeddable tet -> neckpinch/glass signal"
            break
        if prev_cov is not None and abs(prev_cov - cov) < tol and cov < 0.5:
            hist[-1]["note"] = "converged (CoV change < tol)"
            break
        prev_cov = cov
        # normalized Ricci flow step (multiplicative)
        upd = np.where(np.isfinite(K), -(K - Kbar), 0.0)
        ell = ell * np.exp(eta * upd)
        # fix total volume
        Vnow = total_volume(cx, ell)
        if Vnow > 1e-9:
            ell *= (V_target / Vnow) ** (1.0 / 3.0)
    return ell, hist


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    print(f"Stage 1 — normalized discrete Ricci flow on generic S^3 (N={N})")
    print("  *** topology GRANTED (closed S^3 axiom); testing geometry relaxation ***")
    cx = generic_s3(N, seed=20260618)
    chk = verify_sphere(cx, max_link_checks=120)
    sym = symmetry_probe(cx)
    print(f"  substrate: edges={len(cx.edges)} tets={len(cx.tets)} "
          f"chi={chk['euler_char']} closed={chk['closed_no_boundary']} "
          f"links_ok={chk['links_ok']} 2I_regular={sym['looks_2I_regular']}")

    ell0 = np.ones(len(cx.edges))
    K0, _, Vd0, _ = curvature(cx, ell0)
    Kbar0, cov0 = cov_curvature(K0, Vd0)
    print(f"  initial frustrated metric: CoV(K)={cov0:.4f}  Kbar={Kbar0:+.4f}")
    print("  flowing...")
    ell, hist = run_flow(cx, iters=iters)

    cov_final = hist[-1]["cov"]
    drop = cov0 / cov_final if cov_final > 0 else np.inf
    print(f"\n  RESULT: CoV(K) {cov0:.4f} -> {cov_final:.4f}  "
          f"(factor {drop:.1f} reduction)  note={hist[-1].get('note','max-iters')}")

    out = dict(N=N, edges=len(cx.edges), tets=len(cx.tets),
               substrate=chk, symmetry=sym,
               cov_initial=float(cov0), cov_final=float(cov_final),
               Kbar_initial=float(Kbar0), Kbar_final=float(hist[-1]["Kbar"]),
               reduction_factor=float(drop), history=hist)
    with open(f"results_stage1_N{N}.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"  wrote results_stage1_N{N}.json")
    return out


if __name__ == "__main__":
    main()
