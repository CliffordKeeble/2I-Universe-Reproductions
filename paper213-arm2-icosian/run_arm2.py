"""
run_arm2.py  —  Arm 2 production driver (finite-size scaling + rule-space sweep).

  *** 2I IS IMPORTED. Reconstruction test, not ignition. ***

Outputs:
  results.json              full numeric record
  logs/rule_trials.jsonl    every (rho, closure, d_max) variant + (d_s, score)
  figures/arm2_convergence.png   d_s vs N and Arm-2 score vs N with the null band

Bounded for a first pass: N up to 4000; the expensive 1e5-1e6 scaling is gated on
CinC's call once the substrate question (d_s far below 3) is resolved.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from growth import SeedCoupling, IcosianCoupling, grow, topology_fingerprint
from laplacian import scalar_laplacian, connection_laplacian
from spectral import low_eigenvalues, heat_trace_dimension, signature_score
from controls import scramble_null

SEED = 20260618
K_EIG = 250


def one_config(N, d_max, k_children, p_close, closure, rho_mode="parity",
               with_null=False, null_n=10):
    """Grow both arms (same topology under no_constraint) and score."""
    arm1, arm2 = SeedCoupling(), IcosianCoupling()
    if rho_mode == "hash":
        arm2.rho = arm2.rho_hash  # context-hash variant

    G1, h1 = grow(N, arm1, d_max=d_max, k_children=k_children,
                  p_close=p_close, closure=closure, seed=SEED)
    G2, h2 = grow(N, arm2, d_max=d_max, k_children=k_children,
                  p_close=p_close, closure=closure, seed=SEED)

    Ls = scalar_laplacian(G2)
    topo_same = topology_fingerprint(G1) == topology_fingerprint(G2)
    eig_Ls = low_eigenvalues(Ls, k=min(K_EIG, N - 2))
    d_s, d_s_err, _ = heat_trace_dimension(eig_Ls, N)

    Lc1 = connection_laplacian(G1, h1, 2)
    Lc2 = connection_laplacian(G2, h2, 4)
    s1 = signature_score(low_eigenvalues(Lc1, k=min(K_EIG, Lc1.shape[0] - 2)))
    s2 = signature_score(low_eigenvalues(Lc2, k=min(K_EIG, Lc2.shape[0] - 2)))

    rec = dict(N=N, d_max=d_max, k_children=k_children, p_close=p_close,
               closure=closure, rho_mode=rho_mode,
               edges=G2.number_of_edges(),
               topology_identical=bool(topo_same),
               d_s=float(d_s), d_s_err=float(d_s_err),
               arm1_score=float(s1["signature_score"]),
               arm2_score=float(s2["signature_score"]),
               arm2_frac_sup=float(s2["frac_supported"]),
               arm2_frac_forb=float(s2["frac_forbidden"]),
               arm2_ratio=float(s2["ratio_obs"]))
    if with_null:
        nd = scramble_null(G2, arm2._mats, n=null_n, k_eig=K_EIG)
        rec["null_mean"] = float(nd.mean()) if len(nd) else float("nan")
        rec["null_std"] = float(nd.std()) if len(nd) else float("nan")
        rec["null_max"] = float(nd.max()) if len(nd) else float("nan")
        rec["arm2_z"] = (float((s2["signature_score"] - nd.mean()) / nd.std())
                         if len(nd) and nd.std() > 0 else float("nan"))
    return rec


def main():
    results = {"meta": "Paper 213 Arm 2 — 2I IMPORTED — reconstruction test",
               "fss": [], "dmax_sweep": [], "closure_sweep": [], "rho_sweep": []}
    trials = []

    print("== finite-size scaling (d_max=4, no_constraint) ==")
    for N in (1000, 2000, 4000):
        r = one_config(N, 4, 2, 0.5, "no_constraint", with_null=True, null_n=10)
        results["fss"].append(r); trials.append(r)
        print(f"  N={N:5d}  d_s={r['d_s']:.3f}±{r['d_s_err']:.3f}  "
              f"arm2={r['arm2_score']:+.3f} (null {r['null_mean']:+.3f}"
              f"±{r['null_std']:.3f}, z={r['arm2_z']:+.2f})  arm1={r['arm1_score']:+.3f}")

    print("== d_max sweep (N=2000) ==")
    for dm in (3, 4, 5, 6):
        r = one_config(2000, dm, 2, 0.5, "no_constraint")
        results["dmax_sweep"].append(r); trials.append(r)
        print(f"  d_max={dm}  d_s={r['d_s']:.3f}±{r['d_s_err']:.3f}  "
              f"arm2={r['arm2_score']:+.3f}  arm1={r['arm1_score']:+.3f}")

    print("== closure variant sweep (N=2000) ==")
    for cl in ("no_constraint", "flat_triangles", "flat_shortest"):
        r = one_config(2000, 4, 2, 0.5, cl)
        results["closure_sweep"].append(r); trials.append(r)
        print(f"  {cl:16s}  d_s={r['d_s']:.3f}  arm2={r['arm2_score']:+.3f}  "
              f"topo_same={r['topology_identical']}")

    print("== rho variant sweep (N=2000) ==")
    for rm in ("parity", "hash"):
        r = one_config(2000, 4, 2, 0.5, "no_constraint", rho_mode=rm)
        results["rho_sweep"].append(r); trials.append(r)
        print(f"  rho={rm:7s}  arm2={r['arm2_score']:+.3f}")

    # write logs
    with open("logs/rule_trials.jsonl", "w") as f:
        for t in trials:
            f.write(json.dumps(t) + "\n")
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)

    # figure
    fss = results["fss"]
    Ns = [r["N"] for r in fss]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(Ns, [r["d_s"] for r in fss], "o-", label="d_s (shared)")
    ax[0].axhline(3.0, ls="--", c="grey", label="3 (manifold target)")
    ax[0].set_xscale("log"); ax[0].set_xlabel("N"); ax[0].set_ylabel("d_s")
    ax[0].set_title("Spectral dimension vs N\n[2I IMPORTED — Arm 2]"); ax[0].legend()

    a2 = [r["arm2_score"] for r in fss]
    nm = [r["null_mean"] for r in fss]; ns = [r["null_std"] for r in fss]
    ax[1].plot(Ns, a2, "s-", c="C3", label="Arm 2 (2I) score")
    ax[1].plot(Ns, [r["arm1_score"] for r in fss], "^-", c="C0", label="Arm 1 ({I,R})")
    ax[1].fill_between(Ns, np.array(nm) - 2*np.array(ns), np.array(nm) + 2*np.array(ns),
                       color="grey", alpha=0.3, label="scramble null ±2σ")
    ax[1].axhline(1.0, ls="--", c="green", label="reconstruction = 1")
    ax[1].set_xscale("log"); ax[1].set_xlabel("N"); ax[1].set_ylabel("signature score")
    ax[1].set_title("⟨12,20,30⟩ score vs N\n[2I IMPORTED — reconstruction, not ignition]")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/arm2_convergence.png", dpi=120)
    print("\nwrote results.json, logs/rule_trials.jsonl, figures/arm2_convergence.png")


if __name__ == "__main__":
    main()
