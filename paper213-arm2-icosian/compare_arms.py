"""
compare_arms.py  —  THE DECISIVE CONTROL (Arm 2 brief §5).

Run Arm 1 (seed / {I,R} / R^2) and Arm 2 (icosian / 2I / R^4) on the SAME growth
rule, SAME seed, SAME topology (closure='no_constraint', so topology is generated
from the RNG alone and is coupling-independent). Assert L_s is byte-identical, then
compare the connection-Laplacian L_c signature scores.

PREDICTED: the <12,20,30> signature appears in ARM 2 only and is ABSENT in ARM 1.
If so, the honest finding is explicit: the signature is IMPORTED with 2I, not
emergent from golden coupling. That comparison IS the experiment's real result.

  *** 2I IS IMPORTED (Arm 2). This is reconstruction, not ignition. ***

Standalone:  python compare_arms.py            # default N
             python compare_arms.py 5000        # custom N
"""

import sys
import numpy as np

from growth import SeedCoupling, IcosianCoupling, grow, topology_fingerprint
from laplacian import scalar_laplacian, connection_laplacian
from spectral import low_eigenvalues, heat_trace_dimension, signature_score


def run_pair(N, d_max=4, k_children=2, p_close=0.5, seed=20260618, k_eig=800):
    """Grow both arms on identical topology; return their scores + shared d_s."""
    arm1 = SeedCoupling()
    arm2 = IcosianCoupling()

    G1, hol1 = grow(N, arm1, d_max=d_max, k_children=k_children,
                    p_close=p_close, closure="no_constraint", seed=seed)
    G2, hol2 = grow(N, arm2, d_max=d_max, k_children=k_children,
                    p_close=p_close, closure="no_constraint", seed=seed)

    # ---- topology must be identical (the whole point) ----
    fp1, fp2 = topology_fingerprint(G1), topology_fingerprint(G2)
    assert fp1 == fp2, f"TOPOLOGY DIVERGED across arms: {fp1} != {fp2}"
    Ls1 = scalar_laplacian(G1)
    Ls2 = scalar_laplacian(G2)
    assert abs(Ls1 - Ls2).max() == 0.0, "L_s differs across arms — topology not shared"

    # ---- shared spectral dimension (from the common L_s) ----
    eig_Ls = low_eigenvalues(Ls1, k=min(k_eig, N - 2))
    d_s, d_s_err, _ = heat_trace_dimension(eig_Ls, N)

    # ---- connection Laplacians (the holonomy content) ----
    Lc1 = connection_laplacian(G1, hol1, arm1.fibre_dim)
    Lc2 = connection_laplacian(G2, hol2, arm2.fibre_dim)
    eig_Lc1 = low_eigenvalues(Lc1, k=min(k_eig, Lc1.shape[0] - 2))
    eig_Lc2 = low_eigenvalues(Lc2, k=min(k_eig, Lc2.shape[0] - 2))

    # ---- bare-vs-connection control: score L_s too (must show no signature) ----
    sc_Ls = signature_score(eig_Ls)
    sc1 = signature_score(eig_Lc1)
    sc2 = signature_score(eig_Lc2)

    return dict(
        N=N, d_max=d_max, k_children=k_children, p_close=p_close, seed=seed,
        edges=G1.number_of_edges(),
        d_s=float(d_s), d_s_err=float(d_s_err),
        arm1_Lc_score=sc1["signature_score"], arm1_Lc=sc1,
        arm2_Lc_score=sc2["signature_score"], arm2_Lc=sc2,
        Ls_score=sc_Ls["signature_score"], Ls=sc_Ls,
        topology_identical=True,
    )


def _fmt(x):
    return "nan" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:+.3f}"


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    print("=" * 70)
    print("Paper 213 — compare_arms : Arm 1 (seed/{I,R}) vs Arm 2 (icosian/2I)")
    print("  *** 2I IS IMPORTED in Arm 2 — reconstruction test, not ignition ***")
    print("=" * 70)
    r = run_pair(N)
    print(f"N = {r['N']}   edges = {r['edges']}   topology identical = {r['topology_identical']}")
    print(f"shared d_s (from common L_s) = {r['d_s']:.3f} ± {r['d_s_err']:.3f}")
    print("-" * 70)
    print(f"{'object':<22}{'sig score':>12}{'frac_sup':>11}{'frac_forb':>11}{'l2/l1':>9}")
    print(f"{'L_s (bare, shared)':<22}{_fmt(r['Ls_score']):>12}"
          f"{_fmt(r['Ls']['frac_supported']):>11}{_fmt(r['Ls']['frac_forbidden']):>11}"
          f"{_fmt(r['Ls']['ratio_obs']):>9}")
    print(f"{'L_c ARM 1 ({I,R})':<22}{_fmt(r['arm1_Lc_score']):>12}"
          f"{_fmt(r['arm1_Lc']['frac_supported']):>11}{_fmt(r['arm1_Lc']['frac_forbidden']):>11}"
          f"{_fmt(r['arm1_Lc']['ratio_obs']):>9}")
    print(f"{'L_c ARM 2 (2I import)':<22}{_fmt(r['arm2_Lc_score']):>12}"
          f"{_fmt(r['arm2_Lc']['frac_supported']):>11}{_fmt(r['arm2_Lc']['frac_forbidden']):>11}"
          f"{_fmt(r['arm2_Lc']['ratio_obs']):>9}")
    print("-" * 70)
    print(f"target ratio l2/l1 = {r['arm2_Lc']['ratio_target']:.3f} (440/168)")
    print("=" * 70)
    return r


if __name__ == "__main__":
    main()
