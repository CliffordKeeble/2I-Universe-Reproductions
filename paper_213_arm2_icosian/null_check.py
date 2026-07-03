"""
null_check.py  —  does each arm's signature score beat its OWN matched null?

Flag-event diagnostic. The pre-registered score is FROZEN; we only ask whether the
observed score separates from chance, per repo law ("publish the null or don't
publish the match"). Each arm is compared against its own scrambled-coupling null
on the SAME topology:

  Arm 1 : structured {I,R}      vs  random {I,R} per edge        (2-D fibres)
  Arm 2 : structured {s,t}<2I   vs  random icosian per edge      (4-D fibres)

Decides the escalation question: is Arm 1's high score a real <12,20,30> signature
(-> W-108 stop, per pre-registered criteria) or an artefact of the (non-orthogonal)
seed / fibre structure (-> the score is informative but Arm 1 carries structure for
a mundane reason)?
"""

import numpy as np
from growth import SeedCoupling, IcosianCoupling, grow
from laplacian import connection_laplacian
from spectral import low_eigenvalues, signature_score

SEED = 20260618
N = 2000
K_EIG = 200
N_NULL = 20
RNG = np.random.default_rng(12345)


def scrambled_holonomy(G, mats):
    hol = {}
    for (u, v) in G.edges():
        U = mats[RNG.integers(len(mats))]
        hol[(u, v)] = U
        hol[(v, u)] = U.T
    return hol


def score_of(G, hol, f):
    Lc = connection_laplacian(G, hol, f)
    e = low_eigenvalues(Lc, k=K_EIG)
    return signature_score(e)


def null_dist(G, mats, f, n=N_NULL):
    out = []
    for _ in range(n):
        hol = scrambled_holonomy(G, mats)
        out.append(score_of(G, hol, f)["signature_score"])
    return np.array([x for x in out if not np.isnan(x)])


def report(label, obs, nd):
    if len(nd) == 0:
        print(f"\n{label}: null EMPTY (all nan)"); return
    m, sd = nd.mean(), nd.std()
    z = (obs - m) / sd if sd > 0 else float("nan")
    pct = (nd <= obs).mean() * 100
    print(f"\n{label}")
    print(f"   observed score = {obs:+.3f}")
    print(f"   matched null   : mean={m:+.3f} std={sd:.3f} "
          f"min={nd.min():+.3f} max={nd.max():+.3f} (n={len(nd)})")
    print(f"   separation     : z={z:+.2f}, {pct:.0f}th percentile of its own null")


def main():
    arm1, arm2 = SeedCoupling(), IcosianCoupling()
    G1, hol1 = grow(N, arm1, closure="no_constraint", seed=SEED)
    G2, hol2 = grow(N, arm2, closure="no_constraint", seed=SEED)

    s1 = score_of(G1, hol1, 2)
    s2 = score_of(G2, hol2, 4)

    print("OBSERVED (structured coupling):")
    print(f"  Arm1 {{I,R}}  score={s1['signature_score']:+.3f}  "
          f"frac_sup={s1['frac_supported']:+.2f} frac_forb={s1['frac_forbidden']:+.2f} "
          f"l2/l1={s1['ratio_obs']:.3f}")
    print(f"  Arm2 2I     score={s2['signature_score']:+.3f}  "
          f"frac_sup={s2['frac_supported']:+.2f} frac_forb={s2['frac_forbidden']:+.2f} "
          f"l2/l1={s2['ratio_obs']:.3f}")
    print(f"  Arm1 scaled low clusters (->168): "
          f"{[round(x,1) for x in s1['scaled_low']]}")
    print(f"  Arm2 scaled low clusters (->168): "
          f"{[round(x,1) for x in s2['scaled_low']]}")
    print(f"  (target ladder: 168, 440, 624, 960, 1088, 1368, ...)")

    # Arm 1 alphabet {I, R}; Arm 2 alphabet = full 120 icosians
    arm1_alpha = np.array([np.eye(2), arm1.R])
    nd1 = null_dist(G1, arm1_alpha, 2)
    nd2 = null_dist(G2, arm2._mats, 4)

    report("ARM 1  {I,R}  vs random-{I,R} null", s1["signature_score"], nd1)
    report("ARM 2  2I     vs random-icosian null", s2["signature_score"], nd2)

    print("\n" + "=" * 64)
    print("READ: if Arm1 obs sits INSIDE its own {I,R} null -> artefact of the")
    print("2-D / non-orthogonal-seed structure, NOT a real signature (no W-108).")
    print("If Arm1 obs is far ABOVE its own {I,R} null -> real signature -> W-108.")


if __name__ == "__main__":
    main()
