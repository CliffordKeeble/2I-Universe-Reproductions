"""
controls.py  —  Arm 2 controls (pre-registered, Arm 2 brief §5).

  scrambled_2I_null : random icosian per edge on the SAME topology -> the matched
                      null for the 2I-import score (the rigorous "publish the null").
  rho_variants      : default {s,t}-by-parity vs deterministic context-hash into 2I.
  closure_variants  : no_constraint / flat_triangles / flat_shortest.

The Arm-1 head-to-head (the decisive control) lives in compare_arms.py.
"""

import numpy as np
from growth import IcosianCoupling, grow
from laplacian import connection_laplacian
from spectral import low_eigenvalues, signature_score


def scrambled_2I_score(G, mats, rng, k_eig=200):
    hol = {}
    for (u, v) in G.edges():
        U = mats[rng.integers(len(mats))]
        hol[(u, v)] = U
        hol[(v, u)] = U.T
    Lc = connection_laplacian(G, hol, 4)
    e = low_eigenvalues(Lc, k=k_eig)
    return signature_score(e)["signature_score"]


def scramble_null(G, mats, n=12, seed=12345, k_eig=200):
    rng = np.random.default_rng(seed)
    vals = [scrambled_2I_score(G, mats, rng, k_eig) for _ in range(n)]
    vals = np.array([v for v in vals if not np.isnan(v)])
    return vals
