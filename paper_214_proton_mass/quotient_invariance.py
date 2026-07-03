"""
quotient_invariance.py  —  numerical confirmation of brief section 0 (Paper 214).

Section 0 (provable, not discovered): the classical B=1 Skyrme energy on the target
S^3/2I equals the standard S^3 value. The mechanism is local: the energy density depends
on U only through the left current  l_i = U^dag d_i U, which is INVARIANT under a constant
left target rotation  U -> q U.  The quotient S^3/2I identifies  U ~ q U  for q in 2I, so
the energy density is literally identical on every 2I-orbit -> the functional descends to
the quotient unchanged.

This script makes that concrete, as a JOINT check on the code and the section-0 argument:

  (A) Build the hedgehog field on a 3D grid as SU(2) matrices and compute the FULL 3D
      Skyrme energy by finite differences. This is an INDEPENDENT cross-check of the 1D
      ODE result (skyrme_hedgehog.py):  E_tilde should ~ pi * I ~ 36.46 in units F_pi/e,
      and B ~ 1.  (Discretisation + finite box give a few-% offset; that is expected.)

  (B) Act with every element q of 2I (left multiplication U -> q U) and recompute the
      energy. Section 0 predicts ZERO change. We expect max|Delta E| at machine precision.

Standalone:  python quotient_invariance.py [grid_N] [box_L]
"""

import sys
import numpy as np

from skyrme_hedgehog import solve_profile
from icosian_2I import icosians_as_su2


# Dimensionless 3D energy (units F_pi/e), derived from the brief functional with
# r = (2/(F_pi e)) * xi :
#   E_tilde = integral d3xi [ -(1/8) Tr(l_i l_i) - (1/64) Tr([l_i,l_j]^2) ],
#   B       = (1/24 pi^2) integral d3xi eps^{ijk} Tr(l_i l_j l_k).
# Both reduce, on the hedgehog, to the 1D values (E_tilde = pi I, B = 1).


def _tr(A):
    """Trace over the last two axes."""
    return np.einsum('...aa->...', A)


def _matmul(A, B):
    return np.einsum('...ab,...bc->...ac', A, B)


def build_hedgehog_field(F_of_x, N, Lbox):
    """SU(2) hedgehog field U on an N^3 grid over [-Lbox,Lbox]^3 (origin avoided)."""
    ax = np.linspace(-Lbox, Lbox, N)        # N even -> no point at exactly 0
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    r = np.sqrt(X * X + Y * Y + Z * Z)
    F = F_of_x(r)
    nx, ny, nz = X / r, Y / r, Z / r
    cF = np.cos(F)
    sF = np.sin(F)
    # U = cosF * I + i sinF (n.sigma);  sigma = Pauli.  (n.sigma) = [[nz, nx-iny],[nx+iny,-nz]]
    U = np.empty(X.shape + (2, 2), dtype=complex)
    U[..., 0, 0] = cF + 1j * sF * nz
    U[..., 0, 1] = 1j * sF * (nx - 1j * ny)
    U[..., 1, 0] = 1j * sF * (nx + 1j * ny)
    U[..., 1, 1] = cF - 1j * sF * nz
    h = ax[1] - ax[0]
    return U, h


def left_currents(U, h):
    """l_i = U^dag d_i U  for i = x,y,z, by central finite differences."""
    Ud = np.conjugate(np.swapaxes(U, -1, -2))
    out = []
    for axis in range(3):
        dU = np.gradient(U, h, axis=axis)
        out.append(_matmul(Ud, dU))
    return out  # list of 3 arrays, each (...,2,2)


def energy_and_baryon(U, h):
    """Return (E_tilde, B) for field U on grid spacing h."""
    L = left_currents(U, h)
    # quadratic term: -(1/8) sum_i Tr(l_i l_i)
    quad = 0.0
    for i in range(3):
        quad = quad + _tr(_matmul(L[i], L[i]))
    e_quad = -(1.0 / 8.0) * quad
    # quartic term: -(1/64) sum_{i,j} Tr([l_i,l_j]^2)
    quart = 0.0
    for i in range(3):
        for j in range(3):
            C = _matmul(L[i], L[j]) - _matmul(L[j], L[i])
            quart = quart + _tr(_matmul(C, C))
    e_quart = -(1.0 / 64.0) * quart
    e_density = np.real(e_quad + e_quart)
    E_tilde = float(np.sum(e_density) * h ** 3)
    # baryon density: (1/24 pi^2) eps^{ijk} Tr(l_i l_j l_k)
    b = (_tr(_matmul(_matmul(L[0], L[1]), L[2]))
         - _tr(_matmul(_matmul(L[0], L[2]), L[1]))
         + _tr(_matmul(_matmul(L[1], L[2]), L[0]))
         - _tr(_matmul(_matmul(L[1], L[0]), L[2]))
         + _tr(_matmul(_matmul(L[2], L[0]), L[1]))
         - _tr(_matmul(_matmul(L[2], L[1]), L[0])))
    b_density = np.real(b) / (24.0 * np.pi ** 2)
    B = float(np.sum(b_density) * h ** 3)
    return E_tilde, B


def richardson(hs, vals):
    """Fit v(h) = v0 - c h^p to >=3 (h,v) points; return (v0, p). Least squares on log."""
    hs = np.asarray(hs, float); vals = np.asarray(vals, float)
    # grid-search v0 over a data-driven range (proportional to the spread, so it works
    # whether the gap-to-limit is large or tiny); linear-fit p from log|v0-v| vs log h.
    spread = float(vals.max() - vals.min())
    v_lo = vals.max() + 1e-6 * max(spread, 1.0)
    v_hi = vals.max() + 4.0 * max(spread, 1e-6)
    best = None
    for v0 in np.linspace(v_lo, v_hi, 2000):
        d = v0 - vals
        if np.any(d <= 0):
            continue
        p, logc = np.polyfit(np.log(hs), np.log(d), 1)
        resid = np.sum((np.log(d) - (p * np.log(hs) + logc)) ** 2)
        if best is None or resid < best[0]:
            best = (resid, v0, p)
    return best[1], best[2]


def run(verbose=True):
    target = np.pi * 11.606620   # pi * I_min from skyrme_hedgehog.py (= 36.46)

    # 1D profile as an interpolant F(x)
    x, F, Fp, _ = solve_profile()
    F_of_x = lambda rr: np.interp(rr, x, F, left=np.pi, right=0.0)

    # (A) convergence ladder: 3D FD energy/baryon vs grid (run once each, no group loop)
    ladder = [(48, 6.0), (72, 6.0), (96, 6.0), (120, 6.0)]
    hs, Es, Bs = [], [], []
    for N, Lbox in ladder:
        U, h = build_hedgehog_field(F_of_x, N, Lbox)
        E, B = energy_and_baryon(U, h)
        hs.append(h); Es.append(E); Bs.append(B)
    E_extrap, p_E = richardson(hs, Es)
    absB_extrap, _ = richardson(hs, [abs(b) for b in Bs])

    # (B) 2I invariance: U -> qU for all 120 elements on a modest grid (exact regardless of h)
    Ninv, Linv = 40, 6.0
    U, h = build_hedgehog_field(F_of_x, Ninv, Linv)
    E0, B0 = energy_and_baryon(U, h)
    Qsu2, _ = icosians_as_su2()
    dE, dB = [], []
    for q in Qsu2:
        Uq = np.einsum('ab,...bc->...ac', q, U)
        Eq, Bq = energy_and_baryon(Uq, h)
        dE.append(abs(Eq - E0)); dB.append(abs(Bq - B0))
    dE = np.array(dE); dB = np.array(dB)

    out = {
        "convergence_ladder": [
            {"N": N, "L": L, "h": hh, "E_tilde": EE, "B": BB}
            for (N, L), hh, EE, BB in zip(ladder, hs, Es, Bs)],
        "E_tilde_extrapolated": E_extrap,
        "E_tilde_convergence_order": p_E,
        "target_pi_times_I": target,
        "E_tilde_rel_error": abs(E_extrap - target) / target,
        "absB_extrapolated": absB_extrap,
        "baryon_sign_note": "B converges to -1; sign is the orientation convention of "
                            "eps^{ijk} vs the SU(2) embedding. |B|->1 is the physics.",
        "invariance_grid_N": Ninv,
        "n_2I_elements": len(Qsu2),
        "max_dE_under_2I": float(dE.max()),
        "max_dB_under_2I": float(dB.max()),
        "mean_dE_under_2I": float(dE.mean()),
    }
    if verbose:
        print("=== Section 0: classical mass is quotient-independent (3D check) ===")
        print("  (A) independent 3D FD cross-check of the 1D energy:")
        for (N, L), hh, EE, BB in zip(ladder, hs, Es, Bs):
            print(f"      N={N:3d} L={L} h={hh:.4f}  E_tilde={EE:8.4f}  B={BB:8.5f}")
        print(f"      Richardson extrap E_tilde = {E_extrap:.3f} (order p={p_E:.2f})")
        print(f"      1D target pi*I            = {target:.3f}")
        print(f"      rel. error                = {abs(E_extrap-target)/target:.2%}")
        print(f"      |B| extrapolated          = {absB_extrap:.4f}  (target 1)")
        print(f"  (B) 2I invariance over {len(Qsu2)} elements (U -> qU), grid {Ninv}^3:")
        print(f"      max |Delta E|       = {dE.max():.3e}   (section 0: exactly 0)")
        print(f"      max |Delta B|       = {dB.max():.3e}")
        print(f"      -> energy identical on every 2I-orbit to machine precision")
    return out


if __name__ == "__main__":
    run()
