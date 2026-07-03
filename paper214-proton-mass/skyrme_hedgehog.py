"""
skyrme_hedgehog.py  —  B=1 Skyrme soliton, hedgehog ansatz (Paper 214, Phase 1).

Builds the classical B=1 Skyrmion in the hedgehog ansatz and reads off its
dimensionless energy. This is the *validation* layer of the brief: reproduce the
textbook number before anything 2I is attempted. If this file does not reproduce
the known dimensionless soliton energy, the code is wrong and downstream work is
meaningless (Pattern 75: validate before trust).

------------------------------------------------------------------------------
Model (static energy, brief eq.):
    E = integral d3x [ -(F_pi^2/16) Tr(L_i L_i) - (1/32 e^2) Tr([L_i,L_j]^2) ],
    L_i = U^dag d_i U.

Hedgehog ansatz:  U = exp( i tau . nhat F(r) ),  F(0)=pi, F(inf)=0.

Non-dimensionalise with  r = x * (2 / (F_pi e)).  The static energy reduces to

    E = (pi F_pi / e) * I[F],
    I[F] = integral_0^inf dx [ x^2 F'^2 + 2 sin^2 F (1 + F'^2) + sin^4 F / x^2 ].

Split into quadratic (sigma-model) and quartic (Skyrme) parts:
    I2 = integral [ x^2 F'^2 + 2 sin^2 F ]            (the F_pi^2 term)
    I4 = integral [ 2 sin^2 F F'^2 + sin^4 F / x^2 ]  (the 1/e^2 term)
    I  = I2 + I4.

Two convention-independent validation facts the solution MUST satisfy:
  * Derrick / virial theorem for a B=1 Skyrmion:  I2 = I4  (exactly, at the min).
  * Topological baryon number:  B = 1  (independent of the profile, by the BCs).

The Euler-Lagrange profile equation:
    (x^2 + 2 sin^2 F) F'' + 2 x F' + sin 2F (F'^2 - 1) - (sin^2 F sin 2F)/x^2 = 0.

------------------------------------------------------------------------------
INPUTS vs OUTPUTS (brief discipline):
  * F_pi (pion decay constant) and e (Skyrme coupling) are MEASURED/ASSUMED
    scales -- INPUT, not derived. They set the MeV scale only.
  * The geometry's genuine OUTPUT is the dimensionless number I_min (and hence
    M e / F_pi = pi I_min). The dimensionful mass is INPUT-scale x OUTPUT-number.

Standalone:  python skyrme_hedgehog.py
"""

import numpy as np
from scipy.integrate import solve_bvp


# --- domain ----------------------------------------------------------------
# x is the dimensionless radius. x0 small but nonzero (origin is a coordinate
# singularity of the reduced ODE); L large enough that F has decayed to ~0.
X0 = 1.0e-3
L = 30.0
N_MESH = 4000


def _profile_rhs(x, y):
    """RHS for solve_bvp.  y[0]=F, y[1]=F'.  Returns [F', F'']."""
    F, Fp = y
    s = np.sin(F)
    s2f = np.sin(2.0 * F)
    denom = x * x + 2.0 * s * s
    Fpp = (-2.0 * x * Fp - s2f * (Fp * Fp - 1.0) + (s * s * s2f) / (x * x)) / denom
    return np.vstack((Fp, Fpp))


def _profile_bc(ya, yb):
    """Boundary conditions: F(X0) = pi, F(L) = 0."""
    return np.array([ya[0] - np.pi, yb[0] - 0.0])


def solve_profile():
    """Solve the hedgehog BVP. Returns (x_grid, F, Fp) on a fine uniform mesh."""
    x = np.linspace(X0, L, 400)
    # Initial guess: a smooth pi -> 0 kink. Shape is not critical; solve_bvp
    # relaxes it to the true profile.
    F_guess = np.pi * np.exp(-0.8 * (x - X0))
    Fp_guess = np.gradient(F_guess, x)
    sol = solve_bvp(_profile_rhs, _profile_bc, x, np.vstack((F_guess, Fp_guess)),
                    tol=1e-8, max_nodes=200000, verbose=0)
    if not sol.success:
        raise RuntimeError(f"solve_bvp failed: {sol.message}")
    xg = np.linspace(X0, L, N_MESH)
    F, Fp = sol.sol(xg)
    return xg, F, Fp, sol


def energy_integrals(x, F, Fp):
    """Return (I2, I4, I) -- the dimensionless energy pieces, by quadrature."""
    s2 = np.sin(F) ** 2
    i2_dens = x * x * Fp * Fp + 2.0 * s2
    i4_dens = 2.0 * s2 * Fp * Fp + s2 * s2 / (x * x)
    I2 = np.trapezoid(i2_dens, x)
    I4 = np.trapezoid(i4_dens, x)
    return I2, I4, I2 + I4


def baryon_number(x, F, Fp):
    """Topological charge B = -(2/pi) integral sin^2 F F' dx.  Should equal 1."""
    dens = -(2.0 / np.pi) * np.sin(F) ** 2 * Fp
    return np.trapezoid(dens, x)


def bogomolny(x, F, Fp):
    """Topological (Faddeev-Bogomolny) lower bound and the completion residual.

    Identity:  I_density = (x F' + sin^2 F / x)^2 + 2 sin^2 F (1 + F')^2 - 6 sin^2 F F'.
    The two squares are >= 0; the last term integrates to -6 * (-pi/2) = +3 pi.
    Hence  I >= 3 pi,  and  I / (3 pi) -> 1.232 for the true B=1 profile (textbook).
    Returns (I_bound, ratio, completion_check) where completion_check should ~ 0.
    """
    s2 = np.sin(F) ** 2
    squares = (x * Fp + s2 / x) ** 2 + 2.0 * s2 * (1.0 + Fp) ** 2
    topo = -6.0 * s2 * Fp
    I_recon = np.trapezoid(squares + topo, x)
    I2, I4, I = energy_integrals(x, F, Fp)
    I_bound = 3.0 * np.pi
    completion_check = abs(I_recon - I) / I    # identity sanity, ~0
    return I_bound, I / I_bound, completion_check


def el_residual(x, F, Fp):
    """Max |EL residual| over the interior -- a check the ODE is actually solved."""
    Fpp = np.gradient(Fp, x)
    s2f = np.sin(2.0 * F)
    s = np.sin(F)
    res = ((x * x + 2.0 * s * s) * Fpp + 2.0 * x * Fp
           + s2f * (Fp * Fp - 1.0) - (s * s * s2f) / (x * x))
    interior = (x > 0.05) & (x < L - 0.5)
    return float(np.max(np.abs(res[interior])))


# --- INPUT scales (MEASURED/ASSUMED, not derived) --------------------------
# Adkins-Nappi-Witten (1983) best-fit parameter set. These are the canonical
# values used to quote the ~30%-level nucleon mass; they are tuned inputs.
F_PI_ANW_MEV = 129.0   # ANW best-fit pion decay constant (INPUT)
E_ANW = 5.45           # ANW best-fit Skyrme coupling (INPUT, dimensionless)
HBAR_C_MEV_FM = 197.3269804  # MeV.fm (CODATA), only for fm conversion


def classical_mass_mev(I_min, F_pi_mev=F_PI_ANW_MEV, e=E_ANW):
    """Classical soliton mass in MeV = (pi F_pi / e) * I_min.

    F_pi, e are INPUTS. I_min is the geometric OUTPUT.
    """
    return np.pi * F_pi_mev / e * I_min


def run(verbose=True):
    x, F, Fp, sol = solve_profile()
    I2, I4, I = energy_integrals(x, F, Fp)
    B = baryon_number(x, F, Fp)
    res = el_residual(x, F, Fp)
    I_bound, bogo_ratio, completion_check = bogomolny(x, F, Fp)
    M_dimless = np.pi * I            # M e / F_pi
    virial_rel = abs(I2 - I4) / I    # should be ~0
    M_cl = classical_mass_mev(I)

    out = {
        "I2_quadratic": I2,
        "I4_quartic": I4,
        "I_total": I,
        "M_over_Fpi_e": M_dimless,   # dimensionless: M e / F_pi = pi I
        "baryon_number": B,
        "virial_rel_gap": virial_rel,
        "bogomolny_bound_I": I_bound,
        "bogomolny_ratio": bogo_ratio,        # textbook 1.232
        "bogomolny_completion_check": completion_check,
        "el_max_residual": res,
        "Fprime_at_origin": float(Fp[0]),
        "classical_mass_MeV_ANW": M_cl,
        "inputs": {"F_pi_MeV": F_PI_ANW_MEV, "e": E_ANW,
                   "note": "F_pi and e are MEASURED/ASSUMED inputs, not derived"},
    }
    if verbose:
        print("=== B=1 hedgehog Skyrmion (Paper 214 Phase 1) ===")
        print(f"  I2 (quadratic / F_pi^2 term)   = {I2:.6f}")
        print(f"  I4 (quartic  / Skyrme term)    = {I4:.6f}")
        print(f"  I  = I2 + I4                   = {I:.6f}")
        print(f"  virial gap |I2-I4|/I           = {virial_rel:.3e}   (should be ~0)")
        print(f"  M e / F_pi  = pi I             = {M_dimless:.4f}  (ANW lit. 36.5)")
        print(f"  Bogomolny bound I_bound = 3 pi = {I_bound:.4f}")
        print(f"  Bogomolny ratio I/I_bound      = {bogo_ratio:.4f}  (textbook 1.232)")
        print(f"  completion identity check      = {completion_check:.3e}  (should be ~0)")
        print(f"  baryon number B                = {B:.8f}   (should be 1)")
        print(f"  max EL residual (interior)     = {res:.3e}")
        print(f"  F'(0)                          = {Fp[0]:.5f}")
        print(f"  classical mass (ANW inputs)    = {M_cl:.1f} MeV  [scale is INPUT]")
    return out, (x, F, Fp)


if __name__ == "__main__":
    run()
