"""
Independent verification of the three order-of-magnitude G estimates that
retire Paper 63 ("G and the Hierarchy Dilemma").

Three checks, all order-of-magnitude (we chase exponents, not decimals):

  A. Terrestrial envelope  -> kills the lab-variation thesis.
  B. Galactic orbit        -> the discriminator: local model excluded,
                              global/Mach model survives.
  C. Mach (Sciama) estimate -> the surviving order-of-magnitude result.

W-107 / decorrelation discipline: the inputs below were sourced
INDEPENDENTLY by Mr Code (web, June 2026), not lifted from the brief.
Where a value matched the brief's, that is corroboration, not copying.
Sources logged inline. Targets from the brief are printed alongside as
"expected" for comparison; PASS/FAIL is decided on the independent number.

Run:  python verify_G.py
"""

import math

# ----------------------------------------------------------------------
# Units / physical constants (CODATA / standard)
# ----------------------------------------------------------------------
G_meas   = 6.67430e-11        # m^3 kg^-1 s^-2   (CODATA 2018)
c        = 2.99792458e8       # m/s
m_p      = 1.67262192e-27     # kg              (proton mass)
M_sun    = 1.98892e30         # kg
pc       = 3.0857e16          # m               (parsec)
kpc      = 1e3 * pc
yr       = 3.1557e7           # s

# ----------------------------------------------------------------------
# Independently sourced astrophysical inputs (Mr Code, June 2026)
# ----------------------------------------------------------------------
# Oort total local volume density rho0 = 0.100 +/- 0.010 M_sun/pc^3
#   (Holmberg & Flynn / Hipparcos; corroborated by web search 2026-06).
rho0_tot   = 0.100            # M_sun / pc^3
# Baryonic fraction ~0.09 (known matter 0.108; DM ~0.008-0.01 => baryons ~0.09)
rho_b      = 0.09             # M_sun / pc^3
# Local dark matter 0.3 GeV/cm^3 = 0.008 M_sun/pc^3, subdominant & SMOOTH
rho_DM     = 0.008           # M_sun / pc^3
# Thin-disk stellar scale height ~300 pc (canonical; some fits 588 pc)
h_disk_pc  = 300.0
# Solar vertical "bob" amplitude through the plane
z_bob_pc   = 100.0
# Galactocentric orbit
R_orbit    = 8.0 * kpc
M_MW       = 1e11 * M_sun     # enclosed mass ~1e11 M_sun (order of mag)
ecc        = 0.06
T_gal      = 250e6 * yr       # galactic year ~250 Myr
# Universe
N_baryon   = 1e80
R_univ     = 1.3e26          # m  (~13.8 Gly comoving-ish scale used by paper)
# Lab geometry
d_lab      = 9.0e5           # m  (BIPM Paris <-> Florence ~900 km)
L_disk     = h_disk_pc * pc  # m  (~300 pc smallest plausible source scale)
obs_spread = 5e-4            # observed inter-lab G spread (paper Sec 4.2)
# Bounds (independently sourced)
GdotG_LLR  = 1e-12           # /yr  LLR: Gdot/G = (4+/-9)e-13 /yr -> ~1e-12 cap
BBN_bound  = 0.1             # DeltaG/G across BBN epoch
# Cosmological matter fractions for the total-matter cross-check
Omega_b_over_m = 1.0/6.0     # Omega_b/Omega_m ~ 1/6

print("=" * 70)
print("PAPER 63 RETIREMENT -- INDEPENDENT VERIFICATION (3 checks)")
print("=" * 70)

# ----------------------------------------------------------------------
# CHECK A -- terrestrial envelope
# ----------------------------------------------------------------------
print("\n[CHECK A] Terrestrial envelope (kills lab-variation)")
dG_over_G_A = d_lab / L_disk
L_eff = d_lab / obs_spread
print(f"  inputs: d_lab={d_lab:.1e} m, L_disk={L_disk:.2e} m (300 pc),"
      f" obs spread={obs_spread:.0e}")
print(f"  predicted dG/G = d/L          = {dG_over_G_A:.2e}   (expect ~1e-13)")
# at L = 1e18 m
print(f"               (at L=1e18 m)    = {d_lab/1e18:.2e}   (expect ~1e-12)")
print(f"  required source scale L_eff   = {L_eff:.2e} m  (expect ~2e9 m; no source)")
orders_below = math.log10(obs_spread / dG_over_G_A)
passA = dG_over_G_A <= 1e-12
print(f"  predicted is {orders_below:.1f} orders BELOW observed spread")
print(f"  PASS_A (predicted <= ~1e-12, i.e. >=8 orders below obs): {passA}")

# ----------------------------------------------------------------------
# CHECK B -- galactic orbit (the discriminator)
# ----------------------------------------------------------------------
print("\n[CHECK B] Galactic orbit (discriminator: local vs global)")

# --- Local model: G tracks local baryon density ---
# Density change as Sun bobs +/-z through disk of scale height h.
# exp model rho ~ exp(-|z|/h): fractional drop at z_bob.
frac_drop = 1.0 - math.exp(-z_bob_pc / h_disk_pc)
# DM is smooth & ~10% -> damps the *baryonic* swing only slightly.
dm_frac = rho_DM / rho0_tot
local_dG_over_G = frac_drop * (1.0 - dm_frac)   # damped by ~ DM fraction
print(f"  LOCAL model: dG/G ~ d(rho_b)/rho_b")
print(f"    exp-disk fractional density drop over +/-{z_bob_pc:.0f} pc "
      f"(h={h_disk_pc:.0f} pc) = {frac_drop:.2f}")
print(f"    DM fraction {dm_frac:.2f} smooth -> damps swing by ~{dm_frac*100:.0f}%")
print(f"    local dG/G ~ {local_dG_over_G:.2f}   (expect ~0.3, up to ~2)")

# --- Global (Machian) model: G^-1 ~ sum m_i/(r_i c^2) ---
# Term 1: bob slice potential-weighted share.
r_slice    = z_bob_pc * pc
M_slice    = (4.0/3.0) * math.pi * (z_bob_pc * pc)**3 \
             * (rho0_tot * M_sun / pc**3)            # kg
M_univ     = N_baryon * m_p                          # kg
slice_term = (M_slice / r_slice) / (M_univ / R_univ)
print(f"  GLOBAL model: dG/G ~ Machian potential-weighted slice share")
print(f"    M_slice = {M_slice/M_sun:.2e} M_sun (expect ~4e5),"
      f" r_slice={r_slice:.2e} m")
print(f"    bob-slice term = {slice_term:.2e}   (expect ~1e-10)")
# Term 2: eccentric-orbit MW-depth term.
mw_depth   = G_meas * M_MW / (R_orbit * c**2)
ecc_term   = mw_depth * ecc
print(f"    GM_MW/(R c^2) = {mw_depth:.2e}  (expect ~6e-7);"
      f"  x ecc({ecc}) = {ecc_term:.2e} (expect ~1e-8)")
global_dG_over_G = max(slice_term, ecc_term)
print(f"    global total dG/G ~ {global_dG_over_G:.2e}   (expect <= ~1e-8)")

# --- Bounds ---
LLR_per_galyr = GdotG_LLR * (T_gal / yr)
print(f"  BOUNDS: LLR {GdotG_LLR:.0e}/yr x {T_gal/yr:.1e} yr (gal yr)"
      f" = {LLR_per_galyr:.1e} per galactic year")
print(f"          BBN DeltaG/G <~ {BBN_bound} across ~55 orbits")
local_excluded   = local_dG_over_G > LLR_per_galyr * 100   # >~3 orders over
global_consistent = global_dG_over_G < LLR_per_galyr
print(f"  local ({local_dG_over_G:.2f}) exceeds LLR bound "
      f"({LLR_per_galyr:.1e}) by ~{math.log10(local_dG_over_G/LLR_per_galyr):.1f}"
      f" orders -> EXCLUDED: {local_excluded}")
print(f"  global ({global_dG_over_G:.1e}) below LLR bound -> CONSISTENT:"
      f" {global_consistent}")
passB = local_excluded and global_consistent
print(f"  PASS_B (local excluded AND global consistent): {passB}")

# ----------------------------------------------------------------------
# CHECK C -- Mach (Sciama) estimate
# ----------------------------------------------------------------------
print("\n[CHECK C] Mach / Sciama estimate (the surviving result)")
G_mach = R_univ * c**2 / (N_baryon * m_p)
ratio_b = G_mach / G_meas
G_mach_total = G_mach * Omega_b_over_m   # using total matter -> ~6x low
ratio_total = G_mach_total / G_meas
print(f"  G ~ R c^2 / (N m_p) = {G_mach:.3e}  (expect ~7.0e-11)")
print(f"    vs measured {G_meas:.3e}  ->  {abs(ratio_b-1)*100:.1f}% high"
      f"  (factor {ratio_b:.2f})")
print(f"  total-matter (x Omega_b/Omega_m ~1/6) = {G_mach_total:.3e}"
      f"  -> factor {ratio_total:.2f} ({1/ratio_total:.1f}x low)")
passC = (0.5 < ratio_b < 2.0) and (ratio_total < 0.5)
print(f"  PASS_C (baryon est within ~2x AND total-matter ~6x off): {passC}")

# ----------------------------------------------------------------------
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)
print(f"  Check A (terrestrial envelope): {'PASS' if passA else 'FAIL'}")
print(f"  Check B (galactic orbit):       {'PASS' if passB else 'FAIL'}")
print(f"  Check C (Mach estimate):        {'PASS' if passC else 'FAIL'}")
allpass = passA and passB and passC
print(f"\n  ALL THREE: {'PASS -> local model excluded, global/Mach survives' if allpass else 'NOT ALL PASS -- STOP, do not retire'}")
print("  (Order-of-magnitude estimates from standard inputs, NOT 2I derivations.)")
