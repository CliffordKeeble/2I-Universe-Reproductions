"""
controls.py -- Paper 213 Arm 1 controls (brief sec.5), run from the first rule.

  1. Non-golden : sqrt5 -> sqrt2, sqrt3, and a generic det(-1) reflection (a=7).
                  If d_s or any signature score moves materially, it was topology,
                  not gold.
  2. Scrambled  : random {I,R} per edge. Isolates coupling vs growth.
  3. Bare vs connection : L_s vs L_c side by side (always computed in run_arm1).
  4. d_max sweep: confirm d_s is not an artefact of the branching cap.

Each control is a growth spec; run_arm1 executes them and logs (d_s, signature)
so the signal can be read against the control band -- the null distribution for
the signature score (repo rule: publish the null, or don't publish the match).
"""

# Control specs: name -> kwargs override for growth.grow()
CONTROL_SPECS = [
    # 1. non-golden seeds (same topology rule, different reflection)
    {"name": "nongolden_sqrt2", "seed_a": 2.0, "C": "flat-triangles"},
    {"name": "nongolden_sqrt3", "seed_a": 3.0, "C": "flat-triangles"},
    {"name": "nongolden_a7",    "seed_a": 7.0, "C": "flat-triangles"},
    # 2. scrambled coupling
    {"name": "scrambled",       "scramble": True, "C": "flat-triangles"},
]

# 4. d_max sweep values
DMAX_SWEEP = [3, 4, 5, 6, 8]

# Closure-predicate variants (brief sec.3: ship >= 3, log outcomes)
C_VARIANTS = ["flat-triangles", "flat-shortest", "none"]
