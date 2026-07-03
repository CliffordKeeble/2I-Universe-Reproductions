"""
live_cell_recount.py -- Paper 101 v3.1 section 5.3, Mr A's null-arithmetic fix
(FIX D).  Recount the §5.3 null's *denominator* honestly.

Mr A's observation (third cold review).  The §5.3 headline "5,400 cells, one
survivor" flatters the null.  The decision window is 0.01% relative.  Every
candidate value lambda_l +/- c is an INTEGER (integer eigenvalue minus integer
correction).  An integer can only land inside the window [t(1-w), t(1+w)] if the
target t has an integer within t*w of it.  At 0.01%, t*w ~ 0.18 at the proton
scale.  Of the five targets only m_p/m_e sits ~0.15 from an integer; the others
do not sit near an integer at that precision and therefore CANNOT FIRE -- for
reasons nothing to do with the framework.  Their cells were dead before the run.
The honest denominator is the number of *live* cells.

This script does NOT change the decision rule (one expected hit per target) or
the verdict (one survivor).  It corrects the reported denominator.

Two liveness tests are reported, because they DIVERGE for one target and the
divergence is itself the finding:

  (A)  NEAR-INTEGER  (Mr A's stated envelope, NECESSARY condition):
       target t admits an in-window integer iff dist(t, nearest int) <= t*w.

  (B)  REACHABLE  (the FULL condition, definition (2) of the brief):
       some actual candidate (lambda_l +/- c) with l a surviving 2I level and c
       in the correction family lands within t*w of t.  (B) implies (A); the
       converse can fail when the in-window integer is not a producible
       candidate value.

A cell (level l, correction c, sign s, target t) is LIVE iff its target is
REACHABLE -- i.e. the target admits at least one in-window candidate.  All cells
on an unreachable target are dead before the run.  N_live = (# reachable
targets) * (cells per target).

Run on both correction families of §5.3:
  * flat declared-integer set  {1,2,6,11,12,13,18,19,20,30,31,41,42,62,168}
  * orbit-index set            {1,5,6,10,12,15,20,24,30,40,60,120}

Reuses massratio_null (TARGETS, OPERATIONS, survivors_at, CORRECTIONS),
spectrum (surviving_levels), and orbit_index_null (orbit-index enumeration).

Run: python live_cell_recount.py
"""

import csv

from massratio_null import (TARGETS, OPERATIONS, CORRECTIONS as FLAT_CORR,
                            survivors_at)
from spectrum import surviving_levels
from orbit_index_null import all_subgroups, build_2I

PRIMARY = 1e-4   # 0.01% relative -- the locked §5.3 decision window


def realised_orbit_indices():
    """The DERIVED orbit-index set [2I:H] = 120/|H| (gates findings Task 1)."""
    G = build_2I()
    subs = all_subgroups(G)
    orders = {len(H) for H in subs}
    return sorted({120 // o for o in orders})


def nearest_integer_report(t, window):
    """Distance of t to its nearest integer and the window half-width t*window."""
    nint = round(t)
    dist = abs(t - nint)
    halfwidth = t * window
    near = dist <= halfwidth
    return nint, dist, halfwidth, near


def build_combos(corrections):
    levels_raw, _ = surviving_levels(120)
    levels = [(l, lam) for (l, lam, mult) in levels_raw if l > 0]
    combos = []
    for (l, lam) in levels:
        for (_opsym, opsign) in OPERATIONS:
            for c in corrections:
                combos.append((l, lam, _opsym, c, lam + opsign * c))
    return levels, combos


def run_family(name, corrections, window, csv_rows):
    levels, combos = build_combos(corrections)
    n_levels = len(levels)
    n_corr = len(corrections)
    cells_per_target = n_levels * n_corr * len(OPERATIONS)
    full_grid = cells_per_target * len(TARGETS)

    # (B) REACHABLE: which targets have >= 1 candidate inside the window.
    res = survivors_at(combos, window)
    fired = {t: len(res[t]) for t in TARGETS}             # cells that FIRE
    reachable = {t: fired[t] >= 1 for t in TARGETS}

    print("=" * 74)
    print(f"FAMILY: {name}")
    print("=" * 74)
    print(f"surviving levels (l>0)   : {n_levels}")
    print(f"corrections ({n_corr:2d})        : {corrections}")
    print(f"operations               : {len(OPERATIONS)} (+/-)")
    print(f"targets                  : {len(TARGETS)}")
    print(f"cells per target         : {cells_per_target} "
          f"(= {n_levels} x {n_corr} x {len(OPERATIONS)})")
    print(f"FULL GRID                : {full_grid}  "
          f"(= {cells_per_target} x {len(TARGETS)} targets)")

    print(f"\nper-target liveness (window {window*100:.3f}%):")
    print(f"  {'target':10s}  {'value':>13s}  {'near int':>8s}  "
          f"{'dist':>9s}  {'t*w':>8s}  {'(A)near':>7s}  {'(B)reach':>8s}  "
          f"{'fired':>5s}")
    near_targets = []
    reach_targets = []
    for t, tv in TARGETS.items():
        nint, dist, hw, near = nearest_integer_report(tv, window)
        if near:
            near_targets.append(t)
        if reachable[t]:
            reach_targets.append(t)
        print(f"  {t:10s}  {tv:13.6f}  {nint:8d}  {dist:9.5f}  {hw:8.5f}  "
              f"{'YES' if near else 'no':>7s}  "
              f"{'YES' if reachable[t] else 'no':>8s}  {fired[t]:5d}")
        csv_rows.append((name, t, f"{tv:.6f}", nint, f"{dist:.6f}",
                         f"{hw:.6f}", int(near), int(reachable[t]), fired[t],
                         cells_per_target))

    n_live_targets = len(reach_targets)
    n_live_cells = n_live_targets * cells_per_target
    n_dead_cells = full_grid - n_live_cells
    total_fired = sum(fired.values())

    print(f"\n  (A) near-integer targets : {len(near_targets)}  {near_targets}")
    print(f"  (B) reachable targets    : {n_live_targets}  {reach_targets}")
    divergence = sorted(set(near_targets) - set(reach_targets))
    if divergence:
        print(f"  >>> (A)\\(B) divergence   : {divergence}  "
              f"(near an integer, but that integer is NOT a producible "
              f"candidate value -> dead)")

    print(f"\n  LIVE CELLS   N_live      : {n_live_cells}  "
          f"(= {n_live_targets} reachable target(s) x {cells_per_target})")
    print(f"  DEAD CELLS               : {n_dead_cells}  "
          f"(= {full_grid} - {n_live_cells})")
    print(f"  CELLS THAT FIRED         : {total_fired}")

    # The honest headline line for §5.3.
    fired_detail = []
    for t in reach_targets:
        for (l, lam, opsym, c, value, rel) in res[t]:
            fired_detail.append(f"{lam} {opsym} {c} = {value:g} on {t}")
    headline = (
        f"of {n_live_cells} live cells across the reachable target(s) "
        f"{reach_targets}, exactly {total_fired} fired "
        f"({'; '.join(fired_detail)}); the remaining {n_dead_cells} cells "
        f"were on targets no integer-minus-integer value reaches at "
        f"{window*100:.2f}% and were dead before the run.")
    print(f"\n  HEADLINE: {headline}")
    print()

    return {
        "name": name,
        "full_grid": full_grid,
        "cells_per_target": cells_per_target,
        "n_live_targets": n_live_targets,
        "reach_targets": reach_targets,
        "near_targets": near_targets,
        "divergence": divergence,
        "n_live_cells": n_live_cells,
        "n_dead_cells": n_dead_cells,
        "total_fired": total_fired,
        "headline": headline,
    }


def main():
    print("#" * 74)
    print("# Paper 101 v3.1 section 5.3 -- the live-cell recount (FIX D)")
    print("# Mr A's null-arithmetic fix: report the null's TRUE denominator")
    print("#" * 74)
    print()

    # Per-target nearest-integer envelope at the primary window (family-agnostic).
    print("Per-target nearest-integer distances (window 0.01%, t*w half-width):")
    print(f"  {'target':10s}  {'value':>13s}  {'nearest int':>11s}  "
          f"{'dist':>9s}  {'t*w':>9s}  {'integer in window?':>18s}")
    for t, tv in TARGETS.items():
        nint, dist, hw, near = nearest_integer_report(tv, PRIMARY)
        print(f"  {t:10s}  {tv:13.6f}  {nint:11d}  {dist:9.5f}  {hw:9.5f}  "
              f"{'YES' if near else 'no':>18s}")
    print()

    orbit = realised_orbit_indices()
    assert 12 in orbit
    csv_rows = []

    summaries = []
    summaries.append(run_family("flat declared-integer set (15)",
                                FLAT_CORR, PRIMARY, csv_rows))
    summaries.append(run_family("orbit-index set (12)  [headline 5,400]",
                                orbit, PRIMARY, csv_rows))

    print("=" * 74)
    print("SUMMARY -- the recounted denominators")
    print("=" * 74)
    for s in summaries:
        print(f"\n{s['name']}:")
        print(f"  full grid (most cells unreachable) : {s['full_grid']}")
        print(f"  live cells N_live                  : {s['n_live_cells']}")
        print(f"  dead-before-run                    : {s['n_dead_cells']}")
        print(f"  survivors (fired)                  : {s['total_fired']}")
        print(f"  reachable targets                  : {s['reach_targets']}")
        if s["divergence"]:
            print(f"  near-integer-but-dead (flag)       : {s['divergence']}")

    with open("live_cell_recount.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["family", "target", "value", "nearest_int", "dist",
                    "window_halfwidth", "near_integer_A", "reachable_B",
                    "n_fired", "cells_per_target"])
        w.writerows(csv_rows)
    print("\nwrote live_cell_recount.csv")


if __name__ == "__main__":
    main()
