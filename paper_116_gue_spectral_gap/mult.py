"""
Reuse adapter (Gate G5): multiplicities m_Gamma(l) imported from the verified
Paper 118 modules and exposed to arbitrary l_max for the Paper 116 spectral-
statistics pipeline. No multiplicity arithmetic is reimplemented here -- every
value comes from paper_118_eleven_barrier via its exact "route B".

Population is the 118-CODE-1 Part 2 population (90 solvable + 2I):
    59 cyclic     Z_2 .. Z_60
    29 dihedral   2D_2 .. 2D_30
     3 polyhedral 2T, 2O, 2I
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_P118 = os.path.normpath(os.path.join(_HERE, "..", "paper_118_eleven_barrier"))
if _P118 not in sys.path:
    sys.path.insert(0, _P118)

import solvable_control_null as scn   # verified 91-group machinery (commit feeb5da)
import eleven_barrier as eb           # 2I three-route cross-check (commit 96a290f)


def m_exact_series(name, kind, param, classes, order, l_max):
    """m_Gamma(l) for l = 0..l_max via 118's exact route B (never a fork)."""
    if kind == "poly":
        mol = scn.molien(*scn.MOLIEN[name], l_max)
        return [int(mol[l]) for l in range(l_max + 1)]
    if kind == "cyclic":
        return [int(scn.m_cyclic_exact(param, l)) for l in range(l_max + 1)]
    if kind == "dihedral":
        return [int(scn.m_dihedral_exact(param, l)) for l in range(l_max + 1)]
    raise ValueError(f"unknown kind {kind!r}")


def m_char_series(name, kind, param, classes, order, l_max):
    """Independent route A (mpmath character sum) for cross-checks (AT1)."""
    out = []
    for l in range(l_max + 1):
        mr, _res = scn.m_char(classes, order, l)
        out.append(int(mr))
    return out


def population():
    """118 Part-2 population as (name, kind, param, classes, order, p1) tuples."""
    return scn.build_population()


def group_by_name(name):
    for g in population():
        if g[0] == name:
            return g
    raise KeyError(name)


def m_2I_eleven_barrier(l_max):
    """2I multiplicities via the *other* module (eleven_barrier Molien) -- AT1 anchor."""
    return [int(x) for x in eb.molien_coeffs(l_max)]
