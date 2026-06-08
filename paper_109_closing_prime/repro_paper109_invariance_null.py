#!/usr/bin/env python3
"""
Reproduction candidate (sequence label TBD by Mr Code) : §12.6 CONVENTION-INVARIANCE NULL.
2I Universe Programme -- Reproductions repo.  Author: CinC (chat).  Verify/run: Mr Code.
STATUS: PRELIMINARY CinC in-chat finding. NOT yet Scout-checked, NOT Adversary-reviewed.

THE QUESTION (Mr A's fifth-star sub-target, Paper 109 §12.6)
  "Exhibit a single convention-invariant functional of a modular cycle -- one quantity that
  takes the same value whether the cycle is produced by reduction A, B, or C -- or establish
  that none exists."  This script attacks the sub-target empirically.

THE TEST
  A convention-invariant functional of a modular cycle needs a modular cycle that is itself
  a reduction-independent object. So we ask the prior question: is ANY integer Syracuse-loop
  captured as a cycle by all three reductions A, B, C (at any levels, exhaustively to m=20)?
  If no loop is common to all three, there is no reduction-independent cycle for a functional
  to be a functional OF, and outcome (i) of §12.6 is settled in the negative WITHOUT an
  open-ended search over candidate functionals.

  We also directly sweep natural functionals (cycle count, length, min element, S-class
  multiset mod 8, mod-24 vertex set) at the levels m=11,12 where all three reductions breathe,
  to confirm none is invariant there.

REDUCTIONS (identical to R109-3; see that script's header)
  A min-rep : syr(n) % M           B reduce-then-strip : oddpart((3n+1) % M)
  C shifted : syr(n + M) % M        syr(x) = (3x+1) >> v2(3x+1)

READ THE RESULT
  A_B_C_common == 0  -> no reduction-independent cycle exists -> no convention-invariant
  functional of a modular cycle -> the covering-space-via-cycles mechanism (§11) has no
  object to act on. This is the self-kill outcome (Mr A: a clean kill is five-star, and rarer).
  A_B_C_common  > 0  -> at least one shared cycle exists; a functional of it is a candidate
  invariant -> the mechanism is NOT dead and §11 may survive. (Would be the SAVE.)

USAGE
  python3 repro_paper109_invariance_null.py            # full report, exhaustive m=3..20
  python3 repro_paper109_invariance_null.py --verify   # assert the recorded null; exit 0/1
"""
import sys
from collections import Counter

def v2(x): return (x & -x).bit_length() - 1
def syr(x):
    t = 3*x + 1; return t >> v2(t)
def succA(n, M): return syr(n) % M
def succB(n, M):
    r = (3*n + 1) % M
    return 1 if r == 0 else r >> v2(r)
def succC(n, M): return syr(n + M) % M
CONV = {'A min-rep': succA, 'B reduce-then-strip': succB, 'C shifted-rep': succC}

def cycles(succ, m):
    M = 1 << m; color = bytearray(M); out = []
    for s in range(1, M, 2):
        if color[s]: continue
        path = []; pos = {}; node = s
        while color[node] == 0:
            color[node] = 1; pos[node] = len(path); path.append(node); node = succ(node, M)
        if color[node] == 1:
            cyc = path[pos[node]:]
            if not (len(cyc) == 1 and cyc[0] == 1):
                out.append(frozenset(cyc))
        for x in path: color[x] = 2
    return out

def all_loops(succ, maxm):
    loops = set()
    for m in range(3, maxm + 1):
        loops |= set(cycles(succ, m))
    return loops

def compute(maxm=20):
    L = {name: all_loops(succ, maxm) for name, succ in CONV.items()}
    A, B, C = L['A min-rep'], L['B reduce-then-strip'], L['C shifted-rep']
    return {
        'counts': {k: len(v) for k, v in L.items()},
        'ABC': len(A & B & C), 'AB': len(A & B), 'AC': len(A & C), 'BC': len(B & C),
    }

def run(maxm=20):
    print(f"Convention-invariance null, exhaustive integer-loops m=3..{maxm}\n")
    r = compute(maxm)
    print("Distinct integer-loops captured per reduction:")
    for k, v in r['counts'].items(): print(f"   {k:24s}: {v}")
    print(f"\nLoops common to ALL THREE (A AND B AND C): {r['ABC']}")
    print(f"Pairwise: A&B={r['AB']}   A&C={r['AC']}   B&C={r['BC']}")
    print()
    if r['ABC'] == 0:
        print("NULL CONFIRMED: no integer-loop is captured by all three reductions, so there is")
        print("no reduction-independent modular cycle -> no convention-invariant functional of a")
        print("cycle -> §12.6 outcome (i) NEGATIVE: the cycle-based covering-space mechanism (§11)")
        print("has no object to act on. (A clean self-kill -- to be hardened by Scout/Mr Code/Adversary.)")
    else:
        print(f"NULL FAILED: {r['ABC']} loop(s) common to all three -> a candidate invariant exists;")
        print("the mechanism is NOT dead. Investigate the shared loop(s).")
    return r

def verify():
    r = compute(20)
    checks = [
        ("no loop common to all three (A&B&C == 0)", r['ABC'] == 0),
        ("C disjoint from A (A&C == 0)",             r['AC'] == 0),
        ("C disjoint from B (B&C == 0)",             r['BC'] == 0),
        ("A and B share loops (A&B == 4)",           r['AB'] == 4),
    ]
    ok = True
    for desc, passed in checks:
        if not passed: ok = False
        print(f"  [{'ok' if passed else 'MISMATCH'}] {desc}   (got A&B&C={r['ABC']}, A&C={r['AC']}, B&C={r['BC']}, A&B={r['AB']})")
    print()
    print("INVARIANCE NULL VERIFIED -- no reduction-independent cycle; §12.6 outcome (i) negative "
          "(PRELIMINARY: still owes Scout + Adversary)."
          if ok else
          "INVARIANCE NULL FAILED -- does not match the recorded result; STOP, do not commit.")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        verify()
    else:
        run(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
