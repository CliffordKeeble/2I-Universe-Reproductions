#!/usr/bin/env python3
"""
Reproduction R109-3 : CONVENTION-ROBUSTNESS NULL for Paper 109's breathing pattern.
2I Universe Programme -- Reproductions repo.  Author: CinC (chat).  Verify/run: Mr Code.

THE QUESTION (Mr A's fifth-star gate #1)
  T(n) mod 2^m is not well-defined; the breathing pattern (cycles at m=10, 20 = D+M_D,
  2(D+M_D)) is a property of ONE chosen resolution -- the minimal-representative Syracuse
  map. R109-1/R109-2 reproduced THAT map two ways (implementation-correctness). They do
  NOT show the pattern survives a different, equally-natural reduction. If another
  reduction breathes at different m, or everywhere, or nowhere, the 10/20/30 = k*(D+M_D)
  reading is a convention artifact, not structure. This runs the null: same 3n+1 Syracuse
  dynamics, three different natural reductions, compared.

THREE REDUCTIONS (all self-maps of odd residues in [1,2^m); all "equally natural")
  A  minimal-representative : n -> Syr(n) mod 2^m         [strip 2s of the full integer, then reduce]
  B  reduce-then-strip      : r=(3n+1) mod 2^m; n -> oddpart(r)   [stay inside the ring throughout]
  C  shifted-representative : n -> Syr(n + 2^m) mod 2^m   [same residue class, different lift]
  where Syr(x) = (3x+1) >> v2(3x+1) is the odd-to-odd Syracuse step.
  A and B and C agree on the residue class of n but differ on which integer's 2-adic
  valuation is used -- exactly the ambiguity Paper 109 v2.5 concedes in section 9.

READ THE RESULT
  If A, B, C all breathe at m=10 and m=20 (and nowhere else) -> pattern is convention-
  ROBUST; the D+M_D reading earns its keep. If they breathe at different m -> the pattern
  is convention-DEPENDENT and the numerology is an artifact of the chosen lift.

USAGE
  python3 repro_paper109_convention_null.py [MAXM]   exhaustive m=3..MAXM (default 20), all 3 reductions
  python3 repro_paper109_convention_null.py --verify assert the recorded null result; exit 0 (pass) / 1 (stop)
"""
import sys
def v2(x): return (x & -x).bit_length() - 1
def syr(n):
    t = 3*n + 1; return t >> v2(t)

def succA(n, M): return syr(n) % M
def succB(n, M):
    r = (3*n + 1) % M
    if r == 0: return 1
    return r >> v2(r)
def succC(n, M): return syr(n + M) % M

CONV = {'A min-rep': succA, 'B reduce-then-strip': succB, 'C shifted-rep': succC}

def find_cycles(succ, m):
    M = 1 << m
    color = bytearray(M)
    cycles = []
    for s in range(1, M, 2):
        if color[s]: continue
        path, pos, node = [], {}, s
        while color[node] == 0:
            color[node] = 1; pos[node] = len(path); path.append(node)
            node = succ(node, M)
        if color[node] == 1:
            cyc = path[pos[node]:]
            if not (len(cyc) == 1 and cyc[0] == 1):
                cycles.append(cyc)
        for x in path: color[x] = 2
    return cycles

def run(maxm, quiet=False):
    if not quiet:
        print(f"Convention-robustness null, 3n+1 Syracuse, exhaustive m=3..{maxm}\n")
    results = {}
    for name, succ in CONV.items():
        if not quiet:
            print(f"=== Reduction {name} ===")
            print("m  | #cyc | lengths")
        breaths = []
        for m in range(3, maxm+1):
            cs = find_cycles(succ, m)
            lens = sorted(len(c) for c in cs)
            if cs: breaths.append(m)
            if not quiet:
                print(f"{m:2d} | {len(cs):^4d} | {lens}")
        results[name] = breaths
        if not quiet:
            print(f"  -> breathes at m = {breaths}\n")
    A = results['A min-rep']
    robust = all(results[n] == A for n in results)
    if not quiet:
        print("="*52)
        print("SUMMARY: m-levels with non-trivial cycles, per reduction")
        for name, b in results.items():
            print(f"  {name:24s}: {b}")
        print()
        print("CONVENTION-ROBUST: all three reductions breathe at the same levels."
              if robust else
              "CONVENTION-DEPENDENT: the breathing levels DIFFER across reductions "
              "-> the D+M_D reading is an artifact of the chosen lift.")
    return results, robust

def verify():
    """Independent-run gate: assert the convention-dependence reproduces exactly (exhaustive m=3..20).
       NB this is a NULL: the gate confirms the reductions DISAGREE. A 'robust' (agreeing) result
       would mean the published withdrawal was wrong -- so robust => FAIL => STOP."""
    results, robust = run(20, quiet=True)
    expect = {
        'A min-rep':           [10, 11, 12, 20],
        'B reduce-then-strip': [9, 11, 12, 13],
        'C shifted-rep':       list(range(3, 21)),
    }
    ok = True
    for name, exp in expect.items():
        got = results[name]
        flag = "ok" if got == exp else "MISMATCH"
        if got != exp: ok = False
        print(f"  {name:24s}: {got}   expected {exp}   [{flag}]")
    if robust:
        ok = False
        print("  MISMATCH: verdict came out CONVENTION-ROBUST (expected DEPENDENT)")
    print()
    print("R109-3 NULL VERIFIED -- breath levels differ across reductions; the D+M_D structural "
          "reading is convention-dependent and is correctly withdrawn in Paper 109."
          if ok else
          "R109-3 NULL FAILED -- does not match the recorded result; STOP, do not commit.")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        verify()
    else:
        run(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
