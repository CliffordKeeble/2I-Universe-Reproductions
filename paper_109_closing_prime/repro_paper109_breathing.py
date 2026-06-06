#!/usr/bin/env python3
"""
Reproduction R109-1 : the 'breathing pattern' of Paper 109 (The Closing Prime).
2I Universe Programme -- Reproductions repo.  Author: CinC (chat).  Verify: Mr Code.

WHAT THIS TESTS
  Paper 109 Table 5 claims the mod-2^m Syracuse transfer operator has non-trivial
  cycles ONLY at m = 10 (len 26), 11 (len 25), 12 (two: 6 & 7) and 20 (len 22),
  and is otherwise nilpotent (only the fixed point n=1) for 3 <= m <= 32; with the
  m=10 cycle visiting 8/12 vertices mod 24 (missing {3,9,15,21}), S-class split
  S1:7 S3:7 S5:2 S7:10, minimal element 47, and 47's real trajectory reaching 1
  in 38 steps (peak 3077).

OPERATOR DEFINITION  (explicit -- this is what Paper 109 v2.3 left to this listing)
  States : odd integers in [1, 2^m)            (count 2^{m-1})
  Map    : n  |->  Syr(n) mod 2^m              ("minimal-representative Syracuse")
           Syr(n) = (3n+1) >> v2(3n+1)          (odd-to-odd Syracuse step)
  Deterministic (out-degree 1) => cycles are well defined. T does NOT descend to
  Z/2^m (3 != 27 under T); this convention fixes the representative AFTER applying
  the integer map. A different convention may give a different table -- so the
  convention IS part of the result.

SEARCH  : EXHAUSTIVE for every m reported (every odd state visited; no seeding).
USAGE   : python3 repro_paper109_breathing.py [MAXM]      (default 24)
          python3 repro_paper109_breathing.py --verify    (assert claimed values)
"""
import sys
from collections import Counter

def v2(x): return (x & -x).bit_length() - 1
def syr(n):
    m = 3*n + 1
    return m >> v2(m)

def find_cycles(m):
    """All non-trivial cycles of n |-> syr(n) mod 2^m on odd states in [1,2^m)."""
    M = 1 << m
    color = bytearray(M)                 # 0 unseen / 1 on-path / 2 done, indexed by n
    cycles = []
    for s in range(1, M, 2):
        if color[s]:
            continue
        path, pos, node = [], {}, s
        while color[node] == 0:
            color[node] = 1
            pos[node] = len(path)
            path.append(node)
            node = syr(node) % M
        if color[node] == 1:
            cyc = path[pos[node]:]
            if not (len(cyc) == 1 and cyc[0] == 1):
                cycles.append(cyc)
        for x in path:
            color[x] = 2
    return cycles

def sclass(n): return {1:'S1',3:'S3',5:'S5',7:'S7'}[n % 8]

def real_traj(n):
    steps, peak, seq = 0, n, [n]
    while n != 1:
        n = syr(n); steps += 1; peak = max(peak, n); seq.append(n)
    return steps, peak, seq

def run(maxm):
    table = {}
    print("m  | #nontrivial cycles | lengths")
    print("---+--------------------+--------")
    for m in range(3, maxm + 1):
        cs = find_cycles(m)
        table[m] = cs
        print(f"{m:2d} | {len(cs):^18d} | {sorted(len(c) for c in cs)}")
    for m in (10, 20):
        for i, c in enumerate(table.get(m, [])):
            mod24 = sorted(set(x % 24 for x in c))
            missing = [v for v in [1,3,5,7,9,11,13,15,17,19,21,23] if v not in mod24]
            sc = dict(sorted(Counter(sclass(x) for x in c).items()))
            print(f"\n  m={m} cycle #{i}: length={len(c)}  min={min(c)}")
            print(f"    mod-24 vertices visited ({len(mod24)}/12), missing {missing}")
            print(f"    S-class composition: {sc}")
    st, pk, _ = real_traj(47)
    print(f"\n  real Syracuse trajectory of 47: {st} steps, peak {pk}")
    return table

def verify():
    t = run(24)
    def lens(m): return sorted(len(c) for c in t.get(m, []))
    ok = True
    expect = {**{m: [] for m in list(range(3,10))+list(range(13,20))+list(range(21,25))},
              10:[26], 11:[25], 12:[6,7], 20:[22]}
    for m, e in expect.items():
        if lens(m) != e:
            print(f"  MISMATCH m={m}: got {lens(m)}, expected {e}"); ok = False
    c10 = t[10][0]
    if sorted(set(x%24 for x in c10)) != [1,5,7,11,13,17,19,23]: ok=False; print("  MISMATCH m=10 vertices")
    if dict(Counter(sclass(x) for x in c10)) != {'S1':7,'S3':7,'S5':2,'S7':10}: ok=False; print("  MISMATCH m=10 S-class")
    if min(c10) != 47: ok=False; print("  MISMATCH m=10 min element")
    st, pk, _ = real_traj(47)
    if (st, pk) != (38, 3077): ok=False; print(f"  MISMATCH 47 trajectory: {st},{pk}")
    print("\n" + ("REPRODUCTION VERIFIED -- all Paper 109 Table 5 claims reproduce (exhaustive to m=24)."
                  if ok else "REPRODUCTION FAILED -- see mismatches above."))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        verify()
    else:
        run(int(sys.argv[1]) if len(sys.argv) > 1 else 24)
