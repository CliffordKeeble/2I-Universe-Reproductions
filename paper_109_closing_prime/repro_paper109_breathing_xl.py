#!/usr/bin/env python3
"""
Reproduction R109-2 : breathing pattern of Paper 109, pushed to high m.
2I Universe Programme -- Reproductions repo.  Author: CinC (chat).  Verify/run: Mr Code.

PURPOSE
  R109-1 confirmed the breathing pattern EXHAUSTIVELY to m=24 by per-path walking.
  R109-2 pushes the exhaustive boundary toward m=30 -- the level where a THIRD breath
  would first appear (3*(D+M_D)=30) -- and as high as hardware allows, to either
  close or surface "no third breath".

OPERATOR DEFINITION (identical to R109-1; this is the operative convention)
  States : odd integers in [1, 2^m), indexed by i with n = 2i+1, i in [0, 2^{m-1}).
  Map    : n |-> Syr(n) mod 2^m,  Syr(n) = (3n+1) >> v2(3n+1).
  Deterministic out-degree 1 => functional graph => cycles well defined.

ALGORITHM (different from R109-1, so agreement is an independent cross-check)
  Cycle nodes of a functional graph = the nodes that survive iterated removal of
  in-degree-0 nodes (Kahn topological peeling). Fully vectorized; O(N) total.
  succ[] built in chunks (memory-flat). Then peel; survivors (minus the n=1 fixed
  point) are exactly the non-trivial cycle elements; group them by walking succ.

MEMORY (N = 2^{m-1} states): succ uint32 = 4N B, indeg int32 = 4N B, removed = N B.
  m=30 -> N=5.4e8 -> ~9 GB ;  m=31 -> ~18 GB ;  m=32 -> ~20 GB.  (32 GB box: 30/31 safe,
  32 tight -- watch RSS.)  Chunked build keeps transients ~0.5 GB.

USAGE
  python3 repro_paper109_breathing_xl.py [MAXM=32] [MINM=3]
  python3 repro_paper109_breathing_xl.py --verify     # cross-check m<=24 vs R109-1
"""
import sys, time, numpy as np

U64 = np.uint64
CHUNK = 1 << 26          # 64M states per build-chunk -> ~0.5 GB transient

def build_succ(m):
    """succ[i] = odd-index of (Syr(2i+1) mod 2^m), built in memory-flat chunks."""
    N = 1 << (m - 1)
    M = U64(1) << U64(m)
    modmask = M - U64(1)
    succ = np.empty(N, dtype=np.uint32)
    one = U64(1)
    for start in range(0, N, CHUNK):
        end = min(start + CHUNK, N)
        res = np.arange(start, end, dtype=U64)      # i
        res *= U64(6); res += U64(4)                # 3n+1 = 6i+4 (even)
        even = (res & one) == 0                     # strip all factors of 2
        guard = 0
        while even.any():
            np.right_shift(res, one, out=res, where=even)
            even = (res & one) == 0
            guard += 1
            if guard > m + 4:                       # safety; v2(6i+4) <= ~m
                break
        res &= modmask                              # mod 2^m
        res >>= one                                 # odd -> odd-index ((n-1)/2 = n>>1)
        succ[start:end] = res.astype(np.uint32)
        del res, even
    return succ, N

def cycle_nodes(succ, N):
    """Indices on a cycle = survivors of Kahn peeling (in-degree-0 removal)."""
    if N <= (1 << 30):                              # bincount path (fast), int64 transient ok
        indeg = np.bincount(succ, minlength=N).astype(np.int32)
    else:                                           # memory-safe path for very large N
        indeg = np.zeros(N, dtype=np.int32)
        np.add.at(indeg, succ, 1)
    removed = np.zeros(N, dtype=bool)
    frontier = np.flatnonzero(indeg == 0)
    while frontier.size:
        removed[frontier] = True
        tgt = succ[frontier]
        np.add.at(indeg, tgt, -1)                   # decrement successors' in-degree
        cand = tgt[(indeg[tgt] == 0) & (~removed[tgt])]
        frontier = np.unique(cand)
    return np.flatnonzero(~removed)                 # includes the n=1 self-loop (index 0)

def extract_cycles(succ, nodes):
    """Group surviving nodes into individual cycles by walking succ."""
    nodeset = set(int(x) for x in nodes)
    seen, cycles = set(), []
    for s in nodeset:
        if s in seen: continue
        cyc, x = [], s
        while x not in seen:
            seen.add(x); cyc.append(x); x = int(succ[x])
        if x == s:                                  # closed a fresh cycle
            cycles.append(cyc)
    return cycles

def sclass(n): return {1:'S1',3:'S3',5:'S5',7:'S7'}[n % 8]

def analyse(m):
    t0 = time.time()
    succ, N = build_succ(m)
    nodes = cycle_nodes(succ, N)
    cycles = [c for c in extract_cycles(succ, nodes) if not (len(c) == 1 and c[0] == 0)]
    info = []
    for c in cycles:
        ns = [2*i + 1 for i in c]
        info.append((len(c), min(ns), ns))
    info.sort()
    del succ, nodes
    return info, N, time.time() - t0

def detail(length, ns):
    from collections import Counter
    mod24 = sorted(set(x % 24 for x in ns))
    missing = [v for v in [1,3,5,7,9,11,13,15,17,19,21,23] if v not in mod24]
    sc = dict(sorted(Counter(sclass(x) for x in ns).items()))
    return f"len={length} min={min(ns)} verts={len(mod24)}/12 missing={missing} Sclass={sc}"

def run(minm, maxm):
    print(f"R109-2  exhaustive breathing search, m={minm}..{maxm}")
    print("m  | states 2^(m-1) | #cyc | lengths            | time(s)")
    print("---+----------------+------+--------------------+--------")
    breath3 = []
    for m in range(minm, maxm + 1):
        info, N, dt = analyse(m)
        lens = [x[0] for x in info]
        print(f"{m:2d} | {N:>14,} | {len(info):^4d} | {str(lens):<18} | {dt:7.1f}")
        for length, mn, ns in info:
            print(f"      -> {detail(length, ns)}")
            if m >= 25:
                breath3.append((m, length, mn))
    print()
    if breath3:
        print(f"*** POTENTIAL THIRD BREATH: cycles found at m>=25: {breath3}")
    else:
        print(f"*** NO THIRD BREATH found exhaustively through m={maxm} "
              f"(in particular m=30 if reached).")

def verify():
    """Cross-check against R109-1's exhaustively-known values (m<=24)."""
    expect = {**{m: [] for m in list(range(3,10))+list(range(13,20))+list(range(21,25))},
              10:[26], 11:[25], 12:[6,7], 20:[22]}
    ok = True
    for m in range(3, 25):
        info, N, dt = analyse(m)
        lens = sorted(x[0] for x in info)
        exp = expect[m]
        flag = "ok" if lens == exp else "MISMATCH"
        if lens != exp: ok = False
        print(f"m={m:2d}: lengths={lens} expected={exp} [{flag}]  ({dt:.1f}s)")
    print("\n" + ("R109-2 CROSS-CHECK VERIFIED -- peeling algorithm agrees with R109-1 to m=24."
                  if ok else "R109-2 CROSS-CHECK FAILED."))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        verify()
    else:
        maxm = int(sys.argv[1]) if len(sys.argv) > 1 else 32
        minm = int(sys.argv[2]) if len(sys.argv) > 2 else 3
        run(minm, maxm)
