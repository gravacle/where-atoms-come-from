"""Stronger than the local-Clifford search: EXHAUSTIVELY enumerate EVERY CSS code on
n=5 qubits and show none has [[5,1,3]] parameters.  Since any Clifford circuit (of any
depth) preserves n, k and d, this rules out equivalence to a CSS code by ANY such circuit,
not merely by single-qubit Cliffords."""
import sys, itertools
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *

def subspaces(n):
    """all subspaces of F_2^n, as (dim, reduced basis tuple)"""
    vecs = list(range(1, 1<<n))
    out = {}
    for r in range(0, n+1):
        for combo in itertools.combinations(vecs, r):
            rk, _ = rank_gf2(list(combo))
            if rk != r: continue
            # canonical: span as frozenset
            span = set([0])
            for v in combo:
                span |= {s ^ v for s in span}
            key = frozenset(span)
            if key not in out: out[key] = (r, combo)
    return out

def css_params(n, CX, CZ):
    """CSS: X-stabilisers span CX, Z-stabilisers span CZ, need CX . CZ = 0.
       k = n - dim CX - dim CZ.  d = min weight of a nontrivial logical."""
    gens = [mk(v,0,n) for v in CX] + [mk(0,v,n) for v in CZ]
    if not gens: return None
    if all_commute(gens, n) != 0: return None
    k, r = code_k(gens, n)
    if k != 1: return None
    d, m = min_logical_weight(gens, n)
    return k, d

n = 5
print(f"enumerating all subspaces of F_2^{n} ...")
S = subspaces(n)
print(f"  {len(S)} subspaces")
best = 0; found = []
pairs = 0
sslist = [(r, combo) for (r, combo) in S.values()]
for (rx, bx) in sslist:
    for (rz, bz) in sslist:
        if rx + rz != n - 1: continue          # k = n - rx - rz = 1
        ok = True
        for a in bx:
            for b in bz:
                if popcount(a & b) & 1: ok = False; break
            if not ok: break
        if not ok: continue
        pairs += 1
        res = css_params(n, list(bx), list(bz))
        if res is None: continue
        k, d = res
        if d is not None and d > best: best = d
        if d is not None and d >= 3: found.append((bx, bz, d))
print(f"  admissible CSS pairs with k=1 on n=5: {pairs}")
print(f"  BEST DISTANCE ACHIEVED BY ANY CSS CODE WITH n=5, k=1 : d = {best}")
print(f"  any with d >= 3 ? {len(found)}")
print(f"\n  [[5,1,3]] exists as a STABILISER code (exact d=3 verified earlier).")
print(f"  [[5,1,3]] as a CSS code: {'EXISTS' if found else 'DOES NOT EXIST'}")
print("""
  POSITIVE CONTROL -- the SAME machinery on n=7 must FIND a CSS [[7,1,3]] (Steane).
  Restricted to the self-orthogonal family C_X = C_Z = C (which is where Steane lives),
  enumerated over all 3-dimensional subspaces of F_2^7.""")
n = 7
seen = set(); subs3 = []
V = list(range(1, 1<<7))
import random
for a in V:
    for b in V:
        if b <= a: continue
        for c in V:
            if c <= b: continue
            if rank_gf2([a,b,c])[0] != 3: continue
            span = set([0])
            for v in (a,b,c): span |= {t ^ v for t in span}
            key = frozenset(span)
            if key in seen: continue
            seen.add(key); subs3.append((a,b,c))
print(f"  distinct 3-dim subspaces of F_2^7 enumerated: {len(subs3)}  (Gaussian binomial [7 choose 3]_2 = 11811)")
hit = None
for (a,b,c) in subs3:
    basis = (a,b,c)
    ok = True
    for u in basis:
        for w in basis:
            if popcount(u & w) & 1: ok = False; break
        if not ok: break
    if not ok: continue
    gens = [mk(v,0,7) for v in basis] + [mk(0,v,7) for v in basis]
    if all_commute(gens,7): continue
    k,_ = code_k(gens,7)
    if k != 1: continue
    d,_ = min_logical_weight(gens,7)
    if d is not None and d >= 3: hit = (basis, d); break
print(f"  CSS [[7,1,3]] found ? {'YES  basis='+str(hit[0])+'  d='+str(hit[1]) if hit else 'NO'}")
print(f"  SELF-CHECK (the enumerator CAN return d=3 when d=3 exists): {'PASS' if hit else 'FAIL'}")
print("""
  CONCLUSION.  The n=5 sweep is exhaustive over EVERY CSS code on 5 qubits, not merely
  over local-Clifford images of one code.  Since any Clifford circuit of any depth
  preserves [[n,k,d]], no circuit can carry [[5,1,3]] to a CSS code.  The perfect code
  is therefore OUTSIDE the length-2 F_2 chain-complex class, unconditionally.""")
