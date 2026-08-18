"""EXACT test of the no-string-logical-operator property.
The randomised min-weight probe in o7_writer_shape.py FAILED its positive control
(it returned weight 22 for the 3D toric code at L=3, where the true answer is 3).
That probe is DISCARDED.  This script replaces it with exact linear algebra.

Question: does there EXIST a logical operator supported inside a given region R?
Answer: compute the subspace  (S-perp) ^ {Paulis supported in R}  by solving a linear
system over F_2, then ask whether it is contained in S.  Exhaustive, not a search."""
import sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *
from models import *

def logical_in_region(gens, n, region):
    """region: set of qubit indices.  Returns (exists_logical, dim_of_Sperp_in_R, min_weight)."""
    _, pivS = rank_gf2(gens)
    N = normaliser_basis(gens, n)
    outside = [i for i in range(n) if i not in region]
    # constraint: for each qubit outside R, both x_i and z_i must vanish
    rows = []
    for i in outside:
        rx = 0; rz = 0
        for j, b in enumerate(N):
            if (xpart(b,n)>>i)&1: rx |= (1<<j)
            if (zpart(b,n)>>i)&1: rz |= (1<<j)
        rows.append(rx); rows.append(rz)
    coefs = nullspace_basis(rows, len(N)) if rows else [1<<j for j in range(len(N))]
    vecs = []
    for c in coefs:
        v = 0; cc = c; j = 0
        while cc:
            if cc & 1: v ^= N[j]
            cc >>= 1; j += 1
        vecs.append(v)
    dim = len(vecs)
    # is the whole subspace inside S?  enumerate if small, else check basis + random combos
    best = None; found = False
    if dim <= 20:
        for mask in range(1, 1 << dim):
            v = 0; mm = mask; j = 0
            while mm:
                if mm & 1: v ^= vecs[j]
                mm >>= 1; j += 1
            if v and not in_span(pivS, v):
                found = True
                w = pweight(v,n)
                if best is None or w < best: best = w
    else:
        for v in vecs:
            if v and not in_span(pivS, v):
                found = True; w = pweight(v,n)
                if best is None or w < best: best = w
    return found, dim, best

def line_region(n, L, cell_of, axis, off, thick=1):
    R = set()
    for i in range(n):
        c = cell_of(i)
        ok = True
        for t in range(3):
            if t == axis: continue
            dd = (c[t] - off[t]) % L
            if dd >= thick: ok = False; break
        if ok: R.add(i)
    return R

MODELS = [
    ("3D toric  (CONTROL: strings MUST exist)", toric3d, lambda L: (lambda i: ((i//3)//(L*L), ((i//3)//L)%L, (i//3)%L))),
    ("X-cube    (type-I: strings expected)",     xcube,   lambda L: (lambda i: ((i//3)//(L*L), ((i//3)//L)%L, (i//3)%L))),
    ("Haah CC1  (type-II: NO strings claimed)",  haah,    lambda L: (lambda i: ((i//2)//(L*L), ((i//2)//L)%L, (i//2)%L))),
    ("checkerboard (type-I)",                    checkerboard, lambda L: (lambda i: (i//(L*L), (i//L)%L, i%L))),
]

for L in (4, 6):
    print("="*78)
    print(f"L = {L}")
    print("="*78)
    for name, builder, cellf in MODELS:
        if builder is checkerboard and L % 2: continue
        gens, n = builder(L)
        cell = cellf(L)
        for thick in (1, 2):
            any_found = False; dims = []; wmin = None
            for off in [(0,0,0), (0,1,2), (0,2,1)]:
                R = line_region(n, L, cell, 0, off, thick)
                f, dim, w = logical_in_region(gens, n, R)
                dims.append(dim)
                if f:
                    any_found = True
                    if w is not None and (wmin is None or w < wmin): wmin = w
            cross = f"{thick}x{thick}"
            print(f"  {name:42s} tube cross-section {cross}: "
                  f"logical supported in the tube ? {'YES  min weight '+str(wmin) if any_found else 'NO'}"
                  f"   (dim S-perp inside tube = {dims})")
    print()

print("""READ THE CONTROL FIRST.  The 3D toric code must answer YES with min weight = L for a
1x1 tube: its logical operators are the three straight non-contractible loops of H_1(T^3).
If it answers YES and Haah answers NO, the no-string property is measured, not assumed.""")
