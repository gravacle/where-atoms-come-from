"""W-54.  ACCUMULATION: DOES THE WHOLE DETERMINE ANYTHING ITS PARTS DO NOT?

Road item 20. Every deflation in this program was obtained at a SINGLE scale -- one region, one bath,
one constraint -- so none could have detected a collective effect. This asks the collective question.

IT MUST FIRST BE MADE WELL-POSED. On a FIXED finite carrier every quantity is an exact function of
the incidence matrix, so "emergence" in the strong sense (new laws at large scale) cannot arise: that
needs a limit, and these carriers are small and exact. The sharpest version that CAN fail is:

    GLOBAL capacity is m - 1 (T3). Partition the carrier into k parts. How much of that capacity is
    reachable using ONLY records supported inside a single part -- and how much requires records that
    SPAN parts and belong to no part?

If all of it is reachable locally, the whole is exactly its parts and accumulation adds nothing.
If not, then some of a carrier's record capacity lives in records NO PART CAN HOST, and the deficit
is a collective property. That is decidable, and it is not circular: nothing here asks to recover the
incidence, only to compare two capacities computed from it.

PREDICTION STATED IN ADVANCE (so it can be wrong): local capacity should be sum_i (m_i - 1) = m - k,
so the deficit should be exactly k - 1, independent of HOW the carrier is partitioned. If the deficit
depends on the partition's shape, the story is richer than the counting.
"""
import itertools, numpy as np
def rank_gf2(vecs):
    piv=[]; r=0
    for v in vecs:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r
def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)

def cap_restricted(L,PL,allowed):
    """largest independent protectable set drawn ONLY from `allowed` region-subsets"""
    best=0
    for lk in range(L):
        free=[S for S in allowed if lk not in bd(S,PL)]
        best=max(best, rank_gf2([sum(1<<i for i in S) for S in free]))
    return best

def analyse(name,L,PL,parts):
    m=len(PL)
    ALL=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
    glob=cap_restricted(L,PL,ALL)
    LOCAL=[S for S in ALL if any(set(S)<=set(pt) for pt in parts)]
    loc=cap_restricted(L,PL,LOCAL)
    k=len(parts)
    print(f"  {name:>30s} m={m:2d} k={k:2d}  global {glob:3d}  local-only {loc:3d}  "
          f"deficit {glob-loc:3d}   predicted k-1 = {k-1:2d}   {'ok' if glob-loc==k-1 else 'DIFFERS'}")
    return glob,loc,k

def patch(nx,ny):
    vid={(i,j):j*nx+i for j in range(ny) for i in range(nx)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append(('h',i,j))
    for j in range(ny-1):
        for i in range(nx): E.append(('v',i,j))
    ind={e:k for k,e in enumerate(E)}
    PL=[]
    for j in range(ny-1):
        for i in range(nx-1):
            PL.append({ind[('h',i,j)],ind[('v',i+1,j)],ind[('h',i,j+1)],ind[('v',i,j)]})
    return len(E),PL,(nx-1),(ny-1)

print("W-54  ACCUMULATION: is the whole's capacity reachable from inside its parts?")
print()
L,PL,MX,MY=patch(4,3)          # 3x2 = 6 plaquettes
def cells(MX,MY): return [(i,j) for j in range(MY) for i in range(MX)]
def pid(i,j,MX): return j*MX+i
C=cells(MX,MY)
parts_sets={
 "6 singletons":            [[pid(i,j,MX)] for (i,j) in C],
 "3 vertical pairs":        [[pid(i,0,MX),pid(i,1,MX)] for i in range(MX)],
 "2 horizontal rows":       [[pid(i,j,MX) for i in range(MX)] for j in range(MY)],
 "one part (the whole)":    [[pid(i,j,MX) for (i,j) in C]],
 "uneven: 4 + 1 + 1":       [[0,1,2,3],[4],[5]],
}
for nm,pt in parts_sets.items(): analyse(nm,L,PL,pt)

print()
L2,PL2,MX2,MY2=patch(4,4)      # 3x3 = 9 plaquettes
C2=cells(MX2,MY2)
parts2={
 "9 singletons":            [[pid(i,j,MX2)] for (i,j) in C2],
 "3 rows":                  [[pid(i,j,MX2) for i in range(MX2)] for j in range(MY2)],
 "3 columns":               [[pid(i,j,MX2) for j in range(MY2)] for i in range(MX2)],
 "one part (the whole)":    [[pid(i,j,MX2) for (i,j) in C2]],
 "uneven: 5 + 3 + 1":       [[0,1,2,3,4],[5,6,7],[8]],
}
for nm,pt in parts2.items(): analyse(nm,L2,PL2,pt)

print()
print("  READING")
print("    deficit = k - 1 everywhere means: EVERY partition into k parts leaves exactly k-1 of the")
print("    carrier's records unreachable from inside any part. Those records SPAN parts and belong to")
print("    none. The finer the partition, the larger the share of capacity that is collective:")
for nm,pt in parts2.items():
    g,l,k=analyse("   "+nm,L2,PL2,pt)
    print(f"        {nm:>24s}:  {(g-l)}/{g} of the capacity is collective "
          f"({100*(g-l)/max(g,1):.0f}%)")
