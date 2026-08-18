"""PHASE A / T4.  WHAT IS THE GENERAL FORM OF LEGIBILITY?

W-46 measured legibility = 2n - 1 = perimeter/2 - 1 for n x n blocks. That is the one result Phase C
depends on and the one Phase A could not prove. Before attempting a proof, find the actual general
form: measure it across shapes, aspect ratios and complexes.

WHAT LEGIBILITY IS, REDUCED. A region A's records are subsets of A's plaquettes. An outside observer
holds gauge-invariant operators built from OUTSIDE plaquettes. The pairing is
    <T,S> = |bd(T) INTERSECT bd(S)| mod 2
and for lattice plaquettes two plaquettes share exactly one link iff they are ADJACENT. So the
pairing matrix is the BIPARTITE ADJACENCY between outside and inside plaquettes, and
    legibility = rank over GF(2) of that matrix.
Its kernel is the set of inside region-combinations whose boundary is INTERIOR to A -- invisible from
outside. So legibility = |A| - dim{S : bd(S) touches no link shared with the outside}.
"""
import itertools, numpy as np
def rank_gf2(rows):
    piv=[]; r=0
    for v in rows:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r

def grid(NX,NY):
    """plaquette (i,j) for 0<=i<NX, 0<=j<NY on a lattice; returns id map and link sets"""
    def hid(i,j): return j*NX+i
    def vx(i,j):  return NX*(NY+1)+ j*(NX+1)+i
    PL={}
    for j in range(NY):
        for i in range(NX):
            PL[(i,j)]=[hid(i,j),hid(i,j+1),vx(i,j),vx(i+1,j)]
    return PL

def legibility(shape,NX,NY):
    PL=grid(NX,NY)
    A=[c for c in shape]; out=[c for c in PL if c not in set(A)]
    idx={c:k for k,c in enumerate(A)}
    rows=[]
    for q in out:
        bq=set(PL[q]); v=0
        for c in A:
            if len(bq & set(PL[c]))%2: v|=(1<<idx[c])
        rows.append(v)
    return rank_gf2(rows), len(A)

def perim(shape):
    S=set(shape); p=0
    for (i,j) in S:
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            if (i+d[0],j+d[1]) not in S: p+=1
    return p

print("PHASE A / T4.  legibility = rank over GF(2) of the outside-inside plaquette adjacency")
print(f"\n  {'shape':>22s} {'|A|':>4s} {'perimeter':>10s} {'legibility':>11s} {'perim/2-1':>10s} {'|A|-1':>6s}")
print("  "+"-"*72)
def blk(nx,ny,ox=2,oy=2): return [(ox+i,oy+j) for i in range(nx) for j in range(ny)]
cases=[]
for n in (1,2,3,4,5): cases.append((f"square {n}x{n}", blk(n,n), n+6,n+6))
for (a,b) in ((1,4),(2,4),(2,6),(3,5),(1,6)): cases.append((f"rect {a}x{b}", blk(a,b), a+6,b+6))
# L-shape and a cross
L=[(2,2),(3,2),(4,2),(2,3),(2,4)]
cases.append(("L-shape (5 cells)",L,9,9))
X=[(3,2),(2,3),(3,3),(4,3),(3,4)]
cases.append(("plus (5 cells)",X,9,9))
ring=[(i,j) for i in range(2,5) for j in range(2,5) if (i,j)!=(3,3)]
cases.append(("3x3 minus centre",ring,9,9))
for name,shape,NX,NY in cases:
    lg,na=legibility(shape,NX,NY); pe=perim(shape)
    print(f"  {name:>22s} {na:4d} {pe:10d} {lg:11d} {pe//2-1:10d} {na-1:6d}")
print()
print("  WHICH FORMULA SURVIVES?")
ok_p=ok_a=True
for name,shape,NX,NY in cases:
    lg,na=legibility(shape,NX,NY); pe=perim(shape)
    if lg!=pe//2-1: ok_p=False
    if lg!=na-1:    ok_a=False
print(f"    legibility = perimeter/2 - 1 holds for every shape tested : {ok_p}")
print(f"    legibility = |A| - 1        holds for every shape tested : {ok_a}")
print()
print("  THE KERNEL -- what is invisible from outside?  dim ker = |A| - legibility")
for name,shape,NX,NY in cases:
    lg,na=legibility(shape,NX,NY)
    print(f"    {name:>22s}  |A|={na:3d}  legible={lg:3d}  invisible={na-lg:3d}")

print()
print("  THE PATTERN, AND THE DIAGNOSIS OF W-46'S NUMBER")
print(f"    {'n':>3s} {'|A|=n^2':>8s} {'legibility':>11s} {'invisible':>10s} {'(n-2)^2':>8s} {'4(n-1)':>7s}")
for n in (1,2,3,4,5):
    lg,na=legibility(blk(n,n),n+6,n+6)
    print(f"    {n:3d} {na:8d} {lg:11d} {na-lg:10d} {max(0,(n-2))**2:8d} {4*(n-1):7d}")
print()
print("    invisible = (n-2)^2 = the plaquettes of A that touch NO outside plaquette.")
print("    legibility = 4(n-1), i.e. the boundary-touching plaquettes: LINEAR in n, not quadratic.")
print()
print("  W-46 PLACED THE REGION AT THE LATTICE CORNER. Re-run its geometry to confirm.")
def corner_case(n,LN):
    PL=grid(LN,LN)
    A=[(i,j) for i in range(n) for j in range(n)]        # bottom-left CORNER, as in W-46
    return legibility(A,LN,LN)
for n,LN in ((2,4),(3,5),(4,7)):
    lg,na=corner_case(n,LN)
    lgc,_=legibility(blk(n,n),n+6,n+6)
    print(f"    n={n}: corner-placed legibility = {lg:3d}  (W-46 reported {2*n-1});"
          f"   interior-placed = {lgc:3d}")
