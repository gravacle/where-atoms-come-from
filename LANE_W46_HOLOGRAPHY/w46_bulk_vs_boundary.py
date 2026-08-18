"""W-46.  IS THE CAPACITY LAW HOLOGRAPHIC?  BULK STORAGE vs EXTERNAL LEGIBILITY.

The principal proposes that gravity and capacity converge exactly at the Bekenstein/holographic
bound: the maximum information in a region is set by its BOUNDARY AREA. W-42 measured
capacity = area - 1 where "area" is the number of PLAQUETTES -- the EXTENSIVE measure of a 2D
region, its bulk. That is a VOLUME law, the opposite of holographic, and W-42 recorded
cap/area -> 1 while cap/perimeter grows without bound.

BEFORE TREATING THAT AS A DISAGREEMENT, SEPARATE TWO DIFFERENT QUANTITIES THE PROGRAM HAS BEEN
RUNNING TOGETHER:
   STORAGE    -- how many independent records the region can HOLD.                (W-42: bulk)
   LEGIBILITY -- how many of them an observer OUTSIDE the region can DISTINGUISH.  (never measured)
Holography is a claim about what is READABLE FROM OUTSIDE, not about what fits inside. And this
program has already shown these come apart violently: W-36 found a gauge record has NO redundancy
(I(|F|=1)/I(all) = 0.047) and W-37 that reading it requires transport around a CLOSED path.

SO MEASURE LEGIBILITY DIRECTLY. An outside observer holds the algebra of gauge-invariant operators
supported entirely OUTSIDE region A. Ask how many of A's internal records that algebra can separate.
Everything is linear over GF(2): a record is a plaquette subset, its boundary is a link set, and two
internal records are distinguishable from outside iff some outside-supported invariant sees them
differently. Computed by rank, not by simulation.
"""
import itertools, numpy as np

def lattice(nx,ny):
    V=[(i,j) for j in range(ny) for i in range(nx)]; vid={v:k for k,v in enumerate(V)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append((vid[(i,j)],vid[(i+1,j)]))
    H=len(E)
    for j in range(ny-1):
        for i in range(nx): E.append((vid[(i,j)],vid[(i,j+1)]))
    hid=lambda i,j: j*(nx-1)+i
    vx =lambda i,j: H + j*nx + i
    PL=[[hid(i,j),vx(i+1,j),hid(i,j+1),vx(i,j)] for j in range(ny-1) for i in range(nx-1)]
    return len(E),PL,(nx-1),(ny-1)

def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)

def rank_gf2(vecs):
    rows=[v for v in vecs if v]; r=0; piv=[]
    for v in rows:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r

print("W-46  STORAGE vs LEGIBILITY.  region A = a rectangular block of plaquettes inside a lattice.")
print(f"  {'lattice':>9s} {'A (px,py)':>10s} {'|A| plaq':>9s} {'perim(A)':>9s} "
      f"{'STORAGE cap':>12s} {'LEGIBLE from outside':>21s}")
print("  "+"-"*80)
rows=[]
for (nx,ny),(ax,ay) in [((4,4),(1,1)),((5,5),(1,1)),((5,5),(2,2)),((6,6),(2,2)),
                        ((6,6),(3,3)),((7,7),(3,3)),((7,7),(4,4)),((8,8),(4,4))]:
    L,PL,MX,MY=lattice(nx,ny)
    if ax>MX or ay>MY: continue
    A=[j*MX+i for j in range(ay) for i in range(ax)]          # bottom-left block of plaquettes
    Ain=set(A)
    linksA=set()
    for p in A: linksA|=set(PL[p])
    outside=[p for p in range(len(PL)) if p not in Ain]
    # STORAGE: independent records inside A that can be simultaneously protected (W-42 law)
    storage=len(A)-1
    # perimeter of A = links on exactly one plaquette of A, plus lattice-edge links of A
    cnt={}
    for p in A:
        for lk in PL[p]: cnt[lk]=cnt.get(lk,0)+1
    perim=sum(1 for lk,v in cnt.items() if v==1)
    # LEGIBILITY: an outside observer's gauge-invariant operators are Wilson loops built from
    # plaquettes OUTSIDE A. Two internal records differ for that observer iff some outside loop
    # distinguishes them. Over GF(2) an outside loop W_T (T a subset of outside plaquettes) and an
    # internal record W_S commute unless their link supports overlap oddly -- and internal records
    # live on links of A. So the separating power is the rank of the pairing matrix.
    # The pairing <T,S> = |bd(T) INTERSECT bd(S)| mod 2 is BILINEAR over GF(2), so the rank over
    # ALL outside loops equals the rank over a BASIS of them -- the single outside plaquettes.
    # (Enumerating all subsets is both unnecessary and, at these sizes, fatal: the first version of
    #  this lane was killed by the OS at exit 137 building that list.)
    basisA=[(p,) for p in A]
    bdA=[bd(S,PL) for S in basisA]
    vecs=[]
    for q in outside:
        bt=bd((q,),PL)
        v=0
        for i,b in enumerate(bdA):
            if len(bt & b)%2: v|=(1<<i)
        vecs.append(v)
    legible=rank_gf2(vecs)
    rows.append((nx,ny,ax,ay,len(A),perim,storage,legible))
    print(f"  {f'{nx}x{ny}':>9s} {f'{ax}x{ay}':>10s} {len(A):9d} {perim:9d} {storage:12d} {legible:21d}")

print()
print("  READING THE TWO SCALINGS")
print(f"  {'|A|':>5s} {'perim':>6s} {'storage':>8s} {'legible':>8s} {'storage/|A|':>12s} {'legible/perim':>14s}")
for nx,ny,ax,ay,na,pe,st,le in rows:
    print(f"  {na:5d} {pe:6d} {st:8d} {le:8d} {st/max(na,1):12.3f} {le/max(pe,1):14.3f}")
print()
print("  If STORAGE tracks |A| (bulk) and LEGIBLE tracks perimeter, then holography is a statement")
print("  about what escapes the region, not about what fits in it -- and BOTH the principal's point 4")
print("  and W-42's volume law are correct about different quantities.")
print("  If LEGIBLE is a constant independent of size, neither scaling holds and the region is")
print("  externally almost mute, which is a different and stronger claim.")
