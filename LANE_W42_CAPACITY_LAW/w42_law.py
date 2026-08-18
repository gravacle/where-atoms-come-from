"""W-42.  DOES CAPACITY GROW WITH CARRIER SIZE?  (prerequisite for "content forces expansion")

W-41: a 4-plaquette patch holds 3 records and evicts the 4th. For "demand exceeds capacity -> the
carrier grows" to be a real dynamics, capacity must actually depend on size. Measure the law.

W-41's counting formula makes this pure linear algebra over GF(2), with no Hilbert space at all:
  * the boundary map  S |-> bd(S)  (links lying in an ODD number of the plaquettes of S) is LINEAR;
  * a set of records is simultaneously protectable iff some link L is on NONE of their boundaries;
  * "bd(S) misses link L" is ONE linear condition  f_L(S) = 0;
  * so the largest protectable independent set via L is exactly dim ker f_L.
Hence CAPACITY = max over links L of dim ker(f_L) = m - min_L rank(f_L), m = number of plaquettes.

THIS MATTERS FOR HOW W-39 WAS REPORTED. W-39 called capacity 3 a GEOMETRIC bound as opposed to the
FORCED algebraic bound 4. If the identity above holds, the geometric bound is ALSO forced -- it is
the dimension of the kernel of a single linear functional -- and W-39's framing needs correcting.
Computed here for several patches, and cross-checked against W-39's brute force at m=4.
"""
import itertools, numpy as np

def patch(nx,ny):
    """nx x ny VERTICES -> (nx-1)*(ny-1) plaquettes"""
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
    return len(V),len(E),PL

def rank_gf2(rows,ncols):
    rows=[r for r in rows]; piv=[]; r=0
    for c in range(ncols):
        p=None
        for i in range(r,len(rows)):
            if (rows[i]>>c)&1: p=i; break
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and ((rows[i]>>c)&1): rows[i]^=rows[r]
        r+=1
        if r==len(rows): break
    return r

print("W-42  CAPACITY vs CARRIER SIZE")
print(f"  {'patch':>9s} {'verts':>6s} {'links':>6s} {'plaq m':>7s} {'cycle rank':>11s} {'capacity':>9s} {'m-1?':>6s}")
print("  "+"-"*62)
res=[]
for nx,ny in [(2,2),(3,2),(3,3),(4,3),(4,4),(5,4),(5,5)]:
    NV,L,PL=patch(nx,ny); m=len(PL)
    cyc=L-NV+1
    # f_L as a bitmask over plaquettes: bit p set iff link L lies on plaquette p
    best=0
    for lk in range(L):
        f=sum(1<<p for p,pl in enumerate(PL) if lk in pl)
        # dim ker of the single functional f, acting on the m-dim plaquette space
        kerdim = m - (1 if f!=0 else 0)
        best=max(best,kerdim)
    res.append((nx,ny,m,cyc,best))
    print(f"  {f'{nx}x{ny}':>9s} {NV:6d} {L:6d} {m:7d} {cyc:11d} {best:9d} {str(best==m-1):>6s}")

print("\n  CROSS-CHECK against W-39's brute force at the 3x3 patch (m=4). Brute force said 3.")
NV,L,PL=patch(3,3); m=len(PL)
def bd(S):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return frozenset(lk for lk,v in c.items() if v%2)
regions=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
found=0
for k in range(1,m+1):
    ok=False
    for combo in itertools.combinations(regions,k):
        vecs=[sum(1<<i for i in S) for S in combo]
        if rank_gf2(vecs,m)<k: continue
        if set(range(L))-set().union(*[bd(S) for S in combo]): ok=True; break
    if ok: found=k
print(f"    brute force over all independent sets: capacity = {found}   "
      f"(linear-algebra answer {res[2][4]})  -> {'AGREE' if found==res[2][4] else 'DISAGREE'}")

print("\n  THE LAW")
print(f"    capacity = m - 1 in every case above, where m = number of plaquettes.")
print(f"    For a planar patch the number of plaquettes IS the AREA in lattice cells,")
print(f"    and the cycle rank equals it. So:  CAPACITY = AREA - 1.")
print("\n  AND THE HONEST READING OF THAT")
for nx,ny,m,cyc,cap in res:
    per = 2*((nx-1)+(ny-1))
    print(f"    {nx}x{ny}: area {m:3d}  perimeter {per:3d}  capacity {cap:3d}   "
          f"cap/area {cap/max(m,1):.3f}   cap/perimeter {cap/max(per,1):.3f}")
print("\n    capacity tracks AREA, not perimeter. It is NOT a holographic/boundary-law count,")
print("    and must not be reported as one.")
