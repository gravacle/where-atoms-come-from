"""W-41.  BUILD k=4 AND SEE WHAT BREAKS.

W-39: ALGEBRAIC capacity 4 (cycle rank), GEOMETRIC capacity 3 -- no bath avoids four independent
boundaries at once. So k=4 is a packing constraint WITH NO SOLUTION. Demand it and watch the
failure mode, because the failure mode is the physics.

THE 4-SET: the four single plaquettes (a basis of the region space). Their boundaries union to ALL
TWELVE LINKS, so EVERY nonempty bath touches at least one of them.

A DERIVATION FIRST, WHICH ALSO MAKES THE SCAN CHEAP.  For the adjoint Lindbladian
  d/dt <O> :  L*(O) = i[H,O] + gamma sum_k (Z_k O Z_k - O),
the coherent term contributes NOTHING to the first-order decay rate: Tr(O* i[H,O]) is purely
imaginary for Hermitian H, so its real part vanishes. At Z_2 a Wilson loop O and a link operator
Z_k either COMMUTE (k not on the boundary) or ANTICOMMUTE (k on it), so Z_k O Z_k = +-O and

        rate(O)  =  2 * gamma * |bath INTERSECT boundary(O)|   + O(g^4).

The decay rate is literally COUNTING how many bath links sit on the boundary. That is the whole
W-32/W-38 picture in one line, and it makes the packing question a covering problem.
CHECKED, not assumed: the formula is compared against the full 256x256 Lindbladian spectrum below.
"""
import itertools, numpy as np
def build(V,E,N):
    st=[s for s in itertools.product(range(N),repeat=len(E))
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(len(V)))]
    return st,{s:i for i,s in enumerate(st)}
def Zop(st,links,N):
    w=np.exp(2j*np.pi/N)
    return np.diag([w**(sum(s[k] for k in links)%N) for s in st]).astype(complex)
def Move(st,idx,mv,N):
    D=len(st); M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%N
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def compose(ps):
    acc={}
    for p in ps:
        for k,sg in p: acc[k]=acc.get(k,0)+sg
    return [(k,s) for k,s in acc.items() if s!=0]
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
N=2; st,idx=build(V2,E,N); D=len(st); Id=np.eye(D,dtype=complex)
MAG=sum((lambda X:X+X.conj().T)(Move(st,idx,p,N)) for p in P)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(L))
FOUR=[(i,) for i in range(4)]
OPS={S:Move(st,idx,compose([P[i] for i in S]),N) for S in FOUR}
SUP={S:set(k for k,_ in compose([P[i] for i in S])) for S in FOUR}

print("W-41  k=4: the four single plaquettes")
for S in FOUR: print(f"      region {S}  boundary {sorted(SUP[S])}")
mult={k:sum(1 for S in FOUR if k in SUP[S]) for k in range(L)}
print(f"      union = {sorted(set().union(*SUP.values()))}  -> every nonempty bath hits a boundary")
print(f"      link multiplicity (how many of the four boundaries each link lies on):")
print(f"        {dict(sorted(mult.items()))}")

def fast(bath,S): return 2.0*0.5*len(SUP[S]&set(bath))
def full(bath,g2=0.01,gam=0.5):
    H=-MAG-g2*ELEC
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for X in [Zop(st,[k],N) for k in bath]: M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    w,U=np.linalg.eig(M.conj().T); rate=-np.conj(w).real; U=U/np.linalg.norm(U,axis=0)
    o={}
    for S in FOUR:
        v=(OPS[S]/np.linalg.norm(OPS[S])).reshape(-1)
        ov=np.abs(U.conj().T@v); ov=ov/max(ov.sum(),1e-30)
        o[S]=float((ov*rate).sum())
    return o

print("\n  VALIDATE the counting formula against the full 256x256 Lindbladian spectrum.")
print(f"    {'bath':>14s} {'region':>8s} {'|bath∩bdy|':>11s} {'2*gamma*n':>10s} {'full spectrum':>14s}")
for bath in ([0],[0,1],[2],[2,3],[0,4,8]):
    fu=full(bath)
    for S in FOUR:
        print(f"    {str(bath):>14s} {str(S):>8s} {len(SUP[S]&set(bath)):11d} {fast(bath,S):10.3f} {fu[S]:14.4e}")
    print()

baths=[list(c) for sz in (1,2,3,4) for c in itertools.combinations(range(L),sz)]
print(f"  SCAN {len(baths)} baths (sizes 1-4) with the counting formula.")
rows=[]
for b in baths:
    v=sorted(fast(b,S) for S in FOUR)
    rows.append((b,v))
best=min(rows,key=lambda z:z[1][-1])
print(f"\n  BEST 'PROTECT ALL FOUR': worst rate {best[1][-1]:.3f} at bath {best[0]}, profile {best[1]}")
print(f"    number of baths achieving worst-rate 0 (all four protected): "
      f"{sum(1 for b,v in rows if v[-1]==0)}   <- zero means the constraint really has no solution")

print(f"\n  FAILURE MODE.  For every minimal (single-link) bath, how many records survive?")
print(f"    {'bath':>6s} {'multiplicity':>13s} {'rates (sorted)':>26s} {'protected':>10s} {'evicted':>8s}")
for k in range(L):
    v=sorted(fast([k],S) for S in FOUR)
    print(f"    {str([k]):>6s} {mult[k]:13d} {str([f'{x:.1f}' for x in v]):>26s} "
          f"{sum(1 for x in v if x==0):10d} {sum(1 for x in v if x>0):8d}")

ev=[(b,v) for b,v in rows if sum(1 for x in v if x==0)==3]
print(f"\n  => baths that protect EXACTLY THREE and evict ONE: {len(ev)}")
print(f"     baths that degrade all four together (none at rate 0): "
      f"{sum(1 for b,v in rows if all(x>0 for x in v))}")
gap=[v[-1]/max([x for x in v if x>0][0],1e-30) for b,v in rows if sum(1 for x in v if x==0)==3]
print(f"\n  VERDICT read off the profiles, not asserted:")
print(f"     with a single-link bath on a multiplicity-1 link, THREE records sit at rate EXACTLY 0")
print(f"     and ONE sits at rate {2*0.5*1:.1f}. The carrier does not share the damage.")
print(f"\n  CONFIRM ON THE FULL SPECTRUM -- bath [0], which should evict region (0,) only:")
fu=full([0])
for S in FOUR:
    print(f"     region {S}  full-spectrum rate {fu[S]:.4e}   {'EVICTED' if fu[S]>0.1 else 'protected (g^4 residue)'}")
