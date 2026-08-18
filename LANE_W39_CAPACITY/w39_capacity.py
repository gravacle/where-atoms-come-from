"""W-39.  HOW MANY RECORDS CAN ONE CARRIER HOLD AT ONCE?

The principal: dynamical geometry is unlikely to appear at the level of a SINGLE record; it would
emerge where MANY records must share a carrier. A single record is a test particle, and test
particles never source geometry.

W-38: a record survives when its boundary is DISJOINT from where the environment couples. So with
two records, each one's environment must also avoid the OTHER's boundary. That is a packing
problem, and it is the first place in this program where records can constrain EACH OTHER.

FORCED-OR-NOT, FIRST. Two separate bounds, and they must not be conflated:
  (a) ALGEBRAIC, and FORCED: the 15 region-boundaries live in the Z_2 vector space of plaquette
      subsets, dim = cycle rank = 4. So at most 4 INDEPENDENT records exist, and that is just the
      dimension of the physical space -- a counting fact, not a discovery.
  (b) GEOMETRIC, and NOT forced: whether a set of records can be SIMULTANEOUSLY PROTECTED, i.e.
      whether a nonempty bath exists that is disjoint from ALL their boundaries at once.
Only (b) is a result. Both are computed and reported separately.

Then the combinatorics is CHECKED DYNAMICALLY: for a surviving set, run the sieve with that bath and
confirm every member of the set is actually slow.
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

REG={}
for r in range(1,5):
    for S in itertools.combinations(range(4),r):
        mv=compose([P[i] for i in S])
        REG[S]=(frozenset(k for k,_ in mv), Move(st,idx,mv,N))
names=list(REG)

print("W-39  HOW MANY RECORDS CAN THIS CARRIER HOLD AT ONCE?")
print(f"      links {L}, plaquettes 4, cycle rank 4, physical dim {D}")
print()
print("  (a) ALGEBRAIC BOUND -- FORCED, reported so it cannot be mistaken for a result.")
print("      region boundaries span the Z_2 space of plaquette subsets, dim 4.")
print(f"      => at most 4 INDEPENDENT records. This is dim of the physical space, nothing more.")
print()
print("  (b) GEOMETRIC BOUND -- NOT forced. Can a set be SIMULTANEOUSLY PROTECTED?")
print("      A set is simultaneously protectable iff some NONEMPTY bath avoids every boundary in it.")
print()
best={}
for k in range(1,5):
    ok=[]
    for combo in itertools.combinations(names,k):
        # independence over GF(2) on plaquette subsets
        vecs=[sum(1<<i for i in S) for S in combo]
        rank=0; basis=[]
        for v in vecs:
            for b in basis: v=min(v,v^b)
            if v: basis.append(v); basis.sort(reverse=True); rank+=1
        if rank<k: continue                      # not independent records
        union=set().union(*[REG[S][0] for S in combo])
        free=[x for x in range(L) if x not in union]
        if free: ok.append((combo,free))
    best[k]=ok
    if ok:
        combo,free=ok[0]
        print(f"      k={k}: {len(ok):4d} independent sets are simultaneously protectable."
              f"  example {[''.join(map(str,S)) for S in combo]} with bath {free}")
    else:
        print(f"      k={k}: {len(ok):4d} -- NO independent set of {k} records can be protected at once.")

kmax=max(k for k in best if best[k])
print(f"\n  => GEOMETRIC CAPACITY OF THIS CARRIER = {kmax}   (algebraic bound was 4)")
print(f"     The two DISAGREE, so the limit is packing, not dimension.\n"
      if kmax<4 else "     They agree here.\n")

def rates(links,g2=0.01,gam=0.5):
    H=-MAG-g2*ELEC
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for X in [Zop(st,[k],N) for k in links]: M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    w,U=np.linalg.eig(M.conj().T)
    rate=-np.conj(w).real; U=U/np.linalg.norm(U,axis=0)
    out={}
    for S,(sup,O) in REG.items():
        v=(O/np.linalg.norm(O)).reshape(-1)
        ov=np.abs(U.conj().T@v); ov=ov/max(ov.sum(),1e-30)
        out[S]=float((ov*rate).sum())
    return out

print("  DYNAMICAL CHECK -- the combinatorics above is a claim about the dynamics. Test it.")
for k in sorted(best):
    if not best[k]: continue
    combo,free=best[k][0]
    r=rates(free)
    inside=[r[S] for S in combo]
    outside=[v for S,v in r.items() if S not in combo]
    print(f"    k={k}  bath {free}")
    for S in combo:
        print(f"         member {str(S):12s} rate {r[S]:.4e}   boundary {sorted(REG[S][0])}")
    print(f"         slowest NON-member rate {min(outside):.4e}   "
          f"=> members {'ALL slower' if max(inside)<min(outside) else 'NOT all slower'}"
          f" (margin {min(outside)/max(max(inside),1e-30):.1f}x)")

print("\n  CONTROL -- a set that the combinatorics says CANNOT be protected. Its members must not")
print("  all be slow under any bath disjoint from... (there is none, so use the largest partial bath).")
bad=None
for combo in itertools.combinations(names,2):
    union=set().union(*[REG[S][0] for S in combo])
    if not [x for x in range(L) if x not in union]: bad=combo; break
if bad:
    union=set().union(*[REG[S][0] for S in bad])
    print(f"    pair {bad} covers all {len(union)} links -- no bath avoids both.")
    r=rates([0,1])
    print(f"    with an arbitrary bath [0,1]: rates {[f'{r[S]:.3e}' for S in bad]}")
else:
    print("    every pair leaves some link free on this carrier.")
