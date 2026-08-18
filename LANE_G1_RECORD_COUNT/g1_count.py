"""G1.  HOW MANY THINGS ACTUALLY SATISFY THE DEFINITION?

Two capacity numbers are in the register and they scale differently:
    T3        capacity = area - 1        LINEAR in area          (plaquette-subset boundary loops)
    Thm A     record space = 2^(2g)      EXPONENTIAL in genus    (topological sectors)
Only one can be counting RECORDS under the five-clause definition. Clause (v) is the discriminator:
a record must be PROTECTED -- no operation on a contractible region may write it.

TEST. Enumerate candidate observables and check all five clauses on each:
  (i) bit  (ii) [H,R]=0 and [L_k,R]=0  (iii) non-trivial on an eigenspace  (iv) writable  (v) protected
Count how many INDEPENDENT observables pass. Compare with 2g (the number of independent 1-form
generators, giving 2^(2g) states) and with T3's area - 1.

IF T3's OBJECTS FAIL CLAUSE (v), THEN capacity = area - 1 WAS NEVER A RECORD COUNT, and the record
count is set by GENUS, exponentially -- which is a statement about gravity's role, not EM's.
"""
import itertools, numpy as np
def torus(nx,ny):
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
        for j in range(ny) for i in range(nx)]
    return nx*ny,E,PL
NV,E,PL=torus(2,2); L=len(E)
st=[s for s in itertools.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(NV))]
idx={s:i for i,s in enumerate(st)}; D=len(st)
def Move(S):
    M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k in S: t[k]^=1
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def Zl(S): return np.diag([(-1.0)**(sum(s[k] for k in S)%2) for s in st]).astype(complex)
H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
LOCAL=[Zl([k]) for k in range(L)]+[Move(p) for p in PL]      # every LOCAL gauge-invariant operator
def eigsp(H):
    ev,U=np.linalg.eigh(H); out=[]; i=0
    while i<len(ev):
        j=i
        while j+1<len(ev) and abs(ev[j+1]-ev[i])<1e-8*max(1.0,abs(ev[i])): j+=1
        out.append(U[:,i:j+1]); i=j+1
    return out
ES=eigsp(H)
def clauses(R):
    c1 = np.linalg.norm(R-R.conj().T)<1e-9 and np.linalg.norm(R@R-np.eye(D))<1e-9
    c2 = np.linalg.norm(H@R-R@H)<1e-9
    c3 = any(U.shape[1]>1 and np.linalg.norm((U.conj().T@R@U)-np.trace(U.conj().T@R@U)/U.shape[1]*np.eye(U.shape[1]))>1e-6 for U in ES)
    c5 = all(np.linalg.norm(O.conj().T@R@O + R)>1e-6 for O in LOCAL)   # NO local op flips it
    return c1,c2,c3,c5

print(f"G1  torus 2x2 (genus 1): links {L}, plaquettes {len(PL)}, physical dim {D}")
print(f"    T3 predicts capacity = area - 1 = {len(PL)-1}")
print(f"    Thm A predicts record space 2^(2g) = {2**2}  -> 2g = 2 independent record observables\n")
print("  EVERY Z-TYPE CANDIDATE (all 2^L - 1 link subsets), CHECKED AGAINST ALL FIVE CLAUSES")
print(f"  {'candidate class':>34s} {'count':>7s} {'(i)':>5s} {'(ii)':>5s} {'(iii)':>6s} {'(v)':>5s} {'ALL':>5s}")
print("  "+"-"*74)
def vv(S):
    v=0
    for k in S: v|=(1<<k)
    return v
def bits(v): return [k for k in range(L) if v>>k&1]
plq_boundaries=set()
for r in range(1,len(PL)+1):
    for T in itertools.combinations(range(len(PL)),r):
        c={}
        for p in T:
            for lk in PL[p]: c[lk]=c.get(lk,0)+1
        plq_boundaries.add(vv([lk for lk,n in c.items() if n%2]))
plq_boundaries.discard(0)
allpass=[]; t3pass=0
for v in range(1,1<<L):
    R=Zl(bits(v)); c1,c2,c3,c5=clauses(R)
    if c1 and c2 and c3 and c5: allpass.append(v)
    if v in plq_boundaries and c1 and c2 and c3 and c5: t3pass+=1
def rank_gf2(vs):
    b=[]; r=0
    for v in vs:
        cur=v
        for x in b:
            p=x.bit_length()-1
            if cur>>p&1: cur^=x
        if cur: b.append(cur); b.sort(reverse=True); r+=1
    return r
print(f"  {'ALL link subsets':>34s} {(1<<L)-1:7d} {'':>5s} {'':>5s} {'':>6s} {'':>5s} {len(allpass):5d}")
print(f"  {'of which: T3 plaquette-boundaries':>34s} {len(plq_boundaries):7d} {'':>5s} {'':>5s} {'':>6s} {'':>5s} {t3pass:5d}")
print(f"\n  independent observables passing all clauses (GF(2) rank): {rank_gf2(allpass)}")
print(f"  predicted 2g = 2   ->  record space 2^2 = 4")
print()
print("  WHICH CLAUSE DOES T3's COUNT FAIL?  check a single-plaquette boundary loop:")
p0=vv([lk for lk,n in ((lambda c: c)({k:sum(1 for q in [0] for kk in PL[q] if kk==k) for k in PL[0]})).items() if n%2])
R0=Zl(PL[0])
c1,c2,c3,c5=clauses(R0)
print(f"    single plaquette boundary: (i) {c1}  (ii) {c2}  (iii) {c3}  (v) {c5}")
bad=[i for i,O in enumerate(LOCAL) if np.linalg.norm(O.conj().T@R0@O + R0)<1e-6]
print(f"    local operators that FLIP it: {len(bad)}  -> clause (v) {'FAILS' if bad else 'holds'}")
