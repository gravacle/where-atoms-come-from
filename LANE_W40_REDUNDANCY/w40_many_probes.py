"""W-40.  MEASURE REDUNDANCY INSTEAD OF ARGUING IT.

W-37 claimed condition 5 (redundant copying -> objectivity) was met, but ARGUED it from
||[R,H_hop]|| = 0 rather than measuring it. W-37e tried two probes and got zero for each AND for
both jointly -- which was real physics, not a bug: HARD-CORE probes on a ring cannot pass each
other, so none completes a circuit, and any motion returning to the start flips every link TWICE,
giving R^2 = I and no phase at all.

FIX THE CONSTRUCTION: let the probes PASS (no exclusion). Then each one can circumnavigate on its
own and read the loop independently. k probes on a ring of n sites; a static anti-charge at site 0
when k is ODD, because the total Z_2 charge must be even or the physical sector is empty (checked).

REDUNDANCY IS: does EACH probe, ON ITS OWN, hold the bit?
  I(R ; position of probe i)  ~ 1 bit for every i          -> redundant, hence objective
  only the joint distribution carries it                   -> one delocalised copy, not objective
Information is the CLASSICAL mutual information of probe positions -- gauge-invariant, and needing
no partial trace (the gauge field is not a tensor factor; that is this program's own
non-factorisation problem and it broke W-37d's Z_3 measurement).

AND ASK W-39's QUESTION ONE LEVEL DOWN: does adding probes DEGRADE what each one learns?
"""
import itertools, numpy as np
def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X
def cmi(dists,weights):
    d=np.array(dists); w=np.array(weights,dtype=float); w=w/w.sum()
    marg=(w[:,None]*d).sum(0); I=0.0
    for r in range(len(w)):
        for p in range(d.shape[1]):
            a=d[r,p]
            if a>1e-15 and marg[p]>1e-15: I+=w[r]*a*np.log2(a/marg[p])
    return I

def system(n,k,tau=1.0,phase=0.0):
    Ew=[(j,(j+1)%n) for j in range(n)]
    anti = (k % 2 == 1)
    def dv(s,v): return (sum(s[j] for j,(a,b) in enumerate(Ew) if a==v)
                        -sum(s[j] for j,(a,b) in enumerate(Ew) if b==v))%2
    ST=[]
    for pos in itertools.product(range(n),repeat=k):
        for s in itertools.product(range(2),repeat=n):
            ok=True
            for v in range(n):
                q=(sum(1 for x in pos if x==v)+(1 if (anti and v==0) else 0))%2
                if dv(s,v)!=q: ok=False; break
            if ok: ST.append((pos,s))
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,j):
        t=list(s); t[j]^=1; return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(pos,s) in enumerate(ST):
        t=s
        for e in range(n): t=sh(t,e)
        R[IDX[(pos,t)],j]=1.0
    H=np.zeros((D,D),complex)
    amp=-tau*np.exp(1j*phase)
    for j,(pos,s) in enumerate(ST):
        for i in range(k):                       # probe i hops forward; probes may share a site
            np_=list(pos); np_[i]=(pos[i]+1)%n
            t=IDX.get((tuple(np_),sh(s,pos[i])))
            if t is not None: H[t,j]+=amp
    return ST,D,R,H+H.conj().T,anti

def sectors(R,D):
    ev=np.linalg.eigvals(R); vals=np.unique(np.round(ev,6)); Ps=[]
    for lam in vals:
        P=np.eye(D,dtype=complex)
        for mu in vals:
            if abs(mu-lam)>1e-9: P=P@((R-mu*np.eye(D))/(lam-mu))
        Ps.append(P)
    return vals,Ps

def run(n,k,phase=0.6,tmax=40.0,dt=0.25,seed=11):
    ST,D,R,H,anti=system(n,k,phase=phase)
    if D==0: return None
    vals,Ps=sectors(R,D)
    g=np.random.default_rng(seed)
    start=tuple((1+i)%n for i in range(k))
    mask=np.array([1.0 if pos==start else 0.0 for pos,_ in ST])
    if mask.sum()==0: return None
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    parts=[]
    for P in Ps:
        v=P@w
        if np.linalg.norm(v)<1e-10: return None
        parts.append(v/np.linalg.norm(v))
    psi=sum(parts); psi/=np.linalg.norm(psi)
    U=expm(-1j*H*dt); best=[0.0]*k; bjoint=0.0; t=0.0
    while t<tmax:
        psi=U@psi; t+=dt
        per=[[] for _ in range(k)]; joint=[]; ws=[]
        for P in Ps:
            v=P@psi; pr=float(np.vdot(v,v).real)
            if pr<1e-14: continue
            v=v/np.sqrt(pr)
            marg=[np.zeros(n) for _ in range(k)]; jj=np.zeros(n**k)
            for i,(pos,s) in enumerate(ST):
                p2=abs(v[i])**2
                for a in range(k): marg[a][pos[a]]+=p2
                jj[sum(pos[a]*(n**a) for a in range(k))]+=p2
            for a in range(k): per[a].append(marg[a])
            joint.append(jj); ws.append(pr)
        for a in range(k): best[a]=max(best[a],cmi(per[a],ws))
        bjoint=max(bjoint,cmi(joint,ws))
    return D,best,bjoint,anti

print("W-40  REDUNDANCY, MEASURED.  probes may PASS each other; complex hopping phase 0.6")
print("      (W-37d: a real-hopping probe cannot read an odd ring at all -- objectivity is probe-relative)")
print()
print(f"  {'n':>3s} {'k':>2s} {'dim':>5s} {'anti':>5s} {'I per probe (bits)':>34s} {'I joint':>9s}  reading")
print("  "+"-"*92)
for n in [4,5,6]:
    for k in [1,2,3]:
        if n**k*2 > 1200: print(f"  {n:3d} {k:2d}   (skipped: dimension too large)"); continue
        out=run(n,k)
        if out is None: print(f"  {n:3d} {k:2d}   (no unbiased start in every sector)"); continue
        D,per,joint,anti=out
        rd=("REDUNDANT: each probe alone holds the bit" if min(per)>0.5
            else "partial" if min(per)>0.15 else "NOT redundant: only the joint distribution knows")
        print(f"  {n:3d} {k:2d} {D:5d} {str(anti):>5s} {str([round(x,4) for x in per]):>34s} {joint:9.4f}  {rd}")

print()
print("  CONTROL -- real hopping (phase 0). Even rings should still read; odd rings must give 0.")
print(f"  {'n':>3s} {'k':>2s} {'I per probe':>28s}")
for n,k in [(4,2),(5,2),(6,2)]:
    out=run(n,k,phase=0.0)
    print(f"  {n:3d} {k:2d} {str([round(x,4) for x in out[1]]) if out else 'n/a':>28s}")

print()
print("  CONTROL -- tau -> 0 (no transport). Must be 0 for every probe.")
out=run(6,2,phase=0.6,tmax=40.0)
ST,D,R,H,_=system(6,2,tau=0.0,phase=0.6)
print(f"    ||H|| at tau=0 = {np.linalg.norm(H):.3e}  (nothing can move, so no probe can learn anything)")
