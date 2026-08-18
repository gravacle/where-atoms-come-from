"""W-37e.  (1) RESOLVE THE Z_3 NULL.  (2) MEASURE REDUNDANCY DIRECTLY INSTEAD OF ARGUING IT.

(1) W-37d's Z_3 rings returned EXACTLY 0.000000 at every size, including even ones where Z_2 gives
0.999792. Diagnosis to test: the "probe reduced density matrix" was formed by matching gauge configs
s == t ACROSS probe positions -- but physical states at different positions satisfy DIFFERENT Gauss
laws, so no such pair exists and the matrix was always diagonal. Worse, there is no canonical
identification of gauge states across positions at all: the gauge field is NOT a tensor factor.
That is this program's own non-factorisation problem showing up again.
FIX: use the CLASSICAL mutual information of the probe's POSITION, which needs no partial trace and
is gauge-invariant.   I(R;P) = sum_r sum_p p(r) P(p|r) log2[ P(p|r) / P(p) ].

(2) W-37 argued redundancy from ||[R,H_hop]|| = 0 rather than measuring it. Put TWO probes on the
ring at once (at Z_2 two charges need no anti-charge -- total charge is even) and ask what EACH
probe knows on its own. Redundancy means each single probe already holds the bit.
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

def sectors_of(R,D):
    ev=np.linalg.eigvals(R); vals=np.unique(np.round(ev,6))
    Ps=[]
    for lam in vals:
        P=np.eye(D,dtype=complex)
        for mu in vals:
            if abs(mu-lam)>1e-9: P=P@((R-mu*np.eye(D))/(lam-mu))
        Ps.append(P)
    return vals,Ps

def cmi(dists,weights):
    """classical mutual information between the sector label and the position, in bits"""
    dists=np.array(dists); w=np.array(weights); w=w/w.sum()
    marg=(w[:,None]*dists).sum(0)
    I=0.0
    for r in range(len(w)):
        for p in range(dists.shape[1]):
            a=dists[r,p]
            if a>1e-15 and marg[p]>1e-15: I+=w[r]*a*np.log2(a/marg[p])
    return I

# ---------- (1) one probe, Z_N ring, classical MI ----------
def ring1(n,N):
    Ew=[(k,(k+1)%n) for k in range(n)]
    def dv(s,v): return (sum(s[k] for k,(a,b) in enumerate(Ew) if a==v)
                        -sum(s[k] for k,(a,b) in enumerate(Ew) if b==v))%N
    ST=[(p,s) for p in range(n) for s in itertools.product(range(N),repeat=n)
        if all(dv(s,v)==((1 if v==p else 0)-(1 if v==0 else 0))%N for v in range(n))]
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,k,d):
        t=list(s); t[k]=(t[k]+d)%N; return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST):
        t=s
        for k in range(n): t=sh(t,k,1)
        R[IDX[(p,t)],j]=1.0
    H=np.zeros((D,D),complex)
    for k in range(n):
        for j,(p,s) in enumerate(ST):
            if p!=k: continue
            i=IDX.get(((k+1)%n,sh(s,k,1)))
            if i is not None: H[i,j]-=1.0
    return ST,D,R,H+H.conj().T

print("W-37e (1)  Z_3 NULL RESOLVED?  classical MI of the probe POSITION (no partial trace).")
print(f"  {'n':>3s} {'N':>3s} {'dim':>5s} {'sectors':>8s} {'max I(R;position) bits':>23s} {'ceiling log2(N)':>16s}")
for n in [4,5,6,7]:
    for N in [2,3]:
        ST,D,R,H=ring1(n,N); vals,Ps=sectors_of(R,D)
        g=np.random.default_rng(5)
        w=(g.normal(size=D)+1j*g.normal(size=D))*np.array([1.0 if p==1 else 0.0 for p,_ in ST])
        parts=[]
        ok=True
        for P in Ps:
            v=P@w
            if np.linalg.norm(v)<1e-10: ok=False; break
            parts.append(v/np.linalg.norm(v))
        if not ok: print(f"  {n:3d} {N:3d} {D:5d} {len(vals):8d} {'(no unbiased start)':>23s}"); continue
        psi=sum(parts); psi/=np.linalg.norm(psi)
        U=expm(-1j*H*0.25); best=0.0; t=0.0
        while t<30.0:
            psi=U@psi; t+=0.25
            ds=[];ws=[]
            for P in Ps:
                v=P@psi; pr=float(np.vdot(v,v).real)
                if pr<1e-14: continue
                v=v/np.sqrt(pr)
                d=np.zeros(n)
                for i,(p,s) in enumerate(ST): d[p]+=abs(v[i])**2
                ds.append(d); ws.append(pr)
            best=max(best,cmi(ds,ws))
        print(f"  {n:3d} {N:3d} {D:5d} {len(vals):8d} {best:23.6f} {np.log2(N):16.4f}")

# ---------- (2) two probes at once, Z_2 ----------
def ring2(n):
    Ew=[(k,(k+1)%n) for k in range(n)]
    def dv(s,v): return (sum(s[k] for k,(a,b) in enumerate(Ew) if a==v)
                        -sum(s[k] for k,(a,b) in enumerate(Ew) if b==v))%2
    ST=[(p,q,s) for p in range(n) for q in range(n) if p!=q
        for s in itertools.product(range(2),repeat=n)
        if all(dv(s,v)==((1 if v==p else 0)^(1 if v==q else 0)) for v in range(n))]
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,k):
        t=list(s); t[k]^=1; return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(p,q,s) in enumerate(ST):
        t=s
        for k in range(n): t=sh(t,k)
        R[IDX[(p,q,t)],j]=1.0
    H=np.zeros((D,D),complex)
    for j,(p,q,s) in enumerate(ST):
        for (mover,other,isfirst) in ((p,q,True),(q,p,False)):
            np_=(mover+1)%n
            key=(np_,q,sh(s,mover)) if isfirst else (p,np_,sh(s,mover))
            i=IDX.get(key)
            if i is not None: H[i,j]-=1.0
    return ST,D,R,H+H.conj().T

print("\nW-37e (2)  TWO PROBES AT ONCE (Z_2, no anti-charge needed). What does EACH one know alone?")
print(f"  {'n':>3s} {'dim':>5s} {'I(R;probe1)':>13s} {'I(R;probe2)':>13s} {'I(R;both)':>11s}   reading")
for n in [4,6]:
    ST,D,R,H=ring2(n); vals,Ps=sectors_of(R,D)
    g=np.random.default_rng(9)
    w=(g.normal(size=D)+1j*g.normal(size=D))*np.array([1.0 if (p==1 and q==2) else 0.0 for p,q,_ in ST])
    parts=[]
    for P in Ps:
        v=P@w; parts.append(v/np.linalg.norm(v))
    psi=sum(parts); psi/=np.linalg.norm(psi)
    U=expm(-1j*H*0.25); b1=b2=bb=0.0; t=0.0
    while t<30.0:
        psi=U@psi; t+=0.25
        d1=[];d2=[];db=[];ws=[]
        for P in Ps:
            v=P@psi; pr=float(np.vdot(v,v).real)
            if pr<1e-14: continue
            v=v/np.sqrt(pr)
            a=np.zeros(n); b=np.zeros(n); c=np.zeros(n*n)
            for i,(p,q,s) in enumerate(ST):
                w2=abs(v[i])**2; a[p]+=w2; b[q]+=w2; c[p*n+q]+=w2
            d1.append(a); d2.append(b); db.append(c); ws.append(pr)
        b1=max(b1,cmi(d1,ws)); b2=max(b2,cmi(d2,ws)); bb=max(bb,cmi(db,ws))
    rd=("REDUNDANT: each probe alone already holds the bit" if min(b1,b2)>0.5
        else "not redundant: only the pair knows")
    print(f"  {n:3d} {D:5d} {b1:13.6f} {b2:13.6f} {bb:11.6f}   {rd}")
