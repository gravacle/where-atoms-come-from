"""W-43.  INSTALL THE METRIC TERM AND ASK WHETHER IT DOES ANYTHING THE GAUGE FIELD DOES NOT.

CONTAMINATION NOTICE: this lane was designed AFTER reading Gravacle v337. Everything here is
downstream of its Field-Registration Principle (EM enters records through gauge/action phase;
gravity through metric/proper-time action). W-27..W-42 were independent; this is not.

THE CARRIER HAS NEVER HAD A METRIC. Its links have no length, no proper time, no clock. So install
one and find out whether it is DISTINGUISHABLE from the gauge field at the record level.
  If NOT distinguishable -> "registration-level unification" is a degeneracy, and a metric adds
                            nothing to records that a gauge flux does not already add.
  If distinguishable     -> gravity is an independent record-bearing structure here.

THE DISCRIMINATOR, and it is not a matter of taste. In the continuum:
  EM phase     = line integral of a 1-form  -> ODD under path reversal (reverse the path, flip sign)
  proper time  = arc length / elapsed time  -> EVEN under path reversal (reverse it, same value)
So a Hermitian hop carries EM as a LINK PHASE (forward e^{iA}, backward e^{-iA}: odd) and gravity as
a SITE POTENTIAL (a diagonal term; the phase it accrues depends on TIME SPENT, not direction: even).
That is why gravitational redshift is not an Aharonov-Bohm effect.

MEASURED CONSEQUENCE: an EM flux breaks time-reversal and makes the probe CHIRAL -- it drifts one
way round the ring. A potential cannot. So CHIRALITY separates them, and it is a number.
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
def cmi(d,w):
    d=np.array(d); w=np.array(w,float); w=w/w.sum(); marg=(w[:,None]*d).sum(0); I=0.0
    for r in range(len(w)):
        for p in range(d.shape[1]):
            a=d[r,p]
            if a>1e-15 and marg[p]>1e-15: I+=w[r]*a*np.log2(a/marg[p])
    return I

def ring(n, A=None, V=None, tau=1.0):
    """Z_2 gauge ring + probe + static anti-charge at site 0.
       A[k] = EM vector potential on link k  (enters the HOP phase: ODD under reversal)
       V[j] = metric/proper-time potential at site j (DIAGONAL: EVEN under reversal)"""
    A=np.zeros(n) if A is None else np.asarray(A,float)
    V=np.zeros(n) if V is None else np.asarray(V,float)
    Ew=[(j,(j+1)%n) for j in range(n)]
    def dv(s,v): return (sum(s[k] for k,(a,b) in enumerate(Ew) if a==v)
                        -sum(s[k] for k,(a,b) in enumerate(Ew) if b==v))%2
    ST=[(p,s) for p in range(n) for s in itertools.product(range(2),repeat=n)
        if all(dv(s,v)==((1 if v==p else 0)^(1 if v==0 else 0)) for v in range(n))]
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,k):
        t=list(s); t[k]^=1; return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST):
        t=s
        for k in range(n): t=sh(t,k)
        R[IDX[(p,t)],j]=1.0
    H=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST):
        i=IDX.get(((p+1)%n, sh(s,p)))
        if i is not None: H[i,j]+= -tau*np.exp(1j*A[p])
    H=H+H.conj().T
    for j,(p,s) in enumerate(ST): H[j,j]+=V[p]
    return ST,D,R,H

def sectors(R,D):
    ev=np.linalg.eigvals(R); vals=np.unique(np.round(ev,6)); Ps=[]
    for lam in vals:
        P=np.eye(D,dtype=complex)
        for mu in vals:
            if abs(mu-lam)>1e-9: P=P@((R-mu*np.eye(D))/(lam-mu))
        Ps.append(P)
    return vals,Ps

def evolve_stats(n,A=None,V=None,T=25.0,dt=0.25,seed=5):
    """returns (max I(R;probe) over time, chirality = net signed drift round the ring)"""
    ST,D,R,H=ring(n,A,V); vals,Ps=sectors(R,D)
    g=np.random.default_rng(seed)
    mask=np.array([1.0 if p==0 else 0.0 for p,_ in ST])
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    parts=[]
    for P in Ps:
        v=P@w
        if np.linalg.norm(v)<1e-10: return None,None
        parts.append(v/np.linalg.norm(v))
    psi=sum(parts); psi/=np.linalg.norm(psi)
    U=expm(-1j*H*dt); best=0.0; chir=[]; t=0.0
    disp=np.array([((p+n//2)%n)-n//2 for p in range(n)],float)   # signed displacement from site 0
    while t<T:
        psi=U@psi; t+=dt
        ds=[];ws=[]; tot=np.zeros(n)
        for P in Ps:
            v=P@psi; pr=float(np.vdot(v,v).real)
            if pr<1e-14: continue
            v=v/np.sqrt(pr); d=np.zeros(n)
            for i,(p,s) in enumerate(ST): d[p]+=abs(v[i])**2
            ds.append(d); ws.append(pr); tot+=pr*d
        best=max(best,cmi(ds,ws)); chir.append(float(tot@disp))
    return best, float(np.mean(chir))

n=6
print("W-43  METRIC vs GAUGE ON THE SAME CARRIER.  ring n=6, Z_2 gauge + probe + anti-charge")
print("      EM  = link phase A_k  (hop phase; ODD under path reversal)")
print("      GRAV= site potential V_j (diagonal; phase accrues with TIME SPENT; EVEN under reversal)")
print()
print(f"  {'setting':>34s} {'I(R;probe) bits':>16s} {'chirality (signed drift)':>25s}")
print("  "+"-"*80)
cases=[
 ("nothing added",                 None,                              None),
 ("EM uniform A=0.30",             [0.30]*n,                          None),
 ("EM uniform A=0.60",             [0.60]*n,                          None),
 ("EM single link A_0=1.80",       [1.80]+[0.0]*(n-1),                None),
 ("METRIC uniform V=0.30",         None,                              [0.30]*n),
 ("METRIC uniform V=0.60",         None,                              [0.60]*n),
 ("METRIC gradient V=0..1.0",      None,                              list(np.linspace(0,1.0,n))),
 ("METRIC well V at one site",     None,                              [1.5]+[0.0]*(n-1)),
 ("BOTH  A=0.30 and V grad",       [0.30]*n,                          list(np.linspace(0,1.0,n))),
]
res={}
for nm,A,V in cases:
    I,c=evolve_stats(n,A,V)
    res[nm]=(I,c)
    print(f"  {nm:>34s} {I:16.6f} {c:25.6f}")

print()
print("  THE DISCRIMINATOR -- chirality. An EM flux breaks time reversal; a potential cannot.")
em=[abs(res[k][1]) for k in res if k.startswith("EM")]
gr=[abs(res[k][1]) for k in res if k.startswith("METRIC")]
print(f"    |chirality| with EM terms      : {[f'{x:.4f}' for x in em]}")
print(f"    |chirality| with METRIC terms  : {[f'{x:.4f}' for x in gr]}")
print(f"    baseline (nothing added)       : {abs(res['nothing added'][1]):.6f}")
print(f"    -> {'SEPARABLE: EM is chiral, the metric is not' if min(em)>10*max(gr+[1e-12]) else 'NOT separated by this measure'}")

print()
print("  DEGENERACY CHECK -- can a METRIC setting reproduce an EM setting's full signature?")
print(f"    {'pair':>46s} {'dI':>10s} {'d chirality':>13s}")
for a in [k for k in res if k.startswith("EM")]:
    for b in [k for k in res if k.startswith("METRIC")]:
        dI=abs(res[a][0]-res[b][0]); dc=abs(res[a][1]-res[b][1])
        if dI<0.05:
            print(f"    {a+'  vs  '+b:>46s} {dI:10.4f} {dc:13.4f}"
                  f"   {'DEGENERATE' if dc<0.01 else 'separated by chirality alone'}")
