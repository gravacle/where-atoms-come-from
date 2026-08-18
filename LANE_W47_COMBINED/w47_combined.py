"""W-47.  THE COMBINED LAW, AS A FALSIFIABLE PREDICTION.

Every lane so far measured ONE relation. This composes three measured relations into a single
prediction and then tries to break it. That is the difference between a collection of results and
a law.

THE THREE INPUTS, each measured elsewhere in this program, none assumed here:
  EM        the record is a holonomy on a closed boundary of perimeter P             (W-27, W-34)
  ALPHA     it decays at   Gamma = 32 g^4                                            (W-32)
            reading it by transport takes  T_read = 0.30 * P / tau                   (W-37b)
  CAPACITY  only  P/2 - 1  of a region's records are externally distinguishable       (W-46)

COMPOSED PREDICTION. A record is exportable only if it can be read before it decays:
        Gamma * T_read < 1     =>     9.6 * g^4 * P / tau  <  1
Define the composed control parameter
        X = 9.6 * g^4 * P / tau
Then the prediction is a SHARP FUNCTION OF X ALONE -- not of P, g and tau separately. Three
parameters must collapse onto one curve. THAT is the falsifiable content: if the transition happens
at different X for different (P, g, tau), the composition is wrong.

HOW IT COULD FAIL, stated before running:
  * the curves do NOT collapse -> the three relations do not compose; the law is wrong.
  * the transition sits at X far from 1 -> the composition is right in form but the constants are
    wrong, which is a weaker failure and is reported as such.
  * information stays high at large X -> decay does not actually gate the reading.

CONSTRUCTION. A ring of n sites (perimeter P = n) carrying a Z_2 gauge field, a probe charge hopping
at rate tau with a static anti-charge, and a dephasing channel that damages the record at rate
Gamma. The probe's position information about the loop is measured under that dissipation. Gamma is
imposed directly rather than re-derived, since W-32 already measured Gamma = 32 g^4; this lane tests
the COMPOSITION, not the ingredients.
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

def build(n):
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
    return ST,IDX,D,R

def hop(ST,IDX,D,n,tau,phase=0.6):
    H=np.zeros((D,D),complex)
    def sh(s,k):
        t=list(s); t[k]^=1; return tuple(t)
    for j,(p,s) in enumerate(ST):
        i=IDX.get(((p+1)%n, sh(s,p)))
        if i is not None: H[i,j]+= -tau*np.exp(1j*phase)
    return H+H.conj().T

def run(n, tau, Gamma, T=None):
    """max classical mutual information the probe's POSITION holds about the loop,
       while the loop is being dephased at rate Gamma."""
    ST,IDX,D,R=build(n)
    H=hop(ST,IDX,D,n,tau)
    Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2
    # DECAY CHANNEL. The first version used L = P+ - P-, which IS R. W-30a proves a unitary jump
    # operator commuting with the record cannot change <R> at all, so Gamma did nothing and the
    # information came out identical to six decimals across every g^2. Seventh instance of this
    # program's recurring defect. The channel must ANTICOMMUTE with the loop: at Z_2, R is the
    # product of all link shifts and Z_k anticommutes with the shift on link k, so L = Z_k damages
    # the record. W-41's counting formula then gives rate = 2*gamma per bath link on the boundary.
    Zl=np.diag([(-1.0)**s[0] for p_,s in ST]).astype(complex)     # Z on link 0: anticommutes with R
    assert np.linalg.norm(Zl@R + R@Zl) < 1e-9, "channel must anticommute with the record"
    Ld=[Zl]
    Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    gam=Gamma/2.0     # rate = 2*gamma for one bath link on the boundary (W-41)
    for X in Ld: M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    g=np.random.default_rng(5)
    mask=np.array([1.0 if p==1 else 0.0 for p,_ in ST])
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    a=Pp@w; b=Pm@w
    if min(np.linalg.norm(a),np.linalg.norm(b))<1e-12: return None
    a/=np.linalg.norm(a); b/=np.linalg.norm(b); psi=(a+b); psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    T = T if T is not None else max(4.0, 3.0*0.30*n/tau)
    dt=T/60.0; U=expm(M*dt); v=rho.reshape(-1); best=0.0
    for _ in range(60):
        v=U@v; r=v.reshape(D,D); r=r/np.trace(r).real
        ds=[];ws=[]
        for Pr in (Pp,Pm):
            rr=Pr@r@Pr; pr=np.trace(rr).real
            if pr<1e-12: continue
            rr=rr/pr; d=np.zeros(n)
            for i,(p,s) in enumerate(ST): d[p]+=rr[i,i].real
            ds.append(d); ws.append(pr)
        if len(ds)==2: best=max(best,cmi(ds,ws))
    return best

print("W-47 v3  THE COMPOSED LAW.   X = 9.6 * g^4 * P / tau")
print("   inputs, each measured elsewhere: Gamma = 32 g^4 (W-32); T_read = 0.30 P/tau (W-37b);")
print("   externally distinguishable records = P/2 - 1 (W-46). Composition: a record is exportable")
print("   only if it is read before it decays, Gamma*T_read < 1, i.e. X < 1.")
print("")
print("   v1 used a decay channel that COMMUTED with the record, which W-30a proves can do nothing;")
print("   Gamma had no effect and I was identical to six decimals across every g^2. Fixed: the")
print("   channel now ANTICOMMUTES with the loop (asserted in code).")
print("   v2 compared raw I across n, but BARE readability differs by n (the probe reads n=8 poorly")
print("   at this phase whatever the decay). The composed law predicts the SURVIVING FRACTION, so")
print("   the quantity is  I(Gamma) / I(Gamma=0), normalised per (n,tau).")
print("")
print(f"  {'n=P':>4s} {'tau':>6s} {'g^2':>7s} {'X':>9s} {'I':>10s} {'I0 (no decay)':>14s} {'I/I0':>8s}")
print("  " + "-"*66)
rows=[]
for n in (6,8,10,12):
    for tau in (1.0,2.0):
        I0=run(n,tau,0.0)
        if I0 is None or I0<1e-6: continue
        for g2 in (0.05,0.08,0.11,0.14,0.18,0.24,0.32,0.42):
            X=9.6*(g2**2)*n/tau
            I=run(n,tau,32.0*g2**2)
            if I is None: continue
            rows.append((n,tau,g2,X,I,I0,I/I0))
            print(f"  {n:4d} {tau:6.2f} {g2:7.3f} {X:9.3f} {I:10.5f} {I0:14.5f} {I/I0:8.4f}")
print("")
print("  COLLAPSE TEST -- the falsifiable content. Points with similar X but different (P,tau) must")
print("  agree on I/I0. If they do not, the three relations do not compose.")
import collections
by=collections.defaultdict(list)
for n,tau,g2,X,I,I0,f in rows: by[round(np.log(X)*2)/2].append(f)
print(f"  {'X':>9s} {'#pts':>5s} {'I/I0 values':>40s} {'spread':>8s}")
for k in sorted(by):
    v=by[k]
    if len(v)<2: continue
    print(f"  {np.exp(k):9.3f} {len(v):5d} {str([f'{x:.3f}' for x in v]):>40s} {max(v)-min(v):8.3f}")
xs=np.array([r[3] for r in rows]); fs=np.array([r[6] for r in rows])
print("")
print("  IS THE SURVIVING FRACTION A FUNCTION OF X ALONE?")
for lo,hi,lab in [(0,0.5,'X < 0.5   predicted LEGIBLE'),(0.5,1.5,'0.5-1.5   transition'),
                  (1.5,4.0,'1.5-4.0'),(4.0,1e9,'X > 4     predicted ILLEGIBLE')]:
    m=(xs>=lo)&(xs<hi)
    if m.sum(): print(f"    {lab:34s} n={int(m.sum()):3d}  mean I/I0 = {fs[m].mean():.4f}  sd = {fs[m].std():.4f}")
half=xs[np.argsort(np.abs(fs-0.5))[:4]]
print("")
print(f"    X where half the readability is lost (4 nearest points): {sorted(round(float(x),2) for x in half)}")
print(f"    predicted crossing: X = 1")
