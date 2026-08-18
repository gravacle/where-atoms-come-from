"""W-37d.  IS THE ODD-LOOP NULL DEEP, OR IS IT A Z_2 COINCIDENCE?

W-37c found odd rings carry EXACTLY zero information and attributed it to the antiperiodic momenta
being the periodic ones shifted by pi. PROPOSED SHARPER MECHANISM, tested here:
  a Z_2 flux of pi and a SIGN FLIP OF EVERY HOPPING are the same operation when n is odd, because
  negating all n amplitudes shifts the total flux by n*pi = pi (odd) or 0 (even).
  With real hopping the position distribution is time-symmetric, so the two sectors become
  literally indistinguishable BY POSITION.
If that is the mechanism then the null is a property of Z_2 (where pi is the only nontrivial flux,
and it is exactly the sign-flip value) and of what the probe reports -- NOT of records.

PREDICTIONS, stated before running:
  P1  odd n: spectrum(-H(flux 0)) == spectrum(H(flux pi)) exactly. Even n: not.
  P2  at Z_3 an odd ring SHOULD be readable, because 2pi/3 is not a sign flip.
  P3  at Z_2 an odd ring SHOULD become readable to a probe with COMPLEX hopping (which breaks the
      time-symmetry of the position distribution).
Any prediction failing means the mechanism is wrong and something deeper is present.
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
def vn(r):
    ev=np.linalg.eigvalsh((r+r.conj().T)/2); ev=ev[ev>1e-12]
    return float(-(ev*np.log2(ev)).sum())

print("P1 -- is a pi flux the same as flipping every hopping sign?  (single-particle ring)")
print(f"   {'n':>4s} {'spec(-H(0)) == spec(H(pi))?':>30s}")
for n in [4,5,6,7,8,9,10,11]:
    H0=np.zeros((n,n),complex); Hp=np.zeros((n,n),complex)
    for j in range(n):
        H0[(j+1)%n,j]=-1.0
        Hp[(j+1)%n,j]=-1.0*(-1.0 if j==0 else 1.0)   # put the pi flux on one link
    H0=H0+H0.conj().T; Hp=Hp+Hp.conj().T
    a=np.sort(np.round(np.linalg.eigvalsh(-H0),9)); b=np.sort(np.round(np.linalg.eigvalsh(Hp),9))
    print(f"   {n:4d} {('YES  -> sectors indistinguishable by position' if np.allclose(a,b) else 'no'):>30s}")

def ring_ZN(n,N,tau=1.0,phase=0.0):
    """probe hopping a ring of n sites, Z_N gauge field, anti-charge pinned at site 0."""
    Ew=[(k,(k+1)%n) for k in range(n)]
    def dv(s,v): return (sum(s[k] for k,(a,b) in enumerate(Ew) if a==v)
                        -sum(s[k] for k,(a,b) in enumerate(Ew) if b==v))%N
    ST=[]
    for p in range(n):
        for s in itertools.product(range(N),repeat=n):
            ok=all(dv(s,v)==((1 if v==p else 0)-(1 if v==0 else 0))%N for v in range(n))
            if ok: ST.append((p,s))
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,k,d):
        t=list(s); t[k]=(t[k]+d)%N; return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST):
        t=s
        for k in range(n): t=sh(t,k,1)
        R[IDX[(p,t)],j]=1.0
    H=np.zeros((D,D),complex)
    amp=-tau*np.exp(1j*phase)
    for k in range(n):
        for j,(p,s) in enumerate(ST):
            if p!=k: continue
            i=IDX.get(((k+1)%n, sh(s,k,1)))
            if i is not None: H[i,j]+=amp
    return ST,IDX,D,R,H+H.conj().T

def readout(n,N,phase=0.0,tmax=30.0,dt=0.25,seed=5):
    ST,IDX,D,R,H=ring_ZN(n,N,phase=phase)
    ev=np.linalg.eigvals(R); vals=np.unique(np.round(ev,6))
    Ps=[]
    for lam in vals:
        Pk=np.eye(D,dtype=complex)
        for mu in vals:
            if abs(mu-lam)>1e-9: Pk=Pk@((R-mu*np.eye(D))/(lam-mu))
        Ps.append(Pk)
    g=np.random.default_rng(seed)
    w=(g.normal(size=D)+1j*g.normal(size=D))*np.array([1.0 if p==1 else 0.0 for p,_ in ST])
    parts=[]
    for Pk in Ps:
        v=Pk@w
        if np.linalg.norm(v)<1e-10: return None,len(vals),D
        parts.append(v/np.linalg.norm(v))
    psi=sum(parts); psi/=np.linalg.norm(psi)
    U=expm(-1j*H*dt); best=0.0; t=0.0
    while t<tmax:
        psi=U@psi; t+=dt
        br=[]
        for Pk in Ps:
            v=Pk@psi; pr=float(np.vdot(v,v).real)
            if pr<1e-14: br.append((0.0,None)); continue
            v/=np.sqrt(pr)
            M=np.zeros((n,n),complex)
            for i,(p,s) in enumerate(ST):
                for jj,(q,t2) in enumerate(ST):
                    if s==t2: M[p,q]+=v[i]*np.conj(v[jj])
            br.append((pr,M))
        avg=sum(pp*m for pp,m in br if m is not None)
        best=max(best, vn(avg)-sum(pp*vn(m) for pp,m in br if m is not None))
    return best,len(vals),D

print("\nP2 -- Z_3 rings. If the null is a Z_2 coincidence, ODD rings must read fine at Z_3.")
print(f"   {'n':>4s} {'N':>3s} {'dim':>5s} {'sectors':>8s} {'max I(R:probe) bits':>20s}")
for n in [4,5,6,7]:
    for N in [2,3]:
        r,ns,D=readout(n,N)
        print(f"   {n:4d} {N:3d} {D:5d} {ns:8d} {(f'{r:.6f}' if r is not None else 'no unbiased start'):>20s}")

print("\nP3 -- Z_2 ODD ring, probe given a COMPLEX hopping phase (breaks position time-symmetry).")
print(f"   {'n':>4s} {'phase':>8s} {'max I(R:probe) bits':>20s}")
for n in [5,7]:
    for ph in [0.0, 0.3, 0.7854, 1.5708]:
        r,ns,D=readout(n,2,phase=ph)
        print(f"   {n:4d} {ph:8.4f} {(f'{r:.6f}' if r is not None else 'no unbiased start'):>20s}")
