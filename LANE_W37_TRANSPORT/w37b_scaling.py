"""W-37b.  DOES THE READOUT TIME SCALE WITH THE PERIMETER?

W-37 measured one ring and asserted T_read ~ perimeter/tau as an expectation. Measure it.

CARRIER: a bare cycle graph of n vertices and n links -- the minimal object with a closed path and
nothing else. Probe charge hops the ring; static anti-charge pinned at vertex 0 (without it the
total Z_2 charge is odd and the sector is empty -- rechecked at every n).
Gauss: div(s)_v = [v == probe site] + [v == 0]  (mod 2).  Physical dim = 2n.
R = the Wilson loop around the whole ring (flip every link). Verified [R,H_hop] = 0 at every n, so
the transport never disturbs what it reads.

T_read := the FIRST time I(R:probe) reaches 1/2 bit. Reported against n. A straight line through
the origin in (n, T_read) is the perimeter law; the log-log slope is printed rather than assumed.
"""
import itertools, numpy as np

def expm(A):
    nr=np.linalg.norm(A,np.inf)
    n=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**n); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for k in range(1,40):
        T=T@B/k; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(n): X=X@X
    return X
def vn(r):
    ev=np.linalg.eigvalsh((r+r.conj().T)/2); ev=ev[ev>1e-12]
    return float(-(ev*np.log2(ev)).sum())

def ring(n, tau=1.0, cut=False):
    Ework=[(k,(k+1)%n) for k in range(n)]
    def div(s,v):
        return (sum(s[k] for k,(a,b) in enumerate(Ework) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(Ework) if b==v))%2
    ST=[]
    for p in range(n):
        for s in itertools.product(range(2),repeat=n):
            q=lambda v:(1 if v==p else 0)^(1 if v==0 else 0)
            if all(div(s,v)==q(v) for v in range(n)): ST.append((p,s))
    IDX={x:i for i,x in enumerate(ST)}; D=len(ST)
    def sh(s,ls):
        t=list(s)
        for k in ls: t[k]^=1
        return tuple(t)
    R=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST): R[IDX[(p,sh(s,range(n)))],j]=1.0
    H=np.zeros((D,D),complex)
    for k in range(n):
        if cut and k==0: continue
        for j,(p,s) in enumerate(ST):
            if p!=k: continue
            i=IDX.get(((k+1)%n, sh(s,[k])))
            if i is not None: H[i,j]-=tau
    H=H+H.conj().T
    return ST,IDX,D,R,H,n

def measure(n, tau=1.0, cut=False, seed=5, tmax=60.0, dt=0.25):
    ST,IDX,D,R,H,nn=ring(n,tau,cut)
    assert D==2*n, (D,n)
    com=np.linalg.norm(R@H-H@R)
    Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2
    g=np.random.default_rng(seed)
    mask=np.array([1.0 if p==1 else 0.0 for p,_ in ST])
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    a=Pp@w; b=Pm@w
    if np.linalg.norm(a)<1e-12 or np.linalg.norm(b)<1e-12: return None,com,D
    a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    psi0=(a+b); psi0/=np.linalg.norm(psi0)
    def info(psi):
        br=[]
        for Pr in (Pp,Pm):
            v=Pr@psi; pr=float(np.vdot(v,v).real)
            if pr<1e-14: br.append((0.0,None)); continue
            v=v/np.sqrt(pr)
            M=np.zeros((n,n),complex)
            for i,(p,s) in enumerate(ST):
                for jj,(q,t) in enumerate(ST):
                    if s==t: M[p,q]+=v[i]*np.conj(v[jj])
            br.append((pr,M))
        avg=sum(p*m for p,m in br if m is not None)
        return vn(avg)-sum(p*vn(m) for p,m in br if m is not None)
    U=expm(-1j*H*dt); psi=psi0.astype(complex); t=0.0; best=0.0; tread=None
    while t<tmax:
        psi=U@psi; t+=dt
        I=info(psi); best=max(best,I)
        if tread is None and I>=0.5: tread=t
    return (tread,best),com,D

print("W-37b  DOES READOUT TIME SCALE WITH PERIMETER?  bare cycle graphs, tau=1.0")
print(f"  {'n':>4s} {'dim':>5s} {'||[R,H]||':>10s} {'T_read (I=0.5 bit)':>19s} {'max I over run':>15s}")
print("  "+"-"*60)
rows=[]
for n in [4,5,6,7,8,10,12]:
    out,com,D=measure(n)
    if out is None: print(f"  {n:4d} {D:5d} {com:10.1e}   (no unbiased start)"); continue
    tr,best=out
    rows.append((n,tr))
    print(f"  {n:4d} {D:5d} {com:10.1e} {(f'{tr:.2f}' if tr else 'never'):>19s} {best:15.4f}")

good=[(n,t) for n,t in rows if t]
if len(good)>=3:
    x=np.log([n for n,_ in good]); y=np.log([t for _,t in good])
    sl,ic=np.polyfit(x,y,1)
    print(f"\n  d ln(T_read) / d ln(n) = {sl:+.3f}   "
          f"({'PERIMETER LAW: readout time is linear in the loop length' if abs(sl-1)<0.25 else 'NOT linear -- the perimeter expectation is wrong'})")
    print(f"  ratios T_read/n : " + "  ".join(f"n={n}:{t/n:.3f}" for n,t in good))

print("\n  CONTROL -- cut ring at every n. No closed path, so I must never reach 0.5 bit.")
for n in [4,6,8,10]:
    out,com,D=measure(n,cut=True)
    tr,best=out if out else (None,0.0)
    print(f"    n={n:3d}  T_read = {tr if tr else 'never':>7}   max I over the whole run = {best:.2e}")
