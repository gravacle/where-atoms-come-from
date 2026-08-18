"""W-48.  DOES OUR COMPOSED LAW SAY ANYTHING THE STANDARD DECOHERENCE ACCOUNT DOES NOT?

Asked BEFORE touching any experimental data, deliberately. If our law is just exp(-Gamma T) wearing
lattice clothes, then fitting it to C70 or Panda data would be curve-matching a textbook result and
calling it confirmation. The principal's own guard.

STANDARD ACCOUNT: coherence decays exponentially, V = exp(-Gamma T). Readability at the readout time
should then be EXPONENTIAL in X = Gamma * T_read.
OUR MEASUREMENT (W-47): I/I0 vs X looked shallower than exponential at large X.

SO FIT BOTH, over the widest range the carrier allows:
    exponential   I/I0 = exp(-a X)
    power law     I/I0 = c X^(-b)
and report which the data prefers, with residuals. This is decidable.

AND THE OBVIOUS ARTIFACT MUST BE EXCLUDED FIRST. W-47 reported the MAX over time of I(t). A max over
a decaying-but-still-rising signal is an optimisation, and optimising an exponential against a
rising readout can manufacture a power law. So report BOTH:
    I_max   = max over t          (what W-47 used)
    I_fixed = I at the FIXED readout time T_read = 0.30 P/tau, no optimisation
If I_fixed is exponential and only I_max is a power law, the power law is ours, not nature's, and
the honest verdict is that we add nothing to the standard account.
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
    H=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(ST):
        i=IDX.get(((p+1)%n, sh(s,p)))
        if i is not None: H[i,j]+= -1.0*np.exp(1j*0.6)
    return ST,IDX,D,R,H+H.conj().T

def curves(n,tau,Gamma):
    """returns (I_max over t, I at the fixed readout time T_read = 0.30 n/tau)"""
    ST,IDX,D,R,H0=build(n); H=H0*tau
    Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2
    Zl=np.diag([(-1.0)**s[0] for p_,s in ST]).astype(complex)
    assert np.linalg.norm(Zl@R+R@Zl)<1e-9
    Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))+ (Gamma/2.0)*(np.kron(Zl,Zl.conj())-np.kron(Id,Id))
    g=np.random.default_rng(5)
    mask=np.array([1.0 if p==1 else 0.0 for p,_ in ST])
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    a=Pp@w; b=Pm@w
    if min(np.linalg.norm(a),np.linalg.norm(b))<1e-12: return None
    a/=np.linalg.norm(a); b/=np.linalg.norm(b); psi=(a+b); psi/=np.linalg.norm(psi)
    Tread=0.30*n/tau
    Tmax=max(4.0*Tread, 6.0); steps=120; dt=Tmax/steps
    U=expm(M*dt); v=np.outer(psi,psi.conj()).reshape(-1)
    best=0.0; at_fixed=None
    for k in range(1,steps+1):
        v=U@v; r=v.reshape(D,D); r=r/np.trace(r).real
        ds=[];ws=[]
        for Pr in (Pp,Pm):
            rr=Pr@r@Pr; pr=np.trace(rr).real
            if pr<1e-12: continue
            rr=rr/pr; d=np.zeros(n)
            for i,(p,s) in enumerate(ST): d[p]+=rr[i,i].real
            ds.append(d); ws.append(pr)
        if len(ds)==2:
            I=cmi(ds,ws); best=max(best,I)
            if at_fixed is None and k*dt>=Tread: at_fixed=I
    return best,(at_fixed if at_fixed is not None else 0.0)

print("W-48  IS OUR LAW JUST exp(-Gamma T)?")
print(f"  {'n':>3s} {'tau':>5s} {'g^2':>6s} {'X':>8s} {'Imax/I0':>9s} {'Ifix/I0':>9s}")
print("  "+"-"*46)
rows=[]
for n in (6,10,12):
    for tau in (1.0,2.0):
        b=curves(n,tau,0.0)
        if b is None: continue
        I0m,I0f=b
        if I0m<1e-6 or I0f<1e-6: continue
        for g2 in (0.04,0.06,0.09,0.13,0.18,0.25,0.34,0.45):
            X=9.6*(g2**2)*n/tau
            r=curves(n,tau,32.0*g2**2)
            if r is None: continue
            rows.append((n,tau,X,r[0]/I0m,r[1]/I0f))
            print(f"  {n:3d} {tau:5.1f} {g2:6.3f} {X:8.3f} {r[0]/I0m:9.4f} {r[1]/I0f:9.4f}")

X=np.array([r[2] for r in rows]); Fm=np.array([r[3] for r in rows]); Ff=np.array([r[4] for r in rows])
def fits(F,label):
    m=(F>1e-4)&(X>0.02)
    x,f=X[m],F[m]
    ae=np.polyfit(x,np.log(f),1); rese=np.log(f)-np.polyval(ae,x)
    ap=np.polyfit(np.log(x),np.log(f),1); resp=np.log(f)-np.polyval(ap,np.log(x))
    print(f"\n  {label}   n={m.sum()}")
    print(f"    EXPONENTIAL  log(F) = {ae[0]:.3f} X + {ae[1]:.3f}      rms resid {np.sqrt((rese**2).mean()):.4f}")
    print(f"    POWER LAW    log(F) = {ap[0]:.3f} logX + {ap[1]:.3f}   rms resid {np.sqrt((resp**2).mean()):.4f}")
    win = "POWER LAW" if np.sqrt((resp**2).mean()) < np.sqrt((rese**2).mean()) else "EXPONENTIAL"
    print(f"    -> {win} fits better (by {abs(np.sqrt((rese**2).mean())-np.sqrt((resp**2).mean())):.4f} in rms log-residual)")
    return win
w_max = fits(Fm,"I_max  (max over t -- what W-47 reported)")
w_fix = fits(Ff,"I_fixed (at T_read, NO optimisation)")
print()
print("  VERDICT")
if w_fix=="EXPONENTIAL":
    print("    The un-optimised signal IS exponential in X. Standard decoherence accounts for it.")
    print("    Any power law in I_max is OURS -- produced by maximising over time -- not nature's.")
    print("    => our composed law adds NOTHING to the standard account on this axis, and fitting it")
    print("       to C70/Panda data would be curve-matching a textbook result.")
else:
    print("    The un-optimised signal is NOT exponential. That is a genuine departure from the")
    print("    standard account and is the thing worth taking to data.")
