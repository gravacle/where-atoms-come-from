"""W-32b.  IS THE WINDOW REAL?  Three checks the first pass did not make.

 (1) IDENTIFY THE MODE. w32 took the slowest decaying eigenvalue of the WHOLE Lindbladian. That is
     only the record's lifetime if the corresponding eigenmode actually carries <R>. Check the
     overlap |Tr(R^dag * mode)| for the slowest mode, and report the slowest mode that DOES carry R.
 (2) SLOPE IN THE RIGHT REGIME. w32's summary slope used the scan endpoints and straddled the
     non-monotonic region. Fit only the weak-coupling regime, and push to smaller g^2.
 (3) IS IT JUST PERTURBATION THEORY? omega ~ g^2 (first order) and Gamma ~ g^4 (second order) give
     RATIO ~ 1/g^2 by golden-rule counting alone. Test whether the LOCATION of the bath changes the
     exponent -- if it does, the exponent is physics; if every bath gives the same law, it is
     generic weak-coupling counting and must be reported as such.
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
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
CENTER=vid[(1,1)]; CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
PERIM=[k for k in range(len(E)) if k not in CUT]
N=2; st,idx=build(V2,E,N); D=len(st)
R=Move(st,idx,compose(P),N)
MAG=sum((lambda L:L+L.conj().T)(Move(st,idx,p,N)) for p in P)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
Id=np.eye(D,dtype=complex)
Rv=R.reshape(-1); Rv=Rv/np.linalg.norm(Rv)

def spectrum(g2, links, gam=0.5):
    H=-MAG-g2*ELEC
    Ls=[Zop(st,[k],N) for k in links]
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    for L in Ls: M+=gam*(np.kron(L.conj(),L)-np.kron(Id,Id))
    w,V=np.linalg.eig(M)
    # LEFT eigenvectors carry observables: (M^dag) u = conj(lam) u  <=>  d<u|rho>/dt = lam <u|rho>
    wl,U=np.linalg.eig(M.conj().T)
    ov=np.abs(U.conj().T@Rv)          # overlap of each left-eigenmode with the observable R
    ov=ov/np.linalg.norm(U,axis=0)
    return w, np.conj(wl), ov, np.linalg.norm(H@R-R@H)/D

print("W-32b  CHECK 1 -- does the slowest mode actually carry the record?")
print(f"  {'g^2':>8s} {'Gam(slowest ANY)':>17s} {'ov with R':>10s} {'Gam(slowest CARRYING R)':>24s} {'ov':>8s}")
print("  "+"-"*74)
for g2 in [0.0,0.002,0.005,0.01,0.02,0.05,0.1]:
    w,wl,ov,om=spectrum(g2,CUT)
    rate=-wl.real
    m=rate>1e-10
    if not m.any(): print(f"  {g2:8.4f}   no decaying modes"); continue
    i_any=np.where(m)[0][np.argmin(rate[m])]
    carry=m & (ov>1e-6)
    if carry.any():
        i_R=np.where(carry)[0][np.argmin(rate[carry])]
        print(f"  {g2:8.4f} {rate[i_any]:17.6e} {ov[i_any]:10.2e} {rate[i_R]:24.6e} {ov[i_R]:8.3f}")
    else:
        print(f"  {g2:8.4f} {rate[i_any]:17.6e} {ov[i_any]:10.2e} {'NO MODE CARRIES R':>24s} {'--':>8s}")

def ratio_scan(links, gs, label):
    print(f"\n  {label}")
    print(f"  {'g^2':>9s} {'omega':>12s} {'Gamma_R':>14s} {'RATIO':>12s} {'local slope d lnRATIO/d ln g2':>30s}")
    print("  "+"-"*84)
    prev=None; out=[]
    for g2 in gs:
        w,wl,ov,om=spectrum(g2,links)
        rate=-wl.real; m=(rate>1e-10)&(ov>1e-6)
        if not m.any(): print(f"  {g2:9.5f}  no R-carrying decaying mode"); continue
        G=rate[m].min(); r=om/G; out.append((g2,om,G,r))
        sl=""
        if prev: sl=f"{(np.log(r)-np.log(prev[1]))/(np.log(g2)-np.log(prev[0])):30.3f}"
        print(f"  {g2:9.5f} {om:12.6f} {G:14.6e} {r:12.4f} {sl:>30s}")
        prev=(g2,r)
    return out

print()
print("CHECK 2 -- the weak-coupling regime, pushed down two more decades.")
gs=[0.0005,0.001,0.002,0.005,0.01,0.02,0.05]
a=ratio_scan(CUT,gs,"bath on the CUT (disjoint from the rim; [L,R]=0)")
b=ratio_scan(PERIM,gs,"bath on the RIM (the environment sees the record directly; [L,R]!=0)")

print()
print("CHECK 3 -- does the bath LOCATION change the exponent, or is it generic golden-rule counting?")
for nm,o in [("CUT",a),("RIM",b)]:
    if len(o)>=3:
        x=np.log([r[0] for r in o]); y=np.log([r[3] for r in o])
        sl=np.polyfit(x,y,1)[0]
        yg=np.log([r[2] for r in o]); slg=np.polyfit(x,yg,1)[0]
        print(f"  bath on {nm:4s}:  d ln(RATIO)/d ln(g2) = {sl:+.3f}    d ln(Gamma)/d ln(g2) = {slg:+.3f}")
print()
print("  RATIO ~ 1/g^2 (slope -1) with Gamma ~ g^4 (slope +2) is second-order golden-rule protection:")
print("  the bath can only reach the record THROUGH the H-mixing, so decay is one order higher than drive.")
print("  If the two baths give DIFFERENT exponents, the exponent is a fact about where the environment")
print("  couples. If they give the same, it is generic weak-coupling counting and must be said so.")
