"""W-30b.  WHERE DOES g^2 SIT RELATIVE TO THE OBSTRUCTION?  (v2)

W-30a: the no-go has two horns.  g2=0 -> [H,R]=0, d<R>/dt=0 EXACTLY (durable, unwritable).
g2>0 -> R moves (writable, and dissipation gets a handle). g2 is the dial BETWEEN them.
That was not installed. It is where the coupling landed.

ASK: is there a g2 at which a record both MOVES and LASTS?

TWO DEFECTS IN v1, BOTH FIXED HERE.
 (1) v1 read <R>(T) at ONE time. Once g2>0, R does not commute with H, so <R> PRECESSES --
     the gamma=0 control gave 0.587 with no bath at all. A snapshot samples a rotating
     quantity, not a survival envelope. FIX: time-average over the back half of the run,
     and run a gamma=0 control at EVERY g2 so oscillation is separated from decay.
 (2) v1 called g2=5.0 an "interior optimum" when it was the last point of the scan and
     still rising. FIX: report explicitly whether the peak is interior or at the edge.
INTEGRATOR REMOVED. The Lindblad generator is a D^2 x D^2 linear superoperator (256x256
here), so it is exponentiated EXACTLY. No Heun, no step-size error, no purity drift --
this also retires the known w27b_open integrator defect for this lane.
"""
import itertools, numpy as np

def expm(A):
    """Scaling-and-squaring Taylor exponential. scipy is not available here.
    Verified against fine-step integration in CONTROL D below."""
    nrm=np.linalg.norm(A,np.inf)
    n=max(0,int(np.ceil(np.log2(nrm)))+1) if nrm>0 else 0
    B=A/(2.0**n); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for k in range(1,40):
        T=T@B/k; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(n): X=X@X
    return X


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
N=2; st,idx=build(V2,E,N); D=len(st)
R=Move(st,idx,compose(P),N)
Lmag=[Move(st,idx,p,N) for p in P]
Lel=[Zop(st,[k],N) for k in CUT]
I=np.eye(D,dtype=complex); Pplus=(I+R)/2.0
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
MAG =sum(L+L.conj().T for L in Lmag)

def superop(H,Ls,gam):
    Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    for L in Ls: M+=gam*(np.kron(L.conj(),L)-np.kron(Id,Id))
    return M
def track(H,Ls,gam,ts,rho0):
    M=superop(H,Ls,gam); v0=rho0.reshape(-1)
    out=[]
    for t in ts:
        v=expm(M*t)@v0; r=v.reshape(D,D); r=r/np.trace(r).real
        out.append(np.trace(R@r).real)
    return np.array(out)

ts=np.linspace(0,20,81); half=ts>=10.0
print(f"W-30b v2   3x3 patch Z_2, physical dim {D}. R = rim loop. EXACT Lindblad exponential.")
print("           bath = electric on the cut links (the W-29b conjugate), gamma=0.5.")
print("           LASTS = time-average of <R> over t in [10,20]; OSC = same with gamma=0.")
print()
print(f"  {'g^2':>7s} {'MOVES=||[H,R]||/dim':>20s} {'LASTS(gam=.5)':>14s} {'OSC(gam=0)':>12s} "
      f"{'LASTS-OSC':>11s} {'MOVESxLASTS':>12s}")
print("  "+"-"*82)
rows=[]
for g2 in [0.0,0.01,0.03,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0,5.0,8.0,12.0,20.0]:
    H=-MAG-g2*ELEC
    mv=np.linalg.norm(H@R-R@H)/D
    rho0=(Pplus/np.trace(Pplus).real).astype(complex)
    lasts=track(H,Lel,0.5,ts,rho0)[half].mean()
    osc  =track(H,Lel,0.0,ts,rho0)[half].mean()
    rows.append((g2,mv,lasts,osc,mv*abs(lasts)))
    print(f"  {g2:7.2f} {mv:20.4f} {lasts:14.6f} {osc:12.6f} {lasts-osc:11.6f} {mv*abs(lasts):12.6f}")

vals=[r[4] for r in rows]; k=int(np.argmax(vals))
edge = (k==0 or k==len(rows)-1)
print()
print(f"  peak of MOVESxLASTS at g^2 = {rows[k][0]}  value {vals[k]:.6f}  "
      f"-> {'AT THE EDGE OF THE SCAN, not an optimum' if edge else 'INTERIOR OPTIMUM'}")
print()
print("  CONTROL A -- gamma=0 at g2=0 must give exactly 1 (R is conserved, nothing can move it).")
H0=-MAG
print(f"    <R> over [10,20], g2=0, gamma=0   : {track(H0,Lel,0.0,ts,(Pplus/np.trace(Pplus).real).astype(complex))[half].mean():.12f}")
print(f"    <R> over [10,20], g2=0, gamma=0.5 : {track(H0,Lel,0.5,ts,(Pplus/np.trace(Pplus).real).astype(complex))[half].mean():.12f}")
print()
print("  CONTROL B -- magnetic bath (commutes with R) at several g2.")
for g2 in [0.0,0.5,2.0,20.0]:
    H=-MAG-g2*ELEC
    print(f"    g2={g2:6.2f} magnetic bath  <R>[10,20] = "
          f"{track(H,Lmag,0.5,ts,(Pplus/np.trace(Pplus).real).astype(complex))[half].mean():.9f}")
print()
print("  CONTROL C -- large-g2 diagnostic: does the electric bath still act, or has it become")
print("               a symmetry of H?  ||[H, L_el]|| tells us. If it -> 0, the bath stops working.")
for g2 in [0.5,2.0,5.0,20.0]:
    H=-MAG-g2*ELEC
    print(f"    g2={g2:6.2f}  max||[H,L_el]|| = {max(np.linalg.norm(H@L-L@H) for L in Lel):10.4f}"
          f"   ||[H,R]||/dim = {np.linalg.norm(H@R-R@H)/D:9.4f}")

print()
print("  CONTROL D -- is the exponential right? Compare against fine-step RK integration.")
def lind_rhs(rho,H,Ls,gam):
    return -1j*(H@rho-rho@H)+gam*sum(L@rho@L.conj().T-rho for L in Ls)
for g2 in [0.0,0.5,2.0]:
    H=-MAG-g2*ELEC
    rho=(Pplus/np.trace(Pplus).real).astype(complex)
    T=5.0; steps=20000; dt=T/steps
    r=rho.copy()
    for _ in range(steps):
        k1=lind_rhs(r,H,Lel,0.5); k2=lind_rhs(r+0.5*dt*k1,H,Lel,0.5)
        k3=lind_rhs(r+0.5*dt*k2,H,Lel,0.5); k4=lind_rhs(r+dt*k3,H,Lel,0.5)
        r=r+dt*(k1+2*k2+2*k3+k4)/6.0
    rk=np.trace(R@(r/np.trace(r).real)).real
    ex=track(H,Lel,0.5,[T],rho)[0]
    print(f"    g2={g2:5.2f}  RK4(20k steps)   <R>(5) = {rk:.9f}   exact expm = {ex:.9f}   diff = {abs(rk-ex):.2e}")
