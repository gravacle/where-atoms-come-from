"""W-35.  WAS THE INSTRUMENT EVER ABLE TO SEE A RECORD?

Every lane from W-28 to W-34 measured <R> = Tr(R rho). That is an ENSEMBLE AVERAGE. A fair coin has
mean zero and every single flip is definite. So d<R>/dt = 0 is consistent with every individual
run acquiring a perfectly definite record -- and W-29c's 'exactly 0.000000' would then be measuring
the FAIRNESS of the coin, not the absence of a flip.

THE TEST. Unravel the same Lindbladian into individual trajectories and measure the SECOND moment:
    E[<R>_psi]     must stay 0        (the ensemble result, reproduced -- the coin is fair)
    E[<R>_psi ^2]  -> 1 would mean    EVERY RUN ACQUIRES A DEFINITE RECORD
Same generator, same rho(t). If the second moment moves while the first does not, then no ensemble
measurement in this program could ever have detected record formation, in either direction.

CONTROLS THAT MUST FIRE.
  (a) gamma = 0: nothing is monitored, E[<R>^2] must stay at its initial value.
  (b) rho(t) reconstructed from the trajectories must match the EXACT master-equation solution --
      otherwise the unravelling is not of this Lindbladian and proves nothing.
  (c) the JUMP unravelling of the SAME Lindbladian: with unitary L commuting with R every jump
      leaves <R>_psi untouched, so it must give E[<R>^2] = const. If jump and homodyne disagree,
      record formation is a fact about HOW THE ENVIRONMENT IS READ, not about the generator.
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
N=2; st,idx=build(V2,E,N); D=len(st)
R=Move(st,idx,compose(P),N)
MAG=sum((lambda L:L+L.conj().T)(Move(st,idx,p,N)) for p in P)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
Id=np.eye(D,dtype=complex)
Ls=[Zop(st,[k],N) for k in CUT]            # Hermitian AND unitary at N=2
for L in Ls: assert np.linalg.norm(L@L.conj().T-Id)<1e-12 and np.linalg.norm(L-L.conj().T)<1e-12

Pp=(Id+R)/2; Pm=(Id-R)/2

def expm(A):
    nrm=np.linalg.norm(A,np.inf)
    n=max(0,int(np.ceil(np.log2(nrm)))+1) if nrm>0 else 0
    B=A/(2.0**n); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for k in range(1,40):
        T=T@B/k; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(n): X=X@X
    return X

def init_ens(NT,seed=7):
    """NT unbiased pure states, <R> = 0 EXACTLY on every one (equal weight in both sectors)."""
    g=np.random.default_rng(seed)
    W=g.normal(size=(NT,D))+1j*g.normal(size=(NT,D))
    A=W@Pp.T; B=W@Pm.T
    A=A/np.linalg.norm(A,axis=1,keepdims=True); B=B/np.linalg.norm(B,axis=1,keepdims=True)
    Psi=(A+B)/np.sqrt(2.0)
    return Psi/np.linalg.norm(Psi,axis=1,keepdims=True)

def expect(Psi,O):
    return np.real(np.einsum('ij,ij->i',Psi.conj(),Psi@O.T))

def homodyne(g2,gam,T,NT,dt=0.004,seed=11):
    """SPLIT-STEP. The Hamiltonian part is propagated EXACTLY by expm(-iH dt); only the measurement
    part is stochastic. A plain Euler step on -iH psi is NOT norm-preserving: after renormalisation
    it systematically amplifies extremal eigenvectors of H, and those carry definite R, so it
    MANUFACTURES a second moment out of nothing. The gamma=0 control below detects exactly that."""
    g=np.random.default_rng(seed); H=-MAG-g2*ELEC; steps=int(round(T/dt))
    U=expm(-1j*H*dt)
    Psi=init_ens(NT)
    for _ in range(steps):
        Psi=Psi@U.T                                  # exact, unitary, no drift
        if gam>0:
            dP=np.zeros_like(Psi)
            for L in Ls:
                LP=Psi@L.T
                m=np.real(np.einsum('ij,ij->i',Psi.conj(),LP))[:,None]
                dev=LP-m*Psi
                dP=dP-0.5*gam*(dev@L.T-m*dev)*dt+np.sqrt(gam)*dev*g.normal(size=(NT,1))*np.sqrt(dt)
            Psi=Psi+dP
        Psi=Psi/np.linalg.norm(Psi,axis=1,keepdims=True)
    v=expect(Psi,R)
    return v,np.einsum('ij,ik->jk',Psi,Psi.conj())/NT

def jump(g2,gam,T,NT,dt=0.004,seed=23):
    g=np.random.default_rng(seed); H=-MAG-g2*ELEC; steps=int(round(T/dt))
    Psi=init_ens(NT); U=expm(-1j*H*dt)
    for _ in range(steps):
        Psi=Psi@U.T
        fire=g.random(NT)<gam*len(Ls)*dt
        which=g.integers(len(Ls),size=NT)
        for k,L in enumerate(Ls):
            m=fire&(which==k)
            if m.any(): Psi[m]=Psi[m]@L.T
        Psi=Psi/np.linalg.norm(Psi,axis=1,keepdims=True)
    return expect(Psi,R)

def exact_rho(g2,gam,T,rho0):
    H=-MAG-g2*ELEC
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    for L in Ls: M+=gam*(np.kron(L.conj(),L)-np.kron(Id,Id))
    return (expm(M*T)@rho0.reshape(-1)).reshape(D,D)

print("W-35  DOES A RECORD FORM ON INDIVIDUAL TRAJECTORIES?")
print()
print("  GATE 0 -- the integrator must not manufacture the effect. gamma=0 is pure unitary evolution")
print("  and [H,R]=0 at g2=0, so E[<R>^2] MUST be 0. If it is not, nothing below means anything.")
for g2 in [0.0]:
    for dt_ in [0.008,0.004,0.002]:
        v,_=homodyne(g2,0.0,15.0,400,dt=dt_)
        st_="PASS" if abs((v*v).mean())<1e-12 else "FAIL -- integrator is inventing a record"
        print(f"    g2={g2}  gamma=0  dt={dt_:6.4f}   E[<R>^2] = {(v*v).mean():.3e}   {st_}")
print()
print("  GATE 1 -- step-size convergence of the real effect (g2=0, gamma=0.3).")
for dt_ in [0.008,0.004,0.002,0.001]:
    v,_=homodyne(0.0,0.3,15.0,400,dt=dt_)
    print(f"    dt={dt_:6.4f}   E[<R>] = {v.mean():+.5f}   E[<R>^2] = {(v*v).mean():.6f}")

print(f"      dim={D}, bath on the CUT, max||[L,R]|| = {max(np.linalg.norm(L@R-R@L) for L in Ls):.1e}")
print("      every initial state has <R> = 0 EXACTLY, so any second moment is acquired, not present.")
NT=600; T=15.0
print()
print(f"  HOMODYNE unravelling  ({NT} trajectories, T={T})")
print(f"  {'g^2':>7s} {'gamma':>6s} {'E[<R>]':>11s} {'+-SE':>8s} {'E[<R>^2]':>11s} {'|<R>|>0.9':>11s} {'|<R>|>0.99':>11s}")
print("  "+"-"*72)
for g2,gam in [(0.0,0.0),(0.0,0.1),(0.0,0.3),(0.0,1.0),(0.01,0.3),(0.1,0.3),(1.0,0.3)]:
    v,_=homodyne(g2,gam,T,NT)
    print(f"  {g2:7.3f} {gam:6.2f} {v.mean():11.6f} {v.std()/np.sqrt(NT):8.4f} {(v*v).mean():11.6f} "
          f"{(np.abs(v)>0.9).mean():11.3f} {(np.abs(v)>0.99).mean():11.3f}")

print()
print(f"  JUMP unravelling of the SAME Lindbladian  ({NT} trajectories)")
print(f"  {'g^2':>7s} {'gamma':>6s} {'E[<R>]':>11s} {'E[<R>^2]':>11s} {'|<R>|>0.9':>11s}")
print("  "+"-"*54)
for g2,gam in [(0.0,0.3),(0.0,1.0),(0.1,0.3)]:
    v=jump(g2,gam,T,NT)
    print(f"  {g2:7.3f} {gam:6.2f} {v.mean():11.6f} {(v*v).mean():11.6f} {(np.abs(v)>0.9).mean():11.3f}")

print()
print("  CONTROL -- do the trajectories reproduce the EXACT master equation? (same initial ensemble)")
rho0=np.einsum('ij,ik->jk',init_ens(NT),init_ens(NT).conj())/NT
for g2,gam in [(0.0,0.3),(0.1,0.3)]:
    v,rt=homodyne(g2,gam,T,NT)
    ex=exact_rho(g2,gam,T,rho0)
    print(f"    g2={g2:5.2f} gamma={gam}:  ||rho_traj-rho_exact|| = {np.linalg.norm(rt-ex):.4f}"
          f"   <R>_exact = {np.trace(R@ex).real:+.8f}   <R>_traj = {np.trace(R@rt).real:+.6f}")
print()
print("  CONTROL -- growth of E[<R>^2] with time, g2=0, gamma=0.3 (0 at t=0 by construction).")
for T2 in [0.0,1.0,3.0,6.0,10.0,15.0,25.0]:
    v,_=homodyne(0.0,0.3,T2,NT) if T2>0 else (expect(init_ens(NT),R),None)
    print(f"    T={T2:5.1f}   E[<R>] = {v.mean():+.6f}   E[<R>^2] = {(v*v).mean():.6f}")
