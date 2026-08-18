"""W-52.  DOES A RECORD MAKE ANOTHER RECORD MORE LIKELY NEARBY?

Map step 3, in its corrected form. The FIRST version of that step asked whether records ATTRACT,
which imports the classical footprint. What attraction DOES in the classical theory is supply
POSITIVE FEEDBACK -- content making more content more likely in the same place. Ask that instead,
without assuming a sign or a mechanism.

ALL THREE OUTCOMES ARE READ IN ADVANCE AND NONE IS A TEST ANYTHING MUST PASS:
  positive feedback  content concentrates
  negative feedback  what W-39/W-41 already measured (crowding, then eviction). NOT a refutation.
  neither            records indifferent to each other -- the strongest of the three, since it would
                     mean the resource has no dynamics of its own.

CONSTRUCTION. 3x3 patch, Z_2, physical dim 16, four plaquettes in a 2x2 arrangement.
  plaquettes 0,1,2,3.  0-1, 0-2, 1-3, 2-3 SHARE a link.  0-3 and 1-2 are diagonal and share NONE.
So plaquette 1 is a NEIGHBOUR of 0 and plaquette 3 is NOT. That is the only asymmetry in the test.
Record formation is measured as in W-35: homodyne trajectories, E[<R>^2] rising from 0.

THE COMPARISON:
  (a) no record loaded          -> baseline formation on plaquette 1 and on plaquette 3
  (b) a record loaded on 0      -> formation on 1 (neighbour) and on 3 (non-neighbour)
If loading 0 helps 1 more than 3, that is spatially structured positive feedback.
If it hurts, that is crowding, already known.
If neither moves, records are indifferent.

GATE (the W-35 lesson): with gamma = 0 the evolution is unitary and no second moment may be
manufactured. Checked before anything else is read.
"""
import itertools, numpy as np

V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
CENTER=vid[(1,1)]; CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
st=[s for s in itertools.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(9))]
idx={s:i for i,s in enumerate(st)}; D=len(st)
def Move(mv):
    M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%2
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def Zop(links): return np.diag([(-1.0)**(sum(s[k] for k in links)%2) for s in st]).astype(complex)
W=[Move(p) for p in P]
sup=[set(k for k,_ in p) for p in P]
print("W-52  DOES A RECORD HELP ANOTHER FORM NEARBY?")
print(f"      physical dim {D}; plaquette link supports:")
for i in range(4): print(f"        p{i}: {sorted(sup[i])}")
print(f"      shared links: " + ", ".join(f"p{i}-p{j}:{len(sup[i]&sup[j])}"
      for i in range(4) for j in range(i+1,4)))
NEI, FAR = 1, 3     # p1 shares a link with p0; p3 shares none
print(f"      -> NEIGHBOUR of p0 = p{NEI}   NON-NEIGHBOUR = p{FAR}")
for i in range(4):
    u=np.linalg.norm(W[i]@W[i].conj().T-np.eye(D))
    print(f"      p{i}: ||W||={np.linalg.norm(W[i]):.3f}  unitarity {u:.1e}  "
          f"eigenvalues {len(np.unique(np.round(np.linalg.eigvals(W[i]),6)))}")

MAG=sum(w+w.conj().T for w in W)
ELEC=sum(Zop([k]) for k in range(L))
# BATH CHOICE IS DECIDED BY THE PROTECTION CONDITION, NOT BY HABIT.
# v1 used the CUT, which shares links 3 and 7 with p1: the bath does not commute with a single
# plaquette, so it DESTROYS that record instead of monitoring it, and formation was ~0 everywhere.
# Only the rim loop is protected by a cut bath (W-34) -- and one protected record cannot answer a
# two-record question. The links avoided by p0, p1 AND p3 together are exactly {4,9}.
BATH=[4,9]
for i in (0,1,3):
    assert not (sup[i] & set(BATH)), f"bath must avoid p{i}"
Ls=[Zop([k]) for k in BATH]
Pp=[(np.eye(D)+ (W[i]+W[i].conj().T)/2)/2 for i in range(4)]
Pm=[(np.eye(D)- (W[i]+W[i].conj().T)/2)/2 for i in range(4)]
Rop=[(W[i]+W[i].conj().T)/2 for i in range(4)]

def init(NT,load,g,seed=7):
    """NT pure states with <R_NEI> = <R_FAR> = 0 EXACTLY. If load, R_0 is definite as well.
    v1 balanced NEI and then balanced FAR, which BROKE the first balance -- the gamma=0 gate caught
    it (E[<R_NEI>^2] = 2.8e-3 where it must be 0). The plaquette operators COMMUTE, so the four
    joint sectors (+-,+-) are simultaneously well defined: build an equal superposition over all
    four and both first moments vanish by construction."""
    Wv=g.normal(size=(NT,D))+1j*g.normal(size=(NT,D))
    if load: Wv=Wv@Pp[0].T                      # definite record on plaquette 0
    def nz(X):
        n=np.linalg.norm(X,axis=1,keepdims=True); n[n<1e-12]=1.0; return X/n
    Psi=np.zeros_like(Wv)
    for Pa in (Pp[NEI],Pm[NEI]):
        for Pb in (Pp[FAR],Pm[FAR]):
            Psi=Psi+nz(Wv@(Pa@Pb).T)
    return nz(Psi)

def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X

def run(load,g2,gam,T=15.0,NT=600,dt=0.004,seed=11):
    g=np.random.default_rng(seed)
    H=-MAG-g2*ELEC; U=expm(-1j*H*dt); Psi=init(NT,load,np.random.default_rng(seed+1))
    steps=int(round(T/dt))
    for _ in range(steps):
        Psi=Psi@U.T
        if gam>0:
            dP=np.zeros_like(Psi)
            for Lk in Ls:
                LP=Psi@Lk.T
                m=np.real(np.einsum('ij,ij->i',Psi.conj(),LP))[:,None]
                dev=LP-m*Psi
                dP=dP-0.5*gam*(dev@Lk.T-m*dev)*dt+np.sqrt(gam)*dev*g.normal(size=(NT,1))*np.sqrt(dt)
            Psi=Psi+dP
        Psi=Psi/np.linalg.norm(Psi,axis=1,keepdims=True)
    out={}
    for i in (NEI,FAR,0):
        v=np.real(np.einsum('ij,ij->i',Psi.conj(),Psi@Rop[i].T))
        out[i]=( (v*v).mean(), v.mean() )
    return out

print()
print("  GATE -- gamma = 0 is unitary; no second moment may be manufactured.")
for load in (False,True):
    o=run(load,0.0,0.0)
    print(f"    load={str(load):5s}  E[<R_{NEI}>^2] = {o[NEI][0]:.3e}   E[<R_{FAR}>^2] = {o[FAR][0]:.3e}")
print()
print("  FORMATION.  E[<R>^2] rises from 0 as each record becomes definite on individual runs.")
print(f"  {'g^2':>6s} {'gamma':>6s} {'load p0':>8s} {'E[<R_1>^2] NEIGHBOUR':>21s} {'E[<R_3>^2] FAR':>16s} {'N-F':>9s}")
print("  "+"-"*74)
rows={}
for g2 in (0.0,0.05):
    for gam in (0.3,1.0):
        for load in (False,True):
            o=run(load,g2,gam)
            rows[(g2,gam,load)]=o
            print(f"  {g2:6.2f} {gam:6.2f} {str(load):>8s} {o[NEI][0]:21.6f} {o[FAR][0]:16.6f} "
                  f"{o[NEI][0]-o[FAR][0]:9.4f}")
print()
print("  THE COMPARISON -- does loading p0 change formation on its NEIGHBOUR more than on the FAR one?")
print(f"  {'g^2':>6s} {'gamma':>6s} {'d(neighbour)':>13s} {'d(far)':>9s} {'difference':>11s}   reading")
print("  "+"-"*74)
for g2 in (0.0,0.05):
    for gam in (0.3,1.0):
        dn=rows[(g2,gam,True)][NEI][0]-rows[(g2,gam,False)][NEI][0]
        df=rows[(g2,gam,True)][FAR][0]-rows[(g2,gam,False)][FAR][0]
        d=dn-df
        rd=("positive, spatially structured" if d>0.02 else
            "negative (crowding)" if d<-0.02 else "indifferent")
        print(f"  {g2:6.2f} {gam:6.2f} {dn:13.5f} {df:9.5f} {d:11.5f}   {rd}")
