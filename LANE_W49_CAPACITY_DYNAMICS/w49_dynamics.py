"""W-49.  GIVE CAPACITY DYNAMICS. DOES THE CARRIER GROW, OR ONLY EVICT?

W-41: on a FIXED carrier, demand beyond capacity is resolved by EVICTION. W-42: capacity = area - 1.
W-45: capacity passes all four of gravity's functional markers -- but it is static bookkeeping with
no equation of motion. This lane gives it one.

THE TRAP, AND IT IS W-31's. If the rule is "grow when full", the answer is installed and three
adversaries will (rightly) refute it. So the structural variable is gated by COST ALONE and never
mentions the records:
    H = -(W0 + h.c.) - (W1 + h.c.) - Nhat (x) (W2 + h.c.) - g2 sum Z + mu*Nhat - Delta*sigma_x
The occupation Nhat decides whether the third plaquette participates. It carries an energy cost mu
and tunnels at Delta. NOTHING in H references a record, a boundary, or a capacity.
Capacity is then 1 when n=0 (two active plaquettes) and 2 when n=1 (three).

THE FORCED-OR-NOT CHECK, AND IT IS DECISIVE, SO IT RUNS FIRST.
If the Lindbladian has a UNIQUE steady state then <Nhat> at t=infinity is independent of EVERY
feature of the initial state, so "does the carrier grow when more records are present" would be
void by construction -- the same trap W-33 fell into. Count the zero modes BEFORE measuring anything.
If unique, the question must be asked as a CORRELATION inside the steady state, or at finite time,
and the lane says which it used.

THE MEASUREMENT. Load k independent records (definite plaquette fluxes) and ask whether the carrier's
size responds. Reported both ways: <Nhat> versus k, and the connected correlation
<Nhat R_i> - <Nhat><R_i>, which is the content-sources-geometry statement.
"""
import itertools, numpy as np

# --- carrier: a 4x2 vertex ladder -> 8 vertices, 10 links, 3 plaquettes, cycle rank 3 ---
nx,ny=4,2
V=[(i,j) for j in range(ny) for i in range(nx)]; vid={v:k for k,v in enumerate(V)}
E=[]
for j in range(ny):
    for i in range(nx-1): E.append((vid[(i,j)],vid[(i+1,j)]))
NH=len(E)
for j in range(ny-1):
    for i in range(nx): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j: j*(nx-1)+i
vx =lambda i,j: NH + j*nx + i
PL=[[(hid(i,0),+1),(vx(i+1,0),+1),(hid(i,1),-1),(vx(i,0),-1)] for i in range(nx-1)]
st=[s for s in itertools.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(len(V)))]
idx={s:i for i,s in enumerate(st)}; DG=len(st)
print(f"W-49  carrier {nx}x{ny}: {len(V)} vertices, {L} links, {len(PL)} plaquettes, "
      f"cycle rank {L-len(V)+1}, gauge dim {DG}")

def Move(mv):
    M=np.zeros((DG,DG),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%2
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def Zop(links):
    return np.diag([(-1.0)**(sum(s[k] for k in links)%2) for s in st]).astype(complex)
W=[Move(p) for p in PL]
ELEC=sum(Zop([k]) for k in range(L))
IG=np.eye(DG,dtype=complex); I2=np.eye(2,dtype=complex)
sx=np.array([[0,1],[1,0]],complex); Nh=np.array([[0,0],[0,1]],complex)
D=DG*2
def K(a,b): return np.kron(a,b)

def H_of(mu,Delta,g2):
    H =K(-(W[0]+W[0].conj().T),I2)+K(-(W[1]+W[1].conj().T),I2)
    H+=K(-(W[2]+W[2].conj().T),Nh)              # third plaquette participates ONLY if n=1
    H+=K(-g2*ELEC,I2)+K(IG,mu*Nh)-Delta*K(IG,sx)
    return H
def liou(mu,Delta,g2,bath,gam=0.5):
    H=H_of(mu,Delta,g2); Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for k in bath:
        X=K(Zop([k]),I2); M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    return M

BATH=[0]
print()
print("  FORCED-OR-NOT, FIRST. A UNIQUE steady state makes <Nhat>(inf) independent of the initial")
print("  state, so the whole question would be void by construction (the W-33 trap).")
print(f"  {'mu':>6s} {'Delta':>7s} {'g2':>6s} {'#zero modes':>12s} {'next rate':>11s}")
for mu,Delta,g2 in [(0.5,0.3,0.05),(1.5,0.3,0.05),(0.5,0.0,0.05),(0.5,0.3,0.00)]:
    ev=np.linalg.eigvals(liou(mu,Delta,g2,BATH)); rate=-ev.real
    nz=int((rate<1e-9).sum()); nxt=rate[rate>=1e-9].min() if (rate>=1e-9).any() else float('nan')
    print(f"  {mu:6.2f} {Delta:7.2f} {g2:6.2f} {nz:12d} {nxt:11.3e}"
          f"   {'UNIQUE -> <N>(inf) cannot depend on k' if nz==1 else 'DEGENERATE'}")

print()
print("  Because of that, the question is asked as a CORRELATION inside the steady state, and also")
print("  at finite time. Neither is void when the steady state is unique.")

def steady(M):
    w,V_=np.linalg.eig(M); i=np.argmin(np.abs(w))
    r=V_[:,i].reshape(D,D); r=(r+r.conj().T)/2
    tr=np.trace(r).real
    return r/tr if abs(tr)>1e-12 else r

Rops=[K(W[i]+W[i].conj().T,I2)/2.0 for i in range(3)]   # record observables (Hermitian part)
NHAT=K(IG,Nh)
print()
print(f"  {'mu':>6s} {'Delta':>7s} {'g2':>6s} {'<Nhat>':>9s} "
      f"{'corr(N,R0)':>11s} {'corr(N,R1)':>11s} {'corr(N,R2)':>11s}")
print("  "+"-"*70)
for mu in (0.2,0.5,1.0,2.0):
    for Delta in (0.3,):
        for g2 in (0.05,):
            r=steady(liou(mu,Delta,g2,BATH))
            n=np.trace(NHAT@r).real
            cs=[]
            for Ri in Rops:
                c=np.trace(NHAT@Ri@r).real-n*np.trace(Ri@r).real
                cs.append(c)
            print(f"  {mu:6.2f} {Delta:7.2f} {g2:6.2f} {n:9.5f} "
                  f"{cs[0]:11.3e} {cs[1]:11.3e} {cs[2]:11.3e}")

print()
print("  CONTROL -- Delta = 0 freezes the structure; <Nhat> must stay at its initial value and all")
print("  correlations must be trivial.")
for mu in (0.2,2.0):
    r=steady(liou(mu,0.0,0.05,BATH))
    n=np.trace(NHAT@r).real
    print(f"    mu={mu:4.2f} Delta=0   <Nhat> = {n:.6f}   corr(N,R0) = "
          f"{np.trace(NHAT@Rops[0]@r).real-n*np.trace(Rops[0]@r).real:.3e}")
print()
print("  CONTROL -- does mu actually control the structure? <Nhat> must fall as mu rises.")
for mu in (0.0,0.5,1.0,2.0,4.0,8.0):
    r=steady(liou(mu,0.3,0.05,BATH))
    print(f"    mu={mu:4.1f}   <Nhat> = {np.trace(NHAT@r).real:.6f}")
