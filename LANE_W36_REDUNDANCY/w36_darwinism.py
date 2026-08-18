"""W-36.  IS THE RECORD OBJECTIVE?  REDUNDANCY, WITH AN EXPLICIT ENVIRONMENT.

W-35 showed the record forms on every trajectory under homodyne monitoring and never under jump
monitoring -- same generator, same rho(t). Left there, record formation would be a fact about how
somebody chooses to read the environment, which is not good enough.

Zurek's criterion settles it WITHOUT choosing an unravelling: a record is OBJECTIVE when MANY
DISJOINT FRAGMENTS of the environment EACH carry a full copy, so every observer who interrogates
any fragment gets the same value. That is a property of the global system+environment state alone.

So stop tracing the environment out. Build it.
  SYSTEM      : the 3x3 patch, physical dim 16. Record R = rim loop (which W-34 selected).
  ENVIRONMENT : n qubits, each coupling to ONE cut link:  H_int = kappa * sum_k Z_cut(k) (x) sz_k.
                The coupling references the CUT, never the rim, and never R.
  MEASURE     : the Holevo information each fragment holds about R,
                  I(R:F) = S(sum_r p_r rho_F^r) - sum_r p_r S(rho_F^r)   [bits]
                bounded by H(R) = 1 bit. Averaged over ALL fragments of each size.

WHAT THE SHAPES MEAN, decided before running:
  I(R:F) rising to ~1 bit at SMALL f and then FLAT  -> redundancy: many independent copies, objective.
  I(R:F) rising only linearly, reaching 1 bit at f=n -> no redundancy: one copy, spread out, not objective.
CONTROLS: kappa=0 must give exactly 0 at every f; I(R:F) can never exceed 1 bit; f=n must equal the
total information; and the same run with the coupling moved ONTO the rim (where it can disturb R).
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
N=2; st,idx=build(V2,E,N); DS=len(st)
R=Move(st,idx,compose(P),N)
MAG=sum((lambda L:L+L.conj().T)(Move(st,idx,p,N)) for p in P)
IdS=np.eye(DS,dtype=complex)
Pp=(IdS+R)/2; Pm=(IdS-R)/2

def expm(A):
    nrm=np.linalg.norm(A,np.inf)
    n=max(0,int(np.ceil(np.log2(nrm)))+1) if nrm>0 else 0
    B=A/(2.0**n); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for k in range(1,40):
        T=T@B/k; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(n): X=X@X
    return X

sz=np.array([[1,0],[0,-1]],complex)
def kron_list(ops):
    out=np.array([[1]],complex)
    for o in ops: out=np.kron(out,o)
    return out

SPREAD=True
def run(nq, links, kappa, T, seed=3):
    """System (x) nq environment qubits, exact unitary evolution of the JOINT state."""
    DE=2**nq; DT=DS*DE
    H=np.kron(-MAG,np.eye(DE,dtype=complex))
    for k in range(nq):
        Zk=Zop(st,[links[k%len(links)]],N) if SPREAD else Zop(st,[links[0]],N)
        ops=[np.eye(2,dtype=complex)]*nq; ops[k]=sz
        H=H+kappa*np.kron(Zk,kron_list(ops))
    U=expm(-1j*H*T)
    g=np.random.default_rng(seed)
    w=g.normal(size=DS)+1j*g.normal(size=DS)
    a=Pp@w; b=Pm@w; a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    psiS=(a+b)/np.sqrt(2.0); psiS/=np.linalg.norm(psiS)     # <R> = 0 exactly
    plus=np.ones(2,complex)/np.sqrt(2.0)
    psiE=kron_list([plus.reshape(2,1)]*nq).reshape(-1) if nq else np.array([1.0+0j])
    psi=np.kron(psiS,psiE); psi=U@psi
    return psi.reshape((DS,)+(2,)*nq), DE

def vn(rho):
    ev=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    ev=ev[ev>1e-12]
    return float(-(ev*np.log2(ev)).sum())

def holevo(psiT, nq, frag):
    """I(R:F) in bits: Holevo quantity of the ensemble {p_r, rho_F^r} indexed by the R sector."""
    branches=[]
    for Proj in (Pp,Pm):
        v=np.tensordot(Proj,psiT,axes=([1],[0]))
        p=float(np.vdot(v,v).real)
        if p<1e-14: branches.append((0.0,None)); continue
        v=v/np.sqrt(p)
        keep=[0]+[1+i for i in frag]
        tr=[ax for ax in range(1+nq) if ax not in keep]
        M=np.transpose(v,keep+tr)
        d=int(np.prod([M.shape[i] for i in range(len(keep))])); 
        M=M.reshape(d,-1)
        rho=M@M.conj().T
        # trace out the SYSTEM factor, keep only the fragment
        rho=rho.reshape(DS,2**len(frag),DS,2**len(frag))
        rhoF=np.einsum('ijik->jk',rho)
        branches.append((p,rhoF))
    ps=[b[0] for b in branches]
    avg=sum(p*b for p,b in branches if b is not None)
    return vn(avg)-sum(p*vn(b) for p,b in branches if b is not None), ps

def profile(nq, links, kappa, T, tag):
    psiT,DE=run(nq,links,kappa,T)
    print(f"\n  {tag}   (n={nq} env qubits, kappa={kappa}, T={T})")
    print(f"    {'|F|':>4s} {'I(R:F) bits':>13s} {'#fragments avgd':>16s}")
    print("    "+"-"*38)
    vals=[]
    for f in range(nq+1):
        combos=list(itertools.combinations(range(nq),f))
        if len(combos)>20: combos=combos[:20]
        I=np.mean([holevo(psiT,nq,c)[0] for c in combos])
        vals.append(I)
        print(f"    {f:4d} {I:13.6f} {len(combos):16d}")
    return vals

print("W-36  IS THE RECORD OBJECTIVE?  Holevo information about R held by environment fragments.")
print(f"      system dim {DS}; coupling touches the CUT only; H(R) = 1 bit is the ceiling.")
print()
print("  SATURATION FIRST. A partial-information plot means nothing until the WHOLE environment")
print("  holds the full bit. Scan the coupling until I(R:all) saturates, then read the shape.")
print(f"    {'kappa':>7s} {'T':>6s} {'I(R:all 6) bits':>17s}")
print("    "+"-"*32)
best=None
for kap,T in [(0.6,6.0),(1.2,6.0),(2.0,8.0),(3.0,12.0),(5.0,12.0),(8.0,16.0)]:
    psiT,_=run(6,CUT,kap,T)
    tot=holevo(psiT,6,tuple(range(6)))[0]
    print(f"    {kap:7.2f} {T:6.1f} {tot:17.6f}")
    if best is None or tot>best[2]: best=(kap,T,tot)
print(f"    -> using kappa={best[0]}, T={best[1]} (total {best[2]:.4f} bits)")
KAP,TT=best[0],best[1]

v_cut = profile(6, CUT, KAP, TT, "SPREAD: one qubit per cut link -- [L,R]=0, the QND case")
v_off = profile(6, CUT, 0.0, TT, "CONTROL kappa=0 -- must be exactly 0 at every |F|")
SPREAD=False
v_same = profile(6, CUT, KAP, TT, "SAME LINK: all 6 qubits read ONE cut link (Zurek's setup)")
SPREAD=True
v_rim = profile(6, PERIM, KAP, TT, "CONTROL on the RIM -- environment can disturb R")

print()
print("  READING THE SHAPE. Redundancy = a SMALL fragment already holds nearly the whole bit.")
print(f"  {'case':>10s} {'I(|F|=1)':>10s} {'I(|F|=3)':>10s} {'I(all)':>9s} {'I(1)/I(all)':>12s}  verdict")
print("  "+"-"*76)
for nm,v in [("spread",v_cut),("same-link",v_same),("rim",v_rim)]:
    tot=v[-1]
    if tot<1e-9:
        print(f"  {nm:>10s} {'--':>10s} {'--':>10s} {tot:9.4f} {'--':>12s}  no information anywhere")
        continue
    fr=v[1]/tot
    verdict=("REDUNDANT: one fragment nearly suffices -> objective" if fr>0.5 else
             "PARTIAL redundancy" if fr>0.2 else
             "NO redundancy: information is delocalised, only the whole environment knows")
    print(f"  {nm:>10s} {v[1]:10.4f} {v[3]:10.4f} {tot:9.4f} {fr:12.3f}  {verdict}")
