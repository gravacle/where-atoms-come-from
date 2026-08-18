"""W-34.  THE PREDICTABILITY SIEVE -- LET THE DYNAMICS NAME THE RECORD.

Every measurement so far NOMINATED the rim Wilson loop as the record and then asked how it fares.
That installs the answer. Zurek's point is the opposite: the pointer observable is SELECTED by the
dynamics as the one the environment damages least. So do not nominate anything.

The left-eigenmodes of the Liouvillian ARE the observables with definite decay rates:
    d<O_k>/dt = -Gamma_k <O_k>.
So diagonalise the Liouvillian, sort by decay rate, and ASK WHAT THE SLOWEST OBSERVABLES ARE.
If the rim loop comes out on top, it was selected. If something else does, that something else is
what this construction actually records -- and nominating the rim was an imported expectation.

IDENTIFICATION. Physical dim 16 = 2^4: one effective qubit per independent cycle. The 16 magnetic
operators are the products of subsets of the 4 plaquettes (the empty subset is the identity; the
full subset is the RIM LOOP by discrete Stokes). Electric character is measured as the weight of
the mode that is diagonal in the electric basis. Each slow mode is reported with its best match.
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
Lmag=[Move(st,idx,p,N) for p in P]
MAG=sum(L+L.conj().T for L in Lmag)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
Id=np.eye(D,dtype=complex)

NAMES={}; OPS={}
for r in range(5):
    for S in itertools.combinations(range(4),r):
        O=Move(st,idx,compose([P[i] for i in S]),N) if S else Id.copy()
        nm = "identity" if not S else ("RIM LOOP (all 4)" if r==4 else
             f"{'plaquette' if r==1 else str(r)+'-plaquette loop'} {S}")
        NAMES[S]=nm; OPS[S]=O/np.linalg.norm(O)

def report(g2,links,gam=0.5,topk=6):
    H=-MAG-g2*ELEC
    Ls=[Zop(st,[k],N) for k in links]
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    for L in Ls: M+=gam*(np.kron(L.conj(),L)-np.kron(Id,Id))
    w,U=np.linalg.eig(M.conj().T)          # left modes = observables
    rate=-np.conj(w).real
    order=np.argsort(rate)
    print(f"\n  g2={g2}  bath on {'CUT' if links is CUT else 'RIM'}  gamma={gam}")
    print(f"    {'rate':>12s} {'|diag wt|':>10s}  best operator match")
    print("    "+"-"*62)
    shown=0
    for i in order:
        if shown>=topk: break
        O=U[:,i].reshape(D,D)
        n=np.linalg.norm(O)
        if n<1e-12: continue
        O=O/n
        dw=np.linalg.norm(np.diag(np.diag(O)))     # electric-type weight
        best=max(OPS.items(), key=lambda kv: abs(np.vdot(kv[1].reshape(-1),O.reshape(-1))))
        ov=abs(np.vdot(best[1].reshape(-1),O.reshape(-1)))
        tag=NAMES[best[0]] if ov>0.3 else ("electric-type (diagonal)" if dw>0.7 else "mixed / no clean match")
        print(f"    {rate[i]:12.6e} {dw:10.3f}  {tag}   (overlap {ov:.3f})")
        shown+=1

print("W-34  PREDICTABILITY SIEVE. Slowest observables = what the dynamics actually protects.")
for g2 in [0.0,0.001,0.01,0.1]:
    report(g2,CUT)
report(0.0,PERIM); report(0.01,PERIM)

print()
print("  CROSS-CHECK -- rank the 16 magnetic operators by their OWN decay rate, bath on the cut,")
print("  g2=0.01. If the rim loop is genuinely the pointer it must be the slowest of the sixteen.")
H=-MAG-0.01*ELEC; Ls=[Zop(st,[k],N) for k in CUT]
M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
for L in Ls: M+=0.5*(np.kron(L.conj(),L)-np.kron(Id,Id))
w,U=np.linalg.eig(M.conj().T); rate=-np.conj(w).real
Uinv_cols=U/np.linalg.norm(U,axis=0)
rows=[]
for S,O in OPS.items():
    if not S: continue
    ov=np.abs(Uinv_cols.conj().T@O.reshape(-1))
    ov=ov/max(ov.sum(),1e-30)
    rows.append((float((ov*rate).sum()), NAMES[S], len(S)))
for r,nm,k in sorted(rows):
    print(f"    mean decay rate {r:12.6e}   {nm}")
