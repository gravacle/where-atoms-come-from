"""W-33.  THE SELECTION PROBLEM.  W-29c: from a maximally mixed start every physical bath left
<R> at EXACTLY 0.000000. I read that as 'the dynamics does not select'. Test that reading.

THE FIRST QUESTION IS STRUCTURAL, AND IT IS DECISIVE BEFORE ANY DYNAMICS IS RUN.
A Lindbladian with a UNIQUE steady state cannot select anything, ever, from any initial condition:
everything relaxes to the same rho_ss. Selection REQUIRES a degenerate steady-state manifold.
So: count the zero modes of the Liouvillian. That is one eigen-decomposition, and it settles
whether W-29c's null was a fact about the bath or a fact about the whole construction.

THE SECOND QUESTION REFRAMES W-29c.
A maximally mixed initial state carries NO information. A record of nothing is nothing, so <R>=0
is the CORRECT answer there, not a failure. What a record must do is CORRELATE with something.
So drive the initial state through a range of bias and ask what the final <R> remembers:
    PRESERVES   final = initial            -- memory
    AMPLIFIES   |final| > |initial|         -- genuine selection, the tie gets broken
    ERASES      final -> 0                  -- no record
The distinction is read off the slope d<R>(inf)/d<R>(0), not asserted.
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
Lmag=[Move(st,idx,p,N) for p in P]
MAG=sum(L+L.conj().T for L in Lmag)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
Id=np.eye(D,dtype=complex)

def liou(g2,links,gam=0.5,extraLs=None):
    H=-MAG-g2*ELEC
    Ls=[Zop(st,[k],N) for k in links]+(extraLs or [])
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    for L in Ls: M+=gam*(np.kron(L.conj(),L)-np.kron(Id,Id))
    return M

print("W-33  QUESTION 1 -- HOW MANY STEADY STATES? A unique steady state forbids selection outright.")
print(f"  {'g^2':>8s} {'bath':>10s} {'#zero modes':>12s} {'next rate':>13s}   reading")
print("  "+"-"*70)
for links,nm in [(CUT,"cut"),(PERIM,"rim"),(list(range(len(E))),"all links")]:
    for g2 in [0.0,0.001,0.01,0.1,1.0]:
        ev=np.linalg.eigvals(liou(g2,links))
        rate=-ev.real
        nz=int((rate<1e-9).sum())
        nxt=rate[rate>=1e-9].min() if (rate>=1e-9).any() else float('nan')
        rd = "UNIQUE -> selection impossible" if nz==1 else f"DEGENERATE x{nz} -> a record manifold exists"
        print(f"  {g2:8.4f} {nm:>10s} {nz:12d} {nxt:13.3e}   {rd}")

print()
print("W-33  QUESTION 2 -- WHAT DOES THE FINAL STATE REMEMBER?")
print("  Initial state: rho(b) = (I + b*R)/dim, a physical state biased toward the R=+1 sector.")
print("  Evolved to the steady state by projecting onto the Liouvillian's zero modes (t -> inf, exact).")
print()
def steady_from(M, rho0, tol=1e-9):
    w,V=np.linalg.eig(M)
    Vi=np.linalg.inv(V)
    c=Vi@rho0.reshape(-1)
    c[np.abs(w.real)>tol]=0.0            # kill every decaying mode: this IS t -> infinity
    r=(V@c).reshape(D,D)
    tr=np.trace(r).real
    return r/tr if abs(tr)>1e-12 else r

for links,nm in [(CUT,"cut"),(PERIM,"rim")]:
    for g2 in [0.0,0.001,0.01,0.1]:
        M=liou(g2,links)
        xs=[];ys=[]
        for b in [-0.9,-0.5,-0.2,0.0,0.2,0.5,0.9]:
            rho0=(Id+b*R)/D
            r=steady_from(M,rho0.astype(complex))
            xs.append(np.trace(R@rho0).real); ys.append(np.trace(R@r).real)
        sl=np.polyfit(xs,ys,1)[0]
        kind=("ERASES: nothing is remembered" if abs(sl)<1e-6 else
              "PRESERVES exactly: perfect memory" if abs(sl-1)<1e-6 else
              f"AMPLIFIES x{sl:.3f}" if sl>1 else f"PARTIAL retention {sl:.3f}")
        print(f"  bath={nm:3s} g2={g2:7.4f}   slope d<R>(inf)/d<R>(0) = {sl:+.9f}   {kind}")
        if abs(g2)<1e-12 and nm=="cut":
            print("      (point by point:  " +
                  "  ".join(f"{a:+.2f}->{c:+.3f}" for a,c in zip(xs,ys)) + ")")

print()
print("W-33  CONTROL -- W-29c reproduced: an UNBIASED start must give exactly 0 at every setting.")
for links,nm in [(CUT,"cut"),(PERIM,"rim")]:
    for g2 in [0.0,0.01,1.0]:
        r=steady_from(liou(g2,links),(Id/D).astype(complex))
        print(f"    bath={nm:3s} g2={g2:6.3f}  <R>(inf) from maximally mixed = {np.trace(R@r).real:+.12f}")
