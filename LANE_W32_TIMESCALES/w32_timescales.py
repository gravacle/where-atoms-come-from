"""W-32.  THE APPROXIMATE-CONSERVATION ESCAPE.  Does a WINDOW exist, and does the coupling set it?

W-30 concluded 'no value of the coupling escapes the obstruction'. That reading measured survival
at a FIXED ABSOLUTE TIME while the WRITE timescale was itself changing with g^2. Records in nature
are never exactly conserved -- magnetic domains flip, crystals anneal, DNA mutates. What makes them
records is tau_decay >> tau_observation. So the question is not 'does <R> survive to t=10', it is:

        HOW MANY TIMES CAN THE RECORD BE WRITTEN BEFORE IT DECAYS?

That ratio is dimensionless. If it diverges as the coupling goes to zero, a window exists and the
coupling sets its width -- which would be a SECOND role for alpha, distinct from W-30's.

MEASURED EXACTLY, NO TIME-SERIES FITTING.
  Gamma = the slowest nonzero decay rate of the LINDBLADIAN SPECTRUM restricted to the sector that
          can change <R>. Obtained as an eigenvalue of the D^2 x D^2 generator. No fitting window,
          no integrator, no choice of T.
  omega = the coherent rate at which <R> is driven, ||[H,R]||/D.
  RATIO = omega / Gamma  =  coherent write-rotations per record lifetime.

CONTROLS.  (a) at g2=0 omega must be exactly 0.  (b) Gamma must vanish when the bath is off.
(c) the power laws are read off successive log-log slopes, NOT assumed.
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

def analyse(N, bathname, bathlinks, g2list):
    st,idx=build(V2,E,N); D=len(st)
    R=Move(st,idx,compose(P),N)
    Lmag=[Move(st,idx,p,N) for p in P]
    MAG=sum(L+L.conj().T for L in Lmag)
    ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(len(E)))
    Ls=[Zop(st,[k],N) for k in bathlinks]
    Id=np.eye(D,dtype=complex)
    cLR=max(np.linalg.norm(L@R-R@L) for L in Ls)
    print(f"\n  N={N}  dim={D}  bath = electric on {bathname} ({len(Ls)} ops)   "
          f"max||[L,R]|| = {cLR:.3e}")
    print(f"  {'g^2':>8s} {'omega=||[H,R]||/D':>18s} {'Gamma(slowest)':>15s} "
          f"{'RATIO omega/Gamma':>18s} {'d ln(Gam)/d ln(g2)':>19s}")
    print("  "+"-"*84)
    prev=None; rows=[]
    for g2 in g2list:
        H=-MAG-g2*ELEC
        M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
        for L in Ls: M+=0.5*(np.kron(L.conj(),L)-np.kron(Id,Id))
        ev=np.linalg.eigvals(M)
        re=-ev.real
        dec=re[re>1e-10]
        Gam=dec.min() if dec.size else 0.0
        om=np.linalg.norm(H@R-R@H)/D
        slope=""
        if prev and prev[1]>0 and Gam>0 and prev[0]>0 and g2>0:
            slope=f"{(np.log(Gam)-np.log(prev[1]))/(np.log(g2)-np.log(prev[0])):19.3f}"
        rows.append((g2,om,Gam,om/Gam if Gam>0 else float('inf')))
        print(f"  {g2:8.4f} {om:18.6f} {Gam:15.6e} "
              f"{(om/Gam if Gam>0 else float('inf')):18.4f} {slope:>19s}")
        if g2>0: prev=(g2,Gam)
    return rows

print("W-32   3x3 patch. RATIO = coherent write-rotations per record lifetime (dimensionless).")
g2s=[0.0,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0]

print("\n=== BATH ON THE CUT LINKS (disjoint from the rim -- commutes with R) ===")
rows_cut=analyse(2,"the CUT",CUT,g2s)

print("\n=== CONTROL: BATH ON THE RIM ITSELF (the bath now looks directly at the record) ===")
rows_rim=analyse(2,"the RIM",PERIM,g2s)

print("\n=== CONTROL: bath OFF entirely. Gamma must be 0 at every g^2. ===")
st,idx=build(V2,E,2); D=len(st)
R=Move(st,idx,compose(P),2); Lmag=[Move(st,idx,p,2) for p in P]
MAG=sum(L+L.conj().T for L in Lmag)
ELEC=sum(Zop(st,[k],2)+Zop(st,[k],2).conj().T for k in range(len(E)))
Id=np.eye(D,dtype=complex)
for g2 in [0.0,0.05,1.0]:
    H=-MAG-g2*ELEC
    M=-1j*(np.kron(Id,H)-np.kron(H.T,Id))
    ev=np.linalg.eigvals(M); re=-ev.real; dec=re[re>1e-10]
    print(f"    g2={g2:6.3f}  gamma=0   #decaying modes = {dec.size}   "
          f"max|Re(eigenvalue)| = {np.abs(ev.real).max():.3e}")

print("\n=== READING THE SCALING ===")
print("  If Gamma ~ g^4 and omega ~ g^2 then RATIO ~ 1/g^2 and DIVERGES as the coupling -> 0:")
print("  a window exists, and its width is set by the coupling.")
print("  If RATIO is flat or falls as g^2 -> 0, no window: the W-30 obstruction stands.")
for nm,rows in [("cut bath",rows_cut),("rim bath",rows_rim)]:
    fin=[r for r in rows if r[0]>0 and np.isfinite(r[3])]
    if len(fin)>=2:
        a,b=fin[0],fin[-1]
        s=(np.log(b[3])-np.log(a[3]))/(np.log(b[0])-np.log(a[0]))
        print(f"    {nm}: RATIO at g2={a[0]} is {a[3]:.2f}; at g2={b[0]} is {b[3]:.2f};  "
              f"d ln(RATIO)/d ln(g2) = {s:.3f}")
