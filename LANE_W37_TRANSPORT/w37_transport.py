"""W-37.  CAN A LOCAL PROBE READ THE NONLOCAL RECORD?  AHARONOV-BOHM TRANSPORT.

W-36: the record R (a product over 8 rim links) is NONLOCAL, and a local environment cannot copy it,
so it is not objective. But gauge theory has a standard answer -- transport a charge around the loop.
Every interaction is LOCAL at every instant; the nonlocality lives in the PATH and is paid for in
TIME proportional to the perimeter.

CONSTRUCTION. Explicit matter, no shortcuts:
  a PROBE CHARGE that hops between the 8 rim vertices, plus a STATIC ANTI-CHARGE at the centre
  (without it the total Z_2 charge is odd and NO physical state exists at all -- a genuine
  constraint, checked in code rather than assumed).
  GAUSS becomes  div(s)_v = [v is the probe's site] + [v is the centre]  (mod 2),
  so the physical sector MOVES with the probe -- the constraint is matter-dependent.
  HOP  H = -tau * sum_k ( |k+1><k| (x) U_link(k) + h.c. ),  U = the parallel transporter on that
  link. This is the ONLY interaction, it touches ONE link at a time, and it never mentions R.
  The probe spreads both ways around the ring and interferes: that interference is the loop.

MEASURE. Holevo information the PROBE POSITION holds about R, in bits, ceiling 1.
CONTROLS: tau=0 must give 0; a CUT RING (one hop removed, so no closed path) must give 0 -- that is
the decisive control, since it keeps every local interaction and removes only the topology.
"""
import itertools, numpy as np

V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E); NV=9; CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
PERIM=[k for k in range(L) if k not in CUT]
RIMV=[vid[(0,0)],vid[(1,0)],vid[(2,0)],vid[(2,1)],vid[(2,2)],vid[(1,2)],vid[(0,2)],vid[(0,1)]]
def link_between(u,v):
    for k,(a,b) in enumerate(E):
        if (a,b)==(u,v) or (b,a)==(u,v): return k
    raise ValueError((u,v))
RIML=[link_between(RIMV[k],RIMV[(k+1)%8]) for k in range(8)]
assert sorted(RIML)==sorted(PERIM), (RIML,PERIM)
print("W-37  probe hops the 8 rim vertices; static anti-charge at the centre.")
print(f"      rim vertices {RIMV}\n      rim links    {RIML}  (== PERIM, so the ring IS the record's support)")

def div(s,v):
    return (sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2

STATES=[]
for p in range(8):
    for s in itertools.product(range(2),repeat=L):
        q=lambda v: (1 if v==RIMV[p] else 0)^(1 if v==CENTER else 0)
        if all(div(s,v)==q(v) for v in range(NV)): STATES.append((p,s))
IDX={x:i for i,x in enumerate(STATES)}; D=len(STATES)
print(f"      physical dim with matter = {D}   ({D//8} gauge states per probe position)")

# CHECK the constraint is real: drop the anti-charge and the sector must be EMPTY.
empty=[s for s in itertools.product(range(2),repeat=L)
       if all(div(s,v)==(1 if v==RIMV[0] else 0) for v in range(NV))]
print(f"      CONSTRAINT CHECK -- probe with NO compensating anti-charge: {len(empty)} physical states"
      f"  ({'correct: odd total charge is forbidden' if not empty else 'UNEXPECTED'})")

def shift(s,links):
    t=list(s)
    for k in links: t[k]^=1
    return tuple(t)
def Rop():
    M=np.zeros((D,D),complex)
    for j,(p,s) in enumerate(STATES):
        M[IDX[(p,shift(s,PERIM))],j]=1.0
    return M
def Hhop(tau,skip=None):
    M=np.zeros((D,D),complex)
    for k in range(8):
        if skip is not None and k==skip: continue
        kn=(k+1)%8
        for j,(p,s) in enumerate(STATES):
            if p!=k: continue
            t=shift(s,[RIML[k]])
            i=IDX.get((kn,t))
            if i is not None: M[i,j]-=tau
    return M+M.conj().T
R=Rop()
print(f"      R unitary defect {np.linalg.norm(R@R.conj().T-np.eye(D)):.1e}, "
      f"eigenvalues {len(np.unique(np.round(np.linalg.eigvals(R),6)))}, "
      f"||[R,H_hop]|| = {np.linalg.norm(R@Hhop(1.0)-Hhop(1.0)@R):.1e}")
Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2

def expm(A):
    n=max(0,int(np.ceil(np.log2(np.linalg.norm(A,np.inf))))+1) if np.linalg.norm(A,np.inf)>0 else 0
    B=A/(2.0**n); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for k in range(1,40):
        T=T@B/k; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(n): X=X@X
    return X
def vn(r):
    ev=np.linalg.eigvalsh((r+r.conj().T)/2); ev=ev[ev>1e-12]
    return float(-(ev*np.log2(ev)).sum())
def init(seed=5):
    g=np.random.default_rng(seed)
    mask=np.array([1.0 if p==0 else 0.0 for p,_ in STATES])
    w=(g.normal(size=D)+1j*g.normal(size=D))*mask
    a=Pp@w; b=Pm@w; a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    v=(a+b); return v/np.linalg.norm(v)
def probe_rdm(psi):
    """trace out the gauge field, keep the probe's position (8 states)"""
    M=np.zeros((8,8),complex)
    for i,(p,s) in enumerate(STATES):
        for j,(q,t) in enumerate(STATES):
            if s==t: M[p,q]+=psi[i]*np.conj(psi[j])
    return M
def info(psi):
    br=[]
    for Pr in (Pp,Pm):
        v=Pr@psi; pr=float(np.vdot(v,v).real)
        br.append((pr, probe_rdm(v/np.sqrt(pr)) if pr>1e-14 else None))
    avg=sum(p*m for p,m in br if m is not None)
    return vn(avg)-sum(p*vn(m) for p,m in br if m is not None)

psi0=init()
print(f"\n      initial state: probe at site 0, <R> = {np.vdot(psi0,R@psi0).real:+.1e} (unbiased by construction)")
print("\n  I(R : probe position) in bits as the probe goes around.  tau = 1.0")
print(f"  {'T':>7s} {'I(R:probe)':>12s}   {'P(probe back at 0)':>19s}")
print("  "+"-"*46)
Hh=Hhop(1.0)
for T in [0.0,0.5,1.0,2.0,3.0,4.0,6.0,8.0,12.0,20.0]:
    psi=expm(-1j*Hh*T)@psi0
    back=sum(abs(psi[i])**2 for i,(p,_) in enumerate(STATES) if p==0)
    print(f"  {T:7.2f} {info(psi):12.6f}   {back:19.4f}")

print("\n  CONTROL A -- tau = 0 (no transport at all). Must be 0 at every time.")
for T in [4.0,12.0]:
    print(f"    T={T:5.1f}   I = {info(expm(-1j*Hhop(0.0)*T)@psi0):.3e}")
print("\n  CONTROL B -- THE DECISIVE ONE. Cut the ring (remove ONE hop), keeping every local")
print("  interaction and removing only the closed path. If the information is topological it dies.")
Hc=Hhop(1.0,skip=4)
for T in [4.0,8.0,12.0,20.0]:
    psi=expm(-1j*Hc*T)@psi0
    print(f"    T={T:5.1f}   I(cut ring) = {info(psi):.6f}")
