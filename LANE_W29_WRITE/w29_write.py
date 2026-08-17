# W-29 — THE WRITE STEP. W-28 measured PERSISTENCE: the value was put there by hand and stayed.
# FORMATION means the record ACQUIRES a value it did not have, correlated with something that
# happened. Start with the record UNDETERMINED and vary the interior; does the record end up
# carrying the interior's state?
import numpy as np, itertools
from collections import defaultdict
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E); NV=len(V2); N=2; CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]; PERIM=[k for k in range(L) if k not in CUT]
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(NV))]
idx={s:j for j,s in enumerate(st)}; D=len(st)
def hid(i,j): return j*2+i
def vidx(i,j): return 6+j*3+i
PLQ=[[(hid(i,j),+1),(vidx(i+1,j),+1),(hid(i,j+1),-1),(vidx(i,j),-1)] for j in range(2) for i in range(2)]
def loop(moves):
    M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        M[idx[tuple(t)],j]=1.0
    return M
P=[loop(m) for m in PLQ]; W=loop([(k,+1) for k in PERIM])
Zk=[np.diag([(-1.0)**s[k] for s in st]).astype(complex) for k in range(L)]

# THE INTERIOR (what happens) vs THE RECORD (the perimeter flux)
Zin=Zk[CUT[0]]                                    # an interior observable
print(f"  patch 3x3, Z_2, physical dim {D}")
print(f"  INTERIOR observable: Z on cut link {CUT[0]}.   RECORD: perimeter flux W.")
print(f"  || [Z_interior, W] || = {np.linalg.norm(Zin@W-W@Zin):.3e}   "
      f"(nonzero = the interior can drive the record)")

def evolve(rho,H,Ls,gamma,T,steps=3000):
    dt=T/steps
    def d(r):
        x=-1j*(H@r-r@H)
        for Lo in Ls: x+=gamma*(Lo@r@Lo.conj().T-r)
        return x
    for _ in range(steps):
        k1=d(rho); rho=rho+dt*d(rho+0.5*dt*k1); rho=(rho+rho.conj().T)/2; rho/=np.trace(rho).real
    return rho
def proj(Op,val):
    ev,U=np.linalg.eigh(Op.real); sel=U[:,np.abs(ev-val)<1e-9]; return sel@sel.conj().T

print(f"\n  THE WRITE. Record starts UNDETERMINED (maximally mixed over both flux sectors).")
print(f"  Interior prepared in +1 or -1. Evolve, then read the record.")
print(f"  {'g2':>5}{'gamma':>7}{'T':>6}{'<W> | interior +1':>20}{'<W> | interior -1':>20}{'SEPARATION':>13}")
for g2,gamma,T in ((0.0,0.5,10.0),(1.0,0.5,10.0),(1.0,0.0,10.0),(0.5,0.5,20.0),(2.0,0.5,20.0)):
    H=-sum(P)-g2*sum(Zk)
    out=[]
    for val in (+1.0,-1.0):
        Pin=proj(Zin,val)
        rho=Pin/np.trace(Pin).real                 # interior definite, record undetermined within it
        rho=evolve(rho,H,P,gamma,T)
        out.append(float(np.real(np.trace(rho@W))))
    print(f"  {g2:>5.1f}{gamma:>7.2f}{T:>6.1f}{out[0]:>20.6f}{out[1]:>20.6f}{abs(out[0]-out[1]):>13.6f}")
print()
print("  SEPARATION > 0 means the record ENDED UP carrying the interior's state -- it was written.")
print("  SEPARATION = 0 means nothing was written: the record is blind to what happened inside.")
