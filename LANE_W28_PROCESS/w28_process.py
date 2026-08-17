# W-28 — THE PROCESS, ASSEMBLED.
# Carrier: 3x3 patch, perimeter is a cycle AND a separator (W-27). Z_2.
# POINTER: the record must (a) commute with H so the dynamics does not scramble it, and (b) be what
# the bath monitors so the environment einselects rather than erases it. W-27b failed both: it
# dephased in Z while H's plaquette terms rotate out of the electric basis.
# In the MAGNETIC phase the plaquette/loop fluxes commute with H. Couple the bath to THOSE.
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
P=[loop(m) for m in PLQ]; Wper=loop([(k,+1) for k in PERIM])
Zk=[np.diag([(-1.0)**s[k] for s in st]).astype(complex) for k in range(L)]

def run(g2,bath,gamma,T,steps=3000):
    H=-sum(P)-g2*sum(Zk)
    Ls=P if bath=="magnetic" else [Zk[k] for k in CUT]
    # pointer check
    pc=max(np.linalg.norm(H@Wper-Wper@H),0)
    def d(rho):
        r=-1j*(H@rho-rho@H)
        for Lo in Ls: r+=gamma*(Lo@rho@Lo.conj().T-rho)
        return r
    # two writes: different perimeter flux sectors
    ev,U=np.linalg.eigh(Wper.real)
    plus=U[:,ev>0]; minus=U[:,ev<0]
    ra=plus@plus.conj().T/plus.shape[1]; rb=minus@minus.conj().T/minus.shape[1]
    dt=T/steps
    for _ in range(steps):
        for r_ in ("a","b"): pass
        k1=d(ra); ra=ra+dt*d(ra+0.5*dt*k1); ra=(ra+ra.conj().T)/2; ra/=np.trace(ra).real
        k1=d(rb); rb=rb+dt*d(rb+0.5*dt*k1); rb=(rb+rb.conj().T)/2; rb/=np.trace(rb).real
    S=lambda r:(lambda w:-float(np.sum(w*np.log2(w))))(np.array([x for x in np.linalg.eigvalsh(r) if x>1e-12]))
    td=0.5*np.sum(np.abs(np.linalg.eigvalsh(ra-rb)))
    return pc,S(ra),td,float(np.real(np.trace(ra@Wper))),float(np.real(np.trace(rb@Wper)))

print(f"  patch 3x3, Z_2, physical dim {D}. RECORD = the perimeter flux (+1 vs -1 sector).")
print(f"  {'phase':>10}{'bath':>11}{'gamma':>7}{'||[H,W]||':>11}{'S(global)':>11}{'TD':>8}"
      f"{'<W>_A':>9}{'<W>_B':>9}")
for g2,tag in ((0.0,"magnetic"),(1.0,"mixed")):
    for bath in ("magnetic","electric"):
        for gamma in (0.5,):
            pc,Sg,td,wa,wb=run(g2,bath,gamma,10.0)
            print(f"  {tag:>10}{bath:>11}{gamma:>7.2f}{pc:>11.2e}{Sg:>11.4f}{td:>8.4f}{wa:>9.4f}{wb:>9.4f}")
print()
print("  ||[H,W]|| = 0 means the record commutes with the dynamics -- a genuine pointer.")
print("  TD near 1 with S(global) > 0 means the distinction SURVIVED the entropy being paid.")
print("  <W>_A and <W>_B are the written values; they should stay near +1 and -1 if the record holds.")
