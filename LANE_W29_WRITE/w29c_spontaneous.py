# W-29c — DOES THE DYNAMICS FORCE THE RECORD ON ITS OWN?
# W-29/29b treated writing as an external act and hit an obstruction: a record that commutes with H
# survives but cannot be written; one that does not commute can be written but does not survive.
# THE PRINCIPAL'S POINT: if the dynamics of the parts FORCES the record, that is the process.
# So: start MAXIMALLY UNDETERMINED and let it run. Does the flux sector become definite by itself?
import numpy as np, itertools
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
def evolve(rho,H,Ls,gamma,T,steps=4000):
    dt=T/steps
    def d(r):
        x=-1j*(H@r-r@H)
        for Lo in Ls: x+=gamma*(Lo@r@Lo.conj().T-r)
        return x
    for _ in range(steps):
        k1=d(rho); rho=rho+dt*d(rho+0.5*dt*k1); rho=(rho+rho.conj().T)/2; rho/=np.trace(rho).real
    return rho
S=lambda r:(lambda w:-float(np.sum(w*np.log2(w))))(np.array([x for x in np.linalg.eigvalsh(r) if x>1e-12]))
rho0=np.eye(D,dtype=complex)/D                     # MAXIMALLY UNDETERMINED. Nothing written.
print(f"  patch 3x3, Z_2, dim {D}. START: maximally mixed. <W> at t=0 = {float(np.real(np.trace(rho0@W))):.6f}")
print(f"  Does the record become definite ON ITS OWN?  |<W>| -> 1 would mean the dynamics wrote it.")
print(f"\n  {'bath':>12}{'g2':>5}{'gamma':>7}{'T':>6}{'<W>':>12}{'|<W>|':>10}{'S(final)':>11}")
for tag,Ls in (("plaquette",P),("electric",[Zk[k] for k in CUT]),
               ("proj +1",[ (np.eye(D)+W.real)/2 ]),("none",[])):
    for g2,gamma,T in ((1.0,0.5,20.0),(0.0,0.5,20.0)):
        if tag=="none" and gamma>0: gamma=0.0
        r=evolve(rho0.copy(),-sum(P)-g2*sum(Zk),Ls,gamma,T)
        wv=float(np.real(np.trace(r@W)))
        print(f"  {tag:>12}{g2:>5.1f}{gamma:>7.2f}{T:>6.1f}{wv:>12.6f}{abs(wv):>10.6f}{S(r):>11.4f}")
print()
print("  A bath whose Lindblad operator is a PROJECTOR onto a flux sector is the one that could")
print("  select. If even that leaves |<W>| at 0, no dynamics of this form writes a record here.")
