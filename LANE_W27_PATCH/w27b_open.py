# W-27b — OPEN THE SYSTEM. A closed pure system has S = 0 forever and can pay for nothing.
# Couple the ELECTRIC BOUNDARY (the cut at the interior vertex) to a bath that dephases it in the
# electric basis -- the environment "reads" the boundary. Lindblad:
#   drho/dt = -i[H,rho] + gamma * sum_k ( Z_k rho Z_k^dag - rho )
# THE PROCESS TEST: write two different interior configurations, evolve, and ask whether a
# measurement ON THE BOUNDARY still tells them apart after the entropy has been paid.
import numpy as np, itertools
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E); NV=len(V2); N=2; CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
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
P=[loop(m) for m in PLQ]
Zk=[np.diag([(-1.0)**s[k] for s in st]).astype(complex) for k in range(L)]
H=-sum(P)-1.0*sum(Zk)
def lind(rho,gamma):
    d=-1j*(H@rho-rho@H)
    for k in CUT: d+=gamma*(Zk[k]@rho@Zk[k].conj().T-rho)
    return d
def evolve(rho,gamma,T,steps=4000):
    dt=T/steps
    for _ in range(steps):
        k1=lind(rho,gamma); k2=lind(rho+0.5*dt*k1,gamma)
        rho=rho+dt*k2
        rho=(rho+rho.conj().T)/2; rho/=np.trace(rho).real
    return rho
def S(rho):
    w=np.linalg.eigvalsh(rho); w=w[w>1e-12]; return float(-np.sum(w*np.log2(w)))
def red(rho,keep):
    from collections import defaultdict
    g=defaultdict(list)
    for j,s in enumerate(st): g[tuple(s[l] for l in keep)].append(j)
    ks=list(g); rest=[l for l in range(L) if l not in keep]
    rm=[{tuple(st[j][l] for l in rest):j for j in g[k]} for k in ks]
    m=len(ks); R=np.zeros((m,m),dtype=complex)
    for a in range(m):
        for b in range(m):
            R[a,b]=sum(rho[ja,rm[b][rk]] for rk,ja in rm[a].items() if rk in rm[b])
    return R
td=lambda a,b: 0.5*np.sum(np.abs(np.linalg.eigvalsh(a-b)))
# two interior configurations differing in the CUT distribution at fixed total
tot=lambda s: sum(s[l] for l in CUT)%N
pr=[(i,j) for i in range(D) for j in range(i+1,D)
    if tot(st[i])==tot(st[j]) and tuple(st[i][l] for l in CUT)!=tuple(st[j][l] for l in CUT)]
i0,j0=pr[0]
print(f"  patch 3x3, Z_2, physical dim {D}.  electric boundary = links {CUT}")
print(f"  A cut-config {tuple(st[i0][l] for l in CUT)} vs B {tuple(st[j0][l] for l in CUT)}  (same total)")
print(f"\n  {'gamma':>7}{'T':>6}{'S(global)':>12}{'S(boundary)':>13}{'TD on boundary':>16}")
for gamma in (0.0,0.2,1.0):
    for T in (2.0,10.0):
        ra=np.zeros((D,D),dtype=complex); ra[i0,i0]=1
        rb=np.zeros((D,D),dtype=complex); rb[j0,j0]=1
        ra=evolve(ra,gamma,T); rb=evolve(rb,gamma,T)
        print(f"  {gamma:>7.2f}{T:>6.1f}{S(ra):>12.6f}{S(red(ra,CUT)):>13.6f}{td(red(ra,CUT),red(rb,CUT)):>16.6f}")
print("\n  gamma=0 is the closed control: S(global) must stay 0 and nothing is paid.")
print("  gamma>0 pays entropy. If TD stays high while S(global) rises, the boundary keeps the")
print("  distinction THROUGH the dissipation -- which is what a written record does.")
