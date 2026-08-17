# W-26 — THE BOUNDARY TOTAL IS A MIRROR. IS THE BOUNDARY DISTRIBUTION A MEMORY?
# Gauss fixes the SUM of the n boundary links (= the enclosed charge) at every instant, so the total
# can never disagree with the present interior and can never retain a past one. But Gauss says
# NOTHING about how flux is DISTRIBUTED among the n links: n-1 boundary degrees of freedom are
# unconstrained by the interior. Every boundary measurement in this program has been of the total.
# TEST: two states with the SAME boundary total and DIFFERENT distributions. Evolve. Does a boundary
# measurement still tell them apart later?
import numpy as np, itertools
from collections import defaultdict
def wheel(n): return [(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)], n+1
N=3; n=4
E,V=wheel(n); L=len(E); SPK=list(range(n)); RIM=list(range(n,2*n))
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[i] for i,(a,b) in enumerate(E) if a==v)
           -sum(s[i] for i,(a,b) in enumerate(E) if b==v))%N==0 for v in range(V))]
idx={s:j for j,s in enumerate(st)}; D=len(st)
def loop_op(moves):
    M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        M[idx[tuple(t)],j]=1.0
    return M
T=[loop_op([(k,+1),(n+k,+1),((k+1)%n,-1)]) for k in range(n)]
w=np.exp(2j*np.pi/N)
Ee=np.diag([sum(w**s[l] for l in range(L)) for s in st])
H=-(sum(np.exp(1j)*x+np.exp(-1j)*x.conj().T for x in T))-1.0*(Ee+Ee.conj().T)
H=(H+H.conj().T)/2; ev,U=np.linalg.eigh(H)

def red(psi,keep):
    grp=defaultdict(list)
    for j,s in enumerate(st): grp[tuple(s[l] for l in keep)].append(j)
    keys=list(grp); rest=[l for l in range(L) if l not in keep]
    rm=[{tuple(st[j][l] for l in rest):j for j in grp[k]} for k in keys]
    m=len(keys); R=np.zeros((m,m),dtype=complex)
    for a in range(m):
        for b in range(m):
            R[a,b]=sum(psi[ja]*np.conj(psi[rm[b][rk]]) for rk,ja in rm[a].items() if rk in rm[b])
    return R
td=lambda a,b: 0.5*np.sum(np.abs(np.linalg.eigvalsh(a-b)))
tot=lambda s: sum(s[l] for l in SPK)%N

# two basis states: SAME boundary total, DIFFERENT boundary distribution
pairs=[(i,j) for i in range(D) for j in range(i+1,D)
       if tot(st[i])==tot(st[j]) and tuple(st[i][l] for l in SPK)!=tuple(st[j][l] for l in SPK)]
i0,j0=pairs[0]
print(f"  wheel n={n}, Z_{N}, physical dim {D}")
print(f"  state A spokes {tuple(st[i0][l] for l in SPK)}  total {tot(st[i0])}")
print(f"  state B spokes {tuple(st[j0][l] for l in SPK)}  total {tot(st[j0])}   <- SAME total, different distribution")
a=np.zeros(D,dtype=complex); a[i0]=1; b=np.zeros(D,dtype=complex); b[j0]=1
Ua,Ub=U.conj().T@a,U.conj().T@b
print(f"\n  {'t':>7}{'TD on the BOUNDARY (spokes)':>30}{'TD on the RIM (control)':>26}")
for t in (0.0,0.5,1.0,2.0,5.0,20.0,100.0):
    pa=U@(np.exp(-1j*ev*t)*Ua); pb=U@(np.exp(-1j*ev*t)*Ub)
    print(f"  {t:>7.1f}{td(red(pa,SPK),red(pb,SPK)):>30.6f}{td(red(pa,RIM),red(pb,RIM)):>26.6f}")
print("\n  TD = 1 means perfectly distinguishable, 0 means the difference is gone.")
print("  If the boundary TD stays high while the total is pinned by Gauss, the DISTRIBUTION is")
print("  carrying something the total structurally cannot.")
