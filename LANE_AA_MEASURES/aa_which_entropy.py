# WHICH ENTROPY IS THERE, AND WHICH DOES A RECORD NEED?
# Four candidates, computed on the wheel. Z_3, ground state, theta=1.0, g^2=1.0.
import numpy as np, itertools
from collections import defaultdict
def wheel(n): return [(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)], n+1
def physical(E,V,N,L):
    keep=[]
    for s in itertools.product(range(N),repeat=L):
        ok=True
        for v in range(V):
            t=sum(s[i] for i,(a,b) in enumerate(E) if a==v)-sum(s[i] for i,(a,b) in enumerate(E) if b==v)
            if t%N: ok=False; break
        if ok: keep.append(s)
    return keep
def loop_op(st,idx,moves,N):
    D=len(st); M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        M[idx[tuple(t)],j]=1.0
    return M
def Sreg(psi,st,keep,N,L):
    grp=defaultdict(list)
    for j,s in enumerate(st): grp[tuple(s[l] for l in keep)].append(j)
    keys=list(grp); rest=[l for l in range(L) if l not in keep]
    rm=[{tuple(st[j][l] for l in rest):j for j in grp[k]} for k in keys]
    m=len(keys); R=np.zeros((m,m),dtype=complex)
    for a in range(m):
        for b in range(m):
            R[a,b]=sum(psi[ja]*np.conj(rm[b][rk]) if False else psi[ja]*np.conj(psi[rm[b][rk]])
                       for rk,ja in rm[a].items() if rk in rm[b])
    w=np.linalg.eigvalsh(R); w=w[w>1e-12]
    return float(-np.sum(w*np.log(w))/np.log(N))

N=3; n=5
E,V=wheel(n); L=len(E)
st=physical(E,V,N,L); idx={s:j for j,s in enumerate(st)}; D=len(st)
T=[loop_op(st,idx,[(k,+1),(n+k,+1),((k+1)%n,-1)],N) for k in range(n)]
w=np.exp(2j*np.pi/N)
Ee=np.diag([sum(w**s[l] for l in range(L)) for s in st])
Hm=sum(np.exp(1j)*x+np.exp(-1j)*x.conj().T for x in T)
H=-(Hm)-1.0*(Ee+Ee.conj().T); H=(H+H.conj().T)/2
ev,U=np.linalg.eigh(H); psi=U[:,0]

print(f"  wheel n={n}, Z_{N}: V={V} vertices, L={L} links.  All entropies in log base {N}.\n")
print("  1. THERMODYNAMIC / GLOBAL ENTROPY OF THE STATE")
rho_pure=np.outer(psi,psi.conj())
wv=np.linalg.eigvalsh(rho_pure); wv=wv[wv>1e-12]
print(f"     S(whole system) = {float(-np.sum(wv*np.log(wv))/np.log(N)):.12f}")
print("     ZERO, and it stays zero under unitary evolution. The system is closed and pure.")
print("     There is no bath, no coarse-graining, no dissipation. NOTHING CAN BE PAID HERE.\n")
print("  2. ENTANGLEMENT ENTROPY OF A REGION  (what I have been measuring)")
for j in (1,2,3):
    print(f"     S(first {j} rim links) = {Sreg(psi,st,[n+k for k in range(j)],N,L):.6f}")
print("     Linear in the link count. This is a bookkeeping of correlations INSIDE a pure state.")
print("     It is not produced, not dissipated, and not conserved-and-transferred. It is a ledger.\n")
print("  3. CONFIGURATIONAL ENTROPY -- NATIVE, imported from nothing")
print(f"     dim(full) = {N}^{L} = {N**L};  dim(physical) = {D} = {N}^{round(np.log(D)/np.log(N))}")
print(f"     S_config = log_{N} dim(physical) = {np.log(D)/np.log(N):.4f} = the CYCLE RANK of the graph")
print("     This is a CAPACITY: how many states a record could occupy. It is combinatorial,")
print("     needs no state, no Hamiltonian and no coupling, and is fixed by the carrier alone.\n")
print("  4. GAUGE ENTROPY -- the information DISCARDED by imposing gauge invariance")
orbit=N**(V-1)
print(f"     gauge orbit size = {N}^(V-1) = {orbit};  S_gauge = {V-1:.0f}.0000")
print(f"     physical {np.log(D)/np.log(N):.0f} + gauge {V-1} = {np.log(D)/np.log(N)+V-1:.0f} = L = {L}. Exactly half each on a wheel.")
print("     AND PER C-1 THIS IS WHERE THE RECORD WOULD LIVE: the coarse-graining that would")
print("     produce entropy is exactly the one that throws the record away.")
