# IS AN ENTROPY AREA LAW FORCED ON A LATTICE GAUGE THEORY?
# The Bekenstein-shaped expectation is S(region) ~ |boundary|. On a lattice the region's entropy
# carries a term from the BOUNDARY FLUX DISTRIBUTION -- counting the flux values available on the
# boundary links. If S is determined by the boundary link COUNT alone, an area law is bookkeeping.
# TEST: vary the region at FIXED boundary size, and vary the boundary size at fixed region size.
# If S tracks the boundary count and nothing else, it is forced.
import numpy as np, itertools
def wheel(n):
    return [(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)], n+1
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
        k=idx.get(tuple(t))
        if k is None: return None
        M[k,j]=1.0
    return M
def S_of(psi,st,keep,N,L):
    """entropy of the reduced state on a set of LINKS, in the physical basis"""
    from collections import defaultdict
    grp=defaultdict(list)
    for j,s in enumerate(st):
        grp[(tuple(s[l] for l in keep))].append(j)
    keys=list(grp); m=len(keys)
    R=np.zeros((m,m),dtype=complex)
    rest=[l for l in range(L) if l not in keep]
    restmap=defaultdict(dict)
    for a,k in enumerate(keys):
        for j in grp[k]: restmap[a][tuple(st[j][l] for l in rest)]=j
    for a in range(m):
        for b in range(m):
            tot=0
            for rk,ja in restmap[a].items():
                jb=restmap[b].get(rk)
                if jb is not None: tot+=psi[ja]*np.conj(psi[jb])
            R[a,b]=tot
    w=np.linalg.eigvalsh(R); w=w[w>1e-12]
    return float(-np.sum(w*np.log(w))/np.log(N))

N=3
print(f"  Z_{N}. Ground state at theta=1.0, g^2=1.0. Entropy in log base {N}.")
print(f"  {'wheel n':>8}{'region (links)':>28}{'|region|':>10}{'|boundary|':>12}{'S(region)':>12}")
for n in (4,5,6):
    E,V=wheel(n); L=len(E)
    st=physical(E,V,N,L); idx={s:j for j,s in enumerate(st)}
    T=[loop_op(st,idx,[(k,+1),(n+k,+1),((k+1)%n,-1)],N) for k in range(n)]
    w=np.exp(2j*np.pi/N)
    Ee=np.diag([sum(w**s[l] for l in range(L)) for s in st])
    Hm=sum(np.exp(1j*1.0)*x+np.exp(-1j*1.0)*x.conj().T for x in T)
    H=-(Hm)-1.0*(Ee+Ee.conj().T); H=(H+H.conj().T)/2
    psi=np.linalg.eigh(H)[1][:,0]
    # regions: growing arcs of the rim. boundary of an arc of j rim links = the 2 spokes at its ends
    for j in range(1,n):
        region=[n+k for k in range(j)]
        bdry=2                                     # an arc always has two endpoints
        print(f"  {n:>8}{str(region):>28}{j:>10}{bdry:>12}{S_of(psi,st,region,N,L):>12.6f}")
    print()
print("  If S depends only on |boundary| it would be CONSTANT down each block (boundary is always 2).")
print("  If S grows with |region| the entropy is NOT determined by the boundary count, and an")
print("  area law here would be a real statement rather than bookkeeping.")
