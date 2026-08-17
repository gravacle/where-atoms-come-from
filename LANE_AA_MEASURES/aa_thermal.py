# THE ENTROPY WAS REMOVED BY CHOOSING THE GROUND STATE OF A CLOSED SYSTEM, NOT BY THE FIELD.
# Put it back: thermal states at temperature T. Then ask the question a record actually poses --
# is there an ENTROPY DEFICIT? Does anything stay ordered while the rest disorders?
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
def Sred(rho,st,keep,N,L):
    grp=defaultdict(list)
    for j,s in enumerate(st): grp[tuple(s[l] for l in keep)].append(j)
    keys=list(grp); rest=[l for l in range(L) if l not in keep]
    rm=[{tuple(st[j][l] for l in rest):j for j in grp[k]} for k in keys]
    m=len(keys); R=np.zeros((m,m),dtype=complex)
    for a in range(m):
        for b in range(m):
            R[a,b]=sum(rho[ja,rm[b][rk]] for rk,ja in rm[a].items() if rk in rm[b])
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
ev,U=np.linalg.eigh(H)
RIM=[n+k for k in range(n)]; SPK=list(range(n))
Smax=lambda k: k    # log_N of N^k, the maximum for k links
print(f"  wheel n={n}, Z_{N}, {L} links. Entropies in log base {N}. Max for k links = k.")
print(f"  {'T':>8}{'S(global)':>12}{'S(rim)':>10}{'S(spokes)':>12}{'I(rim:spokes)':>15}"
      f"{'rim deficit':>13}{'spoke deficit':>15}")
for Tt in (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 50.0):
    if Tt==0.0:
        p=np.zeros(D); p[0]=1.0
    else:
        b=np.exp(-(ev-ev[0])/Tt); p=b/b.sum()
    rho=(U*p)@U.conj().T
    wv=p[p>1e-12]; Sg=float(-np.sum(wv*np.log(wv))/np.log(N))
    Sr=Sred(rho,st,RIM,N,L); Ss=Sred(rho,st,SPK,N,L)
    print(f"  {Tt:>8.2f}{Sg:>12.6f}{Sr:>10.6f}{Ss:>12.6f}{Sr+Ss-Sg:>15.6f}"
          f"{Smax(n)-Sr:>13.6f}{Smax(n)-Ss:>15.6f}")
print()
print("  DEFICIT = how far a region is from maximally disordered. A record lives in a deficit:")
print("  order held while the surroundings disorder. If rim and spoke deficits fall together the")
print("  carrier has no preferred place to keep one; if they separate, it does.")
