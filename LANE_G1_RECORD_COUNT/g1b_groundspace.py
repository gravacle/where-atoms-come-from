"""G1b.  HOW MANY RECORDS?  -- counted ON THE GROUND SPACE, where a record actually lives.

g1 checked clause (v) on the FULL Hilbert space and every candidate passed trivially: a short string
maps OUT of the ground space, so it cannot flip anything there, and the "protection" test was vacuous.
It also counted operators that commute with H but act as the IDENTITY on the ground space. Both are
fixed by doing the whole thing inside the code space, which is where the record is.

THE COUNT, done properly:
    candidates    Z(S) with S in Z^1 (even overlap with every plaquette) -- these commute with H
    non-trivial   S not in B^1 (not a combination of vertex stars) -- otherwise it is the identity
                  on the ground space
    independent   the GF(2) rank of what survives
  PREDICTION from Theorem A: dim H^1 = 2g = 2, giving a record space of 2^(2g) = 4.

AND THE QUESTION THIS SETTLES: T3's `capacity = area - 1` counts plaquette-subset boundary loops.
Do any of those survive as records? A plaquette boundary is a COBOUNDARY-like object; if it acts as
the identity on the ground space it was never a record, and `area - 1` was never a record count.
"""
import itertools, numpy as np
def torus(nx,ny):
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
        for j in range(ny) for i in range(nx)]
    return nx*ny,E,PL
def analyse(nx,ny,label):
    NV,E,PL=torus(nx,ny); L=len(E)
    st=[s for s in itertools.product(range(2),repeat=L)
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(NV))]
    idx={s:i for i,s in enumerate(st)}; D=len(st)
    def Move(S):
        M=np.zeros((D,D),complex)
        for j,s in enumerate(st):
            t=list(s)
            for k in S: t[k]^=1
            t=tuple(t)
            if t in idx: M[idx[t],j]=1.0
        return M
    def Zl(S): return np.diag([(-1.0)**(sum(s[k] for k in S)%2) for s in st]).astype(complex)
    H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
    ev,U=np.linalg.eigh(H); tol=1e-8*max(1.0,abs(ev).max())
    g=int(sum(1 for e in ev if abs(e-ev[0])<tol)); G=U[:,:g]
    def vv(S):
        v=0
        for k in S: v|=(1<<k)
        return v
    def bits(v): return [k for k in range(L) if v>>k&1]
    def rank_gf2(vs):
        b=[]; r=0
        for v in vs:
            cur=v
            for x in b:
                p=x.bit_length()-1
                if cur>>p&1: cur^=x
            if cur: b.append(cur); b.sort(reverse=True); r+=1
        return r
    # vertex stars span B^1 -- two link sets differing by one of these act IDENTICALLY on the ground
    # space, so the record count must be taken in the QUOTIENT Z^1 / B^1, not in the link space.
    stars=[vv([k for k in range(L) if E[k][0]==x or E[k][1]==x]) for x in range(NV)]
    def rk(vs):
        b=[]; r=0
        for v in vs:
            cur=v
            for x in b:
                p=x.bit_length()-1
                if cur>>p&1: cur^=x
            if cur: b.append(cur); b.sort(reverse=True); r+=1
        return r
    rstars=rk(stars)
    Zcoc=[v for v in range(1,1<<L) if all(bin(v & vv(P)).count('1')%2==0 for P in PL)]
    surv=[]
    for v in Zcoc:
        b=G.conj().T@Zl(bits(v))@G
        if np.linalg.norm(b-np.trace(b)/g*np.eye(g))>1e-8: surv.append(v)   # non-trivial on ground space
    # T3's objects: plaquette-subset boundary loops
    t3=set()
    for r in range(1,len(PL)+1):
        for T in itertools.combinations(range(len(PL)),r):
            c={}
            for p in T:
                for lk in PL[p]: c[lk]=c.get(lk,0)+1
            t3.add(vv([lk for lk,n in c.items() if n%2]))
    t3.discard(0)
    t3surv=[v for v in t3 if v in surv]
    print(f"  {label}: dim {D}, ground degeneracy {g}")
    print(f"    Z-operators commuting with H            : {len(Zcoc)}")
    indep = rk(list(stars)+surv) - rstars          # rank MODULO the vertex stars = dim H^1
    print(f"    of those, NON-TRIVIAL on the ground space: {len(surv)}")
    print(f"    INDEPENDENT modulo vertex stars (= dim H^1): {indep}   [raw link-space rank was {rank_gf2(surv)}]")
    print(f"    predicted 2g = {int(np.log2(g))}   -> record space 2^(2g) = {g}")
    t3indep = rk(list(stars)+t3surv) - rstars if t3surv else 0
    print(f"    T3 plaquette-boundary loops             : {len(t3)}   surviving: {len(t3surv)}   "
          f"INDEPENDENT modulo stars: {t3indep}")
    return indep, g, t3indep
print("G1b  RECORD COUNT ON THE GROUND SPACE\n")
analyse(2,2,"torus 2x2 (genus 1)")
print()
analyse(2,3,"torus 2x3 (genus 1)")
print()
print("  READING: if the independent count equals 2g and NO T3 object survives, then")
print("  `capacity = area - 1` never counted records, and the record count is set by GENUS.")
