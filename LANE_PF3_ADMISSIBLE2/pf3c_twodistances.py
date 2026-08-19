"""Is PF-3B's failure the DEFINITION's or the CARRIER's?

G-12/G-15 already say R3 needs TWO numbers: the code distance is min(systole, cosystole).
G-11's table reported d = 3,3,3 for the bouquet -- but that is the Z-SYSTOLE alone. If the
X-COSYSTOLE is 1, the bouquet has d = 1 and genuinely fails clause (v), and DEF-A is
exonerated. Computed here for the bouquet AND the torus, over GF(2)."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
def rank2(M):
    M=M.copy()%2; rows,cols=M.shape; r=0
    for c in range(cols):
        p=next((i for i in range(r,rows) if M[i,c]),None)
        if p is None: continue
        M[[r,p]]=M[[p,r]]
        for i in range(rows):
            if i!=r and M[i,c]: M[i]^=M[r]
        r+=1
    return r
def inspan(v,B):
    B=[b for b in B]
    for b in B:
        if v==0: break
        h=v.bit_length()-1
        if b.bit_length()-1==h: v^=b
    return v==0
def basis_of(rows):
    piv={}
    for v in rows:
        while v:
            h=v.bit_length()-1
            if h in piv: v^=piv[h]
            else: piv[h]=v; break
    return list(piv.values())
def minwt(space_rows, sub_rows, nL):
    """min Hamming weight of an element of span(space) NOT in span(sub)."""
    S=basis_of(space_rows); B=basis_of(sub_rows)
    best=None
    for m in range(1,1<<len(S)):
        v=0
        for i in range(len(S)):
            if (m>>i)&1: v^=S[i]
        if v and not inspan(v,B):
            w=bin(v).count('1')
            if best is None or w<best: best=w
    return best

def analyse(name,V,E,PLQ):
    nL=len(E)
    d1=np.zeros((V,nL),dtype=np.int8)
    for k,(a,b) in enumerate(E): d1[a,k]^=1; d1[b,k]^=1
    d2=np.zeros((nL,len(PLQ)),dtype=np.int8)
    for k,p in enumerate(PLQ):
        for e in p: d2[e,k]^=1
    def toint(v): return int(''.join(map(str,v)),2)
    # Z-logicals: cycles (ker d1) modulo boundaries (im d2)   -> Z-systole
    def nullspace(M):
        M=M.copy()%2; rows,cols=M.shape; pc=[]; r=0
        for c in range(cols):
            p=next((i for i in range(r,rows) if M[i,c]),None)
            if p is None: continue
            M[[r,p]]=M[[p,r]]
            for i in range(rows):
                if i!=r and M[i,c]: M[i]^=M[r]
            pc.append(c); r+=1
        out=[]
        for fc in [c for c in range(cols) if c not in pc]:
            v=np.zeros(cols,dtype=np.int8); v[fc]=1
            for i,p_ in enumerate(pc): v[p_]=M[i,fc]
            out.append(v)
        return out
    Z1=[toint(v) for v in nullspace(d1)]
    B1=[toint(d2[:,c]) for c in range(d2.shape[1])] or [0]
    syst=minwt(Z1,B1,nL)
    # X-logicals: cocycles (ker d2^T) modulo vertex stars (row space d1) -> cosystole
    Zc=[toint(v) for v in nullspace(d2.T)] if len(PLQ) else [1<<i for i in range(nL)]
    Bc=[toint(d1[r,:]) for r in range(V)]
    cosy=minwt(Zc,Bc,nL)
    d=min(x for x in (syst,cosy) if x is not None)
    say(f"  {name:<32}{nL:>6}{str(syst):>12}{str(cosy):>12}{d:>10}")
    return syst,cosy,d

say("="*88); say("THE TWO DISTANCES -- systole AND cosystole, over GF(2)"); say("="*88)
say(f"  {'carrier':<32}{'links':>6}{'Z-systole':>12}{'X-cosystole':>12}{'d = min':>10}")
def torus(n):
    ind={}; E=[]; k=0
    for j in range(n):
        for i in range(n):
            ind[('h',i,j)]=k; E.append((j*n+i,j*n+(i+1)%n)); k+=1
            ind[('v',i,j)]=k; E.append((j*n+i,((j+1)%n)*n+i)); k+=1
    P=[[ind[('h',i,j)],ind[('v',(i+1)%n,j)],ind[('h',i,(j+1)%n)],ind[('v',i,j)]] for j in range(n) for i in range(n)]
    return n*n,E,P
for n in (2,3):
    V,E,P=torus(n); analyse(f"torus {n}x{n} (manifold)",V,E,P)
for k in (2,3):
    E=[]; v=1
    for _ in range(k):
        a,b=v,v+1; v+=2; E+=[(0,a),(a,b),(b,0)]
    analyse(f"bouquet of {k} triangles (pinch)",v,E,[])
say("")
say("  READ: G-11's table reported d = 3,3,3 for the bouquet. That is the Z-SYSTOLE alone.")
say("        If the X-cosystole is 1 the true code distance is 1, the bouquet fails clause (v)")
say("        on its own merits, and PF-3B's 171 admissible flippers indict the CARRIER, not DEF-A.")
