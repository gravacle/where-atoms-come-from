"""T-8: a SECOND CARRIER FAMILY, and O-14's specific need.

PF-3B failed because the bouquet of triangles has Z-systole 3 but X-COSYSTOLE 1, so its true
code distance is 1 and it cannot test clause (v) at all. O-14 asks for a NON-MANIFOLD with BOTH
distances >= 2. The bouquet had no 2-cells, which is exactly why its cosystole collapsed -- so
the candidates here all HAVE 2-cells and are non-manifold in some other way."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
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
def toint(v): return int(''.join(map(str,v)),2)
def basis_of(rows):
    piv={}
    for v in rows:
        while v:
            h=v.bit_length()-1
            if h in piv: v^=piv[h]
            else: piv[h]=v; break
    return piv
def inspan(v,piv):
    for h in sorted(piv,reverse=True):
        if (v>>h)&1: v^=piv[h]
    return v==0
def minwt(space,subpiv,cap=1<<20):
    S=list(basis_of(space).values())
    if len(S)>20: return None
    best=None
    for mm in range(1,1<<len(S)):
        v=0
        for i in range(len(S)):
            if (mm>>i)&1: v^=S[i]
        if v and not inspan(v,subpiv):
            w=bin(v).count('1')
            if best is None or w<best: best=w
    return best

def analyse(name,V,E,F,note=""):
    nL=len(E)
    d1=np.zeros((V,nL),dtype=np.int8)
    for k,(a,b) in enumerate(E): d1[a,k]^=1; d1[b,k]^=1
    d2=np.zeros((nL,len(F)),dtype=np.int8)
    for k,p in enumerate(F):
        for e in p: d2[e,k]^=1
    assert not ((d1@d2)%2).any(), name+": d1 d2 != 0"
    h1=(nL-rank2(d1))-rank2(d2)
    Z1=[int(toint(v)) for v in nullspace(d1)]
    B1=basis_of([int(toint(d2[:,c])) for c in range(d2.shape[1])] or [0])
    syst=minwt(Z1,B1)
    Zc=[int(toint(v)) for v in nullspace(d2.T)] if len(F) else [1<<i for i in range(nL)]
    Bc=basis_of([int(toint(d1[r,:])) for r in range(V)])
    cosy=minwt(Zc,Bc)
    d = min(x for x in (syst,cosy) if x is not None) if (syst and cosy) else (syst or cosy)
    ok = (syst or 0)>=2 and (cosy or 0)>=2 and h1>=1
    say(f"  {name:<38}{V:>4}{nL:>5}{len(F):>4}{h1:>6}{str(syst):>9}{str(cosy):>10}{str(d):>5}   "
        f"{'QUALIFIES' if ok else ''}  {note}")
    return ok,h1,syst,cosy

say("="*104); say("T-8 / O-14   A NON-MANIFOLD CARRIER WITH BOTH DISTANCES >= 2"); say("="*104)
say(f"  {'candidate':<38}{'V':>4}{'E':>5}{'F':>4}{'H_1':>6}{'systole':>9}{'cosyst':>10}{'d':>5}")

# torus 2x2 -- the incumbent, a manifold, for reference
def torus(n):
    ind={}; E=[]; k=0
    for j in range(n):
        for i in range(n):
            ind[('h',i,j)]=k; E.append((j*n+i,j*n+(i+1)%n)); k+=1
            ind[('v',i,j)]=k; E.append((j*n+i,((j+1)%n)*n+i)); k+=1
    F=[[ind[('h',i,j)],ind[('v',(i+1)%n,j)],ind[('h',i,(j+1)%n)],ind[('v',i,j)]] for j in range(n) for i in range(n)]
    return n*n,E,F
analyse("torus 2x2 (MANIFOLD, incumbent)",*torus(2))
analyse("bouquet of 2 triangles (PF-3B's)",5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],[],"cosystole 1 -- disqualified")

# TORUS + A DUPLICATE PLAQUETTE: those four edges now lie in THREE faces, so the complex is
# NOT a manifold anywhere along them, while H_1 is untouched (the new face bounds nothing new).
V,E,F=torus(2)
analyse("torus 2x2 + duplicated plaquette",V,E,F+[F[0]],"edge in 3 faces -- NON-MANIFOLD")
V,E,F=torus(2)
analyse("torus 2x2 + two duplicated plaquettes",V,E,F+[F[0],F[1]],"NON-MANIFOLD")
V,E,F=torus(3)
analyse("torus 3x3 + duplicated plaquette",V,E,F+[F[0]],"NON-MANIFOLD")

# TORUS + A DISK CAPPING A NON-CONTRACTIBLE CYCLE. links 0 and 2 are h(0,0) and h(1,0): a wrap.
V,E,F=torus(2)
analyse("torus 2x2 + disk capping a wrap",V,E,F+[[0,2]],"non-manifold along the capped cycle")

# BOOK: k square pages sharing one spine edge -> non-manifold along the spine
# BOOK: k square pages sharing one spine edge -> non-manifold along the spine
def book(pages, cyc=4):
    E=[(0,1)]; F=[]                      # the spine
    v=2
    for _ in range(pages):
        prev=1; ids=[0]
        for s in range(cyc-2):
            E.append((prev,v)); ids.append(len(E)-1); prev=v; v+=1
        E.append((prev,0)); ids.append(len(E)-1)
        F.append(ids)
    return v,E,F
for p in (2,3):
    analyse(f"book, {p} square pages on one spine",*book(p),"non-manifold along the spine")

say("")
say("  A TRAP, RECORDED. 'torus + duplicated plaquette' is non-manifold as a COMPLEX and is NOT a")
say("  second carrier: adding a face that is already there leaves rank(d2), B_1, H_1 and the whole")
say("  STABILISER GROUP unchanged. Same code, different drawing. Only a face that changes the")
say("  stabiliser group changes the carrier.")
say("")
say("  THE SECOND CARRIER: torus 2x2 + a disk capping a non-contractible wrap.")
say("    non-manifold along the capped cycle (those edges lie in THREE faces),")
say("    H_1 = 1 where the torus has 2 -- so it is [[8,1,2]] against the toric [[8,2,2]],")
say("    both distances 2, and the same Hilbert dimension 256.")
