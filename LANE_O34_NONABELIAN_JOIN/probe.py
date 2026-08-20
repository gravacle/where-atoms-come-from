"""Probe: can the model enumerate records on D(D_4), dim 64?  O-28 caps at 20 minimal projections."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
def close(gens,d):
    def m(a,b): return tuple(a[b[i]] for i in range(d))
    E=tuple(range(d)); S={E}; fr=[E]
    while fr:
        nf=[]
        for x in fr:
            for g in gens:
                y=m(x,g)
                if y not in S: S.add(y); nf.append(y)
        fr=nf
    return sorted(S)
d=4; G=close([(1,2,3,0),(1,0,3,2)],d); n=len(G); gi={g:i for i,g in enumerate(G)}
def mul(a,b): return tuple(a[b[i]] for i in range(d))
def inv(a):
    r=[0]*d
    for i,x in enumerate(a): r[x]=i
    return tuple(r)
e=tuple(range(d)); D=n*n
def ket(g1,g2): return gi[g1]*n+gi[g2]
A=np.zeros((D,D))
for h in G:
    for g1 in G:
        for g2 in G:
            A[ket(mul(mul(h,g1),inv(h)),mul(mul(h,g2),inv(h))), ket(g1,g2)]+=1.0/n
B=np.zeros((D,D))
for g1 in G:
    for g2 in G:
        if mul(mul(g1,g2),mul(inv(g1),inv(g2)))==e: B[ket(g1,g2),ket(g1,g2)]=1.0
H=-(A+B)
M=RecordModel(H.astype(complex))
say(f"  D(D_4) minimal torus: |G|={n}, dim {D}")
say(f"  eigenvalues {[round(float(v),6) for v,_,_ in M.es]}  multiplicities {[m for _,_,m in M.es]}")
say(f"  ground space {M.ground_space()[1]}")
say(f"  minimal projections in the commutant: {len(M.projs)}   (O-28 caps enumeration at 20)")
say(f"  non-commuting pairs in D_4: {sum(1 for a in G for b in G if mul(a,b)!=mul(b,a))} of {n*n}"
    f"   -> flux transport is {'NON-TRIVIAL' if any(mul(a,b)!=mul(b,a) for a in G for b in G) else 'trivial'}")
