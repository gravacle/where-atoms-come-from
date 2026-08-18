"""What a holonomy is, measured. Z2 gauge field on a 3x3 torus: each link carries a bit.
A loop's HOLONOMY = XOR of the bits on it.

CONTRACTIBLE loop  -> its holonomy EQUALS the XOR of the plaquettes it encloses (discrete
                      Stokes). It is therefore determined by what is inside. Not new information.
NON-CONTRACTIBLE   -> there is no inside. No combination of plaquettes reproduces it.
                      It is information that exists nowhere locally."""
import numpy as np, itertools
L=3; V=L*L; ind={}; k=0; E=[]
for j in range(L):
    for i in range(L):
        ind[('h',i,j)]=k; E.append((j*L+i,j*L+(i+1)%L)); k+=1
        ind[('v',i,j)]=k; E.append((j*L+i,((j+1)%L)*L+i)); k+=1
nE=len(E)
PLAQ=[[ind[('h',i,j)],ind[('v',(i+1)%L,j)],ind[('h',i,(j+1)%L)],ind[('v',i,j)]] for j in range(L) for i in range(L)]
def vec(links):
    v=np.zeros(nE,dtype=np.int8)
    for l in links: v[l]^=1
    return v
P=np.array([vec(p) for p in PLAQ]).T                       # plaquette boundaries, columns

def in_span(target,basis):
    """is target an XOR of some subset of the basis columns? GF(2) solve."""
    A=np.concatenate([basis,target.reshape(-1,1)],axis=1)%2
    rows,cols=A.shape; r=0
    for c in range(cols-1):
        p=next((i for i in range(r,rows) if A[i,c]),None)
        if p is None: continue
        A[[r,p]]=A[[p,r]]
        for i in range(rows):
            if i!=r and A[i,c]: A[i]^=A[r]
        r+=1
    for i in range(rows):
        if A[i,-1] and not A[i,:-1].any(): return False
    return True

print("Z2 gauge field, 3x3 torus:", nE, "links,", len(PLAQ), "plaquettes\n")
loops={
 "one plaquette (contractible)"      : PLAQ[0],
 "2x2 block boundary (contractible)" : [ind[('h',0,0)],ind[('h',1,0)],ind[('v',2,0)],ind[('v',2,1)],
                                        ind[('h',1,2)],ind[('h',0,2)],ind[('v',0,1)],ind[('v',0,0)]],
 "horizontal wrap (NON-contractible)": [ind[('h',i,0)] for i in range(L)],
 "vertical wrap (NON-contractible)"  : [ind[('v',0,j)] for j in range(L)],
}
print(f"  {'loop':<38}{'length':>7}{'closed?':>9}{'= XOR of plaquettes?':>23}")
for nm,lk in loops.items():
    v=vec(lk)
    d=np.zeros(V,dtype=np.int8)
    for l in lk:
        a,b=E[l]; d[a]^=1; d[b]^=1
    print(f"  {nm:<38}{len(lk):>7}{('yes' if not d.any() else 'NO'):>9}"
          f"{('YES - determined inside' if in_span(v,P) else 'NO - nothing determines it'):>23}")

print("\n  How many independent loop-holonomies are NOT determined by any plaquettes?")
Z1=[]
M=np.zeros((V,nE),dtype=np.int8)
for l,(a,b) in enumerate(E): M[a,l]^=1; M[b,l]^=1
A=M.copy()%2; rows,cols=A.shape; pc=[]; r=0
for c in range(cols):
    p=next((i for i in range(r,rows) if A[i,c]),None)
    if p is None: continue
    A[[r,p]]=A[[p,r]]
    for i in range(rows):
        if i!=r and A[i,c]: A[i]^=A[r]
    pc.append(c); r+=1
free=[c for c in range(cols) if c not in pc]
for fc in free:
    v=np.zeros(cols,dtype=np.int8); v[fc]=1
    for i,p_ in enumerate(pc): v[p_]=A[i,fc]
    Z1.append(v)
undet=sum(1 for m in range(1,1<<len(Z1))
          if not in_span(np.bitwise_xor.reduce([Z1[i] for i in range(len(Z1)) if (m>>i)&1]),P))
det=2**len(Z1)-undet
print(f"    closed loops in total                    : 2^{len(Z1)} = {2**len(Z1)}")
print(f"    of those, DETERMINED by plaquettes       : {det} = 2^{int(np.log2(det))}   (these are B_1)")
print(f"    loops NOT determined by any plaquettes   : {undet}")
print(f"    but they fall into CLASSES -- two loops differing by plaquettes carry the SAME holonomy:")
print(f"      classes = {2**len(Z1)} / {det} = {2**len(Z1)//det} = 2^dim H_1,  so dim H_1 = {int(np.log2(2**len(Z1)//det))}")
print(f"\n  dim H_1 = 2.  Two independent holonomies nothing local can fix.  THAT is the record.")
