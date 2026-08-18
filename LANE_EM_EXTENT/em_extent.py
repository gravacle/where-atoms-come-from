"""How BIG is the EM object that carries one record, and how much of EM is record?

Not a philosophical question -- both are counts. Torus LxL, Z2 gauge theory.
  dim Z_1 = E - V + 1   EM's gauge-invariant configuration content (the Gauss sector)
  dim B_1 = F - 1       pure gauge -- contractible, reachable locally
  dim H_1 = 2g = 2      the records
  d       = L           minimum support of a non-trivial record operator (minimal cycle)
Verified against explicit GF(2) rank computation at small L."""
import numpy as np
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
def torus(L):
    V=L*L; ind={}; k=0; E=[]
    for j in range(L):
        for i in range(L):
            ind[('h',i,j)]=k; E.append((j*L+i, j*L+(i+1)%L)); k+=1
            ind[('v',i,j)]=k; E.append((j*L+i, ((j+1)%L)*L+i)); k+=1
    F=[[ind[('h',i,j)],ind[('v',(i+1)%L,j)],ind[('h',i,(j+1)%L)],ind[('v',i,j)]] for j in range(L) for i in range(L)]
    d1=np.zeros((V,len(E)),dtype=np.int8)
    for e,(a,b) in enumerate(E): d1[a,e]^=1; d1[b,e]^=1
    d2=np.zeros((len(E),len(F)),dtype=np.int8)
    for f,pl in enumerate(F):
        for e in pl: d2[e,f]^=1
    return V,len(E),len(F),d1,d2

print("VERIFICATION -- closed forms against explicit GF(2) ranks")
print(f"  {'L':>3}{'links':>7}{'dim Z_1':>10}{'formula':>9}{'dim B_1':>10}{'formula':>9}{'dim H_1':>10}")
for L in (2,3,4,5):
    V,E,F,d1,d2=torus(L); r1,r2=rank2(d1),rank2(d2)
    z,b=E-r1,r2
    print(f"  {L:>3}{E:>7}{z:>10}{L*L+1:>9}{b:>10}{L*L-1:>9}{z-b:>10}"
          + ("   PASS" if (z==L*L+1 and b==L*L-1 and z-b==2) else "   FAIL"))
print()
print("HOW MUCH OF EM IS RECORD, AND HOW BIG IS THE OBJECT THAT CARRIES IT")
print(f"  {'L':>4}{'links':>8}{'dim Z_1':>10}{'dim H_1':>9}{'record fraction':>18}{'d (support)':>13}{'d / links':>11}")
for L in (2,3,4,8,16,64,256):
    E=2*L*L; z=L*L+1
    print(f"  {L:>4}{E:>8}{z:>10}{2:>9}{2/z:>18.6e}{L:>13}{L/E:>11.6f}")
print()
print("  RECORD FRACTION      2/(L^2+1)  -> 0     the record is a VANISHING fraction of EM's content")
print("  ABSOLUTE EXTENT      d = L      -> inf   and the object carrying it grows WITHOUT BOUND")
print("  RELATIVE EXTENT      d/links = 1/(2L) -> 0")
print()
print("  So the record-carrying EM object is neither local nor the whole system:")
print("  unboundedly large in absolute size, vanishing as a fraction, and with")
print("  ZERO local content -- no local operator reads it (Thm C, 1.59e-16).")
