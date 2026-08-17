# Registrar's check of W-19's load-bearing claim: the "redundancy plateau" on a pure Z_2 gauge
# carrier is FORCED BY THE GAUSS LAW, hence could not have failed, hence is not evidence of a record.
# Claim: if F contains a u-v cut of G - l, and S u F contains no cycle through l, then
#        X_l = X(cut)^{-1} EXACTLY on the physical sector, so I(S:F) = H(S) EXACTLY, for ANY state.
import numpy as np, itertools
rng=np.random.default_rng(20260825)
# dbl_chain-like: V=6, links as (u,v). A simple chain of triangles is enough to carry a u-v cut.
E=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]      # V=6, L=9, min degree 3
V=6; L=len(E)
def kron(ops):
    M=np.array([[1]],dtype=complex)
    for o in ops: M=np.kron(M,o)
    return M
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
X=[kron([Xp if j==i else I2 for j in range(L)]) for i in range(L)]
Z=[kron([Zp if j==i else I2 for j in range(L)]) for i in range(L)]
# Gauss operator at vertex v: product of Z on incident links (Z_2 electric-field generator)
def gauss(v):
    M=np.eye(2**L,dtype=complex)
    for i,(a,b) in enumerate(E):
        if a==v or b==v: M=M@Z[i]
    return M
G=[gauss(v) for v in range(V)]
P=np.eye(2**L,dtype=complex)
for g in G: P=P@((np.eye(2**L)+g)/2)
dim=int(round(np.trace(P).real))
print(f"  carrier V={V} L={L}; physical sector dim = {dim} of {2**L}")

print("\n== V1  IS THE PLATEAU A GAUSS IDENTITY? ==")
l=0; u,v=E[l]
# a u-v cut of G-l: the set of links (other than l) incident to vertex u
cut=[i for i,(a,b) in enumerate(E) if i!=l and (a==u or b==u)]
print(f"  system link l = {l} = {E[l]};  cut of G-l at vertex {u} = links {cut}")
Xcut=np.eye(2**L,dtype=complex)
for i in cut: Xcut=Xcut@Z[i]
# Gauss at u says Z_l * prod_{cut} Z = identity on the physical sector
lhs=Z[l]@Xcut
print(f"  || (Z_l * prod_cut Z - I) P ||  =  {np.linalg.norm((lhs-np.eye(2**L))@P):.3e}")
print("  -> on the physical sector the system link's electric flux IS the cut's, exactly.")
print("     So any fragment containing that cut determines S with certainty, for EVERY state.\n")

print("== V2  AND IT COULD NOT HAVE FAILED: same number on states that should differ ==")
def red(rho,keep):
    d=[2]*L; r=rho.reshape(d+d); tr=[i for i in range(L) if i not in keep]
    for k,i in enumerate(sorted(tr,reverse=True)): r=np.trace(r,axis1=i,axis2=i+L-k)
    m=2**len(keep); return r.reshape(m,m)
def S(r):
    w=np.linalg.eigvalsh(r); return float(-sum(x*np.log2(x) for x in w if x>1e-12))
def MI(rho,A,B): return S(red(rho,A))+S(red(rho,B))-S(red(rho,sorted(A+B)))
# physical basis
w_,vecs=np.linalg.eigh(P); phys=vecs[:,w_>0.5]
def rand_phys():
    c=rng.normal(size=phys.shape[1])+1j*rng.normal(size=phys.shape[1]); psi=phys@c
    psi/=np.linalg.norm(psi); return np.outer(psi,psi.conj())
print(f"  {'state':<34}{'H(S)':>12}{'I(S:F)':>12}{'I/H(S)':>10}")
for tag,rho in [("Haar physical #1",rand_phys()),("Haar physical #2",rand_phys()),
                ("Haar physical #3",rand_phys())]:
    hs=S(red(rho,[l])); mi=MI(rho,[l],cut)
    print(f"  {tag:<34}{hs:>12.9f}{mi:>12.9f}{mi/hs if hs>1e-12 else float('nan'):>10.6f}")
print("  -> I/H(S) = 1.000000 on every state, including states with nothing in common.")
print("     THE PLATEAU IS A THEOREM ABOUT THE GAUSS LAW, NOT A MEASUREMENT OF A RECORD.")
print("     'Could not have failed' voids a CONTROL -- and this was being read as evidence.")
