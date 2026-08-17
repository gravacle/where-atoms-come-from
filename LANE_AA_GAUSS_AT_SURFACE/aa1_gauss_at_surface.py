# ASSUMPTION AUDIT 1 — is "impose the Gauss law at every vertex" valid AT THE SURFACE?
#
# W-19's headline: the redundancy plateau is a Gauss identity -- Z_l * prod_cut Z = I on the
# physical sector, so I(S:F) = H(S) for EVERY state. It could not have failed.
# THE UNEXAMINED ASSUMPTION UNDERNEATH IT: Gauss was imposed at EVERY vertex, INCLUDING the vertices
# on the cut. That is a BULK statement. The extended-Hilbert-space construction does the opposite:
# when you cut, the severed links get endpoints whose charge is LEFT FREE, and that freedom IS the
# edge mode. Imposing Gauss everywhere may be the act that deletes what we are looking for.
#
# ISOLATION: same carrier, same links, same region, same fragment, same states. ONE thing moves --
# whether the Gauss constraint is imposed at the cut vertex.
import numpy as np
rng=np.random.default_rng(20260827)
E=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]
V=6; L=len(E)
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
def kr(ops):
    M=np.array([[1]],dtype=complex)
    for o in ops: M=np.kron(M,o)
    return M
Z=[kr([Zp if j==i else I2 for j in range(L)]) for i in range(L)]
def gauss(v):
    M=np.eye(2**L,dtype=complex)
    for i,(a,b) in enumerate(E):
        if a==v or b==v: M=M@Z[i]
    return M
def projector(free_vertices):
    """impose Gauss at every vertex EXCEPT those listed as free (charge allowed there)"""
    P=np.eye(2**L,dtype=complex)
    for v in range(V):
        if v not in free_vertices: P=P@((np.eye(2**L)+gauss(v))/2)
    return P

l=0; u,vv=E[l]
cut=[i for i,(a,b) in enumerate(E) if i!=l and (a==u or b==u)]   # a u-v cut of G-l at vertex u
Xcut=np.eye(2**L,dtype=complex)
for i in cut: Xcut=Xcut@Z[i]

print(f"  carrier V={V} L={L};  system link l={l}={E[l]};  cut at vertex {u} = links {cut}")
print()
print("== ARM A — GAUSS IMPOSED EVERYWHERE (what W-19 did) ==")
PA=projector(set())
print(f"   physical dim = {int(round(np.trace(PA).real))} of {2**L}")
print(f"   || (Z_l * prod_cut Z - I) P ||  =  {np.linalg.norm((Z[l]@Xcut-np.eye(2**L))@PA):.3e}")
print()
print(f"== ARM B — GAUSS NOT IMPOSED AT THE CUT VERTEX {u} (the edge-mode prescription) ==")
PB=projector({u})
print(f"   physical dim = {int(round(np.trace(PB).real))} of {2**L}")
print(f"   || (Z_l * prod_cut Z - I) P ||  =  {np.linalg.norm((Z[l]@Xcut-np.eye(2**L))@PB):.3e}")
print("   ONE VARIABLE MOVED: whether the constraint is imposed at the cut vertex. Nothing else.")
print()

def red(rho,keep):
    d=[2]*L; r=rho.reshape(d+d); tr=[i for i in range(L) if i not in keep]
    for k,i in enumerate(sorted(tr,reverse=True)): r=np.trace(r,axis1=i,axis2=i+L-k)
    m=2**len(keep); return r.reshape(m,m)
def S(r):
    w=np.linalg.eigvalsh(r); return float(-sum(x*np.log2(x) for x in w if x>1e-12))
def MI(rho,A,B): return S(red(rho,A))+S(red(rho,B))-S(red(rho,sorted(A+B)))
def rand_state(P):
    w_,vec=np.linalg.eigh(P); basis=vec[:,w_>0.5]
    c=rng.normal(size=basis.shape[1])+1j*rng.normal(size=basis.shape[1])
    psi=basis@c; psi/=np.linalg.norm(psi); return np.outer(psi,psi.conj())

print("== THE CONSEQUENCE: is I(S:F) STILL FORCED TO H(S)? ==")
print(f"   {'':<28}{'H(S)':>12}{'I(S:F)':>12}{'I/H(S)':>10}")
for tag,P in (("ARM A  Gauss everywhere",PA),("ARM B  free at the cut",PB)):
    for k in range(3):
        rho=rand_state(P); hs=S(red(rho,[l])); mi=MI(rho,[l],cut)
        print(f"   {tag if k==0 else '':<28}{hs:>12.9f}{mi:>12.9f}{(mi/hs if hs>1e-12 else float('nan')):>10.6f}")
print()
print("   ARM A: forced to 1.000000 on every state -- the theorem W-19 found.")
print("   ARM B: if I/H(S) < 1, the identity is GONE and the plateau is no longer forced.")
print("   THE ASSUMPTION 'IMPOSE GAUSS AT THE CUT' WAS DOING THE WORK, NOT THE PHYSICS.")
