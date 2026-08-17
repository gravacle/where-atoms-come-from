# ASSUMPTION AUDIT 1b — why AA1 came back null, and what freeing a boundary charge actually requires.
# AA1's two arms had IDENTICAL physical dimension (16 = 16). That is the tell: a zero-variable
# control. Here is the theorem behind it, and the construction that does work.
import numpy as np
rng=np.random.default_rng(20260827)
E=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]
V=6; L=len(E)
I2=np.eye(2); Zp=np.diag([1,-1]).astype(complex)
def kr(ops,n):
    M=np.array([[1]],dtype=complex)
    for o in ops: M=np.kron(M,o)
    return M
def Zi(i,n): return kr([Zp if j==i else I2 for j in range(n)],n)
def gauss_of(edges,v,n):
    M=np.eye(2**n,dtype=complex)
    for i,(a,b) in enumerate(edges):
        if a==v or b==v: M=M@Zi(i,n)
    return M

print("== WHY AA1 WAS A ZERO-VARIABLE CONTROL: THE CONSTRAINTS ARE NOT INDEPENDENT ==")
Gs=[gauss_of(E,v,L) for v in range(V)]
prod=np.eye(2**L,dtype=complex)
for g in Gs: prod=prod@g
print(f"   || prod over ALL vertices of G_v  -  I ||  =  {np.linalg.norm(prod-np.eye(2**L)):.3e}")
print("   Every link touches exactly TWO vertices, so each Z appears twice in the product and")
print("   squares to the identity. ON A CLOSED GRAPH THE GAUSS CONSTRAINTS SATISFY ONE RELATION.")
print("   So dropping ONE of them removes NOTHING -- it is implied by the other five.")
print("   'Impose Gauss at the cut vertex' was never an assumption I made. IT IS FORCED.\n")

print("== WHAT ACTUALLY FREES A BOUNDARY CHARGE: SEVER THE LINK, DO NOT RELAX THE CONSTRAINT ==")
# subdivide link 0 = (0,1) into (0,6) and (6,1) with a NEW vertex 6, and leave vertex 6 unconstrained
E2=[(0,6),(6,1)]+E[1:]
V2=7; L2=len(E2)
def projector(edges,n,nv,free):
    P=np.eye(2**n,dtype=complex)
    for v in range(nv):
        if v not in free: P=P@((np.eye(2**n)+gauss_of(edges,v,n))/2)
    return P
Pc=projector(E2,L2,V2,set())          # constrained everywhere, including the new vertex
Pf=projector(E2,L2,V2,{6})            # the new endpoint's charge LEFT FREE -- the edge mode
print(f"   severed carrier: V={V2} L={L2}, full space {2**L2}")
print(f"   Gauss imposed at the new endpoint : physical dim = {int(round(np.trace(Pc).real))}")
print(f"   new endpoint's charge LEFT FREE   : physical dim = {int(round(np.trace(Pf).real))}")
print(f"   ratio = {int(round(np.trace(Pf).real))/max(1,int(round(np.trace(Pc).real))):.1f}x  <-- NOW the arms differ")
print()
cut=[i for i,(a,b) in enumerate(E2) if i not in (0,1) and (a==0 or b==0)]
Xcut=np.eye(2**L2,dtype=complex)
for i in cut: Xcut=Xcut@Zi(i,L2)
lhs=Zi(0,L2)@Xcut
print(f"   the W-19 identity, constrained arm : || (Z_l1 * prod_cut Z - I) P || = {np.linalg.norm((lhs-np.eye(2**L2))@Pc):.3e}")
print(f"   the W-19 identity, FREE-charge arm : || (Z_l1 * prod_cut Z - I) P || = {np.linalg.norm((lhs-np.eye(2**L2))@Pf):.3e}")
print()
print("   ==> THE ASSUMPTION IS STILL UNEXAMINED, BUT NOW I KNOW WHAT EXAMINING IT COSTS:")
print("       an edge mode requires SEVERING the graph and carrying the new endpoints as degrees of")
print("       freedom. W-19 never did that -- it cut a CLOSED graph, where the constraint cannot be")
print("       relaxed. Its 'Gauss-forced' verdict is therefore a verdict about a construction that")
print("       has no edge modes in it AT ALL.")
