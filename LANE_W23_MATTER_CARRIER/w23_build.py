# W-23 — THE TWO-DANGLING-LINK CARRIER WITH MATTER. Attempt two, corrected construction.
# WHY: pure gauge freezes the boundary charge ([H,Q] = 0 exactly, verified) because a dangling link
# lies on no cycle and no plaquette term can contain it. Only MATTER can move charge onto a boundary.
# And one dangling link is not enough: the bulk constraints then FORCE the boundary charge.
#
# CARRIER: triangle {0,1,2} in the bulk, two dangling links to boundary vertices 3 and 4.
#   links  0:(0,1) 1:(1,2) 2:(2,0)   [the plaquette]   3:(0,3) 4:(1,4)   [dangling]
#   Z_2 matter at every vertex.  Gauss imposed at BULK vertices only.
#   qubits = 5 links + 5 matter = 10  ->  dim 1024, exactly diagonalisable.
import numpy as np
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
NL,NV=5,5; n=NL+NV                       # qubits 0..4 links, 5..9 matter at vertices 0..4
E=[(0,1),(1,2),(2,0),(0,3),(1,4)]
BULK=[0,1,2]; BDRY=[3,4]
def q(i,P):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
Zl=[q(i,Zp) for i in range(NL)]; Xl=[q(i,Xp) for i in range(NL)]
tau=[q(NL+v,Zp) for v in range(NV)]      # matter charge (diagonal)
mu =[q(NL+v,Xp) for v in range(NV)]      # matter flip
def G(v):
    M=tau[v].copy()
    for i,(a,b) in enumerate(E):
        if a==v or b==v: M=M@Zl[i]
    return M
print(f"  carrier: V={NV} L={NL}, matter at every vertex, {n} qubits, dim {2**n}")
print(f"  Gauss operators are products of the matter charge and the incident electric fields.")

# hopping term: flips matter at both ends and the link between -> gauge invariant, MOVES CHARGE
def hop(i):
    a,b=E[i]; return mu[a]@Xl[i]@mu[b]
print(f"\n  gauge invariance of the hopping term (max commutator with any bulk Gauss):")
print(f"    max || [hop_l, G_v] || = {max(np.linalg.norm(hop(i)@G(v)-G(v)@hop(i)) for i in range(NL) for v in BULK):.3e}")

P=np.eye(2**n,dtype=complex)
for v in BULK: P=P@((np.eye(2**n)+G(v))/2)
dphys=int(round(np.trace(P).real))
print(f"  physical dim (Gauss imposed on the 3 bulk vertices only) = {dphys} of {2**n}")

Q3,Q4=G(3),G(4)
w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]
e3=sorted(set(np.round(np.linalg.eigvalsh(B.conj().T@Q3@B),6)))
e34=sorted(set(np.round(np.linalg.eigvalsh(B.conj().T@(Q3@Q4)@B),6)))
print(f"\n  boundary charge Q_3 on the physical subspace : {e3}")
print(f"  total Q_3*Q_4 on the physical subspace       : {e34}")

def H(g2=1.0,J=0.8,m=0.5):
    plaq=Xl[0]@Xl[1]@Xl[2]
    return -(1.0/g2)*(plaq+plaq.conj().T)/1.0 - g2*sum(Zl) - J*sum(hop(i) for i in range(NL)) - m*sum(tau)
print(f"\n  == THE TEST PURE GAUGE FAILED: DOES THE CHARGE MOVE? ==")
for J in (0.0,0.8):
    h=H(J=J)
    print(f"    J={J}:  || [H, Q_3] || = {np.linalg.norm(h@Q3-Q3@h):.3e}"
          f"   || [H, Q_3*Q_4] || = {np.linalg.norm(h@(Q3@Q4)-(Q3@Q4)@h):.3e}")
print("    J is the matter hopping. At J=0 this is pure gauge and the charge is frozen.")
