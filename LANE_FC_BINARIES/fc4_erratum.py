# ERRATUM against W-22 and against fc3, found by fc3's own printed diagnostic.
# fc3 reported "distinct eigenvalues of Q on the physical subspace : 1". Q is CONSTANT there, so
# "H is not a function of Q" could not have failed. The discriminator is VOID as run.
# And W-22's claim of TWO superselection sectors is wrong ON THE PHYSICAL SUBSPACE. Here is why.
import numpy as np
I2=np.eye(2); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
def build(edges,bulk):
    n=len(edges)
    Z=[op(i,Zp,n) for i in range(n)]
    def gauss(v):
        M=np.eye(2**n,dtype=complex)
        for i,(a,b) in enumerate(edges):
            if a==v or b==v: M=M@Z[i]
        return M
    P=np.eye(2**n,dtype=complex)
    for v in bulk: P=P@((np.eye(2**n)+gauss(v))/2)
    return n,gauss,P

print("== ONE DANGLING LINK: the bulk constraints FORCE the boundary charge ==")
E1=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4),(0,6)]
n1,g1,P1=build(E1,list(range(6)))
prod=np.eye(2**n1,dtype=complex)
for v in range(6): prod=prod@g1(v)
Z9=op(9,Zp,n1)
print(f"   || prod_bulk G_v  -  Z_9 || = {np.linalg.norm(prod-Z9):.3e}")
print("   Each link with BOTH ends in the bulk contributes Z twice and cancels. The dangling link")
print("   contributes ONCE. So prod_bulk G_v = Z_9 = Q exactly, and imposing bulk Gauss = +1")
print("   FORCES Q = +1. ONE SECTOR, NOT TWO.")
Qp=P1@g1(6)@P1
ev=np.round(np.linalg.eigvalsh(Qp+1e-9*(np.eye(2**n1)-P1)),6)
print(f"   physical dim {int(round(np.trace(P1).real))};  distinct eigenvalues of Q there: "
      f"{sorted(set(np.round(np.linalg.eigvalsh((P1@g1(6)@P1)[np.ix_(*[np.argsort(-np.linalg.eigvalsh(P1))[:0]]*2)] if False else P1@g1(6)@P1),6)))-{0.0} if False else 'see below'}")
w_,vec=np.linalg.eigh(P1); B=vec[:,w_>0.5]
print(f"   eigenvalues of Q restricted to the physical subspace: "
      f"{sorted(set(np.round(np.linalg.eigvalsh(B.conj().T@g1(6)@B),6)))}")

print("\n== TWO DANGLING LINKS: now the boundary charge is FREE ==")
E2=E1+[(1,7)]                                            # second dangling link at a new boundary vertex
n2,g2f,P2=build(E2,list(range(6)))
w2,vec2=np.linalg.eigh(P2); B2=vec2[:,w2>0.5]
Qtot=g2f(6)@g2f(7)                                       # total boundary charge
print(f"   physical dim {int(round(np.trace(P2).real))} of {2**n2}")
print(f"   eigenvalues of Q_6 on the physical subspace     : "
      f"{sorted(set(np.round(np.linalg.eigvalsh(B2.conj().T@g2f(6)@B2),6)))}")
print(f"   eigenvalues of Q_6*Q_7 (total) on physical      : "
      f"{sorted(set(np.round(np.linalg.eigvalsh(B2.conj().T@Qtot@B2),6)))}")
print("\n   -> with TWO boundary links the individual boundary charge takes BOTH values while the")
print("      TOTAL is still forced. That is the real structure: boundary charge is free to be")
print("      DISTRIBUTED, and only its sum is fixed by the bulk. ONE dangling link cannot show it.")
