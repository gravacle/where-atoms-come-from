# THE DISCRIMINATOR: is our boundary charge EM-shaped or GRAVITY-shaped?
# Both EM and gravity make a conserved charge a BOUNDARY term. That property does not distinguish
# them and cannot support the identification on its own.
# WHAT DOES: in gravity the bulk Hamiltonian is a CONSTRAINT that vanishes on shell, so the entire
# energy IS the boundary term -- charge and Hamiltonian are the same object. In EM they are
# different: charge is a surface integral, energy is a bulk density.
# TEST: on the physical subspace, is H a function of the boundary charge Q alone?
# If yes -> gravity-shaped. If H has structure Q cannot label -> EM-shaped.
import numpy as np
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
OPEN=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4),(0,6)]
n=len(OPEN); BULK=list(range(6))
Z=[op(i,Zp,n) for i in range(n)]; X=[op(i,Xp,n) for i in range(n)]
def gauss(v):
    M=np.eye(2**n,dtype=complex)
    for i,(a,b) in enumerate(OPEN):
        if a==v or b==v: M=M@Z[i]
    return M
Q=gauss(6)
P=np.eye(2**n,dtype=complex)
for v in BULK: P=P@((np.eye(2**n)+gauss(v))/2)
w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]                   # physical subspace basis
d=B.shape[1]
cycles=[[0,1,2],[3,4,5],[6,7,8]]
def wil(c):
    M=np.eye(2**n,dtype=complex)
    for i in c: M=M@X[i]
    return M
print(f"  physical subspace dimension: {d}")
for g2 in (0.3,1.0,3.0):
    H=-(1.0/g2)*sum(wil(c) for c in cycles)-g2*sum(Z)
    Hp=B.conj().T@H@B; Qp=B.conj().T@Q@B
    eh=np.round(np.linalg.eigvalsh(Hp),9); eq=np.round(np.linalg.eigvalsh(Qp),9)
    print(f"\n  g^2={g2}")
    print(f"    distinct eigenvalues of H on the physical subspace : {len(set(eh))}")
    print(f"    distinct eigenvalues of Q on the physical subspace : {len(set(eq))}")
    # can H be written as a function of Q? only if H is constant on each Q-eigenspace
    spread=[]
    for q in sorted(set(eq)):
        idx=[i for i,x in enumerate(np.round(np.linalg.eigvalsh(Qp),9)) if x==q]
        _,V=np.linalg.eigh(Qp); sub=V[:,idx]
        Hsub=sub.conj().T@Hp@sub
        ev=np.linalg.eigvalsh(Hsub)
        spread.append(ev.max()-ev.min())
    print(f"    spread of H WITHIN each Q-sector: {['%.6f'%s for s in spread]}")
    print(f"    -> H is a function of Q alone?  {all(s<1e-9 for s in spread)}")
print()
print("  READING: if H varies WITHIN a fixed Q-sector, the Hamiltonian carries structure the")
print("  boundary charge cannot label -- energy is a BULK density and charge is a surface term.")
print("  THAT IS THE EM SHAPE, NOT THE GRAVITY SHAPE. In gravity the two coincide on shell.")
