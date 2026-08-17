# CORRECTION, RECORDED NOT PATCHED: w23_build measured [H, G_v] -- the GAUSS OPERATOR at a boundary
# vertex. Every gauge-invariant H commutes with every Gauss operator BY CONSTRUCTION, so that test
# was forced and returned 0.000e+00 no matter what. Same defect class as W-19's Gauss identity.
# THE PHYSICAL OBJECT IS THE MATTER CHARGE IN A REGION, prod_{v in R} tau_v. That is NOT a constraint
# and it is NOT automatically conserved: hopping across the region's boundary flips it.
import numpy as np
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
NL,NV=5,5; n=NL+NV
E=[(0,1),(1,2),(2,0),(0,3),(1,4)]; BULK=[0,1,2]
def q(i,P):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
Zl=[q(i,Zp) for i in range(NL)]; Xl=[q(i,Xp) for i in range(NL)]
tau=[q(NL+v,Zp) for v in range(NV)]; mu=[q(NL+v,Xp) for v in range(NV)]
def G(v):
    M=tau[v].copy()
    for i,(a,b) in enumerate(E):
        if a==v or b==v: M=M@Zl[i]
    return M
def hop(i):
    a,b=E[i]; return mu[a]@Xl[i]@mu[b]
def H(g2=1.0,J=0.8,m=0.5):
    plaq=Xl[0]@Xl[1]@Xl[2]
    return -(1.0/g2)*plaq - g2*sum(Zl) - J*sum(hop(i) for i in range(NL)) - m*sum(tau)
P=np.eye(2**n,dtype=complex)
for v in BULK: P=P@((np.eye(2**n)+G(v))/2)

print("== THE TWO OBJECTS I HAD COLLAPSED ==")
R=[0,3]                                   # a region: bulk vertex 0 and the boundary vertex it feeds
Qgauss=G(3)                               # the CONSTRAINT at a boundary vertex -- always conserved
Qmat=tau[0]@tau[3]                        # the MATTER CHARGE in the region -- the physical object
for J in (0.0,0.8):
    h=H(J=J)
    print(f"  J={J}:  || [H, GAUSS at bdry] || = {np.linalg.norm(h@Qgauss-Qgauss@h):>10.3e}"
          f"   || [H, MATTER CHARGE in R] || = {np.linalg.norm(h@Qmat-Qmat@h):>10.3e}")
print("  -> the Gauss operator is conserved at every J, forced by gauge invariance.")
print("     the MATTER CHARGE moves as soon as the hopping is switched on. THAT is the dynamical object.\n")

print("== AND IS IT DURABLE? evolve and watch the region's charge ==")
w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]
h=B.conj().T@H(J=0.8)@B; Qm=B.conj().T@Qmat@B
ev,U=np.linalg.eigh(h)
rng=np.random.default_rng(20260829)
c=rng.normal(size=B.shape[1])+1j*rng.normal(size=B.shape[1]); c/=np.linalg.norm(c)
print(f"   physical dim {B.shape[1]};  <Q_matter> over time:")
print(f"   {'t':>7}{'<Q_mat>':>14}")
for t in (0.0,0.5,1.0,2.0,5.0,10.0,25.0):
    ct=U@(np.exp(-1j*ev*t)*(U.conj().T@c))
    print(f"   {t:>7.1f}{float(np.real(ct.conj()@(Qm@ct))):>14.6f}")
print("\n   a value that moves under the theory's own dynamics is a value with a history.")
print("   a value pinned by a constraint is a label. THE FIRST ONE IS NOW ON THE TABLE.")
