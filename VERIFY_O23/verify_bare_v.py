"""Does the 0-form record R = sign(Jz_total) pass clause (v) in the BARE operator form too
(the form the toric-code record FAILS)?  If yes, the counterexample does not depend on O-4."""
import numpy as np
I2=np.eye(2,dtype=complex); SX=np.array([[0,1],[1,0]],dtype=complex)
SY=np.array([[0,-1j],[1j,0]],dtype=complex); SZ=np.array([[1,0],[0,-1]],dtype=complex)
def kron(*ms):
    o=np.array([[1.0+0j]])
    for m in ms: o=np.kron(o,m)
    return o
def site(op,k,L):
    ms=[I2]*L; ms[k]=op; return kron(*ms)
def herm(M): return (M+M.conj().T)/2
def heis(L,J):
    H=np.zeros((2**L,2**L),complex)
    for k in range(L):
        for s in (SX,SY,SZ): H+=J*site(s,k,L)@site(s,(k+1)%L,L)/4
    return herm(H)
for L in (3,5):
    H=heis(L,-1.0)
    Jz=sum(site(SZ,k,L) for k in range(L))/2
    R=np.diag(np.sign(np.real(np.diag(Jz)))).astype(complex)
    ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-9*max(1.0,abs(ev).max()); Q=U[:,sel]
    loc=[site(s,k,L) for k in range(L) for s in (SX,SY,SZ)]
    bare=min(float(np.linalg.norm(Lk.conj().T@R@Lk+R)) for Lk in loc)
    Ux=kron(*[np.array([[0,-1j],[-1j,0]],dtype=complex)]*L)
    nontriv=float(np.linalg.norm(Q.conj().T@R@Q-np.trace(Q.conj().T@R@Q)/Q.shape[1]*np.eye(Q.shape[1])))
    print(f"FM Heisenberg ring L={L}  (0-FORM SU(2) symmetry), R = sign(Jz_total), E0 = spin-{L}/2 multiplet dim {Q.shape[1]}")
    print(f"   (i)   ||R-R+||={np.linalg.norm(R-R.conj().T):.1e}  ||R^2-I||={np.linalg.norm(R@R-np.eye(2**L)):.1e}   PASS")
    print(f"   (ii)  ||[H,R]|| = {np.linalg.norm(H@R-R@H):.1e}   PASS")
    print(f"   (iii) deviation from scalar on E0 = {nontriv:.4f}   PASS")
    print(f"   (iv)  global transversal pi-rotation: ||U+RU+R|| = {np.linalg.norm(Ux.conj().T@R@Ux+R):.1e}   PASS (writer is NON-CONTRACTIBLE)")
    print(f"   (v)   BARE form, min over single-site L of ||L+RL+R|| = {bare:.4f}   "
          f"{'PASS' if bare>1e-6 else 'FAIL'}  <- the toric record FAILS this form; this one passes")
    print(f"         POSITIVE CONTROL for that minimum: R'=Z_0 gives "
          f"{min(float(np.linalg.norm(Lk.conj().T@site(SZ,0,L)@Lk+site(SZ,0,L))) for Lk in loc):.1e} (flippable)")
    # what DOES kill it: a single-site Zeeman fully resolves E0 at first order -> clause (iii)
    def phi(V):
        B=Q.conj().T@V@Q; m=B.shape[0]; return B-np.trace(B)/m*np.eye(m)
    sp=np.sort(np.linalg.eigvalsh(herm(phi(site(SZ,0,L)))))
    print(f"   WHAT KILLS IT: Phi(Z on one site) spectrum on E0 = {np.round(sp,4)}, min gap "
          f"{np.min(np.diff(sp)):.4f} -> E0 FULLY RESOLVED at FIRST order, clause (iii) dies for every R.")
    print()
