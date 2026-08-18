"""SUPERSEDED AND WITHDRAWN -- the anticommutation test in this draft had a
null-space bug (it accepted numerically-zero combinations).  Its line
'exists single-site combination ANTICOMMUTING with R = True' is WRONG.
Use verify_dR_0form2.py, which builds the image space explicitly.
Kept only so the error is on the record."""
"""PROBE 7.  The lane proposes to RECLASSIFY C-2 as 'a 0-form symmetry has d_R = 1'.
What T5 actually proves is D1 > 0 (some local term of the CHARGE splits the multiplet at first
order).  By the lane's OWN O-3 result (CASE C/D), splitting != record death, so D1 > 0 does NOT
imply d_R = 1.  Test it directly on a genuine 0-form (SU(2)) multiplet."""
import numpy as np, itertools
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
def phi(V,Q):
    B=Q.conj().T@V@Q; m=B.shape[0]
    return B-np.trace(B)/m*np.eye(m)

for L in (3,4,5):
    H=heis(L,-1.0)                      # ferromagnet: ground multiplet = spin L/2, dim L+1
    ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-9*max(1.0,abs(ev).max()); Q=U[:,sel]; n=Q.shape[1]
    Jz=sum(site(SZ,k,L) for k in range(L))/2
    Jzb=Q.conj().T@Jz@Q; Jzb=herm(Jzb)
    w,Ub=np.linalg.eigh(Jzb)
    if n%2: 
        print(f"  L={L}: multiplet dim n={n} is ODD -> any involution on E0 has unequal +-1 counts,")
        print(f"        so clause (iv) (U+RU=-R preserves trace) is IMPOSSIBLE.  No record on E0 at all.")
        continue
    # R = sign(Jz) in the multiplet basis: balanced +-1, so trace 0 -> clause (iv) is possible
    R0=Ub@np.diag(np.sign(w))@Ub.conj().T
    R0=herm(R0)
    print(f"  L={L}: multiplet dim n={n}, R = sign(Jz) on E0, eigenvalues {np.round(np.linalg.eigvalsh(R0),6)}")
    # the space of Phi(V) reachable by SINGLE-SITE Hermitian V
    basis=[]
    for k in range(L):
        for s in (SX,SY,SZ):
            basis.append(phi(site(s,k,L),Q))
    # (a) does ANY single-site V give [Phi(V),R] != 0  (the lane's criterion (3))
    c_max=max(float(np.linalg.norm(A@R0-R0@A)) for A in basis)
    # (b) does ANY nonzero element of span{Phi(V)} ANTICOMMUTE with R  (the lane's d_R definition)
    rows=[]
    for A in basis:
        M=A@R0+R0@A
        rows.append(np.concatenate([M.real.ravel(),M.imag.ravel()]))
    Mx=np.array(rows)
    # also need the coefficient map: nonzero combination c with sum c_i (A_i R + R A_i) = 0 AND sum c_i A_i != 0
    s=np.linalg.svd(Mx,compute_uv=False)
    ns_dim=int(np.sum(s<1e-9*max(1.0,s.max())))
    # find such a combination if it exists
    _,_,Vt=np.linalg.svd(Mx)
    anti=None
    if ns_dim>0:
        for j in range(len(s)-ns_dim,len(basis)):
            c=Vt[j]
            A=sum(c[i]*basis[i] for i in range(len(basis)))
            if np.linalg.norm(A)>1e-8: anti=A; break
    D1=max(float(np.linalg.norm(A)) for A in basis)
    print(f"        D1 (max ||Phi(single-site V)||)                      = {D1:.4f}   -> T5 holds, first-order splitting")
    print(f"        max ||[Phi(V), R]|| over single-site V               = {c_max:.4f}")
    print(f"        exists single-site combination ANTICOMMUTING with R  = {anti is not None}"
          f"   -> d_R {'= 1' if anti is not None else '> 1'} by the lane's own definition")
    # cross-check with the honest question: does a WEIGHT-1 unitary flip R on E0?
    best=0.0
    rng=np.random.default_rng(3)
    for _ in range(400):
        k=rng.integers(0,L)
        A=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); A=herm(A)
        wv,uv=np.linalg.eigh(A); Uu=uv@np.diag(np.exp(1j*rng.normal()*wv))@uv.conj().T
        Us=site(Uu,k,L)
        Rb=Q.conj().T@Us.conj().T@(Q@R0@Q.conj().T)@Us@Q
        best=max(best,float(np.linalg.norm(Rb+R0)))
    print(f"        400 random single-site UNITARIES: max ||P U+ R U P + R|| = {best:.4f}"
          f"   (a flip would give 0)")
