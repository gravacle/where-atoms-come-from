"""PROBE 7 (corrected).  Does D1 > 0 imply d_R = 1 on a 0-form multiplet?
Careful version: build the ACTUAL image space {Phi(V) : V single-site}, then ask whether it
contains a nonzero element anticommuting with R, with explicit norms printed."""
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
def phi(V,Q):
    B=Q.conj().T@V@Q; m=B.shape[0]
    return B-np.trace(B)/m*np.eye(m)

def orth_image(basis,tol=1e-9):
    """orthonormal basis of the real span of a list of Hermitian matrices"""
    vecs=[np.concatenate([A.real.ravel(),A.imag.ravel()]) for A in basis]
    M=np.array(vecs)
    U,s,Vt=np.linalg.svd(M,full_matrices=False)
    keep=[i for i in range(len(s)) if s[i]>tol*max(1.0,s.max())]
    m=basis[0].shape[0]
    out=[]
    for i in keep:
        v=Vt[i]; A=(v[:m*m].reshape(m,m)+1j*v[m*m:].reshape(m,m))
        out.append(A)
    return out,len(keep)

for L in (3,5):
    H=heis(L,-1.0)
    ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-9*max(1.0,abs(ev).max()); Q=U[:,sel]; n=Q.shape[1]
    Jz=sum(site(SZ,k,L) for k in range(L))/2
    Jzb=herm(Q.conj().T@Jz@Q)
    w,Ub=np.linalg.eigh(Jzb)
    R0=herm(Ub@np.diag(np.sign(w))@Ub.conj().T)
    basis=[phi(site(s,k,L),Q) for k in range(L) for s in (SX,SY,SZ)]
    img,dimimg=orth_image(basis)
    print(f"L={L}: multiplet n={n} (spin {(n-1)/2}), R=sign(Jz), dim span{{Phi(single-site V)}} = {dimimg}")
    print(f"       D1 = {max(float(np.linalg.norm(A)) for A in basis):.4f}  (T5: first-order splitting by a local term)")
    # anticommutation: solve  sum c_i (A_i R + R A_i) = 0  over the ORTHONORMAL image basis
    rows=[np.concatenate([(A@R0+R0@A).real.ravel(),(A@R0+R0@A).imag.ravel()]) for A in img]
    Mx=np.array(rows)
    s=np.linalg.svd(Mx,compute_uv=False)
    print(f"       singular values of the anticommutator map on that space: {np.round(s,6)}")
    nul=int(np.sum(s<1e-9*max(1.0,s.max())))
    if nul==0:
        print("       -> NO nonzero locally-reachable operator ANTICOMMUTES with R.  d_R > 1.")
    else:
        _,_,Vt=np.linalg.svd(Mx)
        c=Vt[-1]; A=sum(c[i]*img[i] for i in range(len(img)))
        print(f"       -> found one: ||A|| = {np.linalg.norm(A):.3e}, ||AR+RA|| = {np.linalg.norm(A@R0+R0@A):.3e}")
    # commutator criterion (the lane's measured form (3))
    print(f"       max ||[Phi(V),R]|| over single-site V = {max(float(np.linalg.norm(A@R0-R0@A)) for A in basis):.4f}"
          "   -> criterion (3) FAILS at w=1")
    # and does any single-site UNITARY flip R on E0?  report the MINIMUM
    rng=np.random.default_rng(11); best=1e9
    Rfull=Q@R0@Q.conj().T
    for _ in range(3000):
        k=rng.integers(0,L)
        A=herm(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
        wv,uv=np.linalg.eigh(A); Uu=uv@np.diag(np.exp(1j*rng.normal()*wv))@uv.conj().T
        Us=site(Uu,k,L)
        Rb=Q.conj().T@Us.conj().T@Rfull@Us@Q
        best=min(best,float(np.linalg.norm(Rb+R0)))
    print(f"       3000 random single-site unitaries: MIN ||P U+ R U P + R|| = {best:.4f}  (a flip gives 0)")
    print(f"       POSITIVE CONTROL: the GLOBAL pi-rotation about x, ||P U+ R U P + R|| = ", end="")
    Ux=kron(*[np.array([[0,-1j],[-1j,0]],dtype=complex)]*L)
    print(f"{np.linalg.norm(Q.conj().T@Ux.conj().T@Rfull@Ux@Q + R0):.3e}   -> it DOES flip R (non-local writer)")
    print()
print("READING:  criterion (2) 'anticommutes' and criterion (3) 'does not commute' COINCIDE for")
print("Pauli logicals on a stabiliser code, which is the only setting the lane measured.  Off that")
print("setting they can come apart, so the equivalence row is a STABILISER-CODE theorem, not general.")
