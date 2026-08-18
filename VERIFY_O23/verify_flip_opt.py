"""Minimise ||P U^dag R U P + R|| over SINGLE-SITE unitaries (random search + local refinement).
POSITIVE CONTROL: the same optimiser on a record that IS single-site flippable."""
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
def u2(p):
    A=p[0]*SX+p[1]*SY+p[2]*SZ
    w,U=np.linalg.eigh(A)
    return U@np.diag(np.exp(1j*w))@U.conj().T
def anneal(f,rng,nstart=40,nstep=600):
    best=1e9;bp=None
    for _ in range(nstart):
        p=rng.normal(size=3)*2.0; v=f(p); step=1.0
        for t in range(nstep):
            q=p+rng.normal(size=3)*step
            vq=f(q)
            if vq<v: p,v=q,vq
            else: step*=0.995
        if v<best: best,bp=v,p
    return best
rng=np.random.default_rng(5)
for L in (3,5):
    H=heis(L,-1.0); ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-9*max(1.0,abs(ev).max()); Q=U[:,sel]
    Jz=sum(site(SZ,k,L) for k in range(L))/2
    Rg=np.diag(np.sign(np.real(np.diag(Jz)))).astype(complex)
    R0=Q.conj().T@Rg@Q
    best=min(anneal(lambda p,k=k: float(np.linalg.norm(
        Q.conj().T@(site(u2(p),k,L).conj().T@Rg@site(u2(p),k,L))@Q+R0)),rng) for k in range(L))
    print(f"L={L}: MIN over single-site unitaries of ||P U+ R U P + R|| = {best:.4f}   (||R0||={np.linalg.norm(R0):.3f})")
    Rc=site(SZ,0,L)
    bc=anneal(lambda p: float(np.linalg.norm(site(u2(p),0,L).conj().T@Rc@site(u2(p),0,L)+Rc)),rng)
    print(f"      POSITIVE CONTROL, R = Z_0 (X_0 flips it), SAME optimiser: {bc:.3e}")
    Ux=kron(*[np.array([[0,-1j],[-1j,0]],dtype=complex)]*L)
    print(f"      the GLOBAL transversal writer reaches: {np.linalg.norm(Q.conj().T@(Ux.conj().T@Rg@Ux)@Q+R0):.3e}")
