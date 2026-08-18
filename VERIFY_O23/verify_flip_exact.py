"""EXACT version of PROBE 7's flip test.  For a single-site unitary U on site k,
P U^dag R U P lies in the span of {P s_a^(k) R s_b^(k) P : a,b in I,X,Y,Z}.
If -R is NOT in that span for any k, then NO single-site unitary flips R: d_R > 1 exactly."""
import numpy as np
I2=np.eye(2,dtype=complex); SX=np.array([[0,1],[1,0]],dtype=complex)
SY=np.array([[0,-1j],[1j,0]],dtype=complex); SZ=np.array([[1,0],[0,-1]],dtype=complex)
PS=[I2,SX,SY,SZ]
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
def inspan(target,ops,tol=1e-8):
    m=target.shape[0]
    M=np.array([np.concatenate([A.real.ravel(),A.imag.ravel()]) for A in ops]).T
    t=np.concatenate([target.real.ravel(),target.imag.ravel()])
    sol,res,rk,sv=np.linalg.lstsq(M,t,rcond=None)
    return float(np.linalg.norm(M@sol-t))
for L in (3,5):
    H=heis(L,-1.0); ev,U=np.linalg.eigh(H)
    sel=np.abs(ev-ev[0])<1e-9*max(1.0,abs(ev).max()); Q=U[:,sel]; n=Q.shape[1]
    Jz=sum(site(SZ,k,L) for k in range(L))/2
    w,Ub=np.linalg.eigh(herm(Q.conj().T@Jz@Q))
    R0=herm(Ub@np.diag(np.sign(w))@Ub.conj().T); Rfull=Q@R0@Q.conj().T
    print(f"L={L}, multiplet n={n}, R = sign(Jz_total) (global operator, [H,R]=0: "
          f"{np.linalg.norm(H@ (np.sign(np.round(2*np.diag(np.real(Jz)))) * np.eye(2**L)) - 0):.0e} placeholder)")
    Rg=np.diag(np.sign(np.real(np.diag(Jz)))).astype(complex)   # global R = sign(Jz_total)
    print(f"   global R=sign(Jz_tot):  ||R-R+||={np.linalg.norm(Rg-Rg.conj().T):.1e}  "
          f"||R^2-I||={np.linalg.norm(Rg@Rg-np.eye(2**L)):.1e}  ||[H,R]||={np.linalg.norm(H@Rg-Rg@H):.1e}")
    Ux=kron(*[np.array([[0,-1j],[-1j,0]],dtype=complex)]*L)
    print(f"   (iv) global transversal pi-rotation: ||U+RU + R|| (FULL space) = "
          f"{np.linalg.norm(Ux.conj().T@Rg@Ux+Rg):.2e}   -> WRITABLE")
    worst=1e9
    for k in range(L):
        ops=[Q.conj().T@(site(PS[a],k,L)@Rfull@site(PS[b],k,L))@Q for a in range(4) for b in range(4)]
        r=inspan(-R0,ops)
        worst=min(worst,r)
    print(f"   (v) EXACT: min over sites of dist(-R, span{{P s_a R s_b P}}) = {worst:.3e}"
          f"   -> {'NO single-site unitary can flip R' if worst>1e-6 else 'a single-site flip EXISTS'}")
    # POSITIVE CONTROL: a carrier where a single-site flip DOES exist -- R = Z on qubit 0, L qubits
    Rc=site(SZ,0,L); Q2=np.eye(2**L,dtype=complex)
    ops=[site(PS[a],0,L)@Rc@site(PS[b],0,L) for a in range(4) for b in range(4)]
    print(f"   POSITIVE CONTROL (R = Z_0, flippable by X_0): dist = {inspan(-Rc,ops):.3e}")
    # does a single-site perturbation FULLY RESOLVE the multiplet (clause (iii) death)?
    def phi(V):
        B=Q.conj().T@V@Q; m=B.shape[0]; return B-np.trace(B)/m*np.eye(m)
    sp=np.linalg.eigvalsh(herm(phi(site(SZ,0,L))))
    print(f"   Phi(Z on site 0) spectrum on E0 = {np.round(sp,4)}  min gap {np.min(np.diff(np.sort(sp))):.4f}"
          f"  -> {'FULLY RESOLVED at first order: clause (iii) dies for every R' if np.min(np.diff(np.sort(sp)))>1e-8 else 'not resolved'}")
    print()
