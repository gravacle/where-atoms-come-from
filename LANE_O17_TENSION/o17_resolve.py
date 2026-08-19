"""O-17: the record is destroyed at order 2 (O-16) and recorded at order 2n* (PF-5).

THE SUSPICION: the two arms used DIFFERENT COUPLING CLASSES.
  PF-5 coupled with sum_l Z_l -- Z-type, COMMUTING with the Z-type record label.
  O-16 jumped with single-site {X,Y,Z}, which includes operators ANTICOMMUTING with it.
Clause (ii) already rules on this: [L_k, R] = 0 for all k.

NOTE ON METHOD: the commutant of a SET equals the commutant of the algebra it generates, so
the algebra never has to be built. (record_model.py builds it and does not need to for this;
that is an efficiency defect in the model, recorded.)"""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
P1={'I':I2,'X':X,'Y':Y,'Z':Z}
def pauli(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,P1[c])
    return M
n=5; S=[pauli(g) for g in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
H0=-sum(S); Zb=pauli('ZZZZZ'); Xb=pauli('XXXXX'); N=2**n
def singles(letters):
    return [pauli(''.join(c if j==q else 'I' for j in range(n))) for q in range(n) for c in letters]
def _gram(gens):
    """M-dagger M for the stacked commutator map, without ever forming M.
       rank(M) = rank(M-dagger M), and the Gram matrix is N^2 x N^2 rather than
       (n_gens * N^2) x N^2 -- the direct SVD was the reason v2 timed out."""
    I=np.eye(N); G=np.zeros((N*N,N*N),dtype=complex)
    for g in gens:
        C=np.kron(I,g.T)-np.kron(g,I)
        G+=C.conj().T@C
    return G
def _null(gens, tol=1e-8):
    G=_gram(gens); w,V=np.linalg.eigh(G)
    keep=w < tol*max(1.0,abs(w).max())
    return V[:,keep]
def commutant_dim(gens):
    return int(_null(gens).shape[1])
def herm_dim(gens):
    """dimension of the HERMITIAN part of the commutant, as a REAL vector space."""
    Ns=_null(gens); B=[Ns[:,i].reshape(N,N) for i in range(Ns.shape[1])]
    cand=[]
    for A in B: cand += [(A+A.conj().T)/2, 1j*(A-A.conj().T)/2]
    if not cand: return 0
    R=np.vstack([np.concatenate([A.real.reshape(-1),A.imag.reshape(-1)]) for A in cand])
    return int(np.linalg.matrix_rank(R, tol=1e-7))

say("="*94); say("O-17  DID THE TWO ARMS USE DIFFERENT COUPLING CLASSES?"); say("="*94)
say(f"  carrier [[5,1,3]], dim {N}, code space 2, d = 3")
say(f"  SELF-CHECK ||[Zbar,H0]|| = {np.linalg.norm(Zb@H0-H0@Zb):.1e}  ||{{Zbar,Xbar}}|| = {np.linalg.norm(Zb@Xb+Xb@Zb):.1e}")
say("")
say("1.  DOES CLAUSE (ii) HOLD?   [L_k, R] = 0 for every Lindblad operator")
say(f"  {'environment':<34}{'max ||[L_k, Zbar]||':>22}{'clause (ii)':>14}")
for nm,letters in (("PF-5's set: single-site {Z}","Z"),("O-16's set: single-site {X,Y,Z}","XYZ")):
    m=max(np.linalg.norm(L@Zb-Zb@L) for L in singles(letters))
    say(f"  {nm:<34}{m:>22.4f}{('HOLDS' if m<1e-9 else 'VIOLATED'):>14}")
say("")
say("2.  DOES A RECORD EXIST AT ALL?   dim of the HERMITIAN commutant; 1 = scalars only = no record")
say(f"  {'environment':<34}{'dim A-prime':>14}{'dim Herm(A-prime)':>20}{'record possible?':>19}")
for nm,letters in (("no environment (bare code)",""),("single-site {Z}","Z"),("single-site {X,Y,Z}","XYZ")):
    g=[H0]+[L for L in singles(letters)] if letters else [H0]
    gens=g+[L.conj().T for L in g]
    cd=commutant_dim(gens); hd=herm_dim(gens)
    say(f"  {nm:<34}{cd:>14}{hd:>20}{('NO -- scalars only' if hd<=1 else 'yes'):>19}")
say("")
say("3.  READ")
say("    If {X,Y,Z} violates clause (ii) AND leaves only scalars in the commutant, then O-16")
say("    measured the decay of an object that is NOT A RECORD for that environment. O-17 is then")
say("    two arms run on two different (H,{L_k}) pairs, only one of which has a record at all.")
