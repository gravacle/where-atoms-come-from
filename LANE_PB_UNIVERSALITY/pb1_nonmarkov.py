"""PHASE B / B1.  DOES T1 SURVIVE A NON-MARKOVIAN ENVIRONMENT?

T1 was proved for the GKSL generator, which assumes a Markovian, weak-coupling, memoryless bath.
Every result in this program inherits that assumption. It is the single most suspect hypothesis in
the whole account, because einselection stories most often fail exactly there.

THE ARGUMENT THAT IT SURVIVES, and it is simpler than the GKSL proof. Model system AND environment
explicitly and unitarily, with
        H_tot = H_sys (x) I + I (x) H_env + H_int
If R acts on the system alone, [H_sys,R]=0, and [H_int, R (x) I]=0, then R (x) I commutes with the
FULL H_tot -- because H_env acts on the other factor. Then <R> is EXACTLY conserved for all time,
with no Markov approximation, no weak coupling, no memorylessness, and no bath spectral assumption.
GKSL is not needed at all: T1 is a CONSERVATION statement.

TEST. An explicit finite environment (hence genuinely non-Markovian: it recurs, it back-acts, it has
memory), strong coupling, long times. Then break each hypothesis and confirm it breaks.
"""
import numpy as np
rng=np.random.default_rng(5)
def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X
def U(D):
    Q,_=np.linalg.qr(rng.normal(size=(D,D))+1j*rng.normal(size=(D,D))); return Q

def run(DS,DE,break_what=None,T=40.0,coupling=3.0):
    Q=U(DS)
    r=np.exp(1j*rng.uniform(0,2*np.pi,DS))
    R=Q@np.diag(r)@Q.conj().T                                  # unitary system observable
    Hs=Q@np.diag(rng.normal(size=DS))@Q.conj().T; Hs=(Hs+Hs.conj().T)/2   # [Hs,R]=0 exactly
    He=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(DE,DE))+1j*rng.normal(size=(DE,DE)))
    # interaction built to commute with R (x) I: system factor diagonal in R's eigenbasis
    Hint=np.zeros((DS*DE,DS*DE),complex)
    for _ in range(3):
        A=Q@np.diag(rng.normal(size=DS))@Q.conj().T
        B=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(DE,DE))+1j*rng.normal(size=(DE,DE)))
        Hint=Hint+coupling*np.kron((A+A.conj().T)/2,B)
    if break_what=='Hs':   Hs=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(DS,DS))+1j*rng.normal(size=(DS,DS)))
    if break_what=='Hint':
        Hint=np.zeros((DS*DE,DS*DE),complex)
        for _ in range(3):
            A=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(DS,DS))+1j*rng.normal(size=(DS,DS)))
            B=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(DE,DE))+1j*rng.normal(size=(DE,DE)))
            Hint=Hint+coupling*np.kron(A,B)
    Htot=np.kron(Hs,np.eye(DE))+np.kron(np.eye(DS),He)+Hint
    RI=np.kron(R,np.eye(DE))
    com=np.linalg.norm(Htot@RI-RI@Htot)
    A=rng.normal(size=(DS*DE,DS*DE))+1j*rng.normal(size=(DS*DE,DS*DE))
    rho=A@A.conj().T; rho=rho/np.trace(rho).real
    vals=[]
    for t in (0.0,1.0,5.0,15.0,T):
        Ut=expm(-1j*Htot*t)
        vals.append(np.trace(RI@(Ut@rho@Ut.conj().T)).real)
    return com,vals

print("PHASE B / B1.  T1 UNDER AN EXPLICIT FINITE ENVIRONMENT -- no Markov approximation anywhere.")
print("  The environment is finite, so it RECURS and has MEMORY. Coupling is strong (3.0).")
print(f"\n  {'case':>34s} {'||[H_tot, R(x)I]||':>19s} {'<R> drift over t=0..40':>24s}")
print("  "+"-"*82)
for DS,DE in ((3,4),(4,6),(5,8)):
    com,v=run(DS,DE)
    print(f"  {f'DS={DS} DE={DE}  all hypotheses hold':>34s} {com:19.3e} {max(v)-min(v):24.3e}")
print()
print("  DROP EACH HYPOTHESIS -- each must break conservation.")
for bw,name in ((None,'all hold'),('Hs','[H_sys,R] != 0'),('Hint','[H_int, R(x)I] != 0')):
    com,v=run(4,6,break_what=bw)
    print(f"    {name:>24s}   ||[H_tot,R]|| = {com:9.3e}   drift = {max(v)-min(v):9.3e}")
print()
print("  READING: T1 needs no Markov approximation, no weak coupling and no memorylessness.")
print("  It is a CONSERVATION statement, and GKSL was never required to obtain it.")
