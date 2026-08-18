"""F-10 v2.  The v1 discriminator was mis-designed: it measured d<Zbar>/dt under a
DAVIES bath, but F-9 already proved no detailed-balance bath forms a record at ANY
weight, so the answer was fixed by the bath class before weight could matter.

Correct structure:
  A. NECESSARY condition -- Knill-Laflamme: minimum weight acting non-trivially on the
     code space at all.  Bath-independent, so it cannot inherit F-9's exclusion.
  B. SUFFICIENCY -- an explicit NON-EQUILIBRIUM (non-detailed-balance) bath, and what
     its coupling must weigh.
  C. CONTROL -- can any weight-1 jump operator, Hermitian or not, do it?

PRE-REGISTERED PREDICTION: threshold == d.  At d=2, weight 1 gives exactly zero."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('f7_davies.py').read().split('say("="*104); say("0.')[0])

E,V=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E-E[0])<1e-9))
Pg=V[:,:gs]@V[:,:gs].conj().T
P1={'X':X,'Y':1j*(X@Z),'Z':Z}
def supp(M,tol=1e-9):
    """number of qubits M acts on non-trivially, by partial-trace test."""
    n=0
    for l in range(L):
        A=op({l:Z},L); B=op({l:X},L)
        if np.linalg.norm(M@A-A@M)>tol or np.linalg.norm(M@B-B@M)>tol: n+=1
    return n

say("="*100); say("F-10 v2   THE FORMATION THRESHOLD"); say("="*100)
say(f"  carrier 2x2 torus, dim {2**L}, ground degeneracy {gs}, code distance d = 2")
say(f"  logical weights: Zbar={LOGINFO[0]}, Xbar={LOGINFO[2]}")

say(""); say("A.  NECESSARY -- minimum weight acting NON-TRIVIALLY on the code space (Knill-Laflamme)")
say(f"    {'weight':>8}{'# couplings':>13}{'max ||P C P - (trC/4) P||':>28}{'verdict':>34}")
thr=None
for w in (1,2,3):
    mx=0.0; arg=None
    for sites in itertools.combinations(range(L),w):
        for letters in itertools.product('XYZ',repeat=w):
            C=op({s:P1[c] for s,c in zip(sites,letters)},L)
            M=Pg@C@Pg; c=np.trace(M)/gs; n=np.linalg.norm(M-c*Pg)
            if n>mx: mx=n; arg=(sites,''.join(letters))
    v='acts as a SCALAR -- cannot form' if mx<1e-10 else f'DISTINGUISHES code states ({arg[1]}@{arg[0]})'
    if mx>1e-10 and thr is None: thr=w
    ncoup=len(list(itertools.combinations(range(L),w)))*3**w
    say(f"    {w:>8}{ncoup:>13}{mx:>28.3e}{v:>34}")
say(f"    -> MINIMUM WEIGHT WITH ANY CODE-SPACE ACTION = {thr}    CODE DISTANCE d = 2    "
    f"{'PREDICTION CONFIRMED' if thr==2 else 'PREDICTION FALSIFIED'}")

say(""); say("B.  SUFFICIENCY -- an explicit NON-EQUILIBRIUM bath (single jump op, no KMS partner)")
Pp=(np.eye(2**L)+Zbar)/2; Pm=(np.eye(2**L)-Zbar)/2
sig=Pp@Xbar@Pm                                      # logical lowering |Zbar=+1><Zbar=-1|
say(f"    logical lowering operator sigma^- = P+ Xbar P- :  support = {supp(sig)} qubits, "
    f"||sigma^-|| = {np.linalg.norm(sig):.4f}")
say(f"    is it Hermitian? ||sig - sig^dag|| = {np.linalg.norm(sig-sig.conj().T):.4f}  "
    f"(a HERMITIAN jump operator can only dephase, never select)")
rho=Pg/gs
say(f"    {'t':>7}{'<Zbar>':>12}{'code weight':>14}{'purity on code':>17}")
def step(r,jumps,dt):
    d=np.zeros_like(r)
    for J in jumps:
        Jd=J.conj().T; d+=J@r@Jd-0.5*(Jd@J@r+r@Jd@J)
    return r+dt*d
J=[sig]; dt=0.02; t=0.0
say(f"    {t:>7.2f}{np.real(np.trace(rho@Zbar)):>12.6f}{np.real(np.trace(Pg@rho)):>14.6f}"
    f"{np.real(np.trace((Pg@rho@Pg)@(Pg@rho@Pg)))/max(np.real(np.trace(Pg@rho))**2,1e-30):>17.6f}")
for n in range(1,301):
    rho=step(rho,J,dt); t+=dt
    if n in (25,50,100,200,300):
        cw=np.real(np.trace(Pg@rho)); rc=Pg@rho@Pg/cw
        say(f"    {t:>7.2f}{np.real(np.trace(rho@Zbar)):>12.6f}{cw:>14.6f}"
            f"{np.real(np.trace(rc@rc)):>17.6f}")

say(""); say("C.  CONTROL -- can ANY weight-1 jump operator select?  (Hermitian and non-Hermitian)")
mx=0.0; argm=None
singles=[(l,c,op({l:P1[c]},L)) for l in range(L) for c in 'XYZ']
for (l1,c1,A) in singles:
    for (l2,c2,B) in singles:
        Lop=(A+1j*B)/2
        r=Pg/gs; Ld=Lop.conj().T
        d=Lop@r@Ld-0.5*(Ld@Lop@r+r@Ld@Lop)
        v=abs(np.real(np.trace(d@Zbar)))
        if v>mx: mx=v; argm=(c1,l1,c2,l2)
say(f"    swept {len(singles)**2} operators of the form (A + iB)/2 with A,B of weight 1")
say(f"    max |d<Zbar>/dt| = {mx:.3e}   at ({argm[0]}_{argm[1]} + i {argm[2]}_{argm[3]})/2")
say(f"    -> {'ZERO: no weight-1 jump operator can form a record' if mx<1e-10 else 'NONZERO -- PREDICTION FALSIFIED'}")
