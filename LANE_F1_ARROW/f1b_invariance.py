"""F-1 part 3, corrected. The v1 control conjugated the observable the WRONG WAY
(Q-dag Zbar Q where covariance needs Q Zbar Q-dag), so its 'falsification' measured a
rotated observable, not a change in chi.

Two distinct questions, both answered here:
  (a) COVARIANCE -- an instrument check: chi(O in r) must equal chi(QOQ-dag in QrQ-dag).
  (b) THE ARROW -- can a system-only unitary DESTROY the correlation? The invariant is the
      full mutual information I(S:B). Local unitaries preserve it exactly. The system can
      MOVE which observable carries the bath's knowledge; it cannot remove it."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
nS=2**L; nB=4
EB=np.array([0.0,0.7,1.3,2.1]); bb=np.array([1.0,0.3,-0.2,-0.9]); HB=np.diag(EB); Bop=np.diag(bb); beta=2.0
def rho_SB(A,lam):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),HB)+lam*np.kron(A,Bop)
    w,U=np.linalg.eigh(Ht); w=w-w.min(); M=(U*np.exp(-beta*w))@U.conj().T
    return M/np.trace(M)
def vN(r):
    e=np.linalg.eigvalsh(r); e=e[e>1e-13]; return float(-(e*np.log2(e)).sum())
def chi(r,O):
    I2=np.eye(nS); out=[]
    for s in (+1,-1):
        P=np.kron((I2+s*O)/2,np.eye(nB)); blk=P@r@P; p=np.real(np.trace(blk))
        if p<1e-12: out.append((0.0,np.zeros((nB,nB)))); continue
        out.append((p,(blk/p).reshape(nS,nB,nS,nB).trace(axis1=0,axis2=2)))
    av=sum(p*rB for p,rB in out)
    return max(vN(av)-sum(p*vN(rB) for p,rB in out),0.0)
def mutual(r):
    rS=r.reshape(nS,nB,nS,nB).trace(axis1=1,axis2=3); rB=r.reshape(nS,nB,nS,nB).trace(axis1=0,axis2=2)
    return vN(rS)+vN(rB)-vN(r)

r=rho_SB(Zbar,0.8); c0=chi(r,Zbar); I0=mutual(r)
say("="*98); say("F-1 part 3, CORRECTED"); say("="*98)
say(f"  chi(Zbar:B) = {c0:.8f} bits     I(S:B) = {I0:.8f} bits")
say("")
say("(a) COVARIANCE -- instrument check, chi(O in r) vs chi(Q O Q-dag in Q r Q-dag)")
rng=np.random.default_rng(5); wc=0.0; wi=0.0; wfixed=0.0
for _ in range(12):
    M=rng.normal(size=(nS,nS))+1j*rng.normal(size=(nS,nS)); Q,_=np.linalg.qr(M)
    Us=np.kron(Q,np.eye(nB)); rp=Us@r@Us.conj().T
    wc=max(wc,abs(chi(rp,Q@Zbar@Q.conj().T)-c0))          # covariant: must be 0
    wi=max(wi,abs(mutual(rp)-I0))                          # I(S:B): must be 0
    wfixed=max(wfixed,abs(chi(rp,Zbar)-c0))                # FIXED observable: may change
say(f"    max |chi(Q O Q-dag in Q r Q-dag) - chi(O in r)| = {wc:.3e}   "
    f"{'PASS -- covariant, so v1 measured a rotated observable' if wc<1e-8 else 'FAIL'}")
say("")
say("(b) THE ARROW -- what a system-only unitary can and cannot do")
say(f"    max |I(S:B) change| over 12 random system-only unitaries = {wi:.3e}   "
    f"{'EXACTLY INVARIANT' if wi<1e-8 else 'CHANGES'}")
say(f"    max |chi about the FIXED label Zbar| change              = {wfixed:.3e}   (may move -- and does)")
say("")
say("    READ: a system-only operation can MOVE which observable the bath knows about.")
say("          It cannot change I(S:B) at all. The correlation is not erasable from inside.")
say("")
say("(c) HOW MUCH IS UNREACHABLE, by coupling weight")
Ze=op({ind[('h',0,0)]:Z},L)
say(f"    {'coupling':<24}{'weight':>8}{'I(S:B) bits':>15}{'chi(Zbar:B) bits':>19}")
for nm,A,w in (("Zbar  (logical)",Zbar,2),("Zbar2 (logical)",Zbar2,2),("Z_e   (single site)",Ze,1),("identity",np.eye(nS),0)):
    rr=rho_SB(A,0.8); say(f"    {nm:<24}{w:>8}{mutual(rr):>15.8f}{chi(rr,Zbar):>19.8f}")
