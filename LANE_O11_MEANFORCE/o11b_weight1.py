"""O-11 follow-up. The weight-1 column in o11_meanforce was NOISE (~1e-15, sign oscillating),
so the fitted slope 0.2741 was a fit to zero and the auto-verdict misread it.

THE REAL QUESTION: is the weight-1 zero due to CODE DISTANCE (Knill-Laflamme), or merely
an accidental symmetry of the single operator picked? A coupling A gives <Zbar> = 0 by
symmetry alone whenever [A,Xbar] = 0. So the test must sweep ALL weight-1 couplings and
separate those that COMMUTE with the writer from those that do not."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
nS=2**L; E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
nB=4; EB=np.array([0.0,0.7,1.3,2.1]); bb=np.array([1.0,0.3,-0.2,-0.9]); HB=np.diag(EB); Bop=np.diag(bb)
def rho_MF(A,lam,beta=2.0):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),HB)+lam*np.kron(A,Bop)
    w,U=np.linalg.eigh(Ht); w=w-w.min(); M=(U*np.exp(-beta*w))@U.conj().T
    R=M.reshape(nS,nB,nS,nB).trace(axis1=1,axis2=3); return R/np.trace(R)
def zb(r):
    rc=Pg@r@Pg; return np.real(np.trace(rc@Zbar))/np.real(np.trace(rc))
P1={'X':X,'Y':1j*(X@Z),'Z':Z}

say("="*100); say("O-11b  IS THE WEIGHT-1 ZERO KNILL-LAFLAMME, OR JUST A SYMMETRY?"); say("="*100)
say(f"  NOISE FLOOR: <Zbar> for A = identity (no coupling at all) = {zb(rho_MF(np.eye(nS),0.8)):.3e}")
say("")
say(f"  ALL 24 WEIGHT-1 COUPLINGS at lambda = 0.8   (Xbar support = the writer)")
say(f"  {'group':<34}{'count':>7}{'max |<Zbar>|':>16}{'argmax':>12}")
comm=[]; anti=[]
for l in range(L):
    for c in 'XYZ':
        A=op({l:P1[c]},L); k=np.linalg.norm(A@Xbar-Xbar@A)
        v=abs(zb(rho_MF(A,0.8)))
        (comm if k<1e-9 else anti).append((f"{c}_{l}",k,v))
for nm,grp in (("weight-1, COMMUTES with Xbar",comm),("weight-1, does NOT commute with Xbar",anti)):
    if not grp: say(f"  {nm:<34}{0:>7}"); continue
    m=max(grp,key=lambda t:t[2])
    say(f"  {nm:<34}{len(grp):>7}{m[2]:>16.3e}{m[0]:>12}")
say("")
say("  WEIGHT-2 COMPARISON (d = 2) -- couplings that are logical operators")
for nm,A in (("Zbar  (weight 2 = d)",Zbar),("Zbar2 (weight 2 = d)",Zbar2)):
    k=np.linalg.norm(A@Xbar-Xbar@A); say(f"    {nm:<24} ||[A,Xbar]|| = {k:>8.3f}   |<Zbar>| = {abs(zb(rho_MF(A,0.8))):.6f}")
say("")
say("  DECIDING TEST -- a weight-1 operator that DOES break the writer symmetry:")
if anti:
    best=max(anti,key=lambda t:t[2])
    say(f"    strongest is {best[0]} with ||[A,Xbar]|| = {best[1]:.3f} and |<Zbar>| = {best[2]:.3e}")
    say(f"    -> {'STILL ZERO: the suppression is KNILL-LAFLAMME, not symmetry' if best[2]<1e-10 else 'NONZERO: weight-1 DOES bias, so clause (2) fails'}")
else:
    say("    NONE EXIST -- every weight-1 Pauli commutes with Xbar, so the zero is SYMMETRY, not distance.")
    say("    In that case clause (2) is NOT tested by this carrier and O-10's threshold argument is the only warrant.")
