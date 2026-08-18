"""O-11: does a STRONG-COUPLING THERMAL environment form a record?

F-13 clause (1) claims the RFP needs dynamics OUTSIDE detailed balance. But a strongly
coupled thermal bath relaxes to the MEAN-FORCE Gibbs state rho_MF = Tr_B[e^{-beta H_tot}]/Z,
which depends on the COUPLING operators, not on H_S alone.

PRE-REGISTERED (REGISTER_V001, commit b2272eb):
  (1) lambda -> 0 gives <Zbar> -> 0 for every coupling  [the control]
  (2) weight-d coupling: <Zbar> != 0, rising at order lambda^1, matching the closed form
        <Zbar> = [Z_B(+1) - Z_B(-1)] / [Z_B(+1) + Z_B(-1)],  Z_B(s) = sum_j e^{-beta(E_j + s lambda b_j)}
  (3) weight-1 coupling: bias SUPPRESSED as lambda^d
Carrier: toric code, 2x2 torus, L=8 qubits, dim 256, d=2. Bath: 4 levels."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])

nS=2**L
E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
say("="*104); say("O-11  MEAN-FORCE GIBBS STATE ON THE CODE SPACE"); say("="*104)
say(f"  system: toric code 2x2, dim {nS}, ground degeneracy {gs}, d = 2")
say(f"  SELF-CHECK  ||[Zbar,H_S]|| = {np.linalg.norm(Zbar@H0-H0@Zbar):.3e}  "
    f"||{{Zbar,Xbar}}|| = {np.linalg.norm(Zbar@Xbar+Xbar@Zbar):.3e}   "
    f"{'PASS' if max(np.linalg.norm(Zbar@H0-H0@Zbar),np.linalg.norm(Zbar@Xbar+Xbar@Zbar))<1e-9 else 'FAIL'}")

# --- bath: 4 levels, DIAGONAL, deliberately ASYMMETRIC so B -> -B is not a symmetry
nB=4
EB=np.array([0.0,0.7,1.3,2.1]); bb=np.array([1.0,0.3,-0.2,-0.9])
HB=np.diag(EB); Bop=np.diag(bb)
say(f"  bath: {nB} levels, H_B = {EB}, B = {bb}   (sum b_j = {bb.sum():+.3f}, so B -> -B is NOT a symmetry)")

def rho_MF(A,lam,beta):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),HB)+lam*np.kron(A,Bop)
    w,U=np.linalg.eigh(Ht); w=w-w.min()
    M=(U*np.exp(-beta*w))@U.conj().T
    R=M.reshape(nS,nB,nS,nB).trace(axis1=1,axis2=3)
    return R/np.trace(R)

def zbar_code(r):
    rc=Pg@r@Pg; t=np.real(np.trace(rc))
    return np.real(np.trace(rc@Zbar))/t, t

beta=2.0
Ze=op({ind[('h',0,0)]:Z},L)                    # weight-1
say(""); say(f"  {'lambda':>9}{'<Zbar> A=Zbar (wt 2=d)':>26}{'closed form':>15}{'<Zbar> A=Z_e (wt 1)':>22}{'code weight':>13}")
rows=[]
for lam in (0.0,0.02,0.05,0.1,0.2,0.4,0.8):
    zd,cw=zbar_code(rho_MF(Zbar,lam,beta))
    ZBp=np.exp(-beta*(EB+lam*bb)).sum(); ZBm=np.exp(-beta*(EB-lam*bb)).sum()
    cf=(ZBp-ZBm)/(ZBp+ZBm)
    z1,_=zbar_code(rho_MF(Ze,lam,beta))
    rows.append((lam,zd,cf,z1))
    say(f"  {lam:>9.3f}{zd:>26.8f}{cf:>15.8f}{z1:>22.3e}{cw:>13.6f}")

say("")
mx=max(abs(r[1]-r[2]) for r in rows)
say(f"  CLOSED FORM AGREEMENT  max |numeric - analytic| = {mx:.3e}   {'PASS' if mx<1e-9 else 'FAIL'}")
say(f"  CONTROL lambda=0        <Zbar> = {rows[0][1]:.3e} (wt d), {rows[0][3]:.3e} (wt 1)   "
    f"{'PASS -- weak coupling recovers ordinary Gibbs' if max(abs(rows[0][1]),abs(rows[0][3]))<1e-12 else 'FAIL'}")

sm=[r for r in rows if 0<r[0]<=0.1]
sd=np.polyfit(np.log([r[0] for r in sm]),np.log([abs(r[1]) for r in sm]),1)[0]
s1=np.polyfit(np.log([r[0] for r in sm]),np.log([abs(r[3]) for r in sm]),1)[0]
say(f"  log-log slope in lambda   weight-d coupling : {sd:.4f}   (predicted 1)")
say(f"  log-log slope in lambda   weight-1 coupling : {s1:.4f}   (predicted d = 2)")

say(""); say("  SYMMETRY CONTROL -- a weight-d coupling that COMMUTES with the writer Xbar must give 0")
for nm,A in (("Xbar (commutes with Xbar)",Xbar),("Zbar (anticommutes)",Zbar)):
    c=np.linalg.norm(A@Xbar-Xbar@A); z,_=zbar_code(rho_MF(A,0.4,beta))
    say(f"    A = {nm:<28} ||[A,Xbar]|| = {c:>8.3f}   <Zbar> = {z:+.3e}   "
        f"{'ZERO as predicted' if abs(z)<1e-10 else 'NONZERO'}")

say(""); say("="*104)
big=abs(rows[-1][1])
say(f"  VERDICT: a GENUINELY THERMAL state of H_tot -- detailed balance w.r.t. the TOTAL --")
say(f"           gives <Zbar> = {rows[-1][1]:+.6f} at lambda = {rows[-1][0]} with a weight-d coupling.")
say(f"           F-13 clause (1) is {'CONTRADICTED' if big>1e-6 else 'UPHELD'}.")
say(f"           F-13 clause (2) {'SURVIVES in scaling form' if s1>sd+0.5 else 'does NOT survive'}: "
    f"weight-d acts at order {sd:.2f}, weight-1 at order {s1:.2f}.")
say("="*104)
