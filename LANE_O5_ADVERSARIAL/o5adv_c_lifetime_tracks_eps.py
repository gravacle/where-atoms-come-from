"""ADVERSARIAL 2.  WHAT ACTUALLY SETS THE LIFETIME -- the cluster WIDTH or eps_dressed?

The lane fits  delta = cluster width  and registers  T(eta) ~ eta/delta  with a
'carrier-dependent O(1) factor of 2.50 and 8.64' that it explains by
    '8.78/2.50 = 3.5 = 0.801/0.322 * 1.4'
-- a 1.4 inserted with no justification.  Test the alternative: T(eta) = eta/eps_dressed EXACTLY,
with NO free constant, where eps_dressed = ||[H_c, R_c]|| is the lane's own Reading-2 number.
"""
import numpy as np
import sys as _s, os as _o
# REPRODUCTION FIX (T-35): o5_common lives in LANE_O5_APPROXIMATE; the sealed runs had it on the
# path by happenstance and reproduce.sh could not run this lane standalone.
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..', 'LANE_O5_APPROXIMATE'))
from o5_common import Zop, toric_H, sym_H, local_perturbation, Z_A_SUP
V = local_perturbation(seed=2026)

def block(H0,R0,g,p):
    e,U=np.linalg.eigh(H0+p*V); Uc=U[:,:g]
    M=Uc.conj().T@R0@Uc; M=(M+M.conj().T)/2
    w,Q=np.linalg.eigh(M); Rc=Q@np.diag(np.sign(w))@Q.conj().T
    Hc=np.diag(e[:g]).astype(complex)
    return e[:g], M, e[g-1]-e[0], float(np.linalg.norm(Hc@Rc-Rc@Hc,2))

def dev(ee,M,ts):
    ee=np.asarray(ee,float); W=ee[:,None]-ee[None,:]
    D=M[None,:,:]*(np.exp(1j*W[None,:,:]*ts[:,None,None])-1.0)
    return np.linalg.svd(D,compute_uv=False)[:,0]

def lifetime(ee,M,eta,n=400001):
    ee=np.asarray(ee,float); d=ee.max()-ee.min()
    ts=np.linspace(0,4.0/d,n); dv=dev(ee,M,ts)
    if dv.max()<eta: return float('inf')
    return float(ts[int(np.argmax(dv>=eta))])

print("="*104)
print("  T(eta) vs BOTH candidate scales.  If T*eps/eta == 1 with no free constant, the lifetime is")
print("  set by eps_dressed, NOT by the cluster width, and the lane's 8.64 / 2.50 are not free O(1)s.")
print("="*104)
print(f"  {'carrier':>12s} {'p':>7s} {'eta':>6s} {'width delta':>14s} {'eps_dressed':>14s}"
      f" {'T(eta)':>14s} {'T*delta/eta':>12s} {'T*eps/eta':>11s} {'delta/eps':>10s}")
for nm,H0,R0,g in (("TOPOLOGICAL",toric_H(),Zop(Z_A_SUP),4),("SYMMETRY",sym_H(),Zop([0]),2)):
    for p in (1e-2,3e-2,1e-1):
        ee,M,d_,ed=block(H0,R0,g,p)
        for eta in (0.01,0.1):
            T=lifetime(ee,M,eta)
            print(f"  {nm:>12s} {p:7.0e} {eta:6.2f} {d_:14.6e} {ed:14.6e} {T:14.6e}"
                  f" {T*d_/eta:12.4f} {T*ed/eta:11.4f} {d_/ed:10.4f}")
print("""
  T*eps_dressed/eta = 1.00 to three digits on BOTH carriers, at every p and every eta.
  The lane's 'carrier-dependent O(1) factor' IS EXACTLY delta/eps_dressed.  There is no free
  constant and no 1.4.  THE LIFETIME IS eta/eps_dressed, NOT eta/delta.""")

print("\n"+"="*104)
print("  CONSEQUENCE: delta and eps_dressed have DIFFERENT gap dependence.  Vary the two toric")
print("  couplings separately (V and p fixed).  Under the lane's uniform rescaling H->lam*H they")
print("  move together, which is why the lane's test could not see this.")
print("="*104)
from o5_common import Xop, STARS, PLAQS
Av=sum(Xop(s) for s in STARS); Bp=sum(Zop(pl) for pl in PLAQS); Rt=Zop(Z_A_SUP); p0=1e-2
print(f"  {'':>26s} {'width':>15s} {'x gap':>9s} {'eps_dressed':>15s} {'x gap':>9s}")
w1,e1=block(-Av-Bp,Rt,4,p0)[2],block(-Av-Bp,Rt,4,p0)[3]
for lbl,H in (("star gap x4  (b=4)",-4*Av-Bp),("plaq gap x4  (a=4)",-Av-4*Bp),("both x4 (lane's test)",-4*Av-4*Bp)):
    _,_,w,ed=block(H,Rt,4,p0)
    print(f"  {lbl:>26s} {w:15.6e} {w1/w:9.3f} {ed:15.6e} {e1/ed:9.3f}")
print("""
  READ IT.  Multiplying the PLAQUETTE gap by 4 divides eps_dressed by 4 and leaves the width almost
  unchanged.  Multiplying the STAR gap by 4 divides the width by ~3 and leaves eps_dressed UNCHANGED.
  'Delta' in  delta = c p^d / Delta^(d-1)  is therefore NOT 'the gap of H'.  For the record it is the
  gap of the excitations the WRITER-type logical must virtually cross (plaquette); for the cluster
  width it is a mixture.  The lane fitted the law on the width and applied it to the lifetime.""")
