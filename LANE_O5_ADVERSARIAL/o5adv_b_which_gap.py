"""ADVERSARIAL 1b.  WHICH gap is 'Delta'?  Vary star and plaquette couplings SEPARATELY."""
import numpy as np
from o5_common import Zop, Xop, toric_H, local_perturbation, Z_A_SUP, X_A_SUP, STARS, PLAQS
V = local_perturbation(seed=2026); p0=1e-2; gt=4
Av = sum(Xop(s) for s in STARS); Bp = sum(Zop(pl) for pl in PLAQS)
Rt = Zop(Z_A_SUP)

def split_along(H, R0, g, p):
    """cluster width, and the splitting RESOLVED ALONG the record: the spread of PHP eigenvalues
    that distinguishes R0=+1 from R0=-1 sectors."""
    e,U = np.linalg.eigh(H+p*V); Uc=U[:,:g]
    Hc = np.diag(e[:g]); M = Uc.conj().T@R0@Uc; M=(M+M.conj().T)/2
    w,Q=np.linalg.eigh(M); Rc=Q@np.diag(np.sign(w))@Q.conj().T
    return e[g-1]-e[0], np.linalg.norm(Hc@Rc-Rc@Hc,2)

print("  VARY STAR COUPLING b:  H = -b*sum A_v - sum B_p  (star gap 4b, plaq gap 4)")
print(f"  {'b':>6s} {'star gap 4b':>12s} {'width':>16s} {'width*b':>16s} {'eps_dressed':>15s} {'eps*b':>15s}")
for b in (0.5,1.0,2.0,4.0,8.0,16.0):
    w_,ed = split_along(-b*Av-Bp, Rt, gt, p0)
    print(f"  {b:6.2f} {4*b:12.2f} {w_:16.8e} {w_*b:16.8e} {ed:15.6e} {ed*b:15.6e}")
print("\n  VARY PLAQUETTE COUPLING a:  H = -sum A_v - a*sum B_p  (star gap 4, plaq gap 4a)")
print(f"  {'a':>6s} {'plaq gap 4a':>12s} {'width':>16s} {'width*a':>16s} {'eps_dressed':>15s} {'eps*a':>15s}")
for a in (0.5,1.0,2.0,4.0,8.0,16.0):
    w_,ed = split_along(-Av-a*Bp, Rt, gt, p0)
    print(f"  {a:6.2f} {4*a:12.2f} {w_:16.8e} {w_*a:16.8e} {ed:15.6e} {ed*a:15.6e}")
print("""
  CONCLUSION.  The Z-type record Z_0 Z_1 is split by a virtual STAR-defect path, so the 'Delta' in
  delta = c p^d/Delta^(d-1) is THE STAR GAP, not 'the gap' of H.  The lane's lam-rescaling test
  moves BOTH gaps at once and therefore cannot identify which one the law refers to; it also cannot
  fail, being algebraically equivalent to the p-sweep it already fitted.""")
