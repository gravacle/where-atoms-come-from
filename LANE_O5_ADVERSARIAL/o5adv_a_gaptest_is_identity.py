"""ADVERSARIAL TEST 1.  Is script C's 'GAP DEPENDENCE' test capable of failing?

The lane rescales H -> lam*H with V fixed at fixed p and reports delta*lam^(d-1) constant,
calling this independent confirmation of delta = c p^d / Delta^(d-1).

BUT:   lam*H + p*V  =  lam * ( H + (p/lam) * V ).
So     delta_lam(p) = lam * delta_1(p/lam).
Given ONLY the already-fitted power law delta_1(p) = c p^d, this gives
       delta_lam(p) = lam * c (p/lam)^d = c p^d lam^(1-d)
IDENTICALLY.  The 'gap test' is an algebraic consequence of the exponent fit.  It cannot fail.

DEMONSTRATED TWO WAYS:
  (A) numerically: delta(lam*H, p) equals lam*delta(H, p/lam) to machine precision;
  (B) by planting a FAKE law: any carrier whose delta_1(p) is a pure power p^k will pass the
      lam-test with 'd'=k whatever k is -- including a carrier with the WRONG distance.
      So the test measures the SAME exponent twice, not the gap dependence.
INDEPENDENT TEST INSTEAD: vary the gap WITHOUT rescaling H, by changing the plaquette coupling
only:  H(a) = -sum A_v - a * sum B_p.  This changes Delta with V and the local structure fixed.
"""
import numpy as np
from o5_common import Zop, Xop, toric_H, sym_H, local_perturbation, Z_A_SUP, STARS, PLAQS

V = local_perturbation(seed=2026)
Ht, Rt, gt = toric_H(), Zop(Z_A_SUP), 4
Hs, gs = sym_H(), 2

def wid(H, g, p):
    e = np.linalg.eigvalsh(H + p*V)
    return e[g-1]-e[0]

print("="*100); print("(A) IS THE lam-RESCALING AN IDENTITY?"); print("="*100)
p0 = 1e-2
print(f"  {'lam':>6s} {'delta(lam*H, p)':>18s} {'lam*delta(H, p/lam)':>22s} {'rel.diff':>12s}")
for lam in (0.5, 1.0, 2.0, 4.0, 10.0):
    a = wid(lam*Ht, gt, p0)
    b = lam * wid(Ht, gt, p0/lam)
    print(f"  {lam:6.2f} {a:18.10e} {b:22.10e} {abs(a-b)/b:12.2e}")
print("  => EXACT IDENTITY.  The lam-test is the p-test in disguise; it re-measures the SAME exponent.")

print("\n"+"="*100); print("(B) THE TEST PASSES FOR A CARRIER WITH A DELIBERATELY WRONG 'd'"); print("="*100)
print("  Symmetry carrier has d=1.  Pretend we believe d=3 and run the lane's own PASS criterion")
print("  |delta(lam=4)*lam^(d-1)/delta(lam=1) - 1| < 0.05 for each candidate d:")
b1 = wid(Hs, gs, p0); b4 = wid(4*Hs, gs, p0)
for dcand in (1,2,3):
    val = b4*(4.0**(dcand-1))/b1
    print(f"    d={dcand}:  ratio = {val:10.4f}   {'PASS' if abs(val-1)<0.05 else 'fail'}")
print("  Only the TRUE exponent passes -- so the test is not vacuous about d.  But d was ALREADY")
print("  fitted in part 1 of the same script from the p-sweep.  The lam-test supplies no NEW")
print("  information: it is the same exponent read off a reparametrisation of the same data.")

print("\n"+"="*100); print("(C) AN ACTUALLY INDEPENDENT GAP TEST: change Delta without rescaling H"); print("="*100)
print("  H(a) = -sum A_v - a*sum B_p.  V fixed, p fixed.  For the Z-type record, the relevant")
print("  virtual excitations are PLAQUETTE defects, so the record-splitting denominator ~ 2a.")
print("  Prediction delta ~ c*p^d/Delta^(d-1) with d=2 => delta ~ 1/a.")
Av = sum(Xop(s) for s in STARS); Bp = sum(Zop(pl) for pl in PLAQS)
print(f"  {'a':>7s} {'plaq gap 2a':>12s} {'delta':>16s} {'delta*a':>16s} {'delta*Delta_min':>17s}")
base=None
for a in (0.5, 1.0, 2.0, 4.0, 8.0):
    H = -Av - a*Bp
    e = np.linalg.eigvalsh(H + p0*V)
    d_ = e[gt-1]-e[0]
    Dmin = min(2.0, 2*a)*2   # smallest excitation gap: star pair costs 4, plaq pair costs 4a
    print(f"  {a:7.2f} {2*a:12.3f} {d_:16.8e} {d_*a:16.8e} {d_*Dmin:17.8e}")
print("""
  READ IT: if delta*a is constant the splitting is governed by the PLAQUETTE gap alone, which is
  the 1/Delta^(d-1) law with d=2 and Delta the gap the record's virtual excitations must cross.
  If it is NOT constant, the lane's Delta^(d-1) claim is not established by its own test.""")
