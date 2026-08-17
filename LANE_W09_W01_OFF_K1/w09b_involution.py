# W-09 leg B — NAME THE OPERATIVE VARIABLE CORRECTLY.
# Leg A shows B1q (spectator, NO vertex in both loops) gives numbers IDENTICAL to K1 (pinch, NO
# vertex in neither). So "three versus four coefficients" is not the name, and "p00 = 0" is not
# the name either -- B1q has p00 > 0 and behaves exactly like K1.
# HYPOTHESIS: they are identical because they are exchanged by the involution the register ALREADY
# CARRIES at W-03 -- "PINCH AND SPECTATOR ARE THE SAME OBJECT, exactly: 00 <-> 11 and 10 <-> 01."
import numpy as np
rng = np.random.default_rng(20260816); N = 200000
INV = {"00":"11", "11":"00", "10":"01", "01":"10"}

def chars(occ, f, c):
    m = {"00": np.ones_like(f, dtype=complex), "10": np.exp(-1j*f),
         "01": np.exp(1j*c),                   "11": np.exp(1j*(c-f))}
    return [m[k] for k in occ]
def zero_in_hull(pts):
    A = np.sort(np.angle(np.stack(pts, axis=0)), axis=0)
    g = np.diff(np.concatenate([A, A[:1] + 2*np.pi], axis=0), axis=0)
    return g.max(axis=0) <= np.pi + 1e-12

f = rng.uniform(-np.pi, np.pi, N); c = rng.uniform(-np.pi, np.pi, N)
K1  = ["01","10","11"]        # pinch only     (v0 in BOTH loops; no vertex in neither)
B1q = ["00","01","10"]        # spectator only (spectator in NEITHER; no vertex in both)
print("== B1  IS B1q's CLASS SET THE INVOLUTION IMAGE OF K1's? ==")
print(f"  K1  occupies {K1}")
print(f"  its involution image {sorted(INV[k] for k in K1)}")
print(f"  B1q occupies {B1q}")
print(f"  EQUAL: {sorted(INV[k] for k in K1) == B1q}")
print("  -> K1 and B1q are exchanged by W-03's registered involution. That, and not the")
print("     coefficient count, is why leg A's two rows are identical to the last digit.\n")

print("== B2  THE OPERATIVE VARIABLE, ISOLATED ==")
print("   One variable: whether the incidence occupies ALL FOUR classes. Everything else fixed.")
print(f"   {'occupied classes':<22} {'#':>2} {'firing region':>14} {'f->-f flips':>12} {'curvature-aware?':>17}")
for occ in [["01","10"],["00","01"],["01","10","11"],["00","01","10"],["00","10","11"],
            ["00","01","11"],["00","01","10","11"]]:
    fire = zero_in_hull(chars(occ, f, c)); flip = int((fire != zero_in_hull(chars(occ,-f,c))).sum())
    print(f"   {str(occ):<22} {len(occ):>2} {fire.mean():>14.6f} {flip:>12} "
          f"{'YES' if flip>0 else 'NO':>17}")
print()
print("   EVERY three-class carrier: firing region 1/4, curvature-aware.")
print("   The ONLY four-class carrier: region 1/2, curvature-BLIND (0 of 200000 flips).")
print("   FOUR CLASSES REQUIRES a vertex in BOTH loops AND a vertex in NEITHER.")
print("   K1 has the first and not the second. B1q has the second and not the first.")
print()
print("   ==> THE OPERATIVE VARIABLE IS NOT 'three versus four coefficients' AND NOT 'p00 = 0'.")
print("       IT IS WHETHER THE INCIDENCE OCCUPIES ALL FOUR CLASSES -- pinch AND spectator")
print("       together. Either alone reproduces K1 exactly, by an involution already of record.")
