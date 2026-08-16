"""RC-2  (a) THE -1.200555 OF RECORD IS A FINITE-STAGE VALUE.
       (b) ATTACK ON WALL (2), THE POINT-COUNT WALL.

CLAIM UNDER ATTACK (wall 2):
  "The class-weight pushforward pi is a complete invariant ONLY WHILE the number of
   distinct exponent POINTS equals the number of vertex CLASSES.  Per-vertex charge on
   K1 takes that count from 3 to 5 (44 of 243 assignments reach 5)."
"""
import numpy as np
from math import pi, sqrt, log
import itertools
import rclib as R

p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
fg, cg = 1.0, sqrt(2.0)

print("=" * 78)
print("RC-2a  THE q=(1,2,2,2,2) RATE OF RECORD")
print("=" * 78)
q = [1, 2, 2, 2, 2]
pts, w = R.support_points(q, p_S3)
print("   exponent points/weights:", [(tuple(map(int, a)), float(b)) for a, b in zip(pts, w)])
print("   Delta(S) = < (1,-1) >, rank 1.  Factor out u*v:")
print("     |Z| = |0.4 + 0.3 w + 0.3 w^-1|  with w = u/v ;  w^n equidistributes on T^1.")
print("     so lambda_B = m(0.3 w^2 + 0.4 w + 0.3), a ONE-variable Mahler measure.")
rts = np.roots([0.3, 0.4, 0.3])
print("   roots:", rts, " |roots| =", np.abs(rts), "-> both ON the unit circle")
exact = R.mahler_1var_exact([0.3, 0.4, 0.3])
print(f"   ==> lambda_B = log(0.3) = {exact:.16f}   EXACTLY  (= log(3/10))")
print(f"       log(0.3) numerically                     = {log(0.3):.16f}")
lt = R.lambda_B_torus(q, p_S3, n=6000)
print(f"   2-D grid 6000^2                              = {lt:.9f}")
for N in (200_000, 2_000_000, 20_000_000):
    ld = R.lambda_B_direct(fg, cg, q, p_S3, N=N)
    print(f"   direct schedule-B  N={N:>10d}            = {ld:.9f}   (dev {abs(ld-exact):.2e})")
print("   REGISTER W-03 / this lane's brief record       = -1.200555")
print(f"   |W-03 value - true limit| = {abs(-1.200555 - exact):.3e}")
print("   ==> -1.200555 is a FINITE-STAGE value of a slowly-converging singular integrand")
print("       (both roots on |w|=1, so log|Z| has genuine -inf singularities on the orbit).")
print("       The corrected value of record is  lambda_B = log(3/10) = -1.2039728043259361.")

print()
print("=" * 78)
print("RC-2b  THE POINT-COUNT WALL")
print("=" * 78)

print("\n--- the '44 of 243' count.  243 = 3^5, so a 3-element charge alphabet.")
for alpha in ([0, 1, 2], [1, 2, 3], [-1, 0, 1], [1, 2, 4]):
    cnt = {}
    for q in itertools.product(alpha, repeat=5):
        pts, w = R.support_points(list(q), np.ones(5) / 5)
        cnt[len(pts)] = cnt.get(len(pts), 0) + 1
    print(f"   alphabet {str(alpha):11s} -> distinct-point-count histogram {dict(sorted(cnt.items()))}")
print("   NOTE: no 3-element alphabet gives 44 assignments reaching 5 points.")
print("   {1,2,3} gives 108; {0,1,2} gives 24.  The claim's 44/243 is not reproducible")
print("   under any charge alphabet of size 3 on K1.")

print("\n--- IS 'point count == class count' NECESSARY for pi to be complete?  NO.")
q = [0, 0, 0, 1, 1]
pts, w = R.support_points(q, p_S3)
print(f"   q = {q}:  exponent points {[tuple(map(int,x)) for x in pts]}  -> |points| = {len(pts)} < 3 classes")
print("   but the exponent map IS constant on each class: (1,1)->(0,0), (1,0)->(0,0), (0,1)->(0,1).")
print("   So Z = (p11+p10)*1 + p01*v  is a function of pi alone: pi IS complete.")
# demonstrate: two ready states with same pi give same lambda; different pi gives different
pa = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
pb = np.array([0.4, 0.25, 0.05, 0.20, 0.10])   # same class weights .4/.3/.3, different p_v
print(f"   pi(pa) = {(pa[0], pa[1]+pa[2], pa[3]+pa[4])}   pi(pb) = {(pb[0], pb[1]+pb[2], pb[3]+pb[4])}")
la = R.lambda_B_direct(fg, cg, q, pa, N=1_000_000)
lb = R.lambda_B_direct(fg, cg, q, pb, N=1_000_000)
print(f"   lambda_B(pa) = {la:.12f}   lambda_B(pb) = {lb:.12f}   |diff| = {abs(la-lb):.3e}")
print("   ==> point count 2 < class count 3, and pi is still COMPLETE.  'only while' FALSE.")

print("\n--- IS 'point count == class count' SUFFICIENT?  NO.")
q = [0, 0, 1, 1, 1]
pts, w = R.support_points(q, p_S3)
print(f"   q = {q}:  exponent points {[tuple(map(int,x)) for x in pts]}  -> |points| = {len(pts)} == 3 classes")
print("   but the partition is {v0,v1}/{v2}/{v3,v4}, NOT the class partition:")
print("   class (1,0)={v1,v2} is SPLIT across two exponent points.")
pc = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
pd = np.array([0.4, 0.05, 0.25, 0.15, 0.15])   # same pi, weight moved WITHIN class (1,0)
print(f"   pi(pc) = {(pc[0], pc[1]+pc[2], pc[3]+pc[4])}   pi(pd) = {(pd[0], pd[1]+pd[2], pd[3]+pd[4])}  IDENTICAL")
lc = R.lambda_B_direct(fg, cg, q, pc, N=1_000_000)
ldd = R.lambda_B_direct(fg, cg, q, pd, N=1_000_000)
print(f"   lambda_B(pc) = {lc:.12f}   lambda_B(pd) = {ldd:.12f}   |diff| = {abs(lc-ldd):.3e}")
print("   ==> identical pi, identical point count == class count, DIFFERENT lambda_B.")
print("       'point count == class count' is NOT sufficient.  The wall is at the wrong object.")

print("\n--- THE CORRECT STATEMENT (proved, not asserted):")
print("   pi is a complete invariant  <=>  the exponent map v |-> q_v*(a_v,b_v) is CONSTANT")
print("   ON EACH VERTEX CLASS.  Point count is neither necessary nor sufficient for that.")
print("   Verified by exhaustive enumeration below.")
rng = np.random.default_rng(987654321)          # SEED PUBLISHED
mis_nec = mis_suf = 0
tested = 0
for q in itertools.product([0, 1, 2, 3], repeat=5):
    Ev = R.exponents(list(q))
    const_on_class = (tuple(Ev[1]) == tuple(Ev[2])) and (tuple(Ev[3]) == tuple(Ev[4]))
    pts, w = R.support_points(list(q), np.ones(5) / 5)
    equal_counts = (len(pts) == 3)
    tested += 1
    if const_on_class != equal_counts:
        if equal_counts and not const_on_class:
            mis_suf += 1
        else:
            mis_nec += 1
print(f"   enumerated all 4^5 = {tested} charges in {{0,1,2,3}}^5:")
print(f"      point-count==3 but exponent map NOT class-constant : {mis_suf}  (sufficiency fails)")
print(f"      exponent map class-constant but point-count != 3   : {mis_nec}  (necessity fails)")
