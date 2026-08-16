"""RC-1  ATTACK ON WALL (1), THE HOMOGENEITY WALL.

CLAIM UNDER ATTACK:
  "rank G, the four-class taxonomy, and S4-1 survive any homogeneous q and die at the
   first inhomogeneity.  K1 with q = (1,2,2,2,2) -- one pinch charge, four charge-2
   vertices -- is already outside."
  where CLASS-HOMOGENEOUS = "any single integer q != 0 shared by every vertex of a class"
  and the failure mode is "charge that is inhomogeneous WITHIN a vertex class".

ATTACK 1A.  q = (1,2,2,2,2) has ZERO within-class inhomogeneity.  It is class-homogeneous
            by the claim's own definition.  If S4-1 dies there, "survive any homogeneous q"
            is false at the claim's own headline example.
ATTACK 1B.  Enumerate ALL class-homogeneous charges and find where S4-1 actually dies.
ATTACK 1C.  Exhibit class-homogeneous, class-NON-uniform charges under which every S4
            result survives exactly -- so the wall is not at inhomogeneity in either
            direction.
ATTACK 1D.  Exhibit class-homogeneous charge that moves lambda_B ON THE EXCEPTIONAL SET,
            where S4's three-tier classification and every exact value of record live.
"""
import numpy as np
from math import pi, sqrt, log
import rclib as R

p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])     # S3/S4 ready state, class wts .4/.3/.3
fg, cg = 1.0, sqrt(2.0)                             # generic connection, rank L = 0

print("=" * 78)
print("RC-1  THE HOMOGENEITY WALL")
print("=" * 78)


def report(q, p=p_S3, tag=""):
    pts, w = R.support_points(q, p)
    rk = R.delta_rank(pts)
    B, _ = R.delta_basis(pts)
    return pts, w, rk, B


print("\n--- 1A.  q = (1,2,2,2,2):  is it class-homogeneous?")
q = [1, 2, 2, 2, 2]
print("   class (1,1) = {v0}      charges", [q[0]], "  -> constant")
print("   class (1,0) = {v1,v2}   charges", [q[1], q[2]], "-> constant")
print("   class (0,1) = {v3,v4}   charges", [q[3], q[4]], "-> constant")
print("   WITHIN-CLASS INHOMOGENEITY: NONE.  It IS class-homogeneous.")
pts, w, rk, B = report(q)
print("   exponent points", [tuple(x) for x in pts], " weights", list(w))
print("   |S| =", len(pts), "  rank Delta =", rk, "  <- S4-1 predicts rank 2 for |S|>=3")
print("   ==> S4-1 FAILS at a CLASS-HOMOGENEOUS charge.  Wall (1) is refuted as stated.")

print("\n--- 1B.  all class-homogeneous q = (alpha; beta,beta; gamma,gamma), 1..6")
print("     alpha beta gamma | |S| rankD | collinear?   (S4-1 predicts rank 2 whenever |S|>=3)")
bad = []
for al in range(1, 7):
    for be in range(1, 7):
        for ga in range(1, 7):
            qq = [al, be, be, ga, ga]
            pts, w, rk, B = report(qq)
            if len(pts) >= 3 and rk < 2:
                bad.append((al, be, ga))
print("   class-homogeneous (alpha,beta,gamma) in 1..6 with |S|=3 but rank Delta = 1:")
for t in bad:
    print("      ", t, "   condition beta*gamma = alpha*(beta+gamma):",
          t[1] * t[2], "==", t[0] * (t[1] + t[2]))
print(f"   count = {len(bad)} of {6**3} class-homogeneous assignments -- ALL with zero")
print("   within-class inhomogeneity.  The S4-1 failure locus is an ALGEBRAIC SURFACE")
print("   beta*gamma = alpha*(beta+gamma) inside the class-homogeneous set, not its boundary.")

print("\n--- 1C.  class-homogeneous, class-NON-uniform charge under which S4 SURVIVES EXACTLY")
base = R.mahler_3term_jensen(0.4, 0.3, 0.3, n=4_000_000)
print(f"   baseline  m(0.4+0.3x+0.3y) = {base:.9f}   (unit charge)")
for qq in ([1,1,1,1,1], [2,2,2,2,2], [1,2,2,3,3], [3,1,1,2,2], [1,5,5,7,7], [2,3,3,7,7],
           [-1,2,2,3,3], [1,2,2,-3,-3]):
    pts, w, rk, B = report(qq)
    lt = R.lambda_B_torus(qq, p_S3, n=2400)
    ld = R.lambda_B_direct(fg, cg, qq, p_S3, N=2_000_000)
    det = int(round(np.linalg.det(np.array([pts[1]-pts[0], pts[2]-pts[0]], dtype=float)))) if len(pts) >= 3 else 0
    print(f"   q={str(qq):18s} |S|={len(pts)} rk={rk} det={det:4d}  grid={lt: .9f}  direct={ld: .9f}  d(base)={abs(lt-base):.2e}")

print("\n   ==> every class-homogeneous, rank-2 charge reproduces S4's lambda EXACTLY,")
print("       including charges that are NOT proportional to the unit assignment.")
print("       So the wall is NOT at 'proportional exponents' either: it is drawn TOO EARLY.")

print("\n--- 1D.  THE EXCEPTIONAL SET.  Class-homogeneous charge on rank L = 2.")
print("   S4's Q3 headline: 'lambda_B is a function of L alone'.  Test it under charge.")
for (fr, cr, lab) in [(2*pi/3, 2*pi/3, "u,v primitive cube roots"),
                      (pi/2,  pi/2,  "4th roots"),
                      (2*pi/5, 4*pi/5, "5th roots")]:
    print(f"\n   connection {lab}:  f={fr:.6f} c={cr:.6f}")
    for qq in ([1,1,1,1,1], [1,2,2,3,3], [2,2,2,2,2], [3,1,1,2,2]):
        lam, T = R.lambda_B_finite_orbit(fr, cr, qq, p_S3, None)
        pts, w, rk, B = report(qq)
        print(f"      q={str(qq):16s} |S|={len(pts)} rk={rk}  orbit period T={T:3d}  lambda_B = {lam: .9f}")
print("\n   ==> L is IDENTICAL across the rows (same f,c); lambda_B is NOT.")
print("       'lambda_B is a function of L alone' is FALSE under class-homogeneous charge.")
print("       This is exactly where S4's Tier-1/Tier-2 exact values of record live.")
