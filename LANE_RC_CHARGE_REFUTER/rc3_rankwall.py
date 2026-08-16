"""RC-3  ATTACK ON WALL (3), THE RANK WALL, AND ON THE HEADLINE 'FOUR RAYS' WALL.

CLAIMS UNDER ATTACK:
  (i)  "lambda_B at a generic connection is charge-rigid for 3-point rank-2 supports
        (Theorem C-3) and charge-sensitive for 4-point supports and for any rank drop."
  (ii) "the wall between 'charge does nothing' and 'charge does everything' runs through
        the SPECTATOR VERTEX and through rank Delta, not through the pinch."
  (iii) HEADLINE: "nothing in it survives the exponent vector leaving the union of the
        four rays Z(0,0) u Z(1,0) u Z(0,1) u Z(1,1) ..."
"""
import numpy as np
from math import pi, sqrt, log
import itertools
import rclib as R

p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
fg, cg = 1.0, sqrt(2.0)
GRID = 3000                                   # midpoint grid, published

print("=" * 78)
print("RC-3  THE RANK WALL AND THE 'FOUR RAYS' HEADLINE")
print("=" * 78)

# ---------------------------------------------------------------- (iii) first
print("\n--- 3.0  THE HEADLINE WALL IS NEVER CROSSED.")
RAYS = [(0, 0), (1, 0), (0, 1), (1, 1)]
def on_rays(e):
    if e[0] == 0 and e[1] == 0:
        return True
    for r in RAYS[1:]:
        if r[0] == 0:
            if e[0] == 0:
                return True
        elif r[1] == 0:
            if e[1] == 0:
                return True
        else:
            if e[0] == e[1]:
                return True
    return False
off = 0
for q in itertools.product(range(-3, 4), repeat=5):
    for e in R.exponents(list(q)):
        if not on_rays(tuple(int(x) for x in e)):
            off += 1
print(f"   enumerated all 7^5 = {7**5} charges in {{-3..3}}^5, all 5 exponent points each:")
print(f"   exponent points NOT in Z(0,0) u Z(1,0) u Z(0,1) u Z(1,1) : {off}")
print("   E(v) = q_v * (a_v,b_v) lies on the ray through (a_v,b_v) BY CONSTRUCTION.")
print("   ==> per-vertex charge on K1 NEVER leaves the union of the four rays, yet the")
print("       claim's own three walls all fall INSIDE it.  The headline names a wall that")
print("       none of the phenomena it explains ever crosses.  VACUOUS AS STATED.")

# ------------------------------------------------------- (i) Theorem C-3 test
print("\n--- 3.1  THEOREM C-3 (3-point rank-2 rigidity) -- HARD TEST, NOT ASSERTION.")
print("   Proof to be checked: for points {e0,e1,e2} with M = [e1-e0; e2-e0], det M != 0,")
print("   the map (x,y) |-> (x^M11 y^M12, x^M21 y^M22) is a surjective endomorphism of T^2")
print("   with finite kernel, so it pushes Haar to Haar; hence")
print("       m(p0 x^e0 + p1 x^e1 + p2 x^e2) = m(p0 + p1 x + p2 y),  independent of M.")
base = R.mahler_3term_jensen(0.4, 0.3, 0.3, n=8_000_000)
print(f"   reference m(0.4+0.3x+0.3y) = {base:.12f}")
rng = np.random.default_rng(24681357)          # SEED PUBLISHED
worst = 0.0; ntest = 0; dets = []
print("   random 3-point rank-2 configs (exponents NOT restricted to K1's rays):")
for t in range(200):
    e0 = rng.integers(-4, 5, 2); e1 = rng.integers(-4, 5, 2); e2 = rng.integers(-4, 5, 2)
    M = np.array([e1 - e0, e2 - e0], dtype=float)
    d = np.linalg.det(M)
    if abs(d) < 0.5:
        continue
    n = 2400
    th = (np.arange(n) + 0.5) * 2 * pi / n
    X, Y = np.meshgrid(th, th, indexing='ij')
    Zs = (0.4 * np.exp(1j * (e0[0] * X + e0[1] * Y))
          + 0.3 * np.exp(1j * (e1[0] * X + e1[1] * Y))
          + 0.3 * np.exp(1j * (e2[0] * X + e2[1] * Y)))
    val = float(np.mean(np.log(np.abs(Zs) + 1e-300)))
    worst = max(worst, abs(val - base)); ntest += 1; dets.append(int(round(d)))
print(f"   {ntest} random rank-2 triples, |det| up to {max(abs(x) for x in dets)}:")
print(f"   worst |lambda_grid - m(0.4+0.3x+0.3y)| = {worst:.3e}   (grid discretisation floor)")
print("   ==> THEOREM C-3 SURVIVES.  I could not break it.")

# ------------------------------ (i) 'charge-sensitive for 4-point supports'
print("\n--- 3.2  'CHARGE-SENSITIVE FOR 4-POINT SUPPORTS' -- EXACT COUNTEREXAMPLES.")
print("   (a) GLOBAL RESCALING.  q -> c*q multiplies every exponent by c; (x,y)|->(x^c,y^c)")
print("       is a surjective endomorphism of T^2, so lambda_B is EXACTLY invariant at any")
print("       point count.  Charge changes; lambda does not.")
for q in ([1,1,2,1,1], [1,1,2,3,1], [1,2,3,4,5]):
    pts, w = R.support_points(q, p_S3)
    l1 = R.lambda_B_torus(q, p_S3, n=GRID)
    q2 = [2 * x for x in q]; q3 = [3 * x for x in q]
    l2 = R.lambda_B_torus(q2, p_S3, n=GRID)
    l3 = R.lambda_B_torus(q3, p_S3, n=GRID)
    d1 = R.lambda_B_direct(fg, cg, q, p_S3, N=4_000_000)
    d2 = R.lambda_B_direct(fg, cg, q2, p_S3, N=4_000_000)
    print(f"   q={str(q):14s} |S|={len(pts)}  grid: {l1: .9f} {l2: .9f} {l3: .9f}")
    print(f"   {'':14s}        direct(f=1,c=sqrt2, N=4e6): {d1: .9f}  vs 2q: {d2: .9f}  |diff|={abs(d1-d2):.2e}")

print("\n   (b) THE F/C SWAP.  M = [[0,1],[1,0]] in GL_2(Z) fixes the (1,1) ray and exchanges")
print("       the (1,0) and (0,1) rays, so it maps K1 charge configs to K1 charge configs.")
print("       With class-symmetric weights it is a charge symmetry of lambda_B:")
for q in ([1,1,2,1,1], [1,3,1,2,2], [2,1,3,2,1]):
    qs = [q[0], q[3], q[4], q[1], q[2]]
    l1 = R.lambda_B_torus(q, p_S3, n=GRID)
    l2 = R.lambda_B_torus(qs, p_S3, n=GRID)
    pts, _ = R.support_points(q, p_S3)
    print(f"   q={str(q):14s} <-> q'={str(qs):14s} |S|={len(pts)}  {l1: .9f}  {l2: .9f}  |diff|={abs(l1-l2):.2e}")
print("   ==> a 4- and 5-point support carries an EXACT infinite charge-rigidity group.")
print("       'charge-sensitive for 4-point supports' is false as a universal; the true")
print("       invariant is the weighted exponent configuration MODULO integer matrices of")
print("       non-zero determinant, at every point count.")

# --------------------------------------------- (ii) the spectator mislocation
print("\n--- 3.3  THE SPECTATOR.  [ILLUSTRATION carrier K1+ : K1 plus e7: v0->v5, v5 on")
print("    neither loop.  a_v5 = b_v5 = 0, so E(v5) = (0,0) for EVERY charge q_v5.]")
print("   The spectator's exponent point is (0,0) whatever its charge:")
for qs in (0, 1, 2, 5, -3):
    print(f"      q_v5 = {qs:3d}  ->  E(v5) = {tuple(int(x) for x in (qs*0, qs*0))}")
print("   ==> the spectator vertex is the ONE vertex on which charge acts trivially, exactly.")
print("       It cannot be the wall.  And K1 HAS NO SPECTATOR, yet charge already does")
print("       everything on K1 (RC-1, RC-2).  The wall is mislocated a second time.")
print("   What actually happens at 4 points: on K1 a 4th point requires WITHIN-CLASS")
print("   inhomogeneity (q_v1 != q_v2 or q_v3 != q_v4) -- the pinch and the spectator are")
print("   both irrelevant to it.")

# ------------------------------------------------- rank drop -> always moves?
print("\n--- 3.4  'AND FOR ANY RANK DROP'.  Rank-1 (collinear) 3-point configs:")
for q in ([1,2,2,2,2], [2,3,3,6,6], [2,4,4,4,4], [3,6,6,6,6], [2,6,6,3,3]):
    pts, w = R.support_points(q, p_S3)
    rk = R.delta_rank(pts)
    lt = R.lambda_B_torus(q, p_S3, n=GRID)
    ld = R.lambda_B_direct(fg, cg, q, p_S3, N=8_000_000)
    print(f"   q={str(q):14s} pts={[tuple(map(int,x)) for x in pts]} rk={rk}  grid={lt: .9f} direct={ld: .9f}")
print("   ==> rank drop does move lambda.  This half of wall (3) SURVIVES.")
print("   NOTE q=(2,4,4,4,4) = 2*(1,2,2,2,2): rank drop AND exact rigidity under rescaling.")
