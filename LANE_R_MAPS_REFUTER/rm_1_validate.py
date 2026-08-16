"""
rm_1_validate.py -- LANE R (MAPS REFUTER) validation gate.
Reproduce, from scratch, every S4 number this refutation touches.
If this file does not reproduce S4, nothing below it is admissible.
"""
import numpy as np, math, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from rm_lib import *

print("=" * 78)
print("RM-1  VALIDATION GATE -- reproducing S4 from scratch")
print("=" * 78)

print("\n--- 1. PUBLISHED INCIDENCE MATRICES (S4 published none) ---")
for key in ["B1", "B1s", "B1p", "B1q", "B0a", "B0b"]:
    K = ALL[key]()
    b0, b1, b2, r1, r2 = K.betti()
    dd = K.d1d2()
    print(f"\n[{key}] {K.name}   V={K.nV} E={K.nE} F={K.nF} chi={K.chi()}")
    print(f"   edges (tail,head): {K.edges}")
    print(f"   faces (edge,sign): {K.faces}")
    print(f"   gamma_F edges = {K.gF}    gamma_C edges = {K.gC}")
    print(f"   gamma_F is a cycle: {K.is_cycle(K.gF)}   gamma_C is a cycle: {K.is_cycle(K.gC)}")
    print(f"   rank d1={r1} rank d2={r2}   b0={b0} b1={b1} b2={b2}   "
          f"max|d1.d2| = {0 if dd.size == 0 else int(np.max(np.abs(dd)))}")
    print(f"   d1 =")
    for row in K.d1():
        print("       ", row)
    if K.nF:
        print(f"   d2^T (rows = faces) =")
        d2 = np.array(K.d2()).T
        for row in d2:
            print("       ", list(map(int, row)))
    cc = K.class_counts()
    print(f"   class counts {{(a,b):n}} = " +
          str({f"{a}{b}": n for (a, b), n in sorted(cc.items())}))
    pi = K.pi_uniform()
    print(f"   SENSE U pi = (p00,p10,p01,p11) = "
          f"({pi[0]:.9f},{pi[1]:.9f},{pi[2]:.9f},{pi[3]:.9f})")

print("\n\n--- 2. S4 TABLE (:574-583) lambda_B SENSE U, REPRODUCED THREE WAYS ---")
S4_U = {"B0a": -0.747659833, "B0b": -0.810930216, "B1": -0.756573586,
        "B1p": -0.693147181, "B1q": -0.741029583, "B1s": -0.724759919}
print(f"{'car':5s} {'closed form':>16s} {'jensen quad':>16s} {'2D quad':>16s} "
      f"{'S4':>14s} {'|dev vs S4|':>12s}")
worst = 0.0
for key in ["B0a", "B0b", "B1", "B1p", "B1q", "B1s"]:
    K = ALL[key]()
    pi = K.pi_uniform()
    lc = lambda_B_closed(pi)
    lj = mahler_jensen(pi, n=400000)
    l2 = mahler_2d(pi, n=4000)
    dev = abs(lc - S4_U[key])
    worst = max(worst, dev)
    print(f"{key:5s} {lc:16.12f} {lj:16.12f} {l2:16.12f} {S4_U[key]:14.9f} {dev:12.3e}")
print(f"worst |closed - S4| over 6 carriers = {worst:.3e}")

print("\n--- 3. DIRECT SCHEDULE-B SIMULATION (independent of any Mahler identity) ---")
f0, c0 = 1.0, math.sqrt(2.0)     # S4's stated check point, :605
print(f"(f,c) = (1.0, sqrt(2));  N = 2e6;  relation lattice rank = "
      f"{relation_lattice_rank(f0, c0)[0]}")
for key in ["B0a", "B0b", "B1", "B1p", "B1q", "B1s"]:
    K = ALL[key]()
    pi = K.pi_uniform()
    d = lambda_B_direct(pi, f0, c0, N=2000000)
    lc = lambda_B_closed(pi)
    print(f"   {key:5s} direct {d:14.9f}   closed {lc:14.9f}   dev {abs(d-lc):.2e}")

print("\n--- 4. THE CONVENTION CHECK (S4 :155 uses cos(k(f+c)); this lane cos(k(f-c))) ---")
rng = np.random.Generator(np.random.PCG64(20260816))   # SEED 20260816
mx = 0.0
for _ in range(2000):
    f, c = rng.uniform(0, 2 * np.pi, 2)
    k = int(rng.integers(1, 50))
    pi = (0.0, 0.4, 0.4, 0.2)
    z_mine = abs(Z_from_pi(pi, f, c, k))
    p0, q, r = pi[3], pi[1], pi[2]
    s4 = math.sqrt(max(0.0, p0*p0 + q*q + r*r + 2*p0*q*math.cos(k*c)
                       + 2*p0*r*math.cos(k*f) + 2*q*r*math.cos(k*(f + c))))
    # S4's formula with c -> -c must equal mine
    s4_flip = math.sqrt(max(0.0, p0*p0 + q*q + r*r + 2*p0*q*math.cos(-k*c)
                            + 2*p0*r*math.cos(k*f) + 2*q*r*math.cos(k*(f - c))))
    mx = max(mx, abs(z_mine - s4_flip))
print(f"   seed 20260816, 2000 samples: max| |Z_k|_mine - S4(:155) with c->-c | = {mx:.3e}")
print("   -> the sign difference is a labelling convention, not an error. Declared.")

print("\n--- 5. Z_k FROM VERTICES == Z_k FROM CLASS PUSHFORWARD (the definition check) ---")
rng = np.random.Generator(np.random.PCG64(777))        # SEED 777
mx = 0.0
for key in ["B0a", "B0b", "B1", "B1p", "B1q", "B1s"]:
    K = ALL[key]()
    p = [1.0 / K.nV] * K.nV
    pi = K.pi_from_p(p)
    for _ in range(200):
        f, c = rng.uniform(0, 2 * np.pi, 2)
        k = int(rng.integers(1, 200))
        mx = max(mx, abs(Z_from_vertices(K, p, f, c, k) - Z_from_pi(pi, f, c, k)))
print(f"   seed 777, 1200 samples: max|Z_vertexwise - Z_classwise| = {mx:.3e}")

print("\n--- 6. THE EXACT ROWS S4 CLAIMS (:590-596) ---")
print(f"   m(2/9 + 3/9 y + 4/9 x)   = {mahler_cassaigne_maillot(2/9,3/9,4/9):.12f}   "
      f"S4 -0.747659833081")
print(f"   m(0.4 + 0.4x + 0.2y)     = {mahler_cassaigne_maillot(0.2,0.4,0.4):.12f}   "
      f"S4 -0.756573585640")
print(f"   m(1/7 + 3/7x + 3/7y)     = {mahler_cassaigne_maillot(1/7,3/7,3/7):.12f}   "
      f"S4 -0.741029582571")
print(f"   m(5/11 + 5/11x + 1/11y)  = {mahler_cassaigne_maillot(1/11,5/11,5/11):.12f}  "
      f"S4 -0.724759919461")
print(f"   m(0.4 + 0.3x + 0.3y)     = {mahler_cassaigne_maillot(0.3,0.3,0.4):.12f}   "
      f"S4 -0.767507880358  (SENSE C, 3 classes)")
print(f"   B1p Jensen log max(1/2,1/2) = {math.log(0.5):.12f}   S4 -0.693147180560")

print("\n--- 7. S4 CONTROL 3 AS S4 STATED IT ---")
b1 = ALL["B1"](); b1s = ALL["B1s"]()
lu_b1 = lambda_B_closed(b1.pi_uniform())
lu_b1s = lambda_B_closed(b1s.pi_uniform())
print(f"   SENSE U   B1 {lu_b1:.12f}   B1s {lu_b1s:.12f}   |diff| = {abs(lu_b1-lu_b1s):.3e}")
print(f"   S4 reported                                         |diff| = 3.181e-02")
lc3 = lambda_B_closed((0.0, 0.4, 0.3, 0.3))
print(f"   SENSE C   both = {lc3:.12f}   |diff| = 0.0e+00     (S4's own next line)")
print("\nVALIDATION GATE: PASSED if worst dev above is <= 1e-8 for closed forms.")
