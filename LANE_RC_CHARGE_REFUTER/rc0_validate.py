"""RC-0  VALIDATE.  Reproduce the corpus before attacking it.
Publishes d1, d2, grid conventions, seeds."""
import numpy as np
from math import pi, sqrt, log
import rclib as R

np.set_printoptions(linewidth=140)
print("=" * 78)
print("RC-0  VALIDATION -- carrier, incidence, closed form, S4 baselines")
print("=" * 78)

d1 = R.d1_matrix(); d2 = R.d2_matrix()
print("\nd1 (5 x 6), rows v0..v4, cols e1..e6:")
print(d1)
print("\nd2 (6 x 1), rows e1..e6, col F:")
print(d2.T, "(transposed for display)")
print("\nd1 @ d2 =", (d1 @ d2).ravel(), "   -> d^2 = 0 :", bool(np.all(d1 @ d2 == 0)))
print("rank d1 =", np.linalg.matrix_rank(d1), "  rank d2 =", np.linalg.matrix_rank(d2))
b0 = 5 - np.linalg.matrix_rank(d1)
b1 = (6 - np.linalg.matrix_rank(d1)) - np.linalg.matrix_rank(d2)
b2 = 1 - np.linalg.matrix_rank(d2)
print(f"b0={b0} b1={b1} b2={b2}   chi = 5-6+1 = {5-6+1}   (S1 sec.2: b0=1 b1=1 b2=0)")

print("\na_v (on gamma_F) =", R.A_INC, "   b_v (on gamma_C) =", R.B_INC)
print("classes: (1,1)={v0}  (1,0)={v1,v2}  (0,1)={v3,v4}   -- S4 sec.2")

# --- closed form vs INDEPENDENT direct matrix action on C^5
print("\n-- Z_k : closed form vs direct 5x5 matrix action, random (f,c,k,q,p)")
rng = np.random.default_rng(11223344)          # SEED PUBLISHED
worst = 0.0
for trial in range(4000):
    f = rng.uniform(0, 2 * pi); c = rng.uniform(0, 2 * pi)
    k = int(rng.integers(1, 12))
    q = rng.integers(-3, 4, size=5)
    p = rng.random(5); p /= p.sum()
    z1 = R.Z_closed(k, f, c, q, p)
    z2 = R.Z_direct(k, f, c, q, p, seed=int(rng.integers(0, 10**9)))
    worst = max(worst, abs(z1 - z2))
print(f"   4000 trials, max |closed - direct| = {worst:.3e}")

# --- S1's own published connection, unit charge, W-01's overlap
f0, c0 = pi, 3 * pi / 2                        # S1 sec.6: W_F=-1, W_C=-i
p_S1 = np.array([0.5, 0, 0, 0.25, 0.25])
z = R.Z_closed(1, f0, c0, [1] * 5, p_S1)
print(f"\n-- W-01 firing on S1's own connection: Z_1 = {z:.12g}  |Z_1| = {abs(z):.3e}")

# --- S3/S4 ready state, generic connection, unit charge
p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
fg, cg = 1.0, sqrt(2.0)                        # GENERIC connection used throughout
rk, gens = R.relation_lattice_rank(fg, cg, mmax=120)
print(f"\n-- generic connection f=1.0, c=sqrt(2):  rank L = {rk}  (searched |m|,|n|<=120)")

m_exact = R.mahler_3term_jensen(0.4, 0.3, 0.3, n=4_000_000)
print(f"   m(0.4 + 0.3x + 0.3y)  Jensen-in-y, 4e6 nodes  = {m_exact:.9f}")
print(f"   register of record                            = -0.767507880")
lg = R.lambda_B_torus([1]*5, p_S3, n=3000)
print(f"   2-D midpoint grid 3000^2                      = {lg:.9f}")
ld = R.lambda_B_direct(fg, cg, [1]*5, p_S3, N=2_000_000)
print(f"   direct schedule-B, N=2e6, f=1,c=sqrt2         = {ld:.9f}")

# --- the resonant headline connection of S3
fr, cr = 2.0, 1.1
print(f"\n-- S3 headline f=2.0 c=1.1 : -11f+20c = {-11*fr+20*cr:.3f}  (resonant, erratum)")
