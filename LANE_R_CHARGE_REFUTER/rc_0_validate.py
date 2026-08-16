"""V0 -- validate my own implementation against the sealed corpus before attacking anything."""
import numpy as np
import rc_lib as R

print("=" * 96)
print("V0.1  CARRIER K1, PUBLISHED (d1, d2, d1 d2 = 0, Betti)")
print("=" * 96)
d1 = R.d1_matrix(5, R.K1_EDGES)
d2 = R.d2_matrix(6, R.K1_FACES)
print("d1 (rows v0..v4, cols e1..e6) =\n", d1)
print("d2^T (cols e1..e6) =", d2.T)
print("d1 @ d2 =", (d1 @ d2).ravel(), " max|.| =", int(np.max(np.abs(d1 @ d2))))
V, E, F = 5, 6, 1
b0 = V - np.linalg.matrix_rank(d1.astype(float))
b2 = F - np.linalg.matrix_rank(d2.astype(float))
b1 = (E - np.linalg.matrix_rank(d1.astype(float))) - np.linalg.matrix_rank(d2.astype(float))
print(f"chi = {V-E+F}   b0 = {b0}   b1 = {b1}   b2 = {b2}   (S1 section 2: 0,1,1,0)")

Ecls = R.class_vectors(5, R.K1_EDGES, R.K1_GAMMA_F, R.K1_GAMMA_C)
print("unit-charge exponent vectors E_v = (a_v,b_v):", [tuple(x) for x in Ecls])

print()
print("=" * 96)
print("V0.2  THE FUNCTIONAL: operator form vs closed form, and gauge invariance")
print("=" * 96)
rng = np.random.default_rng(90210001)
worst_mc, worst_g = 0.0, 0.0
for _ in range(400):
    q = rng.integers(-3, 4, size=5)
    Ev = (Ecls.T * q).T
    ph = rng.uniform(0, 2 * np.pi, 5)
    amp = rng.uniform(0.05, 1.0, 5)
    amp /= np.linalg.norm(amp)
    s = amp * np.exp(1j * ph)
    p = np.abs(s) ** 2
    f, c = rng.uniform(0, 2 * np.pi, 2)
    k = int(rng.integers(1, 8))
    worst_mc = max(worst_mc, abs(R.Zk_matrix(Ev, s, f, c, k) - R.Zk(Ev, p, f, c, k)))
    # gauge: a_e -> a_e + th_tgt - th_src leaves W_F, W_C (hence f,c) fixed; and the
    # fibrewise action s_v -> g_v^{q_v} s_v leaves p_v fixed.  Check both at once.
    th = rng.uniform(0, 2 * np.pi, 5)
    s2 = s * np.exp(1j * q * th)
    worst_g = max(worst_g, abs(R.Zk_matrix(Ev, s2, f, c, k) - R.Zk_matrix(Ev, s, f, c, k)))
print(f"  max |operator form - closed form|                     = {worst_mc:.3e}   (400 samples, seed 90210001)")
print(f"  max |Z_k(gauge-transformed section) - Z_k|            = {worst_g:.3e}")

print()
print("=" * 96)
print("V0.3  REPRODUCE THE CORPUS BEFORE ATTACKING IT")
print("=" * 96)
p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
f0, c0 = 1.0, np.sqrt(2.0)
lamA = np.log(abs(R.Zk(Ecls, p_S3, f0, c0, 1)))
print(f"  lambda_A at (f,c)=(1,sqrt2), unit charge   = {lamA:.12f}   [S4/S5 report -0.493553117]")
lamB = R.lambda_B_generic(Ecls, p_S3, Nx=16384)
print(f"  lambda_B generic, unit charge (Nx=16384)   = {lamB:.12f}   [record -0.767507880]")
lamBd = R.lambda_direct(Ecls, p_S3, f0, c0, 400000)
print(f"  lambda_B direct, N=4e5, (1,sqrt2)          = {lamBd:.12f}")
print(f"  S1 section 6 instance W_F=-1, W_C=-i, p=(1/2,0,0,1/4,1/4):")
p_S1 = np.array([0.5, 0.0, 0.0, 0.25, 0.25])
z = R.Zk(Ecls, p_S1, np.pi, 3 * np.pi / 2, 1)
print(f"      Z_1 = {z:.15f}   |Z_1| = {abs(z):.3e}   [W-01: exactly 0]")

# the exceptional carrier value of record, log(4/9), on B0b's exponent profile
print()
print("  cross-check of an exact rate of record (W-03: B0b = log(4/9)):")
Eb0b = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [1, 0], [1, 0], [1, 1], [1, 1], [0, 1]])
pb0b = np.ones(9) / 9.0
print(f"      B0b SENSE-U lambda_B (Nx=32768) = {R.lambda_B_generic(Eb0b, pb0b, Nx=32768):.12f}"
      f"   log(4/9) = {np.log(4/9):.12f}")
