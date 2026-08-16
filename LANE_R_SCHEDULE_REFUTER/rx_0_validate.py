"""R-0: validate my from-scratch implementation against the SEALED corpus."""
import numpy as np
from rx_lib import (check_complex, Z_matrix, Z_closed, logabsZ, schedule_alpha,
                    lam_of_schedule, orbit_order, lambda_B_finite_orbit, D1, D2)

np.set_printoptions(precision=6, suppress=True)
P_S3 = [0.4, 0.15, 0.15, 0.15, 0.15]

print("=== COMPLEX SANITY ===")
print(check_complex())
print("d1 =\n", D1)
print("d2^T =", D2.ravel())

print("\n=== Z_k: matrix action vs closed form (f=2.0,c=1.1,p=S3), k<=200 ===")
dev = max(abs(Z_matrix(k, 2.0, 1.1, P_S3) - Z_closed(k, 2.0, 1.1, P_S3))
          for k in range(0, 201))
print("max deviation =", dev, "   (S3 sec4.1 reports 1.538e-14)")

print("\n=== S3 sec 4.1 recurrence facts ===")
ks = np.arange(1, 401.0)
az = np.abs(Z_closed(ks, 2.0, 1.1, P_S3))
print("min |Z_k| over k<=400 = %.6f at k=%d   (S3: 0.024654 at k=42)"
      % (az.min(), ks[az.argmin()]))
ks2 = np.arange(1, 4001.0)
az2 = np.abs(Z_closed(ks2, 2.0, 1.1, P_S3))
print("sup |Z_k| over k<=4000 = %.6f at k=%d  (S3: 0.999941 at k=377)"
      % (az2.max(), ks2[az2.argmax()]))
print("count |Z_k|>0.99 among first 4000 =", int((az2 > 0.99).sum()), " (S3: 37)")

print("\n=== S3 sec 4.3 schedule-B finite stages, f=2.0 c=1.1 ===")
for N in (1, 2, 5, 10, 20, 42, 50, 100, 200, 400, 1000, 2000, 4000):
    print("  N=%6d  (1/N)log|Omega_N| = %.6f" % (N, lam_of_schedule(np.arange(1, N + 1), 2.0, 1.1, P_S3)))
print("  N=200000            -> %.6f   (S3: -0.767026)"
      % lam_of_schedule(np.arange(1, 200001), 2.0, 1.1, P_S3))

print("\n=== ERRATUM check: f=2.0,c=1.1 is EXACTLY RESONANT, -11f+20c = %.3e ==="
      % (-11 * 2.0 + 20 * 1.1))
print("  subtorus average (register: -0.767014993):")
# subtorus: (u^k,v^k) = (e^{-2ik}, e^{1.1ik}); relation u^11 v^20 = 1.
# parametrise H = {(t^20, t^-11)} ; k -> t = t0^k with t0 chosen so
# (t0^20,t0^-11) = (u,v).  Just average log|Z| over H directly.
T = 4_000_000
th = (np.arange(T) + 0.5) / T * 2 * np.pi     # midpoint grid on the circle H
x = np.exp(1j * 20 * th)      # plays u
y = np.exp(-1j * 11 * th)     # plays v
val = np.log(np.abs(x * y * 0.4 + x * 0.3 + y * 0.3))
print("  Haar mean over the subtorus H = ker(chi_{11,20}) : %.9f" % val.mean())
print("  generic full-torus value m(0.4+0.3x+0.3y) (register: -0.767507880):")
G = 3000
tt = (np.arange(G) + 0.5) / G * 2 * np.pi
X, Y = np.meshgrid(np.exp(1j * tt), np.exp(1j * tt), indexing="ij")
print("   %.9f" % np.log(np.abs(0.4 * X * Y + 0.3 * X + 0.3 * Y)).mean())

print("\n=== S3 sec 6(f) table, reproduced ===")
for (f, c, ref) in [(2.0, 1.1, -0.767026), (2.0, 2.0, -1.203587),
                    (3.14159, 1.57080, -0.804719), (3.14159, 4.71239, -0.804719)]:
    v = lam_of_schedule(np.arange(1, 200001), f, c, P_S3)
    print("  f=%.5f c=%.5f  lambda_B(N=2e5) = %.6f   S3 says %.6f" % (f, c, v, ref))

print("\n=== EXACT finite-orbit connections (pi and pi/2, pi and 3pi/2) ===")
for (f, c) in [(np.pi, np.pi / 2), (np.pi, 3 * np.pi / 2)]:
    m = orbit_order(f, c)
    print("  f=%.6f c=%.6f  orbit order m=%s  lambda_B(exact cyclic mean)=%.9f"
          % (f, c, m, lambda_B_finite_orbit(m, f, c, P_S3)))
    print("     |Z_k| for k=0..%d :" % (m - 1),
          np.round(np.abs(Z_closed(np.arange(m, dtype=float), f, c, P_S3)), 9))
