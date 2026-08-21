"""
R-1: ATTACK on  "lambda(floor(n^alpha)) = lambda_B for EVERY alpha > 0".

Strategy.  The claim's proof (THEOREM R3-2, Abel summation) is stated for
alpha = 1/q, q a positive integer -- i.e. alpha <= 1, where floor(n^alpha)
is a NON-DECREASING SURJECTION onto N with block lengths c_m = (m+1)^q - m^q.
For alpha > 1 floor(n^alpha) is a strictly increasing SUBSEQUENCE: it SKIPS
integers.  Skipping is invisible when the orbit closure is a connected torus
(Weyl's polynomial theorem covers it) but it is NOT invisible when the closure
is a FINITE cyclic group -- squares/cubes are not equidistributed mod m.

S3 sec 6(f) puts two finite-orbit connections on the record.  Take them.
"""
import numpy as np
from rx_lib import (Z_closed, logabsZ, schedule_alpha, lam_of_schedule,
                    lam_alpha, orbit_order, lambda_B_finite_orbit)

np.seterr(invalid="ignore")  # T-35 NONDET fix: the expected floor(n^alpha)>int64 cast RuntimeWarning is stderr whose position vs stdout varies with capture buffering; seterr changes reporting only, never a computed value

P_S3 = [0.4, 0.15, 0.15, 0.15, 0.15]
GEN = (2.0, 1.1000001)   # a NON-resonant nearby connection, for contrast


def resid_hist(ks, m):
    r = np.bincount(np.asarray(ks) % m, minlength=m)
    return r / r.sum()


print("#" * 78)
print("# EXHIBIT 1 -- S3's own finite-orbit connection  f = pi, c = pi/2")
print("#" * 78)
f, c = np.pi, np.pi / 2
m = orbit_order(f, c)
lamB = lambda_B_finite_orbit(m, f, c, P_S3)
G = logabsZ(np.arange(m, dtype=float), f, c, P_S3)
print("orbit order m = %d ; log|Z_k|, k=0..3 = %s" % (m, np.round(G, 9)))
print("lambda_B (schedule k_n = n, EXACT cyclic mean) = %.9f" % lamB)
print("lambda_A (schedule k_n = 1)                    = %.9f" % G[1])
print()
print(" alpha |   lambda(floor(n^a)) N=1e6  |  N=1e7      | residue histogram mod 4")
for a in [0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]:
    k6 = schedule_alpha(10**6, a)
    k7 = schedule_alpha(10**7, a)
    print("  %4.2f | %22.9f | %11.9f | %s"
          % (a, lam_of_schedule(k6, f, c, P_S3),
             lam_of_schedule(k7, f, c, P_S3),
             np.round(resid_hist(k7, m), 6)))
print()
print("PREDICTION for integer alpha>=2 from elementary number theory:")
print("  n^2 mod 4 = 1 if n odd, 0 if n even  -> density (1/2, 1/2, 0, 0)")
print("  => lambda(n^2) = (log|Z_0| + log|Z_1|)/2 = %.9f" % ((G[0] + G[1]) / 2))
print("  n^3 mod 4 : n even -> 0 ; n=1,5,..->1 ; n=3,7,..->3  -> (1/2,1/4,0,1/4)")
print("  => lambda(n^3) = %.9f" % (0.5 * G[0] + 0.25 * G[1] + 0.25 * G[3]))
print("  n^4 mod 4 = 0 (n even), 1 (n odd) -> same as n^2 : %.9f" % ((G[0] + G[1]) / 2))
print("  lambda_B                                        = %.9f" % lamB)
print("  GAP at alpha=2 : %.9f" % (((G[0] + G[1]) / 2) - lamB))

print()
print("#" * 78)
print("# EXHIBIT 2 -- the same sweep at a NON-resonant connection (control)")
print("#" * 78)
f2, c2 = GEN
print("  f=%.7f c=%.7f   -11f+20c = %.3e (NOT resonant)"
      % (f2, c2, -11 * f2 + 20 * c2))
for a in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    print("  alpha=%4.2f  lambda(N=1e7) = %.9f"
          % (a, lam_of_schedule(schedule_alpha(10**7, a), f2, c2, P_S3)))
print("  full-torus m(0.4+0.3x+0.3y) = -0.767507880 (register)")

print()
print("#" * 78)
print("# EXHIBIT 3 -- AMPLIFIER: a finite orbit with an EXACT ZERO off the")
print("#              range of the alpha>=2 schedules.  lambda_B = -infinity,")
print("#              lambda(n^2) FINITE.")
print("#" * 78)
# u = e^{-i f} = i  and  v = e^{i c} = i   ->  f = -pi/2 (== 3pi/2), c = pi/2
f3, c3 = 3 * np.pi / 2, np.pi / 2
p3 = [0.5, 0.125, 0.125, 0.125, 0.125]      # p0 = 1/2, q = r = 1/4
m3 = orbit_order(f3, c3)
Zs = Z_closed(np.arange(m3, dtype=float), f3, c3, p3)
print("f = 3pi/2, c = pi/2 ; p = %s  (p0=1/2, q=r=1/4)" % p3)
print("u = e^{-if} = %s   v = e^{ic} = %s   orbit order m = %d"
      % (np.round(np.exp(-1j * f3), 12), np.round(np.exp(1j * c3), 12), m3))
print("Z_k, k=0..3 :", np.round(Zs, 12))
print("|Z_k|       :", np.round(np.abs(Zs), 12))
print("  -> Z_2 = 2*p0 - 1 = 0 EXACTLY.  Formation criterion: chi_0=uv=-1,")
print("     chi_F=u=i, chi_C=v=i ; G = <i> != {1}, so formation occurs.")
G3 = logabsZ(np.arange(m3, dtype=float), f3, c3, p3)
print("log|Z_k| :", G3)
print("lambda_B (k_n = n)   = %s   [S3-2: 'immediate exact formation']"
      % lambda_B_finite_orbit(m3, f3, c3, p3))
for a in [1.0, 2.0, 3.0, 4.0]:
    k = schedule_alpha(10**6, a)
    lam = lam_of_schedule(k, f3, c3, p3)
    hits2 = int((k % 4 == 2).sum())
    print("  alpha=%3.1f : lambda(N=1e6) = %-22s  cells landing on k=2 mod 4 : %d"
          % (a, "%.9f" % lam if np.isfinite(lam) else str(lam), hits2))
print("  Omega_N under k_n=n is EXACTLY 0 from N=2 onward;")
print("  under k_n=n^2 it is never 0 and decays at rate log(1/sqrt2)/2 = %.9f"
      % (0.5 * (G3[0] + G3[1])))

print()
print("#" * 78)
print("# EXHIBIT 4 -- where the wall actually is: alpha<=1 survives at the")
print("#              finite orbit; alpha>1 non-integer also survives;")
print("#              INTEGER alpha>=2 is the break.")
print("#" * 78)
f, c = np.pi, np.pi / 2
lamB = lambda_B_finite_orbit(4, f, c, P_S3)
rows = []
for a in [0.05, 0.1, 0.2, 1/3, 0.5, 0.75, 0.9, 0.99, 1.0, 1.01, 1.1, 1.25,
          1.5, 1.75, 1.9, 1.99, 2.0, 2.01, 2.5, 3.0, 3.01, 4.0]:
    k = schedule_alpha(10**7, a)
    rows.append((a, lam_of_schedule(k, f, c, P_S3), k.max()))
print(" alpha   lambda(N=1e7)    dev from lambda_B    k_max")
for a, l, km in rows:
    print("  %5.3f  %13.9f    %+14.9f   %d" % (a, l, l - lamB, km))
print("\n lambda_B = %.9f   lambda_A = %.9f" % (lamB, logabsZ(1.0, f, c, P_S3)))
