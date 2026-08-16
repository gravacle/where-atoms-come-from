"""
R-3: (a) CONTROL -- at connections whose orbit closure is INFINITE the claim
         survives every alpha I can test.  Needs exact reduction of k*f mod 2pi
         because k reaches 1e15+; float64 loses the angle entirely.  Decimal at
         50 digits with pi to 50 places (declared in the IMPORT AUDIT).
     (b) FAMILY SCAN -- the break is not one example.  For every finite orbit
         order m and every alpha = 2, lambda(n^2) != lambda_B whenever the
         square map on Z/m is not measure preserving, i.e. m >= 3.
"""
from decimal import Decimal, getcontext
import numpy as np
from rx_lib import logabsZ, Z_closed

getcontext().prec = 60
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
TWOPI = 2 * PI
P_S3 = [0.4, 0.15, 0.15, 0.15, 0.15]


def lam_exact_angles(N, alpha, f, c, p):
    """(1/N) sum log|Z_{k_n}| with k_n = floor(n^alpha), angles reduced mod 2pi
    at 60 digits before being handed to float."""
    fd, cd = Decimal(repr(f)), Decimal(repr(c))
    p0, q, r = p[0], p[1] + p[2], p[3] + p[4]
    A = Decimal(repr(alpha))
    tot = 0.0
    for n in range(1, N + 1):
        k = n ** int(alpha) if float(alpha).is_integer() else int(Decimal(n) ** A)
        K = Decimal(k)
        th_u = float(-(K * fd) % TWOPI)      # arg u^k ,  u = e^{-if}
        th_v = float((K * cd) % TWOPI)       # arg v^k ,  v = e^{+ic}
        z = (np.exp(1j * (th_u + th_v)) * p0
             + np.exp(1j * th_u) * q
             + np.exp(1j * th_v) * r)
        tot += np.log(abs(z))
    return tot / N


print("=" * 78)
print("(a) CONTROL at INFINITE-orbit connections, exact angle reduction")
print("=" * 78)
for (f, c, tag, ref) in [
        (2.0, 1.1, "S3 headline, EXACTLY RESONANT (-11f+20c=0), subtorus", -0.767014993),
        (2.0, 1.1000001, "generic / non-resonant", -0.767507880)]:
    print("\n f=%.7f c=%.7f   %s" % (f, c, tag))
    print("   reference lambda_B of record = %.9f" % ref)
    for a, N in [(1.0, 200000), (1.5, 200000), (2.0, 200000), (3.0, 100000),
                 (4.0, 60000), (0.5, 200000), (0.25, 200000)]:
        v = lam_exact_angles(N, a, f, c, P_S3)
        print("   alpha=%4.2f  N=%7d  lambda = %.9f   dev = %+.6f"
              % (a, N, v, v - ref))

print()
print("=" * 78)
print("(b) FAMILY SCAN -- finite orbit of order m, u = e^{2pi i s/m}, v = e^{2pi i t/m}")
print("    lambda_B = mean over Z/m ; lambda(n^2) = mean over squares mod m")
print("    p = S3 ready state (0.4, .15,.15,.15,.15).  EXACT arithmetic on Z/m.")
print("=" * 78)
print("  m   s   t    lambda_B      lambda(n^2)    gap        lambda(n^3)   gap")
nbreak = 0
ntot = 0
for m in range(2, 17):
    for s in range(m):
        for t in range(m):
            # u = e^{-i f} = e^{2pi i s/m} -> f = -2pi s/m ; v = e^{i c} = e^{2pi i t/m}
            f = -2 * np.pi * s / m
            c = 2 * np.pi * t / m
            if s == 0 and t == 0:
                continue                       # trivial connection, excluded
            # true order of (u,v)
            from math import gcd
            order = m // gcd(gcd(s, t), m)
            if order != m:
                continue                       # count each order once, at its own m
            G = logabsZ(np.arange(m, dtype=float), f, c, P_S3)
            lamB = float(np.mean(G))
            sq = np.array([pow(n, 2, m) for n in range(1, 20 * m + 1)])
            cu = np.array([pow(n, 3, m) for n in range(1, 20 * m + 1)])
            l2 = float(np.mean(G[sq]))
            l3 = float(np.mean(G[cu]))
            ntot += 1
            if abs(l2 - lamB) > 1e-9 or abs(l3 - lamB) > 1e-9:
                nbreak += 1
            if (s, t) in [(1, 1), (2, 1), (1, 2)] and m <= 8:
                print("  %2d  %2d  %2d  %11.6f  %11.6f  %+9.6f  %11.6f  %+9.6f"
                      % (m, s, t, lamB, l2, l2 - lamB, l3, l3 - lamB))
print("\n  connections of order m in 2..16 scanned : %d" % ntot)
print("  of these, alpha=2 or alpha=3 DIFFERS from lambda_B : %d  (%.1f%%)"
      % (nbreak, 100.0 * nbreak / ntot))

print()
print("  order-2 connections (the ONLY finite orbits where n^2 is uniform mod m):")
for (s, t) in [(1, 0), (0, 1), (1, 1)]:
    f, c = -2 * np.pi * s / 2, 2 * np.pi * t / 2
    G = logabsZ(np.arange(2.0), f, c, P_S3)
    sq = np.array([pow(n, 2, 2) for n in range(1, 41)])
    print("    (s,t)=(%d,%d)  lambda_B=%.9f  lambda(n^2)=%.9f  gap=%.2e"
          % (s, t, float(np.mean(G)), float(np.mean(G[sq])),
             abs(float(np.mean(G[sq])) - float(np.mean(G)))))
