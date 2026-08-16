"""
R-2: EXACT-ARITHMETIC redo.  In rx_1 every row with alpha >= 2.5 was garbage:
floor(n^alpha) overflowed int64 (k_max printed as 2^63-1) and had already left
float64's exact-integer range 2^53.  This file recomputes with exact integers.

At a FINITE-ORBIT connection only  k mod m  matters, and that is exact:
   integer alpha       : k mod m = pow(n, alpha, m)            (exact, small)
   alpha = p/2         : k = isqrt(n^p)                        (exact, big int)
   other alpha         : k = int(Decimal(n)**Decimal(alpha)) at 50 digits
N is chosen per row so that the exact computation finishes: the limit is a
density statement about residues and 1e5 samples already pin the densities to
1e-3, while the integer-alpha rows (the ones that carry the refutation) are
run at 1e6 with modular exponentiation.
"""
from decimal import Decimal, getcontext
from math import isqrt
import numpy as np
from rx_lib import logabsZ, orbit_order

getcontext().prec = 50
P_S3 = [0.4, 0.15, 0.15, 0.15, 0.15]


def kmod(N, alpha, m, half=None):
    out = np.empty(N, dtype=np.int64)
    if float(alpha).is_integer():
        a = int(alpha)
        if a == 0:
            out[:] = 1 % m
            return out
        for n in range(1, N + 1):
            out[n - 1] = pow(n, a, m)
        return out
    if half is not None:                      # alpha = half/2, half odd int
        for n in range(1, N + 1):
            out[n - 1] = isqrt(n ** half) % m
        return out
    A = Decimal(repr(alpha))
    for n in range(1, N + 1):
        out[n - 1] = int(Decimal(n) ** A) % m
    return out


ROWS = [  # (alpha, N, half-exponent or None)
    (0.00, 1_000_000, None),
    (0.25,   200_000, None),
    (0.50, 1_000_000, 1),
    (1.00, 1_000_000, None),
    (1.50, 1_000_000, 3),
    (1.90,   200_000, None),
    (1.99,   200_000, None),
    (2.00, 1_000_000, None),
    (2.01,   200_000, None),
    (2.50,   500_000, 5),
    (3.00, 1_000_000, None),
    (3.50,   300_000, 7),
    (4.00, 1_000_000, None),
    (5.00, 1_000_000, None),
    (6.00, 1_000_000, None),
]


def run(f, c, p, label):
    m = orbit_order(f, c)
    G = logabsZ(np.arange(m, dtype=float), f, c, p)
    lamB = float(np.mean(G))
    print("=" * 78)
    print(label)
    print("  orbit order m = %d   log|Z_k| = %s" % (m, np.round(G, 9)))
    print("  lambda_B (k_n=n, EXACT cyclic mean) = %.9f" % lamB)
    print("  lambda_A (k_n=1)                    = %.9f" % G[1 % m])
    print()
    print("  alpha        N     lambda(floor(n^a))   dev from lambda_B   residue density mod %d" % m)
    for a, N, half in ROWS:
        km = kmod(N, a, m, half)
        lam = float(np.mean(G[km]))
        hist = np.bincount(km, minlength=m) / N
        print("  %5.2f  %9d  %18.9f  %+16.9f   %s"
              % (a, N, lam, lam - lamB, np.round(hist, 6)))
    print()


run(np.pi, np.pi / 2, P_S3,
    "EXHIBIT 1 (exact)  f = pi, c = pi/2   [S3 sec6(f) row 3, 'finite orbit']")

run(np.pi, 3 * np.pi / 2, P_S3,
    "EXHIBIT 1b (exact) f = pi, c = 3pi/2  [S3 sec6(f) row 4, 'order 4']")

run(3 * np.pi / 2, np.pi / 2, [0.5, 0.125, 0.125, 0.125, 0.125],
    "EXHIBIT 3 (exact)  f = 3pi/2, c = pi/2, p0=1/2 q=r=1/4  -- Z_2 = 0 EXACTLY")

print("=" * 78)
print("CLOSED-FORM PREDICTIONS (proved, not fitted), f=pi c=pi/2, m=4")
G = logabsZ(np.arange(4.0), np.pi, np.pi / 2, P_S3)
print("  a>=2 EVEN : n^a mod 4 = 0 (n even) / 1 (n odd), density (1/2,1/2,0,0)")
print("              -> lambda = (G0+G1)/2 = %.9f" % (0.5 * (G[0] + G[1])))
print("  a>=3 ODD  : n^a mod 4 = 0 (n even) / n mod 4 (n odd), (1/2,1/4,0,1/4)")
print("              -> lambda = G0/2 + (G1+G3)/4 = %.9f"
      % (0.5 * G[0] + 0.25 * G[1] + 0.25 * G[3]))
print("  (G1 = G3 here, so both integer families give the SAME broken value)")
print("  lambda_B = %.9f   lambda_A = %.9f" % (float(np.mean(G)), G[1]))
print("  GAP lambda(n^2) - lambda_B = %+.9f" % (0.5 * (G[0] + G[1]) - float(np.mean(G))))
