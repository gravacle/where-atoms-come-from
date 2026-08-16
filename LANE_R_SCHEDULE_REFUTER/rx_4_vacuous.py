"""
R-4: the claim's own verification instrument, tested for the failure it was
supposed to be able to detect.

The claim says: "The Abel limit verified for beta = 1/alpha - 1 = -0.9, -0.5,
0, 1, 3, 9 (i.e. alpha = 10, 2, 1, 0.5, 0.25, 0.1): all -> -0.767508 at M=1e7."

The Abel object is   L(beta,M) = sum_{m<=M} m^beta G(m w) / sum_{m<=M} m^beta,
with G(mw) = log|Z_m|.  It sums over EVERY integer m <= M with a positive
smooth weight.  For alpha <= 1 that is a faithful re-description of the
schedule (floor(n^alpha) hits every m, with multiplicity c_m ~ q m^{q-1}).
For alpha > 1 the schedule floor(n^alpha) is a strictly increasing SUBSEQUENCE
that SKIPS integers -- and no positive-weight sum over all m can see a skip.

So for alpha > 1 the Abel check is VACUOUS BY CONSTRUCTION: it returns the
full-orbit average whatever the schedule does.  Demonstrated below on the
connection where the schedule provably does something else.
"""
import numpy as np
from rx_lib import logabsZ, orbit_order

P_S3 = [0.4, 0.15, 0.15, 0.15, 0.15]
f, c = np.pi, np.pi / 2
m = orbit_order(f, c)
G = logabsZ(np.arange(m, dtype=float), f, c, P_S3)
lamB = float(np.mean(G))

print("connection f=pi c=pi/2, orbit order m=%d, lambda_B = %.9f" % (m, lamB))
print()
print("  alpha  beta=1/alpha-1   Abel L(beta,M=1e7)   TRUE lambda(floor(n^alpha))")
M = 10**7
mm = np.arange(1, M + 1, dtype=np.int64)
Gm = G[mm % m]
for alpha in [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 10.0]:
    beta = 1.0 / alpha - 1.0
    w = mm.astype(float) ** beta
    L = float(np.dot(w, Gm) / w.sum())
    a = int(alpha)
    if float(alpha).is_integer():
        km = np.array([pow(int(n), a, m) for n in range(1, 1000001)])
        true = float(np.mean(G[km]))
    else:
        n = np.arange(1, 10**7 + 1, dtype=float)
        km = (np.floor(n ** alpha).astype(np.int64)) % m
        true = float(np.mean(G[km]))
    print("  %5.2f  %+12.4f   %18.9f   %18.9f%s"
          % (alpha, beta, L, true, "   <== MISMATCH" if abs(L - true) > 1e-6 else ""))

print()
print("READ IT:  the Abel column is flat at lambda_B for every alpha, including")
print("the alphas where the actual schedule is nowhere near lambda_B.  The")
print("instrument returns the answer the claim wanted for reasons that have")
print("nothing to do with the schedule.  It is the S4 Control-1 defect again:")
print("a control that could not have failed.")

print()
print("=" * 78)
print("EXACT ARITHMETIC on the amplifier (fractions, no floating point)")
print("=" * 78)
from fractions import Fraction as F
# u = v = i ; p0 = 1/2, q = r = 1/4.  Z_k = (uv)^k p0 + u^k q + v^k r.
# i^k cycles 1, i, -1, -i.  Work in Z[i] with rational coefficients.
I = [(F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1))]   # i^k, k mod 4
p0, qq, rr = F(1, 2), F(1, 4), F(1, 4)
for k in range(4):
    uv = I[(2 * k) % 4]          # (uv)^k = (i*i)^k = (-1)^k
    uk = I[k % 4]
    vk = I[k % 4]
    re = uv[0] * p0 + uk[0] * qq + vk[0] * rr
    im = uv[1] * p0 + uk[1] * qq + vk[1] * rr
    print("  Z_%d = %s + %s i     |Z_%d|^2 = %s" % (k, re, im, k, re * re + im * im))
print()
print("  Z_2 = 0 EXACTLY (= 2*p0 - 1 with p0 = 1/2).  Therefore:")
print("    schedule k_n = n     : Omega_N = 0 EXACTLY for all N >= 2,  lambda_B = -infinity")
print("    schedule k_n = n^2   : n^2 mod 4 in {0,1} only, so Z=0 is NEVER reached;")
print("                           lambda = (log 1 + log 2^{-1/2})/2 = -(1/4) log 2")
print("                                  = %.12f" % (-0.25 * np.log(2)))
print("  The gap between the two is INFINITE, not numerical.")
