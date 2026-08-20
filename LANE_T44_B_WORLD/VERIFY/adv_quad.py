"""ADV-FLOAT independent instrument: continuum Fourier quadrature for the critical
lattice kernel on Z^3 (an instrument the lane was FORBIDDEN to use -- adversary-side
only), plus independent recomputation of the D=2 / D=1 discriminators.

a(x) = G(0) - G(x) = (1/(2pi)^3) int (1 - cos(k.x)) / (1 - (cos k1 + cos k2 + cos k3)/3)
The integrand is bounded (removable at 0); midpoint tensor grid, two resolutions for an
empirical error estimate.  These floats DECIDE NOTHING alone -- they cross-check the
lane's exact brackets with an entirely different method.
"""
import numpy as np
from math import comb, pi, log2

def a_quad(targets, N):
    th = (np.arange(N) + 0.5) / N * 2 * pi - pi   # midpoint grid avoiding 0
    c1 = np.cos(th)[:, None, None]
    c2 = np.cos(th)[None, :, None]
    c3 = np.cos(th)[None, None, :]
    denom = 1.0 - (c1 + c2 + c3) / 3.0
    out = {}
    k1 = th[:, None, None]; k2 = th[None, :, None]; k3 = th[None, None, :]
    for t in targets:
        num = 1.0 - np.cos(t[0] * k1 + t[1] * k2 + t[2] * k3)
        out[t] = float(np.mean(num / denom))
    return out

TARGETS = [(2,0,0),(4,0,0),(6,0,0),(8,0,0),(12,0,0),(16,0,0),
           (2,2,0),(4,4,0),(8,8,0),(2,2,2),(4,4,4),(8,8,8)]
A1 = a_quad(TARGETS, 160)
A2 = a_quad(TARGETS, 240)
print("quadrature a(x), N=160 vs N=240 (drift = empirical error):")
for t in TARGETS:
    print("  a%s = %.6f  (drift %.2e)" % (t, A2[t], abs(A2[t] - A1[t])))
H = lambda t1, t2: A2[t2] - A2[t1]
print("increment ratios (independent continuum instrument):")
print("  axis 2->4 : %.6f" % (H((4,0,0),(8,0,0)) / H((2,0,0),(4,0,0))))
print("  axis 4->8 : %.6f" % (H((8,0,0),(16,0,0)) / H((4,0,0),(8,0,0))))
print("  fdiag 2->4: %.6f" % (H((4,4,0),(8,8,0)) / H((2,2,0),(4,4,0))))
print("  bdiag 2->4: %.6f" % (H((4,4,4),(8,8,8)) / H((2,2,2),(4,4,4))))
print("  exponent -log2 axis 4->8: %.4f" % (-log2(H((8,0,0),(16,0,0)) / H((4,0,0),(8,0,0)))))
print("  coefficient 2*8*H(8->16): %.6f   (3/(2pi) = %.6f)"
      % (16 * H((8,0,0),(16,0,0)), 3 / (2 * pi)))
WATSON = 1.5163860591519780
print("  G(2) = Watson - a(2) = %.6f ; d*G: 2G(2)=%.4f 4G(4)=%.4f 8G(8)=%.4f 16G(16)=%.4f"
      % (WATSON - A2[(2,0,0)],
         2 * (WATSON - A2[(2,0,0)]), 4 * (WATSON - A2[(4,0,0)]),
         8 * (WATSON - A2[(8,0,0)]), 16 * (WATSON - A2[(16,0,0)])))

# ---------------- D=2 discriminator, independent (rotation bijection, no lane code) ----
def a2d(d, M):
    """a_2D(d,0) = sum_{m<=M} (C(2m,m)^2 - C(2m,m+d/2)^2)/16^m, float log-space."""
    from math import lgamma, exp
    lf = [0.0]
    for i in range(1, 2 * M + 2): lf.append(lf[-1] + np.log(i))
    s = 0.0
    h = d // 2
    for m in range(M + 1):
        l0 = 2 * (lf[2*m] - 2*lf[m]) - 2*m*np.log(16.0)/1.0
        lc = lf[2*m] - lf[m] - lf[m]
        p0 = exp(2*lc - 2*m*np.log(4.0))
        pt = exp(2*(lf[2*m] - lf[m+h] - lf[m-h]) - 2*m*np.log(4.0)) if m >= h else 0.0
        s += p0 - pt
    return s
a2 = {d: a2d(d, 12000) for d in (2, 4, 8, 16)}
i1 = a2[4] - a2[2]; i2 = a2[8] - a2[4]; i3 = a2[16] - a2[8]
print("D=2 independent (M=12000): increments %.6f %.6f %.6f ratios %.4f %.4f  (LOG => ~1)"
      % (i1, i2, i3, i2 / i1, i3 / i2))

# ---------------- D=1 discriminator, independent ---------------------------------------
def a1d(d, M):
    from math import exp
    lf = [0.0]
    for i in range(1, 2 * M + 2): lf.append(lf[-1] + np.log(i))
    s = 0.0
    h = d // 2
    for m in range(M + 1):
        p0 = exp(lf[2*m] - 2*lf[m] - 2*m*np.log(2.0))
        pt = exp(lf[2*m] - lf[m+h] - lf[m-h] - 2*m*np.log(2.0)) if m >= h else 0.0
        s += p0 - pt
    return s
a1 = {d: a1d(d, 200000) for d in (2, 4, 8)}
j1 = a1[4] - a1[2]; j2 = a1[8] - a1[4]
print("D=1 independent (M=200000): increments %.5f %.5f ratio %.4f  (LIN => ~2; exact a(d)=d)"
      % (j1, j2, j2 / j1))

# ---------------- D=2 venue at mu=1/6 (subcritical there) ------------------------------
from fractions import Fraction
def series2d(mu, a, b, K):
    mu = Fraction(mu); S = Fraction(0)
    for k in range(K + 1):
        if (k + a + b) % 2 == 0 and k >= abs(a) + abs(b):
            n = comb(k, (k + a + b)//2) * comb(k, (k + a - b)//2)
            S += n * mu ** k
    return S
S8 = series2d(Fraction(1,6), 8, 0, 260)
S10 = series2d(Fraction(1,6), 10, 0, 260)
print("D=2 at mu=1/6 (exact partial sums, K=260): r = %.6f  (published 0.063870)"
      % float(S10 / S8))
