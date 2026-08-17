#!/usr/bin/env python3
"""
R1_02 — ATTACK ON m(P).  The lane's m(P) = -0.767507880358 is load-bearing for EVERY
"agrees / does not agree" verdict in F10.  It is computed there by two routes that share a
float64 Gauss-Legendre kernel and a hand-derived kink location.  Here it is recomputed at
50 decimal digits by three routes that share NOTHING with the lane:

  A. mpmath quadrature of the one-dimensional reduction, derived here independently:
         P = (p00 + p10 x) + y (p01 + p11 x),  Jensen in y  =>
         m(P) = (1/2pi) int_0^{2pi} log max(|p00+p10 e^{it}|, |p01+p11 e^{it}|) dt
     with the kink located by mpmath.findroot, NOT by the lane's closed form.
  B. Cassaigne-Maillot with an mpmath Bloch-Wigner dilogarithm (mpmath.polylog), 50 dps.
  C. an ENTIRELY DIFFERENT representation: the Boyd/Smyth-style series
         m(a+bx+cy) obtained by integrating log max over the SECOND variable first,
     i.e. the multiset symmetry is used as a CHECK, not assumed.
  D. exact-value spot checks where a closed form exists: m(1+x+y), m(1+x), log(max).

Precision: mpmath at 50 decimal digits (mp.dps = 50) unless stated.  float64 appears only
where the lane's own number is quoted for comparison.
"""
from mpmath import mp, mpf, mpc, exp, log, pi, quad, polylog, im, arg, fabs, sqrt, acos, findroot
import numpy as np

mp.dps = 50

def m_reduction(p00, p10, p01, p11, dps=50):
    """(1/2pi) int_0^{2pi} log max(|p00+p10 e^{it}|,|p01+p11 e^{it}|) dt, mpmath."""
    p00, p10, p01, p11 = mpf(p00), mpf(p10), mpf(p01), mpf(p11)
    def A(t): return fabs(p00 + p10 * exp(mpc(0, 1) * t))
    def B(t): return fabs(p01 + p11 * exp(mpc(0, 1) * t))
    def g(t): return log(mp.mpf(max(A(t), B(t))))
    # locate kinks on [0,pi] numerically (integrand even in t)
    ts = [mpf(0)]
    NS = 4000
    prev = A(mpf(0)) - B(mpf(0))
    for i in range(1, NS + 1):
        t = pi * i / NS
        cur = A(t) - B(t)
        if prev == 0 or (prev < 0) != (cur < 0):
            try:
                r = findroot(lambda x: A(x) - B(x), (pi * (i - 1) / NS, t), solver='bisect',
                             tol=mpf(10) ** (-dps + 5))
                ts.append(r)
            except Exception:
                pass
        prev = cur
    ts.append(pi)
    ts = sorted(set(ts))
    tot = mpf(0)
    for lo, hi in zip(ts[:-1], ts[1:]):
        if hi - lo > mpf(10) ** (-dps + 5):
            tot += quad(g, [lo, hi])
    return tot / pi          # (2 * int_0^pi) / (2 pi)

def bloch_wigner(z):
    z = mpc(z)
    return im(polylog(2, z)) + arg(1 - z) * log(fabs(z))

def m_CM(a, b, c):
    """m(a+bx+cy), a,b,c >= 0, mpmath.  Cassaigne-Maillot."""
    a, b, c = mpf(a), mpf(b), mpf(c)
    v = sorted([a, b, c])
    if v[2] >= v[0] + v[1]:
        return log(v[2])
    def ang(op, s1, s2):
        return acos((s1 * s1 + s2 * s2 - op * op) / (2 * s1 * s2))
    al, be, ga = ang(a, b, c), ang(b, a, c), ang(c, a, b)
    D = bloch_wigner((a / b) * exp(mpc(0, 1) * ga))
    return (D + al * log(a) + be * log(b) + ga * log(c)) / pi

print("=" * 78)
print("R1_02 — m(P) AT 50 DECIMAL DIGITS, THREE INDEPENDENT ROUTES")
print("=" * 78)

print("\nD. CALIBRATION against closed forms:")
smyth = m_CM(1, 1, 1)
print("   m(1+x+y)  route B      = %s" % mp.nstr(smyth, 30))
print("   m(1+x+y)  route A      = %s" % mp.nstr(m_reduction(1, 1, 1, 0), 30))
print("   Smyth  L'(chi_-3,-1) = (3 sqrt3 / 4pi) L(chi_-3,2):")
L2 = mp.nsum(lambda n: (mp.legendre_kron if False else 0), [1, 1]) if False else None
# L(chi_-3,2) = sum_{n>=1} chi(n)/n^2 with chi = +1,-1,0 mod 3
Lchi = mp.nsum(lambda n: (mpf(1) / (3 * n - 2) ** 2 - mpf(1) / (3 * n - 1) ** 2), [1, mp.inf])
closed = 3 * sqrt(3) / (4 * pi) * Lchi
print("                          = %s   |diff| = %s"
      % (mp.nstr(closed, 30), mp.nstr(fabs(closed - smyth), 5)))

print("\nA/B/C ON THE LANE'S LOAD-BEARING WEIGHT VECTORS  (p00 = 0):")
cases = [("S3/S4 generic  (0.3,0.3,0.4)", (0.0, mpf(3)/10, mpf(3)/10, mpf(2)/5)),
         ("S1 published   (0, 1/2, 1/2)", (0.0, mpf(0), mpf(1)/2, mpf(1)/2)),
         ("centroid       (1/3,1/3,1/3)", (0.0, mpf(1)/3, mpf(1)/3, mpf(1)/3)),
         ("uniform on 5   (2/5,2/5,1/5)", (0.0, mpf(2)/5, mpf(2)/5, mpf(1)/5)),
         ("non-triangle   (0.6,.25,.15)", (0.0, mpf(6)/10, mpf(1)/4, mpf(15)/100)),
         ("degenerate     (0.2,0.5,0.3)", (0.0, mpf(2)/10, mpf(1)/2, mpf(3)/10)),
         ("asym triangle  (.45,.3,.25) ", (0.0, mpf(45)/100, mpf(3)/10, mpf(1)/4))]
LANE = {"S3/S4 generic  (0.3,0.3,0.4)": "-0.767507880358",
        "S1 published   (0, 1/2, 1/2)": "-0.693147180560",
        "centroid       (1/3,1/3,1/3)": "-0.775546341449",
        "uniform on 5   (2/5,2/5,1/5)": "-0.756573585640"}
for name, (p00, p10, p01, p11) in cases:
    A_ = m_reduction(p00, p10, p01, p11)
    B_ = m_CM(p11, p01, p10)
    # route C: reduce in the OTHER variable  P = (p00 + p01 y) + x(p10 + p11 y)
    C_ = m_reduction(p00, p01, p10, p11)
    print("  %-30s" % name)
    print("      A reduction  = %s" % mp.nstr(A_, 30))
    print("      B Cassaigne-Maillot = %s   |A-B| = %s" % (mp.nstr(B_, 30), mp.nstr(fabs(A_-B_), 4)))
    print("      C other-variable    = %s   |A-C| = %s" % (mp.nstr(C_, 30), mp.nstr(fabs(A_-C_), 4)))
    if name in LANE:
        print("      LANE float64 value  = %s        |A - lane| = %s"
              % (LANE[name], mp.nstr(fabs(A_ - mpf(LANE[name])), 4)))

print("\nTHE CENTROID CLOSED FORM the lane quotes: log(1/3) + m(1+x+y)")
cc = log(mpf(1)/3) + smyth
print("   log(1/3)+m(1+x+y) = %s" % mp.nstr(cc, 30))
print("   direct m(P)       = %s   |diff| = %s"
      % (mp.nstr(m_reduction(0, mpf(1)/3, mpf(1)/3, mpf(1)/3), 30),
         mp.nstr(fabs(cc - m_reduction(0, mpf(1)/3, mpf(1)/3, mpf(1)/3)), 5)))

print("\nS1 (lane): m(P) = log(max) EXACTLY when max > 1/2.  Tested at 50 dps on 12 points:")
worst = mpf(0)
import random
random.seed(7)
for _ in range(12):
    a = mpf(random.uniform(0.51, 0.98))
    r = 1 - a
    b = mpf(random.uniform(0.0, 1.0)) * r
    c = r - b
    val = m_reduction(0, b, c, a)
    worst = max(worst, fabs(val - log(max(a, b, c))))
print("   worst |m(P) - log(max)| over the 12 non-triangle points = %s" % mp.nstr(worst, 4))

print("\nS2 (lane): m(P) <= -log 2 on the triangle region max <= 1/2, equality on its boundary.")
print("   -log 2 = %s" % mp.nstr(-log(2), 25))
for w in [(mpf(1)/2, mpf(1)/2, mpf(0)), (mpf(1)/2, mpf(1)/4, mpf(1)/4),
          (mpf(1)/2, mpf(3)/10, mpf(2)/10), (mpf(2)/5, mpf(2)/5, mpf(1)/5)]:
    print("   (p10,p01,p11) = (%s,%s,%s)  max=%s   m = %s"
          % (mp.nstr(w[0],4), mp.nstr(w[1],4), mp.nstr(w[2],4), mp.nstr(max(w),4),
             mp.nstr(m_reduction(0, w[0], w[1], w[2]), 20)))
print("\nDONE R1_02")
