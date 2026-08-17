#!/usr/bin/env python3
"""
R2_05 — FINDING F3's NUMERICAL EVIDENCE IS A COMPARISON OF AN ORBIT AVERAGE AGAINST A
QUADRATURE ERROR, AND ONE OF ITS FOUR ROWS IS VACUOUS.

THE TARGET'S F3 EVIDENCE (M1_03 truncation block, Diophantine pair, N=1e6, (0.3,0.3,0.4)):
    eps=1e-1  orbit -0.758503420  vs torus mean -0.758503816   (4.0e-07)
    eps=1e-2  orbit -0.767418102  vs torus mean -0.767418907   (8.1e-07)
    eps=1e-3  orbit -0.767506494  vs torus mean -0.767506966   (4.7e-07)
    eps=1e-4  orbit -0.767507874  vs torus mean -0.767508002   (1.3e-07)
The "torus mean" column is a 3000 x 3000 midpoint grid.  TWO DEFECTS.

DEFECT 1 (an inequality violation that is a PROOF the grid column is wrong).
    f_eps = max(log|P|, log eps) >= log|P| POINTWISE, so int f_eps >= m(P) for EVERY eps.
    m(P) = -0.767507880358.  The eps=1e-4 grid value -0.767508002 is BELOW it.  Impossible.
    So that row's 1.3e-07 "deviation" is the grid's error, not the orbit's.

DEFECT 2 (a vacuous row).  On the Diophantine orbit min_{k<=1e6}|Z_k| = 4.9415e-04 (R2_03),
    which is ABOVE 1e-4.  So at eps = 1e-4 the truncation NEVER FIRES on the orbit and the
    "truncated orbit average" is bit-for-bit the untruncated one (M1_03 case A at N=1e6 is
    -0.767507873838; the eps=1e-4 row prints -0.767507874).  That row tests nothing about
    truncation.  "Could not have failed" voids a CONTROL -- and this is a control.

WHAT I DO INSTEAD.  int f_eps is computed to ~1e-12 in CLOSED FORM plus a graded quadrature:
    |P| = |a + b(th) w| with a = p10, b(th) = |p01 + p11 e^{i th}| and w uniform, so
        int f_eps = m(P) + (1/pi) int_0^pi J(a,b(th),eps) dth,
        J(a,b,eps) = (1/pi) int_0^{s0} ( log eps - (1/2) log( (a-b)^2 + 4ab sin^2(s/2) ) ) ds,
        s0 = 2 arcsin( sqrt( (eps^2-(a-b)^2) / (4ab) ) ),  J = 0 unless |a-b| < eps.
    Both integrals use a GEOMETRICALLY GRADED mesh (60 dyadic levels) with 40-point
    Gauss-Legendre per level, which resolves the log singularity at a = b.  Convergence is
    demonstrated by doubling the mesh.
ISOLATION LEDGER: HELD FIXED -- P, eps, the definition of f_eps.  THE ONE THING THAT MOVES:
    the method used to evaluate int f_eps (3000^2 midpoint grid -> closed form + graded
    quadrature).  If the grid column were accurate the two would agree to 1e-12.
Precision: float64; the closed-form reduction removes the 2-D grid entirely.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W08_M1_IDENTIFICATION")
from M1_02_mahler_machinery import m_R1

P10, P01, P11 = 0.3, 0.3, 0.4
A = P10
mP = m_R1(0.0, P10, P01, P11)
GLN = 40
GX, GW = np.polynomial.legendre.leggauss(GLN)


def graded(f, lo, hi, levels=60, at_lo_singular=True):
    """integrate f on [lo,hi] with a dyadic mesh graded towards lo (or hi)."""
    tot = 0.0
    a, b = lo, hi
    for _ in range(levels):
        mid = 0.5 * (a + b) if at_lo_singular else 0.5 * (a + b)
        if at_lo_singular:
            x = 0.5 * (mid - a) * GX + 0.5 * (mid + a)
            tot += float(np.sum(GW * f(x))) * 0.5 * (mid - a)
            a = mid
        else:
            x = 0.5 * (b - mid) * GX + 0.5 * (b + mid)
            tot += float(np.sum(GW * f(x))) * 0.5 * (b - mid)
            b = mid
    return tot


def bth(th):
    return np.sqrt(P01 ** 2 + P11 ** 2 + 2 * P01 * P11 * np.cos(th))


def J(b, eps):
    """(1/pi) int_0^{s0} ( log eps - 0.5 log( D^2 + 4ab sin^2(s/2) ) ) ds ,  D = a-b."""
    D = A - b
    if abs(D) >= eps:
        return 0.0
    q = (eps ** 2 - D ** 2) / (4.0 * A * b)
    q = min(1.0, max(0.0, q))
    s0 = 2.0 * np.arcsin(np.sqrt(q))

    def g(s):
        return np.log(eps) - 0.5 * np.log(D * D + 4.0 * A * b * np.sin(s / 2.0) ** 2)
    return graded(g, 0.0, s0, levels=70, at_lo_singular=True) / np.pi


def int_feps(eps, levels=70):
    """m(P) + (1/pi) int_0^pi J(b(th),eps) dth, integrated over the band where |a-b|<eps."""
    lo_c = ((A - eps) ** 2 - (P01 ** 2 + P11 ** 2)) / (2 * P01 * P11)
    hi_c = ((A + eps) ** 2 - (P01 ** 2 + P11 ** 2)) / (2 * P01 * P11)
    lo_c = min(1.0, max(-1.0, lo_c)); hi_c = min(1.0, max(-1.0, hi_c))
    th_hi = np.arccos(lo_c)          # cos decreasing on [0,pi]
    th_lo = np.arccos(hi_c)
    thstar = np.arccos((A ** 2 - P01 ** 2 - P11 ** 2) / (2 * P01 * P11))

    def h(thv):
        return np.array([J(bth(t), eps) for t in np.atleast_1d(thv)])
    # graded towards thstar from both sides (J has a log-type kink there)
    left = graded(h, th_lo, thstar, levels=levels, at_lo_singular=False)
    right = graded(h, thstar, th_hi, levels=levels, at_lo_singular=True)
    return mP + (left + right) / np.pi


print("=" * 78)
print("R2_05 — int f_eps, DONE PROPERLY, AGAINST M1_03's 3000x3000 GRID")
print("=" * 78)
print("\nm(P) = %.12f    (the LOWER BOUND every int f_eps must respect)" % mP)
GRID = {1e-1: -0.758503816, 1e-2: -0.767418907, 1e-3: -0.767506966, 1e-4: -0.767508002}
ORB = {1e-1: -0.758503420, 1e-2: -0.767418102, 1e-3: -0.767506494, 1e-4: -0.767507874}
print("\n    eps      M1_03 grid value    this lane (closed form +      grid error     >= m(P)?")
print("                                   graded quad, 2 mesh sizes)")
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    v1 = int_feps(eps, levels=60)
    v2 = int_feps(eps, levels=90)
    ok = "YES" if GRID[eps] >= mP else "NO  <-- IMPOSSIBLE"
    print("   %6.0e   %16.9f    %16.9f (%.1e)   %+.3e   %s"
          % (eps, GRID[eps], v2, abs(v2 - v1), GRID[eps] - v2, ok))
print("\n  So the 'torus mean' column carries errors of 4e-7 .. 1e-6, which is the SAME SIZE as")
print("  the deviations M1_03 reports as evidence that the orbit average converges.  The")
print("  reported agreement measures the GRID, not the orbit.")

print("\n  DEFECT 2, checked: does the eps=1e-4 truncation ever fire on the orbit at N=1e6?")
A62 = 3413011746732233848; B62 = 2708909218571285002; TWO62 = 1 << 62
SH = 31; MK = (1 << SH) - 1
k = np.arange(1, 10 ** 6 + 1, dtype=np.int64)


def ph(Anum):
    hi = ((k * (Anum >> SH)) % (1 << SH)) << SH
    lo = k * (Anum & MK)
    val = (hi + lo) % TWO62
    return (val >> SH).astype(np.float64) / (1 << SH) + (val & MK).astype(np.float64) / float(TWO62)


x = np.exp(2j * np.pi * ph(A62)); y = np.exp(2j * np.pi * ph(B62))
az = np.abs(P10 * x + P01 * y + P11 * x * y)
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    n = int(np.sum(az < eps))
    print("     eps = %6.0e :  #{k <= 1e6 : |Z_k| < eps} = %d" % (eps, n))
print("     min|Z_k| over k <= 1e6 = %.6e" % float(az.min()))
plain = float(np.mean(np.log(az)))
trunc = float(np.mean(np.maximum(np.log(az), np.log(1e-4))))
print("     untruncated orbit average          = %.12f" % plain)
print("     eps=1e-4 truncated orbit average   = %.12f   difference = %.1e"
      % (trunc, abs(plain - trunc)))
print("     -> the eps=1e-4 row of F3's evidence is the untruncated average relabelled.")
print("        It could not have failed.  VOIDS THAT ROW AS A CONTROL (not the theorem F3,")
print("        which is proved from monotone convergence and stands).")
print("\nDONE R2_05")
