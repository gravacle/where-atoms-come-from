#!/usr/bin/env python3
"""
R1_05 — THE RANK-1 HYPOTHESIS GAP, EXHIBITED, AND THE BOYD-LAWTON ACCUMULATION AUDITED.

PART A.  THE (1,1) RESONANCE  u v = 1  (i.e. conj(W_F) W_C = 1, i.e. W_C = W_F: the flat
holonomy equal to the curvature holonomy).  This is a codimension-one family of K1
connections and row 3 of the lane's own M1_07 table.  There
    P|_H = p11 + p10 z + p01 z^{-1},   |Z_k| = |0.4 + 0.6 cos(k theta)|   (u = e^{i theta}),
which VANISHES whenever cos(k theta) = -2/3.  Q = 0.3 z^2 + 0.4 z + 0.3 is self-reciprocal;
both roots lie EXACTLY on |z| = 1 (proved exactly in R1_04 via gcd(Q, Q*) of degree 2).  So on
this rank-1 locus log|Z_k| is UNBOUNDED and its Birkhoff average is subject to exactly the
inhomogeneous-Diophantine failure that T4 confines to H = T^2.

Exhibited: two connections with the IDENTICAL relation lattice L = Z(1,1) and the identical
weights, one converging to m(Q) = log(0.3) and one driven arbitrarily far below it.

PART B.  BOYD-LAWTON ACCUMULATION: the noise floor of the two routes the lane used.
The lane's np.roots branch caps at degree 1200; above that it uses a 2^23-point quadrature.
Here the quadrature is run at several M and its own convergence measured, so the reported
+4.716e-09 at (2584,1597) can be compared with the method's resolution rather than trusted.

Precision: mpmath at the dps stated; float64 where labelled; exact Fractions for phases.
"""
import numpy as np
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, polyroots, fabs, log as mlog, exp as mexp, pi as mpi, cos as mcos

mp.dps = 60
P10, P01, P11 = Fr(3, 10), Fr(3, 10), Fr(2, 5)

print("=" * 78)
print("R1_05 — PART A: TWO CONNECTIONS ON THE SAME RANK-1 LOCUS, DIFFERENT LIMITS")
print("=" * 78)
tstar = mp.acos(-mpf(2) / 3) / (2 * mpi)
print("\nQ(z) = 0.3 z^2 + 0.4 z + 0.3;  roots EXACTLY on |z|=1 at additive phase +- %s turns"
      % mp.nstr(tstar, 25))
print("m(Q) = log(0.3) = %s      [S4 row (1,1) = -1.203972804]" % mp.nstr(mlog(mpf(3) / 10), 20))

def avg_from_rho(rho, N):
    """(1/N) sum_{k<=N} log|0.4+0.6 cos(2 pi k rho)|; phase reduced EXACTLY in big ints,
    only cos/log in float64."""
    num, den = rho.numerator, rho.denominator
    tot = 0.0; ph = 0
    for k in range(1, N + 1):
        ph = (ph + num) % den
        tot += float(np.log(abs(0.4 + 0.6 * np.cos(2 * np.pi * (ph / den)))))
    return tot / N

rho1 = Fr(int(mp.nstr(1 / (20 * mpi) * mpf(10) ** 60, 70).split('.')[0]), 10 ** 60)
print("\nA1  BENIGN POINT OF THE LOCUS.  rho = 1/(20 pi) (the erratum's own rotation number),")
print("    as a 60-digit rational.  Weights and relation lattice as above.")
for N in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
    print("      N = %-9d  (1/N) sum log|Z_k| = %.12f" % (N, avg_from_rho(rho1, N)))
print("      -> converging on log(0.3) = %.12f" % float(np.log(0.3)))

DIG = 220
mp.dps = 400
tstar_hi = mp.acos(-mpf(2) / 3) / (2 * mpi)
T = Fr(int(mp.nstr(tstar_hi * mpf(10) ** DIG, DIG + 20).split('.')[0]), 10 ** DIG)
rho2 = (Fr(3) + T) / 10
tenrho = (10 * rho2) % 1
off = fabs(mpf(tenrho.numerator) / tenrho.denominator - tstar_hi)
print("\nA2  LIOUVILLE POINT OF THE SAME LOCUS.  rho = (3 + T)/10, T = tstar truncated at 1e-%d." % DIG)
print("    (rho becomes irrational on appending M1_06's own sparse tail, which is below")
print("     10^-(10^220) and cannot move k <= 10.)")
print("    distance of the k=10 orbit point to the singular phase = %s turns" % mp.nstr(off, 6))
tot = mpf(0)
for k in range(1, 11):
    f = (k * rho2) % 1
    val = fabs(mpf(2) / 5 + mpf(3) / 5 * mcos(2 * mpi * mpf(f.numerator) / f.denominator))
    if k == 10:
        print("    |Z_10| = %s   log|Z_10| = %s" % (mp.nstr(val, 8), mp.nstr(mlog(val), 12)))
    tot += mlog(val)
mQ = mlog(mpf(3) / 10)
print("    (1/10) sum_{k<=10} log|Z_k| = %s   against m(Q) = %s"
      % (mp.nstr(tot / 10, 14), mp.nstr(mQ, 14)))
print("    DIP BELOW m(Q) = %s" % mp.nstr(tot / 10 - mQ, 10))
mp.dps = 60
print("""
  CONSEQUENCE.  Two connections with IDENTICAL relation lattice L = Z(1,1) -- the lane's
  NAMED OPERATIVE VARIABLE -- and identical weights have different limits; iterating with
  growing gaps drives the second to liminf = -infinity.  Therefore:
   (i)  T4's sentence "What is fragile is only the NAMING of the exponential rate as m(P):
        that needs H = T^2 plus the inhomogeneous Diophantine condition of T2(c)" is FALSE.
        The same fragility lives on a rank-1 locus, on row 3 of the lane's own table.
   (ii) F6/M1_07's 16-row table is presented as a table of LIMITS.  Row (1,1) is not a
        limit; it is a Mahler measure that the orbit average need not converge to.
   (iii)The other 15 rows ARE limits, and freely so -- R1_04 proves exactly, over Q, that
        their Q has no root on |z| = 1, so log|Q| is continuous there and Weyl alone
        suffices.  The lane gives no such argument for any of them.""")

# ------------------------------------------------------------------------------
print("\n" + "=" * 78)
print("R1_05 — PART B: THE ACCUMULATION'S NOISE FLOOR")
print("=" * 78)
MP = mpf("-0.767507880357775871645874051819")

def quad_m(mm, nn, MQ):
    th = np.arange(MQ) * (2 * np.pi / MQ)
    val = np.zeros(MQ, dtype=complex)
    for e, cc in [(nn, 0.3), (-mm, 0.3), (nn - mm, 0.4)]:
        val += cc * np.exp(1j * e * th)
    return float(np.mean(np.log(np.maximum(np.abs(val), 1e-300))))

print("\n  (m,n)          M=2^20            M=2^22            M=2^24        spread    dev vs m(P)")
for (mm, nn) in [(11, 20), (97, 61), (610, 377), (2584, 1597)]:
    vals = [quad_m(mm, nn, 1 << e) for e in (20, 22, 24)]
    print("  (%5d,%5d)  %.12f  %.12f  %.12f  %.1e   %+.4e"
          % (mm, nn, vals[0], vals[1], vals[2], max(vals) - min(vals), vals[2] - float(MP)))
print("\n  lane reported deviations: (11,20) +4.929e-04 ; (97,61) -7.687e-06 ;")
print("                            (610,377) +5.251e-07 ; (2584,1597) +4.716e-09")
print("  -> the M-spread above is the quadrature's own resolution.  Read the last row's")
print("     +4.716e-09 against it before treating it as a measured convergence rate.")
print("\nDONE R1_05")
