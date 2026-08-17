#!/usr/bin/env python3
"""
M1_06 — EQUIDISTRIBUTION IS NOT ENOUGH.  AN EXPLICIT LIOUVILLE PAIR (u,v) WITH NO
MULTIPLICATIVE RELATION FOR WHICH  (1/N) sum_{k<=N} log|Z_k|  DOES NOT CONVERGE TO m(P).

This is the point the brief asks about and the point W-03 recorded as MISSING from the
IMPORT AUDIT.  Weyl equidistribution gives Birkhoff convergence for RIEMANN-INTEGRABLE
observables.  log|P| is not one: it has logarithmic singularities wherever P vanishes on
T^2, and (M1_02, Z1) P DOES vanish on T^2 exactly when max(p10,p01,p11) <= 1/2.  So the
question "does the average converge to m(P)?" is not answered by Weyl, and the answer is:
NOT FOR EVERY EQUIDISTRIBUTING (u,v).

THE CONSTRUCTION.  Weights p10 = p01 = p11 = 1/3, so P = (x + y + xy)/3 and
m(P) = log(1/3) + m(1+x+y) = -0.7755463414...  The zeros of P on T^2 are the two points
(x,y) = (w, conj w) and (conj w, w) with w = exp(2 pi i/3); in additive coordinates
(alpha,beta) in (R/Z)^2 they are (1/3, 2/3) and (2/3, 1/3).

Put  a_1 = 1,  a_{j+1} = a_j + C_j * 10^{a_j}  (so a_{j+1} > 2 a_j), k_j = 10^{a_j}, and
     alpha = 1/3 + sum_j d_j 10^{-a_j},     beta = 2/3 + sum_j c_j 10^{-a_j},
with (d_j, c_j) alternating between (1,1) and (1,2).

(i)  NO MULTIPLICATIVE RELATION.  m*alpha + n*beta in Z forces m*A + n*B in Q where
     A = sum d_j 10^{-a_j}, B = sum c_j 10^{-a_j}.  m*A + n*B = sum (m d_j + n c_j)10^{-a_j}
     is a sparse series with bounded coefficients and gaps a_{j+1} - a_j -> infinity, hence
     (Liouville) irrational unless m d_j + n c_j = 0 for all large j.  The alternation forces
     m + n = 0 and m + 2n = 0, hence m = n = 0.  So L = {0}: the orbit IS equidistributed in
     T^2 by Weyl, and every hypothesis in the "orbit dense in T^2" clause holds.

(ii) 10^a = 1 mod 3 for every a >= 0, so k_j * (1/3) = 1/3 and k_j * (2/3) = 2/3 mod 1, and
     frac(k_j alpha) = 1/3 + delta_j,  frac(k_j beta) = 2/3 + eta_j  with
     delta_j = d_{j+1} 10^{-(a_{j+1}-a_j)} (1+o(1)),  eta_j likewise.  The orbit point at
     k = k_j sits within ~10^{-(a_{j+1}-a_j)} of a ZERO of P.

(iii) Near a simple zero |P| is comparable to the distance, so
        log|Z_{k_j}| = -(a_{j+1} - a_j) log 10 + O(1) = -C_j k_j log 10 + O(1),
      and therefore
        (1/k_j) sum_{k<=k_j} log|Z_k|  <=  -C_j log 10 + O(1).
      Choosing C_j -> infinity gives liminf_N (1/N) sum log|Z_k| = -infinity, while
      m(P) = -0.7755.  Choosing C_j = C constant gives a liminf below m(P) - C log 10.

WHAT IS VERIFIED NUMERICALLY BELOW (float64 for the arithmetic, EXACT Fractions for the
phases): the first dip, with a_2 - a_1 = 300, computed from exact rational phase reduction
and a validated local expansion of P.  The deeper dips are the same computation with larger
gaps and are stated, not run: their offsets are 10^(-10^301) and no machine holds them.
"""
import numpy as np
from fractions import Fraction
from M1_02_mahler_machinery import m_R1

W = np.exp(2j * np.pi / 3.0)
X0, Y0 = W, np.conj(W)                 # a zero of (x+y+xy)/3 on T^2
ZEROS = [(Fraction(1, 3), Fraction(2, 3)), (Fraction(2, 3), Fraction(1, 3))]


def P_direct(fx, fy):
    x = np.exp(2j * np.pi * float(fx))
    y = np.exp(2j * np.pi * float(fy))
    return (x + y + x * y) / 3.0


def P_local(zero_idx, dx, dy):
    """|P| near a zero, from the exact first-order expansion.
       P(x0 e^{2 pi i dx}, y0 e^{2 pi i dy}) = (2 pi i/3)[dx (x0 + x0 y0) + dy (y0 + x0 y0)]
       + O(d^2).  dx, dy are float64 (they may be as small as 1e-300)."""
    x0 = np.exp(2j * np.pi * float(ZEROS[zero_idx][0]))
    y0 = np.exp(2j * np.pi * float(ZEROS[zero_idx][1]))
    return (2j * np.pi / 3.0) * (dx * (x0 + x0 * y0) + dy * (y0 + x0 * y0))


print("=" * 78)
print("M1_06 — EQUIDISTRIBUTION IS NOT SUFFICIENT.  A LIOUVILLE COUNTEREXAMPLE.")
print("=" * 78)
mP = m_R1(0.0, 1 / 3, 1 / 3, 1 / 3)
print("\nweights (1/3,1/3,1/3):  m(P) = %.12f" % mP)

# ---- validate the local expansion against direct evaluation at moderate offsets
print("\nVALIDATION of the local expansion (must agree with direct evaluation):")
for d in (1e-4, 1e-6, 1e-8):
    fx = ZEROS[0][0] + Fraction(d).limit_denominator(10 ** 18)
    fy = ZEROS[0][1] + Fraction(2 * d).limit_denominator(10 ** 18)
    direct = abs(P_direct(fx, fy))
    loc = abs(P_local(0, d, 2 * d))
    print("    offset %8.1e :  direct |P| = %.12e   expansion = %.12e   rel dev %.2e"
          % (d, direct, loc, abs(direct - loc) / direct))

# ---- the construction, exactly
A1, A2 = 1, 301                    # a_1, a_2   (a_3 = 301 + C*10^301, symbolic)
d1, d2 = 1, 1
c1, c2 = 1, 2
alpha = Fraction(1, 3) + Fraction(d1, 10 ** A1) + Fraction(d2, 10 ** A2)
beta = Fraction(2, 3) + Fraction(c1, 10 ** A1) + Fraction(c2, 10 ** A2)
print("\nalpha = 1/3 + 1e-1 + 1e-301   (plus a tail below 10^-(10^301), immaterial to k<=10)")
print("beta  = 2/3 + 1e-1 + 2e-301")

K = 10
logs = []
print("\n   k   frac(k alpha)-1/3      frac(k beta)-2/3      |Z_k|            log|Z_k|")
for k in range(1, K + 1):
    fx = (k * alpha) % 1
    fy = (k * beta) % 1
    # distance to each zero, exactly, as a Fraction on the circle
    best = None
    for zi, (gx, gy) in enumerate(ZEROS):
        dx = (fx - gx) % 1
        dx = dx if dx <= Fraction(1, 2) else dx - 1
        dy = (fy - gy) % 1
        dy = dy if dy <= Fraction(1, 2) else dy - 1
        d = abs(dx) + abs(dy)
        if best is None or d < best[0]:
            best = (d, zi, dx, dy)
    dist, zi, dx, dy = best
    if dist < Fraction(1, 10 ** 6):
        az = abs(P_local(zi, float(dx), float(dy)))
        tag = "  <- local expansion (offset %.1e)" % float(dist)
    else:
        az = abs(P_direct(fx, fy))
        tag = ""
    logs.append(np.log(az))
    print("  %2d   %+18.10e  %+18.10e   %12.5e   %12.4f%s"
          % (k, float((fx - Fraction(1, 3)) % 1 if (fx - Fraction(1, 3)) % 1 <= Fraction(1, 2)
                      else ((fx - Fraction(1, 3)) % 1) - 1),
             float((fy - Fraction(2, 3)) % 1 if (fy - Fraction(2, 3)) % 1 <= Fraction(1, 2)
                   else ((fy - Fraction(2, 3)) % 1) - 1),
             az, np.log(az), tag))

run = float(np.mean(logs))
print("\n  (1/10) sum_{k<=10} log|Z_k| = %.6f        m(P) = %.6f" % (run, mP))
print("  DIP BELOW m(P) = %.4f    [predicted -(a_2-a_1) log 10 / k_1 = %.4f]"
      % (run - mP, -(A2 - A1) * np.log(10) / 10.0))
print("\n  This is one dip.  The same construction with a_{j+1} - a_j = C_j 10^{a_j} puts a dip")
print("  of depth ~ C_j log 10 at N = 10^{a_j} for every j, so with C_j -> infinity")
print("      liminf_N (1/N) sum_{k<=N} log|Z_k| = -infinity   while   m(P) = -0.7755463.")
print("  The pair (u,v) = (e^{2 pi i alpha}, e^{2 pi i beta}) has NO multiplicative relation,")
print("  so its orbit IS dense and equidistributed in T^2.  EQUIDISTRIBUTION IS NOT ENOUGH.")

# ---- the other side: what a Diophantine connection actually does
print("\n" + "-" * 78)
print("CONTRAST: how close does a DIOPHANTINE orbit get to the zeros?  alpha = -2^(1/3),")
print("beta = 4^(1/3) (Schmidt), weights (1/3,1/3,1/3), k <= 10^7, exact int64 phase reduction.")
al = -(2.0 ** (1 / 3)) % 1.0
be = (4.0 ** (1 / 3)) % 1.0
D = 2 ** 39
An = int(np.floor(al * D)); Bn = int(np.floor(be * D))
dA = al - An / D; dB = be - Bn / D
print("        N        min_{k<=N} |Z_k|     min dist to a zero    sqrt(N)*dist")
prev = 0
mn = np.inf
mnd = np.inf
for N in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    k = np.arange(prev + 1, N + 1, dtype=np.int64)
    fa = np.mod(((k * An) % D) / D + k * dA, 1.0)
    fb = np.mod(((k * Bn) % D) / D + k * dB, 1.0)
    x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
    az = np.abs(x + y + x * y) / 3.0
    mn = min(mn, float(az.min()))
    for (gx, gy) in ZEROS:
        ddx = np.abs(((fa - float(gx) + 0.5) % 1) - 0.5)
        ddy = np.abs(((fb - float(gy) + 0.5) % 1) - 0.5)
        mnd = min(mnd, float(np.min(np.hypot(ddx, ddy))))
    prev = N
    print("   %10d   %16.6e   %16.6e   %10.4f" % (N, mn, mnd, np.sqrt(N) * mnd))
print("   The closest approach shrinks like N^{-1/2} — the spacing of N points in T^2 — so")
print("   its contribution to the average is O(log N / N) and dies.  The Liouville pair above")
print("   makes that same quantity e^{-C N} instead, and the contribution is O(1).")
print("   THE HYPOTHESIS THAT DECIDES IT is a lower bound on dist((u^k,v^k), Z(P)), i.e. an")
print("   INHOMOGENEOUS DIOPHANTINE CONDITION relative to the two zeros of P — not merely")
print("   the homogeneous one that gives equidistribution.")
print("\nDONE M1_06")
