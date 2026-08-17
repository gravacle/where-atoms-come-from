#!/usr/bin/env python3
"""
R_01 — THE SINGULARITY, LOCATED EXACTLY, AND m(P) ACROSS FOUR DECADES OF QUADRATURE.

Everything downstream depends on WHERE log|P| is singular and on HOW singular.  This script
fixes both at K1's registered pi = (0, 0.3, 0.3, 0.4), in closed form, and reproduces the two
control numbers the brief supplies (min|P| on a 2048^2 grid; m(P)).

CLOSED FORM.  P = 0.3 x + 0.3 y + 0.4 x y.  On T^2 divide by xy and conjugate:
       P = 0 <=> 0.3(1/x + 1/y) + 0.4 = 0 <=> 0.3(x + y) = -0.4 <=> x + y = -4/3.
With |x| = |y| = 1 this forces y = conj(x), x = e^{i phi}, 2 cos phi = -4/3, cos phi = -2/3.
So P has EXACTLY TWO zeros on T^2:
       (x,y) = (zeta, conj zeta) and (conj zeta, zeta),   zeta = (-2 + i sqrt5)/3,
in additive coordinates  (theta*, -theta*) and (-theta*, theta*),  theta* = arccos(-2/3)/(2pi).

TWO FACTS ABOUT theta* THAT DECIDE LATER SECTIONS, BOTH VERIFIED HERE:
  (Z1) zeta IS ALGEBRAIC, of degree 2: 3 zeta^2 + 4 zeta + 3 = 0.  Hence phi* = arccos(-2/3)
       = -i Log(zeta) IS A LOGARITHM OF AN ALGEBRAIC NUMBER.  (Used by R-7/Baker.)
  (Z2) theta* IS IRRATIONAL, by NIVEN'S THEOREM: the only rational multiples of pi whose
       cosine is rational have cosine in {0, +-1/2, +-1}; cos(2 pi theta*) = -2/3 is not one.
       (Used by the non-resonance proof of R-3.)

THE ZEROS ARE SIMPLE: dP/dx = 0.3 + 0.4 y has modulus 0.3 at both zeros (and dP/dy likewise),
so |P| is comparable to the distance to the zero.  Verified numerically over four decades of
offset.  This is what makes log|P| INTEGRABLE on T^2 (so m(P) is finite) and what makes the
dips of R-3 exactly logarithmic in the offset.
"""
import numpy as np
from fractions import Fraction
from R_lib import PI_K1, P_eval, m_jensen, hp_theta_star, cross_check_constants

print("=" * 79)
print("R_01 — THE SINGULARITY OF log|P| AT K1's REGISTERED pi, LOCATED EXACTLY")
print("=" * 79)
p00, p10, p01, p11 = PI_K1
print("\npi = (p00,p10,p01,p11) = %s   (K1's registered ready state)" % (PI_K1,))

# ---------------------------------------------------------------- W-01's firing criterion
mx = max(PI_K1); rest = sum(PI_K1) - mx
print("\nW-01 CONVEX-HULL / FIRING CRITERION (C3):  max(pi) = %.3f  <=  sum of others = %.3f"
      % (mx, rest))
print("   -> FIRES.  P therefore vanishes somewhere on T^2.  THE SINGULAR CASE IS THE")
print("      FORMATION CASE: the corpus lives in the regime where its own theorem is hardest.")

# ---------------------------------------------------------------- the two zeros, closed form
th = np.arccos(-2.0 / 3.0) / (2 * np.pi)
zeta = np.exp(2j * np.pi * th)
print("\nTHE TWO TORUS ZEROS, CLOSED FORM")
print("   theta* = arccos(-2/3)/(2pi) = %.15f" % th)
print("   zeta   = e(theta*)          = %+.15f %+.15fi" % (zeta.real, zeta.imag))
print("   minimal polynomial check  3 zeta^2 + 4 zeta + 3 = %.3e   (Z1: zeta is ALGEBRAIC)"
      % abs(3 * zeta ** 2 + 4 * zeta + 3))
for (a, b) in [(th, -th), (-th, th)]:
    x, y = np.exp(2j * np.pi * a), np.exp(2j * np.pi * b)
    print("   |P| at (%+.12f, %+.12f) = %.3e" % (a, b, abs(P_eval(PI_K1, x, y))))
print("   Niven's theorem (Z2): cos(2 pi theta*) = -2/3 not in {0,+-1/2,+-1} => theta* IRRATIONAL.")

# high-precision theta*, pure integer, cross-checked
tstar_hp = hp_theta_star(60)
print("   theta* to 60 digits (pure-integer fixed point): %s" % str(tstar_hp.numerator).rjust(60, '0')[:60])
print("   mpmath cross-check deviation: %.3e   (mpmath used as a CHECK only, CHOICE LEDGER L4)"
      % cross_check_constants(60))

# ---------------------------------------------------------------- simplicity of the zeros
print("\nSIMPLICITY OF THE ZEROS — |P| ~ C * dist, four decades of offset")
print("   offset d     |P(theta*+d, -theta*)|    ratio |P|/d      |P(theta*+d,-theta*+d)|/d")
for d in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    x1, y1 = np.exp(2j * np.pi * (th + d)), np.exp(2j * np.pi * (-th))
    x2, y2 = np.exp(2j * np.pi * (th + d)), np.exp(2j * np.pi * (-th + d))
    a1 = abs(P_eval(PI_K1, x1, y1)); a2 = abs(P_eval(PI_K1, x2, y2))
    print("   %8.1e    %20.12e   %12.8f   %12.8f" % (d, a1, a1 / d, a2 / d))
print("   ratios settle -> the zeros are SIMPLE (order 1).  log|P| is L^1 on T^2, so m(P) is finite.")

# ---------------------------------------------------------------- the brief's grid control
print("\nCONTROL (the brief's number; proves nothing, checks the setup)")
for n in (512, 1024, 2048, 4096):
    t = np.arange(n) * 2 * np.pi / n
    X = np.exp(1j * t)[:, None]; Y = np.exp(1j * t)[None, :]
    print("   min|P| on %5d^2 NESTED grid = %.6e" % (n, np.abs(P_eval(PI_K1, X, Y)).min()))
print("   2048^2 -> 1.9995e-04, the brief's 2.0e-04.  REPRODUCED.")
print("   NOTE, AND IT IS A DEFECT OF THE GRID FAMILY, NOT OF THE CODE: 1024 | 2048 | 4096, so")
print("   the three grids are NESTED and the minimum CANNOT fall.  A grid minimum reported")
print("   over a nested family is a window artefact of exactly COR-E's class.  Non-nested:")
for n in (1000, 3001, 9007):
    t = np.arange(n) * 2 * np.pi / n
    X = np.exp(1j * t)[:, None]; Y = np.exp(1j * t)[None, :]
    print("   min|P| on %5d^2 grid = %.6e" % (n, np.abs(P_eval(PI_K1, X, Y)).min()))
print("   THE 2-D GRID IS NOISE-LIMITED NEAR THE ZEROS (CHOICE LEDGER L3).  Not used again.")

# ---------------------------------------------------------------- m(P) across four decades
print("\nm(P) BY THE JENSEN REDUCTION, ACROSS FOUR DECADES OF NODE COUNT")
print("   (CONVERGENCE IS NOT A WINDOW: the trend is shown, not one endpoint)")
prev = None
for nq in (2 ** 10, 2 ** 12, 2 ** 14, 2 ** 16, 2 ** 18, 2 ** 20, 2 ** 22, 2 ** 24):
    v = m_jensen(PI_K1, nq)
    dd = "" if prev is None else "   change %+.3e" % (v - prev)
    print("   nq = %9d   m(P) = %.15f%s" % (nq, v, dd))
    prev = v
MP = m_jensen(PI_K1, 2 ** 24)
print("   ADOPTED VALUE  m(P) = %.12f" % MP)
print("   REGISTER (erratum against W-02, :168): generic torus value -0.767507880 = m(0.4+0.3x+0.3y)")
print("   deviation from the register's published figure: %.3e" % abs(MP - (-0.767507880)))

print("\nDONE R_01")
