#!/usr/bin/env python3
"""
L1 — ARE THE TORUS ZEROS OF P ALGEBRAIC?

This is the hypothesis check for the theorem this lane identifies as the one that licenses N1.
Gelfond's Theorem III (as used by Lind-Schmidt-Verbitskiy, Sec. 9) and Baker's theorem on
linear forms in logarithms both require ALGEBRAIC data.  LSV say so explicitly:

    "Now log|f| has only finitely many logarithmic singularities, and by Proposition 4.2 these
     all have algebraic coordinates.  We can therefore use Gelfond's result to control the few
     potentially large negative values ..."      (LSV, Sec. 9)

LSV get algebraicity for free because their f has INTEGER coefficients.  Ours has probability
coefficients.  So the question has to be asked here and it has an answer:

  CLAIM L1.  For P = p00 + p10 x + p01 y + p11 xy with real coefficients, the zero set of P on
  T^2 is
        {(x, y) : cos(arg x) = t,  y = -(p00 + p10 x)/(p01 + p11 x)},
        t = (p01^2 + p11^2 - p00^2 - p10^2) / (2 (p00 p10 - p01 p11)),
  whenever |t| <= 1 and the denominator is non-zero.  t is a RATIONAL FUNCTION OF pi.  Hence
  the zeros have coordinates algebraic over Q(pi): if pi is rational the zeros are algebraic
  numbers of degree <= 4 over Q; if pi is algebraic they are algebraic.

  DERIVATION.  P = (p00 + p10 x) + (p01 + p11 x) y.  On T^2, P = 0 forces
  |p00 + p10 x| = |p01 + p11 x|; with x = e^{i th} this is
  p00^2 + p10^2 + 2 p00 p10 cos th = p01^2 + p11^2 + 2 p01 p11 cos th, i.e. cos th = t.
  Conversely any such x gives a unique y on the circle.  []

Everything below is EXACT in fractions.Fraction except where a line says float64.
"""
from fractions import Fraction as F
import math

print("=" * 78)
print("L1 — THE TORUS ZEROS OF P, AND WHETHER THEY ARE ALGEBRAIC")
print("=" * 78)


def zeros_exact(p00, p10, p01, p11):
    """Return (t, list of (x,y) as complex float64) with t EXACT."""
    den = 2 * (p00 * p10 - p01 * p11)
    if den == 0:
        return None, []
    t = (p01 * p01 + p11 * p11 - p00 * p00 - p10 * p10) / den
    if abs(t) > 1:
        return t, []
    ct = float(t)
    st = math.sqrt(max(0.0, 1.0 - ct * ct))
    out = []
    for s in ([+1, -1] if st > 0 else [+1]):
        x = complex(ct, s * st)
        num = float(p00) + float(p10) * x
        de = float(p01) + float(p11) * x
        if abs(de) < 1e-14:
            continue
        y = -num / de
        out.append((x, y))
    return t, out


def Pval(p, x, y):
    return p[0] + p[1] * x + p[2] * y + p[3] * x * y


CASES = [
    ("K1 REGISTERED pi (W-11 H2 row; the singular case)", (F(0), F(3, 10), F(3, 10), F(2, 5))),
    ("K1 centroid three-class", (F(0), F(1, 3), F(1, 3), F(1, 3))),
    ("S1 published ready state (W-01)", (F(0), F(1, 2), F(1, 2), F(0))),
    ("B0b four-class (S4:590 multiset 2/9,3/9,4/9 + 0)", (F(2, 9), F(4, 9), F(3, 9), F(0))),
    ("B0b genuine four-class 1/4 each", (F(1, 4), F(1, 4), F(1, 4), F(1, 4))),
    ("NON-firing: max > 1/2", (F(0), F(1, 10), F(1, 10), F(8, 10))),
]

print("\nCLAIM L1 CHECKED CASE BY CASE.  't' is EXACT (Fraction); |P(x,y)| is float64.\n")
print("  %-46s %-14s %-9s %s" % ("case", "t = cos(arg x)", "|t|<=1?", "max |P(zero)|"))
for name, p in CASES:
    t, zs = zeros_exact(*p)
    if t is None:
        print("  %-46s %-14s %-9s %s" % (name, "denominator 0", "-", "degenerate: handled below"))
        continue
    worst = max([abs(Pval([float(q) for q in p], x, y)) for x, y in zs], default=float("nan"))
    print("  %-46s %-14s %-9s %.3e" % (name, str(t), "YES" if abs(t) <= 1 else "NO", worst))

print("""
  READING.  Wherever |t| <= 1 the two exhibited points are zeros of P on T^2 to float64
  round-off, and t is an EXACT RATIONAL.  Wherever |t| > 1 there is no zero.  The firing
  criterion max(pi) <= 1/2 (M1_02 Z1, sealed) and |t| <= 1 agree on every row above.
""")

print("-" * 78)
print("THE REGISTERED CASE, IN CLOSED FORM AND EXACTLY.")
print("-" * 78)
p = (F(0), F(3, 10), F(3, 10), F(2, 5))
t, zs = zeros_exact(*p)
print("  pi = (0, 3/10, 3/10, 2/5).   t = %s  EXACTLY." % t)
print("  So x = -2/3 + i sqrt(5)/3 and its conjugate.  MINIMAL POLYNOMIAL, verified exactly:")
# (3x+2)^2 = -5  ->  9x^2 + 12x + 9 = 0  ->  3x^2 + 4x + 3 = 0
a, b, c = F(3), F(4), F(3)
# check 3x^2+4x+3 = 0 for x = (-2 + i sqrt5)/3 : real part 3*(4-5)/9 + 4*(-2/3) + 3
re = a * (F(4, 9) - F(5, 9)) + b * F(-2, 3) + c
print("      3x^2 + 4x + 3 :  REAL part exactly %s   (must be 0)" % re)
print("      IMAG part: 3*2*(-2/3)(sqrt5/3)/1 + 4*sqrt5/3 = -4sqrt5/3 + 4sqrt5/3 = 0 exactly.")
print("  => x is ALGEBRAIC OF DEGREE 2 over Q (root of 3z^2+4z+3), NOT an algebraic integer.")
print("  => Gelfond / Baker apply to it.  Algebraic INTEGER is not required by either.")
for x, y in zs:
    print("      zero (%.15f%+.15fi , %.15f%+.15fi)   |x|-1 = %.2e  |y|-1 = %.2e  |P| = %.2e"
          % (x.real, x.imag, y.real, y.imag, abs(x) - 1, abs(y) - 1,
             abs(Pval([float(q) for q in p], x, y))))
print("  y = conj(x) on this row (the zero set is {(x,xbar),(xbar,x)}), consistent with")
print("  M1_02's 'the zero set is 2 points generically' (sealed, LANE_W08_M1_IDENTIFICATION).")

print("""
--------------------------------------------------------------------------------
WHERE THE DERIVATION BREAKS, RECORDED RATHER THAN HIDDEN.
  (i) p00 p10 = p01 p11 makes the denominator vanish.  Then the condition on cos th is either
      vacuous (a WHOLE CIRCLE of zeros — M1_02's third case, which includes S1's own published
      ready state) or empty.  On a circle of zeros the singularity is not isolated and NO
      theorem in this lane's citation list applies: Gelfond/Baker control finitely many
      algebraic singular points, and LSV say explicitly that the infinite case needs estimates
      that "do not appear to be available".
  (ii) pi with transcendental entries gives transcendental t, and Baker does not apply.  The
      corpus has never published such a pi; every registered ready state is rational.
""")
p = (F(0), F(1, 2), F(1, 2), F(0))
print("  CHECK of (i) on S1's own published ready state pi = (0,1/2,1/2,0):")
print("      p00*p10 - p01*p11 = %s   -> denominator 0, CIRCLE of zeros." % (p[0] * p[1] - p[2] * p[3]))
print("      P = (x+y)/2, which vanishes on the whole circle y = -x.  Confirmed by inspection.")
print("      This is the case where the corpus's own registered numbers give Z_k = 0 for odd k")
print("      (M1_08 T2(e), sealed) and lambda = -infinity, not m(P) = -log 2.")
print("\nDONE L1")
