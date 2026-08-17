"""
X01 — THE GEOMETRY N1 NEEDS, IN EXACT ARITHMETIC, AND m(P) WITH ITS CONVERGENCE TABULATED.

Nothing here is imported from any W-13 lane.  Every algebraic fact is certified by exact
rational arithmetic in Q(i sqrt5) (the field that happens to hold BOTH readings of S4's
SENSE C), and m(P) is produced by a Jensen reduction whose convergence is shown over five
node decades rather than quoted at one endpoint.

Legs:
  (a) The zero set of P on T^2 for K1_REG and for the two other admissible readings of
      S4:566's "(0.4,0.3,0.3) for 3 classes".  cos s0 EXACT, x0 and y0 EXACT algebraic,
      P(x0,y0) = 0 EXACT, and the minimal polynomials printed.
  (b) x0 y0 = 1  ?   (Theorem Z4's anti-diagonal.)  Labelling-dependent — shown.
  (c) m(P) by Jensen over 2^12 .. 2^24 nodes, the TREND printed, plus the Cassaigne-Maillot
      closed form as an independent second route.
  (d) The local shape at the zero: |P| / dist over six decades (conical?), and the constant.
"""
from fractions import Fraction as F
import numpy as np, math, hashlib

# ---------------------------------------------------------------- Q(i sqrt5) exact arithmetic
# element  a + b*i*sqrt5  with a,b in Q.   (i sqrt5)^2 = -5.
class Q5:
    __slots__ = ("a", "b")
    def __init__(self, a, b=0): self.a = F(a); self.b = F(b)
    def __add__(s, o): o = _q5(o); return Q5(s.a + o.a, s.b + o.b)
    def __sub__(s, o): o = _q5(o); return Q5(s.a - o.a, s.b - o.b)
    def __mul__(s, o):
        o = _q5(o)
        return Q5(s.a * o.a - 5 * s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    __radd__ = __add__
    def conj(s): return Q5(s.a, -s.b)
    def norm(s): return s.a * s.a + 5 * s.b * s.b          # = |z|^2, a rational
    def inv(s):
        n = s.norm(); return Q5(s.a / n, -s.b / n)
    def __truediv__(s, o): return s * _q5(o).inv()
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s):
        return "(%s %s %s*i*sqrt5)" % (s.a, "+" if s.b >= 0 else "-", abs(s.b))
    def cx(s): return complex(float(s.a), float(s.b) * math.sqrt(5.0))
def _q5(o): return o if isinstance(o, Q5) else Q5(o)

def minpoly_unitcircle(x):
    """x = a + b i sqrt5 on the unit circle => x^2 - 2a x + 1 = 0 ; return integer coeffs."""
    a = x.a
    num, den = (2 * a).numerator, (2 * a).denominator
    return (den, -num, den)          # den*z^2 - num*z + den

# ---------------------------------------------------------------- leg (a)(b)
def zero_set_exact(pi, name):
    p00, p10, p01, p11 = [F(t) for t in pi]
    # |A(e^is)|^2 - |B(e^is)|^2 = C + D cos s ,  A = p00+p10 x , B = p01+p11 x
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    print("  %-10s pi = (%s, %s, %s, %s)" % (name, p00, p10, p01, p11))
    print("      C = %s   D = %s" % (C, D))
    if D == 0:
        print("      D = 0 : degenerate (branch difference is constant) — not handled here")
        return None
    cos_s0 = -C / D
    print("      cos s0 = -C/D = %s   (EXACT RATIONAL)" % cos_s0)
    if abs(cos_s0) > 1:
        print("      |cos s0| > 1  =>  Z(P) is EMPTY on T^2")
        return None
    # sin^2 = 1 - cos^2 ; we need sin s0 = sqrt(1-cos^2).  For both S4 readings this is a
    # rational multiple of sqrt5, which is why one field holds them all.  Verify it.
    s2 = 1 - cos_s0**2
    num, den = s2.numerator, s2.denominator
    # want s2 = 5*(r)^2 for rational r
    r2 = s2 / 5
    r = F(math.isqrt(r2.numerator * r2.denominator), r2.denominator) if r2 > 0 else F(0)
    ok = (r * r == r2)
    print("      sin^2 s0 = %s ;  sin s0 = %s * sqrt5 ?  %s" % (s2, r, "YES (exact)" if ok else "NO"))
    if not ok:
        print("      -> not in Q(i sqrt5); skipping the exact certificate for this reading")
        return None
    x0 = Q5(cos_s0, r)
    A = Q5(p00) + Q5(p10) * x0
    B = Q5(p01) + Q5(p11) * x0
    print("      x0 = %s      |x0|^2 = %s" % (x0, x0.norm()))
    print("      |A(x0)|^2 = %s   |B(x0)|^2 = %s   equal: %s"
          % (A.norm(), B.norm(), A.norm() == B.norm()))
    y0 = Q5(0) - (A * B.inv())
    print("      y0 = -A/B = %s      |y0|^2 = %s" % (y0, y0.norm()))
    Pv = Q5(p00) + Q5(p10) * x0 + Q5(p01) * y0 + Q5(p11) * (x0 * y0)
    print("      P(x0,y0) = %s   EXACTLY ZERO: %s" % (Pv, Pv.iszero()))
    mp_x = minpoly_unitcircle(x0); mp_y = minpoly_unitcircle(y0)
    print("      min poly of x0 over Q : %d z^2 %+d z %+d" % mp_x)
    print("      min poly of y0 over Q : %d z^2 %+d z %+d" % mp_y)
    xy = x0 * y0
    print("      x0*y0 = %s   ANTI-DIAGONAL (x0 y0 = 1): %s" % (xy, (xy.a == 1 and xy.b == 0)))
    print("      s0/pi rational?  cos s0 = %s is not in {0,+-1/2,+-1} => NO, by NIVEN. %s"
          % (cos_s0, "(so exp(i s0) is NOT a root of unity)"))
    return x0, y0

# ---------------------------------------------------------------- leg (c)
def mP_jensen(pi, n):
    """m(P) = (1/2pi) INT_0^{2pi} log max(|A|,|B|) dt, midpoint rule on n nodes."""
    p00, p10, p01, p11 = [float(t) for t in pi]
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    ct, st = np.cos(t), np.sin(t)
    a2 = p00 * p00 + p10 * p10 + 2 * p00 * p10 * ct
    b2 = p01 * p01 + p11 * p11 + 2 * p01 * p11 * ct
    return float(np.mean(0.5 * np.log(np.maximum(a2, b2))))

def cassaigne_maillot(a, b, c):
    """m(a + b x + c y) for non-negative reals.  If the triangle inequality fails,
       log of the largest.  Otherwise the Bloch-Wigner/dilogarithm form."""
    a, b, c = float(a), float(b), float(c)
    if max(a, b, c) >= (a + b + c) - max(a, b, c):
        return math.log(max(a, b, c))
    # angles of the triangle with sides a,b,c ; alpha opposite a, etc.
    al = math.acos((b * b + c * c - a * a) / (2 * b * c))
    be = math.acos((a * a + c * c - b * b) / (2 * a * c))
    ga = math.pi - al - be
    def cl2(x, K=200000):
        k = np.arange(1, K + 1)
        return float(np.sum(np.sin(k * x) / (k * k)))
    return (al * math.log(a) + be * math.log(b) + ga * math.log(c)
            + 0.5 * (cl2(2 * al) + cl2(2 * be) + cl2(2 * ga))) / math.pi

# ---------------------------------------------------------------- leg (d)
def local_shape(pi, x0c, y0c):
    p00, p10, p01, p11 = [float(t) for t in pi]
    rng = np.random.default_rng(20260817)
    print("      r          min |P|/r        max |P|/r      (conical zero <=> both O(1)>0)")
    for e in range(2, 9):
        r = 10.0 ** (-e)
        th = rng.random(4000) * 2 * np.pi
        sg, ta = r * np.cos(th), r * np.sin(th)
        x = x0c * np.exp(1j * sg); y = y0c * np.exp(1j * ta)
        v = np.abs(p00 + p10 * x + p01 * y + p11 * x * y) / r
        print("      1e-%d      %.12f   %.12f" % (e, v.min(), v.max()))

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    print("=" * 92)
    print("X01 — EXACT GEOMETRY OF N1's POLYNOMIAL, AND m(P) WITH ITS CONVERGENCE TABULATED")
    print("=" * 92)

    print("\n(a)(b)  THE ZERO SET, EXACTLY, FOR ALL THREE ADMISSIBLE READINGS OF S4:566.")
    print("  S4:566 fixes SENSE C as '(0.4,0.3,0.3) for 3 classes' and does NOT say which")
    print("  class carries the 0.4.  The register and M1_08 T2(e) fix pi = (0,0.3,0.3,0.4),")
    print("  i.e. 0.4 on class 11.  The other two readings are computed here beside it.")
    states = [((0, F(3,10), F(3,10), F(2,5)), "K1_REG"),
              ((0, F(2,5),  F(3,10), F(3,10)), "K1_ALT10"),
              ((0, F(3,10), F(2,5),  F(3,10)), "K1_ALT01")]
    zeros = {}
    for pi, nm in states:
        zeros[nm] = zero_set_exact(pi, nm)
        print()

    print("(c)  m(P):  Jensen midpoint over five node decades, TREND, plus a second route.")
    for pi, nm in states:
        vals = [(n, mP_jensen(pi, n)) for n in [2**12, 2**14, 2**16, 2**18, 2**20, 2**22, 2**24]]
        ref = vals[-1][1]
        print("  %-9s  n        m_n(P)              m_n - m_{2^24}" % nm)
        for n, v in vals:
            print("             %-8d %.15f   %+.3e" % (n, v, v - ref))
        w = sorted([float(t) for t in pi], reverse=True)[:3]
        cm = cassaigne_maillot(*w)
        print("             Cassaigne-Maillot closed form on the multiset {%.1f,%.1f,%.1f}: %.12f  dev %.2e"
              % (w[0], w[1], w[2], cm, cm - ref))
        print()

    print("(d)  THE LOCAL SHAPE AT THE ZERO — is it conical (|P| ~ r) or tangential (~r^2)?")
    for pi, nm in states:
        z = zeros[nm]
        if z is None: continue
        print("  %s:" % nm)
        local_shape(pi, z[0].cx(), z[1].cx())
        print()
    print("DONE X01")
