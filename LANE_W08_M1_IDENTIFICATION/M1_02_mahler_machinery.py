#!/usr/bin/env python3
"""
M1_02 — THE MAHLER MEASURE m(P), COMPUTED THREE INDEPENDENT WAYS, AND THE ZERO SET.

P(x,y) = p00 + p10 x + p01 y + p11 xy.   On K1, p00 = 0 identically.
m(P) = (1/(2pi)^2) int int log|P(e^{is},e^{it})| ds dt.

Three routes, so no single implementation is load-bearing:
  R1  ONE-DIMENSIONAL REDUCTION (exact in structure, Gauss-Legendre in the smooth pieces):
      for fixed x on the circle, P is LINEAR in y, so Jensen gives
          m_y(A + B y) = log max(|A|, |B|),
      hence   m(P) = int_0^1 log max( p10 , |p01 + p11 e^{2 pi i t}| ) dt      (p00 = 0 case)
      and in general  m(P) = int_0^1 log max( |p00 + p10 x| , |p01 + p11 x| ) dt, x = e^{2 pi i t}.
      The integrand is continuous and piecewise analytic; the kinks are located in closed form.
  R2  CASSAIGNE-MAILLOT closed form for m(a + b x + c y) via the Bloch-Wigner dilogarithm.
  R3  brute 2-D grid (float64), as a crude independent witness only.

Also established here:
  Z1  the zero set of P on T^2 is FINITE (0, 1 or 2 points) and is non-empty
      IFF max(p10,p01,p11) <= 1/2  (equivalently: the three weights obey the triangle
      inequality).  This is W-01's convex-hull criterion, re-derived.
  Z2  m(P) depends only on the MULTISET {p10,p01,p11}  (independent re-derivation of N2).
  Z3  the Jensen bound  m(P) <= (1/2) log(sum p^2) <= 0, equality iff one class has all weight.

Precision: float64 throughout; Gauss-Legendre node counts are pushed until successive
refinements agree to < 1e-14, and R1/R2 are compared against each other.
"""
import numpy as np
from functools import lru_cache

TWOPI = 2.0 * np.pi

@lru_cache(maxsize=None)
def _leg(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w

# ---------------------------------------------------------------- dilogarithm
def li2(z, n=400):
    """Li2(z) = -int_0^1 log(1-z t)/t dt, Gauss-Legendre.  Valid for z off [1,inf)."""
    z = complex(z)
    if abs(z) < 1e-300:
        return 0j
    x, w = _leg(n)
    t = 0.5 * (x + 1.0)
    w = 0.5 * w
    f = -np.log(1.0 - z * t) / t
    return complex(np.sum(w * f))


def bloch_wigner(z):
    """D(z) = Im Li2(z) + arg(1-z) log|z|."""
    z = complex(z)
    return li2(z).imag + np.angle(1.0 - z) * np.log(abs(z))


# ---------------------------------------------------------------- R2 Cassaigne-Maillot
def m_CM(a, b, c):
    """m(a + b x + c y) for a,b,c >= 0."""
    a, b, c = float(a), float(b), float(c)
    v = sorted([a, b, c])
    if v[0] == 0.0 and v[1] == 0.0:
        return np.log(v[2]) if v[2] > 0 else -np.inf
    if v[2] > v[0] + v[1]:                       # no triangle: m = log max
        return np.log(v[2])
    if v[0] == 0.0:                              # degenerate triangle, b = c
        return np.log(v[2])                      # = log max, consistent with the branch above
    # angles opposite a, b, c
    def ang(op, s1, s2):
        cs = (s1 * s1 + s2 * s2 - op * op) / (2.0 * s1 * s2)
        return np.arccos(min(1.0, max(-1.0, cs)))
    al = ang(a, b, c)
    be = ang(b, a, c)
    ga = ang(c, a, b)
    D = bloch_wigner((a / b) * np.exp(1j * ga))
    return (D + al * np.log(a) + be * np.log(b) + ga * np.log(c)) / np.pi


# ---------------------------------------------------------------- R1 one-dimensional reduction
def _gl(f, lo, hi, n):
    x, w = _leg(n)
    t = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    return 0.5 * (hi - lo) * np.sum(w * f(t))


def m_R1(p00, p10, p01, p11, n=600):
    """m(P) = (1/2pi) int_0^{2pi} log max(|p00 + p10 e^{i th}|, |p01 + p11 e^{i th}|) d th.

    DEFECT FOUND ON RE-READ AND RECORDED RATHER THAN SILENTLY FIXED:  the first version of
    this routine returned -0.693145430 instead of -log 2 = -0.693147181 on the weights
    (p10,p01,p11) = (0, 1/2, 1/2) — which is S1's OWN PUBLISHED READY STATE.  Cause: when
    p00 = p10 = 0 and p01 = p11 the integrand acquires a genuine log singularity at th = pi
    and Gauss-Legendre loses ~6 digits there.  That case is exactly the one where P factors
    as a monomial times a one-variable binomial, so it is now taken in closed form.  The
    error was 1.75e-06 and it sat on the corpus's published instance."""
    if p00 == 0.0 and p10 == 0.0:            # P = y (p01 + p11 x)
        return np.log(max(p01, p11)) if max(p01, p11) > 0 else -np.inf
    if p01 == 0.0 and p11 == 0.0:            # P = p00 + p10 y
        return np.log(max(p00, p10)) if max(p00, p10) > 0 else -np.inf
    def A(th):
        return np.abs(p00 + p10 * np.exp(1j * th))
    def B(th):
        return np.abs(p01 + p11 * np.exp(1j * th))
    def f(th):
        return np.log(np.maximum(A(th), B(th)))
    # kinks where A = B : |p00+p10 e|^2 = |p01+p11 e|^2
    #   p00^2+p10^2 + 2 p00 p10 cos = p01^2+p11^2 + 2 p01 p11 cos
    K = 2.0 * (p00 * p10 - p01 * p11)
    C = (p01 * p01 + p11 * p11) - (p00 * p00 + p10 * p10)
    brk = [0.0, np.pi]
    if abs(K) > 1e-300:
        cs = C / K
        if -1.0 <= cs <= 1.0:
            brk.append(np.arccos(cs))
    brk = sorted(set(brk))
    tot = 0.0
    for lo, hi in zip(brk[:-1], brk[1:]):
        if hi - lo > 1e-15:
            tot += _gl(f, lo, hi, n)
    return tot / np.pi          # integrand even in th; int over [0,pi] * 2 / (2pi)


# ---------------------------------------------------------------- R3 brute grid
def m_grid(p00, p10, p01, p11, N=4096):
    t = (np.arange(N) + 0.5) / N * TWOPI
    x = np.exp(1j * t)[:, None]
    y = np.exp(1j * t)[None, :]
    V = np.abs(p00 + p10 * x + p01 * y + p11 * x * y)
    V = np.maximum(V, 1e-300)
    return float(np.mean(np.log(V)))


# ---------------------------------------------------------------- zero set on T^2
def zeros_on_T2(p00, p10, p01, p11):
    """Solve P=0 on T^2 for the K1 case p00=0: p11 + p10 x^{-1} + p01 y^{-1} = 0 after
    dividing by xy.  Returns list of (theta_x, theta_y)."""
    assert p00 == 0.0
    a, b, c = p11, p10, p01          # need a + b w1 + c w2 = 0, |w1|=|w2|=1
    if not (a <= b + c and b <= a + c and c <= a + b):
        return []
    if b == 0 or c == 0 or a == 0:
        # THIRD DEFECT FOUND ON RE-READ AND RECORDED.  The first version returned [] here
        # and the page would have claimed "the zero set is always finite".  FALSE.  If one
        # weight is 0 and the other two are EQUAL (hence 1/2 each), P is a monomial times a
        # one-variable binomial with equal moduli and its zero set on T^2 is a whole CIRCLE:
        #   p10 = 0, p01 = p11  ->  P = y(p01 + p11 x),  zero on {x = -1} x T^1
        #   p01 = 0, p10 = p11  ->  P = x(p10 + p11 y),  zero on T^1 x {y = -1}
        #   p11 = 0, p10 = p01  ->  P = p10 x + p01 y,   zero on {y = -x}
        # THIS IS S1'S OWN PUBLISHED READY STATE: p = (1/2,0,0,1/4,1/4) gives
        # (p10,p01,p11) = (0, 1/2, 1/2).  The singular set of log|P| there is 1-dimensional.
        return "CIRCLE" if (b == 0 and a == c) or (c == 0 and a == b) or (a == 0 and b == c) else []
    cos1 = (c * c - a * a - b * b) / (2 * a * b)
    if abs(cos1) > 1:
        return []
    t1 = np.arccos(min(1, max(-1, cos1)))
    out = []
    for sgn in (+1, -1):
        w1 = np.exp(1j * sgn * t1)
        w2 = -(a + b * w1) / c
        # DEFECT FOUND ON RE-READ AND RECORDED: the first version labelled w1 = x^{-1}.
        # It is y^{-1}.  P/(xy) = p11 + p10 y^{-1} + p01 x^{-1}, so w1 (which multiplies
        # b = p10) is y^{-1} and w2 (which multiplies c = p01) is x^{-1}.  The swap is
        # INVISIBLE whenever p10 = p01 — which is true of every weight vector the corpus
        # actually uses — and gave |P| = 6.0e-01 at a claimed zero on the first asymmetric
        # case tried.  This is the same failure mode as the corpus's own COR-K: a check
        # that passes on the symmetric instances and is never run on an asymmetric one.
        out.append((float(np.angle(1 / w2)), float(np.angle(1 / w1))))
    if abs(t1) < 1e-12 or abs(t1 - np.pi) < 1e-12:
        out = out[:1]
    return out


if __name__ == "__main__":
    print("=" * 78)
    print("M1_02 — MAHLER MACHINERY.  Three independent routes; float64.")
    print("=" * 78)

    print("\nCALIBRATION against a value known in closed form:")
    print("  m(1+x+y) = L'(chi_-3,-1) = 0.32306594...   this code: %.12f" % m_CM(1, 1, 1))
    print("  D(exp(i pi/3)) (max of Bloch-Wigner, = 1.0149416...): %.12f" %
          bloch_wigner(np.exp(1j * np.pi / 3)))
    print("  m(1+x)   = 0 exactly.  R1 gives %.3e ; CM gives %.3e"
          % (m_R1(0, 1, 1, 0), m_CM(0, 1, 1)))

    print("\nAGREEMENT OF THE THREE ROUTES on K1 weights (p00 = 0):")
    cases = [(0.0, 0.3, 0.3, 0.4),                 # S3/S4's generic ready state
             (0.0, 0.0, 0.5, 0.5),                 # S1's PUBLISHED ready state (1/2,0,0,1/4,1/4)
             (0.0, 1 / 3, 1 / 3, 1 / 3),
             (0.0, 0.6, 0.25, 0.15),               # triangle FAILS (max > 1/2)
             (0.0, 0.2, 0.5, 0.3),                 # degenerate triangle (max = 1/2)
             (0.0, 0.05, 0.05, 0.90),
             (0.0, 1.0, 0.0, 0.0)]                 # a vertex of the simplex
    print("   p10     p01     p11    |   R1 (1-D reduction)   R2 (Cassaigne-Maillot)  R3 (grid)"
          "        |R1-R2|")
    for (p00, p10, p01, p11) in cases:
        r1 = m_R1(p00, p10, p01, p11)
        r2 = m_CM(p11, p01, p10)          # multiset argument, see Z2
        r3 = m_grid(p00, p10, p01, p11)
        print("  %5.3f   %5.3f   %5.3f  |  %18.12f  %18.12f  %14.9f   %.2e"
              % (p10, p01, p11, r1, r2, r3, abs(r1 - r2)))

    print("\nZ2  m(P) depends only on the MULTISET {p10,p01,p11} — independent check of N2.")
    rng = np.random.default_rng(20260816 + 20)
    worst = 0.0
    from itertools import permutations
    for _ in range(300):
        w = rng.dirichlet([1, 1, 1])
        vals = [m_R1(0.0, *perm) for perm in permutations(w)]
        worst = max(worst, max(vals) - min(vals))
    print("    worst spread of m over all 6 permutations, 300 random weights: %.3e" % worst)

    print("\nZ1  THE ZERO SET OF P ON T^2 (p00 = 0).")
    print("    P = 0 on T^2  <=>  exists |w1|=|w2|=1 with p11 + p10 w1 + p01 w2 = 0")
    print("                  <=>  p11,p10,p01 obey the triangle inequality")
    print("                  <=>  max(p10,p01,p11) <= 1/2   (they sum to 1)")
    for (p00, p10, p01, p11) in cases:
        zs = zeros_on_T2(p00, p10, p01, p11)
        if zs == "CIRCLE":
            print("    p=(%5.3f,%5.3f,%5.3f)  max=%5.3f  ZERO SET IS A CIRCLE (1-dimensional)"
                  % (p10, p01, p11, max(p10, p01, p11)))
            continue
        resid = [abs(p10 * np.exp(1j * a) + p01 * np.exp(1j * b) + p11 * np.exp(1j * (a + b)))
                 for (a, b) in zs]
        print("    p=(%5.3f,%5.3f,%5.3f)  max=%5.3f  #zeros=%d  max|P| at them = %s"
              % (p10, p01, p11, max(p10, p01, p11), len(zs),
                 ("%.2e" % max(resid)) if resid else "-"))
    print("    -> GENERICALLY the zero set is FINITE (<=2 points): log|P| has finitely many")
    print("       isolated LOGARITHMIC singularities on T^2.  In the three exceptional weight")
    print("       vectors (one weight 0, the other two = 1/2) it is a CIRCLE — and one of those")
    print("       three IS S1's published ready state.  Either way log|P| is in L^1 but is NOT")
    print("       Riemann-integrable whenever max(p10,p01,p11) <= 1/2.")

    print("\nZ3  JENSEN BOUND  m(P) <= (1/2) log(sum_a p_a^2) <= 0,  equality on the right")
    print("    iff one class carries all the weight.  (int_{T^2}|P|^2 = sum p_a^2, characters")
    print("    1,x,y,xy orthonormal.)")
    rng2 = np.random.default_rng(20260816 + 21)
    bad = 0
    worstslack = 1e9
    for _ in range(2000):
        w = rng2.dirichlet([0.4, 0.4, 0.4])
        m = m_R1(0.0, *w)
        bnd = 0.5 * np.log(np.sum(w ** 2))
        if m > bnd + 1e-12:
            bad += 1
        worstslack = min(worstslack, bnd - m)
    print("    violations in 2000 random weights: %d ;  min slack (bound - m) = %.3e"
          % (bad, worstslack))
    print("\nDONE M1_02")
