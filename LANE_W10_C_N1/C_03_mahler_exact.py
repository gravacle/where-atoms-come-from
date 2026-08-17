#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 3.  N1's MAHLER LIMIT ON A GENUINELY FOUR-TERM POLYNOMIAL.

S4:599 records for B0b: "genuinely 4-term and does not factor. QUADRATURE ONLY,
-0.810930216216".  W-03 corrected this to "log(4/9) exactly" without giving the mechanism.
This script (a) supplies the mechanism, (b) reproduces the number independently by the Jensen
reduction with convergence PRINTED across n, (c) says how many digits are defensible, and
(d) settles whether either four-class carrier's P has zeros on T^2 — the hypothesis question.

Precision: float64 for the convergence tables; mpmath at mp.dps = 50 for the digit defence;
exact rational arithmetic for the factorisation and branch-ratio tests.

NO 2D GRID IS USED FOR ANY REPORTED VALUE.  A 2D grid is run in section 5 only to exhibit its
failure mode where it fails, and its agreement where it does not.
"""
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import sys

mp.mp.dps = 50

# ---------------------------------------------------------------------------- the reduction
def jensen_integrand(t, p00, p10, p01, p11):
    """log max(|p00 + p10 e^{it}|, |p01 + p11 e^{it}|)  —  Jensen in y."""
    c = np.cos(t)
    A2 = p00 * p00 + p10 * p10 + 2 * p00 * p10 * c
    B2 = p01 * p01 + p11 * p11 + 2 * p01 * p11 * c
    return 0.5 * np.log(np.maximum(A2, B2))

def m_jensen(p, n, offset=0.5):
    """Trapezoid on a periodic integrand, nodes OFFSET by `offset` steps.
    offset = 0.5 (midpoint) is the convention of record here: it avoids landing exactly on
    t = 0 or t = pi, where a degenerate P can vanish and return -inf from one node."""
    p00, p10, p01, p11 = p
    t = 2 * np.pi * (np.arange(n) + offset) / n
    return float(np.mean(jensen_integrand(t, p00, p10, p01, p11)))

def m_jensen_mp(p, n, offset=mp.mpf(1) / 2):
    p00, p10, p01, p11 = [mp.mpf(x) for x in p]
    s = mp.mpf(0)
    for k in range(n):
        t = 2 * mp.pi * (k + offset) / n
        c = mp.cos(t)
        A2 = p00 ** 2 + p10 ** 2 + 2 * p00 * p10 * c
        B2 = p01 ** 2 + p11 ** 2 + 2 * p01 * p11 * c
        s += mp.log(mp.sqrt(mp.mpf(max(A2, B2))))
    return s / n

# ---------------------------------------------------------------------------- exact structure
def factors_exact(p):
    """P = p00 + p10 x + p01 y + p11 xy factors as (a+bx)(c+dy) iff p00*p11 == p10*p01."""
    p00, p10, p01, p11 = [Fr(x) for x in p]
    return p00 * p11 == p10 * p01

def branch_ratio_constant_exact(p):
    """|p00+p10 x| / |p01+p11 x| is constant on |x|=1  iff
       (p00^2+p10^2) p01 p11 == (p01^2+p11^2) p00 p10.
       Two solution branches: PROPORTIONAL (=> P factors) or REVERSED (=> P does not)."""
    p00, p10, p01, p11 = [Fr(x) for x in p]
    lhs = (p00 * p00 + p10 * p10) * p01 * p11
    rhs = (p01 * p01 + p11 * p11) * p00 * p10
    if lhs != rhs:
        return None
    if p00 * p11 == p10 * p01:
        return "PROPORTIONAL (P factors)"
    return "REVERSED  (p00,p10) prop (p11,p01)  — P does NOT factor, m(P) still closed-form"

def zeros_on_T2_exact(p):
    """P has a zero on T^2 iff the two Jensen branch moduli cross.
       |p00+p10 x| runs over [|p00-p10|, p00+p10] and |p01+p11 x| over [|p01-p11|, p01+p11],
       BOTH monotone in cos t, so they cross iff
            ((p00+p10) - (p01+p11)) * (|p00-p10| - |p01-p11|)  <=  0 .
       Exact in Fraction arithmetic."""
    p00, p10, p01, p11 = [Fr(x) for x in p]
    hi = (p00 + p10) - (p01 + p11)
    lo = abs(p00 - p10) - abs(p01 - p11)
    return hi * lo <= 0, hi, lo

def min_absP_numeric(p, n=4096):
    p00, p10, p01, p11 = p
    t = 2 * np.pi * np.arange(n) / n
    A = np.abs(p00 + p10 * np.exp(1j * t))
    B = np.abs(p01 + p11 * np.exp(1j * t))
    # min over y of |A + B y| for |y|=1 is | |A| - |B| |
    return float(np.min(np.abs(A - B)))

# ---------------------------------------------------------------------------- carriers' weights
B0b_U = (Fr(4, 9), Fr(2, 9), Fr(1, 9), Fr(2, 9))      # (p00, p10, p01, p11)
B4_U  = (Fr(1, 6), Fr(1, 6), Fr(1, 6), Fr(3, 6))
SENSE_C4 = (Fr(1, 4), Fr(1, 4), Fr(1, 4), Fr(1, 4))
K1_U  = (Fr(0), Fr(2, 5), Fr(2, 5), Fr(1, 5))

CASES = [("B0b SENSE U  {00:4,10:2,01:1,11:2}/9", B0b_U, "-0.810930216216", mp.log(mp.mpf(4)/9)),
         ("B4  SENSE U  {00:1,10:1,01:1,11:3}/6", B4_U,  "-0.693147181",    -mp.log(2)),
         ("SENSE C 4-class (1,1,1,1)/4",          SENSE_C4, "-1.386294361120", -mp.log(4)),
         ("K1  SENSE U  (0,2,2,1)/5  [3-class]",  K1_U,  "-0.756573586",    None)]

if __name__ == "__main__":
    print("=" * 100)
    print("SECTION 1 — THE MECHANISM S4 MISSED.  Why a NON-FACTORING four-term P can still be exact.")
    print("=" * 100)
    print("""
  RECIPROCITY.  For |x| = 1 and real a,b:   |a + b x| = |x| * |a xbar + b| = |b + a x|.
  So the Jensen branch moduli |p00 + p10 x| and |p01 + p11 x| have CONSTANT ratio iff
      (p00,p10) prop (p01,p11)   [P factors]        or    (p00,p10) prop (p11,p01)   [REVERSAL].
  On the REVERSAL locus (p00,p10) = kappa (p11,p01), P does NOT factor, yet
      m(P) = log max(kappa, 1) + log max(p01, p11)          <-- CLOSED FORM, two terms.
  This is the same reciprocity that supplies the transpositions (00 10) and (01 11) in the
  multiset theorem.  B0b sits ON the reversal locus.  "Does not factor" never implied
  "quadrature only", and S4's inference at :599 is the error, not its arithmetic.
""")
    for name, p, s4, closed in CASES:
        fac = factors_exact(p)
        br = branch_ratio_constant_exact(p)
        z, hi, lo = zeros_on_T2_exact(p)
        print(f"  {name}")
        print(f"     factors over R ?  {fac}      branch-ratio constant ?  {br}")
        print(f"     zeros on T^2 ?    {z}   (hi={hi}, lo={lo}, product={hi*lo})"
              f"   min|P| over T^2 (numeric) = {min_absP_numeric([float(q) for q in p]):.6f}")
    print("""
  B0b: (p00,p10) = (4,2)/9 = 2 * (2,1)/9 = 2 * (p11,p01).  kappa = 2.
       m(P) = log max(2,1) + log max(1/9, 2/9) = log 2 + log(2/9) = log(4/9).   EXACT.
  B4 : not on the reversal locus, but the second branch DOMINATES everywhere
       (|1+3x|^2 - |1+x|^2 = 8 + 4 cos t >= 4 > 0), so m(P) = log 3 - log 6 = -log 2.  EXACT.""")

    print("\n" + "=" * 100)
    print("SECTION 2 — INDEPENDENT REPRODUCTION BY THE JENSEN REDUCTION, WITH CONVERGENCE IN n.")
    print("            (float64 trapezoid; the closed form is the reference, not the target.)")
    print("=" * 100)
    for name, p, s4, closed in CASES:
        pf = [float(q) for q in p]
        print(f"\n  {name}     S4 publishes {s4}")
        ref = float(closed) if closed is not None else None
        prev = None
        for e in range(8, 25):
            n = 2 ** e
            val = m_jensen(pf, n)
            d_ref = f"{abs(val-ref):.3e}" if ref is not None else "   --    "
            d_prev = f"{abs(val-prev):.3e}" if prev is not None else "   --    "
            print(f"     n = 2^{e:<2d} = {n:>9d}   m = {val:.15f}   |m - closed| = {d_ref}"
                  f"   |m - m(n/2)| = {d_prev}")
            prev = val

    print("\n" + "=" * 100)
    print("SECTION 3 — THE DIGIT DEFENCE.  mpmath, mp.dps = %d." % mp.mp.dps)
    print("=" * 100)
    l49 = mp.log(mp.mpf(4) / mp.mpf(9))
    print(f"  log(4/9)  = {mp.nstr(l49, 40)}")
    print(f"  -log 2    = {mp.nstr(-mp.log(2), 40)}")
    print(f"  -log 4    = {mp.nstr(-mp.log(4), 40)}")
    for name, p, s4, closed in CASES[:3]:
        pmp = [mp.mpf(q.numerator) / q.denominator for q in p]
        for n in (256, 1024, 4096):
            v = m_jensen_mp(pmp, n)
            print(f"  {name[:28]:28s} n={n:5d}  mp Jensen = {mp.nstr(v, 30)}"
                  f"   |v - closed| = {mp.nstr(abs(v - closed), 6)}")
    print("""
  S4's published B0b figure is -0.810930216216 (12 decimal places).
  log(4/9) = -0.8109302162163287639560... so ALL TWELVE OF S4's PLACES ARE CORRECT
  (the 12th place rounds 2163|287 -> 216).  Digits defensible: as many as asked for — the
  value is a closed form in elementary constants, not a quadrature result.  The Jensen
  trapezoid reproduces it to float64 epsilon by n = 2^12 and to 25+ places in mpmath,
  because the integrand is real-analytic and periodic (no zero of P on T^2: section 4).""")

    print("\n" + "=" * 100)
    print("SECTION 4 — THE ZERO SET.  THE LOAD-BEARING FACT ABOUT BOTH FOUR-CLASS CARRIERS.")
    print("=" * 100)
    print("""  CRITERION (exact, proved on this page).  For non-negative weights,
      P has a zero on T^2   <=>   ((p00+p10) - (p01+p11)) * (|p00-p10| - |p01-p11|) <= 0.
  PROOF.  P = A(x) + B(x) y with A = p00+p10 x, B = p01+p11 x.  For fixed x on T, min over
  |y|=1 of |A + B y| is ||A| - |B||, so P has a torus zero iff |A(x)| = |B(x)| for some x.
  |A|^2 = p00^2+p10^2+2 p00 p10 cos t and |B|^2 likewise, both AFFINE and INCREASING in cos t
  (weights non-negative), so |A| sweeps [|p00-p10|, p00+p10] and |B| sweeps [|p01-p11|, p01+p11]
  monotonically together; they meet iff the two intervals' endpoints interleave.  []
  It reduces to M1's Z1 on K1: with p00 = 0 it becomes (2 p10 - 1)(1 - 2 p01)(1 - 2 p11) <= 0,
  i.e. max(p10,p01,p11) <= 1/2.  CHECKED below against a numeric minimiser.""")
    rng = np.random.default_rng(20260816 + 3)
    ntest = 20000
    w = rng.dirichlet(np.ones(4), size=ntest)
    # ANALYTIC min over T^2 of |P|, from the criterion's own proof:
    #   g(c) = |A|^2-|B|^2 is AFFINE in c = cos t, so it has at most one root in [-1,1].
    #   crossing  -> min |P| = 0 ; no crossing -> min is at an endpoint c = +-1.
    hi = (w[:, 0] + w[:, 1]) - (w[:, 2] + w[:, 3])          # (|A|-|B|) at t = 0
    lo = np.abs(w[:, 0] - w[:, 1]) - np.abs(w[:, 2] - w[:, 3])  # (|A|-|B|) at t = pi
    pred_zero = (hi * lo <= 0)
    # correct analytic min of ||A|-|B|| over c = cos t in [-1,1] when there is no crossing:
    # h(c) = |sqrt(a+b c) - sqrt(a'+b' c)| has at most ONE interior critical point (the
    # stationarity condition b^2(a'+b'c) = b'^2(a+bc) is linear in c).  Candidates: c = +-1
    # and that point when it lies in range.
    a  = w[:, 0] ** 2 + w[:, 1] ** 2; b  = 2 * w[:, 0] * w[:, 1]
    a2 = w[:, 2] ** 2 + w[:, 3] ** 2; b2 = 2 * w[:, 2] * w[:, 3]
    den = b * b * b2 - b2 * b2 * b
    with np.errstate(divide="ignore", invalid="ignore"):
        ccrit = np.where(np.abs(den) > 1e-14, (b2 * b2 * a - b * b * a2) / den, np.nan)
    def h(c):
        return np.abs(np.sqrt(np.maximum(a + b * c, 0)) - np.sqrt(np.maximum(a2 + b2 * c, 0)))
    cand = np.vstack([h(np.ones(ntest)), h(-np.ones(ntest)),
                      np.where(np.isfinite(ccrit) & (np.abs(ccrit) <= 1),
                               h(np.clip(np.nan_to_num(ccrit), -1, 1)), np.inf)])
    analytic = np.where(pred_zero, 0.0, cand.min(axis=0))
    numeric = np.array([min_absP_numeric(w[i], 8192) for i in range(ntest)])
    nozero = ~pred_zero
    err = np.abs(numeric - analytic)[nozero]
    print(f"  ANALYTIC min|P| over T^2 (from the criterion's own proof) vs NUMERIC min over an")
    print(f"  8192-node x-grid, on the {int(nozero.sum())} of {ntest} draws with NO crossing:")
    print(f"       max |analytic - numeric| = {err.max():.3e}   mean = {err.mean():.3e}")
    # binary test done the only way a grid can do it: does |A|^2-|B|^2 change sign on the grid?
    tt = 2 * np.pi * (np.arange(4096) + 0.5) / 4096
    cc = np.cos(tt)
    gg = (a[:, None] - a2[:, None]) + (b[:, None] - b2[:, None]) * cc[None, :]
    signchange = (gg.min(axis=1) * gg.max(axis=1)) <= 0
    disagree = int(np.sum(signchange != pred_zero))
    print(f"  BINARY test (does |A|^2-|B|^2 change sign on a 4096-node grid?) vs the criterion:"
          f"  {ntest - disagree} / {ntest} agree,  {disagree} disagree")
    # measure of the has-a-zero set
    big = 4_000_000
    w2 = rng.dirichlet(np.ones(4), size=big)
    g1 = (w2[:, 0] + w2[:, 1]) ** 2 - (w2[:, 2] + w2[:, 3]) ** 2
    gm = (w2[:, 0] - w2[:, 1]) ** 2 - (w2[:, 2] - w2[:, 3]) ** 2
    f4 = float(np.mean(g1 * gm <= 0))
    w3 = rng.dirichlet(np.ones(3), size=big)
    f3 = float(np.mean(w3.max(axis=1) <= 0.5))
    print(f"""
  HOW BIG IS THE ZERO-HAVING SET IN THE READY-STATE SIMPLEX?  Uniform (Dirichlet(1,1,1,1)):
     FOUR classes  P(P has a zero on T^2) = {f4:.6f}   over {big} draws
     THREE classes P(P has a zero on T^2) = {f3:.6f}   over {big} draws   [the medial triangle]
  BOTH ARE EXACTLY 1/4, and the four-class value is proved here:
     write s = p00+p10, d1 = p00-p10, d2 = p01-p11.  Uniform on the 3-simplex is uniform on
     {{0<=s<=1, |d1|<=s, |d2|<=1-s}}.  The condition is (2s-1)(|d1|-|d2|) <= 0, whose conditional
     probability is (1-s)/(2s) for s>1/2 and s/(2(1-s)) for s<1/2; integrating against the
     density 4s(1-s)/(2/3) gives (3/2)(1/12 + 1/12) = 1/4.  []
  SO: MOVING FROM THREE CLASSES TO FOUR DOES NOT CHANGE HOW OFTEN log|P| IS SINGULAR.
  It is 1/4 of ready states either way.  The brief's premise — "a four-class P HAS torus zeros
  (W-09: the firing region is exactly 1/2)" — CONFLATES TWO DIFFERENT SETS.  W-09's 1/2 is the
  measure, in the CONNECTION space (f,c), of the set where SOME non-negative weight vector
  annihilates Z_1.  The Lawton question is about a FIXED weight vector's zero set in T^2.  The
  first is a union over p of the second, and a union of measure-zero sets is not measure zero.""")
    print("""
  APPLIED:
     B0b SENSE U   hi = +3/9 > 0,  lo = +1/9 > 0   ->  NO ZERO ON T^2.  min|P| = 1/9 exactly
                   (|A| - |B| = 2|1+2x|/9 - |1+2x|/9 = |1+2x|/9 >= 1/9).
     B4  SENSE U   hi = -2/6 < 0,  lo = -2/6 < 0   ->  NO ZERO ON T^2.  min|P| = sqrt(6)/9
                   = 0.2721655269759087, attained at cos t = -2/3 (an INTERIOR critical point
                   of ||A|-|B||, not an endpoint: 3 sqrt(2+2c) = sqrt(10+6c) gives c = -2/3).
                   [DEFECT D-5, recorded: the first draft of this line asserted 1/3, the
                    endpoint value, without computing the interior critical point.]
     SENSE C 4-cls hi =  0                          ->  A WHOLE CIRCLE of zeros at x = -1.
  BOTH FOUR-CLASS CARRIERS THE CORPUS OWNS ARE ZERO-FREE ON T^2 UNDER SENSE U.
  Therefore log|P| is CONTINUOUS on T^2 for them, Weyl alone gives Birkhoff convergence, and
  the Lawton/Diophantine question the brief calls load-bearing DOES NOT BITE ON EITHER OF THEM.
  Under SENSE C — the corpus's OTHER published four-class column, lambda = -1.386294361120 —
  it bites maximally: the zero set is one-dimensional.""")

    print("\n" + "=" * 100)
    print("SECTION 5 — THE 2D GRID, RUN ONLY TO SHOW WHERE IT FAILS AND WHERE IT DOES NOT.")
    print("=" * 100)
    for name, p, s4, closed in CASES[:3]:
        pf = [float(q) for q in p]
        print(f"\n  {name}   closed form {float(closed):.12f}")
        for n in (256, 1024, 4096):
            t = 2 * np.pi * (np.arange(n) + 0.5) / n
            X, Y = np.meshgrid(np.exp(1j * t), np.exp(1j * t), indexing="ij")
            Pv = pf[0] + pf[1] * X + pf[2] * Y + pf[3] * X * Y
            g = float(np.mean(np.log(np.abs(Pv))))
            print(f"     2D grid {n}x{n}   m = {g:.12f}   |err| = {abs(g-float(closed)):.3e}")
    print("""
  The 2D grid is fine on B0b and B4 (analytic integrand) and is noise-limited only where P has
  torus zeros.  S4's "QUADRATURE ONLY, stable to 12 places" for B0b was therefore CORRECT
  ARITHMETIC ON A WRONG PREMISE: the premise was that no closed form exists.""")
    sys.exit(0)
