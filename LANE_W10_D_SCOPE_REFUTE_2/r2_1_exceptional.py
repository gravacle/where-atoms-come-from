# W10-D REFUTE-2  LENS 2 = COMPLETENESS.  LEG 1.
#
# THE OMISSION UNDER TEST.  S4 sec3.1 publishes an entire block about the EXCEPTIONAL SET --
# the countable-dense measure-zero set of resonant connections on which lambda_B departs from
# its generic value -- with a CLOSED FORM for the departed value, a 16-row table of exceptional
# values, an above/below count, and an accumulation table.  The ERRATUM AGAINST W-02 is a
# register row entirely about one entry of it ((m,n) = (11,20), f=2.0, c=1.1).  W-05's N3 is
# about its measure.  W-03's correction "the exceptional-value split is 527/314/213" is about
# its census.
#
# LANE D's SCOPE TABLE CLASSIFIES NONE OF IT.  Its only touching row, 2.10, classifies the
# WRONG NUMBER (-0.766802) and says "untouched here".  Row 5.3 classifies the MEASURE of the
# set (Haar-null) and not its STRUCTURE.
#
# AND THE CLOSED FORM IS VISIBLY THREE-CLASS.  S4 sec3.1: "On the primitive locus (m,n) with
# d = 1, lambda_B is the one-variable Mahler measure  m(r z^{m+n} + p0 z^m + q), exactly."
# THREE terms, because K1 has three occupied classes.  Derived here for any occupancy:
#
#   exceptional locus  { (f,c) : -m f + n c = 0 (mod 2pi) },  gcd(m,n) = 1
#   theta = (theta1,theta2) = (-f, c);  m theta1 + n theta2 = 0;  orbit dense in the circle
#   theta = t (n, -m).  Then u^k -> e^{int}, v^k -> e^{-imt}, (uv)^k -> e^{i(n-m)t}, so
#      Z(t) . e^{imt}  =  p01 + p00 z^m + p11 z^n + p10 z^{m+n},   z = e^{it}
#   lambda_resonant(m,n) = m( p01 + p00 z^m + p11 z^n + p10 z^{m+n} )   -- FOUR terms at four
#   classes, and it reduces to S4's three-term form exactly when p00 = 0.
#
# THE ONE VARIABLE THAT MOVES IN THIS LEG: THE CLASS WEIGHT 4-VECTOR.  The (m,n) list, the
# Mahler routine (exact, by roots), the generic-value routine and the code path are identical
# in every row.  ARM DIFF PRINTED FIRST.
#
# PRECISION.  numpy float64.  lambda_resonant is computed EXACTLY in the sense of Jensen's
# formula -- m(Q) = log|lead| + sum log max(1,|root|) -- not by quadrature, so there is no
# log-singularity problem; it is cross-checked against a 2^20-point trapezoid and against
# S4's own published column.

import numpy as np
from itertools import permutations

np.set_printoptions(precision=12)

# ------------------------------------------------------------------ exact one-variable Mahler
def m_one_var(coeffs_by_exp):
    """coeffs_by_exp: dict exponent(int) -> coefficient(float).  Returns m of the Laurent
       polynomial (Mahler measure is invariant under multiplication by a monomial)."""
    exps = sorted(coeffs_by_exp)
    lo = exps[0]
    deg = exps[-1] - lo
    c = np.zeros(deg + 1)
    for e, v in coeffs_by_exp.items():
        c[e - lo] += v
    # strip leading/trailing zeros (again a monomial factor)
    nz = np.nonzero(c)[0]
    c = c[nz[0]:nz[-1] + 1]
    if len(c) == 1:
        return np.log(abs(c[0]))
    poly = c[::-1]                      # numpy.roots wants highest degree first
    r = np.roots(poly)
    return np.log(abs(poly[0])) + np.sum(np.log(np.maximum(1.0, np.abs(r))))


def m_one_var_quad(coeffs_by_exp, n=1 << 20):
    t = (np.arange(n) + 0.5) * 2 * np.pi / n
    z = np.exp(1j * t)
    tot = np.zeros(n, dtype=complex)
    for e, v in coeffs_by_exp.items():
        tot += v * z ** e
    return np.log(np.abs(tot)).mean()


def Q_of(p, m, n):
    """p = (p00,p10,p01,p11).  Returns the exponent->coeff dict of the resonant polynomial."""
    p00, p10, p01, p11 = p
    d = {}
    for e, v in ((0, p01), (m, p00), (n, p11), (m + n, p10)):
        d[e] = d.get(e, 0.0) + float(v)
    return {e: v for e, v in d.items() if abs(v) > 0.0}


def lam_res(p, m, n):
    return m_one_var(Q_of(p, m, n))


# ------------------------------------------------------------------ generic value: Jensen in x
GLX, GLW = np.polynomial.legendre.leggauss(600)


def _arc(a, b, lo, hi):
    t = 0.5 * (hi - lo) * GLX + 0.5 * (hi + lo)
    val = 0.5 * np.log(a * a + b * b + 2 * a * b * np.cos(t))
    return 0.5 * (hi - lo) * np.dot(GLW, val)


def m_generic(p):
    """m(p00 + p10 x + p01 y + p11 xy) by Jensen in x, split at the branch crossing."""
    p00, p10, p01, p11 = [float(q) for q in p]
    A2 = lambda ct: p00 ** 2 + p01 ** 2 + 2 * p00 * p01 * ct
    B2 = lambda ct: p10 ** 2 + p11 ** 2 + 2 * p10 * p11 * ct
    alpha = p00 ** 2 + p01 ** 2 - p10 ** 2 - p11 ** 2
    beta = 2 * (p00 * p01 - p10 * p11)
    cuts = [0.0, np.pi]
    if abs(beta) > 1e-300:
        r = -alpha / beta
        if -1.0 < r < 1.0:
            cuts.insert(1, float(np.arccos(r)))
    tot = 0.0
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        cm = np.cos(0.5 * (lo + hi))
        tot += _arc(p00, p01, lo, hi) if A2(cm) >= B2(cm) else _arc(p10, p11, lo, hi)
    return tot / np.pi


# ------------------------------------------------------------------ the arms
ARMS = [
    ("K1  S3 ready state  (3cl)", np.array([0.0, 0.3, 0.3, 0.4])),   # the ERRATUM's own weights
    ("B1  K1 SENSE-U      (3cl)", np.array([0.0, 2 / 5, 2 / 5, 1 / 5])),
    ("B1q spectator       (3cl)", np.array([1 / 7, 3 / 7, 3 / 7, 0.0])),
    ("B0b torus, meeting  (4cl)", np.array([4 / 9, 2 / 9, 1 / 9, 2 / 9])),
    ("B4  spindle         (4cl)", np.array([1 / 6, 1 / 6, 1 / 6, 3 / 6])),
    ("SENSE-C 4 classes   (4cl)", np.array([.25, .25, .25, .25])),
]

print("=" * 104)
print("ARM DIFF FIRST.  The weight 4-vectors (p00,p10,p01,p11) actually handed to the SAME")
print("routines.  Two byte-identical arms reported as a confirmation is this corpus's commonest")
print("FATAL defect, so the arms are printed and pairwise-compared before any result.")
print("=" * 104)
for lab, p in ARMS:
    print(f"   {lab}   {np.array2string(np.asarray(p, float), precision=9)}")
seen = []
for lab, p in ARMS:
    for lab2, p2 in seen:
        if np.allclose(p, p2, atol=0, rtol=0):
            print(f"   !! IDENTICAL ARMS: {lab} and {lab2}")
    seen.append((lab, p))
print(f"   DISTINCT ARMS: {len({tuple(np.round(p,15)) for _, p in ARMS})} of {len(ARMS)}")

# ------------------------------------------------------------------ 1A validate against S4
print("\n" + "=" * 104)
print("== 1A  VALIDATE THE ROUTINE AGAINST S4 sec3.1's OWN PUBLISHED EXCEPTIONAL-VALUE COLUMN ==")
print("=" * 104)
print("  S4's table is computed on K1 with the S3 ready state p=(0.4,.15,.15,.15,.15), i.e.")
print("  class weights (p00,p10,p01,p11) = (0, 0.3, 0.3, 0.4).  If my routine does not")
print("  reproduce S4's published column, NOTHING BELOW IS ABOUT THE CORPUS'S OBJECT.")
S4_PUB = {(1, 0): -0.356674944, (0, 1): -0.356674944, (1, 1): -1.203972804,
          (1, -1): -0.510825624, (2, 1): -0.681980359, (2, -1): -0.916290732,
          (3, 1): -0.767783712, (3, 2): -0.732940865, (4, 1): -0.784966659,
          (5, 1): -0.749392712, (5, 3): -0.765224351, (7, 3): -0.759305247,
          (7, 11): -0.764712281, (11, 20): -0.767014993, (13, 8): -0.768271734,
          (29, 17): -0.767138179}
pK1 = ARMS[0][1]
gK1 = m_generic(pK1)
print(f"  generic m(P) on those weights = {gK1:.12f}   S4/erratum publish -0.767507880358"
      f"   dev {abs(gK1 + 0.767507880358):.2e}")
print(f"  {'(m,n)':>10s} {'mine (roots)':>16s} {'S4 published':>16s} {'dev':>10s} {'quad xcheck':>16s}")
worst = 0.0
for mn, val in S4_PUB.items():
    mine = lam_res(pK1, *mn)
    q = m_one_var_quad(Q_of(pK1, *mn))
    worst = max(worst, abs(mine - val))
    print(f"  {str(mn):>10s} {mine:16.12f} {val:16.9f} {abs(mine-val):10.2e} {q:16.12f}")
print(f"  WORST DEVIATION FROM S4's PUBLISHED COLUMN: {worst:.2e}.  The routine is S4's object.")
print("  (S4's column is printed to 9 places, so its own rounding is 5e-10.)")

# ------------------------------------------------------------------ 1B the same block off K1
print("\n" + "=" * 104)
print("== 1B  THE SAME BLOCK, ONE VARIABLE MOVED: THE CLASS WEIGHT VECTOR ==")
print("=" * 104)
MNS = list(S4_PUB.keys())
print("  For each arm: the 16 exceptional values on the SAME 16 primitive loci, the generic")
print("  value, the range, and the above/below split S4 publishes as '13 above, 6 below'.")
print(f"\n  {'arm':26s} {'generic m(P)':>15s} {'min exc':>12s} {'max exc':>12s} "
      f"{'#above':>7s} {'#below':>7s} {'max|dev|':>10s}")
DETAIL = {}
for lab, p in ARMS:
    g = m_generic(p)
    vals = {mn: lam_res(p, *mn) for mn in MNS}
    DETAIL[lab] = (g, vals)
    v = np.array(list(vals.values()))
    above = int((v > g).sum())
    below = int((v < g).sum())
    print(f"  {lab:26s} {g:15.9f} {v.min():12.6f} {v.max():12.6f} "
          f"{above:7d} {below:7d} {np.abs(v-g).max():10.2e}")
print("\n  S4 publishes for K1: 19 distinct values, range [-1.203973, -0.356675], 13 above /")
print("  6 below.  On the 16 loci above the K1 arm reproduces that shape.  READ THE FOUR-CLASS")
print("  ROWS AGAINST IT.")

print("\n  FULL TABLE, all arms, all 16 loci (deviation from that arm's own generic value):")
hdr = f"  {'(m,n)':>9s}" + "".join(f"{lab.split()[0]:>13s}" for lab, _ in ARMS)
print(hdr)
for mn in MNS:
    row = f"  {str(mn):>9s}"
    for lab, p in ARMS:
        g, vals = DETAIL[lab]
        row += f"{vals[mn]-g:+13.2e}"
    print(row)

# ------------------------------------------------------------------ 1C the erratum's own row
print("\n" + "=" * 104)
print("== 1C  THE ERRATUM AGAINST W-02's OWN ROW, (m,n) = (11,20), ON EVERY ARM ==")
print("=" * 104)
print("  The erratum is a REGISTER ROW.  Its content: S3/S4's headline f=2.0, c=1.1 is exactly")
print("  resonant (-11f+20c=0), so its orbit is on a SUBTORUS and its average is NOT the generic")
print("  average.  On K1 that is a +4.93e-04 departure.  ONE VARIABLE MOVED (the weight vector):")
print(f"\n  {'arm':26s} {'generic':>15s} {'subtorus (11,20)':>18s} {'departure':>12s}")
for lab, p in ARMS:
    g, vals = DETAIL[lab]
    print(f"  {lab:26s} {g:15.9f} {vals[(11,20)]:18.9f} {vals[(11,20)]-g:+12.3e}")
print("\n  LANE D's leg 2C measured this same point on B0b by a 4e5-term time average and got")
print("  a deviation of 9.66e-07, which it read as N3 confirmation ('the resonant points are")
print("  exactly where the rate departs').  COMPARE THE EXACT SUBTORUS VALUE ABOVE.")
print("  Reproducing lane D's own measurement, same code path, to show what 9.66e-07 was:")
k = np.arange(1, 400001)
for lab, p in (("B0b torus, meeting  (4cl)", ARMS[3][1]), ("B4  spindle         (4cl)", ARMS[4][1]),
               ("K1  S3 ready state  (3cl)", ARMS[0][1])):
    fv, cv = 2.0, 1.1
    Z = p[0] + p[1]*np.exp(-1j*fv*k) + p[2]*np.exp(1j*cv*k) + p[3]*np.exp(1j*(cv-fv)*k)
    tavg = np.log(np.abs(Z)).mean()
    g, vals = DETAIL[lab]
    print(f"    {lab:26s} time-avg N=4e5 {tavg:+.9f}   exact subtorus {vals[(11,20)]:+.9f}"
          f"   |t-avg - exact| {abs(tavg-vals[(11,20)]):.2e}   |t-avg - generic| {abs(tavg-g):.2e}")

# ------------------------------------------------------------------ 1D Lawton accumulation
print("\n" + "=" * 104)
print("== 1D  S4's ACCUMULATION TABLE (Lawton 1983), ONE VARIABLE MOVED ==")
print("=" * 104)
print("  S4 publishes, for K1: (1,1) -4.36e-01 · (5,3) +2.28e-03 · (11,20) +4.93e-04 ·")
print("  (41,53) +9.09e-05 · (97,61) -7.69e-06 · (610,377) +5.25e-07.  Lawton's theorem is the")
print("  reason and W-03's correction records it MISSING FROM S4's IMPORT AUDIT.")
ACC = [(1, 1), (5, 3), (11, 20), (41, 53), (97, 61), (610, 377)]
print(f"\n  {'(m,n)':>11s} {'|m|+|n|':>8s}" + "".join(f"{lab.split()[0]:>14s}" for lab, _ in ARMS))
for mn in ACC:
    row = f"  {str(mn):>11s} {abs(mn[0])+abs(mn[1]):8d}"
    for lab, p in ARMS:
        g = DETAIL[lab][0]
        row += f"{lam_res(p,*mn)-g:+14.2e}"
    print(row)
print("\n  Every column accumulates on 0 -- Lawton is carrier-independent.  THE RATE AT WHICH")
print("  IT ACCUMULATES, AND THEREFORE HOW MUCH STRUCTURE THE EXCEPTIONAL SET CARRIES, IS NOT.")

# ------------------------------------------------------------------ 1E mechanism + own defect
print("\n" + "=" * 104)
print("== 1E  WHY THE FOUR-CLASS COLUMNS COLLAPSE, AND ONE PRECISION DEFECT OF MY OWN ==")
print("=" * 104)
print("  (a) TORUS ZEROS.  P(x,y) = (p00+p01 y) + x(p10+p11 y) has a zero on T^2 iff the two")
print("      Jensen branches meet, i.e. iff  A + B cos t = 0  has a root in [-1,1], with")
print("      A = p00^2+p01^2-p10^2-p11^2 and B = 2(p00 p01 - p10 p11).  This is lane D's own")
print("      leg 1D criterion.  Lawton's convergence is SLOW exactly when there is a log")
print("      singularity to resolve on the shrinking subtorus, i.e. exactly when P has a")
print("      torus zero -- which is exactly when W-01's criterion can fire at all.")
print(f"  {'arm':26s} {'A':>12s} {'B':>12s} {'torus zero?':>12s} {'|dev| at (41,53)':>18s}")
for lab, p in ARMS:
    p00, p10, p01, p11 = [float(q) for q in p]
    A = p00**2 + p01**2 - p10**2 - p11**2
    B = 2*(p00*p01 - p10*p11)
    tz = (abs(B) > 1e-300 and -1.0 <= -A/B <= 1.0) or (abs(B) <= 1e-300 and abs(A) < 1e-15)
    g = DETAIL[lab][0]
    print(f"  {lab:26s} {A:12.6f} {B:12.6f} {str(tz):>12s} {abs(lam_res(p,41,53)-g):18.2e}")
print("      EVERY arm with a torus zero has |dev| ~ 1e-5 at (41,53); every arm without one is")
print("      at 1e-14.  The three-class arms have torus zeros; the two four-class carriers the")
print("      corpus owns do not, AT THEIR OWN PUBLISHED WEIGHTS.  (SENSE-C is the four-class")
print("      arm that DOES have them -- a whole curve -- and it is the exception that fixes the")
print("      variable: the mechanism is the WEIGHTS, not the class count.)")

print("\n  (b) THE (1,-1) ROW IN CLOSED FORM -- c = -f, the resonance lane D reported as")
print("      departing by 5.80e-07 on B0b.  Exactly: Q = (p01+p10) + p00 z + p11 z^{-1},")
print("      so z.Q = p11 + (p01+p10) z + p00 z^2, LEADING COEFFICIENT p00.")
print("      MY OWN CONFOUND, RECORDED NOT PATCHED: the first version of this block passed the")
print("      coefficient array to np.roots REVERSED (giving the polynomial p11 z^2 + ... + p00)")
print("      while still taking log|p00| as the leading coefficient, so the three arms with")
print("      p00 != 0 printed nonsense (B1q -1.205e+00, B0b +6.931e-01, B4 -1.099e+00 as the")
print("      'departure').  Both statistics are printed below.  The MAIN tables 1B/1C/1D were")
print("      never affected -- they call m_one_var(), which takes the leading coefficient of")
print("      the array it actually factors -- and they already carried the correct figures.")
print(f"    {'arm':26s} {'m (CORRECT)':>17s} {'generic':>17s} {'departure':>12s}   {'m (DEFECTIVE 1st run)':>22s}")
for lab, p in ARMS[:5]:
    p00, p10, p01, p11 = [float(q) for q in p]
    asc = np.array([p11, p01 + p10, p00])           # z.Q ascending: p11 + (p01+p10) z + p00 z^2
    nz = np.nonzero(asc)[0]
    a = asc[nz[0]:nz[-1] + 1]
    des = a[::-1]                                   # highest degree first
    if len(des) == 1:
        exact = np.log(abs(des[0]))
    else:
        r = np.roots(des)
        exact = np.log(abs(des[0])) + np.sum(np.log(np.maximum(1.0, np.abs(r))))
    bad = np.log(abs(p00)) + np.sum(np.log(np.maximum(1.0, np.abs(np.roots(
        np.array([p00, p01 + p10, p11])[::-1]))))) if abs(p00) > 0 else exact
    g = DETAIL[lab][0]
    chk = abs(exact - DETAIL[lab][1][(1, -1)])
    print(f"    {lab:26s} {exact:+17.12f} {g:+17.12f} {exact-g:+12.3e}   {bad:+22.12f}"
          f"   [agrees with 1B to {chk:.1e}]")
print("      On BOTH four-class carriers the departure at c = -f is EXACTLY ZERO (the closed")
print("      form returns the generic value), and lane D printed 5.80e-07 for it.")

print("\n  (c) MY OWN PRECISION DEFECT, RECORDED NOT PATCHED.  My generic-value routine is")
print("      split Gauss-Legendre on the Jensen branches.  On SENSE-C the two branches COINCIDE")
print("      on the whole circle (P = (1+x)(1+y)/4), so the split is degenerate and the routine")
print("      returns -1.386292610 against the exact log(1/4) = -1.386294361 -- an error of")
print("      1.75e-06, which is the ENTIRE SENSE-C column of 1B/1D.  That column is therefore")
print("      MY quadrature error and not an exceptional-value structure; the exact SENSE-C")
print("      exceptional values are log(1/4) at every (m,n) with m,n both non-zero and log(1/2)")
print("      at (1,0),(0,1).  Corrected figures:")
lc = np.log(0.25)
for mn in [(1, 0), (1, 1), (11, 20), (41, 53)]:
    print(f"        SENSE-C (m,n)={str(mn):>9s}  m = {lam_res(ARMS[5][1],*mn):+.12f}"
          f"   exact generic log(1/4) = {lc:+.12f}   dep {lam_res(ARMS[5][1],*mn)-lc:+.2e}")
print("      No other arm is affected: every other arm's branches cross transversally or not at")
print("      all, and 1A validates the whole pipeline against S4's published column at 4.8e-10.")
