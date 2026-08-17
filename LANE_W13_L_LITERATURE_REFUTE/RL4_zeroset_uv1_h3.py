#!/usr/bin/env python3
"""
RL4 — THE ZERO SET: WHAT L-17 ACTUALLY SHOWS, WHAT THE BRIEF'S BICONDITIONAL DOES NOT,
      AND A REDUCTION OF H3 THAT THE TARGET MISSED AT ITS OWN pi.

Three EXACT computations.  These are not controls and are not scored as evidence for
anything; they are exact facts in fractions.Fraction / Python integers.

(A) THE BRIEF'S HEADLINE BICONDITIONAL IS FALSE OFF THE THREE-CLASS LOCUS.
    The commissioning brief states: "So P has a torus zero IF AND ONLY IF W-01's criterion
    fires", identifying the firing criterion with max(pi) <= sum of the others, i.e.
    max(pi) <= 1/2.  On the THREE-class locus p00 = 0 that is a theorem (M1_02 Z1, sealed,
    and re-derived exactly below).  On the FOUR-class locus it is FALSE, and W-09 already
    registered the reason in another guise: with all four classes occupied, the fourth
    character uv is DETERMINED by u and v, so the four coefficients are not four free unit
    vectors and Wendel/hull reasoning does not transfer.  EXHIBITED below with an exact
    counterexample.  Neither the target lane nor the brief carries this.

(B) L-17 ("W-08's uv=1 quadratic and N1's zero set are ONE object seen twice") IS TRUE AT
    K1's REGISTERED pi AND IS NOT GENERAL.  A torus zero (x,y) lies on {xy = 1} iff x is a
    root of  q(x) = p10 x^2 + (p00+p11) x + p01  of modulus 1.  BOTH roots have modulus 1
    iff |p01/p10| = 1, i.e. p10 = p01 -- the loop-symmetric case, which K1's registered pi
    satisfies.  Otherwise it can still happen degenerately (x = -1 is a root when
    p10 - p00 - p11 + p01 = 0).  Exhibited below on both sides.

(C) THE THING THE TARGET MISSED, AND IT CUTS ITS OWN WAY.  Because at K1's registered pi
    the whole zero set lies inside {xy = 1}, and (x,y) |-> xy is 1-Lipschitz from T^2 (with
    the sup metric) to T,
            dist( (u^k,v^k), Z(P) )  >=  c * | (uv)^k - 1 |.
    So at the corpus's own registered ready state, H3 -- the INHOMOGENEOUS condition the
    target proposes and names the operative variable -- FOLLOWS FROM A HOMOGENEOUS
    condition on the SINGLE character uv = conj(W_F) W_C.  H3 at K1's pi is a one-dimensional
    Diophantine statement about one number, not a two-dimensional inhomogeneous one.
    That STRENGTHENS the target's L-15 and NARROWS its L-16.
"""
from fractions import Fraction as F
import math
import numpy as np

print("=" * 78)
print("RL4 — ZERO SET: THE BICONDITIONAL, L-17's SCOPE, AND A REDUCTION OF H3.")
print("=" * 78)


def has_zero_exact(p):
    """EXACT.  P = (p00+p10 x) + (p01+p11 x) y vanishes for some (x,y) in T^2 iff
       F(c) := (p00^2+p10^2+2 p00 p10 c) - (p01^2+p11^2+2 p01 p11 c) has a root in [-1,1].
       F is affine in c, so: iff F(1) and F(-1) have opposite signs or one is zero."""
    p00, p10, p01, p11 = p
    F1 = (p00 + p10) ** 2 - (p01 + p11) ** 2
    Fm = (p00 - p10) ** 2 - (p01 - p11) ** 2
    return (F1 == 0) or (Fm == 0) or ((F1 > 0) != (Fm > 0))


def t_of(p):
    p00, p10, p01, p11 = p
    den = 2 * (p00 * p10 - p01 * p11)
    if den == 0:
        return None
    return (p01 * p01 + p11 * p11 - p00 * p00 - p10 * p10) / den


print("""
--------------------------------------------------------------------------------
(A) "TORUS ZERO  <=>  max(pi) <= 1/2".  EXACT, ON BOTH LOCI.
--------------------------------------------------------------------------------""")
print("  %-34s %-10s %-14s %-14s %s" % ("pi", "max(pi)", "max<=1/2?", "torus zero?", "agree?"))
rows = [
    ("K1 registered", (F(0), F(3, 10), F(3, 10), F(2, 5))),
    ("three-class centroid", (F(0), F(1, 3), F(1, 3), F(1, 3))),
    ("three-class, max>1/2", (F(0), F(1, 10), F(1, 10), F(8, 10))),
    ("three-class edge", (F(0), F(1, 2), F(1, 4), F(1, 4))),
    ("FOUR-CLASS counterexample", (F(1, 2), F(3, 10), F(1, 10), F(1, 10))),
    ("four-class uniform", (F(1, 4), F(1, 4), F(1, 4), F(1, 4))),
    ("four-class B0b-like", (F(2, 9), F(4, 9), F(3, 9), F(0))),
    ("FOUR-CLASS counterexample 2", (F(9, 20), F(7, 20), F(3, 20), F(1, 20))),
]
bad = []
for name, p in rows:
    mx = max(p)
    le = mx <= F(1, 2)
    hz = has_zero_exact(p)
    ok = (le == hz)
    if not ok:
        bad.append((name, p))
    print("  %-34s %-10s %-14s %-14s %s" % (str(tuple(str(q) for q in p)), str(mx),
                                            "YES" if le else "no", "YES" if hz else "no",
                                            "yes" if ok else "*** NO ***"))
print("\n  DISAGREEMENTS: %d" % len(bad))
for name, p in bad:
    p00, p10, p01, p11 = [float(q) for q in p]
    print("     %s : pi = %s" % (name, tuple(str(q) for q in p)))
    print("       |p00 + p10 x| ranges over [%.4f, %.4f];  |p01 + p11 x| over [%.4f, %.4f]"
          % (abs(p00 - p10), p00 + p10, abs(p01 - p11), p01 + p11))
    print("       both are AFFINE in cos(arg x) with the SAME argument, so equality needs")
    print("       cos(arg x) = t = %s, and |t| = %s > 1.  NO TORUS ZERO, while max(pi) <= 1/2."
          % (t_of(p), abs(t_of(p))))
    # brute-force confirmation on a fine grid, float64
    n = 4096
    th = np.arange(n) * (2 * np.pi / n)
    ex = np.exp(1j * th)
    A = np.abs(p00 + p10 * ex)[:, None]
    B = np.abs(p01 + p11 * ex)[:, None]
    mn = float(np.min(np.abs(A - B)))
    tot = np.abs((p00 + p10 * ex)[:, None] + (p01 + p11 * ex)[:, None] * np.exp(1j * th)[None, :])
    print("       brute force 4096^2 grid: min |P| = %.6e   (min over x of ||a(x)|-|b(x)|| = %.6e)"
          % (float(tot.min()), mn))
print("""
  READING.  On the THREE-class locus p00 = 0 the biconditional is a theorem: p10 <= p01+p11
  is max <= 1/2, and p10 >= |p01-p11| is also exactly max(p01,p11) <= 1/2.  On the FOUR-class
  locus it FAILS, and the failure is not exotic -- it happens at pi with max(pi) EQUAL to the
  registered threshold.  THE BRIEF'S SENTENCE "P has a torus zero if and only if W-01's
  criterion fires" IS THREE-CLASS-SCOPED, exactly like every other W-01 statement W-09
  rescoped.  Nothing in the target lane depends on it -- every leg runs at p00 = 0 -- but the
  target repeats the brief's framing without the qualification, and this is the fifth time in
  this program that a three-class result has been quoted forward unqualified.""")

print("""
--------------------------------------------------------------------------------
(B) L-17's SCOPE.  WHEN DO THE TORUS ZEROS LIE ON THE SUBTORUS {xy = 1}?
    y = 1/x requires q(x) = p10 x^2 + (p00+p11) x + p01 = 0 with |x| = 1.  The product of q's
    roots is p01/p10, so BOTH roots lie on |x|=1 iff p10 = p01.  A degenerate single root at
    x = -1 (when p10 - p00 - p11 + p01 = 0) also lies there -- row 8 below is that case, and
    it is NOT a counterexample to the analysis, only to the naive "iff".  EXACT test:
--------------------------------------------------------------------------------""")
print("  %-36s %-12s %-14s %s" % ("pi", "p10 == p01?", "zeros exist?", "max |x*y - 1| over Z(P)"))
for name, p in rows:
    if not has_zero_exact(p):
        print("  %-36s %-12s %-14s -" % (str(tuple(str(q) for q in p)),
                                         "YES" if p[1] == p[2] else "no", "no"))
        continue
    t = t_of(p)
    if t is None:
        print("  %-36s %-12s %-14s degenerate (zero CIRCLE)" % (str(tuple(str(q) for q in p)),
                                                                "YES" if p[1] == p[2] else "no", "yes"))
        continue
    ct = float(t)
    st = math.sqrt(max(0.0, 1 - ct * ct))
    worst = 0.0
    for s in (+1, -1):
        x = complex(ct, s * st)
        num = float(p[0]) + float(p[1]) * x
        de = float(p[2]) + float(p[3]) * x
        if abs(de) < 1e-13:
            continue
        y = -num / de
        worst = max(worst, abs(x * y - 1))
    print("  %-36s %-12s %-14s %.6e" % (str(tuple(str(q) for q in p)),
                                        "YES" if p[1] == p[2] else "no", "yes", worst))
print("""
  READING.  |xy - 1| = 0 to round-off on every row with p10 = p01, and equals 2 and 0.97 --
  i.e. O(1) -- on rows 4 and 7, which have p10 != p01 and honest two-point zero sets.  Row 8
  has p10 != p01 and t = -1 exactly, a DEGENERATE single zero at (x,y) = (-1,-1), which lies
  on {xy=1} for a different reason; it is recorded rather than hidden.
  So L-17's sentence "P's zeros on T^2 are (x, xbar) with x + xbar = -4/3, which is exactly
  that quadratic" is a fact about pi with p10 = p01 -- which K1's registered state satisfies
  -- and NOT a general fact about P.  L-17's CONCLUSION at K1's pi is correct; its stated
  generality is not.  W-08's uv=1 locus and N1's zero set coincide there because the
  registered ready state is SYMMETRIC IN THE TWO LOOPS, not because the two phenomena are
  one object.  On a carrier whose designated loops carry different weight they separate.""")

print("""
--------------------------------------------------------------------------------
(C) THE REDUCTION THE TARGET MISSED.  H3 AT K1's REGISTERED pi IS HOMOGENEOUS.
--------------------------------------------------------------------------------""")
p = (0.0, 0.3, 0.3, 0.4)
x0 = complex(-2.0 / 3.0, math.sqrt(5.0) / 3.0)
Z0 = [(x0, x0.conjugate()), (x0.conjugate(), x0)]
print("  Z(P) = {(x0, conj x0), (conj x0, x0)},  and x0 * conj(x0) = |x0|^2 = 1 EXACTLY.")
print("  So Z(P) is contained in the subtorus S = {(x,y) in T^2 : xy = 1}.")
print("  phi(x,y) = xy is 1-Lipschitz for the sup-metric on T^2 (|x1y1 - x2y2| <= |x1-x2|+|y1-y2|),")
print("  and phi == 1 on Z(P).  Hence for every k")
print("        |(uv)^k - 1| = |phi(u^k,v^k) - phi(z)| <= 2 * dist_inf((u^k,v^k), Z(P)).")
print("  i.e.  dist((u^k,v^k), Z(P)) >= (1/2)|(uv)^k - 1|.   H3 FOLLOWS FROM A BOUND ON ONE NUMBER.")
print("\n  NUMERICAL CHECK OF THE INEQUALITY (it is an inequality, so it CAN fail; it does not):")
rng_a = [("arm A algebraic uv = (-33+56i)/65",
          (math.atan2(4, 3), math.atan2(12, 5))),
         ("arm B transcendental f=1, c=sqrt2", (-1.0, math.sqrt(2.0)))]
for lab, (fa, fb) in rng_a:
    ks = np.arange(1, 200001, dtype=np.float64)
    tu = np.mod(ks * (fa / (2 * np.pi)), 1.0)
    tv = np.mod(ks * (fb / (2 * np.pi)), 1.0)
    eu = np.exp(2j * np.pi * tu)
    ev = np.exp(2j * np.pi * tv)
    d = np.full(ks.shape, np.inf)
    for (zx, zy) in Z0:
        d = np.minimum(d, np.maximum(np.abs(eu - zx), np.abs(ev - zy)))
    lhs = np.abs(eu * ev - 1.0)
    ratio = lhs / (2 * d)
    print("   %-38s max over k<=2e5 of |(uv)^k-1| / (2 dist) = %.6f   (must be <= 1)  min dist = %.3e"
          % (lab, float(ratio.max()), float(d.min())))
print("""
  CONSEQUENCE, AND IT CUTS BOTH WAYS.
  * IT STRENGTHENS L-15.  At the corpus's own registered ready state, H3 is implied by the
    single homogeneous condition  || k * arg(uv)/2pi || >= c k^{-tau}.  That is a statement
    about ONE character of the connection, checkable, and supplied by Baker whenever uv is
    algebraic and not a root of unity -- which is a weaker hypothesis than "u and v algebraic
    and multiplicatively independent".  THEOREM L's hypotheses can be weakened at this pi.
  * IT NARROWS L-16.  The target names the operative variable "the INHOMOGENEOUS Diophantine
    type of the connection relative to the zero set of P".  That is right in general and is
    NOT the minimal name at the corpus's own registered state, where the inhomogeneous
    condition collapses to a homogeneous condition on the product character uv.  Since every
    number the corpus publishes for N1 is computed at that pi, the operative variable AS
    EXERCISED BY THE CORPUS is one degree simpler than the name the target attached.
  * AND IT DOES NOT RESCUE ARM B.  arg(uv)/2pi = (sqrt2 - 1)/(2 pi) is still a linear form in
    1, sqrt2 and pi, so no effective bound is available there either.  L-4 STANDS.
--------------------------------------------------------------------------------""")
print("\nDONE RL4")
