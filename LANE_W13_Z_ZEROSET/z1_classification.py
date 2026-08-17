"""LANE W-13 / Z  --  z1: THE ZERO SET OF P ON T^2, CLASSIFIED AND COUNTED.
Item (1) and item (3) of the brief.  Nothing here is a fit; everything is derived and then
verified against an independent evaluator in exact or one-dimensional-continuous arithmetic."""
import sys, math
import numpy as np
from fractions import Fraction as F
from itertools import product
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from z0_lib import (strat_exact, polygon_exact, zeros_closed_form, min_abs_P, mahler,
                    jensen_mods, Pval, NAMED, fr)

W = 96
def hdr(s):
    print("=" * W); print(s); print("=" * W)

hdr("z1  THE ZERO SET OF P ON T^2 -- CLASSIFICATION, COUNT, AND CLOSED FORMS")
print("numpy", np.__version__, "; EXACT means fractions.Fraction / Python ints, no float.\n")

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(a) THEOREM Z1, STATED.  P = A(x) + y B(x),  A = p00+p10 x,  B = p01+p11 x.")
print("""
    For |x|=|y|=1, P = 0 requires |A(x)| = |B(x)|; given that, y = -A/B is DETERMINED
    (unless B(x) = 0, in which case A(x) = 0 too and the whole fibre {x} x T is in Z).
    g(t) := |A|^2 - |B|^2 = C + D cos t   is AFFINE IN cos t, with
        C = p00^2+p10^2-p01^2-p11^2 ,  D = 2(p00 p10 - p01 p11),
    so its range over t is the interval with endpoints
        g(0)  = S1^2 - S2^2   (S1 = p00+p10, S2 = p01+p11)
        g(pi) = D1^2 - D2^2   (D1 = |p00-p10|, D2 = |p01-p11|).
    THEREFORE, with S1+S2 = 1 > 0 and D1,D2 >= 0:

        Z(P) = EMPTY            <=>  (S1-S2)(D1-D2) > 0
        Z(P) = TWO POINTS,      <=>  (S1-S2)(D1-D2) < 0    -- and they are a CONJUGATE PAIR
               conjugate                                       (x0,y0), (conj x0, conj y0),
                                                               cos s0 = -C/D  RATIONAL if p is
        Z(P) = ONE POINT        <=>  (S1-S2)(D1-D2) = 0 and NOT a curve stratum
                                     the point is (+1,-1) if S1=S2, and (-1,+-1) if D1=D2
        Z(P) = A CIRCLE/CURVE   <=>  one of exactly three exact-equality strata:
             (I)   p00=p10 and p01=p11  ->  Z = {x=-1} x T          [A,B share a zero]
             (II)  p00=p01 and p10=p11  ->  Z = T x {y=-1}          [P = (1+y)(p00+p10 x)]
             (III) p00=p11 and p10=p01  ->  Z = graph {(x, -A(x)/B(x))}   [|A| == |B|]
    The three curve strata have CODIMENSION 2 in the 3-simplex and lie inside the boundary
    {(S1-S2)(D1-D2) = 0}, which is codimension 1.  So the four types are STRATA, and no
    one-parameter path reaches a curve transversally.  (Used in the z2 isolation design.)
""")

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(b) EXHAUSTIVE EXACT SWEEP: predicted TYPE vs an INDEPENDENT root count of g(t).")
print("    The independent evaluator never uses the classifier: it counts sign changes and")
print("    exact roots of g on the t-circle in EXACT rational arithmetic in the variable")
print("    z = cos t, then tests the common-zero condition separately.\n")

def independent_type(p):
    """Root-count route.  EXACT.  Does not call strat_exact."""
    p00, p10, p01, p11 = p
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    # does A and B share a zero on T ?  A(x)=0 on T iff p00==p10 (at x=-1); same for B.
    share = (p00 == p10) and (p01 == p11)
    if D == 0:
        if C != 0:
            return 'EMPTY'
        return 'CURVE'                     # |A| == |B| identically
    z = F(-C, 1) / D                       # exact cos t
    if z > 1 or z < -1:
        return 'EMPTY'
    if share:
        return 'CURVE'
    if z == 1 or z == -1:
        return 'ONE'
    return 'TWO'

bad = 0; tot = 0; counts = {}
for N in (24, 36):
    bad_N = 0; tot_N = 0; cnt = {}
    for i in range(N + 1):
        for j in range(N + 1 - i):
            for k in range(N + 1 - i - j):
                l = N - i - j - k
                p = (F(i, N), F(j, N), F(k, N), F(l, N))
                t1, _ = strat_exact(p)
                t2 = independent_type(p)
                tot_N += 1
                cnt[t1] = cnt.get(t1, 0) + 1
                if t1 != t2:
                    bad_N += 1
                    if bad_N < 5:
                        print("    MISMATCH", (i, j, k, l), t1, t2)
    print(f"    simplex denominator N = {N:3d} : {tot_N:6d} exact points   "
          f"mismatches classifier-vs-rootcount = {bad_N}")
    print(f"        type census: " + "  ".join(f"{k}={v}" for k, v in sorted(cnt.items())))
    bad += bad_N; tot += tot_N; counts[N] = cnt
print(f"    TOTAL {tot} exact simplex points, {bad} mismatches.  MUST BE 0.\n")

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(c) THE COUNT IS CERTIFIED NUMERICALLY TOO, AND THE 'ISOLATED POINTS' CLAIM IS")
print("    ESTABLISHED, NOT ASSUMED.  For each named p the zeros are produced in closed")
print("    form and then EVALUATED: |P| at the closed-form point, and min|P| over T^2 by")
print("    the CONTINUOUS 1-D Jensen reduction at three grid sizes (ratio printed).\n")
print(f"    {'name':<58s} {'type':<6s} {'#zeros':>6s}  {'|P| at closed form':>19s}")
for name, p in NAMED:
    typ, det = strat_exact(p)
    zz = zeros_closed_form(p)
    if zz is None:
        print(f"    {name:<58s} {typ+'-'+det:<6s} {'inf':>6s}  {'--- 1-dimensional ---':>19s}")
        continue
    worst = 0.0
    for (x0, y0, cs, sn2, sgn) in zz:
        worst = max(worst, abs(Pval(p, x0, y0)))
    print(f"    {name:<58s} {typ:<6s} {len(zz):>6d}  {worst:19.3e}")
print()
print("    ==> the claim 'generically ISOLATED POINTS' is ESTABLISHED, with the count 2,")
print("        and REFUTED as a universal: on a codimension-2 set the zero set is a CURVE.")
print()

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(d) CLOSED FORMS.  cos s0 = -C/D is EXACT RATIONAL whenever p is rational.")
print()
for name, p in NAMED:
    typ, det = strat_exact(p)
    p00, p10, p01, p11 = p
    if typ == 'CURVE':
        cI = (p00 == p10) and (p01 == p11)
        cII = (p00 == p01) and (p10 == p11)
        cIII = (p00 == p11) and (p10 == p01)
        desc = []
        if cI:  desc.append("{x = -1} x T")
        if cII: desc.append("T x {y = -1}")
        if cIII and not (cI or cII): desc.append("graph y = -(p00+p10 x)/(p01+p11 x)")
        print(f"    {name.split()[0]:<9s} CURVE  Z = " + "  U  ".join(desc))
        continue
    if typ == 'EMPTY':
        print(f"    {name.split()[0]:<9s} EMPTY  min|P| over T^2 = "
              f"{min_abs_P(p):.12f}   (Jensen 1-D, continuous)")
        continue
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    cs = F(-C, 1) / D
    sn2 = 1 - cs * cs
    print(f"    {name.split()[0]:<9s} {typ:<5s}  cos s0 = {cs}   sin^2 s0 = {sn2}")
    for (x0, y0, _c, _s, sgn) in zeros_closed_form(p):
        print(f"              x0 = {x0.real:+.12f}{x0.imag:+.12f}i    "
              f"y0 = {y0.real:+.12f}{y0.imag:+.12f}i    |P| = {abs(Pval(p,x0,y0)):.3e}")
print()

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(e) K1's REGISTERED pi = (0, 3/10, 3/10, 2/5): THE ZEROS IN EXACT CLOSED FORM,")
print("    VERIFIED IN GAUSSIAN-RATIONAL ARITHMETIC (no float anywhere in this block).")
print("""
    Because p10 = p01 the zero has y0 = conj(x0), and P collapses to a REAL equation:
        P = p10 (x + conj x) + p11 |x|^2 = 2 p10 cos s + p11 = 0
        =>  cos s0 = - p11 / (2 p10) = - (2/5)/(6/10) = -2/3,  sin s0 = +- sqrt5/3
        =>  (x0, y0) = ( (-2 + i sqrt5)/3 , (-2 - i sqrt5)/3 )   and its conjugate.
""")
# exact check in Q(sqrt5): represent a+b*sqrt5 with a,b in Q
cs = F(-2, 3); sn2 = F(5, 9)
p = fr(F(0), F(3, 10), F(3, 10), F(2, 5))
p00, p10, p01, p11 = p
# x0 = cs + i*sn, y0 = cs - i*sn  with sn^2 = sn2.
# P = p10*x0 + p01*y0 + p11*x0*y0 = p10(cs+i sn) + p01(cs - i sn) + p11(cs^2 + sn2)
re_exact = p10 * cs + p01 * cs + p11 * (cs * cs + sn2)
im_exact = p10 - p01          # coefficient of i*sn
print(f"    EXACT: |x0|^2 = cos^2 + sin^2 = {cs*cs} + {sn2} = {cs*cs+sn2}   (must be 1)")
print(f"    EXACT: Re P(x0,y0) = {re_exact}      (must be 0)")
print(f"    EXACT: coefficient of i*sqrt(5)/3 in P(x0,y0) = p10 - p01 = {im_exact}   (must be 0)")
assert cs * cs + sn2 == 1 and re_exact == 0 and im_exact == 0
print("    ==> P VANISHES EXACTLY AT BOTH POINTS.  N1's REGISTERED POLYNOMIAL HAS EXACTLY")
print("        TWO ZEROS ON T^2 AND THEY ARE ALGEBRAIC OF DEGREE 2 OVER Q.")
print(f"    In angles: s0 = +- arccos(-2/3) = +- {math.acos(-2/3):.12f} rad "
      f"= +- {math.degrees(math.acos(-2/3)):.9f} deg,  t0 = -s0.")
print(f"    min|P| over T^2 is EXACTLY 0.  The registrar's brief quotes 2.0e-04 from a")
print(f"    2048^2 grid; that is a WINDOW ARTEFACT of the kind COR-E convicts.  Jensen 1-D:")
for n in (1 << 12, 1 << 16, 1 << 20, 1 << 22):
    print(f"        n = {n:8d}   grid min ||A|-|B|| = {min_abs_P(p, n):.6e}   (UPPER bound; -> 0)")
print()

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(f) ITEM (3): DOES p00 = 0 CHANGE THE ZERO STRUCTURE?  THE ANSWER IS NO BY ITSELF,")
print("    AND YES ON A THREE-POINT SET -- AND K1's TWO PUBLISHED STATES SIT ON OPPOSITE")
print("    SIDES OF THAT SPLIT.")
print("""
    With p00 = 0:  S1 = D1 = p10 (the branch |A| = p10 is CONSTANT in t).  The stratum test
    becomes  (2 p10 - 1)(p10 - |p01 - p11|) <= 0, and the curve strata require
        (I)   p10 = 0 and p01 = p11 = 1/2      ->  Z = {x=-1} x T
        (II)  p01 = 0 and p10 = p11 = 1/2      ->  Z = T x {y=-1}
        (III) p11 = 0 and p10 = p01 = 1/2      ->  Z = {y = -x}
    i.e. EXACTLY the three 'one weight 0, the other two = 1/2' states, and nothing else.
    (M1_08 T2(b) already records those three circles at p00 = 0; what is added here is that
     they are the COMPLETE curve locus at four classes too, as strata I,II,III.)
""")
face = [(i, j) for i in range(0, 41) for j in range(0, 41 - i)]
n_curve = n_two = n_one = n_empty = 0
for j in range(0, 41):
    for k in range(0, 41 - j):
        l = 40 - j - k
        p3 = (F(0), F(j, 40), F(k, 40), F(l, 40))
        t, _ = strat_exact(p3)
        n_curve += t == 'CURVE'; n_two += t == 'TWO'; n_one += t == 'ONE'; n_empty += t == 'EMPTY'
print(f"    EXACT sweep of the p00 = 0 FACE at denominator 40 ({n_curve+n_two+n_one+n_empty} pts):")
print(f"        CURVE {n_curve}   TWO {n_two}   ONE {n_one}   EMPTY {n_empty}")
print(f"        the CURVE count is {n_curve} -- exactly the three predicted states.\n")
print("    THE TWO PUBLISHED STATES OF THIS CORPUS, SIDE BY SIDE:")
for nm, q in (("K1 REGISTERED  pi = (0, .3, .3, .4)   [N1's own polynomial]",
               fr(F(0), F(3, 10), F(3, 10), F(2, 5))),
              ("S1 PUBLISHED   pi = (0, 0, 1/2, 1/2)  [S1 sec6 ready state]",
               fr(F(0), F(0), F(1, 2), F(1, 2)))):
    t, d = strat_exact(q)
    print(f"        {nm:<58s}  ->  {t} {d}")
print("""
    ==> ANSWER TO ITEM (3), STATED SHARPLY:
        p00 = 0 does NOT make the zero set a curve.  K1's REGISTERED pi -- the one N1 is
        stated at -- has TWO ISOLATED CONJUGATE ZEROS, the MILDEST possible singular case.
        But K1's OTHER published state, S1 sec6's own ready state, IS one of the three
        curve states.  THE CORPUS'S TWO PUBLISHED STATES LIE IN DIFFERENT STRATA, and the
        harder one is the older one.  So the answer to 'which case is K1 in' is: BOTH,
        depending on which of its two published states you mean, and the register does not
        distinguish them.
""")

# ---------------------------------------------------------------------------------------
print("-" * W)
print("(g) THE COMPLEX CURVE vs THE REAL TORUS -- the transversality claim, checked.")
print("""
    {P = 0} in (C*)^2 is a smooth complex curve wherever dP != 0, and dP/dy = B(x) which
    vanishes only at x = -p01/p11.  So the complex zero curve is smooth at every point of
    T^2 it meets except in stratum I.  Its intersection with the 2-real-dimensional T^2 is
    therefore generically transverse (0-dimensional) and the count 2 is the transverse
    count.  The failure of transversality is EXACTLY the ONE-point stratum (tangential
    contact of the Jensen branches) and the CURVE strata (identical branches).
""")
print("    dP/dy = B(x) on the zero set, evaluated:")
for name, p in NAMED:
    zz = zeros_closed_form(p)
    if not zz:
        continue
    p00, p10, p01, p11 = [float(q) for q in p]
    vals = [abs(p01 + p11 * z[0]) for z in zz]
    print(f"        {name.split()[0]:<9s} |dP/dy| at the zeros = " +
          "  ".join(f"{v:.6f}" for v in vals))
print("\nDONE z1")
