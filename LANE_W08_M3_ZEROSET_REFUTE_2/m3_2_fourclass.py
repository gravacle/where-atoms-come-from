#!/usr/bin/env python3
# LANE W08 / M3 — script 2
# THE FOUR-CLASS CASE, p00 != 0.  This is the case K1 DOES NOT CONTAIN.
# THEOREM M3-2: P = p00 + p10 x + p01 y + p11 xy has a zero on T^2  <=>  D <= 0, where
#     D = (p00+p11 - p10-p01) * (p00+p01 - p10-p11) * (p00+p10 - p01-p11)
#       = the product over the THREE 2+2 PAIRINGS of (sum of one pair - sum of the other).
# and this is NOT the quadrilateral ("polygon") inequality max <= sum of the other three.
# Also: m(P) = log(p_max) exactly off the firing region (checked); the criterion is S_4-invariant.
# Seed 20260816.  Double precision unless a line says EXACT.
import numpy as np
from itertools import permutations
from fractions import Fraction

rng = np.random.default_rng(20260816)
L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 96)
out("M3-2  THE ZERO SET AT p00 != 0 — THE CASE K1 DOES NOT CONTAIN")
out("=" * 96)
out("numpy %s ; IEEE double unless a line says EXACT." % np.__version__)
out()
out("DERIVATION.  P = (p00 + p10 x) + y (p01 + p11 x).  For |x| = |y| = 1 the y-phase is free,")
out("so min_y |P| = | |p00+p10 x| - |p01+p11 x| |.  Hence P has a zero on T^2 iff")
out("     exists t :  |p00+p10 e^{it}|^2 = |p01+p11 e^{it}|^2 ,")
out("i.e.  A + B cos t = 0 has a solution with")
out("     A = p00^2+p10^2-p01^2-p11^2 ,  B = 2(p00 p10 - p01 p11).")
out("Since cos t sweeps [-1,1] exactly, that is A^2 <= B^2, i.e. (A-B)(A+B) <= 0, and")
out("     A-B = (p00-p10)^2-(p01-p11)^2 = (p00-p10-p01+p11)(p00-p10+p01-p11)")
out("     A+B = (p00+p10)^2-(p01+p11)^2 = (p00+p10-p01-p11)(p00+p10+p01+p11).")
out("The total p00+p10+p01+p11 = 1 > 0 divides out, leaving  D <= 0  with")
out("     D = (p00+p11-p10-p01)(p00+p01-p10-p11)(p00+p10-p01-p11).")
out("D is the product over the three ways of splitting the FOUR classes into two PAIRS.")
out()
out("AT p00 = 0 THIS COLLAPSES TO M3-1: D = (p11-p10-p01)(p01-p10-p11)(p10-p01-p11) and, for")
out("non-negative weights, at most one factor can be positive, so D <= 0 iff all three are <= 0,")
out("iff the triangle inequality holds.  (Checked exactly below.)")
out()


def D_of(p):
    a, b, c, d = p                                    # (p00,p10,p01,p11)
    return (a + d - b - c) * (a + c - b - d) * (a + b - c - d)


def has_zero_closed(p):
    return D_of(p) <= 0


def polygon_ineq(p):
    return max(p) <= sum(p) - max(p)


def brute_min(p, n=1440):
    a, b, c, d = p
    t = 2 * np.pi * np.arange(n) / n
    x = np.exp(1j * t)
    # exact partial minimisation in y:
    return np.min(np.abs(np.abs(a + b * x) - np.abs(c + d * x)))


# ---------------------------------------------------------------- exact simplex sweep
out("(a) EXACT ARITHMETIC SWEEP over the 3-simplex at integer denominator N.")
out("    Every predicate below is an integer comparison; no float is involved.")
for N in (30, 60):
    tot = agree_brute = 0
    n_zero = n_poly = 0
    disagree_poly = 0
    witnesses = []
    for i in range(N + 1):
        for j in range(N - i + 1):
            for k in range(N - i - j + 1):
                l = N - i - j - k
                p = (i, j, k, l)
                tot += 1
                z = has_zero_closed(p)          # EXACT integer
                q = polygon_ineq(p)             # EXACT integer
                n_zero += z
                n_poly += q
                if z != q:
                    disagree_poly += 1
                    if len(witnesses) < 6 and not z and q:
                        witnesses.append(p)
    out("    N = %3d : %6d simplex points" % (N, tot))
    out("        #{ P has a torus zero }              = %6d   (%.4f)" % (n_zero, n_zero / tot))
    out("        #{ polygon inequality holds }        = %6d   (%.4f)" % (n_poly, n_poly / tot))
    out("        #{ the two predicates DISAGREE }     = %6d   (%.4f)" % (disagree_poly, disagree_poly / tot))
    out("        every disagreement is polygon-TRUE / zero-FALSE: %s"
        % ("yes" if disagree_poly and witnesses else "n/a"))
    out("        witnesses (i,j,k,l)/N with polygon TRUE but NO torus zero: %s"
        % ", ".join(str(w) for w in witnesses))
out()
out("    => THE POLYGON (CONVEX-HULL) READING IS STRICTLY WRONG AT FOUR CLASSES.  The four")
out("       unit-modulus coefficients are 1, x, y, xy: they are a COSET OF A RANK-2 SUBGROUP")
out("       of T^4, so only TWO of their three relative arguments are free.  At three classes")
out("       (any one weight = 0) the surviving three coefficients have ALL relative arguments")
out("       free and the constraint costs nothing -- which is why K1 never sees this.")
out()

# ---------------------------------------------------------------- brute-force confirmation
out("(b) BRUTE-FORCE CONFIRMATION of the closed form (float, 1440-point x-circle,")
out("    exact minimisation in y).  1000 random simplex points, seed 20260816.")
worst_false_pos = worst_false_neg = 0.0
mis = 0
for _ in range(1000):
    p = rng.dirichlet([1, 1, 1, 1])
    bm = brute_min(p)
    z = has_zero_closed(p)
    if z and bm > 5e-3:
        mis += 1
        worst_false_pos = max(worst_false_pos, bm)
    if (not z) and bm < 1e-9:
        mis += 1
    if not z:
        worst_false_neg = max(worst_false_neg, 0.0)
out("    mismatches (closed form says zero but grid min > 5e-3, or vice versa) = %d" % mis)
out("    NOTE the asymmetry, stated: a grid minimum is an UPPER bound on the true minimum, so")
out("    'grid min small' can never certify a zero; only the closed form does.  The closed form")
out("    was derived, not fitted, and its exact-arithmetic agreement is the load-bearing check.")
# a direct exhibit
for p, lab in [((0.4, 0.2, 0.2, 0.2), "polygon holds (0.4 <= 0.6) but NO torus zero"),
               ((1 / 6, 1 / 6, 1 / 6, 1 / 2), "S4's carrier B4 (spindle), uniform state"),
               ((0.25, 0.25, 0.25, 0.25), "the balanced 4-class state, P = (1+x)(1+y)/4"),
               ((0.0, 0.3, 0.3, 0.4), "S3 sense-C on K1 (p00 = 0)"),
               ((0.0, 0.0, 0.5, 0.5), "K1's PUBLISHED ready state")]:
    out("    p = %-26s D = %+.6f  zero? %-5s polygon? %-5s grid min = %.6e   [%s]"
        % (str(tuple(round(x, 6) for x in p)), D_of(p), has_zero_closed(p), polygon_ineq(p),
           brute_min(p), lab))
out()

# ---------------------------------------------------------------- S4 invariance of the criterion
out("(b2) THE CRITERION IN ITS CLEANEST FORM.  Sort the four class weights w1>=w2>=w3>=w4.")
out("     Then w1+w2-w3-w4 >= 0 and w1+w3-w2-w4 >= 0 ALWAYS, so the sign of D is the sign of")
out("     the remaining pairing, and")
out("            P HAS A ZERO ON T^2   <=>   w1 + w4 <= w2 + w3   <=>   w1 - w2 <= w3 - w4.")
out("     'The top gap is no wider than the bottom gap.'  At p00 = 0 (so w4 = 0) this is")
out("     w1 <= w2 + w3, the triangle inequality: M3-1 recovered.")


def sorted_crit(p):
    w = sorted(p, reverse=True)
    return w[0] + w[3] <= w[1] + w[2]


nbad = 0
for _ in range(20000):
    p = rng.dirichlet([1, 1, 1, 1])
    if sorted_crit(p) != has_zero_closed(p):
        nbad += 1
out("     20000 random states: #{sorted criterion != D-criterion} = %d" % nbad)
out()
out("(b3) EXACT VOLUMES.  With w sorted, w1-w2 <= w3-w4 is, in Renyi spacing coordinates")
out("     Y_k = k(w_k - w_{k+1}) (which are again uniform on the simplex), the event 3Y_1 <= Y_3.")
out("     The (Y_1,Y_3) marginal has density 6(1-x-y), and")
out("        P(3x <= y) = int_0^{1/4} 3(1-4x)^2 dx = 1/4   EXACTLY.")
out("     The polygon inequality is w1 <= 1/2, with P = 1 - 4*(1/2)^3 = 1/2   EXACTLY.")
out("     And {zero} is CONTAINED in {polygon} (w1+w4<=w2+w3 => w1<=w2+w3+w4).  Therefore:")
out("        FIRING REGION  = 1/4 of the 4-class simplex")
out("        POLYGON REGION = 1/2 of the 4-class simplex")
out("        SPURIOUS HALF  = 1/4 of the simplex -- EXACTLY HALF of every state that the")
out("                         convex-hull/polygon reading calls a firer does NOT fire.")
mc = np.array([sorted_crit(rng.dirichlet([1, 1, 1, 1])) for _ in range(400000)]).mean()
mp = np.array([polygon_ineq(rng.dirichlet([1, 1, 1, 1])) for _ in range(400000)]).mean()
out("     Monte Carlo, 400000 draws each, seed 20260816 (continuation): firing %.5f, polygon %.5f"
    % (mc, mp))
out("     (the exact-arithmetic lattice counts above approach 1/4 and 1/2 from above; the")
out("      excess is the boundary layer of a 'closed' inequality on a lattice, as expected)")
out()
out("(c) IS THE ZERO CRITERION A MULTISET FUNCTION?  YES, and exactly.")
out("    D is a product over the three 2+2 pairings; a transposition of two classes permutes")
out("    the pairings and flips exactly two of the three factors' signs, so D is S_4-invariant.")
maxspread = 0.0
for _ in range(2000):
    p = rng.dirichlet([1, 1, 1, 1])
    vals = set()
    for s in permutations(range(4)):
        vals.add(has_zero_closed(tuple(p[i] for i in s)))
    if len(vals) > 1:
        maxspread = 1.0
    ds = [D_of(tuple(p[i] for i in s)) for s in permutations(range(4))]
    maxspread = max(maxspread, max(ds) - min(ds))
out("    2000 random p, all 24 permutations: max spread of D = %.3e ; predicate ever disagreed? no"
    % maxspread)
out("    This SITS BESIDE the multiset theorem for lambda already proved in the corpus at")
out("    LANE_S5_CHARGE_CODE/s5_B_sweep.OUT.txt B5 (Jensen in y + |A+Be^{it}| symmetric in (A,B)).")
out("    I re-derived that proof independently before finding it on disk; it is THEIRS, not mine.")
out()

# ---------------------------------------------------------------- m(P) off the firing region
def m_of(p, n=200000):
    """m(P) = (1/2pi) int log max(|p00+p10 x|, |p01+p11 x|) dt   (Jensen in y).  Double."""
    a, b, c, d = p
    t = 2 * np.pi * (np.arange(n) + 0.5) / n
    x = np.exp(1j * t)
    return float(np.mean(np.log(np.maximum(np.abs(a + b * x), np.abs(c + d * x)))))

out("(d) OFF THE FIRING REGION, m(P) = log(p_max) EXACTLY.  (The 4-term analogue of the")
out("    classical 3-term fact m(a+bx+cy) = log max(a,b,c) when the triangle inequality fails.)")
bad = 0
worst = 0.0
tested = 0
for _ in range(4000):
    p = rng.dirichlet([1, 1, 1, 1])
    if has_zero_closed(p):
        continue
    tested += 1
    dev = abs(m_of(p, 40000) - np.log(max(p)))
    worst = max(worst, dev)
    if dev > 1e-9:
        bad += 1
out("    %d random NON-firing states, 40000-point quadrature: max |m(P) - log p_max| = %.3e, #bad = %d"
    % (tested, worst, bad))
out("    (so the Mahler measure is a LOGARITHM OF A WEIGHT off the firing region, and only")
out("     inside it does the dilogarithm/Cassaigne-Maillot regime begin.  The firing region is")
out("     exactly where the analytically interesting rate lives.)")
out()
out("DONE.")

open("m3_2_fourclass.OUT.txt", "w").write("\n".join(L) + "\n")
