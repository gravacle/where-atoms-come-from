#!/usr/bin/env python3
# LANE W08 / M3 — script 1
# (a) The identity  Z_k = P(u^k, v^k)  checked against direct matrix action on K1.
# (b) THEOREM M3-1: with p00 = 0, P has a zero on T^2  <=>  triangle inequality on (p10,p01,p11).
#     Proved on paper; here checked (i) in EXACT INTEGER ARITHMETIC on an integer simplex grid
#     against the exact closed-form reduction, and (ii) numerically by brute-force torus search.
# (c) The QUANTIFIER: W-01/Theorem A is the CONNECTION-side projection; M3-1 is the STATE-side
#     projection. Exhibit showing the REGISTER's literal wording (quantifier dropped) is false.
# Seed: 20260816 (master).  Double precision unless stated.
import numpy as np
from fractions import Fraction

rng = np.random.default_rng(20260816)
L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 96)
out("M3-1  THE IDENTITY  Z_k = P(u^k, v^k)  ON K1, AND THE ZERO-SET THEOREM AT p00 = 0")
out("=" * 96)
out("numpy %s ; IEEE double unless a line says EXACT." % np.__version__)
out()

# ---------------------------------------------------------------- (a) the identity
FACE_V = [0, 1, 2]
CYC_V = [0, 3, 4]

def Z_direct(s, f, c, k):
    """<M_dF^k s, M_c^k s> by literal matrix action on C^5.  <z,w> = conj(z) w."""
    WF, WC = np.exp(1j * f), np.exp(1j * c)
    a = s.copy()
    b = s.copy()
    for v in FACE_V:
        a[v] *= WF ** k
    for v in CYC_V:
        b[v] *= WC ** k
    return np.vdot(a, b)          # numpy vdot conjugates the FIRST argument

def P_poly(p00, p10, p01, p11, x, y):
    return p00 + p10 * x + p01 * y + p11 * x * y

out("(a) IDENTITY CHECK.  s random complex normalised, (f,c) random, k = 1..12.")
worst = 0.0
for trial in range(400):
    s = rng.normal(size=5) + 1j * rng.normal(size=5)
    s /= np.linalg.norm(s)
    f, c = rng.uniform(0, 2 * np.pi, 2)
    p = np.abs(s) ** 2
    p11, p10, p01, p00 = p[0], p[1] + p[2], p[3] + p[4], 0.0
    u, vv = np.exp(-1j * f), np.exp(1j * c)
    for k in range(1, 13):
        d = abs(Z_direct(s, f, c, k) - P_poly(p00, p10, p01, p11, u ** k, vv ** k))
        worst = max(worst, d)
out("    400 states x 12 k-values, max |Z_k(direct) - P(u^k,v^k)| = %.3e" % worst)
out("    => Z_k IS the Mahler polynomial evaluated at the k-th orbit point.  Same object twice.")
out()

# ---------------------------------------------------------------- (b) exact reduction
# With p00 = 0:  P = x*(p10 + p11*y) + p01*y.  |x| = 1 is a free phase multiplying the whole
# first group, so min_{T^2}|P| = min_{|y|=1| abs( |p10 + p11 y| - p01 ).
# |p10+p11 y| sweeps the closed interval [ |p10-p11| , p10+p11 ].  Hence
#     zero exists  <=>  |p10-p11| <= p01 <= p10+p11   <=>  triangle inequality on the three.
out("(b) THEOREM M3-1 (p00 = 0).  P = x(p10 + p11 y) + p01 y.")
out("    min_{T^2}|P| = dist( p01 , [ |p10-p11| , p10+p11 ] )  -- an EXACT closed form.")
out("    Zero exists  <=>  |p10-p11| <= p01 <= p10+p11  <=>  each of the three <= sum of others.")
out("    With p10+p01+p11 = 1 this is  max(p10,p01,p11) <= 1/2.")
out()

def closed_form_min(p10, p01, p11):
    lo, hi = abs(p10 - p11), p10 + p11
    if p01 < lo:
        return lo - p01
    if p01 > hi:
        return p01 - hi
    return 0.0

def brute_min(p00, p10, p01, p11, n=720):
    t = 2 * np.pi * np.arange(n) / n
    x = np.exp(1j * t)[:, None]
    y = np.exp(1j * t)[None, :]
    return np.abs(P_poly(p00, p10, p01, p11, x, y)).min()

# EXACT check on an integer simplex grid: weights (i,j,k)/N, all predicates integer comparisons.
N = 200
mismatch_exact = 0
n_tri = 0
tot = 0
maxdev_closed = 0.0
for i in range(N + 1):
    for j in range(N - i + 1):
        kk = N - i - j
        tot += 1
        # EXACT integer predicate (no float anywhere in this line)
        tri = (i <= j + kk) and (j <= i + kk) and (kk <= i + j)
        # EXACT integer form of the closed-form reduction, (p10,p01,p11) = (i,j,kk)/N
        red = (abs(i - kk) <= j <= i + kk)
        if tri != red:
            mismatch_exact += 1
        n_tri += tri
out("    EXACT ARITHMETIC, integer simplex denominator N = %d, %d grid points:" % (N, tot))
out("      #{triangle inequality holds}                      = %d" % n_tri)
out("      #{triangle predicate != closed-form predicate}     = %d   <-- must be 0" % mismatch_exact)
out()

# numerical brute force, on a coarser simplex, with the undecidable band reported not hidden
N2 = 60
tol = 3e-3          # a 720x720 torus grid cannot resolve a conical zero below ~this
band = 0
bad = 0
pts = 0
maxdev = 0.0
for i in range(N2 + 1):
    for j in range(N2 - i + 1):
        kk = N2 - i - j
        p10, p01, p11 = i / N2, j / N2, kk / N2
        pts += 1
        cf = closed_form_min(p10, p01, p11)
        bm = brute_min(0.0, p10, p01, p11, n=360)
        maxdev = max(maxdev, abs(bm - cf) if cf > 0 else 0.0)
        tri = (i <= j + kk) and (j <= i + kk) and (kk <= i + j)
        if bm < tol and not tri:
            bad += 1
        if (not tri) and cf < tol:
            band += 1
out("    BRUTE FORCE, 360x360 torus grid, simplex denominator N = %d, %d points:" % (N2, pts))
out("      max | brute-grid min  -  closed-form min |  over the NON-firing region = %.3e" % maxdev)
out("      (the grid minimum is an UPPER bound on the true min; on the firing region")
out("       the true min is exactly 0 and the grid returns O(1/n) -- reported, not hidden)")
out("      #{grid min < %.0e while triangle inequality FAILS} = %d" % (tol, bad))
out("      #{triangle FAILS but true min < %.0e (undecidable band, |p_max-1/2| tiny)} = %d"
    % (tol, band))
out()

# how well does the grid minimum behave INSIDE the firing region (it should go to 0 like 1/n)
out("    Inside the firing region the zero is generically CONICAL; grid min ~ C/n:")
out("    *** CONFOUND RECORDED, NOT SILENTLY FIXED.  My first pass used n = 90,180,360,720,1440.")
out("    *** Every one is a multiple of 90, so the SAME lattice point recurred and the grid")
out("    *** minimum was CONSTANT at 1.478364e-03 for four successive n -- which reads as")
out("    *** 'the minimum is not going to zero' and is a pure grid artefact.  The refinement")
out("    *** ratio, not the value, is the observable.  Both ladders are printed.")
for n in (90, 180, 360, 720, 1440):
    out("      [multiples of 90]  p=(0.30,0.30,0.40)  n=%5d  grid min |P| = %.6e"
        % (n, brute_min(0, .3, .3, .4, n)))
for n in (97, 199, 401, 797, 1601):
    out("      [near-coprime   ]  p=(0.30,0.30,0.40)  n=%5d  grid min |P| = %.6e"
        % (n, brute_min(0, .3, .3, .4, n)))
out()

# ---------------------------------------------------------------- (c) the quantifier
out("(c) THE QUANTIFIER.  Two DIFFERENT projections of the one incidence variety")
out("      {(p,(f,c)) in Delta x T^2 : Z_1 = 0}.")
out("    W-01/Theorem A (S2 audit :246-252) fixes the CONNECTION and quantifies over STATES:")
out("      exists p in Delta with Z_1 = 0  <=>  0 in conv{ e^{i(c-f)}, e^{-if}, e^{ic} }")
out("                                       <=>  every consecutive angular gap <= pi.")
out("    THEOREM M3-1 fixes the STATE and quantifies over CONNECTIONS:")
out("      exists (f,c) with Z_1 = 0  <=>  P has a zero on T^2  <=>  triangle inequality on p.")
out("    THESE ARE NOT THE SAME PREDICATE; they have different arguments.")
out()
out("    REGISTER_V001.md:43 drops the quantifier: '...which vanishes IFF 0 lies in the convex")
out("    hull of three unit-modulus coefficients'.  Read literally (fixed p AND fixed (f,c))")
out("    that iff is FALSE.  Counterexample on K1's OWN published connection:")
f0, c0 = np.pi, 3 * np.pi / 2                     # W_F = -1, W_C = -i
u0, v0 = np.exp(-1j * f0), np.exp(1j * c0)
for name, (p10, p01, p11) in [("published p=(1/2,0,0,1/4,1/4)", (0.0, 0.5, 0.5)),
                              ("all weight at the root      ", (0.0, 0.0, 1.0)),
                              ("S3 sense-C p=(.4,.15,.15,.15,.15)", (0.3, 0.3, 0.4))]:
    z = P_poly(0.0, p10, p01, p11, u0, v0)
    out("      %-34s  Z_1 = %+.6f%+.6fi   |Z_1| = %.3e" % (name, z.real, z.imag, abs(z)))
out("      the three unit-modulus coefficients here are {uv, u, v} = {i, -1, -i}: angles")
out("      90, 180, 270 deg, gaps 90/90/180, all <= 180, so 0 IS in their convex hull.")
out("      Yet Z_1 != 0 for the root state.  The register's sentence needs '...CAN BE MADE TO")
out("      vanish by some ready state iff...'.  The SOURCE (S2 audit Theorem A) has it right.")
out()

# ---------------------------------------------------------------- corollaries
out("(d) COROLLARIES OF M3-1, each checked:")
out("    (i)  'the root can never fire' (W-01) = the extreme non-triangle p=(0,0,1)  ->  P = xy,")
out("         |P| == 1 on T^2, no zero.  min over 720x720 grid = %.6f" % brute_min(0, 0, 0, 1))
out("    (ii) firing region = {max <= 1/2} has simplex-area fraction 1/4 exactly")
out("         (the medial triangle).  EXACT count at denominator N=%d: %d of %d = %.6f"
    % (N, n_tri, tot, n_tri / tot))
out("         (-> 1/4 as N -> inf; the finite-N count exceeds 1/4 by the boundary layer)")
out("    (iii) DIMENSION OF THE ZERO SET.  |P| = 0 is two real equations in two real unknowns,")
out("         so the zero set is generically 0-dimensional.  It is 1-DIMENSIONAL exactly when P")
out("         FACTORS, which at p00 = 0 happens iff ONE weight is 0 and the other two are equal")
out("         (hence both = 1/2):  (0,1/2,1/2) -> P = (y/2)(1+x) ; (1/2,0,1/2) -> P = (x/2)(1+y) ;")
out("         (1/2,1/2,0) -> P = (1/2)(x+y).  Those are the three VERTICES of the medial")
out("         triangle, i.e. the three CORNERS of the firing region.  Everywhere else on the")
out("         boundary the two circles |p10+p11 y| and p01 are TANGENT: one degenerate zero point.")
out("         *** CORRECTION TO MY OWN FIRST PASS: I wrote '1-dimensional on the boundary of the")
out("         *** firing region'.  That is FALSE -- only at the three corners.  Recorded.")
for lab, (p10, p01, p11) in [("interior  (.30,.30,.40)", (.3, .3, .4)),
                             ("CORNER    (.00,.50,.50)", (0., .5, .5)),
                             ("bdy edge  (.25,.25,.50)", (.25, .25, .5)),
                             ("exterior  (.10,.20,.70)", (.1, .2, .7))]:
    n = 2000
    t = 2 * np.pi * np.arange(n) / n
    g = np.abs(p10 + p11 * np.exp(1j * t)) - p01     # zero in x exists at this y iff g = 0
    sign_changes = int(np.sum(np.sign(g[:-1]) * np.sign(g[1:]) < 0))
    zeros_exact = int(np.sum(np.abs(g) < 1e-14))
    out("         %s : sign changes of |p10+p11 y|-p01 = %4d ; |g|<1e-14 on %4d of %d grid y"
        % (lab, sign_changes, zeros_exact, n))
out("         (2 sign changes = 2 transverse conical zeros; g == 0 identically = the corner,")
out("          a whole circle of zeros -- its 'sign changes' count is pure float noise on 0;")
out("          1 touching point = the tangential boundary case; 0 and 0 = no zero at all)")
out()
out("DONE.")

open("m3_1_identity_and_threeclass.OUT.txt", "w").write("\n".join(L) + "\n")
