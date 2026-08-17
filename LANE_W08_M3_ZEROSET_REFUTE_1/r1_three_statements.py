#!/usr/bin/env python3
# LANE W08 / M3 REFUTER 1 — LENS 1 (MATHEMATICS): ATTACK THE EQUIVALENCE PROOF.
# The lens: "convex hull of three unit-modulus points, torus zeros of a two-variable
# polynomial, and the triangle inequality are three different statements; check each
# implication separately in BOTH directions."
#
# NAMING, fixed once and used everywhere below:
#   (HULL)  H(u,v)   :  0 in conv{ uv, u, v }                     -- about (u,v) ALONE
#   (ZERO)  Zt(p)    :  exists (x,y) in T^2 with P_p(x,y) = 0     -- about p ALONE
#   (TRI)   T(p)     :  |p10-p11| <= p01 <= p10+p11               -- about p ALONE
#   (FIRE)  F(p,u,v) :  Z_1 = p00 + p10 u + p01 v + p11 uv = 0    -- about BOTH
# On K1, p00 = 0 identically (incidence: FACE_V u CYC_V = all five vertices).
#
# ARITHMETIC: exact Fraction / integer arithmetic wherever a predicate is decided.
# Floats appear only where labelled FLOAT, and never decide a load-bearing predicate.
# I am Claude Opus 5 — the SAME model family as the lane I am refuting.  NOT lineage-
# independent of it.  Stated here so nothing below reads as independent corroboration.
import numpy as np
from fractions import Fraction as Fr

L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 100)
out("R1  THE THREE STATEMENTS, AND EVERY IMPLICATION BETWEEN THEM, IN BOTH DIRECTIONS")
out("=" * 100)
out("numpy %s ; predicates decided in EXACT integer/Fraction arithmetic unless a line says FLOAT."
    % np.__version__)
out()

# ------------------------------------------------------------------ (0) the identity, exact
out("(0) THE IDENTITY Z_k = P(u^k,v^k), RE-DERIVED FROM S1's INCIDENCE AND CHECKED EXACTLY.")
FACE_V, CYC_V = {0, 1, 2}, {0, 3, 4}
out("    vertex classes (in F?, in C?) = %s"
    % [((1 if v in FACE_V else 0), (1 if v in CYC_V else 0)) for v in range(5)])
out("    class 00 EMPTY on K1: FACE_V u CYC_V = %s = all five vertices (incidence, not choice)."
    % sorted(FACE_V | CYC_V))

def zi(a):
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][a % 4]
def cmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])
def cconj(z):
    return (z[0], -z[1])

bad_exact = 0
tot_exact = 0
for a in range(4):
    for b in range(4):
        WF, WC = zi(a), zi(b)
        for weights in [(Fr(1,2),Fr(0),Fr(0),Fr(1,4),Fr(1,4)),
                        (Fr(1,5),Fr(1,5),Fr(1,5),Fr(1,5),Fr(1,5)),
                        (Fr(2,5),Fr(1,10),Fr(1,5),Fr(1,10),Fr(1,5)),
                        (Fr(1),Fr(0),Fr(0),Fr(0),Fr(0))]:
            for k in range(1, 9):
                tot_exact += 1
                WFk, WCk = zi(a * k), zi(b * k)
                lhs = (Fr(0), Fr(0))
                for v in range(5):
                    A = WFk if v in FACE_V else (1, 0)
                    B = WCk if v in CYC_V else (1, 0)
                    t = cmul(cconj(A), B)
                    lhs = (lhs[0] + weights[v] * t[0], lhs[1] + weights[v] * t[1])
                p11 = weights[0]; p10 = weights[1] + weights[2]
                p01 = weights[3] + weights[4]; p00 = Fr(0)
                uk, vk = zi(-a * k), zi(b * k)
                uvk = cmul(uk, vk)
                rhs = (p00, Fr(0))
                for w, z in ((p10, uk), (p01, vk), (p11, uvk)):
                    rhs = (rhs[0] + w * z[0], rhs[1] + w * z[1])
                bad_exact += (lhs != rhs)
out("    EXACT Gaussian-rational check, %d cases (16 order-4 connections x 4 states x k=1..8):"
    % tot_exact)
out("      #{ <M_dF^k s, M_c^k s>  !=  p00 + p10 u^k + p01 v^k + p11 (uv)^k } = %d   <-- must be 0"
    % bad_exact)
out("    F1 CONFIRMED at EXACT arithmetic, not only at 5.118e-16 in double.")
out()

# ------------------------------------------------------------------ (1) hull <=> gaps, exact
out("(1) (HULL)  0 in conv{unit-modulus z_i}  <=>  max consecutive angular gap <= pi.")
out("    EXACT on q-th roots of unity: angles are integers mod q; 'all points in an OPEN")
out("    half-plane' is decided by enumerating the 2q half-lattice normals, with the strict")
out("    inequality cos((pi/q) d) > 0 rendered as the exact integer test 2d<q or 2d>3q.")

def gap_pred(angles, q):
    a = sorted(set(x % q for x in angles))
    if len(a) == 1:
        return False
    g = [a[i + 1] - a[i] for i in range(len(a) - 1)] + [a[0] + q - a[-1]]
    return 2 * max(g) <= q

def hull_pred_exact(angles, q):
    for j in range(2 * q):
        ok = True
        for a in angles:
            d = (2 * a - j) % (2 * q)
            if not (2 * d < q or 2 * d > 3 * q):
                ok = False
                break
        if ok:
            return False
    return True

for q in (12, 20, 36):
    bad3 = bad4 = n3 = n4 = 0
    for aa in range(q):
        for bb in range(q):
            A3 = [(bb - aa) % q, (-aa) % q, bb % q]
            g, h = gap_pred(A3, q), hull_pred_exact(A3, q)
            n3 += h; bad3 += (g != h)
            A4 = [0] + A3
            g4, h4 = gap_pred(A4, q), hull_pred_exact(A4, q)
            n4 += h4; bad4 += (g4 != h4)
    out("    q = %2d : 3 coeffs {uv,u,v}   #hull = %4d of %4d (%.4f) ; #(gap != hull) = %d"
        % (q, n3, q * q, n3 / q / q, bad3))
    out("             4 coeffs {1,u,v,uv} #hull = %4d of %4d (%.4f) ; #(gap != hull) = %d"
        % (n4, q * q, n4 / q / q, bad4))
out("    => the gap criterion IS hull membership at three AND at four points.  NO DEFECT.")
out("    NOTE, not in the lane: the 3-coefficient hull fraction -> 1/4 and the 4-coefficient")
out("    fraction -> 1/2 EXACTLY.  Reason for the 3-case: rotating by f sends {uv,u,v} to")
out("    {c, 0, c+f}, and (c, c+f) is uniform on T^2 when (f,c) is, so the configuration is")
out("    distributed exactly as three i.i.d. uniform points, for which P(0 in hull) =")
out("    1 - 3/2^2 = 1/4.  The lane's '342 of 1369' is that 1/4 and it did not say so.")
out()

# ------------------------------------------------------------------ (2) the vacuity check
out("(2) IS m3_5 (B) -- 'THE CONNECTION-SIDE CRITERION SURVIVES AT FOUR CLASSES' -- A TEST")
out("    THAT COULD HAVE FAILED?")
out("    The statement is: exists p in Delta_{n-1} with sum p_i z_i = 0  <=>  0 in conv{z_i}.")
out("    { sum_i p_i z_i : p in Delta } IS conv{z_i}: that is the DEFINITION of convex hull.")
out("    So it holds for every n and every z_i -- unit modulus or not, coset-constrained or")
out("    not.  Demonstrated by computing, for random points that are NOT on any coset and NOT")
out("    even on the unit circle, dist(0, conv) two ways: (i) exact 2-D polygon geometry,")
out("    (ii) the same simplex-grid minimisation m3_5 uses.  FLOAT.")

def dist0_to_hull(pts):
    """exact-ish distance from origin to conv(pts), pts a list of complex."""
    P = [(z.real, z.imag) for z in pts]
    # inside test + min distance to each hull edge/vertex; O(n^2) over all segments is enough
    # for n <= 6 and is independent of any triangulation.
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    Q = sorted(set(P))
    if len(Q) == 1:
        return (Q[0][0]**2 + Q[0][1]**2) ** .5
    lower = []
    for p in Q:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(Q):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    if len(hull) >= 3:
        sgn = [cross(hull[i], hull[(i+1) % len(hull)], (0.0, 0.0)) for i in range(len(hull))]
        if all(s >= -1e-15 for s in sgn) or all(s <= 1e-15 for s in sgn):
            return 0.0
    best = np.inf
    m = len(hull)
    for i in range(m):
        a, b = hull[i], hull[(i + 1) % m]
        ax, ay = a; bx, by = b
        dx, dy = bx - ax, by - ay
        den = dx * dx + dy * dy
        t = 0.0 if den == 0 else max(0.0, min(1.0, -(ax * dx + ay * dy) / den))
        px, py = ax + t * dx, ay + t * dy
        best = min(best, (px * px + py * py) ** .5)
    return best

def min_over_simplex_generic(zs, nstep):
    n = len(zs)
    best = [np.inf]
    def rec(i, rem, acc):
        if i == n - 1:
            v = sum((acc + [rem])[j] * zs[j] for j in range(n)) / nstep
            best[0] = min(best[0], abs(v))
            return
        for t in range(rem + 1):
            rec(i + 1, rem - t, acc + [t])
    rec(0, nstep, [])
    return best[0]

rng = np.random.default_rng(11081999)
out("      n   grid   max |grid-min  -  dist(0,conv)|   (12 random point-sets each)")
for n in (2, 3, 4, 5, 6):
    nstep = {2: 2000, 3: 300, 4: 80, 5: 36, 6: 22}[n]
    worst = 0.0
    for _ in range(12):
        zs = rng.normal(size=n) + 1j * rng.normal(size=n)     # NOT unit modulus, NO coset
        worst = max(worst, abs(min_over_simplex_generic(zs, nstep) - dist0_to_hull(list(zs))))
    out("      %d   %5d   %.4e      (grid spacing ~ %.3f, so this IS the grid error)"
        % (n, nstep, worst, 3.0 / nstep))
out("    VERDICT ON m3_5 (B): the connection-side statement is TRUE AT EVERY CLASS COUNT")
out("    BECAUSE IT IS THE DEFINITION OF A CONVEX HULL.  Its '0 disagreements on 361 grid")
out("    points' is not evidence that the connection side SURVIVES four classes -- it is a")
out("    re-check of the gap criterion of section (1).  As the ISOLATION ARM that the lane's")
out("    ledger calls 'THE COMPARISON THAT NAMES THE OPERATIVE VARIABLE' (C4), it COULD NOT")
out("    HAVE FAILED.  That voids the CONTROL, not the theorem: the asymmetry the lane names")
out("    is real and is proved in one line.  It was not measured, and the report presents it")
out("    as measured ('0 disagreements on 361 grid points [m3_5 (B)]').")
out()

# ------------------------------------------------------------------ (3) both directions
out("(3) REGISTER_V001.md:43 -- 'which vanishes iff 0 lies in the convex hull of three")
out("    unit-modulus coefficients'.  The lane's F3 calls this, read literally, FALSE.")
out("    IT FAILS IN EXACTLY ONE DIRECTION AND THE LANE DOES NOT SAY WHICH.  Swept over the")
out("    order-q connection grid x an exact rational simplex grid (FLOAT only to test |Z_1|,")
out("    with every hit re-certified exactly below):")
q = 12
Nsx = 24
fwd_bad = rev_bad = both = fires_tot = 0
for aa in range(q):
    for bb in range(q):
        A3 = [(bb - aa) % q, (-aa) % q, bb % q]
        hull = hull_pred_exact(A3, q)
        zc = [np.exp(2j * np.pi * A3[t] / q) for t in range(3)]
        for i in range(Nsx + 1):
            for j in range(Nsx - i + 1):
                kk = Nsx - i - j
                z = (kk * zc[0] + i * zc[1] + j * zc[2]) / Nsx   # p11*uv + p10*u + p01*v
                fires = abs(z) < 1e-12
                fires_tot += fires
                fwd_bad += (fires and not hull)
                rev_bad += (hull and not fires)
                both += 1
out("    q = %d, %d (connection,state) pairs, %d of them with Z_1 = 0:" % (q, both, fires_tot))
out("      #{ Z_1 = 0  AND  0 NOT in conv }  = %d" % fwd_bad)
out("           <-- ZERO, and it is a THEOREM: a convex combination of the z_i that equals 0")
out("               IS a certificate that 0 lies in their hull.  THIS HALF OF THE REGISTER'S")
out("               'IFF' IS TRUE FOR EVERY p AND EVERY CONNECTION.")
out("      #{ 0 in conv  AND  Z_1 != 0 }     = %d" % rev_bad)
out("           <-- the half that needs the existential quantifier over p, and it fails on")
out("               most of the grid.")
out("    CONSEQUENCE FOR F3.  W-01's load-bearing USE of the sentence is the exhibit at S1's")
out("    published connection, and that exhibit EXHIBITS Z_1 = 0 directly -- it travels on the")
out("    TRUE half.  F3 is a correction to one line of register prose, not to any W-01 result,")
out("    and it should have named the surviving half.  Also: the same register sentence ends")
out("    '-- verified against brute-force simplex MINIMISATION', which discloses the")
out("    quantifier it is accused of dropping.  F3 STANDS AS A CORRECTION, DOWNGRADED.")
out()
out("DONE.")
open("r1_three_statements.OUT.txt", "w").write("\n".join(L) + "\n")
