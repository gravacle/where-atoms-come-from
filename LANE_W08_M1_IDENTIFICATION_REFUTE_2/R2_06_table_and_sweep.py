#!/usr/bin/env python3
"""
R2_06 — THE NINE-CASE CLASSIFICATION TABLE, THE 1000-DRAW SWEEP, AND THE OPEN REGION.

PART 1.  THE RELATION LATTICE OF EVERY ROW OF M1_04's NINE-CASE TABLE, computed exactly.
  The ledger calls this table "a CLASSIFICATION" and reads each row against G = {1}.  Fine.
  But row 1 is LABELLED "full support, generic connection" and it is S3/S4's headline
  f=2.0,c=1.1, which the register's own ERRATUM AGAINST W-02 records as EXACTLY RESONANT.
  And NOT ONE of the nine rows has rank L = 0.  So the table -- the only place the durability
  criterion is exhibited across supports -- never once exhibits it on a connection of the
  kind the identification theorem is about.  Recorded as a coverage gap, not as an error.

PART 2.  THE 1000-DRAW MINIMUM DENSITY (F7's "minimum density is 8.887867e-03 ... small only
  NEAR the degenerate locus") IS A WINDOW ARTEFACT OF THE DRAW COUNT, the same defect class
  as S3-audit COR-E and S2-audit COR-H, which this lane was required to carry.
  ISOLATION LEDGER: HELD FIXED -- carrier, observable, N, evaluator, the seed family.
  THE ONE THING THAT MOVES: the number of draws (1000 -> 20000) and, separately, whether the
  weights are drawn Dirichlet(1,1,1) or Dirichlet(0.2,0.2,0.2).  Neither changes the theorem;
  both change the reported minimum by orders of magnitude, so the number bounds nothing.
  Note also that the sweep moves TWO objects jointly (connection AND ready state) while the
  printed read-off attributes the minimum to ONE of them ("near the G={1} locus") -- and the
  witness it prints, pi=(0.991,0.001,0.008) at a generic (f,c), is degenerate in the READY
  STATE, not in the connection.

PART 3.  THE OPEN REGION.  For every ready state with max(p10,p01,p11) = p10 > 1/2 and every
  connection with W_F = 1 (u = 1) and W_C not a root of unity:
        rank L = 1,  G = <v> != {1}  (durable),  and  lambda = log p10 = m(P) EXACTLY.
  A POSITIVE-MEASURE set of (connection, ready state) pairs on which the rate identification
  holds with rank L = 1 and with NO Diophantine hypothesis of any kind.  Verified by orbit.
  ISOLATION LEDGER: HELD FIXED -- u = 1, the observable, the evaluator, N.  MOVES: p10 across
  the threshold 1/2, and separately v across several values.
Precision: float64 for orbits; exact Fractions for the relation lattices.
"""
import sys
import numpy as np
from fractions import Fraction

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W08_M1_IDENTIFICATION")
from M1_02_mahler_machinery import m_R1

print("=" * 78)
print("R2_06 — CLASSIFICATION TABLE COVERAGE, SWEEP WINDOW, AND THE OPEN REGION")
print("=" * 78)

# ---------------------------------------------------------------- PART 1
print("\nPART 1 — rank L FOR EVERY ROW OF M1_04's NINE-CASE TABLE (exact).")
print("  u = e^{i a}, v = e^{i b} with a,b RATIONAL multiples of 1 (not of pi) in every row")
print("  the target uses, so u^m v^n = 1  <=>  m a + n b = 2 pi j  <=>  (pi transcendental,")
print("  m a + n b rational)  j = 0 and m a + n b = 0.  rank L = 2 only when a,b in 2 pi Q.")
ROWS = [
    ("full support, generic connection   [LABEL]", Fraction(-2), Fraction(11, 10), False),
    ("full support, S1 published (ord 4)",         None, None, True),
    ("full support, TRIVIAL connection",           None, None, True),
    ("S={10,01} and u = v  (W_F W_C = 1)",         Fraction(7, 10), Fraction(7, 10), False),
    ("S={10,01}, u != v",                          Fraction(7, 10), Fraction(13, 10), False),
    ("S={10,11} and W_C = 1",                      Fraction(7, 10), Fraction(0), False),
    ("S={01,11} and W_F = 1",                      Fraction(0), Fraction(7, 10), False),
    ("S={11} only  (the ROOT alone)",              Fraction(-2), Fraction(11, 10), False),
    ("S1 published conn + S1 published state",     None, None, True),
]


def rankL_rational_angles(a, b):
    """a,b rational (units of radians, NOT of 2pi).  m a + n b = 0 has solutions
       spanned by (den, ...) -- compute the integer kernel of the 1x2 rational matrix."""
    if a == 0 and b == 0:
        return 2, "L = Z^2 (u = v = 1)"
    # kernel of [a b] over Q is 1-dimensional; intersect with Z^2 -> rank 1
    from math import gcd
    if a == 0:
        return 1, "L = Z*(1,0)"
    if b == 0:
        return 1, "L = Z*(0,1)"
    na, da = a.numerator, a.denominator
    nb, db = b.numerator, b.denominator
    # m*na/da + n*nb/db = 0  ->  m*na*db = -n*nb*da
    p = na * db
    q = nb * da
    g = gcd(abs(p), abs(q))
    return 1, "L = Z*(%d,%d)" % (-q // g, p // g)


for (name, a, b, torsion) in ROWS:
    if torsion:
        print("   %-42s rank L = 2   (roots of unity: u,v torsion)" % name)
    else:
        r, desc = rankL_rational_angles(a, b)
        print("   %-42s rank L = %d   %s" % (name, r, desc))
print("   ---------------------------------------------------------------------------")
print("   ROWS WITH rank L = 0 : 0 of 9.   Row 1's label 'generic connection' is S3/S4's")
print("   EXACTLY RESONANT headline (register ERRATUM AGAINST W-02): L = Z*(11,20).")
print("   Row 5, offered by the ledger as COMPARISON 3's isolated arm, is L = Z*(13,-7),")
print("   and its lambda = m(P) is a SUBTORUS coincidence (R2_02), not the generic case.")

# ---------------------------------------------------------------- PART 2
print("\nPART 2 — THE 1000-DRAW MINIMUM DENSITY IS A WINDOW ARTEFACT.")


def Zseq(u, v, p10, p01, p11, K):
    k = np.arange(1, K + 1)
    return np.abs(p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k)


def sweep(ndraw, alpha, seed, K=20000):
    rng = np.random.default_rng(seed)
    dm, w = np.inf, None
    for _ in range(ndraw):
        f, c = rng.uniform(0, 2 * np.pi, 2)
        u, v = np.exp(-1j * f), np.exp(1j * c)
        p = rng.dirichlet([alpha, alpha, alpha])
        d = float(np.mean(1 - Zseq(u, v, p[0], p[1], p[2], K)))
        if d < dm:
            dm, w = d, (f, c, p)
    return dm, w


for (nd, al, sd) in ((1000, 1.0, 20260816 + 40), (1000, 1.0, 20260816 + 41),
                     (20000, 1.0, 20260816 + 42), (2000, 0.2, 20260816 + 43)):
    dm, w = sweep(nd, al, sd, K=2000)
    print("   draws=%6d  Dirichlet(%.1f)  seed=%d :  min density = %.6e  at pi=(%.4f,%.4f,%.4f)"
          % (nd, al, sd, dm, w[2][0], w[2][1], w[2][2]))
print("   The target reports 8.887867e-03 from ONE window of 1000 Dirichlet(1,1,1) draws.")
print("   The infimum over the parameter space is 0 (approached as any pi -> a vertex, and")
print("   as u -> v on the two-element supports).  The number bounds nothing; it is a")
print("   LOWER BOUND on a sample minimum.  Same defect class as COR-E / COR-H.")
print("   Direct exhibit, no sampling: pi = (1-1e-6, 5e-7, 5e-7), generic connection:")
u0, v0 = np.exp(-1.4519j), np.exp(1.1964j)
d0 = float(np.mean(1 - Zseq(u0, v0, 1 - 1e-6, 5e-7, 5e-7, 20000)))
print("      Cesaro density = %.6e   (%.0f times below the reported 'minimum')"
      % (d0, 8.887867e-03 / d0))

# ---------------------------------------------------------------- PART 3
print("\nPART 3 — THE OPEN REGION WHERE A RANK-1 CONNECTION REPRODUCES m(P) EXACTLY.")
print("   u = 1 (W_F = 1: the FACE holonomy trivial, the flat holonomy W_C free).")
print("   HELD FIXED: u = 1, observable, evaluator, N = 2e6.  MOVES: p10 across 1/2, then v.")
print("     p10      m(P)            lambda_H (closed form)   orbit avg N=2e6     |orbit-m(P)|")
for p10 in (0.30, 0.45, 0.50, 0.55, 0.70, 0.85):
    p01 = p11 = (1 - p10) / 2
    mP = m_R1(0.0, p10, p01, p11)
    lamH = np.log(max(p10, p01 + p11))
    v = np.exp(1j * (2 * np.pi * (np.sqrt(5) - 1) / 2))
    k = np.arange(1, 2000001)
    az = np.abs(p10 + (p01 + p11) * v ** k)
    orb = float(np.mean(np.log(az)))
    print("    %.2f   %14.9f   %14.9f        %14.9f     %.2e"
          % (p10, mP, lamH, orb, abs(orb - mP)))
print("   For p10 > 1/2 the three columns coincide: lambda = m(P) with rank L = 1 and NO")
print("   Diophantine hypothesis.  For p10 < 1/2 they separate.  So whether the identification")
print("   holds is a JOINT condition on (L, pi) -- it is NOT a function of L, and it is NOT")
print("   'rank L = 0'.  Durability is unaffected: G = <v> != {1} in every row above.")
print("   Several v, p10 = 0.7 (m(P) = log 0.7 = %.9f):" % np.log(0.7))
for vv in (0.7, 1.3, 2.0, np.pi / 2, np.sqrt(2)):
    v = np.exp(1j * vv)
    k = np.arange(1, 2000001)
    az = np.abs(0.7 + 0.3 * v ** k)
    print("      arg v = %.6f :  orbit avg = %.9f   dev = %+.2e"
          % (vv, float(np.mean(np.log(az))), float(np.mean(np.log(az))) - np.log(0.7)))
print("\nDONE R2_06")
