"""
X05 — TWO CHECKS AGAINST THE REGISTRAR'S OWN BRIEF, AND THE CARRIER CENSUS.

 (a) THE BRIEF'S BICONDITIONAL.  My brief for this round asserted "P has a torus zero IF AND
     ONLY IF W-01's criterion fires".  The L-refuter says that is THREE-CLASS-SCOPED and
     exhibits pi = (1/2, 3/10, 1/10, 1/10).  Checked here exactly, both directions, on the
     whole four-class simplex at three denominators, with the counterexample count printed.
 (b) THE CARRIER CENSUS.  The stratum of every SENSE-U row of S4:575 and of SENSE C, from
     S4's own class counts read at the bytes — with B0b coded AS S4 WRITES IT.
 (c) THE MEASURES.  The zero region as a fraction of the four-class and three-class simplices,
     by exact counting on a fine simplex lattice (not by Monte Carlo).
"""
from fractions import Fraction as F
import itertools, math
import numpy as np

def stratum(pi):
    p00, p10, p01, p11 = pi
    S1 = p00+p10; S2 = p01+p11; D1 = abs(p00-p10); D2 = abs(p01-p11)
    curve = ((p00 == p10 and p01 == p11) or (p00 == p01 and p10 == p11)
             or (p00 == p11 and p10 == p01))
    q = (S1-S2)*(D1-D2)
    if curve and q == 0: return "CURVE"
    return "EMPTY" if q > 0 else ("TWO" if q < 0 else "ONE")

def hull_fires(pi):
    """W-01's criterion as the register states it: 0 in the convex hull of the occupied
       unit-modulus coefficients <=> max_a p_a <= sum of the others."""
    return max(pi) <= sum(pi) - max(pi)

def mP(pi, n=1 << 20):
    p00, p10, p01, p11 = [float(t) for t in pi]
    t = (np.arange(n)+0.5)*(2*np.pi/n); ct = np.cos(t)
    a2 = p00*p00+p10*p10+2*p00*p10*ct; b2 = p01*p01+p11*p11+2*p01*p11*ct
    return float(np.mean(0.5*np.log(np.maximum(a2, b2))))

if __name__ == "__main__":
    print("="*96)
    print("X05 — THE BRIEF'S BICONDITIONAL, THE CARRIER CENSUS, AND THE EXACT MEASURES")
    print("="*96)

    print("\n(a) THE BRIEF'S BICONDITIONAL, TESTED EXACTLY ON THE FOUR-CLASS SIMPLEX.")
    print("    Named counterexample first (L-refuter's):")
    ce = (F(1,2), F(3,10), F(1,10), F(1,10))
    print("      pi = (1/2, 3/10, 1/10, 1/10):  max = %s <= %s = sum of the others -> HULL FIRES: %s"
          % (max(ce), sum(ce)-max(ce), hull_fires(ce)))
    print("      (S1-S2)(D1-D2) = (%s)(%s) = %s  ->  stratum %s"
          % ((ce[0]+ce[1])-(ce[2]+ce[3]), abs(ce[0]-ce[1])-abs(ce[2]-ce[3]),
             ((ce[0]+ce[1])-(ce[2]+ce[3]))*(abs(ce[0]-ce[1])-abs(ce[2]-ce[3])), stratum(ce)))
    print("      CONFIRMED: the hull fires and P HAS NO TORUS ZERO.  The brief is wrong as")
    print("      stated, and wrong in the same way W-09 convicted W-01's advertised virtue.")
    print()
    for D in (12, 18, 24):
        n4 = f4 = 0; n3 = f3 = 0
        for a in range(D+1):
            for b in range(D-a+1):
                for c in range(D-a-b+1):
                    d = D-a-b-c
                    pi = (F(a,D), F(b,D), F(c,D), F(d,D))
                    hz = stratum(pi) != "EMPTY"
                    hf = hull_fires(pi)
                    n4 += 1
                    if hz != hf: f4 += 1
                    if a == 0:
                        n3 += 1
                        if hz != hf: f3 += 1
        print("    denominator %2d :  four-class points %6d, DISAGREEMENTS %5d (%.1f%%) ;"
              "  three-class (p00=0) points %5d, DISAGREEMENTS %d"
              % (D, n4, f4, 100.0*f4/n4, n3, f3))
    print("    READ-OFF.  On the p00 = 0 face the biconditional is EXACT (0 disagreements at")
    print("    every denominator).  On the full four-class simplex it fails on a positive")
    print("    fraction.  The one direction that survives everywhere is TORUS ZERO => HULL.")

    print("\n(b) THE CARRIER CENSUS, from S4:575's own class counts, B0b AS S4 WRITES IT.")
    rows = [
        ("B0a ring torus, disjoint", {"00":2, "01":3, "10":4, "11":0}),
        ("B0b ring torus, meeting ", {"00":4, "01":1, "10":2, "11":2}),
        ("B3  horn torus          ", {"00":0, "01":2, "10":2, "11":1}),
        ("B1  K1                  ", {"00":0, "01":2, "10":2, "11":1}),
        ("B4  spindle             ", {"00":1, "01":1, "10":1, "11":3}),
        ("B2  K1 both filled      ", {"00":0, "01":2, "10":2, "11":1}),
        ("B1p K1-bridged          ", {"00":0, "01":3, "10":3, "11":0}),
        ("B1q K1-bridged+spectator", {"00":1, "01":3, "10":3, "11":0}),
        ("B1s K1 subdivided       ", {"00":0, "01":5, "10":5, "11":1}),
    ]
    cnt = {}
    for nm, cc in rows:
        V = sum(cc.values())
        pi = (F(cc["00"], V), F(cc["10"], V), F(cc["01"], V), F(cc["11"], V))
        st = stratum(pi); cnt[st] = cnt.get(st, 0)+1
        print("    %s  SENSE U pi = (%s,%s,%s,%s)  stratum %-6s  m(P) = %.12f"
              % (nm, pi[0], pi[1], pi[2], pi[3], st, mP(pi)))
    print("    SENSE U census over the nine runnable rows: %s" % cnt)
    print("    SENSE C rows (S4:566): 3 classes -> (0,3/10,3/10,2/5) stratum %s ;"
          % stratum((F(0), F(3,10), F(3,10), F(2,5))))
    print("                            4 classes -> (1/4,1/4,1/4,1/4) stratum %s ;"
          % stratum((F(1,4), F(1,4), F(1,4), F(1,4))))
    print("                            2 classes -> (0,1/2,1/2,0)     stratum %s"
          % stratum((F(0), F(1,2), F(1,2), F(0))))
    print("    S1's OWN PUBLISHED READY STATE (0,0,1/2,1/2) stratum %s"
          % stratum((F(0), F(0), F(1,2), F(1,2))))

    print("\n(c) THE EXACT MEASURES, by counting on the simplex lattice (not Monte Carlo).")
    for D in (60, 120, 240):
        tot4 = sing4 = 0; tot3 = sing3 = 0
        for a in range(D+1):
            for b in range(D-a+1):
                for c in range(D-a-b+1):
                    d = D-a-b-c
                    S1 = a+b; S2 = c+d; D1 = abs(a-b); D2 = abs(c-d)
                    sing = (S1-S2)*(D1-D2) <= 0
                    tot4 += 1; sing4 += sing
                    if a == 0:
                        tot3 += 1; sing3 += sing
        print("    D=%3d  four-class: singular fraction %.6f (target 1/4 = 0.250000)"
              "   three-class: %.6f" % (D, sing4/tot4, sing3/tot3))
    print("\nDONE X05")
