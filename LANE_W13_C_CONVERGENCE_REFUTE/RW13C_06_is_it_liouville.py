#!/usr/bin/env python3
"""
RW13C_06 -- IS THE "LIOUVILLE CONSTRUCTION" ACTUALLY LIOUVILLE?

The brief for this refuter names this check explicitly.  The target lane calls its pair
"THE EXPLICIT LIOUVILLE COUNTEREXAMPLE" (W13C_05 C6), names its script
W13C_03_liouville.py, and its conventions page says "THE LIOUVILLE CONSTRUCTION".  It is
built as

        theta_1 = n_1 / DEN  +  10^-M * sqrt(2),      DEN = k_1 * 10^L,  M = 10^6
        theta_2 = n_2 / DEN  +  10^-M * sqrt(3),

with n_1, n_2, DEN exact integers.

THE POINT.  theta_1 lies in Q(sqrt 2) and theta_2 in Q(sqrt 3): each is a QUADRATIC
IRRATIONAL.  By Lagrange every quadratic irrational has an eventually periodic continued
fraction, hence BOUNDED partial quotients, hence irrationality measure EXACTLY 2 -- it is
BADLY APPROXIMABLE, which is the exact opposite of Liouville.  A Liouville number needs
UNBOUNDED irrationality measure.

WHAT IS ACTUALLY TRUE OF THE CONSTRUCTION.  It is a SINGLE ENGINEERED NEAR-RESONANCE: one
rational approximation, of enormous height DEN ~ 10^(L+4), is unusually good (error
1.4*10^-M).  That produces exactly ONE deep dip, at k = k_1, and nothing after it.  It does
NOT give liminf A_N = -infinity for the theta it constructs.

AND THE COMPARISON WITH M1_06 RUNS THE OTHER WAY FROM WHAT C6 SAYS.  M1_06's
alpha = 1/3 + SUM_j d_j 10^{-a_j} with a_{j+1} - a_j -> infinity IS a Liouville number
(unbounded partial quotients) and its asserted construction DOES give liminf = -infinity for
ONE theta.  The target lane traded a genuinely-Liouville-but-unverified construction for a
verified-but-not-Liouville SINGLE DIP.  The real advance of the round is the Baire-category
PROOF (C6), which needs no construction at all.

THREE LEGS.
  (A) The exact integer minimal polynomial of theta_1 is exhibited and verified to vanish at
      theta_1 to hundreds of digits.  Quadratic => Lagrange => bounded partial quotients.
  (B) A SCALED-DOWN replica (M and L small enough to run the continued fraction to the end)
      is built by the target lane's OWN build_theta and its partial quotients are printed:
      ONE huge quotient at the engineered scale, BOUNDED ones afterwards forever.
      ONE VARIABLE: the tail exponent M.
  (C) The same continued fraction for a GENUINE Liouville number of M1_06's shape, which
      shows unbounded quotients at every scale.  Same code path, same qmax.
"""
import numpy as np, mpmath as mp
from fractions import Fraction

def cf(x, nterms, dps):
    mp.mp.dps = dps
    y = mp.mpf(x); out = []
    for _ in range(nterms):
        a = int(mp.floor(y)); out.append(a)
        fr = y - a
        if fr == 0: break
        y = 1/fr
    return out

if __name__ == "__main__":
    print("="*78)
    print("RW13C_06 -- IS THE 'LIOUVILLE CONSTRUCTION' ACTUALLY LIOUVILLE?")
    print("="*78)
    import sys
    sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W13_C_CONVERGENCE")
    from W13C_03_liouville import build_theta, M_TAIL, K1

    print("\n(A) THE MINIMAL POLYNOMIAL OF theta_1.  EXACT INTEGERS.")
    L = 500
    n1, n2, DEN, Z1, Z2, _ = build_theta(L)
    print("   the target lane's own build_theta(L=%d):  DEN = k1 * 10^L has %d digits" % (L, len(str(DEN))))
    print("   n1 has %d digits;  M = 10^%d" % (len(str(n1)), int(np.log10(M_TAIL))))
    print("   theta_1 = n1/DEN + 10^-M sqrt(2)   =>   (theta_1 - n1/DEN)^2 = 2 * 10^-2M")
    print("   =>  DEN^2 * 10^2M * theta_1^2  -  2 n1 DEN 10^2M * theta_1  +  (n1^2 10^2M - 2 DEN^2) = 0")
    print("   ALL THREE COEFFICIENTS ARE INTEGERS, so theta_1 is a QUADRATIC IRRATIONAL in")
    print("   Q(sqrt 2).  Likewise theta_2 in Q(sqrt 3).  BY LAGRANGE ITS CONTINUED FRACTION")
    print("   IS EVENTUALLY PERIODIC, ITS PARTIAL QUOTIENTS ARE BOUNDED, AND ITS")
    print("   IRRATIONALITY MEASURE IS EXACTLY 2.  A LIOUVILLE NUMBER HAS UNBOUNDED MEASURE.")
    print("   ==> THE PAIR IS BADLY APPROXIMABLE.  'LIOUVILLE' IS A MISNOMER.")
    # verify the quadratic vanishes, at 400 dps, on a scaled replica so mpmath can hold it
    Mv = 300
    mp.mp.dps = 900
    th1 = mp.mpf(n1)/DEN + mp.mpf(10)**(-Mv)*mp.sqrt(2)
    A_ = mp.mpf(DEN)**2 * mp.mpf(10)**(2*Mv)
    B_ = -2*mp.mpf(n1)*mp.mpf(DEN)*mp.mpf(10)**(2*Mv)
    C_ = mp.mpf(n1)**2*mp.mpf(10)**(2*Mv) - 2*mp.mpf(DEN)**2
    print("   check (with the tail exponent scaled to M=%d so mpmath can hold it):" % Mv)
    print("      |A theta_1^2 + B theta_1 + C| / A = %s   (relative machine zero)"
          % mp.nstr(abs(A_*th1*th1 + B_*th1 + C_)/A_, 6))

    print("\n(B) A SCALED-DOWN REPLICA, CONTINUED FRACTION RUN PAST THE ENGINEERED SCALE.")
    print("    ONE VARIABLE: the tail exponent M.  Everything else is build_theta's output.")
    Lr = 20
    n1r, n2r, DENr, _, _, _ = build_theta(Lr, k1=100)
    print("    replica: L = %d, k1 = 100, DEN = %d  (%d digits)" % (Lr, DENr, len(str(DENr))))
    for Mr in (60, 90):
        mp.mp.dps = 400
        th = mp.mpf(n1r)/DENr + mp.mpf(10)**(-Mr)*mp.sqrt(2)
        q = cf(th, 45, 400)
        print("    M = %-4d partial quotients: %s" % (Mr, q[:6]))
        big = max(q[1:]); ib = q[1:].index(big)+1
        print("             largest quotient a_%d has %d digits (the engineered approximation)"
              % (ib, len(str(abs(big)))))
        print("             quotients AFTER it: %s ..." % q[ib+1:ib+16])
        print("             max quotient after the engineered one: %d" % max(q[ib+1:]))
    print("    READ-OFF: ONE huge partial quotient, then ordinary bounded ones.  That is a")
    print("    single well-approximable rational of huge height, not Liouville behaviour.")

    print("\n(C) THE SAME CONTINUED FRACTION FOR M1_06's ACTUAL LIOUVILLE NUMBER.")
    print("    alpha = 1/3 + SUM_j 10^{-a_j}, a = 1, 8, 40, 130, 400  (M1_06's shape, scaled")
    print("    so mpmath can hold it).  Same code path, same nterms.")
    mp.mp.dps = 900
    al = mp.mpf(1)/3
    for a in (1, 8, 40, 130, 400):
        al += mp.mpf(10)**(-a)
    q2 = cf(al, 60, 900)
    sizes = [len(str(abs(t))) for t in q2]
    print("    partial-quotient DIGIT COUNTS: %s" % sizes[:40])
    print("    max digit count = %d, and a large quotient RECURS at every scale a_j -- the" % max(sizes))
    print("    signature of a Liouville number.  Contrast (B), where it happens exactly once.")
    print("""
    VERDICT.  MISNOMER, NOT A BROKEN RESULT.
      * C6's Baire-category THEOREM is untouched: it constructs nothing and needs no
        Diophantine class.  liminf A_N = -infinity on a comeager set STANDS.
      * C5's exhibit stands as what it is -- ONE engineered dip at k_1 with a free depth,
        which is exactly the claim the lane makes for it ("for any target T there is a theta
        with H2 and A_{k1} < T") and no more.
      * WHAT FALLS is the LABEL, and with it the sentence "M1_06 is confirmed and enlarged":
        M1_06's alpha IS Liouville and its (asserted) construction DOES give liminf =
        -infinity for a single theta; this lane's theta does neither.  In that one respect
        the new exhibit is WEAKER than the one it says it built past.  In every other
        respect (K1's registered pi, exact integers, seven decades, and a proof) it is
        stronger.  This is the program's standing defect -- SIX consecutive layers misnamed
        the operative variable (W-10 N-7, W-11) -- committed once more, in the file name.""")
    print("\nDONE RW13C_06")
