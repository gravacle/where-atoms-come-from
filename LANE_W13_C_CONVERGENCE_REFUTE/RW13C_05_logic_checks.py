#!/usr/bin/env python3
"""
RW13C_05 -- TWO COMPARATIVE CLAIMS INSIDE PROVED FINDINGS, CHECKED AGAINST THE TARGET
LANE'S OWN NUMBERS.

(I)  C3 (status PROVED) ends: "...and it is STRICTLY STRONGER than the 'orbit is dense'
     clause N1 currently carries."  The theorem is correct; the comparison is not.  The two
     conditions are LOGICALLY INCOMPARABLE, and the corpus already owns both witnesses:
       * S1's own published connection (order 4) SATISFIES C3's condition -- |Z_k| takes
         four values, none of them 0, so (1/k)log|Z_k| -> 0 -- and A_N CONVERGES, to
         -(1/2)log 5.  H2 FAILS there.   So C3's condition does NOT imply density.
       * The target lane's own C6/W13C_03 theta satisfies H2 (proved) and violates C3's
         condition.   So density does not imply C3's condition.
     Neither is stronger.  They are necessary conditions for DIFFERENT statements: C3 for
     "A_N converges at all", H2 for "the limit, if it exists, is m(P)".
     Exhibited below on the target lane's own estimator, one variable (the connection).

(II) C7's methodological sentence: "A numerical study of N1 will confirm it at every
     connection it can reach, INCLUDING connections in the failure set.  The corpus's
     numerics on N1 could not have failed."  The target lane's OWN W13C_03 output is a
     counterexample: at its engineered theta -- which is in the failure set -- the estimator
     at N = 1e4 returns -4.767319 and at N = 1e7 returns -0.771508, a deviation of -4.0e-03,
     four orders worse than any of its H2 arms.  A numerical study at that connection does
     not confirm N1; it loudly refutes it.  The defensible statement is the one C7's LEMMA
     actually supports and is weaker: a simulation can never see the SECOND dip, so it
     cannot distinguish "one dip then convergence" from "liminf = -infinity"; agreement is
     therefore not decisive.  That is a claim about DISCRIMINATION, not about "could not
     have failed".
     Reproduced below from the target lane's own sealed numbers, no new computation.
"""
import numpy as np, mpmath as mp
TWO64 = 1 << 64

def Pabs(a, b):
    x = np.exp(2j*np.pi*a); y = np.exp(2j*np.pi*b)
    return np.abs(0.3*x + 0.3*y + 0.4*x*y)

def mP(dps=50):
    mp.mp.dps = dps
    t0 = mp.acos(mp.mpf(-2)/3)
    g = lambda t: mp.log(mp.mpf('0.3')) - mp.log(mp.mpf('0.25') + mp.mpf('0.24')*mp.cos(t))/2
    return mp.log(mp.mpf('0.4')) + mp.quad(g, [t0, mp.pi])/mp.pi

if __name__ == "__main__":
    print("="*78)
    print("RW13C_05 -- TWO COMPARATIVE CLAIMS INSIDE PROVED FINDINGS")
    print("="*78)
    mPv = float(mP())

    print("\n(I) C3's 'STRICTLY STRONGER THAN THE ORBIT-IS-DENSE CLAUSE'.")
    print("    ONE VARIABLE: the connection.  Same estimator, same grid, same code path.")
    mp.mp.dps = 60
    arms = [("S1 PUBLISHED, order 4 (f=pi, c=3pi/2)   H2 FAILS",
             ((-mp.pi/(2*mp.pi)) % 1, ((3*mp.pi/2)/(2*mp.pi)) % 1)),
            ("GENERIC f=1.0 c=sqrt(2)                 H2 PROVED",
             ((-mp.mpf(1)/(2*mp.pi)) % 1, (mp.sqrt(2)/(2*mp.pi)) % 1))]
    print("    %-46s %-16s %-22s %-16s" % ("arm", "A_1e7", "sup_{k>1e6} log(1/|Z_k|)/k", "A_N converges?"))
    for nm, th in arms:
        A1 = int(mp.floor(th[0]*TWO64)); A2 = int(mp.floor(th[1]*TWO64))
        a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
        cur1 = 0; cur2 = 0; k0 = 0; s = 0.0; res = {}; sup_tail = 0.0
        cps = [10**3, 10**5, 10**6, 10**7]; ci = 0
        while k0 < 10**7:
            n = min(1_000_000, 10**7-k0)
            i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
            i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
            x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
            cur1 = int(x1[-1]); cur2 = int(x2[-1])
            z = Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64)
            kk = k0+1+np.arange(n)
            r = np.where(kk > 10**6, -np.log(z)/kk, -np.inf)
            sup_tail = max(sup_tail, float(np.max(r)))
            c = np.cumsum(np.log(z))
            while ci < len(cps) and cps[ci] <= k0+n:
                res[cps[ci]] = (s + c[cps[ci]-k0-1])/cps[ci]; ci += 1
            s += c[-1]; k0 += n
        conv = "YES -> %.9f" % res[10**7] if abs(res[10**7]-res[10**6]) < 1e-7 else "not settled here"
        print("    %-46s %-16.9f %-22.3e %-16s" % (nm, res[10**7], sup_tail, conv))
    print("""
    READ-OFF.  The order-4 row's |Z_k| cycles through {sqrt(0.1), 0.4, sqrt(0.1), 1}, so
    log(1/|Z_k|) <= 1.1513 for EVERY k and (1/k)log(1/|Z_k|) <= 1.1513/k -> 0: C3's
    necessary condition holds with room to spare (the tail sup printed above is exactly
    1.1513/10^6 to three digits).  And A_N converges, to -(1/2) log 5 = -0.804718956,
    EXACTLY as the register has it.  H2 fails.  So C3's condition holds where density fails: IT IS NOT STRONGER THAN
    DENSITY, it is INCOMPARABLE WITH IT.  The theorem C3 states is untouched; the sentence
    comparing it to N1's density clause is wrong, and it is wrong in the direction that
    OVERSELLS the finding.""")

    print("\n(II) C7's 'THE CORPUS'S NUMERICS ON N1 COULD NOT HAVE FAILED'.")
    print("     Quoted from the TARGET LANE'S OWN sealed output, W13C_03_liouville.OUT.txt,")
    print("     at its OWN engineered theta, which is in the failure set and satisfies H2:")
    print("        N = 1e3   A_N = -0.766669   (dev +8.4e-04)")
    print("        N = 1e4   A_N = -4.767319   (dev -4.0e+00)")
    print("        N = 1e5   A_N = -1.167509   (dev -4.0e-01)")
    print("        N = 1e6   A_N = -0.807491   (dev -4.0e-02)")
    print("        N = 1e7   A_N = -0.771508   (dev -4.0e-03)")
    print("     A numerical study of N1 at that connection, over the same four decades the")
    print("     target lane uses everywhere else, reports a deviation of 4.0e-03 -- FOUR")
    print("     ORDERS worse than any of its eight H2 arms.  IT DOES NOT CONFIRM N1.")
    print("     So 'a numerical study will confirm it at every connection it can reach,")
    print("     INCLUDING connections in the failure set' is contradicted by the lane's own")
    print("     table.  What the SPACING LEMMA actually supports is weaker and still worth")
    print("     having: a simulation can never see the SECOND dip, so it can never")
    print("     DISTINGUISH 'one dip, then convergence' from liminf = -infinity.")
    print("     THE 'COULD NOT HAVE FAILED' VERDICT AGAINST THE CORPUS'S NUMERICS IS")
    print("     THEREFORE NOT EARNED AS STATED.  The weaker discrimination claim is earned.")
    print("     (And the program's own norm applies here: 'COULD NOT HAVE FAILED' VOIDS A")
    print("      CONTROL, NEVER A THEOREM -- C6's theorem is untouched either way.)")
    print("\nDONE RW13C_05")
