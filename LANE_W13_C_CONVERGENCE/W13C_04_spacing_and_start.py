#!/usr/bin/env python3
"""
W13C_04 -- TWO THINGS THE COUNTEREXAMPLE FORCES.

(A) THE ZERO SET IN CLOSED FORM, AND THE SPACING LEMMA.
    An exact closed form for Z(P) at K1's registered pi is derived and checked; it reduces
    the inhomogeneous Diophantine problem to a ONE-dimensional one.  Then:

    SPACING LEMMA.  If dist(k1 theta, z) <= d1 and dist(k2 theta, z) <= d2 with k1 < k2,
    then  || (k2 - k1) z ||  <=  k1 d2 + k2 d1     (distance to Z^2, EXACT identity below).
    So if z satisfies ||q z|| >= c q^-s, a SECOND dip at least as deep as the first forces
        k2  >~  (c / (2 d1))^(1/(s+1))  =  (c e^D / 2)^(1/(s+1)),  D = depth in nats.
    CONSEQUENCE, AND IT IS THE METHODOLOGICAL POINT OF THIS LANE:
    NO COMPUTATION OF ANY FEASIBLE SIZE CAN EXHIBIT TWO DEEP DIPS IN ONE ORBIT, SO NO
    NUMERICAL STUDY CAN EVER SEE liminf A_N = -infinity DIRECTLY.  It has to be proved.
    That is why every numerical check in this corpus "confirmed" N1.

(B) THE STARTING POINT.  Item (4) of the brief.  The orbit is {k theta}, i.e. the orbit of
    the SINGLE point 0 under rotation by theta.  Birkhoff's theorem fixes theta and gives
    a.e. STARTING POINT; 0 is one point, of measure zero, so Birkhoff does not apply.
    EXHIBITED: at the SAME Liouville theta of W13C_03, with the SAME N-grid, the SAME
    estimator and the SAME code path, moving ONLY the starting point makes the failure
    disappear.  One variable, and the verdict flips.
"""
import numpy as np, mpmath as mp
from W13C_01_central import TWO64, mahler_jensen, zeros_on_torus, Pabs
from W13C_03_liouville import build_theta, dip_value, K1

def cf_record_denoms(x, qmax, dps=60):
    """continued-fraction convergent denominators of x up to qmax, with ||q x||."""
    mp.mp.dps = dps
    a = []; y = mp.mpf(x); out = []
    p0, q0, p1, q1 = 0, 1, 1, 0
    for _ in range(200):
        ai = int(mp.floor(y)); a.append(ai)
        p0, p1 = p1, ai*p1 + p0
        q0, q1 = q1, ai*q1 + q0
        if q1 > qmax: break
        nrm = abs(q1*mp.mpf(x) - p1)
        out.append((q1, float(nrm), ai))
        fr = y - ai
        if fr == 0: break
        y = 1/fr
    return out, a

if __name__ == "__main__":
    print("="*78)
    print("W13C_04 (A) -- THE ZERO SET IN CLOSED FORM, AND THE SPACING LEMMA")
    print("="*78)
    mp.mp.dps = 50
    phi = mp.acos(mp.mpf(-1)/9)/(2*mp.pi)
    z1c = mp.mpf(1)/2 - phi/2
    z2c = mp.mpf(1)/2 + phi/2
    zs = zeros_on_torus(40)
    print("""
DERIVATION.  P = 0.3x + 0.3y + 0.4xy = 0 on T^2, x=e(a), y=e(b).  Divide by x:
   0.3 + 0.3 e(b-a) + 0.4 e(b) = 0.   Put phi = b-a, psi = b.
   |0.4| = 0.3|1+e(phi)| => 2+2cos(2 pi phi) = 16/9 => cos(2 pi phi) = -1/9.
   1+e(phi) = 2 cos(pi phi) e(phi/2), and cos(pi phi) = sqrt((1-1/9)/2) = 2/3, so
   0.4 e(psi) = -0.3 * 2 * (2/3) e(phi/2) = -0.4 e(phi/2)  =>  psi = 1/2 + phi/2.
   Hence   b = 1/2 + phi/2,  a = psi - phi = 1/2 - phi/2,  and  a + b = 1 EXACTLY.
   THE TWO ZEROS ARE  (1/2 -+ phi/2, 1/2 +- phi/2),  phi = arccos(-1/9)/(2 pi).""")
    print("   closed form      z = (%s, %s)" % (mp.nstr(z1c,25), mp.nstr(z2c,25)))
    print("   numerical solve  z = (%s, %s)" % (mp.nstr(zs[0][0],25), mp.nstr(zs[0][1],25)))
    print("   |difference| = %.3e   ;  z1 + z2 - 1 = %.3e" %
          (float(abs(z1c-zs[0][0])), float(abs(z1c+z2c-1))))
    print("   |P| at the closed-form point (float64) = %.3e" % Pabs(float(z1c), float(z2c)))
    print("""
   CONSEQUENCE 1.  z2 = 1 - z1, so ||q z|| (sup norm to Z^2) = ||q z1|| for every q:
   THE INHOMOGENEOUS 2-D PROBLEM COLLAPSES TO THE 1-D IRRATIONALITY OF
       z1 = 1/2 - arccos(-1/9)/(4 pi).
   CONSEQUENCE 2.  z1 is irrational: cos(2 pi phi) = -1/9 with phi rational would violate
   Niven's theorem (the only rational cosines of rational angles are 0, +-1/2, +-1).
   So ||q z|| > 0 for every q != 0 and the lemma below is never vacuous.""")

    print("\n   ||q z1|| along the continued-fraction convergents of z1 (record approximations):")
    conv, aa = cf_record_denoms(z1c, 10**9)
    print("   partial quotients: %s ..." % aa[:18])
    print("   %-14s %-16s %-14s" % ("q", "||q z1||", "q^2 ||q z1||"))
    for q, nrm, ai in conv[:20]:
        print("   %-14d %-16.6e %-14.4f" % (q, nrm, q*q*nrm))
    cmin = min(q*q*nrm for q, nrm, ai in conv)
    qmaxr = max(q for q, nrm, ai in conv)
    print("   measured  c := min over convergents q <= %d  of  q^2 ||q z1|| = %.5f" % (qmaxr, cmin))
    print("   (record approximations are always convergents, so this IS the true min over q <= %d)" % qmaxr)
    print("   NOT PROVED beyond that range: nothing is known about the irrationality measure of")
    print("   z1 = 1/2 - arccos(-1/9)/(4 pi).  The lemma below is therefore stated CONDITIONALLY,")
    print("   and an UNCONDITIONAL form is given after it.")
    print("   the partial quotients are small and unremarkable: z1 shows NO Liouville behaviour")
    print("   in this range.  s = 2 is the exponent to use in the lemma below on this evidence.")

    print("""
   THE IDENTITY BEHIND THE LEMMA, EXACT.  Write k_i theta = z + h_i + n_i, n_i in Z^2.
   Then  k1 h2 - k2 h1 = k1(k2 theta - z - n2) - k2(k1 theta - z - n1)
                       = (k2 - k1) z + (k2 n1 - k1 n2),
   so  ||(k2-k1) z||  =  ||k1 h2 - k2 h1||  <=  k1|h2| + k2|h1|.   No hypothesis at all.""")
    rng = np.random.default_rng(20260817)
    worst = 0.0
    for _ in range(20000):
        th = rng.random(2); k1 = int(rng.integers(1, 5000)); k2 = int(rng.integers(k1+1, 20000))
        z = np.array([float(z1c), float(z2c)])
        h1 = (k1*th - z + 0.5) % 1 - 0.5
        h2 = (k2*th - z + 0.5) % 1 - 0.5
        lhs = ((k2-k1)*z + 0.5) % 1 - 0.5
        rhs = (k1*h2 - k2*h1 + 0.5) % 1 - 0.5
        worst = max(worst, float(np.max(np.abs(((lhs-rhs)+0.5) % 1 - 0.5))))
    print("   identity verified on 20000 random (theta,k1,k2): worst residual mod 1 = %.3e" % worst)
    print("   (float64 roundoff on this check is ~k1*k2*2^-52 ~ 2e-8, so the residual IS the roundoff)")

    Lg = 1.026056   # global linear bound |P| >= Lg * dist(x,Z), measured in W13C_01
    print("\n   THE BOUND, EVALUATED.  |P| >= %.6f * dist (W13C_01), so a dip of depth D nats" % Lg)
    print("   means dist <= e^-D / %.6f.  With ||q z|| >= c q^-2, c = %.5f:" % (Lg, cmin))
    print("   %-12s %-24s %-28s" % ("depth D", "dist = e^-D/L", "k2 >= (c/(2 dist))^(1/3)"))
    for D in (1150.7, 4604.6, 18420.5, 39999.9):
        d1 = mp.e**(-mp.mpf(D))/Lg
        k2 = (mp.mpf(cmin)/(2*d1))**(mp.mpf(1)/3)
        print("   %-12.1f %-24s %-28s" % (D, "10^(%.1f)" % float(mp.log10(d1)), "10^(%.1f)" % float(mp.log10(k2))))
    print("""
   READ-OFF.  For the deepest arm of W13C_03 (D = 40000 nats) a second dip of the same
   depth cannot occur before k2 ~ 10^5790.  The number of atoms in the observable universe
   is about 10^80.  NO COMPUTATION WILL EVER SEE THE SECOND DIP.  The divergence
   liminf A_N = -infinity is real, and it is FORMALLY UNOBSERVABLE by simulation.
   UNCONDITIONAL FORM, needing nothing about z.  From ||(k2-k1)z|| <= (k1+k2) d1 <= 2 k2 d1:
   a second dip of depth >= D anywhere below k2 = 10^80 (more indices than there are atoms
   in the observable universe) would force  ||q z1|| <= 2 * 10^80 * e^-D / L  for some
   q <= 10^80, i.e. for D = 40000 an approximation of z1 to 10^(-17290) by a rational of
   denominator <= 10^80 -- an irrationality measure >= 216 in that range.  The continued
   fraction above (partial quotients 0,2,1,2,1,2,1,1,2,1,7,6,1,1,5,4,9,1,... up to q=1.8e7)
   gives no hint of any such thing.  EITHER WAY THE SECOND DIP IS NOT COMPUTABLE.
   AND THE UN-ENGINEERED SPACING IS WORSE STILL: dips of depth D occur along a generic
   orbit with density (2 pi/|det J|) e^{-2D} = 1.779 e^{-2D} per k, so the expected gap
   between depth-D dips is ~ e^{2D}/1.779 -- 10^%d for D = 40000.""" % int(2*39999.9/np.log(10)))

    # ---------------------------------------------------------------- (B) start point
    print("\n" + "="*78)
    print("W13C_04 (B) -- THE STARTING POINT.  ONE VARIABLE: x_0.  SAME theta, SAME EVERYTHING ELSE.")
    print("="*78)
    L = 17372
    n1, n2, DEN, Z1, Z2, _ = build_theta(L)
    pv = dip_value(Z1, Z2, L); logpv = float(mp.log(pv))
    A1 = (n1*TWO64)//DEN; A2 = (n2*TWO64)//DEN
    mPv = float(mahler_jensen(50))
    print("theta = the L=%d Liouville pair of W13C_03; A = (%d, %d)" % (L, A1, A2))
    print("m(P) = %.12f ;  the engineered dip is at k = %d with log|Z_k| = %.1f nats" % (mPv, K1, logpv))
    print("\nSTRUCTURAL FACT FIRST: |P| <= 1 on T^2 with equality ONLY at x = 0 (all three")
    print("occupied characters equal 1).  So the corpus's orbit STARTS AT THE GLOBAL MAXIMUM")
    print("of |P| -- f(0) = 0 exactly -- which is the most favourable start there is, and it")
    print("is still not enough.  |P(0,0)| = %.15f" % Pabs(0.0, 0.0))

    NMAX = 2_000_000
    CPS = [1000, K1, 100000, 1000000, 2000000]
    starts = [("x0 = 0  -- THE CORPUS'S OBJECT (Z_k = P(u^k,v^k), no offset)", (0, 0))]
    rr = np.random.default_rng(20260817)
    for i in range(4):
        w = rr.integers(0, TWO64, size=2, dtype=np.uint64)
        starts.append(("x0 = Haar-random draw %d (Birkhoff's a.e. starting point)" % i, (int(w[0]), int(w[1]))))

    print("\n%-58s %s" % ("starting point", "  ".join("A_%d" % c for c in CPS)))
    for name, (s1, s2) in starts:
        a1 = np.uint64(A1); a2 = np.uint64(A2)
        cur1, cur2 = s1, s2; k0 = 0; ssum = 0.0; res = {}; ci = 0; cps = sorted(CPS)
        while k0 < NMAX:
            n = min(1_000_000, NMAX-k0)
            i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
            i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
            x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
            cur1 = int(x1[-1]); cur2 = int(x2[-1])
            lz = np.log(Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64))
            if (s1, s2) == (0, 0) and k0 < K1 <= k0+n:
                lz[K1-k0-1] = logpv          # spliced true value; see W13C_03
            cs = np.cumsum(lz)
            while ci < len(cps) and cps[ci] <= k0+n:
                res[cps[ci]] = (ssum + cs[cps[ci]-k0-1])/cps[ci]; ci += 1
            ssum += cs[-1]; k0 += n
        print("%-58s %s" % (name[:58], "  ".join("%12.6f" % res[c] for c in CPS)))
    print("\nREAD-OFF.  Same theta.  Same estimator.  Same grid.  Same code path.  The ONLY")
    print("thing that moves is x_0, and the corpus's x_0 = 0 is the ONLY row that fails.")
    print("SO: BIRKHOFF'S THEOREM (a.e. x, theta fixed) IS NOT AVAILABLE TO N1, AND THE")
    print("DIFFERENCE IS NOT COSMETIC -- IT IS THE WHOLE OF THE FAILURE.")
    print("The a.e. result that IS available runs the other way: a.e. THETA, x_0 = 0 fixed.")
    print("It is proved in W13C_05, and it needs Borel-Cantelli in theta, not Birkhoff in x.")
    print("\nDONE W13C_04")
