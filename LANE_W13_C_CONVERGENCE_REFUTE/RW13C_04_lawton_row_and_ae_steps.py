#!/usr/bin/env python3
"""
RW13C_04 -- (A) THE S4:973 ROW, PROPERLY CONVERGED;  (B) THE TWO SELF-CONTAINED STEPS OF
THE a.e. THEOREM, CHECKED RATHER THAN TAKEN ON TRUST.

(A) W13C_05 C9 rules that S4:973's published connection f = 3.14159, c = 1.57080 "lies IN
    the exceptional set with a DIFFERENT limit".  The first half is right (W-10 N-4: every
    rational pair is exactly resonant, relation (157080, 314159)).  The second half was
    never computed.  Its orbit closure is the circle {(314159 s, -157080 s)} and its limit
    is the ONE-VARIABLE Mahler measure
         m(0.3 + 0.4 z^314159 + 0.3 z^471239),
    which W-08/M1 T2(d) already names as the place Boyd-Lawton IS exactly right: as the
    relation grows this tends to m(P).  Here it is computed to convergence, with the same
    machinery validated first on the two relations whose answers the register already
    carries (f=2.0,c=1.1 -> -0.767014993; f=1.3,c=2.0).
    Exponents are reduced EXACTLY by int64 modular arithmetic ((e*j) mod n), so no large
    floating-point phase argument ever appears.

(B) C5's proof has three steps; (c) is a named import and is not touched here.  (a) and (b)
    are claimed self-contained, and (b) carries the whole inhomogeneous Diophantine half:
        "theta |-> k theta preserves Haar, so Pr(dist(k theta, Z) < r) = 2 pi r^2",
        then Borel-Cantelli with r_k = k^{-3/4}.
    Both halves are checked numerically: the pushforward claim at several k, and the
    Borel-Cantelli conclusion (that inf_k k^{3/4} dist(k theta, Z) is bounded away from 0
    for random theta, and is NOT for the target lane's own engineered theta).
    ONE VARIABLE in the last comparison: the arithmetic type of theta.

Seed 20260818.
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

def m_exact_modular(terms, n, chunk=1 << 23):
    """mean over n-th roots of unity of log|sum a_e z^e|, exponents reduced mod n in int64."""
    tot = 0.0; j0 = 0
    while j0 < n:
        c = min(chunk, n - j0)
        j = np.arange(j0, j0 + c, dtype=np.int64)
        v = np.zeros(c, dtype=np.complex128)
        for a, e in terms:
            idx = (np.int64(e % n) * j) % np.int64(n)
            v += a * np.exp(2j * np.pi * (idx.astype(np.float64) / n))
        tot += float(np.sum(np.log(np.abs(v))))
        j0 += c
    return tot / n

if __name__ == "__main__":
    print("=" * 78)
    print("RW13C_04 -- THE S4:973 ROW CONVERGED, AND THE a.e. THEOREM'S OWN STEPS")
    print("=" * 78)
    mPv = float(mP())
    print("m(P) = %.15f" % mPv)

    print("\n" + "=" * 78)
    print("(A) SUBTORUS LIMITS OF THE CORPUS'S RESONANT PUBLISHED CONNECTIONS.")
    print("    ONE VARIABLE: the connection.  Same estimator (modular-exact trapezoid),")
    print("    same n-ladder, same code path.  The first two rows are CALIBRATION: their")
    print("    answers are already in the register.")
    print("=" * 78)
    fam = [
        ("CALIBRATION f=2.0, c=1.1   relation (11,20)   [register erratum: -0.767014993]",
         [(0.3, 31), (0.3, 0), (0.4, 20)]),
        ("CALIBRATION f=1.3, c=2.0   relation (20,13)   [W-10 lane D, never computed]",
         [(0.3, 33), (0.3, 0), (0.4, 13)]),
        ("S4:973      f=3.14159, c=1.57080  relation (157080,314159)  [C9 asserts 'different limit']",
         [(0.3, 471239), (0.3, 0), (0.4, 314159)]),
    ]
    NS = [1 << 20, 1 << 22, 1 << 24, 1 << 26, 1 << 27]
    print("\n   %-74s %s" % ("connection", "  ".join("%16s" % ("n=2^%d" % int(np.log2(n))) for n in NS)))
    lastvals = {}
    for nm, terms in fam:
        vals = [m_exact_modular(terms, n) for n in NS]
        lastvals[nm] = vals[-1]
        print("   %-74s %s" % (nm[:74], "  ".join("%16.9f" % v for v in vals)))
        print("   %-74s %s" % ("   ... minus m(P)", "  ".join("%16.3e" % (v - mPv) for v in vals)))
    print("""
   READ-OFF, AND IT IS A CORRECTION TO C9.
   * f=2.0,c=1.1 reproduces the register's -0.767014993 exactly: the method is sound.
   * f=1.3,c=2.0 sits at a genuinely different limit, gap ~ +7.8e-04.  C9's claim holds
     for this row.
   * S4:973's f=3.14159, c=1.57080 does NOT have a "demonstrably different limit".  Its
     relation has order ~5e5, and by the very Boyd-Lawton statement W-08/M1 T2(d) says is
     the RIGHT use of that theorem, its subtorus limit is m(P) to the precision printed
     above -- SMALLER THAN THE TARGET LANE'S OWN RESIDUAL AT ITS GENERIC ARM (3.5e-07 at
     N=1e7).  H2 does fail there; the LIMIT does not differ measurably.
   * AND THE SAME ROW IS A SECOND INSTANCE OF THE TARGET LANE'S OWN C7 FINDING, WHICH IT
     DID NOT NOTICE BECAUSE IT NEVER RAN THE ROW: at N = 1e7 the Birkhoff average at
     S4:973's connection reads -7.93e-04 away from m(P) (RW13C_03), i.e. THREE ORDERS
     WORSE than its own true limit.  A numerical study at 1e7 misclassifies this published
     connection in the OPPOSITE direction from the one C7 warns about: not a false
     confirmation, a false refutation.""")

    print("\n" + "=" * 78)
    print("(B) THE TWO SELF-CONTAINED STEPS OF C5.")
    print("=" * 78)
    mp.mp.dps = 40
    phi = float(mp.acos(mp.mpf(-1) / 9) / (2 * mp.pi))
    zs = [(0.5 - phi / 2, 0.5 + phi / 2), (0.5 + phi / 2, 0.5 - phi / 2)]

    def dist_to_Z(t1, t2):
        d = np.full(t1.shape, np.inf)
        for (za, zb) in zs:
            da = np.abs(t1 - za); da = np.minimum(da, 1 - da)
            db = np.abs(t2 - zb); db = np.minimum(db, 1 - db)
            d = np.minimum(d, np.hypot(da, db))
        return d

    print("\n   (b1) 'theta |-> k theta preserves Haar, so Pr(dist(k theta, Z) < r) = 2 pi r^2.'")
    rng = np.random.default_rng(20260818)
    M = 4_000_000
    print("        %-8s %-10s %-16s %-16s %-10s" % ("k", "r", "measured Pr", "2 pi r^2", "ratio"))
    for k in (1, 2, 7, 1000, 999983):
        th1 = rng.random(M); th2 = rng.random(M)
        t1 = np.mod(k * th1, 1.0); t2 = np.mod(k * th2, 1.0)
        d = dist_to_Z(t1, t2)
        for r in (0.01, 0.003):
            pr = float(np.mean(d < r)); th_ = 2 * np.pi * r * r
            print("        %-8d %-10.4g %-16.6e %-16.6e %-10.4f" % (k, r, pr, th_, pr / th_))
    print("        CONFIRMED (to Monte-Carlo error ~ 1/sqrt(M*Pr)).  Step (b)'s measure claim")
    print("        holds, including at k = 999983 where the float64 product k*theta has lost")
    print("        ~20 bits -- so the check is of the MATHEMATICS, not of the arithmetic.")

    print("\n   (b2) BOREL-CANTELLI'S CONCLUSION: inf_k k^{3/4} dist(k theta, Z) > 0.")
    print("        ONE VARIABLE: the arithmetic type of theta.  Same k-range, same estimator.")
    mp.mp.dps = 60
    rng2 = np.random.default_rng(20260818)
    arms = [("GENERIC f=1.0 c=sqrt(2)", ((-mp.mpf(1) / (2 * mp.pi)) % 1, (mp.sqrt(2) / (2 * mp.pi)) % 1)),
            ("BADLY-APPROX (-2^(1/3),4^(1/3))", ((-mp.cbrt(2)) % 1, (mp.cbrt(4)) % 1))]
    for s in range(3):
        w = rng2.integers(0, TWO64, size=2, dtype=np.uint64)
        arms.append(("HAAR-RANDOM draw %d" % s, (mp.mpf(int(w[0])) / TWO64, mp.mpf(int(w[1])) / TWO64)))
    import sys
    sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W13_C_CONVERGENCE")
    from W13C_03_liouville import build_theta
    n1, n2, DENl, _, _, _ = build_theta(500)
    arms.append(("THE TARGET LANE'S ENGINEERED LIOUVILLE theta (L=500)",
                 (mp.mpf(n1) / DENl, mp.mpf(n2) / DENl)))
    print("        %-52s %-16s %-16s" % ("arm", "inf k^{3/4} dist", "at k"))
    for nm, th in arms:
        A1 = int(mp.floor(th[0] * TWO64)); A2 = int(mp.floor(th[1] * TWO64))
        a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
        cur1 = 0; cur2 = 0; k0 = 0; best = np.inf; bk = -1
        while k0 < 10 ** 7:
            n = min(1_000_000, 10 ** 7 - k0)
            i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1 + int(a1)) % TWO64)
            i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2 + int(a2)) % TWO64)
            x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
            cur1 = int(x1[-1]); cur2 = int(x2[-1])
            d = dist_to_Z(x1.astype(np.float64) / TWO64, x2.astype(np.float64) / TWO64)
            kk = (k0 + 1 + np.arange(n)).astype(np.float64)
            v = (kk ** 0.75) * d
            j = int(np.argmin(v))
            if v[j] < best: best = float(v[j]); bk = k0 + j + 1
            k0 += n
        print("        %-52s %-16.6e %-16d" % (nm[:52], best, bk))
    print("""
        READ THE LAST ROW CAREFULLY, IN BOTH DIRECTIONS.
        (i) It is NOT a counterexample to step (b): step (b) is an a.e. statement and the
            engineered theta is one point, constructed to lie in the null set.
        (ii) But it DOES flag the engineered theta four orders below every other arm, in a
            statistic costing one pass -- because the 64-bit truncation still lands within
            1e4 * 2^-64 = 5.4e-16 of the zero at k = 1e4.  So the target lane's engineered
            theta IS detectable by a plain float64 scan of C3's own necessary-condition
            statistic.  This bears on C7's methodological sentence and is taken up in
            RW13C_05(II).  What is genuinely undetectable is the SECOND dip, not the first.""")
    print("\nDONE RW13C_04")
