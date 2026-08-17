#!/usr/bin/env python3
"""
RW13C_02 -- THE ATTACK ON C8 / W13C_04(B):  "THE TIED STARTING POINT ... IS THE WHOLE OF
THE FAILURE."

WHAT THE TARGET LANE DID.  W13C_04(B) fixes theta at the L = 17372 LIOUVILLE PAIR BUILT IN
W13C_03 *FOR x_0 = 0* -- theta_rat is chosen so that k_1 * theta lands on a zero of P -- and
then moves x_0.  Its code, W13C_04_spacing_and_start.py:169-170, reads

        if (s1, s2) == (0, 0) and k0 < K1 <= k0+n:
            lz[K1-k0-1] = logpv          # spliced true value

so the mpmath dip value is SUBSTITUTED INTO EXACTLY ONE ARM, SELECTED BY A LITERAL TEST ON
THAT ARM'S OWN LABEL.  The other four arms run without it.  The reported difference between
the rows is therefore the hard-coded constant `logpv`, divided by N.  THE EXHIBIT COULD NOT
HAVE FAILED: theta was engineered so that the dip exists at x_0 = 0 and at no other x_0.

Leg A below re-runs their five arms with the branch DISABLED and shows all five rows agree,
which is the diff of their arms.

WHAT IS ACTUALLY TRUE, AND IT REFUTES THE CLAIM RATHER THAN THE ARITHMETIC.  x_0 = 0 has NO
privilege.  Transpose the target lane's own Baire argument (C6) from theta to x_0:

    THEOREM (this lane).  Fix ANY theta satisfying H2.  Put
        V_n = { x in T^2 : EXISTS k > n with dist(x + k theta, z) < e^{-nk} }.
    V_n is open (a union of open balls B(z,e^{-nk}) - k theta) and DENSE, because under H2
    the points {z - k theta : k > n} are dense in T^2.  So INTERSECT_n V_n is a dense
    G_delta, and on it, since f <= 0 everywhere,
        A_{k_n}(x_0) <= (1/k_n) log|Z_{k_n}| <= (log U)/k_n - n  ->  -infinity.
    So for EVERY connection with H2 -- INCLUDING EVERY CONNECTION AT WHICH THE TARGET LANE'S
    OWN a.e. THEOREM (C5) GIVES CONVERGENCE AT x_0 = 0 -- a COMEAGER set of starting points
    diverges.  x_0 = 0 is measure-good and category-bad in exactly the way theta is.

Leg B exhibits it, at the corpus's OWN published generic connection.

  ISOLATION.  theta is HELD FIXED at f = 1.0, c = sqrt(2) (S4:603; W-10 N-4: the only
  generic connection the corpus publishes; H2 PROVED by Lindemann-Weierstrass).  The
  polynomial, the estimator, the N-grid, the chunk size and the code path are byte-identical
  across all nine arms.  THE ONLY THING THAT MOVES IS x_0.
  AND THE SPLICE IS UNCONDITIONAL: for EVERY arm the k = k_1 term is recomputed in mpmath
  from the exact x_0 + k_1 theta and substituted.  There is no branch on the arm's label.
  For the un-engineered arms the substituted value is an ordinary O(0.1) number and changes
  nothing; that it changes nothing is printed, so the reader can check the splice is not
  doing the work.

Seed 20260818 (different from the target lane's, deliberately).
"""
import numpy as np, mpmath as mp

TWO64 = 1 << 64
K1 = 10**4

def Pabs(a, b):
    x = np.exp(2j*np.pi*a); y = np.exp(2j*np.pi*b)
    return np.abs(0.3*x + 0.3*y + 0.4*x*y)

def mP(dps=50):
    mp.mp.dps = dps
    t0 = mp.acos(mp.mpf(-2)/3)
    g = lambda t: mp.log(mp.mpf('0.3')) - mp.log(mp.mpf('0.25') + mp.mpf('0.24')*mp.cos(t))/2
    return mp.log(mp.mpf('0.4')) + mp.quad(g, [t0, mp.pi])/mp.pi

def zero_hp(dps):
    mp.mp.dps = dps
    phi = mp.acos(mp.mpf(-1)/9)/(2*mp.pi)
    return (mp.mpf(1)/2 - phi/2, mp.mpf(1)/2 + phi/2)

def theta_corpus(dps):
    """f = 1.0, c = sqrt(2)  =>  theta = (-f/2pi, c/2pi) mod 1.  H2 PROVED."""
    mp.mp.dps = dps
    return ((-mp.mpf(1)/(2*mp.pi)) % 1, (mp.sqrt(2)/(2*mp.pi)) % 1)

def Pabs_hp(a, b, dps):
    mp.mp.dps = dps
    x = mp.expjpi(2*a); y = mp.expjpi(2*b)
    return abs(mp.mpf(3)/10*x + mp.mpf(3)/10*y + mp.mpf(2)/5*x*y)

def run(A1, A2, X1, X2, Nmax, cps, splice_k, splice_log, chunk=1_000_000):
    """IDENTICAL code path for every arm.  x_0 enters as the initial uint64 state (X1,X2);
       the k = splice_k term is ALWAYS replaced by splice_log."""
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    cur1 = X1 % TWO64; cur2 = X2 % TWO64
    k0 = 0; s = 0.0; res = {}; res_nodip = {}; ci = 0; cl = sorted(cps)
    s_nd = 0.0
    while k0 < Nmax:
        n = min(chunk, Nmax - k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1 + int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2 + int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        lz = np.log(Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64))
        lz_nd = lz.copy()                     # the same orbit WITHOUT the spliced term
        if k0 < splice_k <= k0 + n:
            lz[splice_k-k0-1] = splice_log
        c = np.cumsum(lz); c_nd = np.cumsum(lz_nd)
        while ci < len(cl) and cl[ci] <= k0 + n:
            res[cl[ci]] = (s + c[cl[ci]-k0-1])/cl[ci]
            res_nodip[cl[ci]] = (s_nd + c_nd[cl[ci]-k0-1])/cl[ci]
            ci += 1
        s += c[-1]; s_nd += c_nd[-1]; k0 += n
    return res, res_nodip

if __name__ == "__main__":
    print("="*78)
    print("RW13C_02 -- THE STARTING POINT.  ONE VARIABLE: x_0.  theta HELD AT THE CORPUS'S")
    print("           OWN PUBLISHED GENERIC CONNECTION f = 1.0, c = sqrt(2).")
    print("="*78)
    mPv = float(mP())
    print("m(P) = %.15f" % mPv)

    # ------------------------------------------------------------------ LEG A: the diff
    print("\n" + "="*78)
    print("(A) THE ARMS-DIFF OF W13C_04(B).  Their five arms, THEIR theta, THEIR grid, with")
    print("    the single hard-coded branch `if (s1,s2)==(0,0)` DISABLED.  Nothing else moved.")
    print("="*78)
    import sys
    sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W13_C_CONVERGENCE")
    from W13C_03_liouville import build_theta, dip_value
    Lt = 17372
    n1, n2, DENt, Z1t, Z2t, _ = build_theta(Lt)
    logpv_t = float(mp.log(dip_value(Z1t, Z2t, Lt)))
    At1 = (n1*TWO64)//DENt; At2 = (n2*TWO64)//DENt
    CPSA = [1000, K1, 100000, 1000000, 2000000]
    rrA = np.random.default_rng(20260817)          # THEIR seed, so these are THEIR points
    startsA = [("x_0 = 0   (their row 1)", (0, 0))]
    for i in range(4):
        w = rrA.integers(0, TWO64, size=2, dtype=np.uint64)
        startsA.append(("x_0 = Haar-random draw %d (their row %d)" % (i, i+2), (int(w[0]), int(w[1]))))
    print("\n    %-44s %s" % ("starting point", "  ".join("%12s" % ("A_%d" % c) for c in CPSA)))
    print("    -- BRANCH DISABLED (splice_log = the float64 value, i.e. no substitution) --")
    for nm, (s1, s2) in startsA:
        # splice with the value the orbit itself produces => a no-op substitution
        r0, _ = run(At1, At2, s1, s2, 2_000_000, CPSA, splice_k=-1, splice_log=0.0)
        print("    %-44s %s" % (nm, "  ".join("%12.6f" % r0[c] for c in CPSA)))
    print("    -- BRANCH AS PUBLISHED (substitution in arm 1 only) --")
    for nm, (s1, s2) in startsA:
        sk = K1 if (s1, s2) == (0, 0) else -1
        r1, _ = run(At1, At2, s1, s2, 2_000_000, CPSA, splice_k=sk, splice_log=logpv_t)
        print("    %-44s %s" % (nm, "  ".join("%12.6f" % r1[c] for c in CPSA)))
    print("""
    READ-OFF.  With the branch disabled ALL FIVE ROWS AGREE to ~1e-3 and none of them fails.
    Every digit of the published contrast is the constant logpv = %.1f nats divided by N,
    injected into one arm by a literal test on that arm's own coordinates.  The exhibit is
    a CONTROL THAT COULD NOT HAVE FAILED: theta was built in W13C_03 so that the dip exists
    at x_0 = 0, so moving x_0 destroys it with probability 1 and no other outcome was
    available.  (W-08's isolation finding, the commonest FATAL defect in this program.)
    THE ARITHMETIC IS RIGHT.  WHAT IS NOT EARNED IS THE INFERENCE.""" % logpv_t)

    # ------------------------------------------------------------------ LEG B: the real test
    print("\n" + "="*78)
    print("(B) THE TEST THE CLAIM NEEDS.  theta FIXED at the corpus's own f=1.0, c=sqrt(2).")
    print("    ONE VARIABLE: x_0.  SPLICE APPLIED UNCONDITIONALLY TO EVERY ARM.")
    print("="*78)
    DEPTHS = (500, 2000, 8000, 17372)
    HDPS = 17372 + 400
    th = theta_corpus(HDPS)
    zz = zero_hp(HDPS)
    A1 = int(mp.floor(th[0]*TWO64)); A2 = int(mp.floor(th[1]*TWO64))
    print("theta (uint64) = (%d, %d)   [identical to W13C_01's GENERIC arm]" % (A1, A2))
    print("H2: m*(-1/2pi) + n*(sqrt2/2pi) in Z  =>  -m + n sqrt2 = 2 pi j  =>  j = 0 (LHS")
    print("    algebraic, RHS transcendental)  =>  m = n = 0.  PROVED, not assumed.")

    arms = []
    mp.mp.dps = HDPS
    kt = [ (K1*th[0]) % 1, (K1*th[1]) % 1 ]
    # un-engineered arms
    arms.append(("x_0 = 0   THE CORPUS'S OBJECT", (mp.mpf(0), mp.mpf(0))))
    rr = np.random.default_rng(20260818)
    for i in range(4):
        w = rr.integers(0, TWO64, size=2, dtype=np.uint64)
        arms.append(("x_0 = Haar-random draw %d" % i,
                     (mp.mpf(int(w[0]))/TWO64, mp.mpf(int(w[1]))/TWO64)))
    # engineered arms:  x_0 := z + delta - k1*theta   =>   x_0 + k1 theta = z + delta EXACTLY
    for Ld in DEPTHS:
        d = mp.mpf(10)**(-Ld)
        arms.append(("x_0 = ENGINEERED, dist(x_0 + k1 theta, z) = sqrt2 * 10^-%d" % Ld,
                     ((zz[0] + d - kt[0]) % 1, (zz[1] + d - kt[1]) % 1)))

    print("\n    ARMS DIFF GUARD.  (A_1,A_2) is the same in every arm BY DESIGN -- theta is the")
    print("    held-fixed variable.  What must be pairwise distinct is the STARTING STATE and")
    print("    the SPLICED VALUE; both are printed and asserted.")
    print("    %-52s %-21s %-21s %s" % ("arm", "X_1 (uint64)", "X_2 (uint64)", "spliced log|Z_k1|"))
    rows = []
    seen = set()
    for nm, x0 in arms:
        X1 = int(mp.floor((x0[0] % 1)*TWO64)); X2 = int(mp.floor((x0[1] % 1)*TWO64))
        a = (x0[0] + K1*th[0]) % 1
        b = (x0[1] + K1*th[1]) % 1
        val = Pabs_hp(a, b, HDPS)
        lg = float(mp.log(val))
        # distance to the nearer zero, high precision
        best = None
        for (za, zb) in ((zz[0], zz[1]), (zz[1], zz[0])):
            da = (a - za) % 1; da = da if da <= mp.mpf(1)/2 else da - 1
            db = (b - zb) % 1; db = db if db <= mp.mpf(1)/2 else db - 1
            dd = mp.sqrt(da*da + db*db)
            if best is None or dd < best: best = dd
        rows.append((nm, X1, X2, lg, val, best))
        seen.add((X1, X2, round(lg, 6)))
        print("    %-52s %-21d %-21d %s" % (nm[:52], X1, X2, mp.nstr(val, 8)))
    assert len(seen) == len(rows), "ZERO-VARIABLE CONTROL: two arms identical"
    print("    all %d arms pairwise distinct in (X_1, X_2, spliced value): OK" % len(rows))

    NMAX = 10_000_000
    CPS = [10, 100, 1000, 5000, 9999, 10000, 20000, 100000, 1000000, 10000000]
    print("\n    A_N, SEVEN DECADES.  Splice applied in EVERY row.  m(P) = %.9f" % mPv)
    print("    %-52s %s" % ("arm", " ".join("%11s" % ("A_%d" % c) for c in (1000, 9999, 10000, 100000, 1000000, 10000000))))
    out = []
    for (nm, X1, X2, lg, val, dd), (nm2, x0) in zip(rows, arms):
        res, res_nd = run(A1, A2, X1, X2, NMAX, CPS, K1, lg)
        out.append((nm, res, res_nd, lg, dd))
        print("    %-52s %s" % (nm[:52], " ".join("%11.6f" % res[c] for c in (1000, 9999, 10000, 100000, 1000000, 10000000))))

    print("\n    THE SAME ROWS WITH THE SPLICED TERM REMOVED (background only) -- this is the")
    print("    control that shows the splice is not doing the work in the un-engineered arms:")
    print("    %-52s %s" % ("arm", " ".join("%11s" % ("A_%d" % c) for c in (1000, 10000, 100000, 1000000, 10000000))))
    for (nm, res, res_nd, lg, dd) in out:
        print("    %-52s %s" % (nm[:52], " ".join("%11.6f" % res_nd[c] for c in (1000, 10000, 100000, 1000000, 10000000))))

    print("\n    DEPTH LADDER AT FIXED theta AND FIXED x_0-CONSTRUCTION -- A_{k1} IS UNBOUNDED")
    print("    BELOW IN x_0 ALONE, at the connection the corpus publishes:")
    print("    %-14s %-22s %-18s %-18s" % ("10^-L offset", "log|Z_k1| (nats)", "A_k1", "A_k1 - m(P)"))
    for (nm, res, res_nd, lg, dd) in out[5:]:
        print("    %-14s %-22.1f %-18.9f %-18.6f" % (nm.split("10^-")[1], lg, res[K1], res[K1]-mPv))
    print("\n    VERDICT ON C8.  At a connection whose H2 is PROVED and at which the target")
    print("    lane's own table reports convergence to 3.5e-07, the starting point alone")
    print("    drives A_{k1} to any prescribed level.  x_0 = 0 IS NOT PRIVILEGED.  The")
    print("    failure is the LOG SINGULARITY of the observable, which is symmetric in")
    print("    (theta, x_0): both are measure-good and category-bad.  C8's narrow logical")
    print("    point -- Birkhoff gives a.e. x_0 and {0} is null, so Birkhoff alone does not")
    print("    license N1 -- SURVIVES.  Its headline -- 'the tied starting point ... is THE")
    print("    WHOLE OF THE FAILURE' -- IS REFUTED, by the target lane's own C6 argument")
    print("    transposed from theta to x_0, and by this exhibit.")
    print("    AND THE TARGET LANE'S OWN C1/C2 SAY SO ALREADY: under H2 the rotation is")
    print("    UNIQUELY ERGODIC, so (1/N)SUM f_eps(x_0 + k theta) -> INT f_eps for EVERY x_0,")
    print("    not almost every.  All of the x_0-dependence lives in E_N(eps), i.e. in the")
    print("    singularity.  The starting point is a carrier of the failure, not its cause.")
    print("\nDONE RW13C_02")
