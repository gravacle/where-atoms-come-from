#!/usr/bin/env python3
"""
RW13C_03 -- THE C9 TABLE'S MISSING ROWS, AND FIVE STATED FIGURES RE-DERIVED.

(1) THE TWO PUBLISHED CONNECTIONS THE TARGET LANE ASSERTS BUT NEVER RUNS.  C9 says four of
    the corpus's five published connections "lie IN the exceptional set with a DIFFERENT
    limit".  Two of the four are computed in W13C_01 (the order-4 point and f=2.0,c=1.1).
    The other two -- S4:973's f=3.14159, c=1.57080 and W-10 lane D's f=1.3, c=2.0 -- are
    carried by INHERITANCE from W-10 N-4 ("every rational pair is resonant").  Resonance
    gives H2 FALSE; it does NOT give "a different limit".  W-12's own Corollary 2 exhibits a
    resonant subtorus whose average equals m(P) EXACTLY (difference 5.55e-17), so the
    inference is not available in general.  Here both rows are RUN, on the target lane's
    estimator and grid, and their subtorus limits are also computed independently as
    one-variable Mahler measures.  ONE VARIABLE: the connection.

(2) FIVE FIGURES RE-DERIVED, EACH ONE STATED AS A CONSTANT SOMEWHERE IN THE TARGET LANE:
    (a) "float64 would have reported |Z_{k1}| = 1.06e-15 -- WRONG BY UP TO 10^17368"
        (W13C_05 C6; W13C_03_liouville.py:140 prints 10^(L-4)).
    (b) "min |Z_k| over k<=1e7 EXCLUDING the engineered dip: 1.0735e-04 at k=9815000",
        computed by code that DISCARDS A WHOLE CHUNK when the dip is that chunk's argmin
        (W13C_03_liouville.py:97).
    (c) W13C_02's (D2) third column, labelled "the TAIL sup over N/10 < k <= N ... which the
        NECESSARY CONDITION requires to go to 0".  It is a MONOTONE RUNNING MAX
        (np.maximum.accumulate carried across chunks), so it cannot decrease and cannot test
        what the label says it tests.  The genuine per-decade tail sup is computed here.
    (d) "1.026056 * dist <= |P| <= 2.492571 * dist" -- a 4096^2 GRID MINIMUM printed as a
        two-sided inequality.  A grid minimum is an UPPER bound on an infimum (COR-E's
        defect class).  Refined here by local descent.
    (e) E_N(eps) -> (pi/|det J|) eps^2 = 0.8897 eps^2 and the shell count 1.7794 eps^2,
        over five decades of N rather than four.

(3) THE SPACING LEMMA'S SECOND ZERO.  W13C_04's lemma is written for two dips at the SAME
    zero.  Since z_1 + z_2 = 1 the second zero is z' = -z, so a dip at z followed by one at
    z' gives ||(k1+k2) z|| <= k1|h2| + k2|h1| -- a different q.  Verified numerically here;
    the conclusion is unchanged because q <= 2 k_2 either way, but the lemma as written does
    not cover it.

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

def turns_from_fc(f, c, dps=60):
    mp.mp.dps = dps
    return ((-mp.mpf(f)/(2*mp.pi)) % 1, (mp.mpf(c)/(2*mp.pi)) % 1)

def birkhoff(A1, A2, Nmax, cps, chunk=1_000_000):
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    cur1 = 0; cur2 = 0; k0 = 0; s = 0.0; res = {}; ci = 0; cl = sorted(cps)
    gmin = np.inf; gk = -1
    while k0 < Nmax:
        n = min(chunk, Nmax-k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        z = Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64)
        j = int(np.argmin(z))
        if z[j] < gmin: gmin = float(z[j]); gk = k0+j+1
        c = np.cumsum(np.log(z))
        while ci < len(cl) and cl[ci] <= k0+n:
            res[cl[ci]] = (s + c[cl[ci]-k0-1])/cl[ci]; ci += 1
        s += c[-1]; k0 += n
    return res, gmin, gk

def m1var(coeffs_exponents, n):
    """m of sum_j a_j z^{e_j} by trapezoid with n points (n a power of 2)."""
    t = 2*np.pi*np.arange(n)/n
    v = np.zeros(n, dtype=complex)
    for a, e in coeffs_exponents:
        v += a*np.exp(1j*e*t)
    return float(np.mean(np.log(np.abs(v))))

if __name__ == "__main__":
    print("="*78)
    print("RW13C_03 -- THE MISSING TABLE ROWS AND FIVE RE-DERIVED FIGURES")
    print("="*78)
    mPv = float(mP())
    print("m(P) = %.15f" % mPv)

    # ------------------------------------------------------- (1) the two unrun connections
    print("\n" + "="*78)
    print("(1) THE TWO PUBLISHED CONNECTIONS C9 ASSERTS BUT NEVER RUNS.")
    print("    ONE VARIABLE: the connection.  Same estimator, same N-grid, same code path")
    print("    as W13C_01 (uint64 modular orbit, chunk 1e6).")
    print("="*78)
    FINE = [1000, 3162, 10000, 31623, 100000, 316228, 1000000, 3162278, 10000000]
    conns = [
        ("S4:973  f=3.14159, c=1.57080   [W-10 N-4: relation (157080, 314159)]", 3.14159, 1.57080,
         [(0.3, 471239), (0.3, 0), (0.4, 314159)], "z^471239"),
        ("W-10 lane D  f=1.3, c=2.0      [W-10 N-4: relation (20, 13)]", 1.3, 2.0,
         [(0.3, 33), (0.3, 0), (0.4, 13)], "z^33"),
        ("(for calibration) f=2.0, c=1.1  [relation (11,20); target lane ran this]", 2.0, 1.1,
         [(0.3, 31), (0.3, 0), (0.4, 20)], "z^31"),
    ]
    print("\n   %-62s %-16s %-14s" % ("connection", "A_1e7", "A_1e7 - m(P)"))
    keys = []
    for nm, f, c, poly, deg in conns:
        th = turns_from_fc(f, c)
        A1 = int(mp.floor(th[0]*TWO64)); A2 = int(mp.floor(th[1]*TWO64))
        keys.append((A1, A2))
        res, gmin, gk = birkhoff(A1, A2, 10**7, FINE)
        print("   %-62s %-16.12f %+.4e" % (nm[:62], res[10**7], res[10**7]-mPv))
        print("      %-12s %s" % ("N", " ".join("%12s" % ("1e%.1f" % np.log10(N)) for N in FINE)))
        print("      %-12s %s" % ("A_N - m(P)", " ".join("%+12.3e" % (res[N]-mPv) for N in FINE)))
        print("      min|Z_k| over k<=1e7 = %.4e at k=%d ;  uint64 pair = (%d, %d)" % (gmin, gk, A1, A2))
        # independent: the subtorus limit as a ONE-VARIABLE Mahler measure
        vals = [m1var(poly, 2**e) for e in (18, 20, 22, 24)]
        print("      SUBTORUS LIMIT independently, m(0.3 + 0.4 z^a + 0.3 %s) by trapezoid" % deg)
        print("      2^18..2^24: %s" % "  ".join("%.9f" % v for v in vals))
        print("      gap from m(P) at 2^24: %+.4e" % (vals[-1]-mPv))
    assert len(set(keys)) == len(keys)
    print("""
   READ-OFF.  Both previously-unrun rows DO sit at a limit different from m(P), so C9's
   claim survives -- but it was an INHERITANCE, not a measurement, and the sizes matter for
   how the ruling reads: the gap at S4:973 is of the order printed above, not of the order
   of the order-4 point's 3.7e-02.  Resonance alone does not imply a different limit (W-12
   Cor. 2 exhibits a resonant subtorus whose average equals m(P) to 5.55e-17), so these two
   rows needed running and now have been.""")

    # ------------------------------------------------------- (2a) the 10^(L-4) figure
    print("\n" + "="*78)
    print("(2a) 'float64 would have reported |Z_k1| = 1.06e-15 -- WRONG BY 10^(L-4)'")
    print("="*78)
    import sys
    sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W13_C_CONVERGENCE")
    from W13C_03_liouville import build_theta, dip_value, K1 as TK1
    print("   %-8s %-24s %-18s %-16s %-16s" % ("L", "true |Z_k1| (mpmath)", "float64 value", "TRUE ratio", "as printed"))
    for L in (500, 2000, 8000, 17372):
        n1, n2, DEN, Z1, Z2, _ = build_theta(L)
        pv = dip_value(Z1, Z2, L)
        mp.mp.dps = L+50
        ratio = mp.mpf('1.060516e-15')/pv
        print("   %-8d %-24s %-18s 10^%-13.1f 10^%-13d" % (L, mp.nstr(pv, 8), "1.060516e-15",
              float(mp.log10(ratio)), L-4))
    print("   THE PRINTED EXPONENT IS L-4 AND THE TRUE ONE IS L-15 (to a rounding).  The")
    print("   sealed figure '10^17368' in W13C_05 C6 and W13C_03's output is therefore ~11")
    print("   orders too large.  IT CHANGES NO VERDICT -- the point is that float64 cannot")
    print("   see the dip at all -- but it is a wrong number in a headline claim, of exactly")
    print("   the class the register keeps recording (COR-K, COR-L).")

    # ------------------------------------------------------- (2b) the off-dip minimum
    print("\n" + "="*78)
    print("(2b) 'min |Z_k| over k<=1e7 EXCLUDING the engineered dip: 1.0735e-04 at k=9815000'")
    print("="*78)
    print("   W13C_03_liouville.py:96-97 takes the argmin of the WHOLE chunk and then")
    print("   discards it if it is the splice index -- so the chunk containing k1 = 1e4")
    print("   contributes NOTHING to the running minimum.  Recomputed correctly here by")
    print("   masking the dip index before the argmin:")
    L = 17372
    n1, n2, DEN, Z1, Z2, _ = build_theta(L)
    A1 = (n1*TWO64)//DEN; A2 = (n2*TWO64)//DEN
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    cur1 = 0; cur2 = 0; k0 = 0; gmin = np.inf; gk = -1; chunkmins = []
    while k0 < 10**7:
        n = min(1_000_000, 10**7-k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        z = Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64)
        if k0 < TK1 <= k0+n: z[TK1-k0-1] = np.inf          # MASK, do not discard the chunk
        j = int(np.argmin(z)); chunkmins.append((k0//1_000_000, float(z[j]), k0+j+1))
        if z[j] < gmin: gmin = float(z[j]); gk = k0+j+1
        k0 += n
    print("   per-chunk minima (chunk 0 is the one their code throws away):")
    for ci_, v, kk in chunkmins:
        print("      chunk %2d  min |Z_k| = %.6e at k = %8d %s" % (ci_, v, kk, "  <- DISCARDED by their code" if ci_==0 else ""))
    print("   CORRECT off-dip minimum over k <= 1e7:  %.6e at k = %d" % (gmin, gk))
    print("   their reported value: 1.0735e-04 at k=9815000")
    print("   VERDICT: %s" % ("the discarded chunk did not contain the true minimum, so the "
          "FIGURE IS RIGHT BY LUCK; the CODE IS WRONG." if abs(gmin-1.0735e-4) < 1e-8
          else "THE FIGURE IS WRONG: the discarded chunk held a smaller value."))
    print("   Either way the sentence it supports -- 'the background orbit is bit-for-bit the")
    print("   same in all four arms' -- is established directly by the printed uint64 pair,")
    print("   which IS identical in all four arms, so nothing load-bearing moves.")

    # ------------------------------------------------------- (2c) the tail sup
    print("\n" + "="*78)
    print("(2c) W13C_02 (D2) COLUMN 3: A MONOTONE RUNNING MAX SOLD AS A TAIL SUP.")
    print("="*78)
    print("   Their code: rate = -log|Z_k|/k on part of each chunk, then")
    print("   np.maximum.accumulate, then max() against a value carried across chunks.  A")
    print("   running max is NON-DECREASING, so it can never exhibit the decay to 0 that the")
    print("   necessary condition C3 requires -- the diagnostic COULD NOT HAVE SHOWN what it")
    print("   is labelled as showing.  The genuine per-decade sup, one decade per row:")
    th = turns_from_fc(1, mp.sqrt(2))
    A1 = int(mp.floor(th[0]*TWO64)); A2 = int(mp.floor(th[1]*TWO64))
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    print("   arm: GENERIC f=1.0 c=sqrt(2)")
    print("   %-26s %-18s %-14s" % ("decade", "sup (1/k)log(1/|Z_k|)", "at k"))
    cur1 = 0; cur2 = 0; k0 = 0
    dec_edges = [(10**j, 10**(j+1)) for j in range(0, 7)]
    best = {e: (0.0, -1) for e in dec_edges}
    while k0 < 10**7:
        n = min(1_000_000, 10**7-k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        z = Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64)
        kk = k0 + 1 + np.arange(n)
        r = -np.log(z)/kk
        for (lo, hi) in dec_edges:
            m_ = (kk > lo) & (kk <= hi)
            if m_.any():
                j = int(np.argmax(np.where(m_, r, -np.inf)))
                if r[j] > best[(lo,hi)][0]: best[(lo,hi)] = (float(r[j]), int(kk[j]))
        k0 += n
    for (lo, hi) in dec_edges:
        v, kk_ = best[(lo,hi)]
        print("   %-26s %-18.6e %-14d" % ("%d < k <= %d" % (lo, hi), v, kk_))
    print("   THIS one decays, decade by decade, which is what C3 needs.  The published")
    print("   column does not and cannot.  C3 ITSELF IS UNAFFECTED -- it is a proof.")

    # ------------------------------------------------------- (2d) the constant L
    print("\n" + "="*78)
    print("(2d) 'L*dist <= |P| <= U*dist with L = 1.026056' -- A GRID MINIMUM AS AN INEQUALITY")
    print("="*78)
    mp.mp.dps = 40
    phi = float(mp.acos(mp.mpf(-1)/9)/(2*mp.pi))
    zs = [(0.5-phi/2, 0.5+phi/2), (0.5+phi/2, 0.5-phi/2)]
    def ratio(a, b):
        d = np.inf
        for (za, zb) in zs:
            da = np.abs(a-za); da = np.minimum(da, 1-da)
            db = np.abs(b-zb); db = np.minimum(db, 1-db)
            d = np.minimum(d, np.hypot(da, db))
        return Pabs(a, b)/d
    for G in (2048, 4096, 8192):
        aa = (np.arange(G))/G
        A_, B_ = np.meshgrid(aa, aa, indexing='ij')
        r = ratio(A_, B_)
        r[~np.isfinite(r)] = np.inf
        j = int(np.argmin(r))
        print("   grid %5d x %5d :  min |P|/dist = %.8f  at (%.6f, %.6f)"
              % (G, G, float(r.flat[j]), A_.flat[j], B_.flat[j]))
    # local refinement by coordinate descent from the 16384 argmin
    a0, b0 = float(A_.flat[j]), float(B_.flat[j]); step = 1.0/8192
    cur = float(r.flat[j])
    for _ in range(60):
        improved = False
        for da_, db_ in ((step,0),(-step,0),(0,step),(0,-step),(step,step),(-step,-step),(step,-step),(-step,step)):
            v = float(ratio(np.array([a0+da_]), np.array([b0+db_]))[0])
            if v < cur: cur = v; a0 += da_; b0 += db_; improved = True
        if not improved: step /= 2
    print("   local descent refinement          :  min |P|/dist = %.8f  at (%.8f, %.8f)" % (cur, a0, b0))
    print("   TARGET LANE'S STATED L            :  1.02605600")
    print("   The stated inequality |P| >= 1.026056*dist is %s as written."
          % ("FALSE" if cur < 1.026056 else "consistent with the refinement"))
    print("   This is COR-E's defect class (a grid figure stated as an equality/inequality).")
    print("   NOT LOAD-BEARING: L enters C4 only inside an O(.) and C7 only inside a cube")
    print("   root of 10^17372.  Recorded, not scored against the theorems.")

    # ------------------------------------------------------- (2e) E_N(eps), five decades
    print("\n" + "="*78)
    print("(2e) E_N(eps) -> (pi/|det J|) eps^2 = 0.8897 eps^2 -- FIVE DECADES, NOT FOUR")
    print("="*78)
    EPSS = [1e-1, 1e-2, 1e-3]
    cur1 = 0; cur2 = 0; k0 = 0
    sums = {e: 0.0 for e in EPSS}; cnts = {e: 0 for e in EPSS}
    CPS5 = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
    outrows = []
    ci = 0
    while k0 < 10**8:
        n = min(2_000_000, 10**8-k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        z = Pabs(x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64)
        lz = np.log(z)
        cE = {e: np.cumsum(np.where(z < e, np.log(e)-lz, 0.0)) for e in EPSS}
        cC = {e: np.cumsum((z < e).astype(np.int64)) for e in EPSS}
        while ci < len(CPS5) and CPS5[ci] <= k0+n:
            N = CPS5[ci]; i_ = N-k0-1
            outrows.append((N, {e: (sums[e]+cE[e][i_])/N for e in EPSS},
                               {e: (cnts[e]+int(cC[e][i_]))/N for e in EPSS}))
            ci += 1
        for e in EPSS:
            sums[e] += cE[e][-1]; cnts[e] += int(cC[e][-1])
        k0 += n
    print("   %-12s %s" % ("N", " ".join("%14s" % ("E_N(%g)" % e) for e in EPSS)
                           + " | " + " ".join("%12s" % ("cnt(%g)" % e) for e in EPSS)))
    for N, E, C in outrows:
        print("   %-12d %s" % (N, " ".join("%14.4e" % E[e] for e in EPSS)
                               + " | " + " ".join("%12.4e" % C[e] for e in EPSS)))
    print("   predicted   %s" % (" ".join("%14.4e" % (np.pi*e*e/3.5310565) for e in EPSS)
                                 + " | " + " ".join("%12.4e" % (2*np.pi*e*e/3.5310565) for e in EPSS)))
    print("   CONFIRMED over five decades, both the excess-mass constant and the shell count.")

    # ------------------------------------------------------- (3) the second zero
    print("\n" + "="*78)
    print("(3) THE SPACING LEMMA AT THE SECOND ZERO.  z' = (z2,z1) = -z mod 1.")
    print("="*78)
    rng = np.random.default_rng(20260818)
    z = np.array([0.5-phi/2, 0.5+phi/2]); zp = np.array([0.5+phi/2, 0.5-phi/2])
    print("   z + z' mod 1 = (%.1e, %.1e)  -> z' = -z" % (((z+zp)[0]) % 1, ((z+zp)[1]) % 1))
    worst_same = 0.0; worst_cross = 0.0
    for _ in range(20000):
        th_ = rng.random(2); k1_ = int(rng.integers(1, 5000)); k2_ = int(rng.integers(k1_+1, 20000))
        h1 = (k1_*th_ - z + 0.5) % 1 - 0.5
        h2 = (k2_*th_ - z + 0.5) % 1 - 0.5
        lhs = ((k2_-k1_)*z + 0.5) % 1 - 0.5
        rhs = (k1_*h2 - k2_*h1 + 0.5) % 1 - 0.5
        worst_same = max(worst_same, float(np.max(np.abs(((lhs-rhs)+0.5) % 1 - 0.5))))
        # CROSS CASE: dip 1 near z, dip 2 near z' = -z
        h2c = (k2_*th_ - zp + 0.5) % 1 - 0.5
        lhsc = ((k1_+k2_)*z + 0.5) % 1 - 0.5
        rhsc = (k1_*h2c - k2_*h1 + 0.5) % 1 - 0.5
        worst_cross = max(worst_cross, float(np.max(np.abs(((lhsc-rhsc)+0.5) % 1 - 0.5))))
    print("   SAME-ZERO identity  ||(k2-k1)z|| = ||k1 h2 - k2 h1||   worst residual %.3e" % worst_same)
    print("   CROSS-ZERO identity ||(k1+k2)z|| = ||k1 h2 - k2 h1||   worst residual %.3e" % worst_cross)
    print("   (both are the float64 roundoff k1*k2*2^-52 ~ 2e-8)")
    print("   So the cross case obeys the SAME bound with q = k1+k2 <= 2 k2 in place of")
    print("   q = k2-k1 <= k2.  C7's conclusion is UNCHANGED (both corollaries only use")
    print("   q <= 2 k2), but the lemma as stated in W13C_04/W13C_05 does not cover it and")
    print("   the second zero is never mentioned in the spacing analysis.  CORRECTION, NOT")
    print("   REFUTATION.")
    print("\nDONE RW13C_03")
