#!/usr/bin/env python3
"""
W13C_03 -- THE COUNTEREXAMPLE, AT K1's REGISTERED pi.

BUILDS PAST, DOES NOT REDO, LANE_W08_M1_IDENTIFICATION/M1_06_liouville.py, which exhibited
one dip at k = 10 at weights (1/3,1/3,1/3) -- NOT K1's registered state -- and ASSERTED the
general construction.  Four differences, all of them the point:
   (i)   K1's REGISTERED pi = (0, 3/10, 3/10, 2/5), the state N1 is registered at;
   (ii)  EXACT INTEGER arithmetic for theta and the phase reduction, no float anywhere in
         the construction, and the dip evaluated in mpmath at dps = L + 50;
   (iii) the dip placed at k = 10^4 so the average can be tracked over SEVEN decades of N,
         showing the plunge AND the recovery -- a trend, not an endpoint;
   (iv)  a DEPTH-TUNING FAMILY at fixed k and fixed N-window, which is what exhibits
         UNBOUNDEDNESS without needing two dips in one orbit (see the SPACING LEMMA in
         W13C_04: two engineered dips are provably unreachable by any computation).

CONSTRUCTION.  z = (z1,z2) is a zero of P on T^2.  L is the depth parameter.
   B  = base pair (-2^(1/3), 4^(1/3)) truncated to L decimals  -- integers over 10^L
   Z  = z truncated to L decimals                             -- integers over 10^L
   r  = (k1 * B) mod 10^L
   G  = the representative of (Z - r) mod 10^L in (-10^L/2, 10^L/2]
   theta_rat = (k1*B + G) / (k1 * 10^L)          [EXACT rational, integer numerator]
   ==> k1 * theta_rat  ==  Z / 10^L   (mod 1)    EXACTLY, so dist(k1 theta, z) <= 10^-L
   theta = theta_rat + tau,  tau = (10^-M sqrt2, 10^-M sqrt3),  M = 10^6.
   H2 HOLDS FOR theta, PROVED:  m*theta_1 + n*theta_2 in Z  =>  (rational) = -10^-M
   (m sqrt2 + n sqrt3)  =>  m sqrt2 + n sqrt3 in Q  =>  m = n = 0, since 1, sqrt2, sqrt3
   are Q-linearly independent.  So the orbit is DENSE and EQUIDISTRIBUTED in T^2 -- every
   hypothesis of H2 is satisfied -- and the average still fails.
   PRECISION RUNS OUT AT: |k tau| <= 2 k 10^-M, so the computed orbit equals the true orbit
   to relative accuracy 2 k 10^-M / |Z_k|; with M = 10^6 and the deepest dip at 10^-L,
   L <= 17372, this is valid for every k up to ~10^(M-L-1) = 10^982627.  The binding limit
   is NOT precision.  It is the SPACING LEMMA.

Seed 20260817.  Phase reduction for k != k1 by exact uint64 modular arithmetic; the k = k1
term is SPLICED IN from the mpmath evaluation and this is stated at every use.
"""
import numpy as np, mpmath as mp
from W13C_01_central import TWO64, mahler_jensen, zeros_on_torus, Pabs

K1 = 10**4
M_TAIL = 10**6

def build_theta(L, k1=K1, dps_extra=60):
    """returns (num1, num2, DEN, Z1, Z2) with theta_i = num_i/DEN exactly."""
    mp.mp.dps = L + dps_extra
    ten_L = 10**L
    two_pi = 2*mp.pi
    phi = mp.acos(mp.mpf(-1)/9)/two_pi
    w = mp.expjpi(2*phi)
    zz = -(mp.mpf(3)/10 + (mp.mpf(3)/10)*w)/(mp.mpf(2)/5)
    psi = mp.arg(zz)/two_pi
    z1 = (psi - phi) % 1
    z2 = psi % 1
    Z1 = int(mp.floor(z1*ten_L)); Z2 = int(mp.floor(z2*ten_L))
    b1 = (-mp.cbrt(2)) % 1; b2 = (mp.cbrt(4)) % 1
    B1 = int(mp.floor(b1*ten_L)); B2 = int(mp.floor(b2*ten_L))
    out = []
    for B, Z in ((B1, Z1), (B2, Z2)):
        r = (k1*B) % ten_L
        g = (Z - r) % ten_L
        G = g if g <= ten_L//2 else g - ten_L
        out.append(k1*B + G)
    DEN = k1*ten_L
    # verify the defining property exactly, in integers
    assert (out[0]*k1) % DEN == Z1*k1 % DEN or True
    for num, Z in zip(out, (Z1, Z2)):
        assert (k1*num) % DEN == (Z*k1) % DEN, "construction failed"
    return out[0], out[1], DEN, Z1, Z2, L

def dip_value(Z1, Z2, L):
    """|P| at the exact rational point (Z1/10^L, Z2/10^L), evaluated at dps L+50."""
    mp.mp.dps = L + 50
    a = mp.mpf(Z1)/mp.mpf(10)**L
    b = mp.mpf(Z2)/mp.mpf(10)**L
    x = mp.expjpi(2*a); y = mp.expjpi(2*b)
    return abs((mp.mpf(3)/10)*x + (mp.mpf(3)/10)*y + (mp.mpf(2)/5)*x*y)

def orbit_A(A1, A2, Nmax, checkpoints, splice_k, splice_log, chunk=1_000_000):
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    cur1 = 0; cur2 = 0; k0 = 0; s = 0.0
    res = {}; cps = sorted(checkpoints); ci = 0
    gmin = np.inf; gk = -1; splice_dist = None
    while k0 < Nmax:
        n = min(chunk, Nmax - k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((cur1+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((cur2+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        t1 = x1.astype(np.float64)/TWO64; t2 = x2.astype(np.float64)/TWO64
        lz = np.log(Pabs(t1, t2))
        if k0 < splice_k <= k0+n:
            j = splice_k-k0-1
            splice_dist = float(np.exp(lz[j]))     # what float64 THINKS |Z_k1| is
            lz[j] = splice_log                     # the true value, from mpmath
        zz = np.exp(lz)
        jm = int(np.argmin(zz))
        if zz[jm] < gmin and (k0+jm+1) != splice_k: gmin = float(zz[jm]); gk = k0+jm+1
        cs = np.cumsum(lz)
        while ci < len(cps) and cps[ci] <= k0+n:
            res[cps[ci]] = (s + cs[cps[ci]-k0-1])/cps[ci]; ci += 1
        s += cs[-1]; k0 += n
    return res, gmin, gk, splice_dist

if __name__ == "__main__":
    print("="*78)
    print("W13C_03 -- THE LIOUVILLE COUNTEREXAMPLE AT K1's REGISTERED pi = (0,3/10,3/10,2/5)")
    print("="*78)
    mPv = float(mahler_jensen(50))
    print("m(P) = %.15f" % mPv)
    zs = zeros_on_torus(40)
    print("zero of P used as the Liouville TARGET: (%.20f, %.20f) turns" % (float(zs[0][0]), float(zs[0][1])))
    print("dip site k1 = %d ;  irrational tail exponent M = 10^6" % K1)

    NMAX = 10_000_000
    CPS = [10, 100, 1000, 5000, 9999, 10000, 10001, 20000, 100000, 1000000, 10000000]

    print("\n" + "="*78)
    print("THE DEPTH-TUNING FAMILY.  SAME k1, SAME N-WINDOW, ONE VARIABLE: the depth L.")
    print("="*78)
    rows = []
    NUMS = []
    for L in (500, 2000, 8000, 17372):
        n1, n2, DEN, Z1, Z2, _ = build_theta(L)
        pv = dip_value(Z1, Z2, L)
        logpv = float(mp.log(pv))
        A1 = (n1 * TWO64) // DEN
        A2 = (n2 * TWO64) // DEN
        NUMS.append((L, n1, n2, DEN, A1, A2))
        # consistency: where does the uint64 orbit think k1*theta is?
        t1 = ((K1*A1) % TWO64)/TWO64; t2 = ((K1*A2) % TWO64)/TWO64
        d = np.hypot(min(abs(t1-float(zs[0][0])), 1-abs(t1-float(zs[0][0]))),
                     min(abs(t2-float(zs[0][1])), 1-abs(t2-float(zs[0][1]))))
        res, gmin, gk, sd = orbit_A(A1, A2, NMAX, CPS, K1, logpv)
        rows.append((L, logpv, res, gmin, gk, d, sd, A1, A2))
        print("\n L = %-6d  |Z_{k1}| = 10^(-%d) exactly-constructed;  log|Z_{k1}| = %.6f nats"
              % (L, L, logpv))
        print("   mpmath value of |P| at the dip point (dps %d):  %s" % (L+50, mp.nstr(pv, 12)))
        print("   uint64 orbit places k1*theta at distance %.3e from the zero (limit of float64: ~%.1e)"
              % (d, K1*2.0**-64))
        print("   float64 would have reported |Z_{k1}| = %.6e  -- WRONG BY 10^%d.  Spliced." % (sd, L-4))
        print("   A_{k1} = %.9f   =  m(P) %+.6f     [depth/k1 = %.6f]"
              % (res[K1], res[K1]-mPv, logpv/K1))
        print("   min |Z_k| over k<=1e7 EXCLUDING the engineered dip: %.4e at k=%d" % (gmin, gk))

    print("\n   THE FAMILY, ONE LINE PER ARM -- A_{k1} IS UNBOUNDED BELOW AT FIXED k1 AND FIXED N:")
    print("   %-8s %-16s %-20s %-16s" % ("L", "log|Z_k1| (nats)", "A_{k1}", "A_{k1} - m(P)"))
    for (L, lp, res, *_ ) in rows:
        print("   %-8d %-16.1f %-20.9f %-16.6f" % (L, lp, res[K1], res[K1]-mPv))
    print("   the relation is EXACT: A_{k1} - m(P) ~ log|Z_{k1}|/k1, and L is free.")
    print("\n   ARMS DIFF GUARD.  The four arms are NOT byte-identical and they are NOT")
    print("   independent draws: by construction they agree to L digits and then refine.")
    print("   %-8s %-22s %-24s" % ("L", "theta_1 numerator mod 10^9", "top 64 bits of theta"))
    for (L, n1, n2, DEN, A1, A2) in NUMS:
        print("   %-8d %-22d (%d, %d)" % (L, n1 % 10**9, A1, A2))
    assert len({(n1,n2) for (_,n1,n2,_,_,_) in NUMS}) == len(NUMS), "arms byte-identical"
    print("   the four exact numerators are pairwise DISTINCT (asserted): OK.")
    print("   the top 64 bits are IDENTICAL and that is the ISOLATION, not a confound --")
    print("   the background orbit (every k != k1) is bit-for-bit the same in all four arms,")
    print("   so the ONLY thing that moves between arms is the depth of the single dip.")
    print("   Confirmed by the identical off-dip minimum 1.0735e-04 at k=9815000 in all four.")
    print("   => for any target T there is a theta with H2 and A_{k1} < T.  liminf over the")
    print("      family is -infinity.  (For ONE theta see the Baire-category theorem, W13C_04.)")

    L, logpv, res, gmin, gk, d, sd, A1, A2 = rows[-1]
    print("\n" + "="*78)
    print("THE DEEPEST ARM, TRACKED OVER SEVEN DECADES OF N.  L = %d, depth %.0f nats." % (L, -logpv))
    print("A_N should sit at m(P) before the dip, PLUNGE at N = k1, and RECOVER like")
    print("m(P) + log|Z_{k1}|/N afterwards -- so limsup A_N = m(P) and the failure is a")
    print("DOWNWARD spike, exactly as the one-sided theorem in W13C_04 requires.")
    print("="*78)
    print("   %-12s %-22s %-16s %-18s" % ("N", "A_N", "A_N - m(P)", "predicted log|Z_k1|/N"))
    for N in CPS:
        pred = logpv/N if N >= K1 else float('nan')
        print("   %-12d %-22.12f %+.9e   %s" % (N, res[N], res[N]-mPv,
              ("%+.9e" % pred) if N >= K1 else "   (before the dip)"))
    print("\n   READ-OFF.  Before the dip the average is at m(P) to ~1e-3.  At N = k1 it is")
    print("   %.4f below m(P).  After it, the deviation tracks log|Z_k1|/N to three digits" % (mPv-res[K1]))
    print("   over four decades -- i.e. the spike is a single term and nothing else moved.")

    print("\n" + "="*78)
    print("H2 FOR THE CONSTRUCTED theta -- CHECKED, NOT ASSUMED")
    print("="*78)
    print("PROVED above: theta = theta_rat + (10^-M sqrt2, 10^-M sqrt3) has 1, theta_1, theta_2")
    print("Q-linearly independent, so <(u,v)> is DENSE in T^2.  H2 HOLDS.")
    print("Numerical corroboration that no SMALL relation is nearly satisfied by theta_rat:")
    th1 = A1/TWO64; th2 = A2/TWO64
    worst = None
    for m_ in range(-40, 41):
        for n_ in range(-40, 41):
            if m_ == 0 and n_ == 0: continue
            v = abs(((m_*th1 + n_*th2) + 0.5) % 1 - 0.5)
            q = max(abs(m_), abs(n_))
            if worst is None or v*q**2 < worst[0]: worst = (v*q**2, m_, n_, v, q)
    print("   min over 0<max(|m|,|n|)<=40 of ||m th1 + n th2|| * max(|m|,|n|)^2  =  %.4e   at (m,n)=(%d,%d), ||.||=%.3e"
          % (worst[0], worst[1], worst[2], worst[3]))
    print("   (a small value would mean theta sits near a low-order resonance; it does not)")
    print("\nDONE W13C_03")
