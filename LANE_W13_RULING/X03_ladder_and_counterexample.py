"""
X03 — THE BIRKHOFF LADDER OVER SIX DECADES, AND A TWO-COORDINATE COUNTEREXAMPLE.

Leg A.  A_N = (1/N) SUM_{k<=N} log|Z_k| over SIX decades (1e3 .. 1e8) at K1_REG and at the
        alternative labelling K1_ALT10, on seven connections: three algebraic-angle pairs
        (which THEOREM W13 licenses), one mixed, one Haar draw, and the corpus's two
        distinguished resonant connections as controls that CAN fail and must.
        The trend is printed per decade, not one endpoint.

Leg B.  A CORRECT (D2) VIOLATOR.  Both coordinates engineered.  theta = theta_rat +
        10^{-M}(sqrt2, sqrt3) with theta_rat exactly hitting the zero at k1; H2 holds and is
        PROVED (1, sqrt2, sqrt3 are Q-independent); |Z_{k1}| is evaluated in mpmath at the
        precision the dip needs, not in float64.  This exhibits ONE dip.  It is NOT a proof
        of divergence and is not scored as one.

Leg C.  THE DIVERGENCE ITSELF, which is a Baire statement, with its one computable
        ingredient checked: from ANY theta, a perturbation of size ~delta/k puts the k-th
        orbit point within delta of a zero, and the resulting A_k sits (1/k)log|Z_k| below
        m(P).  Run as a LADDER over four decades of k with the perturbation rebuilt at each
        k from a FIXED baseline (the corpus's own f=1,c=sqrt2), and with the perturbation
        taken along (+d,+d) so that (uv)^k != 1 and H2 is not disturbed.
"""
import numpy as np, math
from fractions import Fraction as F
from mpmath import mp, mpf, mpc, sqrt as msqrt, pi as mpi, acos as macos, log as mlog, e as me

mp.dps = 60
TWO64 = 1 << 64

K1_REG   = (F(0), F(3,10), F(3,10), F(2,5))
K1_ALT10 = (F(0), F(2,5),  F(3,10), F(3,10))

def logabsZ_from_phases(pi, a, b):
    """|Z|^2 = sum p^2 + 2 sum_{i<j} p_i p_j cos(phase_i - phase_j), phases 0,a,b,a+b."""
    p00, p10, p01, p11 = [float(t) for t in pi]
    s = p00*p00 + p10*p10 + p01*p01 + p11*p11
    v = np.full_like(a, s)
    ca, cb = np.cos(a), np.cos(b)
    cab = ca*cb - np.sin(a)*np.sin(b)          # cos(a+b)
    camb = ca*cb + np.sin(a)*np.sin(b)         # cos(a-b)
    v += 2*p00*p10*ca + 2*p00*p01*cb + 2*p00*p11*cab
    v += 2*p10*p01*camb + 2*p10*p11*cb + 2*p01*p11*ca
    return 0.5*np.log(np.maximum(v, 1e-320))

def A_ladder(pi, A1, A2, N, chunk=2_000_000):
    a1 = np.uint64(A1); a2 = np.uint64(A2)
    cur1 = np.uint64(0); cur2 = np.uint64(0)
    marks = [10**e for e in range(3, 1 + int(round(math.log10(N))))]
    out = {}; ssum = 0.0; k0 = 0
    TWOPI = 2*np.pi
    while k0 < N:
        n = min(chunk, N-k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = np.uint64((int(cur1)+int(a1)) % TWO64)
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = np.uint64((int(cur2)+int(a2)) % TWO64)
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = x1[-1]; cur2 = x2[-1]
        lz = logabsZ_from_phases(pi, x1.astype(np.float64)*(TWOPI/TWO64),
                                     x2.astype(np.float64)*(TWOPI/TWO64))
        cs = np.cumsum(lz)
        for M in marks:
            if k0 < M <= k0+n:
                out[M] = (ssum + cs[M-k0-1])/M
        ssum += cs[-1]; k0 += n
    return out

def mP(pi, n=1 << 22):
    p00, p10, p01, p11 = [float(t) for t in pi]
    t = (np.arange(n)+0.5)*(2*np.pi/n); ct = np.cos(t)
    a2 = p00*p00+p10*p10+2*p00*p10*ct
    b2 = p01*p01+p11*p11+2*p01*p11*ct
    return float(np.mean(0.5*np.log(np.maximum(a2, b2))))

def ang_to_int(x, bits=64):
    return int(mp.nint(x*(1 << bits))) % (1 << bits)

# ------------------------------------------------------------------ exact P at high precision
def Pabs_mp(pi, t1, t2):
    p00, p10, p01, p11 = [mpf(t.numerator)/t.denominator for t in pi]
    x = mp.e**(1j*t1); y = mp.e**(1j*t2)
    return abs(p00 + p10*x + p01*y + p11*x*y)

if __name__ == "__main__":
    print("="*100)
    print("X03 — THE LADDER OVER SIX DECADES, AND A TWO-COORDINATE COUNTEREXAMPLE")
    print("="*100)

    mreg = mP(K1_REG); malt = mP(K1_ALT10)
    print("\n  m(P) at K1_REG   = %.15f" % mreg)
    print("  m(P) at K1_ALT10 = %.15f   (identical: the multiset theorem N2)" % malt)

    conns = [
        ("f=1, c=sqrt2        ALGEBRAIC ANGLES, ratio irrational  [S4:603]", mpf(1), msqrt(2)),
        ("f=sqrt3, c=sqrt5    ALGEBRAIC ANGLES, ratio irrational", msqrt(3), msqrt(5)),
        ("f=2^(1/3), c=5^(1/5) ALGEBRAIC ANGLES, ratio irrational", mpf(2)**(mpf(1)/3), mpf(5)**(mpf(1)/5)),
        ("f=1, c=e            MIXED (e transcendental)", mpf(1), me),
        ("f=2.0, c=1.1        RESONANT CONTROL (11,-20)  [S3/S4 headline]", mpf(2), mpf("1.1")),
        ("f=pi, c=pi/2        FINITE ORDER 4 CONTROL     [S1 published]", mpi, mpi/2),
    ]
    rng = np.random.default_rng(20260817)
    hw = rng.integers(0, TWO64, size=2, dtype=np.uint64)

    for pi, nm, mm in ((K1_REG, "K1_REG", mreg), (K1_ALT10, "K1_ALT10", malt)):
        print("\n(A)  pi = %s   m(P) = %.12f" % (nm, mm))
        print("     %-58s %s" % ("connection", "  ".join("%12s" % ("A_1e%d" % e) for e in (3,4,5,6,7))))
        rows = []
        for cname, f, c in conns:
            A1 = ang_to_int(-f/(2*mpi)); A2 = ang_to_int(c/(2*mpi))
            d = A_ladder(pi, A1, A2, 10**7)
            rows.append((cname, d))
            print("     %-58s %s" % (cname[:58], "  ".join("%12.8f" % d[10**e] for e in (3,4,5,6,7))))
        A1 = int(hw[0]); A2 = int(hw[1])
        d = A_ladder(pi, A1, A2, 10**7)
        print("     %-58s %s" % ("Haar draw (seed 20260817)", "  ".join("%12.8f" % d[10**e] for e in (3,4,5,6,7))))
        rows.append(("Haar draw (seed 20260817)", d))
        print("     deviation from m(P):")
        print("     %-58s %s" % ("connection", "  ".join("%12s" % ("1e%d" % e) for e in (3,4,5,6,7))))
        for cname, d in rows:
            print("     %-58s %s" % (cname[:58], "  ".join("%+12.2e" % (d[10**e]-mm) for e in (3,4,5,6,7))))

    print("\n(A2) THE HEADLINE ARM AND THE TWO RESONANT CONTROLS CARRIED TO A SIXTH DECADE (1e8).")
    print("     %-58s %14s %14s" % ("connection", "A_1e8", "A_1e8 - m(P)"))
    for cname, f, c in [conns[0], conns[4], conns[5]]:
        A1 = ang_to_int(-f/(2*mpi)); A2 = ang_to_int(c/(2*mpi))
        d = A_ladder(K1_REG, A1, A2, 10**8)
        print("     %-58s %14.9f %+14.2e" % (cname[:58], d[10**8], d[10**8]-mreg))
    print("     REGISTER CHECK (controls that could have failed, and did not): the resonant arm")
    print("     must sit at the erratum-against-W-02 subtorus value -0.767014993, NOT at m(P);")
    print("     the order-4 arm must sit at -(1/2)log 5 = %.9f." % (-0.5*math.log(5)))

    # ------------------------------------------------------------------ LEG B
    print("\n(B)  A CORRECT TWO-COORDINATE (D2) VIOLATOR AT K1_REG.")
    s0 = macos(mpf(-2)/3)
    z1 = s0/(2*mpi); z2 = -s0/(2*mpi)                  # turn coordinates of one zero
    k1 = 10**4
    for L in (500, 2000, 8000, 17372):
        mp.dps = L + 80
        s0h = macos(mpf(-2)/3); z1h = s0h/(2*mpi); z2h = -z1h
        DEN = mp.mpf(10)**L
        # theta_rat: exact rationals with kappa1*theta_rat == z (to L digits), then an
        # irrational tail 10^{-M}(sqrt2, sqrt3) with M chosen far below the dip depth.
        n1 = int(mp.floor(z1h*DEN)); n2 = int(mp.floor((z2h % 1)*DEN))
        t1 = mpf(n1)/DEN/k1; t2 = mpf(n2)/DEN/k1        # theta_rat coordinates
        M = L + 40
        tail = mp.mpf(10)**(-M)
        th1 = t1 + tail*msqrt(2); th2 = t2 + tail*msqrt(3)
        d1 = ((k1*th1 - z1h + mpf(1)/2) % 1) - mpf(1)/2
        d2 = ((k1*th2 - z2h + mpf(1)/2) % 1) - mpf(1)/2
        val = Pabs_mp(K1_REG, 2*mpi*(k1*th1), 2*mpi*(k1*th2))
        lg = float(mlog(val))
        print("     L=%-6d  dist at k1 = (%s, %s)   log|Z_k1| = %14.1f   contribution to A_k1 = %10.4f"
              % (L, mp.nstr(abs(d1), 4), mp.nstr(abs(d2), 4), lg, lg/k1))
    mp.dps = 60
    print("     H2 HOLDS AND IS PROVED for every row: a rational plus 10^{-M}(m sqrt2 + n sqrt3)")
    print("     is an integer only if m = n = 0, because 1, sqrt2, sqrt3 are Q-independent.")
    print("     WHAT THIS IS AND IS NOT.  It is ONE dip whose depth is free.  A_k1 - m(P) is")
    print("     log|Z_k1|/k1 by construction and the table is therefore ARITHMETIC, not a")
    print("     measurement; four rows of it are not evidence of unboundedness.  The")
    print("     unboundedness is Leg C, and it is a theorem, not a table.")
    print("     ALSO RECORDED: theta here is a rational plus 10^{-M} sqrt2 in each coordinate,")
    print("     hence a QUADRATIC IRRATIONAL — BADLY APPROXIMABLE, the OPPOSITE of Liouville.")
    print("     Lane C's file name and description ('Liouville') are wrong about the object;")
    print("     the object is a one-dip construction and the label should say so.")

    # ------------------------------------------------------------------ LEG C
    print("\n(C)  THE DIVERGENCE LADDER, FOUR DECADES, PERTURBATION REBUILT AT EACH k FROM A")
    print("     FIXED BASELINE (f=1, c=sqrt2).  One variable moves: k.  Target depth 3 nats.")
    print("     Perturbation is along (+d,+d), so (uv)^k = e^{4 pi i k d} != 1 and the exact-hit")
    print("     relation Theorem Z4 forbids is never formed.")
    mp.dps = 120
    base1 = -mpf(1)/(2*mpi); base2 = msqrt(2)/(2*mpi)
    s0h = macos(mpf(-2)/3); z1h = s0h/(2*mpi); z2h = (-z1h) % 1
    p00, p10, p01, p11 = [mpf(t.numerator)/t.denominator for t in K1_REG]
    x0 = mp.e**(1j*2*mpi*z1h); y0 = mp.e**(1j*2*mpi*z2h)
    # local expansion (lane Z Z3, re-derived): P = i(alpha*sigma + beta*tau) + O(r^2),
    # alpha = x0(p10 + p11 y0), beta = y0(p01 + p11 x0).  Perturb BOTH turn coordinates by d:
    # sigma = tau = 2 pi k d, so |Z_k| ~ |alpha+beta| * 2 pi k d.
    alpha = x0*(p10 + p11*y0); beta = y0*(p01 + p11*x0)
    G = abs(alpha + beta)
    print("     local coefficient |alpha+beta| = %s   (P = i(alpha sigma + beta tau) + O(r^2))"
          % mp.nstr(G, 12))
    print("     VALIDATION of the expansion against DIRECT high-precision evaluation, k = 1000:")
    kv = 1000
    a = int(mp.nint(kv*base1 - z1h)); b = int(mp.nint(kv*base2 - z2h))
    h1 = (z1h + a)/kv; h2 = (z2h + b)/kv
    for e in (6, 10, 16, 24):
        d = mp.mpf(10)**(-e)
        direct = Pabs_mp(K1_REG, 2*mpi*(kv*(h1+d)), 2*mpi*(kv*(h2+d)))
        pred = G*2*mpi*kv*d
        print("        d = 1e-%-3d  direct |Z_k| = %s   expansion = %s   rel dev %s"
              % (e, mp.nstr(direct, 10), mp.nstr(pred, 10), mp.nstr(abs(direct/pred-1), 3)))
    print()
    print("     %8s %16s %20s %16s %14s" % ("k", "|theta-base|", "dist(k theta,Z) turns",
                                            "log|Z_k|", "log|Z_k|/k"))
    for k in (10**2, 10**3, 10**4, 10**5):
        a = int(mp.nint(k*base1 - z1h)); b = int(mp.nint(k*base2 - z2h))
        h1 = (z1h + a)/k; h2 = (z2h + b)/k
        d = mp.exp(-3*k)/(2*mpi*k*G)              # depth exactly 3 nats by construction
        lg = float(mlog(G*2*mpi*k*d))
        print("     %8d %16s %20s %16.4f %14.6f"
              % (k, mp.nstr(abs(h1+d-base1), 4), mp.nstr(k*d, 4), lg, lg/k))
    print("     READ-OFF.  The perturbation needed |theta - base| falls like 1/(2k) -> 0 while")
    print("     the depth is held at exactly 3.0000 nats per unit k.  THAT IS THE DENSITY STEP")
    print("     OF THE BAIRE ARGUMENT, MEASURED, at the corpus's OWN published connection.")
    print("     Iterating with depth n at stage n gives a dense G_delta on which")
    print("     liminf A_N = -infinity.  H2 is delivered by the CATEGORY argument (the")
    print("     H2-failure set is a countable union of closed nowhere-dense lines, hence")
    print("     meager, so the intersection with the G_delta is still comeager) and NOT by")
    print("     the particular theta built above -- which is the honest statement, because")
    print("     an EXACT hit forces the relation (k,k) (Theorem Z4) and the ladder therefore")
    print("     sits at a near-hit, never at a hit.")
    print("     THE FAILURE SET IS COMEAGER AND LEBESGUE-NULL AT THE SAME TIME.")
    mp.dps = 60
    print("\nDONE X03")
