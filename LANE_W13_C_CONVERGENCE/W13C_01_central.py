#!/usr/bin/env python3
"""
W13C_01 -- THE CENTRAL COMPUTATION.
Does (1/N) SUM_{k<=N} log|P(u^k,v^k)| reach m(P) at K1's REGISTERED pi = (0,3/10,3/10,2/5),
where P DOES vanish on T^2?

ONE VARIABLE: theta = (arg u, arg v)/2pi.  Everything else byte-identical across arms.
Seed 20260817.  Phase reduction exact (uint64 modular).  Four+ decades of N, trend reported.
"""
import numpy as np, mpmath as mp
from fractions import Fraction

TWO64 = 1 << 64

# ----------------------------------------------------------------------------- pi and m(P)
P00, P10, P01, P11 = Fraction(0), Fraction(3,10), Fraction(3,10), Fraction(2,5)
COEF = np.array([0.0, 0.3, 0.3, 0.4])

def mahler_jensen(dps=50):
    """m(P) = (1/2pi) INT log max(|p00+p10 e^it|,|p01+p11 e^it|) dt, split at the exact
    crossing so each piece is smooth.  Here p00=0 so branch1 == 0.3 identically and
    branch2 = |0.3+0.4 e^it|, |.|^2 = 0.25+0.24 cos t; they cross at cos t = -2/3."""
    mp.mp.dps = dps
    t0 = mp.acos(mp.mpf(-2)/3)
    # full-circle integral of log|0.3+0.4 e^it| is 2pi log(0.4)  (Jensen, |0.4|>|0.3|)
    g = lambda t: mp.log(mp.mpf('0.3')) - mp.log(mp.mpf('0.25') + mp.mpf('0.24')*mp.cos(t))/2
    corr = mp.quad(g, [t0, mp.pi])            # region where 0.3 is the max
    val = mp.log(mp.mpf('0.4')) + corr/mp.pi  # (1/2pi)*2*corr = corr/pi
    return val

def mahler_grid(n):
    """independent check: plain trapezoid on the (continuous, non-smooth-at-2-points) max."""
    t = 2*np.pi*np.arange(n)/n
    b1 = np.abs(P00.__float__() + float(P10)*np.exp(1j*t))
    b2 = np.abs(float(P01) + float(P11)*np.exp(1j*t))
    return float(np.mean(np.log(np.maximum(b1, b2))))

# ----------------------------------------------------------------------------- zero set
def zeros_on_torus(dps=40):
    """P = e(a)[p10 + p11 e(b) + p01 e(b-a)] with pi=(0,.3,.3,.4):
       0.3 + 0.3 w + 0.4 z = 0, |w|=|z|=1  =>  |0.4|=|0.3(1+w)| => cos(2 pi phi) = -1/9.
       Returns the two zeros (a,b) in turns."""
    mp.mp.dps = dps
    two_pi = 2*mp.pi
    phi = mp.acos(mp.mpf(-1)/9)/two_pi            # in turns, phi = (b-a)
    out = []
    for s in (+1, -1):
        ph = s*phi
        w = mp.e**(1j*two_pi*ph)
        zz = -(mp.mpf('0.3') + mp.mpf('0.3')*w)/mp.mpf('0.4')
        psi = mp.arg(zz)/two_pi                    # = b
        a = (psi - ph) % 1
        b = psi % 1
        out.append((a, b))
    return out

def Pabs(a, b):
    x = np.exp(2j*np.pi*a); y = np.exp(2j*np.pi*b)
    return np.abs(0.3*x + 0.3*y + 0.4*x*y)

# ----------------------------------------------------------------------------- theta arms
def turns_from_fc(f, c, dps=60):
    """u = conj(W_F) = e^{-i f}, v = W_C = e^{i c}  =>  theta = (-f/2pi, c/2pi) mod 1."""
    mp.mp.dps = dps
    return ((-mp.mpf(f)/(2*mp.pi)) % 1, (mp.mpf(c)/(2*mp.pi)) % 1)

def to_u64(x, dps=60):
    mp.mp.dps = dps
    return int(mp.floor(mp.mpf(x) * TWO64)) % TWO64

def cubic_root(dps=60):
    """real root of x^3 - x - 1 (plastic number); Q(rho) is cubic with a complex embedding,
       so (rho, rho^2) is badly approximable (Cassels/Davenport)."""
    mp.mp.dps = dps
    return mp.findroot(lambda x: x**3 - x - 1, mp.mpf('1.3247'))

def cubic_root2(dps=60):
    """real root of x^3 + x - 1 (~0.6823); another cubic field with a complex embedding."""
    mp.mp.dps = dps
    return mp.findroot(lambda x: x**3 + x - 1, mp.mpf('0.68'))

# ----------------------------------------------------------------------------- the run
def birkhoff(A1, A2, Nmax, checkpoints, chunk=1_000_000):
    """exact uint64 modular orbit; returns dict N -> A_N, plus min |Z_k| and argmin."""
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    scale = 1.0/TWO64
    s = 0.0; k0 = 0
    cur1 = np.uint64(0); cur2 = np.uint64(0)
    res = {}; cps = sorted(checkpoints)
    ci = 0
    gmin = np.inf; gargmin = -1
    while k0 < Nmax:
        n = min(chunk, Nmax - k0)
        inc1 = np.full(n, a1, dtype=np.uint64); inc1[0] = np.uint64((int(cur1) + int(a1)) % TWO64)
        inc2 = np.full(n, a2, dtype=np.uint64); inc2[0] = np.uint64((int(cur2) + int(a2)) % TWO64)
        x1 = np.cumsum(inc1, dtype=np.uint64)
        x2 = np.cumsum(inc2, dtype=np.uint64)
        cur1 = x1[-1]; cur2 = x2[-1]
        t1 = x1.astype(np.float64)*scale
        t2 = x2.astype(np.float64)*scale
        ex = np.exp(2j*np.pi*t1); ey = np.exp(2j*np.pi*t2)
        z = np.abs(0.3*ex + 0.3*ey + 0.4*ex*ey)
        j = int(np.argmin(z))
        if z[j] < gmin: gmin = float(z[j]); gargmin = k0 + j + 1
        lz = np.log(z)
        cs = np.cumsum(lz)
        while ci < len(cps) and cps[ci] <= k0 + n:
            idx = cps[ci] - k0 - 1
            res[cps[ci]] = (s + cs[idx]) / cps[ci]
            ci += 1
        s += cs[-1]
        k0 += n
    return res, gmin, gargmin

def exact_check(theta_mp, ks, A1, A2, dps=80):
    """256-bit-equivalent exact check of the uint64 reduction at selected k."""
    mp.mp.dps = dps
    out = []
    for k in ks:
        t1_true = (k*theta_mp[0]) % 1
        t2_true = (k*theta_mp[1]) % 1
        t1_u64 = ((k*A1) % TWO64)/TWO64
        t2_u64 = ((k*A2) % TWO64)/TWO64
        d = max(abs(float(t1_true) - t1_u64), abs(float(t2_true) - t2_u64))
        out.append((k, d))
    return out

if __name__ == "__main__":
    print("="*78)
    print("W13C_01 -- THE CENTRAL COMPUTATION AT K1's REGISTERED pi = (0, 3/10, 3/10, 2/5)")
    print("="*78)

    mP = mahler_jensen(50)
    print("\nm(P) by the Jensen reduction, split at the exact crossing cos t = -2/3:")
    print("   m(P) = %s" % mp.nstr(mP, 20))
    print("   corpus value (REGISTER erratum vs W-02 / M1_07): -0.767507880358")
    print("   |difference| = %.3e" % abs(float(mP) - (-0.767507880358)))
    for n in (2**12, 2**16, 2**20, 2**24):
        print("   plain trapezoid n=%-9d  %.15f   dev %.3e" % (n, mahler_grid(n), mahler_grid(n)-float(mP)))

    print("\nDOES P VANISH ON T^2 AT THIS pi?  W-01's convex-hull criterion: max(0.4) <= 0.6 = sum of others -> YES.")
    zs = zeros_on_torus(40)
    for i,(a,b) in enumerate(zs):
        print("   zero %d at (a,b) = (%s, %s) turns   |P| there = %.3e"
              % (i+1, mp.nstr(a,25), mp.nstr(b,25), Pabs(float(a), float(b))))
    # simplicity of the zeros: gradient magnitude
    for i,(a,b) in enumerate(zs):
        a=float(a); b=float(b)
        x = np.exp(2j*np.pi*a); y = np.exp(2j*np.pi*b)
        # P = 0.3x + 0.3y + 0.4xy ; dP/da = 2pi i (0.3x + 0.4xy), dP/db = 2pi i (0.3y + 0.4xy)
        dPa = 2j*np.pi*(0.3*x + 0.4*x*y); dPb = 2j*np.pi*(0.3*y + 0.4*x*y)
        J = np.array([[dPa.real, dPb.real],[dPa.imag, dPb.imag]])
        print("   zero %d: real Jacobian of P, |det| = %.6f, smallest singular value = %.6f"
              % (i+1, abs(np.linalg.det(J)), np.linalg.svd(J, compute_uv=False)[-1]))
        print("           (nonzero => SIMPLE zero => |P| ~ (const)*dist locally)")
    # global linear lower bound |P| >= L*dist  (needed by the sufficiency proof)
    rng = np.random.default_rng(20260817)
    G = 4096
    aa, bb = np.meshgrid(np.arange(G)/G, np.arange(G)/G, indexing='ij')
    vals = Pabs(aa, bb)
    d = np.full_like(vals, np.inf)
    for (a,b) in zs:
        da = np.abs(aa-float(a)); da = np.minimum(da, 1-da)
        db = np.abs(bb-float(b)); db = np.minimum(db, 1-db)
        d = np.minimum(d, np.hypot(da, db))
    m = d > 1e-9
    L = float(np.min(vals[m]/d[m])); U = float(np.max(vals[m]/d[m]))
    print("   GLOBAL TWO-SIDED BOUND on a %dx%d grid:  L*dist <= |P(x)| <= U*dist,  L = %.6f, U = %.6f"
          % (G, G, L, U))
    print("   analytic upper bound: Lip(P) <= 2 pi * sqrt((0.3+0.4)^2+(0.3+0.4)^2) = %.6f  (>= U, OK: %s)"
          % (2*np.pi*np.hypot(0.7,0.7), 2*np.pi*np.hypot(0.7,0.7) >= U))
    print("   (L>0 confirms both zeros are SIMPLE and is the constant the sufficiency proof needs;")
    print("    U<inf is the constant the necessity and Baire-category proofs need)")

    # ------------------------------------------------------------------- the arms
    mp.mp.dps = 60
    rho = cubic_root()
    arms = []
    arms.append(("GENERIC (S4:603, W-10 N-4: the ONLY generic connection the corpus publishes) f=1.0 c=sqrt(2)",
                 turns_from_fc(1, mp.sqrt(2)), "H2 PROVED (Lindemann-Weierstrass)"))
    arms.append(("BADLY-APPROX A  theta = (-2^(1/3), 4^(1/3)) mod 1   [Cassels/Davenport cubic]",
                 ((-mp.cbrt(2)) % 1, (mp.cbrt(4)) % 1), "H2 PROVED (cubic field)"))
    arms.append(("BADLY-APPROX B  theta = (rho, rho^2) mod 1, rho^3=rho+1  [Cassels/Davenport cubic]",
                 (rho % 1, (rho**2) % 1), "H2 PROVED (cubic field)"))
    sig = cubic_root2()
    arms.append(("BADLY-APPROX C  theta = (sig, sig^2) mod 1, sig^3+sig=1  [Cassels/Davenport cubic]",
                 (sig % 1, (sig**2) % 1), "H2 PROVED (cubic field)"))
    rng2 = np.random.default_rng(20260817)
    for s in range(4):
        # Haar-random on T^2 with FULL 64-bit entropy in each coordinate
        w = rng2.integers(0, TWO64, size=2, dtype=np.uint64)
        arms.append(("HAAR-RANDOM seed 20260817 draw %d" % s,
                     (mp.mpf(int(w[0]))/TWO64, mp.mpf(int(w[1]))/TWO64), "H2 a.s."))
    arms.append(("CONTROL: S1 PUBLISHED  W_F=-1, W_C=-i  (f=pi, c=3pi/2; finite order 4)",
                 turns_from_fc(mp.pi, 3*mp.pi/2), "H2 FAILS (order 4)"))
    arms.append(("CONTROL: S3/S4 HEADLINE f=2.0, c=1.1  (exactly resonant, relation (11,20))",
                 turns_from_fc(2, mp.mpf('1.1')), "H2 FAILS (rank-1 subtorus)"))

    print("\nARMS DIFF GUARD -- the uint64 pair actually fed to the estimator, per arm:")
    keys = []
    for name, th, hh in arms:
        A1, A2 = to_u64(th[0]), to_u64(th[1])
        keys.append((A1, A2))
        print("   %-96s  A = (%20d, %20d)" % (name[:96], A1, A2))
    assert len(set(keys)) == len(keys), "TWO ARMS ARE BYTE-IDENTICAL -- ZERO-VARIABLE CONTROL"
    print("   all %d arms pairwise distinct: OK  (no zero-variable control)" % len(keys))

    NMAX = 10_000_000
    CPS = [10**3, 10**4, 10**5, 10**6, 10**7]
    FINE = [1000, 3162, 10000, 31623, 100000, 316228, 1000000, 3162278, 10000000]

    print("\n" + "="*78)
    print("A_N = (1/N) SUM_{k<=N} log|Z_k|      m(P) = %.12f" % float(mP))
    print("FOUR DECADES OF N (10^3 .. 10^7).  TREND, not an endpoint.")
    print("="*78)
    for (name, th, hh), (A1, A2) in zip(arms, keys):
        res, gmin, gk = birkhoff(A1, A2, NMAX, sorted(set(CPS+FINE)))
        print("\n%s\n   %s" % (name, hh))
        print("   %-12s %-22s %-14s %-14s" % ("N", "A_N", "A_N - m(P)", "|dev|*sqrt(N)"))
        for N in FINE:
            dev = res[N] - float(mP)
            print("   %-12d %-22.15f %+.6e   %.4f" % (N, res[N], dev, abs(dev)*np.sqrt(N)))
        # trend statement
        d3, d7 = abs(res[10**3]-float(mP)), abs(res[10**7]-float(mP))
        # decay exponent from a least-squares fit of log|dev| on log N over the four decades
        xs = np.log10(np.array(FINE, dtype=float))
        ys = np.log10(np.array([max(abs(res[N]-float(mP)), 1e-300) for N in FINE]))
        slope = float(np.polyfit(xs, ys, 1)[0])
        print("   TREND over the four decades 10^3 -> 10^7: |dev| %.3e -> %.3e   (ratio %.1f)"
              % (d3, d7, d3/max(d7,1e-300)))
        print("   TREND as an exponent: |A_N - m(P)| ~ N^(%+.3f)  (least squares over the 9 checkpoints)" % slope)
        print("   min_{k<=1e7} |Z_k| = %.6e at k = %d" % (gmin, gk))
        ec = exact_check(th, [10**3, 10**5, 10**7], A1, A2)
        print("   uint64-vs-exact phase deviation: " + "  ".join("k=%d: %.2e" % (k,d) for k,d in ec))

    print("\nDONE W13C_01")
