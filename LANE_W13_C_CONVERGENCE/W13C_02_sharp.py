#!/usr/bin/env python3
"""
W13C_02 -- THE SHARP CONDITION, NUMERICALLY.
Three diagnostics, all at K1's REGISTERED pi = (0,3/10,3/10,2/5):

 (D1) THE EXACT DECOMPOSITION.  For every eps>0,
        A_N  =  (1/N) SUM f_eps(k theta)  -  E_N(eps),
        f_eps := max(log|P|, log eps)  CONTINUOUS,
        E_N(eps) := (1/N) SUM_{k<=N, |Z_k|<eps} log(eps/|Z_k|)  >= 0.
      Under H2 the first term -> INT f_eps -> m(P).  So convergence holds IFF
        lim_{eps->0} limsup_N E_N(eps) = 0.
      This is an IDENTITY, verified to machine precision, not a model.

 (D2) THE INHOMOGENEOUS APPROACH.  min_{k<=N} dist(k theta, Z(P)) as a function of N.
      For the sufficiency proof one wants a polynomial lower bound; for the necessity
      direction one wants max_k (1/k) log(1/|Z_k|) -> 0.  Both reported over 4+ decades.

 (D3) THE SHELL COUNT.  #{k<=N : dist_k < eps}/N against the area 2*pi*eps^2 -- the
      ingredient the sufficiency proof needs from equidistribution.

Seed 20260817.  Exact uint64 phase reduction.  Same arms, same code path, as W13C_01.
"""
import numpy as np, mpmath as mp
from W13C_01_central import (TWO64, mahler_jensen, zeros_on_torus, turns_from_fc,
                             to_u64, cubic_root, Pabs)

def orbit_chunks(A1, A2, Nmax, chunk=1_000_000):
    a1 = np.uint64(A1 % TWO64); a2 = np.uint64(A2 % TWO64)
    cur1 = 0; cur2 = 0; k0 = 0
    while k0 < Nmax:
        n = min(chunk, Nmax - k0)
        inc1 = np.full(n, a1, dtype=np.uint64); inc1[0] = np.uint64((cur1 + int(a1)) % TWO64)
        inc2 = np.full(n, a2, dtype=np.uint64); inc2[0] = np.uint64((cur2 + int(a2)) % TWO64)
        x1 = np.cumsum(inc1, dtype=np.uint64); x2 = np.cumsum(inc2, dtype=np.uint64)
        cur1 = int(x1[-1]); cur2 = int(x2[-1])
        yield k0, x1.astype(np.float64)/TWO64, x2.astype(np.float64)/TWO64
        k0 += n

def dist_to_Z(t1, t2, zs):
    d = np.full(t1.shape, np.inf)
    for (a, b) in zs:
        da = np.abs(t1 - float(a)); da = np.minimum(da, 1-da)
        db = np.abs(t2 - float(b)); db = np.minimum(db, 1-db)
        d = np.minimum(d, np.hypot(da, db))
    return d

if __name__ == "__main__":
    print("="*78)
    print("W13C_02 -- THE SHARP CONDITION.  pi = (0, 3/10, 3/10, 2/5).")
    print("="*78)
    mPv = float(mahler_jensen(50))
    zs = zeros_on_torus(40)
    print("m(P) = %.15f ;  Z(P) on T^2 = two simple zeros at" % mPv)
    for a,b in zs: print("     (%.15f, %.15f)" % (float(a), float(b)))

    mp.mp.dps = 60
    rho = cubic_root()
    arms = [
        ("GENERIC  f=1.0 c=sqrt(2)   [S4:603, the corpus's only generic connection]", turns_from_fc(1, mp.sqrt(2))),
        ("BADLY-APPROX A  (-2^(1/3), 4^(1/3))", ((-mp.cbrt(2)) % 1, (mp.cbrt(4)) % 1)),
        ("BADLY-APPROX B  (rho, rho^2), rho^3=rho+1", (rho % 1, (rho**2) % 1)),
    ]
    rng2 = np.random.default_rng(20260817)
    for s in range(2):
        w = rng2.integers(0, TWO64, size=2, dtype=np.uint64)
        arms.append(("HAAR-RANDOM draw %d" % s, (mp.mpf(int(w[0]))/TWO64, mp.mpf(int(w[1]))/TWO64)))
    arms.append(("CONTROL: S3/S4 HEADLINE f=2.0 c=1.1 (resonant, H2 FAILS)", turns_from_fc(2, mp.mpf('1.1'))))
    arms.append(("CONTROL: S1 PUBLISHED (order 4, H2 FAILS)", turns_from_fc(mp.pi, 3*mp.pi/2)))

    NMAX = 10_000_000
    CPS  = [10**3, 10**4, 10**5, 10**6, 10**7]
    EPSS = [1e-1, 1e-2, 1e-3, 1e-4]

    for name, th in arms:
        A1, A2 = to_u64(th[0]), to_u64(th[1])
        print("\n" + "-"*78)
        print("%s\n   theta(uint64) = (%d, %d)" % (name, A1, A2))
        # accumulators
        sum_f = 0.0
        sum_feps = {e: 0.0 for e in EPSS}
        sum_E    = {e: 0.0 for e in EPSS}
        cnt      = {e: 0   for e in EPSS}
        run_min_d = np.inf; run_max_rate = -np.inf; argmax_rate = -1
        rows = []
        ci = 0
        for k0, t1, t2 in orbit_chunks(A1, A2, NMAX):
            z = Pabs(t1, t2)
            d = dist_to_Z(t1, t2, zs)
            lz = np.log(z)
            # cumulative within chunk
            c_f = np.cumsum(lz)
            c_fe = {e: np.cumsum(np.maximum(lz, np.log(e))) for e in EPSS}
            c_E  = {e: np.cumsum(np.where(z < e, np.log(e) - lz, 0.0)) for e in EPSS}
            c_c  = {e: np.cumsum((z < e).astype(np.int64)) for e in EPSS}
            c_md = np.minimum.accumulate(d)
            kk = k0 + 1 + np.arange(len(lz))
            rate = np.where(kk*10 >= (k0+len(lz)), -lz/kk, -np.inf)  # TAIL sup: k > N/10
            c_mr = np.maximum.accumulate(rate)
            while ci < len(CPS) and CPS[ci] <= k0 + len(lz):
                N = CPS[ci]; i = N - k0 - 1
                row = dict(N=N,
                           A=(sum_f + c_f[i])/N,
                           feps={e: (sum_feps[e] + c_fe[e][i])/N for e in EPSS},
                           E={e: (sum_E[e] + c_E[e][i])/N for e in EPSS},
                           cnt={e: (cnt[e] + int(c_c[e][i]))/N for e in EPSS},
                           mind=min(run_min_d, float(c_md[i])),
                           maxrate=max(run_max_rate, float(c_mr[i])))
                rows.append(row); ci += 1
            sum_f += c_f[-1]
            for e in EPSS:
                sum_feps[e] += c_fe[e][-1]; sum_E[e] += c_E[e][-1]; cnt[e] += int(c_c[e][-1])
            run_min_d = min(run_min_d, float(c_md[-1]))
            j = int(np.argmax(rate))
            if rate[j] > run_max_rate: run_max_rate = float(rate[j]); argmax_rate = k0+j+1

        print("\n   (D1) THE IDENTITY  A_N = (1/N)SUM f_eps - E_N(eps).  max |residual| over all rows:")
        resid = max(abs(r['A'] - (r['feps'][e] - r['E'][e])) for r in rows for e in EPSS)
        print("        %.3e   (machine zero => the decomposition is an identity, not a model)" % resid)

        print("\n   (D1) E_N(eps) -- THE ONLY THING THAT CAN BREAK CONVERGENCE.  Must -> 0 as eps -> 0.")
        print("        %-10s %-14s %-14s %-14s %-14s" % ("N", "E_N(1e-1)", "E_N(1e-2)", "E_N(1e-3)", "E_N(1e-4)"))
        for r in rows:
            print("        %-10d %-14.3e %-14.3e %-14.3e %-14.3e"
                  % (r['N'], r['E'][1e-1], r['E'][1e-2], r['E'][1e-3], r['E'][1e-4]))

        print("\n   (D1) (1/N) SUM f_eps(k theta)  -- the CONTINUOUS part; converges with NO Diophantine input")
        print("        %-10s %-16s %-16s %-16s %-16s" % ("N", "eps=1e-1", "eps=1e-2", "eps=1e-3", "eps=1e-4"))
        for r in rows:
            print("        %-10d %-16.12f %-16.12f %-16.12f %-16.12f"
                  % (r['N'], r['feps'][1e-1], r['feps'][1e-2], r['feps'][1e-3], r['feps'][1e-4]))

        print("\n   (D2) min_{k<=N} dist(k theta, Z)  and  max_{k<=N} (1/k) log(1/|Z_k|)")
        print("        %-10s %-16s %-14s %-24s" % ("N", "min dist", "sqrt(N)*mindist", "tailsup (1/k)log(1/|Z_k|)"))
        for r in rows:
            print("        %-10d %-16.6e %-14.4f %-22.6e" % (r['N'], r['mind'], np.sqrt(r['N'])*r['mind'], r['maxrate']))
        print("        the third column is the TAIL sup over N/10 < k <= N, which is the quantity")
        print("        the NECESSARY CONDITION requires to go to 0:  limsup_k (1/k)log(1/|Z_k|) = 0.")

        print("\n   (D3) #{k<=N : |Z_k| < eps}/N  against the exact area |{|P|<eps}| = (2 pi/|det J|) eps^2")
        print("        %-10s %-14s %-14s %-14s %-14s" % ("N", "eps=1e-1", "eps=1e-2", "eps=1e-3", "eps=1e-4"))
        for r in rows:
            print("        %-10d %-14.4e %-14.4e %-14.4e %-14.4e"
                  % (r['N'], r['cnt'][1e-1], r['cnt'][1e-2], r['cnt'][1e-3], r['cnt'][1e-4]))
        print("        area (2pi/3.531056)eps^2: %-14.4e %-14.4e %-14.4e %-14.4e"
              % tuple(2*np.pi*e**2/3.5310565 for e in EPSS))
        print("        PREDICTED E_N(eps) -> (pi/|det J|) eps^2 = %.4e %.4e %.4e %.4e"
              % tuple(np.pi*e**2/3.5310565 for e in EPSS))
    print("\nDONE W13C_02")
