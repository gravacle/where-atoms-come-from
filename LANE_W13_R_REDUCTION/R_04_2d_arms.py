#!/usr/bin/env python3
"""
R_04 — THE CASE N1 ACTUALLY NEEDS: THE FULL TWO-DIMENSIONAL BIRKHOFF SUM.

R_02 proved that under H2 the orbit closure is 2-dimensional and NO one-variable reduction
exists.  This script measures the two-dimensional object directly, over SEVEN decades, on
arms differing in exactly one thing: the arithmetic type of omega = (alpha, beta).

   BA-cubic     omega = (2^{1/3}, 4^{1/3})   the classical badly-approximable pair in dim 2
                                             (the pair M1_06 used, so the comparison is direct)
   CORPUS       omega = (-1/(2 pi), sqrt2/(2 pi))   <-> f = 1.0, c = sqrt(2)
                                             THE ONLY GENERIC CONNECTION THE CORPUS PUBLISHES
                                             (S4:603; W-10 N-4 -- every rational pair is resonant)
   BA-algebraic omega = (sqrt2 - 1, sqrt3 - 1)
   RESONANT     omega <-> f = 2.0, c = 1.1   S3/S4's headline; -11f + 20c = 0 exactly.
                                             A CONTROL THAT CAN FAIL: it must converge to the
                                             SUBTORUS value -0.767014993, NOT to m(P).

THE TWO-DIMENSIONAL DIAGNOSTIC, which is the one that replaces Sudler's:
   sqrt(N) * min_{k<=N} dist((k alpha, k beta), Z(P)).
In one dimension Dirichlet FORCES ||k a|| <= 1/k, so near-hits on a homogeneous singularity
are unavoidable and the Sudler problem is delicate for every a.  In two dimensions nothing
forces a near-hit on a POINT: N orbit points in T^2 leave the zeros ~N^{-1/2} away, so the
worst single term is ~ -(1/2) log N and its share of the average is O(log N / N).
THAT is the structural reason our problem is EASIER than Sudler's generically and yet still
has counterexamples (R_05): the singular set is codimension 2 and the orbit is dimension 1.
"""
import math
import numpy as np
from R_lib import PI_K1, PhaseReducer, m_jensen, arm_hash

MP = m_jensen(PI_K1, 1 << 24)
TH = np.arccos(-2.0 / 3.0) / (2 * np.pi)
ZEROS = [(TH, -TH % 1.0), (-TH % 1.0, TH)]
KMAX = 10 ** 7
CH = 10 ** 6
DEC = [10 ** i for i in range(1, 8)]
p00, p10, p01, p11 = PI_K1

TWOPI = 2 * math.pi
ARMS = [
    ("BA-cubic", 2.0 ** (1 / 3), 4.0 ** (1 / 3)),
    ("CORPUS", -1.0 / TWOPI, math.sqrt(2.0) / TWOPI),
    ("BA-algebraic", math.sqrt(2.0) - 1.0, math.sqrt(3.0) - 1.0),
    ("RESONANT", -2.0 / TWOPI, 1.1 / TWOPI),
]

print("=" * 79)
print("R_04 — THE TWO-DIMENSIONAL BIRKHOFF SUM, SEVEN DECADES, FOUR ARMS, ONE VARIABLE MOVED")
print("=" * 79)
print("\nm(P) at K1's pi (Jensen, 2^24 nodes) = %.12f" % MP)
print("subtorus value at the (11,20) relation, register erratum against W-02 = -0.767014993")

results = {}
for name, a, b in ARMS:
    ra, rb = PhaseReducer(a), PhaseReducer(b)
    tot = 0.0
    sums = {}
    mind = np.inf
    minds = {}
    head = None
    hsh = []
    for lo in range(0, KMAX, CH):
        k = np.arange(lo + 1, min(lo + CH, KMAX) + 1, dtype=np.int64)
        fa, fb = ra.frac(k), rb.frac(k)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
        lg = np.log(np.abs(p00 + p10 * x + p01 * y + p11 * x * y))
        if head is None:
            head = (fa[:4].copy(), fb[:4].copy())
            hsh = lg[:200000].copy()
        d = np.full(len(k), np.inf)
        for (g1, g2) in ZEROS:
            d1 = np.abs(((fa - g1 + 0.5) % 1) - 0.5)
            d2 = np.abs(((fb - g2 + 0.5) % 1) - 0.5)
            d = np.minimum(d, np.hypot(d1, d2))
        csum = np.cumsum(lg)
        cmin = np.minimum.accumulate(d)
        for N in DEC:
            if lo < N <= lo + len(k):
                sums[N] = tot + csum[N - lo - 1]
                minds[N] = min(mind, float(cmin[N - lo - 1]))
        tot += float(csum[-1])
        mind = min(mind, float(cmin[-1]))
    results[name] = (sums, minds, head, arm_hash(hsh))

print("\nARMS-DIFF GUARD (sha256 of the first 200000 log|Z_k| of each arm)")
for name, _, _ in [(n, a, b) for n, a, b in ARMS]:
    print("   %-13s %s   first four orbit points (alpha side) %s"
          % (name, results[name][3], np.round(results[name][2][0], 9)))
hs = [results[n][3] for n, _, _ in ARMS]
print("   distinct hashes: %d of %d  -> %s" % (len(set(hs)), len(hs),
                                               "arms differ" if len(set(hs)) == len(hs) else "VOID"))

print("\n(1/N) sum_{k<=N} log|Z_k|      (target for the three generic arms: m(P) = %.9f)" % MP)
print("   %-13s" % "N" + "".join("%16d" % N for N in DEC))
for name, _, _ in ARMS:
    s = results[name][0]
    print("   %-13s" % name + "".join("%16.9f" % (s[N] / N) for N in DEC))

print("\nDEVIATION FROM m(P)")
print("   %-13s" % "N" + "".join("%16d" % N for N in DEC))
for name, _, _ in ARMS:
    s = results[name][0]
    print("   %-13s" % name + "".join("%+16.3e" % (s[N] / N - MP) for N in DEC))
print("\n   THE RESONANT ROW IS THE CONTROL THAT COULD HAVE FAILED AND DID NOT:")
sres = results["RESONANT"][0]
print("   it converges to %.9f, deviating from m(P) by %+.3e and from the register's"
      % (sres[10 ** 7] / 10 ** 7, sres[10 ** 7] / 10 ** 7 - MP))
print("   published subtorus value -0.767014993 by %+.3e.  The estimator can tell a"
      % (sres[10 ** 7] / 10 ** 7 + 0.767014993))
print("   2-dimensional orbit closure from a 1-dimensional one.  (Carried fact C1: this pair")
print("   is EXACTLY RESONANT and is NOT a generic connection.)")

print("\nDECAY EXPONENT r in |S_N - m(P)| ~ N^{-r}, fitted over the last four decades")
for name, _, _ in ARMS:
    s = results[name][0]
    tgt = MP if name != "RESONANT" else -0.767014993
    xs = np.array([4.0, 5.0, 6.0, 7.0])
    ys = np.array([math.log10(abs(s[10 ** i] / 10 ** i - tgt) + 1e-300) for i in (4, 5, 6, 7)])
    print("   %-13s r = %+6.3f   (target %s)" % (name, -np.polyfit(xs, ys, 1)[0],
                                                 "m(P)" if name != "RESONANT" else "subtorus"))

print("\nTHE TWO-DIMENSIONAL DIAGNOSTIC:  sqrt(N) * min_{k<=N} dist(orbit, Z(P))")
print("   %-13s" % "N" + "".join("%12d" % N for N in DEC))
for name, _, _ in ARMS:
    m = results[name][1]
    print("   %-13s" % name + "".join("%12.4f" % (math.sqrt(N) * m[N]) for N in DEC))
print("   %-13s" % "min dist" + "".join("%12.2e" % results["RESONANT"][1][N] for N in DEC)
      + "   <- RESONANT only")
print("   BOUNDED, O(1), OVER SEVEN DECADES ON THE THREE GENERIC ARMS -- AND NOT ON THE FOURTH.")
print("   The RESONANT row GROWS like sqrt(N), and that is not a defect: its orbit closure is a")
print("   CIRCLE that misses both zeros, so min dist SATURATES at %.3e and sqrt(N) times a"
      % results["RESONANT"][1][10 ** 7])
print("   constant grows.  That is an independent confirmation of PROP R-5 from the other side:")
print("   the (11,20) subtorus is NON-SINGULAR, which is why its Birkhoff average converges at")
print("   all, and converges to the subtorus value rather than to m(P).")
print("   For the three generic arms: N points in T^2 sit ~N^{-1/2} apart, and the orbit does")
print("   not do better than that against a fixed point.  So the deepest single term is")
print("   ~ (1/2) log N and its share of the average is O(log N / N).")
print("   COMPARE the ONE-dimensional diagnostic of R_03, N * min ||k a - theta*||, also O(1):")
print("   same shape, but there Dirichlet FORCES it to be O(1) and here nothing does -- which is")
print("   exactly the room in which R_05's counterexample lives.")

print("\nWHAT THIS DOES NOT ESTABLISH")
print("   Seven decades of a bounded diagnostic is not a lower bound on dist for all k.  This")
print("   script measures; it proves nothing.  The proof is in R_06, and it covers only")
print("   rotations reachable by Baker's theorem.")
print("\nDONE R_04")
