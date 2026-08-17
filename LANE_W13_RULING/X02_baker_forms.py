"""
X02 — THE SHARP HYPOTHESIS, AND WHY IT IS SUPPLIED BY A PUBLISHED THEOREM AT THE CORPUS'S
      OWN PUBLISHED GENERIC CONNECTION.

The two hypotheses that decide N1 beyond H2 are
  (D1)  || m alpha + n beta ||  >=  c1 max(|m|,|n|)^{-tau1}      for all (m,n) != (0,0)
  (D2)  dist( (u^k, v^k), Z(P) )  >=  c2 k^{-tau2}               for all k >= 1
with alpha = -f/2pi, beta = c/2pi.  This script does three things:

 (a) WRITES BOTH OF THEM AS INHOMOGENEOUS LINEAR FORMS IN LOGARITHMS OF ALGEBRAIC NUMBERS
     with ALGEBRAIC coefficients, at f = 1.0, c = sqrt(2), pi = K1_REG — and certifies every
     ingredient the form needs (x0 algebraic: X01; Log(x0) and 2pi i Q-independent: Niven;
     the constant term non-zero).  The transcendence of u = e^{-i} is IRRELEVANT to the form:
     what enters is the ANGLE f, and f = 1 is algebraic.
 (b) MEASURES both quantities over enough decades to show a TREND, at f=1,c=sqrt2 and at
     three controls, one of which is engineered to violate (D2) and must be separated.
 (c) SWEEPS the claim that the licence is a property of the ANGLES being algebraic, by
     running four algebraic-angle pairs the corpus never published and one transcendental
     pair, and reporting the same diagnostics.

Phases are advanced in EXACT integer arithmetic modulo 2^64 (A = round(angle/2pi * 2^64));
the truncation is bounded and printed, and a two-limb 2^128 path is run on a subsample as a
cross-check that the single-limb figures are not an artefact of the surrogate rational.
"""
import numpy as np, math
from fractions import Fraction as F
from mpmath import mp, mpf, sqrt as msqrt, pi as mpi, cos as mcos, acos as macos

mp.dps = 60
TWO64 = 1 << 64
TWO128 = 1 << 128

def ang_to_int(x, bits=64):
    """x a real in units of turns (mpf); return round(x * 2^bits) mod 2^bits."""
    return int(mp.nint(x * (1 << bits))) % (1 << bits)

def orbit_min_dist(A1, A2, z1, z2, N, chunk=2_000_000):
    """min over k<=N of the sup-metric distance from (k*alpha, k*beta) to (z1,z2), where
       A1,A2,z1,z2 are turn coordinates as uint64 fixed point.  Returns per-decade minima."""
    a1 = np.uint64(A1); a2 = np.uint64(A2)
    cur1 = np.uint64(0); cur2 = np.uint64(0)
    out = {}
    marks = [10**e for e in range(1, 1 + int(round(math.log10(N))))]
    k0 = 0; best = np.inf
    while k0 < N:
        n = min(chunk, N - k0)
        i1 = np.full(n, a1, dtype=np.uint64); i1[0] = cur1 + a1
        i2 = np.full(n, a2, dtype=np.uint64); i2[0] = cur2 + a2
        x1 = np.cumsum(i1, dtype=np.uint64); x2 = np.cumsum(i2, dtype=np.uint64)
        cur1 = x1[-1]; cur2 = x2[-1]
        d1 = np.abs((x1.astype(np.float64) - float(z1)) / TWO64)
        d1 = np.minimum(d1, 1.0 - d1)
        d2 = np.abs((x2.astype(np.float64) - float(z2)) / TWO64)
        d2 = np.minimum(d2, 1.0 - d2)
        d = np.maximum(d1, d2)
        run = np.minimum.accumulate(d)
        for M in marks:
            if k0 < M <= k0 + n:
                out[M] = min(best, float(run[M - k0 - 1]))
        best = min(best, float(run[-1]))
        k0 += n
    return out

def hom_min(alpha, beta, M):
    """min over 0 < max(|m|,|n|) <= M of || m alpha + n beta ||, computed in float128-ish
       via mpmath-derived exact rationals scaled to 2^96."""
    S = 1 << 96
    A = int(mp.nint(alpha * S)); B = int(mp.nint(beta * S))
    best = None; arg = None
    for m in range(-M, M + 1):
        if m == 0:
            ns = [n for n in range(1, M + 1)]
        else:
            ns = range(-M, M + 1)
        for n in ns:
            if m == 0 and n == 0:
                continue
            v = (m * A + n * B) % S
            d = min(v, S - v)
            if best is None or d < best:
                best = d; arg = (m, n)
    return best / S, arg

# ---------------------------------------------------------------------------- the forms
def print_forms():
    s0 = macos(mpf(-2) / 3)
    print("  pi = K1_REG = (0, 3/10, 3/10, 2/5).  Z(P) = {(x0, conj x0), (conj x0, x0)},")
    print("  x0 = (-2 + i sqrt5)/3, a root of 3 z^2 + 4 z + 3 (X01, exact).  s0 = arg x0 =")
    print("  arccos(-2/3) = %s" % mp.nstr(s0, 25))
    print()
    print("  CONNECTION f = 1.0, c = sqrt(2) (S4:603; W-10 N-4: the ONLY generic connection")
    print("  this corpus publishes).  u = e^{-i f}, v = e^{i c}.  BOTH ARE TRANSCENDENTAL")
    print("  (Lindemann-Weierstrass).  THAT IS NOT THE QUESTION.  The forms are:")
    print()
    print("  (D2), first coordinate.   dist over the first circle from arg(u^k) = -k f to +-s0:")
    print("        Lambda = -k f -+ s0 + 2 pi n .        Multiply by i:")
    print("        i Lambda = (-i k f)  -+  Log(x0)  +  n (2 pi i)")
    print("                   \\_______/     \\______/       \\______/")
    print("                    beta_0        log of an      log of 1")
    print("                    ALGEBRAIC     ALGEBRAIC      (a period)")
    print("        beta_0 = -i k f is ALGEBRAIC because f = 1 IS ALGEBRAIC, and non-zero for")
    print("        k >= 1; height ~ k.  Log(x0) and 2 pi i are Q-LINEARLY INDEPENDENT because")
    print("        s0/2pi is irrational, which is NIVEN's theorem applied to cos s0 = -2/3.")
    print("        => BAKER's theorem on inhomogeneous linear forms gives |Lambda| > C H^{-K},")
    print("        H = max height ~ k, hence (D2) with tau2 = K.")
    print()
    print("  (D2), second coordinate is the same form with f replaced by c = sqrt(2), which is")
    print("        also ALGEBRAIC, and with the second zero coordinate y0 = conj(x0).")
    print()
    print("  (D1).  || m alpha + n beta || = |(-m + n sqrt2) - 2 pi q| / 2pi.  Multiply by i:")
    print("        i Lambda = i(-m + n sqrt2)  -  q (2 pi i) ,   beta_0 = i(-m+n sqrt2) != 0")
    print("        ALGEBRAIC of degree <= 4 and height << max(|m|,|n|)^2; one logarithm, 2 pi i.")
    print("        Same theorem.  => (D1) with some effective tau1.")
    print("        SANITY CHECK ON THE SHAPE:  at n = 0 this is |m - 2 pi q| and the bound is")
    print("        the classical EFFECTIVE IRRATIONALITY MEASURE OF pi (Mahler 1953: exponent")
    print("        42; Salikhov 2008: 7.6063...), which is exactly a bound of this shape.")
    print()
    print("  WHAT THE FORM NEEDS, AND WHAT IT DOES NOT.  It needs the two loop ANGLES f, c to")
    print("  be algebraic reals and the zero's coordinates to be algebraic (X01: they are, for")
    print("  every rational pi).  It does NOT need u, v to be algebraic — indeed if f and c")
    print("  are algebraic and non-zero then u and v are transcendental, and that is the")
    print("  regime in which beta_0 is algebraic and non-zero.")

# ---------------------------------------------------------------------------- main
if __name__ == "__main__":
    print("=" * 96)
    print("X02 — (D1) AND (D2) AS INHOMOGENEOUS LINEAR FORMS IN LOGARITHMS; AND THE MEASUREMENTS")
    print("=" * 96)

    print("\n(a) THE FORMS, WRITTEN OUT.")
    print_forms()

    # --------------------------------------------------------------- arms
    s0 = macos(mpf(-2) / 3)
    z1a = s0 / (2 * mpi); z2a = -s0 / (2 * mpi)          # zero  (x0, conj x0) angles / 2pi
    Z1 = ang_to_int(z1a); Z2 = ang_to_int(z2a)

    arms = []
    def add(name, f, c, note):
        al = -f / (2 * mpi); be = c / (2 * mpi)
        arms.append((name, al, be, note))
    add("f=1, c=sqrt2      [S4:603, the corpus's only generic pair]", mpf(1), msqrt(2), "ALGEBRAIC ANGLES")
    add("f=sqrt3, c=sqrt5  [algebraic, never published]", msqrt(3), msqrt(5), "ALGEBRAIC ANGLES")
    add("f=2^(1/3), c=5^(1/5) [algebraic, never published]", mpf(2)**(mpf(1)/3), mpf(5)**(mpf(1)/5), "ALGEBRAIC ANGLES")
    add("f=1, c=e          [e transcendental]", mpf(1), mp.e, "MIXED")
    add("f=pi, c=pi/2      [S4:973 shape: rational multiple of pi]", mpi, mpi/2, "RESONANT (H2 FAILS)")
    add("f=2.0, c=1.1      [S3/S4 headline]", mpf(2), mpf("1.1"), "RESONANT (H2 FAILS)")

    # engineered violator of (D2): choose alpha so that k1*alpha == z1 exactly for k1 = 10^4,
    # then perturb by an irrational tail so H2 still holds.
    k1 = 10**4
    Aviol = (Z1 // k1) + 1                      # integer fixed-point step making k1*A ~= Z1
    alpha_v = mpf(Aviol) / TWO64
    beta_v = msqrt(2) / (2 * mpi)
    arms.append(("ENGINEERED first-coordinate-only near-hit at k=1e4  [MY OWN VOID CONTROL]",
                 alpha_v, beta_v, "VOID: see the note under the table"))

    print("\n(b) MEASURED (D2):  min_{k<=N} dist((u^k,v^k), Z(P)) , sup metric on T^2, SEVEN")
    print("    DECADES.  For a POINT target in T^2 the equidistributed rate is ~N^{-1/2} (N")
    print("    points fill T^2 at spacing N^{-1/2}); a (D2) violation shows as a crash below it.")
    print("    [CORRECTION, RECORDED NOT PATCHED OUT: the first version of this script labelled")
    print("     the reference rate N^{-1}, which is the ONE-dimensional rate.  The numbers were")
    print("     unaffected; only the label was wrong.  M1_06 and lane C both state N^{-1/2}.]")
    print("    Phase step is exact mod 2^64; the surrogate-rational error at N=1e7 is")
    print("    <= 1e7 * 2^-65 = 2.7e-13 turns, four orders below the smallest figure printed.")
    N = 10**7
    print("\n    %-52s %s" % ("arm", "  ".join("%9s" % ("1e%d" % e) for e in range(1, 8))))
    res = {}
    for name, al, be, note in arms:
        A1 = ang_to_int(al); A2 = ang_to_int(be)
        d = orbit_min_dist(A1, A2, Z1, Z2, N)
        res[name] = d
        print("    %-52s %s" % (name[:52], "  ".join("%9.2e" % d[10**e] for e in range(1, 8))))
    print("\n    SAME TABLE AS sqrt(N) * mindist  (a BOUNDED row is the equidistributed rate; a")
    print("    row falling to 0 is a faster-than-equidistributed approach, i.e. toward a (D2)")
    print("    violation; a row GROWING is an orbit whose closure misses the zero altogether):")
    print("    %-52s %s" % ("arm", "  ".join("%9s" % ("1e%d" % e) for e in range(1, 8))))
    for name, al, be, note in arms:
        d = res[name]
        print("    %-52s %s" % (name[:52], "  ".join("%9.3f" % (math.sqrt(10**e) * d[10**e]) for e in range(1, 8))))

    print("\n    MY OWN CONTROL, RECORDED VOID RATHER THAN PATCHED OUT.  The last arm was built to")
    print("    violate (D2) at k = 1e4 and DOES NOT: I engineered only the FIRST coordinate, and")
    print("    dist is the SUP over both, so the second coordinate holds the distance up.  Its")
    print("    row is therefore indistinguishable from the good arms and PROVES NOTHING.  A")
    print("    correct two-coordinate violator is built in X03 and does crash.  The defect is")
    print("    exactly the class W-08's isolation audit names: an arm that could not have shown")
    print("    what its label claims.")
    print("\n(c) MEASURED (D1):  min over 0 < max(|m|,|n|) <= M of ||m alpha + n beta||, and the")
    print("    same scaled by M^2 (the heuristic exponent: three integers m,n,q of size M give")
    print("    ~M^3 values spread over a range ~M, hence a typical minimum ~M^{-2}).")
    print("    %-52s %s" % ("arm", "  ".join("%12s" % ("M=%d" % M) for M in (10, 30, 100, 300, 1000))))
    for name, al, be, note in arms:
        row = []
        for M in (10, 30, 100, 300, 1000):
            v, arg = hom_min(al, be, M)
            row.append(v)
        print("    %-52s %s" % (name[:52], "  ".join("%12.3e" % v for v in row)))
    print("    scaled by M^2:")
    for name, al, be, note in arms:
        row = []
        for M in (10, 30, 100, 300, 1000):
            v, arg = hom_min(al, be, M)
            row.append(v * M * M)
        print("    %-52s %s" % (name[:52], "  ".join("%12.3f" % v for v in row)))
    print("\n    A ZERO (or a value collapsing to the 2^-96 floor) in this row IS H2 failing:")
    for name, al, be, note in arms:
        v, arg = hom_min(al, be, 400)
        print("    %-52s min=%9.2e at (m,n)=%-12s H2: %s"
              % (name[:52], v, str(arg), "FAILS" if v < 1e-25 else "holds up to M=400"))

    print("\n(d) CROSS-CHECK OF THE PHASE PATH.  Single-limb (2^64) against exact Python ints")
    print("    and against a 2^128 two-limb representation, on the f=1,c=sqrt2 arm.")
    al = -mpf(1) / (2 * mpi); be = msqrt(2) / (2 * mpi)
    A64 = ang_to_int(al, 64); A128 = ang_to_int(al, 128)
    worst_int = 0.0; worst_128 = 0.0
    a1 = np.uint64(A64)
    x = np.cumsum(np.full(200000, a1, dtype=np.uint64), dtype=np.uint64)
    for k in (1, 2, 7, 1000, 99991, 200000):
        exact = (k * A64) % TWO64
        worst_int = max(worst_int, abs(int(x[k - 1]) - exact) / TWO64)
        v128 = (k * A128) % TWO128
        worst_128 = max(worst_128, abs(v128 / TWO128 - exact / TWO64))
    print("    max |cumsum - (k*A) mod 2^64| over k<=2e5      : %.3e turns" % worst_int)
    print("    max |2^128 path - 2^64 path| over the same k    : %.3e turns" % worst_128)
    print("    => the single-limb path IS the exact orbit of the surrogate rational, and the")
    print("       surrogate is within 2.7e-13 turns of the true angle at k = 1e7.")
    print("\nDONE X02")
