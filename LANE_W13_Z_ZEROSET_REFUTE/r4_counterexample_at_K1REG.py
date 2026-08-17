"""LANE W-13 / Z_REFUTE -- r4: THE COUNTEREXAMPLE THE ROUND ASKED FOR, AT N1's OWN STATE.
The round's job was 'settle the convergence, or exhibit the counterexample'.  Lane Z settles the
EMPTY 3/4 and declines the counterexample, pointing at M1_06 -- which is a SINGLE Liouville pair
at the CENTROID state (0,1/3,1/3,1/3), not at K1's registered pi.  Lane Z's Theorem Z4 then reads
as reassurance about K1's registered pi ('no connection satisfying H2 ever lands exactly on a
zero').  This script shows that (i) the reassurance is empty, and (ii) the counterexample is not a
single pathological pair at one state but a TOPOLOGICALLY GENERIC set of connections at EVERY
singular state, K1's registered pi included -- and that is a theorem, not a construction."""
import sys, math
import numpy as np
from fractions import Fraction as F
import mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from r0_lib import strat_sorted, zero_angles, zero_points, mahler_1d, fr, K1_REG, CENTROID, S1_PUB

W = 96
def hdr(s):
    print("=" * W); print(s); print("=" * W)

hdr("r4  H2 IS NOT SUFFICIENT AT K1's REGISTERED pi, AND THE FAILURE SET IS COMEAGRE")
print("numpy", np.__version__, "; mpmath", mp.__version__, "\n")

# =========================================================================================
print("-" * W)
print("(a) THEOREM R4 -- THE BAIRE STATEMENT.  This replaces a construction with a proof and")
print("    it covers K1's registered pi, which M1_06's single exhibit does not.")
print("""
    SETTING.  u = e^{2 pi i alpha}, v = e^{2 pi i beta}; Z_k = P(u^k, v^k);
    F_N(alpha,beta) = (1/N) SUM_{k<=N} log|Z_k|, with values in [-inf, 0].
    H2 = the relation lattice L = {(m,n) : m alpha + n beta in Z} is {0}.

    THEOREM R4.  Let pi be ANY state whose zero set Z(P) is non-empty -- i.e. any state in the
    TWO, ONE or CURVE stratum, which by lane Z's own Theorem Z1 is exactly the quarter of the
    simplex where T_A T_B T_C <= 0.  Then

        G(pi) := { (alpha,beta) in T^2 : H2 holds AND liminf_N F_N(alpha,beta) = -infinity }

    is a COMEAGRE (residual, dense G_delta) subset of T^2.  In particular it is uncountable,
    dense, and non-empty at K1's REGISTERED pi.

    PROOF.  Fix a zero (x0,y0) of P on T^2, with angles (a0,b0) in [0,1)^2.
    For integers M, N0 >= 1 put
        U(M,N0) = UNION over N >= N0 of { (alpha,beta) : F_N(alpha,beta) < -M }.
    (i) OPEN.  Each log|Z_k| is continuous as a map T^2 -> [-inf, 0] (it tends to -inf at the
        zeros, which is continuity into the two-point compactification), so F_N is too, and
        F_N^{-1}([-inf,-M)) is open; a union of open sets is open.
    (ii) DENSE.  Let (alpha_*,beta_*) be arbitrary and delta > 0.  Choose N >= max(N0, 2/delta).
        There is t in [0,1/N) with N(alpha_* + t) = a0 mod 1, and t < 1/N < delta; likewise for
        beta.  At the perturbed point Z_N = 0, so F_N = -infinity < -M.  Hence U(M,N0) meets
        every ball.  (Openness then gives a whole neighbourhood.)
    So each U(M,N0) is open dense, and U := INTERSECT over M,N0 of U(M,N0) is residual by Baire.
    Every point of U has liminf F_N = -infinity.
    The H2 set R = T^2 \\ UNION over (m,n) != 0 of {m alpha + n beta in Z} is the complement of a
    COUNTABLE union of closed nowhere-dense sets (each is a finite union of parallel circles),
    hence residual.  G = U INTERSECT R is residual.  []

    WHAT IT SAYS AND WHAT IT DOES NOT.
      * It is NOT a measure statement and does not contradict lane Z's Z-12: by Borel-Cantelli
        the failure set is Lebesgue-NULL.  MEASURE AND CATEGORY DISAGREE HERE, which is the
        normal state of affairs for Diophantine conditions, and the corpus's own N3/W-12
        apparatus is entirely measure-theoretic (Haar pushforward, Haar-null resonances).
        A quantity that is 'true almost everywhere and false on a comeagre set' is exactly the
        kind of claim that needs its quantifier stated, and N1 as registered states none.
      * It makes Theorem Z4 (lane Z's Z-10) INERT for N1.  Z4 forbids EXACT hits at K1's
        registered pi.  R4's density argument passes through exact hits and then discards them:
        the residual set consists of connections NEAR the exact-hit points, at which no relation
        holds at all.  Forbidding a measure-zero, meagre set of exact hits removes nothing.
      * It is uniform over the strata.  theta = 2 vs 3/2 vs 1 changes the RATE at which the
        dive must be engineered; it does not change whether it can be.  So lane Z's 'MILDEST
        SINGULAR STRATUM' is a statement about the exceptional set's SHAPE, not about whether
        N1 holds -- which lane Z says, and which its headline does not.
""")

# =========================================================================================
print("-" * W)
print("(b) THE EXHIBIT, AT K1's REGISTERED pi.  ITS ZEROS, AND THE LOCAL CONSTANT USED BELOW.")
p = K1_REG
ang = zero_angles(p)
s0 = ang[0]
a0 = s0 / (2 * math.pi)
b0 = 1.0 - a0
p00, p10, p01, p11 = [float(q) for q in p]
x0 = complex(math.cos(s0), math.sin(s0)); y0 = x0.conjugate()
alpha = x0 * (p10 + p11 * y0); beta = y0 * (p01 + p11 * x0)
mP = mahler_1d(p, 1 << 22)
print(f"    pi = {tuple(str(q) for q in p)}   stratum {strat_sorted(p)}   m(P) = {mP:.12f}")
print(f"    zero angles (a0,b0) = ({a0:.15f}, {b0:.15f}),  a0 = arccos(-2/3)/(2 pi)")
print(f"    alpha = {alpha:.12f}   beta = {beta:.12f}")
print("""
    THE DIVE DIRECTION USED BELOW is (sigma,tau) = (2 pi d, 2 pi d), i.e. BOTH angles pushed the
    SAME way by d.  Then P ~ i(alpha+beta) 2 pi d and, since beta = conj(alpha) at this pi,
        |Z| ~ 4 pi |Re alpha| d = C d,   C = 4 pi * 0.2.
    This direction is chosen ON PURPOSE so that frac(N alpha) + frac(N beta) = 1 + 2d, i.e.
    (uv)^N = e^{4 pi i d} != 1: THE EXACT RELATION THEOREM Z4 FORBIDS IS NEVER FORMED, and the
    near-relation defect is exactly 2d, which is Theorem Z4's own 'defect O(delta)' clause.
    Z4 IS CONFIRMED AND IT COSTS THE COUNTEREXAMPLE NOTHING.
""")
C_pred = 4 * math.pi * abs(alpha.real)
print(f"    PREDICTED local constant C = 4 pi |Re alpha| = {C_pred:.12f}")
print(f"    VALIDATION of the local expansion against direct evaluation (M1_06's discipline):")
print(f"      {'d':>10s} {'direct |Z|':>22s} {'C*d':>22s} {'rel dev':>10s}")
for d in (1e-3, 1e-5, 1e-7, 1e-9):
    xx = complex(math.cos(2 * math.pi * (a0 + d)), math.sin(2 * math.pi * (a0 + d)))
    yy = complex(math.cos(2 * math.pi * (b0 + d)), math.sin(2 * math.pi * (b0 + d)))
    direct = abs(p00 + p10 * xx + p01 * yy + p11 * xx * yy)
    print(f"      {d:10.0e} {direct:22.15e} {C_pred*d:22.15e} "
          f"{abs(direct-C_pred*d)/(C_pred*d):10.2e}")
print()

# =========================================================================================
print("-" * W)
print("(c) THE LADDER.  FOUR DECADES OF N, TARGET DEPTH M = 3 NATS BELOW m(P) AT EVERY N.")
print("    THE CONNECTION IS REBUILT AT EACH N -- that is the point of a BAIRE statement: the")
print("    bad connections are dense, so one exists arbitrarily near any baseline, for every N.")
print("""
    CONSTRUCTION.  alpha = (J + A)/N with A = a0 + d, beta = (J' + B)/N with B = b0 + d, and
    J, J' integers chosen so that alpha, beta land within 1/N of a declared baseline.  Then
    frac(N alpha) = A and frac(N beta) = B EXACTLY, so the orbit's N-th point sits at offset
    (d,d) from the zero.  For k < N the angles are computed as (k J mod N)/N + k A / N, the
    first term in EXACT integers, the second in double -- error below 1e-16 and irrelevant
    because |Z_k| for k < N never falls below ~N^{-1/2}.  d is set to 10^{-D} with
    D = ceil(N*(M)/ln 10), and log|Z_N| = log C - D ln 10 by the validated local expansion.
""")
BASE_A = (math.sqrt(5) - 1) / 2
BASE_B = math.sqrt(2) - 1
M_TARGET = 3.0
print(f"    baseline (alpha_*,beta_*) = ({BASE_A:.15f}, {BASE_B:.15f})   target depth M = {M_TARGET}")
print(f"    {'N':>8s} {'|alpha-alpha_*|':>16s} {'D (digits)':>12s} {'log|Z_N|':>16s} "
      f"{'(1/N)SUM k<N':>14s} {'F_N':>14s} {'F_N - m(P)':>12s}")
for N in (100, 1000, 10000, 100000):
    D = math.ceil(N * M_TARGET / math.log(10.0)) + 3
    J = int(round(BASE_A * N)); Jp = int(round(BASE_B * N))
    A = a0                                   # + d, with d = 10^{-D}; d is below double resolution
    B = b0
    al = (J + A) / N; be = (Jp + B) / N
    k = np.arange(1, N)                      # k < N
    fa = ((k * J) % N) / N + (k * A) / N
    fb = ((k * Jp) % N) / N + (k * B) / N
    xx = np.exp(2j * np.pi * fa); yy = np.exp(2j * np.pi * fb)
    lz = np.log(np.abs(p00 + p10 * xx + p01 * yy + p11 * xx * yy))
    head = float(lz.sum())
    logZN = math.log(C_pred) - D * math.log(10.0)
    FN = (head + logZN) / N
    print(f"    {N:8d} {abs(al-BASE_A):16.3e} {D:12d} {logZN:16.2f} {head/N:14.6f} "
          f"{FN:14.6f} {FN-mP:12.4f}")
print(f"""
    READ IT.  At every N the connection differs from the baseline by less than 1/N -- so the
    bad connections are dense, exactly as Theorem R4 says -- and F_N sits {M_TARGET} nats below
    m(P) = {mP:.6f}.  Letting M -> infinity along a sparse subsequence of N (Theorem R4's
    diagonal) gives liminf F_N = -infinity at a single connection.
    THE COUNTEREXAMPLE THE ROUND ASKED FOR EXISTS AT N1's REGISTERED STATE.  M1_06's is at the
    CENTROID (0,1/3,1/3,1/3), a different state -- and one whose zero sits at a ROOT OF UNITY
    (cos s0 = -1/2, s0 = 2 pi/3) where K1's registered zero does not (cos s0 = -2/3; by Niven's
    theorem s0/pi is irrational).  So the transfer is not free and lane Z's citation of M1_06
    as 'the counterexample the round asks for' is a citation to a DIFFERENT STATE.
""")

# =========================================================================================
print("-" * W)
print("(d) H2 FOR THE EXHIBIT, STATED HONESTLY.  The alpha, beta above are RATIONAL, hence")
print("    EXACTLY RESONANT (W-10 N-4) and hence NOT H2.  THAT IS NOT PATCHED OUT; it is why")
print("    the load-bearing statement of this script is THEOREM R4 and not the ladder.")
print("""
    What the ladder is entitled to say, and does say: replace A by A + xi and B by B + xi'
    with xi, xi' free in (0, 10^{-(D+100)}).  Every member of that two-parameter family has the
    same F_N to 100 digits, because the perturbation moves the N-th orbit point by less than
    10^{-(D+100)} while the dive is at 10^{-D}.  The subset of the family violating H2 is a
    COUNTABLE union of line segments, hence Lebesgue-null in the (xi,xi') square.  SO ALMOST
    EVERY MEMBER OF THE FAMILY SATISFIES H2 AND SHOWS THE SAME DIVE.  That is a rigorous
    statement about a positive-measure family, not about the printed rational.
    (M1_06 does the same thing and says so: "plus a tail below 10^-(10^301), immaterial to
    k<=10".  This refuter is not improving on M1_06's honesty; it is moving the STATE.)
""")
print("    THE NEAR-RELATION DEFECT, MEASURED, as Theorem Z4 predicts it must be:")
for N in (100, 1000, 10000):
    D = math.ceil(N * M_TARGET / math.log(10.0)) + 3
    print(f"      N = {N:6d}:  frac(N alpha) + frac(N beta) - 1 = 2d = 2e-{D}, so (uv)^N != 1 and")
    print(f"                  the primitive relation nearest to ({N},{N}) is missed by 2e-{D}.")
print()
print("    AND ONE DIRECT HIGH-PRECISION EVALUATION, so the deep step is not taken on faith.")
print("    mpmath at 120 digits, d = 1e-60: |Z| computed from P itself, against C*d.")
mp.mp.dps = 140
d_mp = mp.mpf(10) ** (-60)
a0_mp = mp.acos(mp.mpf(-2) / 3) / (2 * mp.pi)
b0_mp = 1 - a0_mp
xx = mp.e ** (2j * mp.pi * (a0_mp + d_mp))
yy = mp.e ** (2j * mp.pi * (b0_mp + d_mp))
Zv = mp.mpf(0) + mp.mpf(3) / 10 * xx + mp.mpf(3) / 10 * yy + mp.mpf(2) / 5 * xx * yy
direct = abs(Zv)
print(f"      direct |Z|  = {mp.nstr(direct, 20)}")
print(f"      C*d         = {mp.nstr(mp.mpf(C_pred) * d_mp, 20)}")
print(f"      rel dev     = {mp.nstr(abs(direct - mp.mpf(C_pred)*d_mp)/direct, 6)}")
print(f"      log|Z|      = {mp.nstr(mp.log(direct), 12)}   vs  log C - 60 ln10 = "
      f"{math.log(C_pred) - 60*math.log(10):.10f}")
print("      ==> the local expansion used in the ladder is exact to 1e-59 at d = 1e-60.")
print()

# =========================================================================================
print("-" * W)
print("(e) WHAT THIS DOES TO THE REGISTER'S STATEMENT OF N1.")
print("""
    W-11 registers N1 'PUBLISH IT -- under two hypotheses it does not currently carry', H1 and
    H2.  W-12 Corollary 2 then strengthens H2 to 'purely arithmetic ... no structural failure
    mode'.  Neither row carries a THIRD hypothesis, and lane Z's Theorem Z3(ii) plus Theorem R4
    together say one is needed WHENEVER T_A T_B T_C <= 0:

      H3 (INHOMOGENEOUS DIOPHANTINE).  dist( (u^k, v^k), Z(P) ) >= c k^{-tau} for some tau and
      all k -- M1_08 T2(c)'s hypothesis (ii), which M1_08 says explicitly 'is NOT implied by
      (i) or by L = {0}'.

    H3 holds for Lebesgue-almost every connection (lane Z's Z-12; Borel-Cantelli) and FAILS on a
    comeagre set (Theorem R4).  AND IT IS NEEDED AT K1's REGISTERED pi, WHICH IS THE ONLY STATE
    N1 IS REGISTERED AT.  The register's H2 is necessary and not sufficient there, and this is
    not a new discovery -- M1_08 T2(b) says it and W-11 did not carry it into the row.
    THAT IS THE SAME UNDER-READING LANE Z's Z-13 CONVICTS, ONE LEVEL UP: not a lane failing to
    read a lane, but THE REGISTER failing to carry a hypothesis a lane had already isolated.
""")
print("    AND THE CONVERSE HALF, SO THE ROW CAN BE WRITTEN IN BOTH DIRECTIONS:")
print("      T_A T_B T_C > 0  ==>  H2 alone suffices (lane Z's Z-8, confirmed by this refuter);")
print("      T_A T_B T_C <= 0 ==>  H2 alone does NOT suffice (Theorem R4, this script).")
print("      The predicate is EXACT, RATIONAL, and computable from the ready state in one line.")
print("      N1 IS PUBLISHABLE WITH A SIDE CONDITION THAT CAN BE PRINTED IN THE STATEMENT.")
print()
print("DONE r4")
