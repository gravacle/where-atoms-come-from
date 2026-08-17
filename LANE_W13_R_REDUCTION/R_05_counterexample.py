#!/usr/bin/env python3
"""
R_05 — THE COUNTEREXAMPLE, AT K1's OWN REGISTERED READY STATE.

WHAT THIS ADDS TO M1_06, WHICH IS NOT REDONE.  M1_06 (sealed, W-08) exhibited a Liouville pair
defeating equidistribution at pi = (0,1/3,1/3,1/3), where P's torus zeros sit at the RATIONAL
points (1/3,2/3) and (2/3,1/3) and a sparse-DIGIT construction reaches them.  At K1's OWN
registered pi = (0,0.3,0.3,0.4) the zeros are at (theta*, -theta*) with theta* IRRATIONAL
(R_01 Z2), and no digit construction reaches an irrational target.  This lane's construction
is different in kind -- a MONOTONE LADDER (R_ladder) -- and its non-resonance proof is
different in kind too: it turns on the STRICT POSITIVITY of the ladder's offsets, not on
sparse digits.

THEOREM R-3.  At pi = (0, 0.3, 0.3, 0.4) there is an explicit omega = (alpha, beta) in T^2 with
    L(omega) = {(m,n) in Z^2 : m alpha + n beta in Z} = {0}      (so H2 HOLDS, and the orbit is
                                                                  dense and equidistributed)
and
    liminf_N  (1/N) sum_{k<=N} log|Z_k|  =  -infinity,       while m(P) = -0.7675078804.
Hence H2 IS NOT SUFFICIENT FOR N1, at the corpus's own published ready state.

PROOF.
(a) Let gamma1 = theta*, gamma2 = 1 - theta*, so (e(gamma1), e(gamma2)) is a zero of P (R_01).
    Note gamma1 + gamma2 = 1, an integer.
(b) Build alpha from gamma1 and beta from gamma2 on a COMMON rung sequence k_1 | k_2 | ...
    by R_ladder's monotone construction.  Then for every j,
        frac(k_j alpha) = gamma1 + eps_j,   frac(k_j beta) = gamma2 + eps'_j,
        0 < eps_j, eps'_j < 2 k_j / k_{j+1}.
(c) DIVERGENCE.  |Z_k| <= 1 for every k (carried fact C5, W-08), so every term of the sum is
    <= 0 and the running average is bounded above by any single term over N:
        S_{k_j} <= (1/k_j) log|Z_{k_j}|.
    |dP/dtheta1| = 2 pi |p10 + p11 y| <= 2 pi(p10+p11) and |dP/dtheta2| <= 2 pi(p01+p11),
    both <= 2 pi at K1's pi, so
        |Z_{k_j}| <= 2 pi (eps_j + eps'_j) <= 2 pi * 4 k_j / k_{j+1},
    whence S_{k_j} <= -(1/k_j) log( k_{j+1} / (8 pi k_j) ).
    Choosing k_{j+1} = k_j * M_j with M_j >= exp(D_j k_j) gives S_{k_j} <= -D_j + O(1/k_j).
    Take D_j -> infinity.
(d) NON-RESONANCE.  Suppose m alpha + n beta = p in Z with (m,n) != (0,0).  Multiplying by k_j,
        m(gamma1 + eps_j) + n(gamma2 + eps'_j)  in  Z    for every j.
    Using gamma2 = 1 - gamma1 this is (m-n) gamma1 + m eps_j + n eps'_j in Z.  The integers so
    obtained are bounded, hence eventually constant = I, and m eps_j + n eps'_j = I -
    (m-n)gamma1 is then constant in j and tends to 0, so it is 0 and (m-n) gamma1 = I in Z.
    theta* is IRRATIONAL (Niven; R_01 Z2), so m = n.  Then m(eps_j + eps'_j) = 0 with
    eps_j, eps'_j > 0 STRICTLY, so m = 0 and n = 0.  Contradiction.  Hence L(omega) = {0}. []

BY W-12's COROLLARY 1 (carried fact C4) the map from connections to (W_F,W_C) is ONTO T^2 on
every carrier, so this omega IS a connection on K1 -- the counterexample is inside the corpus's
own object, not beside it.

WHAT IS COMPUTED BELOW, AND WHAT IS ONLY STATED.  Computed: the first two rungs, exactly, and
the running average across seven decades against three control arms.  Stated, not run: rungs
three and beyond, whose offsets are below 10^-777 and whose rung positions are past 10^780.
Sparseness is intrinsic (R_ladder): a dip of depth D at k needs the next rung at k e^{Dk}.
"""
import math
import numpy as np
from fractions import Fraction
from R_lib import PI_K1, HighPhaseReducer, m_jensen, arm_hash, diff_arms
from R_ladder import build_ladder, log_fraction

PREC = 900
GAM1 = hp = None
from R_lib import hp_theta_star
GAM1 = hp_theta_star(PREC)
GAM2 = 1 - GAM1
MP = m_jensen(PI_K1, 1 << 24)
p00, p10, p01, p11 = PI_K1
KMAX = 10 ** 7
DEC = [10 ** i for i in range(1, 8)]

print("=" * 79)
print("R_05 — THEOREM R-3: H2 HOLDS, THE LIMIT DIVERGES, AT K1's REGISTERED pi")
print("=" * 79)
print("\npi = %s ;  m(P) = %.12f ;  zeros at (theta*, 1-theta*) and (1-theta*, theta*)"
      % (PI_K1, MP))

# ------------------------------------------------------------------ the local expansion
x0 = np.exp(2j * np.pi * float(GAM1)); y0 = np.exp(2j * np.pi * float(GAM2))
A = (p10 + p11 * y0) * x0
B = (p01 + p11 * x0) * y0
print("\nTHE LOCAL EXPANSION AT THE ZERO, AND WHY THE TWO OFFSETS CANNOT CANCEL")
print("   P(x0 e(d1), y0 e(d2)) = 2 pi i [ A d1 + B d2 ] + O(d^2),  A = (p10+p11 y0) x0,")
print("   B = (p01+p11 x0) y0.   |A| = %.12f  |B| = %.12f   B - conj(A) = %.3e"
      % (abs(A), abs(B), abs(B - np.conj(A))))
psi = math.atan2(A.imag, A.real)
print("   B = conj(A) exactly (real coefficients, y0 = conj(x0)), so with d1,d2 > 0 REAL")
print("   |A d1 + B d2| = |A| sqrt((d1+d2)^2 cos^2 psi + (d1-d2)^2 sin^2 psi) >= |A||cos psi|(d1+d2)")
print("   psi = arg A = %.9f rad, cos psi = %.9f != 0  -> NO CANCELLATION." % (psi, math.cos(psi)))
print("   VALIDATION against direct evaluation, four decades of offset:")
for d in (1e-3, 1e-5, 1e-7, 1e-9):
    d1, d2 = d, 0.37 * d
    direct = abs(p00 + p10 * x0 * np.exp(2j * np.pi * d1) + p01 * y0 * np.exp(2j * np.pi * d2)
                 + p11 * x0 * y0 * np.exp(2j * np.pi * (d1 + d2)))
    loc = abs(2 * np.pi * (A * d1 + B * d2))
    print("      d = %8.1e   direct %.12e   local %.12e   rel dev %.2e"
          % (d, direct, loc, abs(direct - loc) / direct))

# ------------------------------------------------------------------ the ladders
D1, D2 = 0.5, 1.2
k1 = 10
M1 = int(math.ceil(math.exp(D1 * k1)))
k2 = k1 * M1
M2 = 1 << int(math.ceil(D2 * k2 / math.log(2.0)))
ksA, nsA, alA, alpha, epsA = build_ladder(GAM1, k1, 3, [M1, M2])
ksB, nsB, alB, beta, epsB = build_ladder(GAM2, k1, 7, [M1, M2])
assert ksA == ksB
print("\nTHE TWO LADDERS, ON A COMMON RUNG SEQUENCE")
print("   rungs: k_1 = %d, k_2 = %d, k_3 = k_2 * 2^%d  (%d digits)"
      % (ksA[0], ksA[1], M2.bit_length() - 1, len(str(ksA[2]))))
for j in (0, 1):
    la, lb = log_fraction(epsA[j]), log_fraction(epsB[j])
    bnd = log_fraction(Fraction(2 * ksA[j], ksA[j + 1]))
    print("   j=%d   log eps_j  = %11.3f   log eps'_j = %11.3f   proved bound log(2k_j/k_{j+1}) = %11.3f"
          % (j + 1, la, lb, bnd))
    assert epsA[j] > 0 and epsB[j] > 0, "POSITIVITY VIOLATED — the non-resonance proof needs it"
    assert la < bnd and lb < bnd, "LADDER BOUND VIOLATED"
print("   POSITIVITY OF ALL FOUR OFFSETS ASSERTED IN CODE (not claimed in prose).  This is the")
print("   step that makes the pair provably non-resonant: proof step (d) reduces the whole")
print("   lattice L(omega) to the single case m = n with eps_j + eps'_j = 0, and that is excluded")
print("   by a STRICT INEQUALITY the construction guarantees.")
print("   gamma1 + gamma2 = %s (exactly 1, as proof step (d) requires)" % (GAM1 + GAM2))

# ------------------------------------------------------------------ the truncation, disclosed
print("\nTHE TRUNCATION, DISCLOSED BEFORE ANY NUMBER IS REPORTED")
print("   The THEOREM is about the infinite ladder with EXACT theta*.  The COMPUTATION uses the")
print("   ladder truncated at rung 3 with theta* replaced by its %d-digit rational truncation." % PREC)
print("   Consequences, both bounded: (i) alpha differs from the infinite ladder's alpha by less")
print("   than 2/k_3 < 10^-777, so over k <= 10^7 the two orbits are indistinguishable; (ii) the")
print("   truncated (alpha,beta) is a pair of RATIONALS and therefore technically resonant --")
print("   its relation lattice is generated by vectors of size ~k_3 ~ 10^780, so by Boyd-Lawton")
print("   (M1_07, carried fact C6) its subtorus average differs from m(P) far below anything")
print("   visible here.  Neither effect touches the divergence, which is caused by TWO TERMS.")
smallest = None
for mm in range(-60, 61):
    for nn in range(-60, 61):
        if (mm, nn) == (0, 0):
            continue
        v = (mm * alpha + nn * beta) % 1
        v = min(v, 1 - v)
        if smallest is None or v < smallest[0]:
            smallest = (v, mm, nn)
print("   CHECK THAT COULD HAVE FAILED: min over 0 < |m|,|n| <= 60 of ||m alpha + n beta|| for the")
print("   TRUNCATED pair = %.6e at (m,n) = (%d,%d).  No small relation inside the window."
      % (float(smallest[0]), smallest[1], smallest[2]))

# ------------------------------------------------------------------ the arms
TWOPI = 2 * math.pi
ARMS = [
    ("BA-cubic", Fraction(2.0 ** (1 / 3)), Fraction(4.0 ** (1 / 3))),
    ("CORPUS", Fraction(-1.0 / TWOPI) % 1, Fraction(math.sqrt(2.0) / TWOPI)),
    ("LIOUVILLE", alpha, beta),
]
SPECIAL = {ksA[0]: (log_fraction(epsA[0]), log_fraction(epsB[0])),
           ksA[1]: (log_fraction(epsA[1]), log_fraction(epsB[1]))}
logA, logB = math.log(abs(A)), math.log(abs(B))


def run(a, b, special):
    ra, rb = HighPhaseReducer(a), HighPhaseReducer(b)
    out = np.empty(KMAX)
    CH = 10 ** 6
    for lo in range(0, KMAX, CH):
        k = np.arange(lo + 1, min(lo + CH, KMAX) + 1, dtype=np.int64)
        fa, fb = ra.frac(k), rb.frac(k)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
        out[lo:lo + len(k)] = np.log(np.abs(p00 + p10 * x + p01 * y + p11 * x * y))
    for kk, (le, lep) in special.items():
        # |Z| = 2 pi |A eps + B eps'| ; eps' / eps = exp(lep - le) is representable
        r = math.exp(lep - le) if lep - le > -700 else 0.0
        out[kk - 1] = math.log(2 * math.pi) + le + math.log(abs(A + B * r))
    return out


res = {n: run(a, b, SPECIAL if n == "LIOUVILLE" else {}) for n, a, b in ARMS}

print("\nARMS-DIFF GUARD (hashes the OUTPUT vectors, per W-10 N-6)")
nm = [n for n, _, _ in ARMS]
ok = True
for i in range(len(nm) - 1):
    ok &= diff_arms(nm[i], res[nm[i]][:200000], nm[i + 1], res[nm[i + 1]][:200000])
print("   all consecutive pairs differ: %s" % ok)

print("\n(1/N) sum_{k<=N} log|Z_k|   MINUS   m(P) = %.12f    SEVEN DECADES" % MP)
print("   %-11s" % "N" + "".join("%14d" % N for N in DEC))
for n in nm:
    print("   %-11s" % n + "".join("%+14.6f" % (float(np.mean(res[n][:N])) - MP) for N in DEC))

print("\nSTRADDLING THE RUNGS")
GR = [5, 9, 10, 11, 50, 500, 1489, 1490, 1491, 3000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
print("   %-11s" % "N" + "".join("%10d" % N for N in GR))
for n in nm:
    print("   %-11s" % n + "".join("%+10.4f" % (float(np.mean(res[n][:N])) - MP) for N in GR))

print("\n   log|Z_k| AT THE TWO RUNGS (the whole of the effect, two terms out of ten million)")
for kk in (ksA[0], ksA[1]):
    print("      k = %6d   log|Z_k| = %14.4f   contribution to S_k = %+9.4f"
          % (kk, res["LIOUVILLE"][kk - 1], res["LIOUVILLE"][kk - 1] / kk))
print("      predicted from the ladder:   -D_1 = %.2f at k_1,   -D_2 = %.2f at k_2" % (-D1, -D2))

print("\nTHE VERDICT, AND WHAT IT DOES AND DOES NOT SAY")
print("   H2 -- '(conj W_F, W_C) generates a dense subgroup of T^2' -- HOLDS for omega by proof")
print("   step (d), and the Birkhoff average does NOT converge to m(P).  H2 IS NOT SUFFICIENT.")
print("   This is the registrar's verified claim, now EXHIBITED at K1's registered pi rather")
print("   than at the auxiliary pi = (0,1/3,1/3,1/3) of M1_06.")
print("   IT DOES NOT SAY N1 IS FALSE.  It says N1 as registered carries a hypothesis that does")
print("   not do the work its wording implies.  R_06 supplies one that does.")
print("   A NULL READS TWO WAYS AND THIS IS NOT A NULL: the divergence is exhibited, not absent.")
print("\nDONE R_05")
