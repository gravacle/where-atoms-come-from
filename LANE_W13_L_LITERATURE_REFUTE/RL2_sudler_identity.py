#!/usr/bin/env python3
"""
RL2 — THE TARGET'S FINDING L-11 CARRIES A FALSE IDENTITY, AND THE CORRECT ONE IS A RATIO.

TARGET CLAIM (LANE_W13_L_LITERATURE/THEOREMS_AND_CITATIONS.txt Sec. 3.3, and finding L-11):

    "Take P = (1 + x)/2 ... Then |P(u^k)| = |cos(pi k alpha)| and
         SUM_{k<=N} log|Z_k| = log PROD_{k<=N} |2 sin(pi k alpha)| - N log 2 + O(1),
     i.e. N1 in one variable IS THE SUDLER PRODUCT."

The first half is right.  The second half is FALSE, and it is false for a reason that is
visible in the definition the target itself cites: Lubinsky's object, and the Sudler product
of Aistleitner-Technau-Zafeiropoulos, is

    P_N(alpha) = PROD_{r=1}^{N} |2 sin(pi r alpha)| = PROD_{r=1}^{N} |1 - q^r|,  q = e^{2 pi i alpha}

-- a product of (1 - q^r).  The corpus's two-class object is a product of (1 + q^r), because
the coefficients of P are PROBABILITIES and cannot carry a minus sign.  And

    PROD (1 + q^r) = PROD (1 - q^{2r}) / PROD (1 - q^r) = P_N(2 alpha) / P_N(alpha).

EXACT IDENTITY (proved, not fitted):
    SUM_{k<=N} log|Z_k| = log P_N(2 alpha) - log P_N(alpha) - N log 2.

So the one-variable shadow of N1 is a DIFFERENCE OF TWO SUDLER SUMS AT alpha AND 2 alpha,
not one Sudler sum "up to N log 2".  The remainder the target calls O(1) is
log P_N(2 alpha) - 2 log P_N(alpha), and this leg shows over five decades that it is not
bounded.

ONE VARIABLE MOVES: alpha (three arithmetic types).  Everything else -- N schedule, summation
order, estimator, code path -- is identical, built by the same function.
"""
import math
import hashlib
import numpy as np
from mpmath import mp, mpf, sqrt as msqrt, pi as mpi

mp.dps = 60
DEC = [10**2, 10**3, 10**4, 10**5, 10**6]

print("=" * 78)
print("RL2 — 'N1 AT d=1 IS THE SUDLER PRODUCT' IS FALSE AS AN IDENTITY.  IT IS A RATIO.")
print("=" * 78)


def split3(x_mp):
    def chop26(v):
        if v == 0.0:
            return 0.0
        e = math.frexp(v)[1]
        q = 2.0 ** (e - 26)
        return math.floor(v / q) * q
    a1 = chop26(float(x_mp))
    r = x_mp - mpf(a1)
    a2 = chop26(float(r))
    a3 = float(r - mpf(a2))
    return a1, a2, a3


def turns(k, s):
    a1, a2, a3 = s
    return (np.mod(k * a1, 1.0) + np.mod(k * a2, 1.0) + k * a3) % 1.0


def sums(alpha_mp, N):
    """Return (S_Z, S_sud(alpha), S_sud(2 alpha)) as running sums to N, one pass, float64.
       S_Z    = sum log|cos(pi k alpha)|            (= sum log|Z_k| for P=(1+x)/2)
       S_sud  = sum log|2 sin(pi k alpha)|          (log of the SUDLER product)
    """
    sa = split3(alpha_mp % 1)
    sb = split3((2 * alpha_mp) % 1)
    SZ = 0.0
    SS1 = 0.0
    SS2 = 0.0
    CH = 10 ** 6
    done = 0
    while done < N:
        n = min(CH, N - done)
        k = np.arange(done + 1, done + n + 1, dtype=np.float64)
        ta = turns(k, sa)
        tb = turns(k, sb)
        SZ += float(np.sum(np.log(np.abs(np.cos(np.pi * ta)))))
        SS1 += float(np.sum(np.log(np.abs(2 * np.sin(np.pi * ta)))))
        SS2 += float(np.sum(np.log(np.abs(2 * np.sin(np.pi * tb)))))
        done += n
    return SZ, SS1, SS2


ARMS = {
    "sqrt2-1  (bounded p.q. [2,2,2,...])": msqrt(mpf(2)) - 1,
    "golden   (p.q. all 1)              ": (msqrt(mpf(5)) - 1) / 2,
    "1/(2pi)  (generic)                 ": 1 / (2 * mpi),
}

print("""
IDENTITY UNDER TEST, EXACT:   sum log|Z_k|  =  S_sud(2a) - S_sud(a) - N log 2
TARGET'S CLAIM:               sum log|Z_k|  =  S_sud(a)              - N log 2 + O(1)
so the target's remainder is  R_N = S_sud(2a) - 2 S_sud(a),  claimed O(1).
""")

tables = {}
for name, a in ARMS.items():
    print("-" * 78)
    print("ARM alpha = %s   (%s)" % (mp.nstr(a % 1, 25), name.strip()))
    print("   %-10s %-22s %-22s %-16s %s" % ("N", "sum log|Z_k|", "exact identity RHS", "|residual|", "R_N (claimed O(1))"))
    rows = []
    for N in DEC:
        SZ, S1, S2 = sums(a, N)
        rhs = S2 - S1 - N * math.log(2.0)
        R = S2 - 2 * S1
        rows.append((N, SZ, rhs, abs(SZ - rhs), R))
        print("   %-10d %-22.10f %-22.10f %-16.3e %+.6f" % (N, SZ, rhs, abs(SZ - rhs), R))
    tables[name] = rows
    Rs = [r[4] for r in rows]
    grow = abs(Rs[-1]) / max(abs(Rs[0]), 1e-12)
    print("   R_N over five decades: %s" % "  ".join("%+.3f" % r for r in Rs))
    print("   |R_N| grew by a factor %.1f from N=1e2 to N=1e6." % grow)

print("""
--------------------------------------------------------------------------------
ARM DIFF ON OUTPUTS (not only on inputs).""")
hs = {}
for name, rows in tables.items():
    hs[name] = hashlib.sha256(repr([("%.12f" % c if isinstance(c, float) else c) for r in rows for c in r]).encode()).hexdigest()
    print("   sha256 %s  %s" % (hs[name][:24], name.strip()))
print("   ALL THREE TABLES DISTINCT: %s" % (len(set(hs.values())) == 3))
assert len(set(hs.values())) == 3

print("""
--------------------------------------------------------------------------------
READING.

1. THE EXACT IDENTITY HOLDS TO ROUND-OFF IN EVERY ROW.  sum log|Z_k| is
   S_sud(2 alpha) - S_sud(alpha) - N log 2, EXACTLY, for every alpha.  This is elementary:
   (1+q^r) = (1-q^{2r})/(1-q^r).

2. THE TARGET'S REMAINDER R_N = S_sud(2 alpha) - 2 S_sud(alpha) IS NOT O(1).  It is a
   difference of two Sudler sums whose own literature (Lubinsky 1999; Aistleitner-Technau-
   Zafeiropoulos, Amer. J. Math. 145 (2023) 721-764) is precisely about their unbounded
   fluctuation.  Reported above across five decades as a TREND, not an endpoint.

3. WHAT SURVIVES OF L-11, AND IT IS MOST OF IT.  The one-variable shadow of N1 does live
   inside the Sudler-product literature: it is a DIFFERENCE of two Sudler sums, both terms
   are that literature's object, and the dependence on the Diophantine type of alpha is that
   literature's central fact.  L-11's CONCLUSION stands.  Its stated IDENTITY does not, and
   a note that printed it would be printing a false formula in the paragraph where it
   concedes priority.

4. AND THE d=1 OBJECT IS NAMED IN THAT LITERATURE TOO: PROD |2 cos(pi k alpha)| is the
   "cosine product" / (-q;q)_n, treated by exactly the same authors.  L-11 cites the
   (q;q)_n papers for a (-q;q)_n object.
--------------------------------------------------------------------------------""")

print("=" * 78)
print("RL2(B) — R_N IS NOT O(1): AN EXPLICIT alpha, EXACT PHASES, NO FLOAT REDUCTION.")
print("=" * 78)
print("""
The numerics above show R_N wandering; wandering over five decades is not a proof of
unboundedness and is NOT scored as one.  Here is the proof, by exhibit, using the corpus's
own device (M1_06's digit-block construction, sealed):

    alpha = 1/20 + 10^-301 (+ an irrational tail below 10^-(10^301), immaterial to k <= 20).

Then 10*alpha = 1/2 + 10^-300, so |cos(pi*10*alpha)| ~ pi*10^-300 -- a dip of depth ~689 in
sum log|Z_k| at N = 10 -- while |2 sin(pi*10*alpha)| ~ 2, i.e. the SUDLER sum at alpha has
NO dip there at all.  The dip is carried entirely by the S_sud(2 alpha) term of the exact
identity.  Hence S_Z(10) - (S_sud(alpha,10) - 10 log 2) is about -689, not O(1).
""")
from fractions import Fraction as FR
import math as _m

alpha_f = FR(1, 20) + FR(1, 10 ** 301)


def frac(x):
    return x - (x.numerator // x.denominator)


def log_abs_cos_pi(t):
    """|cos(pi t)| for t a Fraction in [0,1), exact reduction then one float sin."""
    d = t - FR(1, 2)
    d = d - round(d)            # nearest-integer reduction, exact in Fraction
    return _m.log(abs(_m.sin(_m.pi * float(d)))) if d != 0 else float("-inf")


def log_abs_2sin_pi(t):
    d = t - round(t)
    return _m.log(2 * abs(_m.sin(_m.pi * float(d)))) if d != 0 else float("-inf")


SZ = 0.0
S1 = 0.0
S2 = 0.0
print("   %-4s %-26s %-16s %-16s" % ("k", "frac(k a) - 1/2", "log|Z_k|", "log|2 sin(pi k a)|"))
for k in range(1, 21):
    t = frac(k * alpha_f)
    t2 = frac(2 * k * alpha_f)
    lz = log_abs_cos_pi(t)
    l1 = log_abs_2sin_pi(t)
    l2 = log_abs_2sin_pi(t2)
    SZ += lz
    S1 += l1
    S2 += l2
    dd = t - FR(1, 2)
    dd = dd - round(dd)
    if k in (9, 10, 11, 20):
        print("   %-4d %-26s %-16.4f %-16.4f" % (k, "%.3e" % float(dd), lz, l1))
    if k in (10, 20):
        rhs = S2 - S1 - k * _m.log(2.0)
        print("      N=%-3d  sum log|Z_k| = %-16.4f  exact identity RHS = %-16.4f  |resid| = %.2e"
              % (k, SZ, rhs, abs(SZ - rhs)))
        print("      N=%-3d  TARGET'S claimed RHS = S_sud(a) - N log2 = %-16.4f   R_N = %.4f"
              % (k, S1 - k * _m.log(2.0), S2 - 2 * S1))
print("""
   READING.  At N = 10 the target's claimed right-hand side misses the left-hand side by
   about 689 nats.  R_N = S_sud(2a) - 2 S_sud(a) is therefore NOT O(1), and the constant
   cannot be absorbed: the same construction with a deeper block makes it arbitrarily large.
   THE IDENTITY IN L-11's EVIDENCE FIELD IS FALSE.  L-11's conclusion -- that the d=1 case
   of N1 sits inside the Sudler-product literature and that the Diophantine-type dependence
   is that literature's central fact -- SURVIVES, because the correct identity is a
   difference of two Sudler sums and both terms are that literature's object.
""")
print("DONE RL2")
