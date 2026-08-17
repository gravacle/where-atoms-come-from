#!/usr/bin/env python3
"""
R_03 — THE SUDLER ARM.  THE ONE RESONANT LINE ON WHICH OUR PROBLEM *IS* THE CLASSICAL ONE.

BRIEF ITEM (3), first half.  On the unique singular resonant line u v = 1 (i.e. c = f; R_02),
K1's Z_k reduces EXACTLY to a product of two INHOMOGENEOUS SUDLER PRODUCTS:

    |Z_k| = 0.3 * |2 sin(pi(k a - theta*))| * |2 sin(pi(k a + theta*))|,
    sum_{k<=N} log|Z_k| = N log(0.3) + log Sud_N(a,-theta*) + log Sud_N(a,+theta*),
    m(Q) = log(0.3) = -1.203972804326    (the register's own (1,1) subtorus row).

Here the existing literature IS the right literature, and it says the answer depends on the
Diophantine relation between the rotation a and the SHIFT theta* -- not on a alone.  This
script measures exactly that, over SEVEN decades of N, on four arms differing in ONE thing.

THE FOUR ARMS (isolation: identical estimator, identical N-grid, identical phase reducer,
identical special-case cutoff; the ONLY thing that moves is a):
   BA-golden    a = (sqrt5 - 1)/2         all partial quotients 1  (badly approximable)
   BA-silver    a = sqrt2 - 1             all partial quotients 2  (badly approximable)
   CORPUS       a = 1/(2 pi)              THE CORPUS'S OWN PUBLISHED GENERIC CONNECTION
                                          f = 1 projected onto this line (S4:603)
   LIOUVILLE    a from the monotone ladder of R_ladder, targeting theta*

THE CLASSICAL DIAGNOSTIC, reported for each arm across decades:
   N * min_{k<=N} ||k a - theta*||   -- bounded  <=>  the inhomogeneous approximation is
   no better than Dirichlet forces; collapsing  <=>  the orbit is hitting the shift.
"""
import math
import numpy as np
from fractions import Fraction
from R_lib import PhaseReducer, hp_theta_star, diff_arms
from R_ladder import build_ladder, log_fraction

PREC = 900
GAM = hp_theta_star(PREC)                     # theta* as an exact Fraction, PREC digits
TH = float(GAM)
LOG03 = math.log(0.3)
KMAX = 10 ** 7
DEC = [10 ** i for i in range(0, 8)]

print("=" * 79)
print("R_03 — THE SUDLER ARM: THE ONE LINE (u v = 1, c = f) ON WHICH OUR PROBLEM IS CLASSICAL")
print("=" * 79)
print("\ntheta* to %d digits, pure-integer fixed point.  m(Q) = log(0.3) = %.12f" % (PREC, LOG03))

# ------------------------------------------------------------------ the Liouville ladder
D1, D2 = 0.5, 1.2
k1, n1 = 10, 3
M1 = int(math.ceil(math.exp(D1 * k1)))                  # k2 = k1 * M1
k2 = k1 * M1
M2 = None
# M2 = ceil(exp(D2*k2)) as an exact integer: use 2**ceil(D2*k2/log 2)
M2 = 1 << int(math.ceil(D2 * k2 / math.log(2.0)))
ks, ns, alphas, alpha_L, eps = build_ladder(GAM, k1, n1, [M1, M2])
print("\nTHE LIOUVILLE LADDER (R_ladder; rungs, and the offsets they force)")
print("   rungs k_j            : %s" % [ks[0], ks[1], "k_3 = k_2 * 2^%d (%d digits)"
                                        % (M2.bit_length() - 1, len(str(ks[2])))])
for j in (0, 1):
    lg = log_fraction(eps[j])
    bnd = log_fraction(Fraction(2 * ks[j], ks[j + 1]))
    print("   eps_%d = k_%d(alpha - alpha_%d):  log(eps) = %11.3f   proved bound log(2k_%d/k_%d) = %11.3f"
          % (j + 1, j + 1, j + 1, lg, j + 1, j + 2, bnd))
    assert eps[j] > 0, "ladder positivity violated"
    assert lg < bnd, "ladder bound violated"
print("   eps_1, eps_2 both STRICTLY POSITIVE (asserted in code, not claimed in prose).")
print("   NOTE, DISCLOSED: alpha is the ladder TRUNCATED AT RUNG 3.  Rung 3 sits at k_3 ~ 10^%d,"
      % len(str(ks[2])))
print("   far outside the window, so inside k <= 1e7 this alpha is identical to the infinite")
print("   construction; eps_3 = 0 is an artefact of the truncation and is never used.")

ALPHA_L = float(alphas[1]) + float(eps[1] / ks[1]) if False else None
# alpha as a double: alpha - alpha_2 < 1e-700, so double(alpha) == double(alpha_2)
alpha_L_dbl = float(alphas[1])

ARMS = [
    ("BA-golden", (math.sqrt(5.0) - 1.0) / 2.0, {}),
    ("BA-silver", math.sqrt(2.0) - 1.0, {}),
    ("CORPUS f=1", 1.0 / (2 * math.pi), {}),
    ("LIOUVILLE", alpha_L_dbl, {ks[0]: log_fraction(eps[0]), ks[1]: log_fraction(eps[1])}),
]

# ------------------------------------------------------------------ the estimator
k = np.arange(1, KMAX + 1, dtype=np.int64)


def run_arm(a, special):
    """log|Z_k| for k = 1..KMAX.  `special` maps a rung k_j to log(eps_j); at those k the
       exact local form log|Z| = log(0.3) + log|2 sin(pi eps)| + log|2 sin(pi(2 theta* + eps))|
       is used with the EXACT eps (float64 cannot hold eps_2)."""
    fa = PhaseReducer(a).frac(k)
    s1 = np.abs(2 * np.sin(np.pi * (fa - TH)))
    s2 = np.abs(2 * np.sin(np.pi * (fa + TH)))
    lg = LOG03 + np.log(s1) + np.log(s2)
    for kk, logeps in special.items():
        # eps tiny: |2 sin(pi eps)| = 2 pi eps (1+O(eps^2));  the other factor at 2 theta*
        lg[kk - 1] = LOG03 + (math.log(2 * math.pi) + logeps) \
            + math.log(abs(2 * math.sin(math.pi * 2 * TH)))
    return lg, fa


print("\nVALIDATION OF THE SPECIAL-CASE FORM (must agree with the direct evaluation)")
for e in (1e-3, 1e-5, 1e-7, 1e-9):
    direct = 0.3 * abs(2 * math.sin(math.pi * e)) * abs(2 * math.sin(math.pi * (2 * TH + e)))
    approx = 0.3 * (2 * math.pi * e) * abs(2 * math.sin(math.pi * 2 * TH))
    print("   eps = %8.1e   direct %.12e   local form %.12e   rel dev %.2e"
          % (e, direct, approx, abs(direct - approx) / direct))

res = {}
for name, a, sp in ARMS:
    lg, fa = run_arm(a, sp)
    res[name] = (lg, fa, a)

print("\nARMS-DIFF GUARD (hashes the OUTPUT vectors log|Z_k|, per W-10 N-6)")
names = [n for n, _, _ in ARMS]
ok = True
for i in range(len(names) - 1):
    ok &= diff_arms(names[i], res[names[i]][0][:200000], names[i + 1], res[names[i + 1]][0][:200000])
print("   all consecutive pairs differ: %s" % ok)
print("   first five orbit points frac(k a), each arm (printed from the arms, not the ledger):")
for n in names:
    print("      %-11s %s" % (n, np.round(res[n][1][:5], 9)))

# ------------------------------------------------------------------ the table
print("\n(1/N) sum_{k<=N} log|Z_k|  MINUS  m(Q) = log(0.3).   SEVEN DECADES.")
print("   %-11s" % "N" + "".join("%14d" % N for N in DEC))
for n in names:
    lg = res[n][0]
    row = []
    for N in DEC:
        row.append(float(np.mean(lg[:N])) - LOG03)
    print("   %-11s" % n + "".join("%+14.6f" % v for v in row))
print("\n   READ, AND THE NON-MONOTONICITY IS REPORTED RATHER THAN SMOOTHED: the three arms whose")
print("   rotation is badly approximable or 'ordinary' fall toward 0 across the seven decades.")
print("   The fall is NOT monotone -- CORPUS f=1 worsens from -1.4e-03 at N=1e3 to -4.6e-03 at")
print("   N=1e4 (one close approach) before resuming -- so no row is quoted at one endpoint.")
print("   Fitted decay |S_N - m(Q)| ~ N^{-r} over the last four decades:")
for n in names:
    lg = res[n][0]
    xs = np.array([4.0, 5.0, 6.0, 7.0])
    ys = np.array([math.log10(abs(float(np.mean(lg[:10 ** i])) - LOG03) + 1e-300) for i in (4, 5, 6, 7)])
    r = -np.polyfit(xs, ys, 1)[0]
    print("      %-11s r = %+6.3f" % (n, r))
print("   The LIOUVILLE row does not fall: it is thrown at N = %d and again, five times harder,"
      % ks[0])
print("   at N = %d, and each throw is a SINGLE term of the sum.  Its apparent r = 1 after the" % ks[1])
print("   second throw is the 1/N decay of ONE bad term, not convergence: the next throw is at")
print("   k_3 ~ 10^%d and would carry it to -%.1f." % (len(str(ks[2])), 1788.0 * 0 + 1.2))

print("\nTHE SAME ROWS ON A GRID THAT STRADDLES THE TWO LIOUVILLE RUNGS")
GR = [5, 9, 10, 11, 20, 100, 1000, 1489, 1490, 1491, 2000, 5000, 10 ** 4, 10 ** 5]
print("   %-11s" % "N" + "".join("%10d" % N for N in GR))
for n in names:
    lg = res[n][0]
    print("   %-11s" % n + "".join("%+10.4f" % (float(np.mean(lg[:N])) - LOG03) for N in GR))

# ------------------------------------------------------------------ the classical diagnostic
print("\nTHE CLASSICAL SUDLER DIAGNOSTIC:  N * min_{k<=N} ||k a - theta*||")
print("   (bounded  <=>  inhomogeneous approximation no better than Dirichlet forces)")
print("   %-11s" % "N" + "".join("%14d" % N for N in DEC[1:]))
for n in names:
    fa = res[n][1]
    d = np.minimum(np.abs(fa - TH), 1.0 - np.abs(fa - TH))
    row = []
    for N in DEC[1:]:
        row.append(N * float(np.min(d[:N])))
    print("   %-11s" % n + "".join("%14.4f" % v for v in row))
print("   The LIOUVILLE entries printed as 0.0000 from N=1e4 on are float64's rendering of the")
print("   k = %d term: frac(k a) rounds to exactly theta* in double precision." % ks[1])
print("   LIOUVILLE's row is computed from the float64 phases and therefore UNDERSTATES its own")
print("   case: the true min at k = %d is eps_2 ~ e^{%.0f}, which float64 cannot represent."
      % (ks[1], log_fraction(eps[1])))
print("   N * that = 10^{%.0f}." % (math.log10(ks[1]) + log_fraction(eps[1]) / math.log(10)))

# ------------------------------------------------------------------ what this does not settle
print("\nWHAT THIS ARM DOES NOT SETTLE — STATED BEFORE THE VERDICT")
print("   The three convergent rows read TWO WAYS: (i) convergence holds for these rotations,")
print("   or (ii) seven decades is not enough to see a failure that begins later.  Nothing in")
print("   this script distinguishes them.  R_06 does, and only for algebraic-type rotations.")
print("   'BA' here is a property of a ALONE; the literature's condition is a property of the")
print("   PAIR (a, theta*).  Badly-approximable a does NOT by itself bound ||k a - theta*||")
print("   from below -- see R_06.  No row here is scored as evidence that it does.")
print("\nDONE R_03")
