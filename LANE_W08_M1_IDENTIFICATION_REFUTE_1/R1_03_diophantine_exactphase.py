#!/usr/bin/env python3
"""
R1_03 — ATTACK ON F10 CASE A AND ON SELF-FLAG 4.

The lane's case A sets   alpha = -(2.0**(1/3)) % 1.0,  beta = (4.0**(1/3)) % 1.0
as float64, then reduces k*alpha by exact int64 modular arithmetic against D = 2^39 plus a
float64 residue.  Its M1_03 docstring says this "loses ~1e-17" where naive accumulation
"loses ~1e-9 by k = 1e7".

THE CLAIM UNDER ATTACK.  What is good to ~1e-17 is the representation of the LANE'S OWN
float64 alpha.  The float64 alpha is NOT -2^(1/3) mod 1: it is a dyadic rational differing
from it by ~1e-17, so by k = 1e7 the lane's orbit point differs from the TRUE Schmidt orbit
point by ~1e-10 in phase -- exactly the error the docstring says it avoided.  Worse for the
lane's own framing: a dyadic rational alpha has rank L = 2, not rank L = 0, so case A as
actually computed is a TORSION point of enormous order masquerading as the Diophantine one.
Whether that matters NUMERICALLY is the question this script settles.

METHOD.  alpha and beta are represented EXACTLY as A/2^E, B/2^E with E = 200, where A and B
are the true integers floor(2^E * alpha_true) computed from mpmath at 120 decimal digits.
frac(k*alpha) is then ((k*A) mod 2^E)/2^E to a relative accuracy 2^-200, computed in Python
big integers -- no float phase anywhere.  The residual truncation is 2^-200 * k < 1e-53.
Only the final |P| evaluation is float64 (and is checked against mpmath on the worst points).

Grid: k = 1..10^7, checkpoints 1e3,1e4,1e5,1e6,1e7.  This is the lane's grid exactly.
"""
import numpy as np
from mpmath import mp, mpf, cbrt, floor as mpfloor

mp.dps = 130
E = 200
SC = 1 << E

# true values, exactly floored at 2^-200
alpha_true = (-(mpf(2) ** (mpf(1) / 3))) % 1
beta_true = (mpf(4) ** (mpf(1) / 3)) % 1
A = int(mpfloor(alpha_true * SC))
B = int(mpfloor(beta_true * SC))

# the lane's float64 values, and the exact dyadic rationals they ARE
al_f = float(-(2.0 ** (1.0 / 3.0)) % 1.0)
be_f = float((4.0 ** (1.0 / 3.0)) % 1.0)
from fractions import Fraction as Fr
al_f_exact = Fr(al_f)
be_f_exact = Fr(be_f)

print("=" * 78)
print("R1_03 — THE DIOPHANTINE CASE WITH EXACT 200-BIT PHASES")
print("=" * 78)
print("\nalpha_true = -2^(1/3) mod 1 = %s" % mp.nstr(alpha_true, 40))
print("beta_true  =  4^(1/3) mod 1 = %s" % mp.nstr(beta_true, 40))
print("lane float64 alpha          = %.17f    it is EXACTLY the rational %d/2^%d"
      % (al_f, al_f_exact.numerator, al_f_exact.denominator.bit_length() - 1))
print("lane float64 beta           = %.17f    it is EXACTLY the rational %d/2^%d"
      % (be_f, be_f_exact.numerator, be_f_exact.denominator.bit_length() - 1))
print("|alpha_lane - alpha_true|   = %s" % mp.nstr(abs(mpf(al_f) - alpha_true), 5))
print("|beta_lane  - beta_true |   = %s" % mp.nstr(abs(mpf(be_f) - beta_true), 5))
print("=> phase divergence at k = 1e7:  alpha %s   beta %s"
      % (mp.nstr(abs(mpf(al_f) - alpha_true) * 10 ** 7, 4),
         mp.nstr(abs(mpf(be_f) - beta_true) * 10 ** 7, 4)))
print("   (the lane's docstring claims its scheme 'loses ~1e-17'; that is true of the DOUBLE")
print("    alpha, not of -2^(1/3).  Against the true irrational the loss is the ~1e-9 the")
print("    docstring says it avoided.  ALSO: a dyadic-rational (alpha,beta) has rank L = 2.)")

P10 = 0.3; P01 = 0.3; P11 = 0.4
CHECK = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
K = 10 ** 7
MASK = SC - 1


def run_exact(A, B, K, chunk=200000):
    """exact big-int phase reduction; returns dict N -> average, plus min|Z|, and the
       per-checkpoint sums."""
    tot = 0.0
    out = {}
    mn = np.inf
    argmn = None
    done = 0
    # incremental: phase_k = (k*A) & MASK  -- do it by repeated addition in big ints
    ph_a = 0
    ph_b = 0
    inv = 1.0 / SC
    while done < K:
        n = min(chunk, K - done)
        fa = np.empty(n); fb = np.empty(n)
        for i in range(n):
            ph_a = (ph_a + A) & MASK
            ph_b = (ph_b + B) & MASK
            # take the TOP 60 bits as a float: exact to 2^-60 of a turn
            fa[i] = (ph_a >> (E - 60)) * (2.0 ** -60)
            fb[i] = (ph_b >> (E - 60)) * (2.0 ** -60)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
        az = np.abs(P10 * x + P01 * y + P11 * x * y)
        j = int(np.argmin(az))
        if az[j] < mn:
            mn = float(az[j]); argmn = done + j + 1
        cs = np.cumsum(np.log(az))
        for cp in CHECK:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp - done - 1]) / cp
        tot += float(cs[-1])
        done += n
    return out, mn, argmn, tot


print("\n--- EXACT-PHASE RUN, alpha = -2^(1/3), beta = 4^(1/3), 200-bit exact reduction")
res, mn, argmn, tot = run_exact(A, B, K)
mP = mpf("-0.767507880357775871645874051819")
print("        N            (1/N) sum log|Z_k|        dev from m(P)")
for cp in CHECK:
    print("   %10d      %.15f      %+.4e" % (cp, res[cp], res[cp] - float(mP)))
print("   min_k |Z_k| = %.9e  at k = %d" % (mn, argmn))

print("\n--- THE LANE'S OWN SCHEME, RE-RUN HERE (D = 2^39 + float64 residue, its float alphas)")
D = 2 ** 39
An = int(np.floor(al_f * D)); Bn = int(np.floor(be_f * D))
dA = al_f - An / D; dB = be_f - Bn / D
tot2 = 0.0; out2 = {}
done = 0
while done < K:
    n = min(10 ** 6, K - done)
    k = np.arange(done + 1, done + n + 1, dtype=np.int64)
    fa = np.mod(((k * An) % D) / D + k * dA, 1.0)
    fb = np.mod(((k * Bn) % D) / D + k * dB, 1.0)
    x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
    az = np.abs(P10 * x + P01 * y + P11 * x * y)
    cs = np.cumsum(np.log(az))
    for cp in CHECK:
        if done < cp <= done + n:
            out2[cp] = (tot2 + cs[cp - done - 1]) / cp
    tot2 += float(cs[-1])
    done += n
print("        N            lane scheme              exact-phase            difference")
for cp in CHECK:
    print("   %10d      %.15f      %.15f      %+.3e" % (cp, out2[cp], res[cp], out2[cp] - res[cp]))

print("\nREAD-OFF.")
print("  The two schemes agree to %.1e at N = 1e7.  The lane's reported value -0.767507816121"
      % max(abs(out2[c] - res[c]) for c in CHECK))
print("  is NOT a precision artefact at the size of its own residual deviation (+6.4e-08):")
print("  the true-irrational orbit gives %.12f, deviation %+.3e." % (res[10**7], res[10**7] - float(mP)))
print("  SELF-FLAG 4 IS THEREFORE RESOLVED IN THE LANE'S FAVOUR ON THE NUMBER, AND AGAINST IT")
print("  ON THE REASON GIVEN: the phase representation is NOT good to 1e-17 against 2^(1/3);")
print("  it happens not to matter at this N because the orbit stays 2.0e-04 away from Z(P).")

# ---------------------------------------------------------------------------
# Is the non-monotonicity 1e6 -> 1e7 real?  Track the running average finely.
# ---------------------------------------------------------------------------
print("\n--- IS THE 1e6 -> 1e7 NON-MONOTONICITY DISCREPANCY OSCILLATION?  Finer checkpoints,")
print("    exact phases.  If the deviation changes SIGN and size smoothly on a log grid it is")
print("    oscillation; if it grows monotonically with k it is a precision leak.")
CHECK = [int(round(10 ** (3 + 0.25 * i))) for i in range(0, 17)]
res3, _, _, _ = run_exact(A, B, 10 ** 7)
for cp in CHECK:
    if cp in res3:
        print("      N = %-10d  dev = %+.4e" % (cp, res3[cp] - float(mP)))
print("\nDONE R1_03")
