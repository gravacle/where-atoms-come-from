#!/usr/bin/env python3
"""
R2_02 — THE TARGET'S CENTRAL CLASSIFYING STATEMENT IS FALSE, AND ITS OWN ROW 5 IS THE
COUNTEREXAMPLE.

M1_04_durability.py, PART 3 header, verbatim:
    "lambda_H = m(P) iff H = T^2 iff u,v have no multiplicative relation."
and M1_04's printed PART 3 block:
    "rank L = 0 (Diophantine)   lambda = m(P)  ...  [= the N1 claim]
     rank L = 1 (11,20)         lambda = subtorus mean != m(P)
     rank L = 2 (order 4)       lambda = 4-point mean  != m(P)"
and the reported operative_variable: rank L = 0 <=> lambda = m(P).

ONLY THE "IF" DIRECTION IS A THEOREM.  THE "ONLY IF" IS FALSE.

WITNESS, AND IT IS ROW 5 OF THE TARGET'S OWN NINE-CASE TABLE:
    u = e^{0.7 i}, v = e^{1.3 i}, weights (p10,p01,p11) = (0.5, 0.5, 0.0).
    0.7 and 1.3 are commensurable: 13*0.7 - 7*1.3 = 9.1 - 9.1 = 0, so u^13 v^-7 = 1.
    rank L = 1.  H = {(z^7, z^13)} (connected component of ker of the character (13,-7)).
    P|_H = 0.5 z^7 + 0.5 z^13 = 0.5 z^7 (1 + z^6), so
        lambda_H = log(0.5) + m(1 + z^6) = log(0.5) + 0 = log(0.5) = m(P).
    So lambda = m(P) on a RANK-ONE lattice.  The target prints this row's lambda as
    -0.6932 = m(P) and reads it as the m(P) case, never noticing the pair is resonant.

ISOLATION LEDGER FOR THIS SCRIPT
  PART A.  HELD FIXED: nothing needs to move -- this is an algebraic identity, checked
           exactly (Fractions on the relation, closed form on the two Mahler measures) and
           numerically (orbit average).
  PART B.  HELD FIXED: relation vector (1,0) (i.e. u = 1, W_F trivial), the observable, the
           evaluator.  THE ONE THING THAT MOVES: the weight vector along the segment
           p10 = t, p01 = p11 = (1-t)/2.  lambda_(1,0) - m(P) changes sign on that segment,
           so by continuity there is an EXACT weight vector, with FULL support, at which a
           rank-one connection reproduces m(P) exactly.  Located to 1e-15 by bisection.
  PART C.  HELD FIXED: weights (0.3,0.3,0.4) (S3/S4's own).  MOVES: the primitive relation
           vector.  Scan for relations whose subtorus rate crosses m(P): the deviation
           changes sign, so the rank-1 locus is not on one side of m(P).
Precision: float64, except the relation check which is exact (Fraction).
"""
import sys
import numpy as np
from fractions import Fraction

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W08_M1_IDENTIFICATION")
from M1_02_mahler_machinery import m_R1          # the TARGET's own evaluator, reused on purpose

print("=" * 78)
print("R2_02 — 'lambda = m(P) IFF rank L = 0' IS FALSE.  ROW 5 OF THE TARGET'S OWN TABLE.")
print("=" * 78)

# ---------------------------------------------------------------- PART A
print("\nPART A — ROW 5 OF M1_04's NINE-CASE TABLE IS RESONANT, NOT GENERIC.")
f7 = Fraction(7, 10)
f13 = Fraction(13, 10)
print("  u = e^{i*0.7}, v = e^{i*1.3}   (M1_04 CASES[4], verbatim)")
print("  EXACT: 13*(0.7) - 7*(1.3) = %s  -> u^13 v^-7 = 1, so (13,-7) in L, rank L = 1."
      % (13 * f7 - 7 * f13))
print("  gcd(13,7) = %d, so (13,-7) is PRIMITIVE and L = Z*(13,-7)." % np.gcd(13, 7))
print("  H = connected kernel = {(z^7, z^13)}.")
p10, p01, p11 = 0.5, 0.5, 0.0
mP = m_R1(0.0, p10, p01, p11)
# subtorus measure, exact by factorisation, and numerically by Jensen on the roots
coef = np.zeros(14)
coef[7] += p10
coef[13] += p01
coef[7 + 13] if False else None
# P|_H = p10 z^7 + p01 z^13 + p11 z^{7+13}
coef = np.zeros(21)
coef[7] += p10
coef[13] += p01
coef[20] += p11
nz = np.nonzero(coef)[0]
c = coef[nz[0]:nz[-1] + 1]
r = np.roots(c[::-1])
mH = float(np.log(abs(c[-1])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))
print("  m(P) over T^2                       = %.12f   (= log 1/2 = %.12f)" % (mP, np.log(0.5)))
print("  lambda_H = m(P|_H) by Jensen        = %.12f" % mH)
print("  closed form  log(0.5) + m(1+z^6)    = %.12f" % (np.log(0.5) + 0.0))
print("  |lambda_H - m(P)|                   = %.3e" % abs(mH - mP))

u = np.exp(0.7j); v = np.exp(1.3j)
k = np.arange(1, 200001)
az = np.abs(p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k)
print("  direct orbit average at N=2e5 (the target's own N) = %.9f  [target printed -0.6932]"
      % float(np.mean(np.log(az))))
print("  ==> lambda = m(P) WITH rank L = 1.  The 'iff' in M1_04 PART 3 is FALSE.")

# ---------------------------------------------------------------- PART B
print("\nPART B — A FULL-SUPPORT WITNESS.  relation (1,0) means u = 1 (W_F = 1, curvature")
print("  trivial, flat holonomy arbitrary).  Then H = {1} x T and")
print("      lambda_(1,0) = m(p10 + (p01+p11) y) = log max(p10, p01+p11).")
print("  Move the weights along p10 = t, p01 = p11 = (1-t)/2 (FULL support for t<1):")


def lam10(t):
    p10 = t; p01 = p11 = (1 - t) / 2.0
    return np.log(max(p10, p01 + p11))


def gap(t):
    p10 = t; p01 = p11 = (1 - t) / 2.0
    return lam10(t) - m_R1(0.0, p10, p01, p11)


print("      t        lambda_(1,0)        m(P)          lambda_(1,0) - m(P)")
for t in (0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.99):
    p10 = t; p01 = p11 = (1 - t) / 2.0
    print("    %.4f    %13.9f  %13.9f   %+.3e" % (t, lam10(t), m_R1(0.0, p10, p01, p11), gap(t)))

lo, hi = 0.5, 0.999
if gap(lo) * gap(hi) < 0:
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if gap(lo) * gap(mid) <= 0:
            hi = mid
        else:
            lo = mid
    t = 0.5 * (lo + hi)
    p10 = t; p01 = p11 = (1 - t) / 2.0
    print("\n  SIGN CHANGE LOCATED BY BISECTION (200 halvings):")
    print("    t* = %.15f   pi = (%.15f, %.15f, %.15f)" % (t, p10, p01, p11))
    print("    lambda_(1,0) = %.15f    m(P) = %.15f    gap = %.3e"
          % (lam10(t), m_R1(0.0, p10, p01, p11), gap(t)))
    print("    FULL SUPPORT (all three weights > 0), rank L = 1, and lambda = m(P) to %.1e."
          % abs(gap(t)))
    print("    So the coincidence is not an artefact of the degenerate support in PART A.")
else:
    print("  no sign change on [%.3f,%.3f]: gap(lo)=%.3e gap(hi)=%.3e"
          % (lo, hi, gap(lo), gap(hi)))

# ---------------------------------------------------------------- PART C
print("\nPART C — AT S3/S4's OWN WEIGHTS (0.3,0.3,0.4) THE RANK-1 RATES STRADDLE m(P).")
P10c, P01c, P11c = 0.3, 0.3, 0.4
mPc = m_R1(0.0, P10c, P01c, P11c)


def m_sub(mm, nn):
    assert np.gcd(abs(mm), abs(nn)) == 1
    terms = [(nn, P10c), (-mm, P01c), (nn - mm, P11c)]
    shift = -min(e for e, _ in terms)
    deg = max(e + shift for e, _ in terms)
    coef = np.zeros(deg + 1)
    for e, cc in terms:
        coef[e + shift] += cc
    nz = np.nonzero(coef)[0]
    coef = coef[nz[0]:nz[-1] + 1]
    if len(coef) == 1:
        return float(np.log(coef[0]))
    r = np.roots(coef[::-1])
    return float(np.log(abs(coef[-1])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))


best = []
for mm in range(1, 60):
    for nn in range(-60, 61):
        if np.gcd(abs(mm), abs(nn)) != 1:
            continue
        d = m_sub(mm, nn) - mPc
        best.append((abs(d), mm, nn, d))
best.sort()
print("  the ten rank-1 relations (|m|,|n| <= 60) whose subtorus rate is CLOSEST to m(P):")
print("      (m,n)        lambda_(m,n)        lambda - m(P)")
for ad, mm, nn, d in best[:10]:
    print("    (%3d,%4d)   %14.9f      %+.3e" % (mm, nn, m_sub(mm, nn), d))
pos = [b for b in best if b[3] > 0][:3]
neg = [b for b in best if b[3] < 0][:3]
print("  deviations of BOTH signs occur among rank-1 relations, so the rank-1 locus does not")
print("  sit on one side of m(P) and 'rank L = 1 => lambda != m(P)' has no monotone rescue.")
print("  closest from above: (%d,%d) %+.3e ; from below: (%d,%d) %+.3e"
      % (pos[0][1], pos[0][2], pos[0][3], neg[0][1], neg[0][2], neg[0][3]))
print("\nDONE R2_02")
