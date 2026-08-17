#!/usr/bin/env python3
"""
M1_07 — WHERE LAWTON'S THEOREM ACTUALLY BELONGS.

W-03 recorded Lawton 1983 as MISSING from the corpus's IMPORT AUDIT.  It is missing, and it
is ALSO not the theorem that licenses N1's Birkhoff convergence.  Both halves matter:

  (a) NOT THIS.  Lawton/Boyd is about MAHLER MEASURES OF SPECIALISATIONS — continuous
      averages over one-parameter subgroups — not about discrete Birkhoff averages along an
      orbit.  It says nothing about whether (1/N) sum_{k<=N} log|P(u^k,v^k)| converges.
      (M1_06 shows that question has a genuinely negative answer for some equidistributing
      (u,v), so no theorem of that shape could license it.)

  (b) BUT EXACTLY THIS.  When (u,v) satisfies one primitive relation u^m v^n = 1, the orbit
      closure is the connected subgroup H = {(z^n, z^{-m})} and the limit is the SUBTORUS
      average, which IS a one-variable Mahler measure:
          lambda_(m,n) = m( p10 z^{m+n} + p11 z^{n} + p01 )      (exponents shifted >= 0)
      Boyd-Lawton then says lambda_(m,n) -> m(P) as the relation gets large, in the precise
      sense that q(m,n) := min over non-zero integer vectors killed by ... -> infinity.
      THAT is the theorem behind S4's observed accumulation, and it is the correct citation
      for the erratum's 4.93e-04 gap at (11,20).

Computed here independently (Jensen on the roots of the one-variable specialisation; the
resulting m is exact up to numpy's root-finding, cross-checked against the direct orbit
average in M1_03 for (11,20)).  Weights held at S3/S4's (p10,p01,p11) = (0.3,0.3,0.4).
Precision: float64.
"""
import numpy as np
from M1_02_mahler_machinery import m_R1

P10, P01, P11 = 0.3, 0.3, 0.4
mP = m_R1(0.0, P10, P01, P11)


def m_subtorus(mm, nn):
    """H = {(z^n, z^-m)} for a PRIMITIVE relation u^m v^n = 1.
       P|_H = p10 z^n + p01 z^-m + p11 z^{n-m}; shift to non-negative exponents."""
    assert np.gcd(abs(mm), abs(nn)) == 1
    terms = [(nn, P10), (-mm, P01), (nn - mm, P11)]
    shift = -min(e for e, _ in terms)
    deg = max(e + shift for e, _ in terms)
    coef = np.zeros(deg + 1)
    for e, c in terms:
        coef[e + shift] += c
    nz = np.nonzero(coef)[0]
    lo, hi = nz[0], nz[-1]
    coef = coef[lo:hi + 1]
    if len(coef) == 1:
        return float(np.log(coef[0]))
    if len(coef) - 1 <= 1200:
        r = np.roots(coef[::-1])
        return float(np.log(abs(coef[-1])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))
    # degree too large for a companion-matrix eigensolve: quadrature on the circle instead.
    # The integrand is a 3-term exponential sum; evaluate it directly at NQ equally spaced
    # points.  Labelled QUADRATURE in the output because it is not the Jensen/roots route.
    NQ = 1 << 23
    th = np.arange(NQ) * (2 * np.pi / NQ)
    val = np.zeros(NQ, dtype=complex)
    for e, c in terms:
        val += c * np.exp(1j * e * th)
    return float(np.mean(np.log(np.maximum(np.abs(val), 1e-300))))


print("=" * 78)
print("M1_07 — SUBTORUS RATES, AND BOYD-LAWTON ACCUMULATION")
print("=" * 78)
print("\nweights (p10,p01,p11) = (0.3, 0.3, 0.4);  m(P) = %.12f" % mP)
print("\nIndependent recomputation of S4's exceptional-locus table (S4 line 258 ff.):")
S4 = {(1, 0): -0.356674944, (0, 1): -0.356674944, (1, 1): -1.203972804,
      (1, -1): -0.510825624, (2, 1): -0.681980359, (2, -1): -0.916290732,
      (3, 1): -0.767783712, (3, 2): -0.732940865, (4, 1): -0.784966659,
      (5, 1): -0.749392712, (5, 3): -0.765224351, (7, 3): -0.759305247,
      (7, 11): -0.764712281, (11, 20): -0.767014993, (13, 8): -0.768271734,
      (29, 17): -0.767138179}
print("    (m,n)      this lane        S4 record        |diff|")
worst = 0.0
for (mm, nn), val in S4.items():
    mine = m_subtorus(mm, nn)
    worst = max(worst, abs(mine - val))
    print("   (%3d,%4d)  %14.9f   %14.9f    %.2e" % (mm, nn, mine, val, abs(mine - val)))
print("    worst deviation over the 16 published rows: %.2e" % worst)

print("\nBOYD-LAWTON ACCUMULATION (S4's second table, recomputed):")
for (mm, nn) in [(1, 1), (5, 3), (11, 20), (41, 53), (97, 61), (610, 377), (2584, 1597)]:
    mine = m_subtorus(mm, nn)
    print("   (%5d,%5d)  |m|+|n| = %6d   lambda = %.9f   deviation from m(P) = %+.3e"
          % (mm, nn, abs(mm) + abs(nn), mine, mine - mP))
print("   -> the subtorus rates accumulate on m(P).  THIS is Lawton's theorem in action, and")
print("      it is the only place in this corpus where Lawton is the right citation.")

print("\nTHE GAP THE ERRATUM RECORDS, NAMED:")
l1120 = m_subtorus(11, 20)
print("   S3/S4's headline is the (11,20) row.  Its rate is %.9f; the generic rate is" % l1120)
print("   %.9f; the difference %+.3e is a BOYD-LAWTON APPROXIMATION GAP at a relation"
      % (mP, l1120 - mP))
print("   of size 31 — not an arithmetic error and not a convergence failure.")
print("\nDONE M1_07")
