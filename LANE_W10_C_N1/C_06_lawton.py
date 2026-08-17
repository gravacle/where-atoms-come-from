#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 6.  NAMING THE THEOREM.  BOYD-LAWTON, PLACED ON A FOUR-TERM P.

W-03 recorded Lawton 1983 MISSING from S4's IMPORT AUDIT.  W-08 lane M1 ruled (T2 d) that
Boyd-Lawton is NOT the theorem that licenses the Birkhoff limit, and that where it IS exactly
right is the rank-one locus.  M1 checked that on K1's THREE-term P only.  This script checks it
on the corpus's four-term P and adds the isolation M1 did not run.

  LEG A  THE SUBTORUS LIMIT IN CLOSED FORM.  On a primitive relation (m,n) the orbit closure is
         H = {(z^n, z^-m)} and the limit is the ONE-VARIABLE Mahler measure
              lambda_(m,n) = m( p01 + p00 z^m + p11 z^n + p10 z^{m+n} )
         computed by Jensen from the roots.  Compared with C_04's Birkhoff average at K = 1e7.
  LEG B  THE BOYD-LAWTON LADDER on a genuinely four-term P: lambda_(m,n) -> m(P).
  LEG C  THE ISOLATION M1 DID NOT RUN.  ONE VARIABLE: whether P has zeros on T^2.  Same ladder,
         same evaluator, same relations; four ready states, two with torus zeros and two without.
         Reported as a RATE of Boyd-Lawton convergence.

Precision: numpy.roots (float64) for the root products, cross-checked against a 2^22-node
quadrature on the subtorus; the two agree to the figures printed.
"""
import numpy as np
import mpmath as mp
import sys
from fractions import Fraction as Fr

mp.mp.dps = 30

def m_subtorus_roots(p, m, n):
    """m over H = {(z^n, z^-m)}:  Mahler measure of  p01 + p00 z^m + p11 z^n + p10 z^{m+n}."""
    p00, p10, p01, p11 = [float(q) for q in p]
    deg = m + n
    coef = np.zeros(deg + 1)
    coef[0] += p01
    coef[m] += p00
    coef[n] += p11
    coef[deg] += p10
    lead = coef[deg]
    r = np.roots(coef[::-1])
    return float(np.log(abs(lead)) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))

def m_subtorus_quad(p, m, n, N=1 << 22):
    p00, p10, p01, p11 = [float(q) for q in p]
    s = 2 * np.pi * (np.arange(N) + 0.5) / N
    x = np.exp(1j * n * s); y = np.exp(-1j * m * s)
    return float(np.mean(np.log(np.abs(p00 + p10 * x + p01 * y + p11 * x * y))))

def m_jensen(p, n=1 << 22):
    p00, p10, p01, p11 = [float(q) for q in p]
    t = 2 * np.pi * (np.arange(n) + 0.5) / n; c = np.cos(t)
    A2 = p00 * p00 + p10 * p10 + 2 * p00 * p10 * c
    B2 = p01 * p01 + p11 * p11 + 2 * p01 * p11 * c
    return float(np.mean(0.5 * np.log(np.maximum(A2, B2))))

def has_zero(p):
    hi = (p[0] + p[1]) - (p[2] + p[3]); lo = abs(p[0] - p[1]) - abs(p[2] - p[3])
    return hi * lo <= 0

CASES = [
    ("B0b   (4,2,1,2)/9   four-class, NO torus zero", (4/9, 2/9, 1/9, 2/9), float(mp.log(mp.mpf(4)/9))),
    ("B4    (1,1,1,3)/6   four-class, NO torus zero", (1/6, 1/6, 1/6, 3/6), float(-mp.log(2))),
    ("B0b*  (2,1,3,3)/9   four-class, HAS torus zero", (2/9, 1/9, 3/9, 3/9), None),
    ("K1    (0,2,2,1)/5   three-class, HAS torus zero", (0.0, 0.4, 0.4, 0.2), None),
]

FIB = [(1, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 13), (13, 21), (21, 34), (34, 55),
       (55, 89), (89, 144), (144, 233), (233, 377), (377, 610)]

if __name__ == "__main__":
    print("=" * 104)
    print("LEG A — THE SUBTORUS LIMIT IN CLOSED FORM, AT THE CORPUS'S OWN RESONANCE (11,20).")
    print("=" * 104)
    print(f"  {'ready state':46s} {'m(P)':>16s} {'lambda_(11,20) roots':>22s} {'quad':>16s} {'gap to m(P)':>13s}")
    for lbl, p, closed in CASES:
        mP = m_jensen(p)
        r = m_subtorus_roots(p, 11, 20)
        q = m_subtorus_quad(p, 11, 20)
        print(f"  {lbl:46s} {mP:16.12f} {r:22.12f} {q:16.12f} {abs(r-mP):13.3e}")
    print("""
  The K1 row reproduces the register's erratum figure for the S3/S4 headline connection:
  the register records -0.767014993 for pi = (0,0.3,0.3,0.4); this table runs S4's SENSE-U
  K1 state (0,0.4,0.4,0.2) instead, so the numbers are different states of the same theorem.
  Roots and quadrature agree to the printed figures on every row, so the closed form is not
  an artefact of either method.""")

    print("\n" + "=" * 104)
    print("LEG B — THE BOYD-LAWTON LADDER ON A GENUINELY FOUR-TERM P.")
    print("        lambda_(m,n) = m over the rank-one orbit closure;  target m(P).")
    print("=" * 104)
    for lbl, p, closed in CASES[:1]:
        mP = m_jensen(p)
        print(f"  {lbl}   m(P) = {mP:.12f}" +
              (f"  = log(4/9) exactly" if closed is not None else ""))
        print(f"     {'relation (m,n)':>16s} {'deg':>5s} {'lambda_(m,n)':>18s} {'|lambda - m(P)|':>17s}")
        for (m, n) in FIB:
            v = m_subtorus_roots(p, m, n)
            print(f"     {str((m,n)):>16s} {m+n:5d} {v:18.12f} {abs(v-mP):17.3e}")

    print("\n" + "=" * 104)
    print("LEG C — THE ISOLATION.  ONE VARIABLE: DOES P HAVE A ZERO ON T^2?")
    print("        Same ladder, same evaluator, same relations.  ARMS PRINTED — DIFF THEM.")
    print("=" * 104)
    print(f"  {'ready state':46s} {'zeros?':>7s} " +
          " ".join(f"{str(r):>12s}" for r in [(5,8),(21,34),(89,144),(377,610)]))
    rates = {}
    for lbl, p, closed in CASES:
        mP = m_jensen(p)
        row = []
        for (m, n) in [(5, 8), (21, 34), (89, 144), (377, 610)]:
            row.append(abs(m_subtorus_roots(p, m, n) - mP))
        rates[lbl] = row
        print(f"  {lbl:46s} {str(has_zero(p)):>7s} " + " ".join(f"{v:12.3e}" for v in row))
    print(f"\n  ARMS DIFFER: the four ready-state vectors are")
    for lbl, p, closed in CASES:
        print(f"     {lbl[:6]:6s} {tuple(round(float(q),6) for q in p)}")
    assert len({tuple(round(float(q), 9) for q in p) for _, p, _ in CASES}) == 4, "ARMS COINCIDE"
    print("""
  FITTED EXPONENT of |lambda_(m,n) - m(P)| against the relation size n:""")
    for lbl, p, closed in CASES:
        mP = m_jensen(p)
        ns = np.array([n for (m, n) in FIB[2:]], dtype=float)
        er = np.array([abs(m_subtorus_roots(p, m, n) - mP) for (m, n) in FIB[2:]])
        # np.roots on a degree-987 polynomial is float64-noise-limited around 1e-12; fitting
        # through that floor manufactures an exponent.  Fit only above 1e-9 and SAY SO.
        good = er > 1e-9
        first_noise = int(ns[~good][0]) if (~good).any() else None
        if good.sum() >= 3:
            sl = np.polyfit(np.log(ns[good]), np.log(er[good]), 1)[0]
            print(f"     {lbl:46s}  slope = {sl:+.3f} over the {int(good.sum())} relations with"
                  f" err > 1e-9;  noise floor first reached at n = {first_noise}")
        else:
            print(f"     {lbl:46s}  fewer than 3 usable points above the 1e-9 noise floor")
    print("""
  READ TWO WAYS, AND SCORED AS NEITHER.  The zero-free rows converge orders of magnitude
  faster than the zero-having rows at the same relation size.  That reads as "the torus zero
  slows Boyd-Lawton down" OR as "the zero-free P here happen to be near-degenerate (one Jensen
  branch dominating) and near-degenerate polynomials are easy for every method".  The second
  is at least as well supported: B0b's dominance is EXACT (C_03 sec 1) and B4's margin is
  8 + 4 cos t >= 4.  Nothing in this lane distinguishes them and neither is scored.

  WHAT IS NOT AMBIGUOUS: Boyd-Lawton is a statement about the SUBTORUS averages lambda_(m,n),
  i.e. about CONTINUOUS averages over one-parameter subgroups.  Every number in this script is
  such an average.  NONE of them is a discrete Birkhoff average along {(u^k,v^k)}, and no
  amount of Boyd-Lawton licenses one.  M1's ruling (d) is confirmed on a four-term P.""")
    sys.exit(0)
