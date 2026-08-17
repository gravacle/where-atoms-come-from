#!/usr/bin/env python3
"""
X_05 — THE CORRECTED CRITERION, AND THE SECOND BOYD-LAWTON ACCUMULATION THE CORPUS DOES NOT HAVE.

CORRECTED STATEMENT OF RECORD (replacing R-5's 'the one singular line is c = f'):
  Let X_omega = closure{k omega} = L(omega)^perp.  At K1's pi, Z(P) = {(e(th*),e(-th*)),
  (e(-th*),e(th*))}, both on the subgroup {xy = 1}, and th* is irrational.
    X_omega meets Z(P)   <=>   every (m,n) in L(omega) has m = n
                         <=>   L(omega) subset Z.(1,1)
                         <=>   rank L = 0, or L = Z.(d,d) for some d >= 1
                         <=>   H2 holds, or u.v is a ROOT OF UNITY.
  In connection coordinates: the SINGULAR RESONANT SET is  (c - f)/2pi in Q  --  a countable
  DENSE family of circles parallel to the diagonal, of which R-5's c = f is the member d = 1.

  On L = Z.(d,d) the closure is d parallel circles, the Birkhoff sum interleaves d branches
  by k mod d, and the limit is
        lambda_d = (1/d) sum_{r<d} m(0.3 w^2 + 0.4 zeta^r w + 0.3 zeta^r),   zeta = e(j/d).
  Exactly ONE branch (r = 0) is the singular inhomogeneous-Sudler pair of R-5; the other
  d-1 are non-singular.  So even on the singular set, the object is a Sudler pair only for
  d = 1; for d >= 2 it is a d-fold interleaving one of whose branches is one.

AND lambda_d -> m(P) AS d -> infinity: a Boyd-Lawton-type accumulation in the SECOND
parameter, running orthogonally to M1_07's accumulation in |m|+|n|.  M1_07's table
accumulates on m(P) along growing PRIMITIVE relations; this one accumulates along growing
INDEX at a FIXED primitive relation.  Neither the register's 16-row table nor M1_07's
accumulation table contains a single row of it.
"""
import numpy as np
from X_lib import PI_K1, m_maxform, m_one_var

p00, p10, p01, p11 = PI_K1
MP = m_maxform(PI_K1, 1 << 24)
print("=" * 79); print("X_05 — THE SECOND ACCUMULATION"); print("=" * 79)
print("\nm(P) = %.12f" % MP)
print("\n   d      lambda_d              lambda_d - m(P)     singular branches (r with unimodular roots)")
for d in [1, 2, 3, 4, 5, 8, 13, 21, 34, 55, 100, 200, 500, 1000]:
    tot, singb = 0.0, 0
    for r in range(d):
        z = np.exp(2j * np.pi * r / d)
        c = np.array([p01 * z, p11 * z, p10], dtype=complex)
        tot += m_one_var(c)
        rr = np.roots(c[::-1])
        if np.min(np.abs(np.abs(rr) - 1.0)) < 1e-12:
            singb += 1
    lam = tot / d
    print("   %-6d %.12f   %+.3e         %d of %d" % (d, lam, lam - MP, singb, d))
print("""
   READ.  lambda_1 = log(0.3) is the register's (1,1) row.  lambda_2, lambda_3, ... appear
   nowhere in this corpus.  Exactly one branch is singular at every d, so the singularity
   never goes away -- the limit merely dilutes it by 1/d.  lambda_d -> m(P) like O(1/d).
   CONSEQUENCE FOR N1: an exceptional-set statement attached to N1 must exclude
   (c-f)/2pi rational, not merely c = f; the excluded set is still Haar-null, so N1's
   scope is unchanged -- what changes is the description of where it fails.
""")
print("DONE X_05")
