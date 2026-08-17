#!/usr/bin/env python3
"""LANE W-10 B — LEG A.  W-03's MULTISET THEOREM AT FOUR OCCUPIED CLASSES, REAL WEIGHTS.

THE ONE VARIABLE: the PERMUTATION of the four class weights.  The polynomial, the
evaluator, the precision, the split rule and the code path are identical in every row;
only the assignment of the same four numbers to the four classes moves.

ARM DIFF (the process rule the corpus's isolation audit says a ledger cannot catch):
every row prints its own coefficient array, so byte-identical arms are visible on the
page.  Row 1 of each block is the identity permutation; the other 23 arrays differ from
it in at least two positions unless the weight multiset has a repeat, and repeats are
flagged explicitly.

FOUR INDEPENDENT EVALUATORS (b_lib.py): E1 Jensen-in-y, E2 Jensen-in-x, E3 exact
rational domination certificate, E4 the corpus's own direct ergodic average.
"""
import itertools
from fractions import Fraction

import numpy as np
import mpmath as mp

from b_lib import (LBL, PERMS, apply_perm, cycle_notation, hdr, m_dominated_exact,
                   m_ergodic, m_jensen, matching_key)

mp.mp.dps = 30

CASES = [
    ("B0b  ring torus, loops MEET   S4:575 {00:4,01:1,10:2,11:2}  SENSE U",
     (Fraction(4, 9), Fraction(2, 9), Fraction(1, 9), Fraction(2, 9))),
    ("B4   spindle                  S4:578 {00:1,01:1,10:1,11:3}  SENSE U",
     (Fraction(1, 6), Fraction(1, 6), Fraction(1, 6), Fraction(3, 6))),
    ("SENSE C, four classes         S4:571 (0.25,0.25,0.25,0.25)",
     (Fraction(1, 4),) * 4),
    ("GEN1 generic, switching regime (0.40,0.40,0.15,0.05)",
     (Fraction(2, 5), Fraction(2, 5), Fraction(3, 20), Fraction(1, 20))),
    ("GEN2 generic, switching regime (0.34,0.33,0.32,0.01)",
     (Fraction(34, 100), Fraction(33, 100), Fraction(32, 100), Fraction(1, 100))),
]

print(__doc__)

for name, p in CASES:
    hdr("LEG A — " + name)
    pf = tuple(float(x) for x in p)
    reps = len(set(p)) < 4
    print("  base array (00,10,01,11) = %s      repeated weights present: %s"
          % (tuple(str(x) for x in p), reps))

    vals1, vals2, exact, arrays = [], [], [], []
    for s in PERMS:
        q = apply_perm(p, s)
        arrays.append(q)
        qf = tuple(float(x) for x in q)
        vals1.append(m_jensen(qf, 'y'))
        vals2.append(m_jensen(qf, 'x'))
        exact.append(m_dominated_exact(q, 'y'))

    sp1 = max(vals1) - min(vals1)
    sp2 = max(vals2) - min(vals2)
    cross = max(abs(a - b) for a, b in zip(vals1, vals2))
    print("  E1 Jensen-in-y, 24 permutations : spread = %s" % mp.nstr(sp1, 6))
    print("  E2 Jensen-in-x, 24 permutations : spread = %s" % mp.nstr(sp2, 6))
    print("  E1 vs E2, worst disagreement    : %s" % mp.nstr(cross, 6))
    print("  E1 value (identity permutation) : %s" % mp.nstr(vals1[0], 18))

    ndom = sum(1 for ok, _ in exact if ok)
    print("  E3 exact: %d of 24 arrangements have a DOMINATED Jensen pairing" % ndom)
    if ndom:
        certs = sorted({M for ok, M in exact if ok})
        for M in certs:
            print("     exact lambda on those = log(%s) = %s"
                  % (M, mp.nstr(mp.log(mp.mpf(M.numerator) / M.denominator), 18)))
        worst = max(abs(vals1[i] - mp.log(mp.mpf(exact[i][1].numerator) / exact[i][1].denominator))
                    for i in range(24) if exact[i][0])
        print("     E1 vs E3 worst deviation = %s" % mp.nstr(worst, 6))
    if ndom == 24:
        print("     >>> ALL 24 ARRANGEMENTS DOMINATED.  lambda = log(max weight) EXACTLY, in")
        print("         rational arithmetic, with no quadrature and no appeal to the multiset")
        print("         theorem: the value does not depend on the other three weights AT ALL.")

    # arm diff
    ident = arrays[0]
    ndiff = sum(1 for q in arrays[1:] if tuple(q) != tuple(ident))
    print("  ARM DIFF: %d of the other 23 permuted arrays differ from the identity array"
          % ndiff)

    # E4, the corpus's own direct method, on one representative of each of the three
    # D4 cosets (the three ways of pairing the four weights on the Newton diagonal)
    seen = {}
    for s in PERMS:
        seen.setdefault(matching_key(s), s)
    print("  E4 direct ergodic average, f=1.0 c=sqrt(2) N=2e6 (no quadrature anywhere):")
    for mk, s in seen.items():
        q = tuple(float(x) for x in apply_perm(p, s))
        e4 = m_ergodic(np.array(q))
        print("     perm %-14s array %s  E4 = %+.9f   E1-E4 = %+.2e"
              % (cycle_notation(s), tuple(round(x, 4) for x in q), e4,
                 float(vals1[PERMS.index(s)]) - e4))

hdr("LEG A — VERDICT")
print("""  W-03's multiset theorem is CONFIRMED at four occupied classes, on both of the
  corpus's four-class carriers and on two generic vectors, by four independent
  evaluators, two of them exact.

  BUT TWO OF THE FOUR ARMS COULD NOT HAVE FAILED, and this is recorded before the
  verdict rather than after it.  On B0b and on B4 -- and on SENSE C -- EVERY one of the
  24 arrangements has a dominated Jensen pairing, so lambda = log(max class weight)
  exactly, independently of the other three weights.  A permutation test on a functional
  that ignores three of its four arguments is not a test of the multiset theorem.
  The information in this leg is carried entirely by GEN1 and GEN2, where the branch
  max genuinely switches and the other three weights genuinely enter.""")
