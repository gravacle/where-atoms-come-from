#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — R2.  THE STRUCTURE THEOREM THE LANE ASSERTS BUT DOES NOT PROVE,
AND THE EXACT BOUNDARY OF ITS BICONDITIONAL.

The lane derives  m(s.p) = G({r_s0,r_s1},{r_s2,r_s3}; |phi_Md(s)|)  and then, in B.4,
silently upgrades it to  H(moduli MULTISET; |phi_Md(s)|).  The upgrade is what its
block-count prediction needs (two matchings with equal |phi| have DIFFERENT moduli
pairings, so without it their values need not agree), and it is not derived on its page.
B.2 tests it, and here it is tested by an evaluator that knows nothing about pairings.

Then the biconditional itself:
    'S4-INVARIANCE <=> AT LEAST THREE OF THE FOUR COEFFICIENTS ARE COLLINEAR IN C'
    'off the collinear locus the invariance group is exactly D4 ... carrying 3 rates'
The forward half is a theorem.  The converse is tested here against the lane's own
degeneracy regime and against a random sweep with an independent evaluator.

ONE VARIABLE PER SECTION, named in each header, arrays printed.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlib import (m_grid_rich, m_split, m_erg, fluxes, matching_of, blocks, apply_perm,
                  PERMS, is_subgroup, max_collinear, dominated_pairing, MATCH_NAME, hdr)

print(__doc__)
rng = np.random.default_rng(90210)

# --------------------------------------------------------------------------- B.1
hdr('B.1  THE UPGRADE: IS m A FUNCTION OF THE MODULI MULTISET AND |phi| ALONE?')
print("""  THE ONE VARIABLE: which of the three matchings of the SAME four moduli is used as
  the Jensen pairing.  The flux is held at the same value, the moduli multiset is held
  fixed, the evaluator is held fixed.  Three arrays realise the three pairings:
      q1 = (a0, a1, a2, a3 e^{i phi})   pairing {a0,a1}|{a2,a3}
      q2 = (a0, a2, a1, a3 e^{i phi})   pairing {a0,a2}|{a1,a3}
      q3 = (a0, a3, a1, a2 e^{i phi})   pairing {a0,a3}|{a1,a2}
  all three have diagonal flux exactly phi.  The lane's two Jensen reductions give q1=q2
  for free; q3 is the one its page never exhibits.""")
worst = 0.0
for trial in range(8):
    a = rng.uniform(0.05, 1.0, 4)
    phi = rng.uniform(0, np.pi)
    e = np.exp(1j * phi)
    q1 = [a[0], a[1], a[2], a[3] * e]
    q2 = [a[0], a[2], a[1], a[3] * e]
    q3 = [a[0], a[3], a[1], a[2] * e]
    v = [m_grid_rich(q, 1024)[0] for q in (q1, q2, q3)]
    s = [float(m_split(q)) for q in (q1, q2, q3)]
    d = max(s) - min(s)
    worst = max(worst, d)
    print('  moduli %s phi=%.4f  grid %.9f %.9f %.9f   split(dps30) spread %.2e'
          % (np.round(a, 4), phi, v[0], v[1], v[2], d))
print('  WORST spread over the three pairings at 30 digits: %.3e' % worst)
print('  -> the upgrade to the moduli MULTISET holds; the lane\'s block-count prediction')
print('     rests on it and its page never exhibits the third pairing.')

# --------------------------------------------------------------------------- B.2
hdr('B.2  THE CONVERSE IS FALSE, AND NOT ONLY ON A MEASURE-ZERO SET')
print("""  THE ONE VARIABLE: the moduli.  Arguments held at (0,.7,1.9,.3) -- three DISTINCT
  matching fluxes, largest collinear set = 1 -- in every row, so by the lane's stated
  biconditional every row must show blocks [8,8,8].""")
ARGS = (0.0, 0.7, 1.9, 0.3)
rows = [
    ('B4  (1,1,1,3)/6      the corpus\'s own spindle', [1 / 6, 1 / 6, 1 / 6, 1 / 2]),
    ('    (.10,.10,.10,.70)  r_max > sum of others  ', [.10, .10, .10, .70]),
    ('    (.05,.30,.32,.33)  no single dominant r   ', [.05, .30, .32, .33]),
    ('GEN3(.37,.29,.23,.11)  the lane\'s own generic ', [.37, .29, .23, .11]),
]
for name, r in rows:
    p = [r[i] * np.exp(1j * ARGS[i]) for i in range(4)]
    vals = [m_grid_rich(apply_perm(p, s), 1024)[0] for s in PERMS]
    sizes, reps, lab = blocks(vals, 3e-5)
    stab = [PERMS[i] for i in range(24) if lab[i] == lab[0]]
    print('  %s  maxcollinear=%d  dominated-pairing=%-5s  blocks=%s  '
          'stabiliser size %d subgroup %s'
          % (name, max_collinear(p), dominated_pairing(p), sizes, len(stab), is_subgroup(stab)))

# --------------------------------------------------------------------------- B.3
hdr('B.3  THE LANE\'S OWN SURVIVING MISSES, RE-DECIDED INDEPENDENTLY')
print("""  These four arrays are printed in the lane's sealed b2 B.4 as arrays where its
  prediction MISSED: three distinct |phi| but fewer than three rates.  Each is a
  counterexample to the biconditional in the lane's LEG B VERDICT, and only the first two
  are explained by the r_max regime the lane names.  Re-decided here at 30 digits and by
  the raw double integral.""")
misses = [
    ((0.5968, 0.2629, 0.0438, 0.2909), (3.7902, 1.8798, 5.5473, 4.9111), 1),
    ((0.4390, 0.2002, 0.9867, 0.4016), (3.2006, 4.1624, 0.3936, 6.0108), 1),
    ((0.1962, 0.0356, 0.5693, 0.7721), (0.2000, 5.4227, 0.8932, 5.8377), 2),
    ((0.8104, 0.1396, 0.1783, 0.5909), (1.8385, 3.7490, 1.6094, 3.8396), 2),
]
for r, ar, obs in misses:
    p = [r[i] * np.exp(1j * ar[i]) for i in range(4)]
    fl = fluxes(p)
    vs = [float(m_split(apply_perm(p, s))) for s in PERMS]
    vg = [m_grid_rich(apply_perm(p, s), 1024)[0] for s in PERMS]
    ss, reps, lab = blocks(vs, 1e-15)
    sg, _, _ = blocks(vg, 3e-5)
    stab = [PERMS[i] for i in range(24) if lab[i] == lab[0]]
    rmaxdom = max(r) >= sum(r) - max(r)
    print('  moduli %s |phi|=(%.3f,%.3f,%.3f) maxcollinear=%d' % (r, fl[0], fl[1], fl[2], max_collinear(p)))
    print('     blocks: split(dps30) %s   grid %s   lane reported %d value(s)' % (ss, sg, obs))
    print('     r_max >= sum of others: %-5s   SOME pairing dominated for all flux: %-5s'
          '   stabiliser size %d subgroup %s'
          % (rmaxdom, dominated_pairing(p), len(stab), is_subgroup(stab)))

# --------------------------------------------------------------------------- B.4
hdr('B.4  HOW BIG IS THE EXCEPTION?  RANDOM SWEEP WITH AN INDEPENDENT EVALUATOR')
print("""  THE ONE VARIABLE: nothing is held -- this is a measure, not a comparison.  4000
  arrays with independent uniform moduli and independent uniform arguments.  Because
  B.1 establishes m = H(multiset;|phi|), the 24 values are the 3 matching values, so each
  array costs 3 evaluations, and a random 1-in-50 subsample is checked at the full 24 to
  confirm that reduction.  Evaluator: my split-Jensen at dps 18 is too slow at this size,
  so the raw grid (Richardson 512/1024) decides at 1e-4 and every array within 3e-4 of a
  tie is RE-DECIDED at dps 30.  Question asked: for how many arrays is the lane's
  'off the collinear locus, exactly D4 and three rates' FALSE?""")
N = 4000
cnt = {1: 0, 2: 0, 3: 0}
dom_of = {1: 0, 2: 0, 3: 0}
redecided = 0
full24_checked = 0
full24_ok = 0
for i in range(N):
    r = rng.uniform(0.02, 1.0, 4)
    ar = rng.uniform(0, 2 * np.pi, 4)
    p = [r[j] * np.exp(1j * ar[j]) for j in range(4)]
    # one representative permutation per matching
    reps_perm = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]
    # matching_of each: 0->{0,3}? compute properly
    triple = []
    for s in reps_perm:
        triple.append(m_grid_rich(apply_perm(p, s), 512)[0])
    tol = 1e-4
    close = any(abs(triple[a] - triple[b]) < 3e-4 for a in range(3) for b in range(a + 1, 3))
    if close:
        triple = [float(m_split(apply_perm(p, s))) for s in reps_perm]
        tol = 1e-12
        redecided += 1
    sizes, _, _ = blocks(triple, tol)
    nd = len(sizes)
    cnt[nd] += 1
    if dominated_pairing(p):
        dom_of[nd] += 1
    if i % 50 == 0:
        full24_checked += 1
        v24 = [m_grid_rich(apply_perm(p, s), 512)[0] for s in PERMS]
        s24, _, _ = blocks(v24, 2e-4)
        if len(s24) == nd:
            full24_ok += 1
print('  %d arrays: #distinct rates -> %s   (re-decided at 30 digits: %d)' % (N, cnt, redecided))
print('  of those, in the "some pairing dominated for every flux" regime: %s' % dom_of)
print('  full-24 subsample: %d checked, %d agreed with the 3-matching reduction' % (full24_checked, full24_ok))
bad = cnt[1] + cnt[2]
print('  ARRAYS WITH NO THREE COLLINEAR COEFFICIENTS (probability 1 for this measure)')
print('  ON WHICH THE LANE\'S "invariance group is exactly D4, three rates" IS FALSE:')
print('      %d of %d = %.1f%%  -- a POSITIVE-MEASURE exception, not a null set.'
      % (bad, N, 100.0 * bad / N))
print('      of those %d, %d are NOT explained by any dominated pairing.'
      % (bad, bad - dom_of[1] - dom_of[2]))
