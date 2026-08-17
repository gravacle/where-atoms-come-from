#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — R5.  THE EXCEPTION, NAMED EXACTLY AND MEASURED.

TWO DEFECTS OF MY OWN, RECORDED HERE RATHER THAN PATCHED AWAY:
  (1) r2's B.4 used three 'representative' permutations (0,1,2,3), (0,2,1,3), (0,3,1,2).
      The first two induce the SAME diagonal matching {00,11}|{10,01}, so two of my three
      arms were the same arm and the sweep could not produce three values.  Its printed
      counts {1:720, 2:3280, 3:0} are that bug and nothing else.  r2's B.1, B.2 and B.3 do
      not use those representatives and stand.
  (2) r2's B.3 labelled the lane's four surviving misses 'NOT explained by any dominated
      pairing'.  That used the FLUX-INDEPENDENT domination test |a-b| >= c+d.  The correct
      test is per-flux, and by it all four ARE explained.  The corrected test is derived
      and used below.

THE EXACT TEST.  For a Jensen pairing {i,j}|{k,l} at flux phi,
    SA(t) - SB(t) = C0 + Re(K e^{it}),   C0 = a_i^2+a_j^2-a_k^2-a_l^2,
    |K(phi)|^2 = 4 a_i^2 a_j^2 + 4 a_k^2 a_l^2 - 8 a_i a_j a_k a_l cos(phi).
One branch dominates pointwise iff |C0| >= |K(phi)|, and then
    lambda = log max(a_i,a_j)   (if C0 >= 0)   or   log max(a_k,a_l)   (if C0 <= 0),
EXACTLY, with no quadrature and with no dependence on phi.  So lambda is LOCALLY CONSTANT
IN THE FLUX on an open set of arrays, which is precisely why the lane's biconditional fails
on a set of positive measure rather than on a null set.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlib import (m_grid_rich, m_split, blocks, apply_perm, PERMS, fluxes, max_collinear, hdr)

print(__doc__)

MATCH_REP = {2: (0, 1, 2, 3),    # diagonal {00,11}|{10,01}
             0: (0, 2, 3, 1),    # {00,10}|{01,11}
             1: (0, 1, 3, 2)}    # {00,01}|{10,11}
for m, s in MATCH_REP.items():
    assert len({frozenset((s[0], s[3]))}) == 1
assert len({frozenset((s[0], s[3])) for s in MATCH_REP.values()}) == 3, 'reps must induce 3 different matchings'


def dominated_value(q, pairing):
    """EXACT lambda by the per-flux domination test, or None.  q is the (permuted) array,
    pairing is 'y' -> {0,1}|{2,3} or 'x' -> {0,2}|{1,3}."""
    idx = (0, 1, 2, 3) if pairing == 'y' else (0, 2, 1, 3)
    a = [abs(complex(q[i])) for i in idx]
    A = [np.angle(complex(q[i])) for i in idx]
    C0 = a[0] ** 2 + a[1] ** 2 - a[2] ** 2 - a[3] ** 2
    phi = (A[3] - A[2]) - (A[1] - A[0])
    K2 = 4 * a[0] ** 2 * a[1] ** 2 + 4 * a[2] ** 2 * a[3] ** 2 - 8 * a[0] * a[1] * a[2] * a[3] * np.cos(phi)
    K = np.sqrt(max(K2, 0.0))
    if C0 >= K:
        return float(np.log(max(a[0], a[1])))
    if -C0 >= K:
        return float(np.log(max(a[2], a[3])))
    return None


def lam(q):
    """(value, 'exact'|'quad')"""
    for pr in ('y', 'x'):
        v = dominated_value(q, pr)
        if v is not None:
            return v, 'exact'
    return m_grid_rich(q, 512)[0], 'quad'


# --------------------------------------------------------------------------- E.1
hdr('E.1  THE LANE\'S FOUR SURVIVING MISSES, EXPLAINED EXACTLY')
misses = [((0.5968, 0.2629, 0.0438, 0.2909), (3.7902, 1.8798, 5.5473, 4.9111)),
          ((0.4390, 0.2002, 0.9867, 0.4016), (3.2006, 4.1624, 0.3936, 6.0108)),
          ((0.1962, 0.0356, 0.5693, 0.7721), (0.2000, 5.4227, 0.8932, 5.8377)),
          ((0.8104, 0.1396, 0.1783, 0.5909), (1.8385, 3.7490, 1.6094, 3.8396))]
for r, ar in misses:
    p = [r[i] * np.exp(1j * ar[i]) for i in range(4)]
    fl = fluxes(p)
    line = []
    for m in (0, 1, 2):
        q = apply_perm(p, MATCH_REP[m])
        v, how = lam(q)
        vs = float(m_split(q))
        line.append('%s %.12f (%s, dps30 %.12f)' % ('abc'[m], v, how, vs))
    print('  moduli %s  |phi|=(%.3f,%.3f,%.3f)' % (r, fl[0], fl[1], fl[2]))
    for l in line:
        print('     lambda_%s' % l)
print("""  Every one of the four is in the per-flux domination regime on at least the two
  matchings that share a value: lambda there is log(max modulus of the dominant branch),
  a number in which neither the flux nor three of the four moduli appear.  So the
  exception is not an accident of the moduli -- it is an OPEN REGION of coefficient
  space on which lambda is constant in the flux.""")

# --------------------------------------------------------------------------- E.2
hdr('E.2  THE MEASURE OF THE EXCEPTION, WITH THE SWEEP BUG FIXED')
print("""  THE ONE VARIABLE: none -- this is a measure.  1500 arrays, moduli and arguments
  independent uniform, three matching representatives that DO induce three different
  diagonal matchings (asserted at the top of this file, which is how the earlier bug
  would have been caught).  Each matching's rate is EXACT where the per-flux domination
  test fires and by the raw grid otherwise; any pair within 3e-4 is re-decided at dps 30.""")
rng = np.random.default_rng(90210)
N = 1500
cnt = {1: 0, 2: 0, 3: 0}
expl = {'all three fluxes in a domination regime': 0, 'other': 0}
redec = 0
for _ in range(N):
    r = rng.uniform(0.02, 1.0, 4)
    ar = rng.uniform(0, 2 * np.pi, 4)
    p = [r[j] * np.exp(1j * ar[j]) for j in range(4)]
    vals, hows = [], []
    for m in (0, 1, 2):
        v, how = lam(apply_perm(p, MATCH_REP[m]))
        vals.append(v)
        hows.append(how)
    if any(abs(vals[a] - vals[b]) < 3e-4 for a in range(3) for b in range(a + 1, 3)) and 'quad' in hows:
        vals = [float(m_split(apply_perm(p, MATCH_REP[m]))) for m in (0, 1, 2)]
        redec += 1
        tol = 1e-12
    else:
        tol = 1e-5
    sizes, _, _ = blocks(vals, tol)
    nd = len(sizes)
    cnt[nd] += 1
    if nd < 3:
        expl['all three fluxes in a domination regime' if all(h == 'exact' for h in hows) else 'other'] += 1
print('  %d arrays: #distinct rates -> %s   (re-decided at dps 30: %d)' % (N, cnt, redec))
print('  of the %d arrays with FEWER than three rates, the classification is %s'
      % (cnt[1] + cnt[2], expl))
print("""  THE LANE'S SENTENCE UNDER TEST: 'off [the three-collinear] locus the invariance group
  is exactly D4 with three rates indexed by the Newton square's diagonal matching.'
  Every array in this sweep is off that locus with probability one.  It is FALSE on
  %d of %d = %.1f%% of them.  That is a positive-measure exception, and the lane's own
  sealed B.4 sweep prints the same thing (its counts {3:304, 1:78, 2:18} on 400 arrays)
  one section above the verdict that states the biconditional."""
      % (cnt[1] + cnt[2], N, 100.0 * (cnt[1] + cnt[2]) / N))
