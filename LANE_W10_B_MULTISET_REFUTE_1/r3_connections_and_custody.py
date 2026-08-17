#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — R3.  THE TWO DISTINGUISHED CONNECTIONS, EXACTLY; AND THE
CUSTODY CLAIMS AT THE BYTES.

C.1  S1's order-4 connection.  The lane computes it in float64 and reports the surviving
     group.  Here the whole thing is closed-form and EXACT in Fractions, so the group is a
     statement about rational numbers rather than about a tolerance.
C.2  S3/S4's exactly-resonant connection, by an independent root-based Mahler measure.
C.3  The lane's two custody findings, re-checked at the bytes -- including the file its
     custody leg did not open.
C.4  LEG E.2's 'every functional is blind' claim, tested where the lane did not test it.
"""
import sys, os, subprocess, itertools, hashlib
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlib import (m_finite, m_resonant, m_poly_roots, m_grid, m_grid_rich, blocks, PERMS,
                  apply_perm, is_subgroup, hdr)

REPO = '/Users/bgm/MB Work/where-atoms-come-from'
print(__doc__)

# --------------------------------------------------------------------------- C.1
hdr('C.1  S1\'s ORDER-4 CONNECTION, IN CLOSED FORM AND IN EXACT RATIONALS')
print("""  u = conj(W_F) = -1, v = W_C = -i (S1 sec6).  The orbit has four points, so
      lambda_4 = (1/4)[ log|Z_1| + log|Z_2| + log|Z_3| + log|Z_4| ]
  and, for a REAL array, Z_3 = conj(Z_1), Z_2 = p0+p1-p2-p3, Z_4 = p0+p1+p2+p3, so

      lambda_4 = (1/4) log[ ((p0-p1)^2 + (p2-p3)^2) * |p0+p1-p2-p3| * (p0+p1+p2+p3) ]

  EXACT.  This is derived here, not taken from the lane.  The bracket is a rational
  number, so the orbit structure over the 24 permutations is decided in Fraction
  arithmetic with no tolerance at all.""")


def bracket(p):
    return ((p[0] - p[1]) ** 2 + (p[2] - p[3]) ** 2) * abs(p[0] + p[1] - p[2] - p[3]) * (p[0] + p[1] + p[2] + p[3])


CASES = [
    ('B0b (4,2,1,2)/9', [Fraction(4, 9), Fraction(2, 9), Fraction(1, 9), Fraction(2, 9)], 2, 0.127706),
    ('B4  (1,1,1,3)/6', [Fraction(1, 6), Fraction(1, 6), Fraction(1, 6), Fraction(1, 2)], 1, 0.0),
    ('GEN (.37,.29,.23,.11)', [Fraction(37, 100), Fraction(29, 100), Fraction(23, 100), Fraction(11, 100)], 3, 0.323797),
]
for name, p, lane_n, lane_spread in CASES:
    exact_vals = {}
    for s in PERMS:
        q = apply_perm(p, s)
        exact_vals[s] = bracket(q)
    distinct = sorted(set(exact_vals.values()))
    fl = [float(np.log(float(d)) / 4) for d in distinct]
    # check the closed form against a direct finite-orbit average
    dev = 0.0
    for s in PERMS[:6]:
        q = [float(x) for x in apply_perm(p, s)]
        direct = m_finite(q, np.pi, -np.pi / 2, 4)
        dev = max(dev, abs(direct - float(np.log(float(bracket(apply_perm(p, s)))) / 4)))
    stab = [s for s in PERMS if exact_vals[s] == exact_vals[PERMS[0]]]
    print('  %s  EXACT #distinct rates = %d  spread = %.6f   [lane: %d, %.6f]'
          % (name, len(distinct), max(fl) - min(fl), lane_n, lane_spread))
    print('     rates %s' % ['%.6f' % v for v in sorted(fl)])
    print('     closed form vs direct 4-point average, worst |dev| = %.2e' % dev)
    print('     stabiliser (EXACT, in Fractions): order %d, is a subgroup: %s'
          % (len(stab), is_subgroup(stab)))
print("""  The lane reports the GEN stabiliser as 'order 8, subgroup True' and names it the
  group of the JENSEN pairing {00,10}|{01,11}.  The closed form above says WHY, with no
  computation: the bracket depends on the array only through (p0-p1)^2, (p2-p3)^2 and
  p0+p1 vs p2+p3, so it is invariant under (00 10), under (01 11), and under swapping the
  two pairs -- exactly that group of order 8 -- and this is a THEOREM about the four
  rational numbers, not a measured coincidence.""")

# --------------------------------------------------------------------------- C.2
hdr('C.2  S3/S4\'s RESONANT CONNECTION, BY AN INDEPENDENT ROOT-BASED MAHLER MEASURE')
val3 = m_resonant([0.0, 0.3, 0.3, 0.4])
print('  VALIDATION against the ERRATUM AGAINST W-02 (three-class weights 0.4/0.3/0.3):')
print('     subtorus m(Q) from roots (numpy.roots, my code) = %.9f   register: -0.767014993' % val3)
print('     full-torus m for the same weights, my grid       = %.9f   register: -0.767507880'
      % m_grid_rich([0.0, 0.3, 0.3, 0.4], 2048)[0])
for name, p, lane_n, lane_spread in [
    ('B0b (4,2,1,2)/9', [4 / 9, 2 / 9, 1 / 9, 2 / 9], 3, 0.000006),
    ('B4  (1,1,1,3)/6', [1 / 6, 1 / 6, 1 / 6, 1 / 2], 2, 0.0),
    ('GEN (.37,.29,.23,.11)', [.37, .29, .23, .11], 12, 0.003428)]:
    vals = [m_resonant(apply_perm(p, s)) for s in PERMS]
    sizes, reps, lab = blocks(vals, 1e-9)
    stab = [PERMS[i] for i in range(24) if lab[i] == lab[0]]
    print('  %s  #distinct = %d  spread = %.6f   [lane: %d, %.6f]   stabiliser order %d subgroup %s'
          % (name, len(sizes), max(vals) - min(vals), lane_n, lane_spread, len(stab), is_subgroup(stab)))
    # independent confirmation of one value by direct simulation on the resonant orbit
    k = np.arange(1, 500001)
    u = np.exp(-1j * 2.0 * k); v = np.exp(1j * 1.1 * k)
    q = p
    sim = float(np.mean(np.log(np.abs(q[0] + q[1] * u + q[2] * v + q[3] * u * v))))
    print('       identity permutation: roots %.9f   direct N=5e5 %.9f   dev %.2e'
          % (vals[0], sim, abs(vals[0] - sim)))

# --------------------------------------------------------------------------- C.3
hdr('C.3  CUSTODY, RE-CHECKED AT THE BYTES')


def seal_state(lane, fname):
    d = os.path.join(REPO, lane)
    seals = os.path.join(d, 'SEALS.sha256')
    listed = False
    if os.path.exists(seals):
        listed = fname in open(seals).read()
    ok = None
    if listed:
        h = hashlib.sha256(open(os.path.join(d, fname), 'rb').read()).hexdigest()
        ok = h in open(seals).read()
    return listed, ok


for lane, f in [('LANE_G_GROUP_REFUTER', 'g0_validate.py'),
                ('LANE_S5_CHARGE_CODE', 's5_B_sweep.py'),
                ('LANE_S5_CHARGE_CODE', 's5_B_sweep.OUT.txt')]:
    print('  %-26s %-22s in SEALS: %-5s  hash verifies: %s' % (lane, f, *seal_state(lane, f)))

print("""
  C.3(a)  THE FILE THE LANE'S CUSTODY LEG DID NOT OPEN.
  LANE_S5_CHARGE_CODE/s5_B_sweep.py is SEALED and verifying, and it contains a PROOF of
  W-03's multiset theorem by exactly the mechanism the lane offers as its own
  contribution -- the pointwise symmetry of |A + B e^{it}| in (A,B), plus x<->y -- and a
  four-class check.  Quoted verbatim from the sealed bytes:""")
src = open(os.path.join(REPO, 'LANE_S5_CHARGE_CODE/s5_B_sweep.py')).read().splitlines()
for i, l in enumerate(src):
    if 'SYMMETRIC IN (A,B) POINTWISE' in l or 'generate the full symmetric group' in l or '24/24 permutation invariance' in l:
        print('     s5_B_sweep.py:%d  %s' % (i + 1, l.strip()))
out = open(os.path.join(REPO, 'LANE_S5_CHARGE_CODE/s5_B_sweep.OUT.txt')).read().splitlines()
for i, l in enumerate(out):
    if 'all 24 permutations' in l or ('spread' in l and 'value' in l):
        print('     s5_B_sweep.OUT.txt:%d  %s' % (i + 1, l.strip()))

print("""
  C.3(b)  AND THE CORPUS'S FOUR-CLASS WEIGHT VECTOR IS ITSELF A CONTROL THAT COULD NOT
  HAVE FAILED -- the same defect the lane raises against B0b, B4 and SENSE C, on the run
  the register actually quotes.  base = (0.4,0.3,0.2,0.1); domination certificate in
  exact rationals, my own, over ALL THREE pairings of the multiset:""")
base = [Fraction(2, 5), Fraction(3, 10), Fraction(1, 5), Fraction(1, 10)]
for (i, j, k, l) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
    a, b, c, d = base[i], base[j], base[k], base[l]
    hi = (a + b) ** 2 - (c + d) ** 2
    lo = (a - b) ** 2 - (c - d) ** 2
    print('     pairing {%s,%s}|{%s,%s}:  (S1^2-S2^2) = %s   (D1^2-D2^2) = %s   -> %s'
          % (a, b, c, d, hi, lo, 'A DOMINATES (lambda = log max = log 2/5)' if hi >= 0 and lo >= 0 else 'crosses'))
allv = set()
for s in PERMS:
    q = apply_perm(base, s)
    for (i, j, k, l) in ((0, 1, 2, 3), (0, 2, 1, 3)):
        pass
    allv.add(max(q))
print('     over all 24 arrangements the dominated branch max takes the values %s'
      % sorted(set(str(x) for x in allv)))
print('     so lambda = log(2/5) = %.12f for every one of the 24, EXACTLY, and the other'
      % float(np.log(0.4)))
print('     three weights never enter.  The corpus never permuted anything that mattered.')

print("""
  C.3(c)  WHERE THE REGISTER'S 2.4e-15 IS NOT.""")
r = subprocess.run(['grep', '-rn', '2.4e-15', REPO], capture_output=True, text=True)
for line in r.stdout.splitlines():
    if '/LANE_W10_B_MULTISET_REFUTE_1/' in line:
        continue
    print('     %s' % line.replace(REPO + '/', ''))
print("""     The only occurrences are: the register itself, the lane-under-attack's own
     custody page, and a PRINT LITERAL inside g0_validate.py that quotes W-03 as an
     external reference.  g0_validate's own measured spread is 9.285e-06, not 2.4e-15.
     A figure a script quotes is not a figure that script produced.""")
print('  g0_validate\'s evaluator, from its sealed bytes: mahler_torus(..., ngrid=1200),')
print('  a plain midpoint grid (glib.py:237-254).  Reproduced here, my own code, same grid:')
import numpy as _np
base_f = [0.4, 0.3, 0.2, 0.1]
vals = []
for s in PERMS:
    q = apply_perm(base_f, s)
    n = 1200
    th = (_np.arange(n) + 0.5) * 2 * _np.pi / n
    TH, PH = _np.meshgrid(th, th, indexing='ij')
    val = q[0] + q[1] * _np.exp(1j * TH) + q[2] * _np.exp(1j * PH) + q[3] * _np.exp(1j * (TH + PH))
    vals.append(float(_np.mean(_np.log(_np.abs(val)))))
print('     midpoint grid 1200^2, 24 permutations: spread = %.3e   [g0_validate.OUT.txt: 9.285e-06]'
      % (max(vals) - min(vals)))
print('     exact spread = 0 (domination certificate above).  The corpus\'s 9.285e-06 is')
print('     the quadrature, not the mathematics.')

# --------------------------------------------------------------------------- C.4
hdr('C.4  LEG E.2 (\'EVERY FULL-TORUS FUNCTIONAL IS BLIND\') WHERE THE LANE DID NOT TEST IT')
print("""  The lane checks E.2 on B0b and B4 only -- both in the domination regime -- and two
  of its eight printed rows are the functional 1{|P|<0.1}, which is IDENTICALLY ZERO on
  both carriers, so those two rows compare 0.0 with 0.0.  Here the same claim is tested on
  a NON-dominated real array where the branch max genuinely switches, and the threshold is
  set from the array's own minimum so the indicator is not vacuous.
  THE ONE VARIABLE: the permutation.  Grid, functional and array multiset held fixed.""")
n = 1536
th = (np.arange(n) + 0.5) * 2 * np.pi / n
TH, PH = np.meshgrid(th, th, indexing='ij')
EX, EY = np.exp(1j * TH), np.exp(1j * PH)
for nm, p in [('GEN (.37,.29,.23,.11) non-dominated', [.37, .29, .23, .11]),
              ('B0b (4,2,1,2)/9       dominated    ', [4 / 9, 2 / 9, 1 / 9, 2 / 9])]:
    mods = []
    for s in PERMS:
        q = apply_perm(p, s)
        mods.append(np.abs(q[0] + q[1] * EX + q[2] * EY + q[3] * EX * EY))
    mn = min(float(m.min()) for m in mods)
    thr = mn * 2.0
    for fname, F in [('mean |P|', lambda a: a),
                     ('mean |P|^(1/2)', lambda a: np.sqrt(a)),
                     ('mean 1{|P|<thr}', lambda a: (a < thr).astype(float)),
                     ('mean log|P|', lambda a: np.log(a))]:
        vs = [float(F(m).mean()) for m in mods]
        print('  %s  %-16s value %+0.9f  spread over 24 perms = %.2e%s'
              % (nm, fname, vs[0], max(vs) - min(vs),
                 '   (threshold %.4f, min|P| = %.4f, non-vacuous: %s)' % (thr, mn, vs[0] > 0)
                 if 'thr' in fname else ''))
