#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — R1.  THE DECISIVE COMPLEX EXPERIMENT, RE-RUN WITH EVALUATORS
THAT DO NOT ASSUME THE JENSEN REDUCTION.

WHY THIS EXISTS.  The lane's entire attribution finding (B-02, B-03) lives in LEG B, and
every number in LEG B comes from m_jensen / m_fast, i.e. from
    m(P) = (1/2pi) INT log max(|p0+p1 x|, |p2+p3 x|) dt.
Its exact evaluator E3 covers only real non-negative dominated arrays and its ergodic
evaluator E4 was run only in LEG A (real arrays).  So the complex half of the lane is
single-evaluator.  Here it is re-run with the raw double integral and with the corpus's
own ergodic method, neither of which knows what a Jensen pairing is.

THE ONE VARIABLE in A.1: the permutation of the four complex coefficients.  The moduli,
the arguments, the evaluator, the grid and the code path are identical in every row of a
block; only the assignment moves, and every array is printed.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlib import (m_grid, m_grid_rich, m_erg, m_split, fluxes, matching_of, blocks,
                  apply_perm, PERMS, is_subgroup, max_collinear, dominated_pairing,
                  MATCH_NAME, hdr)

np.set_printoptions(precision=6, suppress=True)

print(__doc__)

# --------------------------------------------------------------------------- A.0
hdr('A.0  MY EVALUATORS AGAINST VALUES THAT ARE KNOWN EXACTLY (real arrays)')
exact = [
    ('B0b   (4,2,1,2)/9  SENSE U', [4 / 9, 2 / 9, 1 / 9, 2 / 9], np.log(4 / 9)),
    ('B4    (1,1,1,3)/6  SENSE U', [1 / 6, 1 / 6, 1 / 6, 1 / 2], np.log(1 / 2)),
    ('SENSE C (1,1,1,1)/4       ', [.25, .25, .25, .25], np.log(1 / 4)),
    ('corpus (0.4,0.3,0.2,0.1)  ', [.4, .3, .2, .1], np.log(0.4)),
]
for name, p, ex in exact:
    r, g1, g2 = m_grid_rich(p, 1024)
    e = m_erg(p)
    s = float(m_split(p))
    print('  %s  exact %.12f   grid(1024/2048->rich) %.9f   ergodic %.9f   split-mpmath %.15f' %
          (name, ex, r, e, s))
    print('     %s  deviations: grid %.2e   ergodic %.2e   split %.2e' % (' ' * len(name), r - ex, e - ex, s - ex))

# --------------------------------------------------------------------------- A.1
hdr('A.1  THE LANE\'S OWN B.3 ARRAYS, ALL 24 PERMUTATIONS, BY THE RAW DOUBLE INTEGRAL')
ARGS = (0.0, 0.7, 1.9, 0.3)
CASES = [
    ('B0b  (4,2,1,2)/9', [4 / 9, 2 / 9, 1 / 9, 2 / 9],
     {2: -0.767461183337109, 0: -0.801115558404581, 1: -0.810930216216329}),
    ('B4   (1,1,1,3)/6 ', [1 / 6, 1 / 6, 1 / 6, 1 / 2],
     {2: -0.693147180559945, 0: -0.693147180559945, 1: -0.693147180559945}),
    ('GEN2 (.34,.33,.32,.01)', [.34, .33, .32, .01],
     {2: -0.779742059445662, 0: -0.785797264441102, 1: -0.790472715034835}),
    ('GEN3 (.37,.29,.23,.11)', [.37, .29, .23, .11],
     {2: -0.803695834252179, 0: -0.857030888035461, 1: -0.913019920897394}),
]
print('  args = %s in every case.  Lane block values quoted from its sealed'
      ' b2_legB_complex.OUT.txt B.3.' % (ARGS,))
print('  Matching index: 0 = {00,10}|{01,11}, 1 = {00,01}|{10,11}, 2 = {00,11}|{10,01} (diagonal).')
for name, r, lane in CASES:
    p = [r[i] * np.exp(1j * ARGS[i]) for i in range(4)]
    fl = fluxes(p)
    print()
    print('  %s   |phi_a|,|phi_b|,|phi_c| = %.4f %.4f %.4f   maxcollinear=%d   '
          'some-pairing-dominated=%s' % (name, fl[0], fl[1], fl[2], max_collinear(p),
                                         dominated_pairing(p)))
    vals_g, vals_s = [], []
    for s in PERMS:
        q = apply_perm(p, s)
        vals_g.append(m_grid_rich(q, 1024)[0])
    tol = 3e-5
    sizes, reps, lab = blocks(vals_g, tol)
    print('     GRID (no Jensen anywhere): blocks %s at tol %.0e' % (sizes, tol))
    # is each block exactly one diagonal matching?
    okmatch = True
    for b in range(len(reps)):
        ms = set(matching_of(PERMS[i]) for i in range(24) if lab[i] == b)
        if len(ms) != 1:
            okmatch = False
    print('     each block is exactly one diagonal matching: %s' % okmatch)
    for b in range(len(reps)):
        idx = [i for i in range(24) if lab[i] == b]
        ms = sorted(set(matching_of(PERMS[i]) for i in idx))
        mrep = MATCH_NAME[ms[0]] if len(ms) == 1 else 'MIXED'
        lv = lane.get(ms[0]) if len(ms) == 1 else None
        spread = max(vals_g[i] for i in idx) - min(vals_g[i] for i in idx)
        # one ergodic evaluation per block, on the block's first permutation
        e = m_erg(apply_perm(p, PERMS[idx[0]]), N=4_000_000)
        sp = float(m_split(apply_perm(p, PERMS[idx[0]])))
        print('       block %d  n=%2d  %-16s grid %.9f  within-block spread %.2e'
              % (b, len(idx), mrep, reps[b], spread))
        print('                          ergodic N=4e6 %.9f   my-split(dps30) %.15f' % (e, sp))
        if lv is not None:
            print('                          LANE VALUE   %.15f   |grid-lane| %.2e   '
                  '|ergodic-lane| %.2e   |split-lane| %.2e'
                  % (lv, abs(reps[b] - lv), abs(e - lv), abs(sp - lv)))
    if len(reps) > 1:
        gaps = sorted(reps)
        print('     between-block gaps: %s' % ['%.5f' % (gaps[i + 1] - gaps[i]) for i in range(len(gaps) - 1)])
    stab = [PERMS[i] for i in range(24) if lab[i] == lab[0]]
    print('     stabiliser of the identity value: size %d, is a subgroup: %s'
          % (len(stab), is_subgroup(stab)))
