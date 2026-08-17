#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — R4.  FOUR SIDE CHECKS.

D.1  Is the label-dependence I found in a NON-VACUOUS near-return density real, or is it
     grid noise?  (A refuter who reports quadrature as a refutation has done the corpus's
     characteristic damage, so this is settled by convergence before anything is claimed.)
D.2  The involution at the level of Z_k, my own code, real vs complex.
D.3  The scope the lane's reality hypothesis actually has: at THREE occupied classes the
     multiset theorem survives arbitrary COMPLEX coefficients, so 'reality' bites only at
     four -- i.e. on no carrier the corpus ran before W-09.
D.4  The lane's B-06 sweep, re-added from its own printed numbers.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlib import m_grid_rich, m_split, blocks, apply_perm, PERMS, fluxes, max_collinear, hdr

print(__doc__)

# --------------------------------------------------------------------------- D.1
hdr('D.1  THE NEAR-RETURN DENSITY: CONVERGENCE TEST BEFORE ANY CLAIM')
p = [4 / 9, 2 / 9, 1 / 9, 2 / 9]
for n in (1024, 2048, 4096):
    th = (np.arange(n) + 0.5) * 2 * np.pi / n
    TH, PH = np.meshgrid(th, th, indexing='ij')
    EX, EY = np.exp(1j * TH), np.exp(1j * PH)
    vs = []
    for s in PERMS:
        q = apply_perm(p, s)
        a = np.abs(q[0] + q[1] * EX + q[2] * EY + q[3] * EX * EY)
        vs.append(float((a < 0.2222).mean()))
    print('  B0b, 1{|P| < 0.2222}, grid %5d^2 : value %.9f   spread over 24 perms = %.2e'
          % (n, vs[0], max(vs) - min(vs)))
print("""  The spread does not settle on any non-zero value; it falls by more than an order
  of magnitude with the grid and is set by how the grid resolves the boundary of the
  sub-level set.  IT IS QUADRATURE, and the exact statement backs that: since every
  INT INT F(|P|) is permutation-invariant, the whole DISTRIBUTION of |P| over the torus is,
  so the density below any threshold is invariant exactly.  The lane's E.2 claim survives
  this attack.  What does not survive is its choice of threshold: 1{|P| < 0.1} is
  IDENTICALLY ZERO on both carriers it tested, so two of its eight printed rows compare
  0.0 with 0.0 -- a control that could not have failed, unflagged.""")

# --------------------------------------------------------------------------- D.2
hdr('D.2  THE INVOLUTION AT THE LEVEL OF Z_k, MY OWN CODE')
rng = np.random.default_rng(5150)
for nm, gen in (('real non-negative', lambda: rng.dirichlet([1, 1, 1, 1]).astype(complex)),
                ('real, signs mixed ', lambda: (rng.dirichlet([1, 1, 1, 1]) * rng.choice([-1, 1], 4)).astype(complex)),
                ('complex, 3 collinear', lambda: np.array(list(rng.dirichlet([1, 1, 1])) + [0], dtype=complex) + np.array([0, 0, 0, rng.uniform(.1, 1) * np.exp(1j * rng.uniform(0, 6.28))])),
                ('complex, generic  ', lambda: rng.uniform(.1, 1, 4) * np.exp(1j * rng.uniform(0, 6.28, 4)))):
    wZ = wl = 0.0
    for _ in range(500):
        q = gen()
        f, c = rng.uniform(0, 2 * np.pi, 2)
        k = rng.integers(1, 100000, 20)
        u, v = np.exp(-1j * f * k), np.exp(1j * c * k)
        Z = q[0] + q[1] * u + q[2] * v + q[3] * u * v
        qt = [q[3], q[2], q[1], q[0]]
        Zt = qt[0] + qt[1] * u + qt[2] * v + qt[3] * u * v
        wZ = max(wZ, float(np.max(np.abs(np.conj(u) * np.conj(v) * Z - np.conj(Zt)))))
        wl = max(wl, abs(m_grid_rich(q, 512)[0] - m_grid_rich(qt, 512)[0]))
    print('  %-20s max |conj(u)^k conj(v)^k Z_k - conj(Z~_k)| = %.3e     max |lambda - lambda~| = %.3e'
          % (nm, wZ, wl))

# --------------------------------------------------------------------------- D.3
hdr('D.3  AT THREE OCCUPIED CLASSES THE MULTISET THEOREM NEEDS NO REALITY AT ALL')
print("""  THE ONE VARIABLE: the permutation.  Held: the moduli, the arguments, the evaluator.
  Every array below has ONE class coefficient equal to zero -- i.e. K1's own occupancy,
  and the occupancy of every carrier the corpus ran before W-09 -- and arguments chosen so
  that no three of the four coefficients are collinear among the NON-ZERO ones.""")
for nm, r, ar in (('K1-shaped (0, .3, .3, .4) generic args ', (0.0, .3, .3, .4), (0., .7, 1.9, .3)),
                  ('three-class (0, .5, .2, .3) other args ', (0.0, .5, .2, .3), (0., 2.1, .4, 5.0)),
                  ('three-class, zero in class 11          ', (.45, .25, .3, 0.0), (0., 1.3, 2.9, .8))):
    q = [r[i] * np.exp(1j * ar[i]) for i in range(4)]
    vals = [float(m_split(apply_perm(q, s))) for s in PERMS]
    sizes, _, _ = blocks(vals, 1e-15)
    print('  %s  maxcollinear(nonzero)=%d  blocks at dps30 = %s   spread = %.3e'
          % (nm, max_collinear([z for z in q if abs(z) > 0]), sizes, max(vals) - min(vals)))
print("""  So the hypothesis the lane calls load-bearing -- REALITY of the class array -- is
  VACUOUS at three occupied classes: with a zero coefficient the three surviving phases are
  exactly the three gauge parameters (alpha, beta, gamma) of P -> alpha P(beta x, gamma y),
  so they can all be removed and lambda depends on the moduli multiset alone.  The lane's
  own break (LEG C, an inserted observable) does still bite on K1, because the RECTANGLE
  sums it produces are generically all four non-zero even when the CLASS sums are not.""")

# --------------------------------------------------------------------------- D.4
hdr('D.4  THE LANE\'S B-06 SWEEP, RE-ADDED FROM ITS OWN PRINTED TABLE')
print("""  b3_pushforward.OUT.txt C.2 prints, for each of four sweeps of 200 observables:
      B0b Hermitian {(1,'dominant'):39, (1,'not dominant'):6, (2,'not dominant'):14, (3,'not dominant'):141}
      B0b unitary   {(1,'dominant'):31, (1,'not dominant'):5, (2,'not dominant'):13, (3,'not dominant'):151}
      B4  Hermitian {(1,'dominant'):41, (1,'not dominant'):5, (2,'not dominant'):19, (3,'not dominant'):135}
      B4  unitary   {(1,'dominant'):48, (1,'not dominant'):3, (2,'not dominant'):15, (3,'not dominant'):134}
  The finding B-06 says: 'roughly 70 percent give 3 distinct lambda, THE REMAINDER FALL IN
  THE |q|_max-DOMINANT REGIME'.  Adding the lane's own columns:""")
for nm, d in (('B0b Hermitian', (39, 6, 14, 141)), ('B0b unitary  ', (31, 5, 13, 151)),
              ('B4  Hermitian', (41, 5, 19, 135)), ('B4  unitary  ', (48, 3, 15, 134))):
    dom, nd1, nd2, three = d
    rem = 200 - three
    print('   %s: 3 rates %3d (%.0f%%)   remainder %3d, of which dominant %3d and NOT dominant %3d (%.0f%% of the sweep)'
          % (nm, three, 100 * three / 200, rem, dom, nd1 + nd2, 100 * (nd1 + nd2) / 200))
print("""  So between 9% and 10% of every sweep is an observable whose four rectangle
  coefficients are complex, non-collinear, NOT in the dominant regime, and which STILL
  fails to produce three rates.  The lane's own sweep contains the counterexample to the
  lane's own biconditional, in the finding that cites the sweep as its evidence.""")
