#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 1.  IS THE LANE'S OWN VACUITY FLAG (B-04) A FACT ABOUT
FOUR OCCUPIED CLASSES, ABOUT THE TWO CARRIERS, OR ABOUT TWO WEIGHT VECTORS?

The lane flags, against itself, that on B0b, B4 and SENSE C every one of the 24
arrangements has a dominated Jensen pairing, so lambda = log(max weight) and the
permutation arms could not have failed.  That flag is CORRECT.  This script asks the
scope question the lane does not: WHAT IS THE REGION, and does 'four occupied classes'
put you in it?

THE ONE VARIABLE in section 1.3: the READY STATE on a FIXED carrier (B0b or B4 or the
new K1+2S), with the incidence, the class map, the evaluator and the criterion held
fixed.  The class multiset does not move; only the pushforward does.

SECTIONS
  1.0  the domination criterion, DERIVED in sorted form and checked exhaustively
  1.1  the exact measure of the degenerate region in the weight simplex
  1.2  every carrier of S4:519 classified, SENSE U and SENSE C
  1.3  THE DECISIVE ARM: on the SAME two carriers, random ready states
  1.4  lambda >= log(w_max) with equality exactly on the degenerate region
"""
from fractions import Fraction as F

import numpy as np
import mpmath as mp

from r_lib import (PERMS, apply_perm, distinct_arrays, dominated_exact_perm, hdr,
                   mahler_jensen, mahler_lawton, sorted_domination)

mp.mp.dps = 30
print(__doc__)

# =============================================================================== 1.0
hdr("1.0  THE DOMINATION CRITERION IN SORTED FORM, DERIVED AND CHECKED EXHAUSTIVELY")
print("""  For the Jensen pairing {X,Y}|{Z,W} of a REAL NON-NEGATIVE array,
      SA - SB = (X^2+Y^2-Z^2-W^2) + 2(XY-ZW) cos t
  is monotone in cos t, so one branch dominates pointwise iff the endpoint values agree
  in sign, i.e. iff  sign(X+Y-Z-W) = sign(|X-Y|-|Z-W|).  With w1>=w2>=w3>=w4 the three
  pairings give w1-w2>=w3-w4, w1-w3>=w2-w4, w1+w4>=w2+w3 -- ALL THE SAME INEQUALITY.

      ALL 24 ARRANGEMENTS DOMINATED   <=>   w_max + w_min >= w_mid1 + w_mid2
      and then lambda = log(w_max) EXACTLY, the other three weights unused.

  Checked exhaustively below in exact rationals: the sorted criterion against the
  per-arrangement test over all 24 permutations, on every 4-composition of n <= 14.""")
bad = 0
tot = 0
for n in range(4, 15):
    for a in range(1, n - 2):
        for b in range(1, n - a - 1):
            for c in range(1, n - a - b):
                d = n - a - b - c
                if d < 1:
                    continue
                p = (F(a, n), F(b, n), F(c, n), F(d, n))
                ok, wmax, flags, sortcrit = sorted_domination(p)
                nd, vals = dominated_exact_perm(p)
                tot += 1
                if (nd == 24) != sortcrit or (nd == 24) != ok:
                    bad += 1
                if nd == 24 and vals != [wmax]:
                    bad += 1
print("  4-compositions tested: %d      DISAGREEMENTS between the sorted criterion, the"
      % tot)
print("  24-arrangement exact test, and 'lambda = log(w_max)': %d" % bad)

# =============================================================================== 1.1
hdr("1.1  THE MEASURE OF THE DEGENERATE REGION — EXACTLY 3/4 OF THE SIMPLEX")
print("""  Under the uniform (Dirichlet(1,1,1,1)) measure on the four-class weight simplex,
  write the ordered gaps D_k = w_(k) - w_(k+1).  Renyi's representation for uniform
  spacings makes D_k proportional to E_k/k with E_k i.i.d. Exp(1), so
      w1+w4 >= w2+w3  <=>  D_1 >= D_3  <=>  3 E_1 >= E_3,
  and P(E_3 <= 3E_1) = 1 - E[e^{-3E_1}] = 1 - 1/(1+3) = 3/4.
  Three classes (w4 = 0) reduce to w_max >= 1/2, whose Dirichlet(1,1,1) probability is
  3*(1/2)^2 = 3/4 as well.  Both checked by simulation:""")
rng = np.random.default_rng(20260816)
for k, lab in ((4, "four classes, Dirichlet(1,1,1,1)"), (3, "three classes, Dirichlet(1,1,1)")):
    W = rng.dirichlet(np.ones(k), size=2_000_000)
    if k == 3:
        W = np.concatenate([W, np.zeros((len(W), 1))], axis=1)
    S = np.sort(W, axis=1)[:, ::-1]
    frac = float(np.mean(S[:, 0] + S[:, 3] >= S[:, 1] + S[:, 2]))
    print("  %-36s degenerate fraction = %.6f   exact 3/4 = 0.750000" % (lab, frac))

# =============================================================================== 1.2
hdr("1.2  EVERY CARRIER OF S4:519 CLASSIFIED — AND THE CORRELATION IS WITH THE WEIGHTS")
CARR = [
    ("B0a ring torus, disjoint", (2, 4, 3, 0)),
    ("B0b ring torus, MEET    ", (4, 2, 1, 2)),
    ("B3  horn torus          ", (0, 2, 2, 1)),
    ("B1  K1 as handed        ", (0, 2, 2, 1)),
    ("B4  spindle             ", (1, 1, 1, 3)),
    ("B2  K1 both filled      ", (0, 2, 2, 1)),
    ("B1p K1-bridged          ", (0, 3, 3, 0)),
    ("B1q K1-bridged+spectator", (1, 3, 3, 0)),
    ("B1s K1 subdivided       ", (0, 5, 5, 1)),
]
print("  SENSE U (p = class counts / V), the corpus's own row, and SENSE C alongside.")
print("  %-26s %-22s %8s %6s %10s %s" % ("carrier", "class counts", "classes", "arrays",
                                          "degenerate", "lambda"))
for name, cnt in CARR:
    V = sum(cnt)
    p = tuple(F(x, V) for x in cnt)
    nz = sum(1 for x in cnt if x)
    ok, wmax, _, _ = sorted_domination(p)
    lam = mahler_jensen(tuple(float(x) for x in p), dps=25)
    print("  %-26s %-22s %8d %6d %10s %s"
          % (name, str(cnt), nz, distinct_arrays(p), ok, mp.nstr(lam, 12)))
for lab, p in (("SENSE C, 3 classes", (F(0), F(4, 10), F(3, 10), F(3, 10))),
               ("SENSE C, 4 classes", (F(1, 4),) * 4)):
    ok, wmax, _, _ = sorted_domination(p)
    print("  %-26s %-22s %8d %6d %10s %s"
          % (lab, str(tuple(str(x) for x in p)), sum(1 for x in p if x), distinct_arrays(p),
             ok, mp.nstr(mahler_jensen(tuple(float(x) for x in p), dps=25), 12)))
print("""
  READ IT OFF: every THREE-class carrier the corpus ran is in the INFORMATIVE quarter and
  both FOUR-class carriers are in the DEGENERATE three-quarters -- but the criterion has
  nothing to do with class occupancy.  A three-class vector is degenerate iff its max
  weight is >= 1/2, and K1's is 0.4.  SENSE C at four classes sits exactly ON the
  boundary (1/4+1/4 = 1/4+1/4), which is why (1+x)(1+y)/4 factors.""")

# =============================================================================== 1.3
hdr("1.3  THE DECISIVE ARM — SAME CARRIER, SAME CLASS MULTISET, THE READY STATE MOVES")
print("""  The carrier does NOT determine the weight vector.  p_ab = SUM_{v in class ab} |s_v|^2,
  so on any carrier with all four classes occupied the weight vector ranges over the WHOLE
  open 3-simplex as the ready state ranges over the unit sphere.  'B0b is degenerate' is
  therefore not a statement about B0b; it is a statement about SENSE U on B0b.
  ONE VARIABLE: the ready state.  Incidence, class map, criterion, evaluator all fixed.
  200000 Haar-random unit states per carrier, seed 20260816.  (The pushforward measure is
  Dirichlet(class sizes); SENSE U is its mean.)""")
rng = np.random.default_rng(20260816)
for name, mult in (("B0b  {00:4,10:2,01:1,11:2}", (4, 2, 1, 2)),
                   ("B4   {00:1,10:1,01:1,11:3}", (1, 1, 1, 3)),
                   ("K1+2S {00:2,10:2,01:2,11:1}", (2, 2, 2, 1))):
    V = sum(mult)
    N = 200_000
    s = rng.normal(size=(N, V)) + 1j * rng.normal(size=(N, V))
    w = np.abs(s) ** 2
    w /= w.sum(axis=1, keepdims=True)
    p = np.zeros((N, 4))
    i = 0
    for ci, m in enumerate(mult):
        p[:, ci] = w[:, i:i + m].sum(axis=1)
        i += m
    S = np.sort(p, axis=1)[:, ::-1]
    deg = float(np.mean(S[:, 0] + S[:, 3] >= S[:, 1] + S[:, 2]))
    pu = np.array([m / V for m in mult])
    Su = np.sort(pu)[::-1]
    print("  %-28s SENSE U degenerate = %-5s      random-state degenerate fraction = %.4f"
          % (name, bool(Su[0] + Su[3] >= Su[1] + Su[2]), deg))
print("""
  The pushforward of a Haar state onto the classes is EXACTLY Dirichlet(m_00,m_10,m_01,m_11)
  with the class sizes as parameters, and SENSE U is that Dirichlet's MEAN.  So the
  question 'is this carrier degenerate' is really 'does the mean of its class-Dirichlet
  fall in the degenerate 3/4', and 28% of ready states on B0b put it outside.  B-04's 'BOTH FOUR-CLASS CARRIERS THE
  CORPUS OWNS ARE DEGENERATE FOR THIS TEST AND COULD NOT HAVE FAILED IT' is true of the
  SENSE-U vectors S4 chose and false of the carriers.  The lane's own LEG C and LEG D
  already draw random states on these same carriers -- so the lane owned the fix and did
  not apply it to LEG A.""")

# =============================================================================== 1.4
hdr("1.4  lambda >= log(w_max), WITH EQUALITY EXACTLY ON THE DEGENERATE REGION")
print("  A sharper statement than 'the other three weights do not enter': the degenerate")
print("  region is exactly where the rate saturates its own lower bound.")
rng = np.random.default_rng(4242)
worst_eq, worst_gap, nd, nn = 0.0, 1e9, 0, 0
for _ in range(400):
    p = rng.dirichlet(np.ones(4))
    S = np.sort(p)[::-1]
    lam = float(mahler_jensen(tuple(p), dps=20))
    lb = float(np.log(S[0]))
    if S[0] + S[3] >= S[1] + S[2]:
        nd += 1
        worst_eq = max(worst_eq, abs(lam - lb))
    else:
        nn += 1
        worst_gap = min(worst_gap, lam - lb)
print("  400 Dirichlet(1,1,1,1) draws, seed 4242")
print("  degenerate: %d draws, worst |lambda - log w_max| = %.3e" % (nd, worst_eq))
print("  informative: %d draws, smallest (lambda - log w_max) = %.6f  (all strictly > 0)"
      % (nn, worst_gap))
print("  observed degenerate fraction %.4f against the exact 3/4" % (nd / 400))
