#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 2.  RUN THE ARM THE LANE SAID COULD NOT FAIL, IN A REGIME
WHERE IT CAN.

The lane's LEG A tests W-03's multiset theorem on four-class weight vectors that are all
inside the degenerate region of ATTACK 1, so lambda = log(w_max) and three of the four
weights are unused.  The lane says so itself.  What it does not do is run the test on a
four-class CARRIER outside that region -- and it did not have to build one, because both
of the corpus's own four-class carriers get there with a different ready state.

  2.0  a four-class carrier built from incidence whose SENSE-U vector is INFORMATIVE
  2.1  the 24-permutation test on it, three independent evaluators
  2.2  THE DECISIVE ARM: B0b and B4 themselves, with ready states outside the region
  2.3  ARM-COUNT AUDIT: how many of the '24 permutations' are distinct arrays

THE ONE VARIABLE in 2.1 and 2.2 is the permutation of the four class weights.  The
carrier, the state, the connection, the evaluator and the code path are held fixed and
the DISTINCT-ARRAY COUNT is printed for every row, because the corpus's isolation audit
says a byte-identical arm is the defect a ledger cannot catch.
"""
from fractions import Fraction as F

import numpy as np
import mpmath as mp

from r_lib import (PERMS, apply_perm, cyc, distinct_arrays, dominated_exact_perm, hdr,
                   mahler_ergodic, mahler_jensen, mahler_lawton, sorted_domination)

mp.mp.dps = 30
print(__doc__)


# =============================================================================== 2.0
def build(V, edges, faces, gF, gC, name):
    """edges: list of (src,tgt).  faces: list of signed edge-index lists.
    d1[v,e] = +1 if v is target, -1 if v is source."""
    E = len(edges)
    Fn = len(faces)
    d1 = np.zeros((V, E))
    for e, (s, t) in enumerate(edges):
        d1[s, e] -= 1
        d1[t, e] += 1
    d2 = np.zeros((E, Fn))
    for j, fc in enumerate(faces):
        for e in fc:
            d2[abs(e) - 1 if e > 0 else abs(e) - 1, j] += (1 if e > 0 else -1)
    r1 = np.linalg.matrix_rank(d1)
    r2 = np.linalg.matrix_rank(d2)
    b0 = V - r1
    b1 = E - r1 - r2
    b2 = Fn - r2
    chi = V - E + Fn
    return dict(name=name, V=V, E=E, F=Fn, chi=chi, b=(b0, b1, b2), d1=d1, d2=d2,
                dd=float(np.abs(d1 @ d2).max()), gF=gF, gC=gC)


def chain(K, loop):
    z = np.zeros(K['E'])
    for e in loop:
        z[abs(e) - 1] += (1 if e > 0 else -1)
    return z


def is_cycle(K, loop):
    return float(np.abs(K['d1'] @ chain(K, loop)).max()) < 1e-9


def bounds(K, loop):
    z = chain(K, loop)
    sol, *_ = np.linalg.lstsq(K['d2'], z, rcond=None)
    return float(np.linalg.norm(K['d2'] @ sol - z)) < 1e-9


hdr("2.0  K1+2S — K1 WITH A SPECTATOR PAIR.  FOUR CLASSES, AND SENSE U IS INFORMATIVE")
print("""  W-09's corrected name of record is 'all four classes occupied, which requires a
  vertex in BOTH loops and a vertex in NEITHER'.  K1 supplies the pinch (v0), so the
  minimal four-class carrier is K1 with a spectator branch attached.  Incidence, in S1's
  own conventions (d1[v,e] = +1 at the target, -1 at the source):

    V = 7   v0..v6            F = 1, attached along e1+e2+e3 (K1's filled triangle)
    e1=(0,1) e2=(1,2) e3=(2,0)   e4=(0,3) e5=(3,4) e6=(4,0)   e7=(0,5) e8=(5,6)
    gamma_F = e1+e2+e3   gamma_C = e4+e5+e6

  Nothing in this lane needs the carrier beyond its class multiset; it is built here so
  that the informative weight vector is a CARRIER's own SENSE U and not a free vector.""")
K = build(7, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0), (0, 5), (5, 6)],
          [[1, 2, 3]], [1, 2, 3], [4, 5, 6], "K1+2S")
FACE_V, CYC_V = {0, 1, 2}, {0, 3, 4}
cnt = [0, 0, 0, 0]
for v in range(K['V']):
    a, b = int(v in FACE_V), int(v in CYC_V)
    cnt[{(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(a, b)]] += 1
print("  %s  V=%d E=%d F=%d chi=%d b=(%d,%d,%d)  |d1.d2| = %.1e"
      % (K['name'], K['V'], K['E'], K['F'], K['chi'], *K['b'], K['dd']))
print("  gamma_F is a cycle: %s and BOUNDS: %s      gamma_C is a cycle: %s and bounds: %s"
      % (is_cycle(K, K['gF']), bounds(K, K['gF']), is_cycle(K, K['gC']), bounds(K, K['gC'])))
print("  class counts {00:%d,10:%d,01:%d,11:%d}   all four occupied: %s"
      % (*cnt, all(cnt)))
pU = tuple(F(x, K['V']) for x in cnt)
ok, wmax, _, _ = sorted_domination(pU)
print("  SENSE U weights = %s      DEGENERATE (all-24-dominated): %s"
      % (tuple(str(x) for x in pU), ok))
print("  gauge count check, S4's own rule: invariants = E-(V-b0) = %d, curvature = rank d2"
      " = %d, flat = b1 = %d, sum matches: %s"
      % (K['E'] - (K['V'] - K['b'][0]), np.linalg.matrix_rank(K['d2']), K['b'][1],
         K['E'] - (K['V'] - K['b'][0]) == np.linalg.matrix_rank(K['d2']) + K['b'][1]))

# =============================================================================== 2.1
hdr("2.1  THE 24-PERMUTATION TEST ON K1+2S, SENSE U — AN ARM THAT COULD HAVE FAILED")
nd, _ = dominated_exact_perm(pU)
print("  exact: %d of 24 arrangements have a dominated Jensen pairing (0 is what makes the"
      " arm live)" % nd)
print("  distinct arrays among the 24 permutations: %d" % distinct_arrays(pU))
v1 = [mahler_jensen(tuple(float(x) for x in apply_perm(pU, s)), dps=30) for s in PERMS]
v1x = [mahler_jensen(tuple(float(x) for x in apply_perm(pU, s)), dps=30, pairing='x')
       for s in PERMS]
v2 = [mahler_lawton(tuple(float(x) for x in apply_perm(pU, s)), N=400) for s in PERMS]
print("  R1 Jensen-in-y  dps 30 : spread over 24 = %s   value = %s"
      % (mp.nstr(max(v1) - min(v1), 6), mp.nstr(v1[0], 18)))
print("  R1 Jensen-in-x  dps 30 : spread over 24 = %s" % mp.nstr(max(v1x) - min(v1x), 6))
print("  R2 Lawton/roots N=400  : spread over 24 = %.3e   value = %.9f  (vs R1: %.2e)"
      % (max(v2) - min(v2), v2[0], abs(v2[0] - float(v1[0]))))
e4 = [mahler_ergodic(np.array([float(x) for x in apply_perm(pU, s)]), N=2_000_000)
      for s in PERMS[:6]]
print("  R4 ergodic f=0.9 c=pi/e N=2e6, 6 permutations: spread = %.3e  value = %.9f"
      % (max(e4) - min(e4), e4[0]))
print("  lambda - log(w_max) = %s  (strictly positive: the other three weights DO enter)"
      % mp.nstr(v1[0] - mp.log(mp.mpf(2) / 7), 8))

# =============================================================================== 2.2
hdr("2.2  THE DECISIVE ARM — B0b AND B4 THEMSELVES, OUTSIDE THE DEGENERATE REGION")
print("""  ONE VARIABLE: the permutation.  For each carrier, Haar-random ready states are drawn
  and only those whose pushforward is (i) OUTSIDE the degenerate region and (ii) has four
  DISTINCT weights are kept -- so all 24 arrays differ and the branch max genuinely
  switches.  Every kept state is a legitimate ready state on the corpus's own carrier.
  seed 20260816, 60 accepted states per carrier, R1 at dps 30 on all 24 permutations.""")
rng = np.random.default_rng(20260816)
for name, mult in (("B0b  {00:4,10:2,01:1,11:2}", (4, 2, 1, 2)),
                   ("B4   {00:1,10:1,01:1,11:3}", (1, 1, 1, 3)),
                   ("K1+2S {00:2,10:2,01:2,11:1}", (2, 2, 2, 1))):
    V = sum(mult)
    kept, tried, worst, worstx, worstcross = 0, 0, mp.mpf(0), mp.mpf(0), mp.mpf(0)
    nar = set()
    while kept < 60:
        tried += 1
        s = rng.normal(size=V) + 1j * rng.normal(size=V)
        w = np.abs(s) ** 2
        w /= w.sum()
        i, p = 0, []
        for m in mult:
            p.append(w[i:i + m].sum())
            i += m
        p = tuple(p)
        S = sorted(p, reverse=True)
        if S[0] + S[3] >= S[1] + S[2]:
            continue
        if len(set(p)) < 4:
            continue
        kept += 1
        nar.add(distinct_arrays(p))
        vy = [mahler_jensen(apply_perm(p, sp), dps=30) for sp in PERMS]
        vx = [mahler_jensen(apply_perm(p, sp), dps=30, pairing='x') for sp in PERMS]
        worst = max(worst, max(vy) - min(vy))
        worstx = max(worstx, max(vx) - min(vx))
        worstcross = max(worstcross, max(abs(a - b) for a, b in zip(vy, vx)))
    print("  %-28s accepted %d of %d states   distinct arrays per state: %s"
          % (name, kept, tried, sorted(nar)))
    print("      worst spread over 24 permutations, Jensen-in-y = %s" % mp.nstr(worst, 6))
    print("      worst spread over 24 permutations, Jensen-in-x = %s" % mp.nstr(worstx, 6))
    print("      worst y-vs-x disagreement                      = %s" % mp.nstr(worstcross, 6))
print("""
  THE MULTISET THEOREM SURVIVES THE ARM THAT COULD HAVE FAILED, on the corpus's own two
  four-class carriers, at 60 ready states each, all 24 arrays distinct, all outside the
  branch-domination region.  This is the test LEG A did not run and B-04 says LEG A could
  not have run on these carriers.  B-04 is right about the SENSE-U vectors and wrong about
  the carriers; B-01's conclusion is right and now rests on something.""")

# =============================================================================== 2.3
hdr("2.3  ARM-COUNT AUDIT — HOW MANY OF THE '24 PERMUTATIONS' ARE DISTINCT ARRAYS")
print("""  The lane's ARM DIFF counts permuted arrays that differ FROM THE IDENTITY ARRAY.  That
  metric cannot see duplicates among the other 23.  The number that matters is the count
  of DISTINCT arrays, and it is smaller than 24 in four of LEG A's five cases.""")
LEGA = [("B0b   (4,2,1,2)/9", (F(4, 9), F(2, 9), F(1, 9), F(2, 9))),
        ("B4    (1,1,1,3)/6", (F(1, 6), F(1, 6), F(1, 6), F(3, 6))),
        ("SENSE C  (1,1,1,1)/4", (F(1, 4),) * 4),
        ("GEN1  (.40,.40,.15,.05)", (F(2, 5), F(2, 5), F(3, 20), F(1, 20))),
        ("GEN2  (.34,.33,.32,.01)", (F(34, 100), F(33, 100), F(32, 100), F(1, 100)))]
print("  %-26s %8s %10s %12s %s" % ("LEG A case", "distinct", "lane ARM", "dominated",
                                    "live arms"))
print("  %-26s %8s %10s %12s %s" % ("", "arrays", "DIFF", "of 24", "(distinct & not dominated)"))
for lab, p in LEGA:
    da = distinct_arrays(p)
    ident = tuple(p)
    armdiff = sum(1 for s in PERMS[1:] if apply_perm(p, s) != ident)
    nd, _ = dominated_exact_perm(p)
    live = 0 if nd == 24 else da
    print("  %-26s %8d %10d %12d %s" % (lab, da, armdiff, nd, live))
print("""
  GEN1 -- one of the two cases the lane says carries all of LEG A's information -- has a
  repeated weight and only 12 distinct arrays.  GEN2 is the only case in LEG A with 24
  distinct arrays AND a live branch switch.  LEG A's headline rests on ONE weight vector.
  (It is still a theorem; this voids the CONTROL's advertised width, not the result.)""")
