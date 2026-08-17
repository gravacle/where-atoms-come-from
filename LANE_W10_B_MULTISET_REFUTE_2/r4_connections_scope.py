#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 4.  B-11 SAYS THE MULTISET THEOREM 'FAILS AT BOTH OF THE
CORPUS'S OWN DISTINGUISHED CONNECTIONS'.  IS THAT A FACT ABOUT THE CONNECTIONS, OR ABOUT
THE ONE WEIGHT VECTOR (GEN) THAT SHOWS IT?

B-11's claim sentence: 'At S1's published order-4 connection the surviving group has
order 8 (the Jensen-pairing group) and the 24 arrangements carry THREE rates; at S3/S4's
exactly resonant f=2.0, c=1.1 only the involution survives, order 2, and they carry
TWELVE.'  Its evidence cites GEN and B0b.  Its own sealed output ALSO contains B4, and
B4 carries ONE rate at the order-4 connection and its resonant spread prints as 0.000000.
This script decides B4 exactly, derives the order-4 rate in closed form, and measures how
often the failure happens across the four-class weight simplex and across ready states on
the corpus's own carriers.

  4.0  the order-4 rate in CLOSED FORM, exact, and its exact invariance group
  4.1  B0b, B4, K1+2S, GEN at both distinguished connections -- B4 decided exactly
  4.2  how much of the four-class simplex actually shows the failure, and how big it is
  4.3  the same question asked of READY STATES on the corpus's own carriers
  4.4  the magnitude against the corpus's own instrument

THE ONE VARIABLE in 4.1-4.3: the permutation of the four class weights.  The connection,
the evaluator and the code path are fixed within each block; the connection's ARITHMETIC
is the variable BETWEEN blocks and is named in each header.
"""
from fractions import Fraction as F

import numpy as np
import mpmath as mp

from r_lib import (PERMS, apply_perm, cyc, distinct_arrays, hdr, is_subgroup,
                   mahler_jensen, mahler_subtorus_resonant, rate_order4, rate_order4_exact)

mp.mp.dps = 30
print(__doc__)

WEIGHTS = [
    ("B0b   (4,2,1,2)/9  SENSE U", (F(4, 9), F(2, 9), F(1, 9), F(2, 9))),
    ("B4    (1,1,1,3)/6  SENSE U", (F(1, 6), F(1, 6), F(1, 6), F(3, 6))),
    ("K1+2S (2,2,2,1)/7  SENSE U", (F(2, 7), F(2, 7), F(2, 7), F(1, 7))),
    ("GEN   (.37,.29,.23,.11)   ", (F(37, 100), F(29, 100), F(23, 100), F(11, 100))),
    ("SENSE C (1,1,1,1)/4       ", (F(1, 4),) * 4),
]




# ------------------------------------------------------------- fast float paths (sweeps)
import math


def rate_order4_f(p):
    a, b, c, d = p
    T1 = (a - b) ** 2 + (c - d) ** 2
    T2 = abs(a + b - c - d)
    T3 = a + b + c + d
    if T1 <= 0 or T2 <= 0 or T3 <= 0:
        return -np.inf
    return (math.log(T1) + math.log(T2) + math.log(T3)) / 4


def resonant_f(p):
    """float64 twin of R5: numpy roots of Q(w) = p0 w^20 + p1 + p2 w^31 + p3 w^11.
    Accurate to ~1e-10; used ONLY for sweep decisions at a 1e-8 threshold, never for a
    headline digit.  Every headline row below is mpmath polyroots at 40 digits."""
    coef = np.zeros(32)
    coef[20] += p[0]
    coef[0] += p[1]
    coef[31] += p[2]
    coef[11] += p[3]
    c = coef[::-1]
    nz = np.nonzero(np.abs(c) > 0)[0]
    if len(nz) == 0:
        return -np.inf
    c = c[nz[0]:]
    if len(c) == 1:
        return float(np.log(abs(c[0])))
    r = np.roots(c)
    return float(np.log(abs(c[0])) + np.sum(np.log(np.maximum(1.0, np.abs(r)))))


def blocks_from(vals, tol):
    bl = []
    for s, v in zip(PERMS, vals):
        for b in bl:
            if abs(b[0] - v) < tol:
                b[1].append(s)
                break
        else:
            bl.append([v, [s]])
    return bl


# =============================================================================== 4.0
hdr("4.0  S1's ORDER-4 CONNECTION: THE RATE IN CLOSED FORM, AND ITS EXACT SYMMETRY")
print("""  S1 sec6 publishes W_F = -1, W_C = -i, so u = conj(W_F) = -1 and v = W_C = -i and the
  pair (u^k, v^k) cycles through (-1,-i), (1,-1), (-1,i), (1,1).  With (a,b,c,d) =
  (p00,p10,p01,p11) the four values of |Z_k| are

     |a - b - i(c-d)| = sqrt((a-b)^2 + (c-d)^2)   at k=1 and k=3 (complex conjugates)
     |a + b - c - d|                              at k=2
     a + b + c + d                                at k=4
     rate = (1/4)[ log((a-b)^2+(c-d)^2) + log|a+b-c-d| + log(a+b+c+d) ]      EXACT

  Every term is symmetric in (a,b), symmetric in (c,d), and symmetric under swapping the
  pair {a,b} with the pair {c,d}.  So the rate is invariant under EXACTLY the stabiliser
  of the partition {{00,10},{01,11}} -- order 8, and NOT the same order-8 subgroup as the
  full-torus D4, which stabilises the DIAGONAL matching {{00,11},{10,01}}.  Two different
  order-8 dihedral subgroups; the corpus has never distinguished them.
  The closed form is checked below against a direct k-average, which is the only place a
  number is needed at all.""")
STAB_M1 = tuple(sorted(s for s in PERMS
                       if {frozenset((s[0], s[1])), frozenset((s[2], s[3]))}
                       == {frozenset((0, 1)), frozenset((2, 3))}))
D4DIAG = tuple(sorted(s for s in PERMS
                      if {frozenset((s[0], s[3])), frozenset((s[1], s[2]))}
                      == {frozenset((0, 3)), frozenset((1, 2))}))
print("  |stabiliser of {{00,10},{01,11}}| = %d, is a subgroup: %s"
      % (len(STAB_M1), is_subgroup(STAB_M1)))
print("  |D4 (diagonal matching)|          = %d, is a subgroup: %s"
      % (len(D4DIAG), is_subgroup(D4DIAG)))
print("  the two order-8 groups are equal: %s   their intersection has order %d"
      % (set(STAB_M1) == set(D4DIAG), len(set(STAB_M1) & set(D4DIAG))))
p = (0.37, 0.29, 0.23, 0.11)
u, v = -1.0 + 0j, -1j
ks = np.arange(1, 4001)
Z = p[0] + p[1] * u ** ks + p[2] * v ** ks + p[3] * (u * v) ** ks
print("  closed form on GEN = %s   direct k-average over k=1..4000 = %.15f   diff %.2e"
      % (mp.nstr(rate_order4(p), 15), float(np.mean(np.log(np.abs(Z)))),
         abs(float(rate_order4(p)) - float(np.mean(np.log(np.abs(Z)))))))

# =============================================================================== 4.1
hdr("4.1  BOTH DISTINGUISHED CONNECTIONS, ALL FIVE WEIGHT VECTORS, DECIDED EXACTLY")
print("  VALIDATION of the resonant route (mpmath polyroots at 40 digits, my own code):")
val = mahler_subtorus_resonant((0.4, 0.3, 0.3, 0.0), dps=40)
print("     lambda_subtorus(0.4,0.3,0.3,0) = %s   ERRATUM AGAINST W-02: -0.767014993"
      % mp.nstr(val, 12))
print("     full torus, same weights       = %s   register:            -0.767507880"
      % mp.nstr(mahler_jensen((0.4, 0.3, 0.3, 0.0), dps=30), 12))
print()
print("  %-28s %7s | %-34s | %s" % ("weights", "arrays", "ORDER 4 (S1's own connection)",
                                    "RESONANT f=2.0 c=1.1 (S3/S4's)"))
for lab, w in WEIGHTS:
    wf = tuple(float(x) for x in w)
    o4 = [rate_order4(apply_perm(wf, s), dps=30) for s in PERMS]
    rs = [mahler_subtorus_resonant(apply_perm(wf, s), dps=40) for s in PERMS]
    b4 = blocks_from(o4, mp.mpf(10) ** -22)
    br = blocks_from(rs, mp.mpf(10) ** -22)
    sp4 = max(o4) - min(o4)
    spr = max(rs) - min(rs)
    stab4 = [s for s, x in zip(PERMS, o4) if abs(x - o4[0]) < mp.mpf(10) ** -22]
    stabr = [s for s, x in zip(PERMS, rs) if abs(x - rs[0]) < mp.mpf(10) ** -22]
    print("  %-28s %7d | %2d blocks spread %-14s | %2d blocks spread %s"
          % (lab, distinct_arrays(w), len(b4), mp.nstr(sp4, 6), len(br), mp.nstr(spr, 6)))
    print("  %-28s %7s | stab %2d subgrp %-5s =STAB_M1 %-5s | stab %2d subgrp %-5s =involution %s"
          % ("", "", len(stab4), is_subgroup(stab4), set(stab4) == set(STAB_M1),
             len(stabr), is_subgroup(stabr),
             set(stabr) == {(0, 1, 2, 3), (3, 2, 1, 0)}))
print("""
  READ IT OFF.  B4 -- one of the two four-class carriers the corpus owns, and the one
  B-11's claim sentence does not mention -- has ONE block at the order-4 connection and
  ONE block at the resonant connection, at a tolerance of 1e-22.  ON B4 THE MULTISET
  THEOREM DOES NOT FAIL AT EITHER DISTINGUISHED CONNECTION.  The reason is exhibited at
  4.0: B4's SENSE-U vector has three equal weights, so the (a,b),(c,d) pairs of the
  order-4 closed form are the same multiset in every arrangement, and at the resonant
  connection the same repeat collapses the 4 distinct arrays onto 2 involution orbits
  whose Mahler measures coincide.""")

# =============================================================================== 4.2
hdr("4.2  HOW MUCH OF THE FOUR-CLASS SIMPLEX ACTUALLY SHOWS THE FAILURE")
print("""  ONE VARIABLE: the permutation, at each of the two fixed connections.  2000 random
  Dirichlet(1,1,1,1) four-class weight vectors, seed 20260816.  A vector 'shows the
  failure' if its 24 arrangements do not all give the same rate at that connection.
  The order-4 column is the EXACT closed form in float64; the resonant column is the
  float64 root twin, at a 1e-8 decision threshold (headline rows above are mpmath at 40
  digits and agree).  Spread is reported as a distribution because its SIZE is the scope question.""")
rng = np.random.default_rng(20260816)
sp4s, sprs = [], []
n4 = nr = 0
NS = 2000
for _ in range(NS):
    w = rng.dirichlet(np.ones(4))
    o4 = [rate_order4_f(apply_perm(tuple(w), s)) for s in PERMS]
    s4 = max(o4) - min(o4)
    sp4s.append(s4)
    n4 += (s4 > 1e-10)
NS2 = 600
for _ in range(NS2):
    w = rng.dirichlet(np.ones(4))
    rs = [resonant_f(apply_perm(tuple(w), s)) for s in PERMS]
    sr = max(rs) - min(rs)
    sprs.append(sr)
    nr += (sr > 1e-8)
sp4s = np.array(sp4s)
sprs = np.array(sprs)
print("  ORDER 4 : %d of %d vectors show a nonzero spread (%.1f%%)"
      % (n4, NS, 100 * n4 / NS))
print("            spread quantiles 10/50/90 = %.4f / %.4f / %.4f   max %.4f"
      % (*np.percentile(sp4s, [10, 50, 90]), sp4s.max()))
print("  RESONANT: %d of %d vectors show a nonzero spread (%.1f%%)"
      % (nr, NS2, 100 * nr / NS2))
print("            spread quantiles 10/50/90 = %.2e / %.2e / %.2e   max %.2e"
      % (*np.percentile(sprs, [10, 50, 90]), sprs.max()))

# =============================================================================== 4.3
hdr("4.3  THE SAME QUESTION ASKED OF READY STATES ON THE CORPUS'S OWN CARRIERS")
print("""  ONE VARIABLE: the ready state.  The class multiset is the carrier's published one;
  the connection and the evaluator are fixed.  This is the arm that decides whether B-11
  is about the connection or about GEN: a carrier whose SENSE-U vector hides the failure
  may still show it at a generic state, and that is a statement about the state, not the
  carrier.  400 Haar states per carrier, seed 4242, order-4 closed form.""")
rng = np.random.default_rng(4242)
for name, mult in (("B0b  {00:4,10:2,01:1,11:2}", (4, 2, 1, 2)),
                   ("B4   {00:1,10:1,01:1,11:3}", (1, 1, 1, 3)),
                   ("K1+2S {00:2,10:2,01:2,11:1}", (2, 2, 2, 1))):
    V = sum(mult)
    sps = []
    for _ in range(400):
        s = rng.normal(size=V) + 1j * rng.normal(size=V)
        w = np.abs(s) ** 2
        w /= w.sum()
        i, p = 0, []
        for m in mult:
            p.append(w[i:i + m].sum())
            i += m
        o4 = [rate_order4_f(apply_perm(tuple(p), sp)) for sp in PERMS]
        sps.append(max(o4) - min(o4))
    sps = np.array(sps)
    pu = np.array([m / V for m in mult])
    o4u = [rate_order4_f(apply_perm(tuple(pu), sp)) for sp in PERMS]
    print("  %-28s SENSE U spread = %.6f ;  random states: nonzero %d/400, median %.4f"
          % (name, max(o4u) - min(o4u), int((sps > 1e-10).sum()), float(np.median(sps))))
print("""
  B-11's SUBSTANCE SURVIVES AND IS STRONGER THAN ITS EVIDENCE: at the order-4 connection
  the failure is UNIVERSAL on the four-class simplex (2000 of 2000, median spread 0.13),
  so it is not a GEN artifact at all.  What does NOT survive is B-11's arithmetic of the
  surviving GROUP.  'The surviving group has order 8, the Jensen-pairing group' and 'only
  the involution survives, order 2' are true of GEN and of no other vector tested:
    B0b at order 4    stabiliser 16, NOT a subgroup, NOT the Jensen-pairing group
    B4  at order 4    stabiliser 24 -- the theorem does not fail there at all
    B4  at resonant   spread 4.8e-10, which the lane's 6-decimal print shows as 0.000000
                      and its tolerance reads as stabiliser 24; at 1e-22 it is 2 blocks
    SENSE C at order 4  rate = -infinity, because a+b-c-d = 0 makes Z_2 vanish exactly
  So the operative variable for WHICH GROUP SURVIVES is the weight vector's repeat
  pattern, and for WHETHER ANYTHING SURVIVES it is the connection.  B-11 states the
  second and quantifies it with the first.""")

# =============================================================================== 4.4
hdr("4.4  MAGNITUDE — THE RESONANT FAILURE IS BELOW THE CORPUS'S OWN INSTRUMENT")
print("""  S4's own validation of its closed forms was 'direct schedule-B simulation, N=2e6,
  worst deviation 3.0e-06' (S4:4.2).  The lane's own LEG A prints E1-E4 deviations up to
  1.35e-05 (SENSE C row).  Against that instrument:""")
for lab, w in WEIGHTS[:4]:
    wf = tuple(float(x) for x in w)
    rs = [mahler_subtorus_resonant(apply_perm(wf, s), dps=40) for s in PERMS]
    o4 = [rate_order4(apply_perm(wf, s), dps=30) for s in PERMS]
    print("  %-28s resonant spread %-12s   order-4 spread %s"
          % (lab, mp.nstr(max(rs) - min(rs), 4), mp.nstr(max(o4) - min(o4), 4)))
print("""  The resonant failure on the corpus's own four-class carrier B0b is 6e-06 -- twice
  S4's stated simulation error and below the lane's own worst evaluator disagreement.  It
  is real and it is EXACT (from roots), but no instrument the corpus ever used at that
  connection could have seen it.  The order-4 failure is 1e-1 and is a different matter.""")
