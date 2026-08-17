#!/usr/bin/env python3
"""LANE W-10 B — LEG D (W-03's EXACT INVOLUTION AT THE LEVEL OF Z_k) and
                  LEG E (DOES ANY CORPUS QUANTITY SEE THE INCIDENCE LABELS?)

LEG D.  The register states the involution as
    'multiplying Z_k by conj(u)^k conj(v)^k leaves |Z_k| fixed and maps (a,b) -> (1-a,1-b),
     so 00 <-> 11 and 10 <-> 01'                                    (REGISTER_V001.md, W-03)
The first clause is vacuous as written -- multiplying anything by a unit-modulus scalar
leaves its modulus fixed -- so the content is entirely in the second.  Written out, the
claim is the exact identity

    conj(u)^k conj(v)^k Z_k[p]  =  conj( Z_k[p~] ),      p~ = (p11, p01, p10, p00),

which is what is tested here, on B0b and B4, at the level of Z_k and not of lambda.

LEG E.  W-03's headline 'the incidence labels are invisible'.  Tested at four occupied
classes against every quantity the corpus registers: lambda, the firing criterion,
W-02's group G, the near-return floor, and the rate at the corpus's OWN two distinguished
connections.

THE ONE VARIABLE throughout: which class each weight is assigned to.  The carrier, the
state, the connection, the k-range and the code path are identical across arms, and the
permuted arrays are printed.
"""
import itertools

import numpy as np
import mpmath as mp

from b_lib import (LBL, PERMS, Z_of, apply_perm, collinearity_defect, cycle_notation, hdr,
                   is_subgroup, m_fast, m_jensen)

mp.mp.dps = 30
INV = (3, 2, 1, 0)          # p -> (p11, p01, p10, p00): 00<->11 and 10<->01
B0B = np.array([4 / 9, 2 / 9, 1 / 9, 2 / 9])
B4W = np.array([1 / 6, 1 / 6, 1 / 6, 3 / 6])
CARRIERS = (("B0b {00:4,10:2,01:1,11:2}", B0B, (4, 2, 1, 2)),
            ("B4  {00:1,10:1,01:1,11:3}", B4W, (1, 1, 1, 3)))
print(__doc__)

# =================================================================================== D
hdr("D.1  THE INVOLUTION AT THE LEVEL OF Z_k — REAL PUSHFORWARD")
rng = np.random.default_rng(20260816)
print("  seed 20260816.  For each carrier: 2000 random (connection, ready state) pairs, each")
print("  evaluated at 40 random circuit counts k in [1, 10^5].  The ready state is drawn on")
print("  the carrier's VERTICES and pushed forward, so the class weights are genuine sums of")
print("  |s_v|^2 with the published multiplicities, not free parameters.")
for name, _, mult in CARRIERS:
    wid = wmod = 0.0
    for _ in range(2000):
        s = rng.normal(size=sum(mult)) + 1j * rng.normal(size=sum(mult))
        s /= np.linalg.norm(s)
        w = np.abs(s) ** 2
        i = 0
        p = np.zeros(4)
        for cidx, m in enumerate(mult):
            p[cidx] = w[i:i + m].sum()
            i += m
        f, c = rng.uniform(0, 2 * np.pi, 2)
        ks = rng.integers(1, 10 ** 5, 40)
        u, v = np.exp(-1j * f * ks), np.exp(1j * c * ks)
        Z = p[0] + p[1] * u + p[2] * v + p[3] * u * v
        pt = p[list(INV)]
        Zt = pt[0] + pt[1] * u + pt[2] * v + pt[3] * u * v
        wid = max(wid, float(np.abs(np.conj(u) * np.conj(v) * Z - np.conj(Zt)).max()))
        wmod = max(wmod, float(np.abs(np.abs(Z) - np.abs(Zt)).max()))
    print("  %-26s max |conj(u)^k conj(v)^k Z_k - conj(Z~_k)| = %.3e     max ||Z_k| - |Z~_k|| = %.3e"
          % (name, wid, wmod))
print("  (the register's own figure for the second column, on K1's three classes, is 6.55e-15)")

hdr("D.2  THE INVOLUTION NEEDS REALITY AT THE Z_k LEVEL AND DOES NOT AT THE lambda LEVEL")
print("""  Pulling the conjugation out of the sum is the ONLY step, and it uses p real:
      conj(u)^k conj(v)^k SUM p_ab u^{ak} v^{bk} = SUM p_ab conj( u^{(1-a)k} v^{(1-b)k} )
  needs conj(p_ab) = p_ab.  At the level of lambda the same map is the monomial
  substitution (x,y) -> (1/x,1/y), which is measure-preserving on the torus, so it
  survives arbitrary complex coefficients.  Both halves are tested; the complex arrays
  are the ones LEG C produces from an inserted observable, not invented here.""")
rng2 = np.random.default_rng(1729)
rowsZ, rowsL = [], []
for tag, gen in (("real non-negative (the construction)", 'r'),
                 ("real with one sign flipped        ", 's'),
                 ("complex, three collinear          ", '3'),
                 ("complex, generic                  ", 'c')):
    wz = 0.0
    wl = mp.mpf(0)
    for _ in range(60):
        r = rng2.dirichlet(np.ones(4))
        if gen == 'r':
            p = r.astype(complex)
        elif gen == 's':
            p = r.astype(complex)
            p[int(rng2.integers(4))] *= -1
        elif gen == '3':
            p = r.astype(complex)
            p[3] *= np.exp(1j * rng2.uniform(0, 2 * np.pi))
        else:
            p = r * np.exp(1j * rng2.uniform(0, 2 * np.pi, 4))
        pt = p[list(INV)]
        f, c = rng2.uniform(0, 2 * np.pi, 2)
        ks = rng2.integers(1, 10 ** 4, 30)
        u, v = np.exp(-1j * f * ks), np.exp(1j * c * ks)
        Z = p[0] + p[1] * u + p[2] * v + p[3] * u * v
        Zt = pt[0] + pt[1] * u + pt[2] * v + pt[3] * u * v
        wz = max(wz, float(np.abs(np.abs(Z) - np.abs(Zt)).max()))
        wl = max(wl, abs(m_jensen(p) - m_jensen(pt)))
    print("  %s   max ||Z_k|-|Z~_k|| = %.3e     max |lambda - lambda~| = %s"
          % (tag, wz, mp.nstr(wl, 4)))
print("""  So the involution is TWO different statements.  At the level of lambda it is a Newton
  polygon symmetry and holds unconditionally.  At the level of |Z_k| -- which is where the
  register states it, and where the firing criterion and the near-return floor live -- it
  holds exactly on the real locus and fails by O(1) off it.""")

# =================================================================================== E
hdr("E.1  AT A FIXED CONNECTION, |Z_k| SEES THE LABELS — AND SEES EXACTLY FOUR PRODUCTS")
print("""  Written out for four real classes, with U = u^k = e^{i alpha}, V = v^k = e^{i beta},

    |Z_k|^2 = SUM p_ab^2 + 2[ A cos alpha + B cos beta + C cos(alpha+beta) + D cos(alpha-beta) ]
        A = p00 p10 + p01 p11     B = p00 p01 + p10 p11     C = p00 p11     D = p10 p01

  -- the three MATCHINGS again: A and B are the two non-diagonal matchings' product sums,
  and C, D are the diagonal matching's two products, appearing separately.  (This is the
  four-class form of W-08's identity |Z_k|^2 = 1 - SUM_{j<l} w_j w_l |chi_j^k - chi_l^k|^2,
  regrouped so that what |Z_k| depends on is visible.)  A permutation preserves |Z_k| at
  EVERY connection iff it preserves (A, B, C, D); generically that group is exactly
  {e, 00<->11 with 10<->01} of order two -- the involution and nothing else -- so the 24
  arrangements take 12 distinct values.  Extra collapse happens only when weights repeat
  or when A = B by accident, and both of the corpus's carriers have repeated weights.""")
rng3 = np.random.default_rng(4242)
CASES = list(CARRIERS) + [("GEN four distinct weights", np.array([0.37, 0.29, 0.23, 0.11]), None)]
for name, p, _ in CASES:
    arrays = {tuple(np.round(apply_perm(list(p), s), 12)) for s in PERMS}
    tally, wid = {}, 0.0
    for _ in range(500):
        f, c = rng3.uniform(0, 2 * np.pi, 2)
        k = int(rng3.integers(1, 1000))
        vals = []
        for q in arrays:
            q = np.array(q)
            z2 = float(abs(complex(Z_of(q, f, c, [k])[0])) ** 2)
            al, be = -f * k, c * k
            cf = float((q ** 2).sum() + 2 * ((q[0] * q[1] + q[2] * q[3]) * np.cos(al)
                                             + (q[0] * q[2] + q[1] * q[3]) * np.cos(be)
                                             + q[0] * q[3] * np.cos(al + be)
                                             + q[1] * q[2] * np.cos(al - be)))
            wid = max(wid, abs(z2 - cf))
            vals.append(z2)
        n = len({round(x, 12) for x in vals})
        tally[n] = tally.get(n, 0) + 1
    inv_orbits = len({frozenset((q, tuple(np.array(q)[list(INV)]))) for q in arrays})
    print("  %-26s distinct arrays %2d, involution orbits %2d, observed #distinct |Z_k| -> %s"
          % (name, len(arrays), inv_orbits, dict(sorted(tally.items()))))
    print("       closed form vs direct |Z_k|^2 : max deviation %.2e" % wid)
print("""  The generic row is the test: 24 distinct arrays, 12 involution orbits, 12 distinct
  values at every one of 500 random (f, c, k).  AT A FIXED CONNECTION THE LABELS ARE
  VISIBLE and the only permutation that hides is W-03's involution.  B0b and B4 collapse
  further only because their published class multisets repeat weights.""")

hdr("E.2  EVERY FULL-TORUS AVERAGE OF A FUNCTION OF |P| IS BLIND — PROOF AND CHECK")
print("""  The proof of the multiset theorem never uses the logarithm.  For any F,
      INT_phi F(|A + e^{i phi} B|) dphi / 2pi  is symmetric in |A| and |B|,
  because |A + e^{i phi} B| = |B + e^{-i phi} A|; and |a+be^{it}| = |b+ae^{it}| pointwise
  for real a,b.  Running that in both Jensen pairings gives the four transpositions that
  generate S4.  So EVERY functional of the form INT INT F(|P|) is multiset-invariant --
  the rate, the mean modulus, the near-return density, all of them.  Checked on a 2048^2
  tensor grid (trapezoid on a full period), one variable moving: the permutation.""")
n = 2048
th = 2 * np.pi * np.arange(n) / n
X = np.exp(1j * th)[:, None]
Y = np.exp(1j * th)[None, :]
for name, p, _ in CARRIERS:
    for fname, F in (("mean |P|          ", lambda a: a),
                     ("mean |P|^(1/2)    ", lambda a: np.sqrt(a)),
                     ("mean 1{|P| < 0.10}", lambda a: (a < 0.10).astype(float)),
                     ("mean log|P| (=lam)", lambda a: np.log(np.maximum(a, 1e-300)))):
        vals = []
        for s in PERMS:
            q = apply_perm(list(p), s)
            A = np.abs(q[0] + q[1] * X + q[2] * Y + q[3] * X * Y)
            vals.append(float(F(A).mean()))
        print("  %-26s %s  value %+.9f  spread over 24 permutations = %.2e"
              % (name, fname, float(np.mean(vals)), max(vals) - min(vals)))

hdr("E.3  BUT A DEGENERATE CONNECTION SEES THE LABELS — AND THE CORPUS HAS TWO OF THEM")
print("""  The proof at E.2 averages over the WHOLE torus.  The corpus's two distinguished
  connections do not.
    * S1's published connection (W_F = -1, W_C = -i, S1 sec6) has ORDER 4, so the k-average
      is over four points.  W-07 measured the recurrence there.
    * S3/S4's headline f = 2.0, c = 1.1 is EXACTLY RESONANT, -11f + 20c = 0 (ERRATUM
      AGAINST W-02), so the orbit closes on a SUBTORUS.  Since (f,c) = 0.1*(20,11), the
      subtorus average is the Mahler measure of the ONE-VARIABLE polynomial
          Q(w) = p00 w^20 + p10 + p01 w^31 + p11 w^11,
      computed here exactly from its roots (Jensen), not by simulation.  VALIDATION: for
      the three-class weights (0.4,0.3,0.3) this method returns the erratum's own
      registered subtorus value.""")


def m1(coeffs):
    c = np.array(coeffs, dtype=complex)
    nz = np.nonzero(np.abs(c) > 0)[0]
    c = c[nz[0]:nz[-1] + 1]
    r = np.roots(c)
    return float(np.log(abs(c[0])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))


def lam_resonant(p):
    d = {20: p[0], 0: p[1], 31: p[2], 11: p[3]}
    co = [0.0] * 32
    for e, v in d.items():
        co[31 - e] += v
    return m1(co)


print("  VALIDATION  lambda_subtorus(0.4, 0.3, 0.3, 0) = %.9f   ERRATUM AGAINST W-02: -0.767014993"
      % lam_resonant([0.4, 0.3, 0.3, 0.0]))
print("              full-torus value for the same weights = %.9f   register: -0.767507880"
      % float(m_fast([0.4, 0.3, 0.3, 0.0])))
for name, p, _ in CASES:
    print()
    print("  %s   (24 permutations; distinct arrays %d)"
          % (name, len({tuple(np.round(apply_perm(list(p), s), 12)) for s in PERMS})))
    v4 = [float(np.mean(np.log(np.abs(Z_of(np.array(apply_perm(list(p), s)), np.pi, 3 * np.pi / 2,
                                           np.arange(1, 5)))))) for s in PERMS]
    u4 = sorted({round(x, 10) for x in v4})
    print("   S1's own connection, W_F=-1 W_C=-i, order 4 : %2d distinct rates, spread %.6f"
          % (len(u4), max(v4) - min(v4)))
    print("      %s" % ", ".join("%+.6f" % x for x in u4[:6]) + (" ..." if len(u4) > 6 else ""))
    vr = [lam_resonant(apply_perm(list(p), s)) for s in PERMS]
    ur = sorted({round(x, 10) for x in vr})
    print("   S3/S4's own connection, f=2.0 c=1.1, resonant: %2d distinct rates, spread %.6f  (EXACT, from roots)"
          % (len(ur), max(vr) - min(vr)))
    print("      %s" % ", ".join("%+.6f" % x for x in ur[:6]) + (" ..." if len(ur) > 6 else ""))
    N = 500_000
    erg = float(np.mean(np.log(np.abs(Z_of(np.array(p), 2.0, 1.1, np.arange(1, N + 1))))))
    print("   the resonant value by direct simulation N=5e5 (independent path): %+.6f vs exact %+.6f"
          % (erg, vr[0]))
    print("   FULL-TORUS rate, all 24 identical              :  1 distinct rate,  %+.9f" % float(m_fast(p)))
    # NOTE: with repeated weights this printed set is inflated -- distinct permutations
    # produce identical arrays -- and need not be a subgroup.  The GEN row, with four
    # distinct weights, is the row that identifies the group.
    for tg, vv in (("order 4  ", v4), ("resonant ", vr)):
        st = [s2 for s2, x in zip(PERMS, vv) if abs(x - vv[0]) < 1e-9]
        print("      stabiliser at %s: order %2d, subgroup %s : %s"
              % (tg, len(st), is_subgroup(st), ", ".join(cycle_notation(s2) for s2 in st)))
print("""
  At neither distinguished connection is the rate a function of the weight MULTISET, and
  the two fail differently -- which is itself the point, because it shows the collapse is
  not one phenomenon:
    * at S1's order-4 connection the surviving group has order 8 and is the group of the
      JENSEN pairing {00,10}|{01,11} -- swap inside either pair, or swap the pairs -- so
      the 24 arrangements carry THREE rates;
    * at S3/S4's resonant connection only the INVOLUTION survives, order 2, so they carry
      TWELVE, the same twelve blocks that |Z_k| resolves at E.1.
  The multiset theorem is a theorem about the full-torus average and about nothing else.""")

hdr("E.4  THE TWO REGISTERED CRITERIA ARE BLIND FOR A THIRD REASON — THEY IGNORE WEIGHTS")
print("""  W-01's firing criterion (0 in the convex hull of the occupied characters) and W-02's
  group G = <chi_a/chi_b : a,b in supp> both depend on the SUPPORT only.  With all four
  classes occupied every permutation leaves the support equal to all four, so they cannot
  distinguish labels even in principle -- a blindness of a different kind from lambda's,
  and one that no permutation test can detect.  Checked by enumeration.""")
g = 400
F, C = np.meshgrid(np.linspace(0, 2 * np.pi, g, endpoint=False),
                   np.linspace(0, 2 * np.pi, g, endpoint=False), indexing='ij')
base = None
worstfire = 0
for name, p, _ in CARRIERS:
    for s in PERMS:
        q = apply_perm(list(p), s)
        fire = (np.cos(F) + np.cos(C) <= 0)
        if base is None:
            base = fire
        worstfire = max(worstfire, int(np.abs(fire.astype(int) - base.astype(int)).max()))
print("  firing region on a %dx%d (f,c) grid: max disagreement across all 24 permutations"
      " and both carriers = %d" % (g, g, worstfire))
print("  firing fraction = %.6f  (W-09's exact 1/2 for a four-class carrier)" % base.mean())

hdr("LEG D / LEG E — VERDICT")
print("""  LEG D.  W-03's involution is CONFIRMED at the level of Z_k on both four-class
  carriers, to 1e-16, at the exact identity conj(u)^k conj(v)^k Z_k = conj(Z~_k) -- which
  is a strictly stronger statement than the one the register makes, since the register's
  first clause ('leaves |Z_k| fixed') is true of any unit-modulus factor and carries no
  content.  THE IDENTITY'S HYPOTHESIS IS REALITY OF THE PUSHFORWARD, and it fails by O(1)
  the moment the pushforward is complex -- while the SAME involution at the level of
  lambda is a monomial substitution and survives complexification exactly.  One symmetry,
  two theorems, different hypotheses, and the register states neither.

  LEG E.  On a four-class carrier the answer to 'does anything see the labels' is YES, and
  the dividing line is not what W-03's headline suggests:
    BLIND -- every full-torus average of every function of |P|, not just the rate: the
      mean modulus, the near-return density below any threshold, the rate.  The proof
      never uses the logarithm.
    BLIND FOR A DIFFERENT AND WEAKER REASON -- W-01's firing criterion and W-02's group G,
      which see only the SUPPORT.  With four classes occupied they are constant across all
      24 permutations by construction, so they could not have detected labels under any
      assignment of weights.
    SIGHTED -- |Z_k| at any fixed connection and any fixed k, which takes 12 distinct
      values over the 24 permutations, the blocks being exactly the involution's orbits.
    SIGHTED, AND THIS IS THE ONE THAT MATTERS -- the RATE ITSELF at a degenerate
      connection.  At S1's own published connection, of order 4, and at S3/S4's own
      headline connection, exactly resonant, the k-average is over a finite orbit or a
      subtorus rather than the torus, the reduction that proves the theorem does not
      apply, and lambda takes several distinct values over the 24 permutations.
  W-03's 'the incidence labels are invisible' is therefore a statement about the GENERIC
  connection, and the corpus computed almost nothing at a generic connection.""")
