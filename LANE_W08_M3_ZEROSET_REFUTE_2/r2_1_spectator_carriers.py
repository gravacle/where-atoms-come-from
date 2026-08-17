#!/usr/bin/env python3
# LANE W08 / M3 — REFUTER 2 — LENS: SCOPE.
# script 1: BUILD ACTUAL CARRIERS WITH A VERTEX IN NEITHER LOOP.
#
# The lane under attack reasons about "the four-class case K1 does not contain" entirely at the
# level of abstract weight vectors p = (p00,p10,p01,p11) in the 3-simplex.  IT NEVER BUILDS A
# CARRIER.  This script builds two, checks their topology and their gauge-invariant content from
# d1/d2 by integer rank, and runs the literal matrix action on C^V.
#
#   K1S  = K1 + a pendant vertex v5 attached to the root by e7 = (v0 -> v5).
#          V=6 E=7 F=1.  FACE_V={0,1,2}, CYC_V={0,3,4}, v5 in NEITHER loop.
#          ALL FOUR CLASSES OCCUPIED.
#   BARB = "barbell": filled triangle {a,b,c} + unfilled triangle {d,e,g}, joined by a path
#          a -> x -> d.  x is in NEITHER loop; NO vertex is in BOTH.
#          CLASSES 00,10,01 OCCUPIED; 11 EMPTY BY INCIDENCE.
#
# PRECISION: IEEE double for the matrix action (identity checks only, deviations reported);
# EXACT integer arithmetic (Fraction / int) for every predicate that decides an inclusion.
import numpy as np
from fractions import Fraction
from itertools import permutations

rng = np.random.default_rng(20260816)
L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 100)
out("R2-1  CARRIERS WITH A VERTEX IN NEITHER LOOP — BUILT, NOT ASSUMED")
out("=" * 100)
out("numpy %s ; IEEE double for the matrix action, EXACT integer arithmetic for every predicate."
    % np.__version__)
out()


# ------------------------------------------------------------------ generic carrier machinery
def topology(V, edges, faces):
    """d1 (V x E), d2 (E x F).  Integer ranks via exact Fraction Gaussian elimination."""
    E = len(edges)
    F = len(faces)
    d1 = [[0] * E for _ in range(V)]
    for ei, (s, t) in enumerate(edges):
        d1[s][ei] -= 1
        d1[t][ei] += 1
    d2 = [[0] * F for _ in range(E)]
    for fi, cyc in enumerate(faces):
        for ei in cyc:
            d2[ei][fi] += 1

    def rank(M):
        M = [[Fraction(x) for x in row] for row in M]
        rows, cols = len(M), (len(M[0]) if M else 0)
        r = 0
        for c in range(cols):
            piv = None
            for i in range(r, rows):
                if M[i][c] != 0:
                    piv = i
                    break
            if piv is None:
                continue
            M[r], M[piv] = M[piv], M[r]
            pv = M[r][c]
            M[r] = [x / pv for x in M[r]]
            for i in range(rows):
                if i != r and M[i][c] != 0:
                    f = M[i][c]
                    M[i] = [a - f * b for a, b in zip(M[i], M[r])]
            r += 1
            if r == rows:
                break
        return r

    r1, r2 = rank(d1), rank(d2)
    b0 = V - r1
    b1 = (E - r1) - r2
    b2 = F - r2
    return dict(V=V, E=E, F=F, chi=V - E + F, rank_d1=r1, rank_d2=r2, b0=b0, b1=b1, b2=b2,
                gauge_invariants=E - (V - 1))


def classes(V, FACE_V, CYC_V):
    cls = []
    for v in range(V):
        cls.append((1 if v in FACE_V else 0, 1 if v in CYC_V else 0))
    return cls


def push(pv, cls):
    d = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
    for w, c in zip(pv, cls):
        d[c] += w
    return (d[(0, 0)], d[(1, 0)], d[(0, 1)], d[(1, 1)])


def Z_direct(s, cls, f, c, k):
    """<M_dF^k s, M_c^k s> by LITERAL matrix action on C^V.  <z,w> = conj(z) w  (S1 sec3)."""
    WF, WC = np.exp(1j * f), np.exp(1j * c)
    a = np.array(s, dtype=complex)
    b = np.array(s, dtype=complex)
    for v, (inF, inC) in enumerate(cls):
        if inF:
            a[v] *= WF ** k
        if inC:
            b[v] *= WC ** k
    return np.vdot(a, b)


def P_poly(p, x, y):
    p00, p10, p01, p11 = p
    return p00 + p10 * x + p01 * y + p11 * x * y


# ------------------------------------------------------------------ the two carriers
# K1 itself, for reference
K1 = dict(name="K1 (S1 sec1)", V=5,
          edges=[(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)],
          faces=[[0, 1, 2]], FACE_V={0, 1, 2}, CYC_V={0, 3, 4})

# K1S: K1 plus ONE pendant vertex on the root.  e7 = (v0 -> v5).
K1S = dict(name="K1S = K1 + pendant v5 on the root", V=6,
           edges=[(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0), (0, 5)],
           faces=[[0, 1, 2]], FACE_V={0, 1, 2}, CYC_V={0, 3, 4})

# BARBELL: v0,v1,v2 filled triangle ; v4,v5,v6 unfilled triangle ; path v0 -> v3 -> v4.
BARB = dict(name="BARB = filled triangle --x-- unfilled triangle", V=7,
            edges=[(0, 1), (1, 2), (2, 0),          # e1 e2 e3  filled triangle
                   (4, 5), (5, 6), (6, 4),          # e4 e5 e6  unfilled triangle
                   (0, 3), (3, 4)],                 # e7 e8     the connecting path
            faces=[[0, 1, 2]], FACE_V={0, 1, 2}, CYC_V={4, 5, 6})

for K in (K1, K1S, BARB):
    t = topology(K["V"], K["edges"], K["faces"])
    cls = classes(K["V"], K["FACE_V"], K["CYC_V"])
    occ = sorted(set(cls))
    K["cls"] = cls
    K["topo"] = t
    out("CARRIER %s" % K["name"])
    out("   V=%d E=%d F=%d  chi=%d  rank d1=%d rank d2=%d  b0=%d b1=%d b2=%d"
        % (t["V"], t["E"], t["F"], t["chi"], t["rank_d1"], t["rank_d2"], t["b0"], t["b1"], t["b2"]))
    out("   gauge-invariant real parameters E-(V-1) = %d   (W_F and W_C exactly, as on K1)"
        % t["gauge_invariants"])
    out("   vertex classes (inF,inC) = %s" % (cls,))
    out("   OCCUPIED CLASSES: %s      EMPTY BY INCIDENCE: %s"
        % (occ, sorted({(0, 0), (1, 0), (0, 1), (1, 1)} - set(cls))))
    out()

out("*** THE FIRST SCOPE FACT, AND IT IS NOT TOPOLOGICAL. ***")
out("K1 and K1S have IDENTICAL chi, b0, b1, b2, one face, one independent cycle, and the SAME")
out("two gauge invariants W_F, W_C.  They differ only in that K1S has one vertex on neither loop.")
out("So 'p00 = 0' is NOT a topological property of the carrier class the corpus works in -- it is")
out("an accident of K1's MINIMALITY (S1 sec2: 'the smallest complex carrying one face and one")
out("independent cycle').  Every K1 result whose proof uses p00 = 0 is a statement about a")
out("codimension-1 FACE of the four-class simplex, i.e. a measure-zero family, on a carrier that")
out("is not distinguished from K1S by any invariant the corpus computes.")
out()

# ------------------------------------------------------------------ (a) identity on K1S / BARB
out("(a) DOES THE IDENTITY  Z_k = P(u^k,v^k)  SURVIVE A SPECTATOR VERTEX?  (M3-F1's scope)")
worst = {}
for K in (K1S, BARB):
    w = 0.0
    for _ in range(400):
        s = rng.normal(size=K["V"]) + 1j * rng.normal(size=K["V"])
        s /= np.linalg.norm(s)
        f, c = rng.uniform(0, 2 * np.pi, 2)
        p = push(np.abs(s) ** 2, K["cls"])
        u, v = np.exp(-1j * f), np.exp(1j * c)
        for k in range(1, 13):
            w = max(w, abs(Z_direct(s, K["cls"], f, c, k) - P_poly(p, u ** k, v ** k)))
    worst[K["name"]] = w
    out("    %-46s 400 states x k=1..12, max dev = %.3e" % (K["name"][:46], w))
out("    => M3-F1 SURVIVES.  The identity is a statement about the class map, not about K1.")
out("       (It could have failed: the class map on K1S sends a vertex to the trivial character,")
out("        which is the one character the K1 derivation never had to carry.)")
out()

# ------------------------------------------------------------------ (b) the headline sentence
out("(b) THE HEADLINE SENTENCE UNDER TEST.  The lane's headline reads:")
out('      "...but that equivalence is an artefact of p00 = 0, IT FAILS ON ANY CARRIER WITH A')
out('       SPECTATOR VERTEX..."')
out("    BARB HAS A SPECTATOR VERTEX (v3, in neither loop) AND THE EQUIVALENCE HOLDS ON IT AT")
out("    EVERY STATE, because BARB has NO vertex in both loops: p11 = 0 identically, by incidence.")
out("    So on BARB the occupied coefficients are {1, u, v} -- three of them, all ratios free --")
out("    and the state-side criterion is the triangle inequality again.  EXACT CHECK:")


def D_of(p):
    a, b, c, d = p
    return (a + d - b - c) * (a + c - b - d) * (a + b - c - d)


N = 200
mismatch = 0
tot = 0
fire = 0
for i in range(N + 1):
    for j in range(N - i + 1):
        k = N - i - j
        p = (i, j, k, 0)                       # (p00,p10,p01,p11) with p11 = 0  -- BARB's family
        tot += 1
        tri = (i <= j + k) and (j <= i + k) and (k <= i + j)
        zero = D_of(p) <= 0                    # EXACT integer
        fire += zero
        if tri != zero:
            mismatch += 1
out("      BARB's family (p11 = 0), EXACT integer simplex denominator N=%d, %d points:" % (N, tot))
out("        #{triangle predicate != torus-zero predicate} = %d   <-- must be 0" % mismatch)
out("        firing fraction %d/%d = %.6f  (-> 1/4, the medial triangle, exactly as on K1)"
    % (fire, tot, fire / tot))
out()
out("    AND THE CORPUS ALREADY OWNS THIS CARRIER.  S4's control list, line 519:")
out("        'B1q  K1-bridged + spectator vertex   V=7 E=8 ...'")
out("    and S4:582 gives its class multiset: {00:1, 01:3, 10:3}.  A SPECTATOR VERTEX AND NO")
out("    SHARED VERTEX -- exactly BARB.  Its characters are listed on the same line as '1 v u'.")
out("    B1q's own uniform state is p = (1/7, 3/7, 3/7, 0): triangle inequality HOLDS (3/7 <= 1/2)")
out("    and the torus zero EXISTS (sorted w1+w4 = 3/7 <= 4/7 = w2+w3).  The two predicates agree,")
out("    at that state and at every other state of B1q.")
out("    *** AND THE LANE PRINTED IT ITSELF: m3_4_threeclass_vs_fourclass.py:34 / .OUT.txt:13 reads")
out("    *** 'B1q K1-bridged+spectator {00:1,01:3,10:3}.  So the p00 > 0 regime is NOT vacuous'.")
out("    *** The lane cited B1q as evidence that p00 > 0 is realisable, read the multiset that")
out("    *** shows 11 is EMPTY on it, and still wrote 'it fails on any carrier with a spectator")
out("    *** vertex' in the headline.  This is the corpus's signature failure mode -- UNDER-READ,")
out("    *** not under-adversarial -- committed by the lane against its own output file.")
out()
out("    *** CORRECTION TO THE LANE'S HEADLINE.  A spectator vertex is NECESSARY for the break")
out("    *** and NOT SUFFICIENT.  What is needed is ALL FOUR CLASSES NON-EMPTY.  The lane's own")
out("    *** m3_5 (C) proves this ('ANY one of the four slots being zero restores the triangle")
out("    *** inequality') -- so the headline contradicts the lane's own isolation C5.")
out()

# ------------------------------------------------------------------ (c) a REALISED break on K1S
out("(c) THE BREAK, REALISED ON AN ACTUAL CARRIER (not a weight vector).  K1S ready state:")
sv = np.zeros(6)
sv[5] = np.sqrt(0.4)                                   # class 00  (the spectator)
sv[1] = sv[2] = np.sqrt(0.1)                           # class 10  -> 0.2
sv[3] = sv[4] = np.sqrt(0.1)                           # class 01  -> 0.2
sv[0] = np.sqrt(0.2)                                   # class 11  -> 0.2
p = push(sv ** 2, K1S["cls"])
out("      |s_v|^2 = %s   (norm %.12f)" % (np.round(sv ** 2, 6).tolist(), float(np.sum(sv ** 2))))
out("      pushforward (p00,p10,p01,p11) = %s" % (tuple(round(float(x), 6) for x in p),))
out("      polygon / convex-hull reading : max = %.4f <= 1/2  ->  SAYS IT FIRES" % max(p))
out("      D = %+.6f > 0                                     ->  NO ZERO ON T^2, NEVER FIRES"
    % D_of(p))
out("      sorted criterion w1+w4 = %.4f  >  w2+w3 = %.4f     ->  NO ZERO"
    % (sorted(p)[-1] + sorted(p)[0], sorted(p)[1] + sorted(p)[2]))
# certify by direct matrix action over a fine torus grid AND by the exact closed form
n = 2000
t = 2 * np.pi * np.arange(n) / n
x = np.exp(1j * t)
closed_min = float(np.min(np.abs(np.abs(p[0] + p[1] * x) - np.abs(p[2] + p[3] * x))))
best = 1e9
for i in range(0, 400):
    for j in range(0, 400):
        f = 2 * np.pi * i / 400
        c = 2 * np.pi * j / 400
        best = min(best, abs(Z_direct(sv.astype(complex), K1S["cls"], f, c, 1)))
out("      LITERAL 6x6 MATRIX ACTION over a 400x400 (f,c) grid: min |Z_1| = %.6f" % best)
out("      exact partial-minimisation closed form              : min |Z_1| = %.6f" % closed_min)
out("      (a grid minimum is an UPPER bound -- COR-E discipline -- so the grid CANNOT certify")
out("       a zero; it agrees with the closed form, and the closed form is what decides.)")
out("      => on K1S this state can never fire at any connection, while W-01's convex-hull")
out("         sentence read on the state side says it can.  THE BREAK IS REALISED ON A CARRIER.")
out()

# ------------------------------------------------------------------ (d) how big is the break
out("(d) HOW MUCH OF K1S's STATE SPACE.  The lane reports the break as a simplex volume (1/2 vs")
out("    1/4).  That is the volume in the FOUR-CLASS SIMPLEX under a flat Dirichlet(1,1,1,1)")
out("    prior on the CLASS WEIGHTS.  On an actual carrier the natural prior is on the STATE, and")
out("    K1S's class map is 1:2:2:1, not 1:1:1:1 -- the pushforward of a uniform state is NOT")
out("    uniform on the simplex.  Reported because the lane's 1/4 and 1/2 are prior-dependent:")
for lab, alpha in [("Dirichlet(1,1,1,1) on the CLASSES (the lane's prior)", None),
                   ("uniform |s_v|^2 on K1S's SIX vertices, pushed forward", K1S)]:
    hits_poly = hits_zero = 0
    M = 200000
    if alpha is None:
        ps = rng.dirichlet([1, 1, 1, 1], size=M)
    else:
        pv = rng.dirichlet([1] * K1S["V"], size=M)
        ps = np.stack([pv[:, 5], pv[:, 1] + pv[:, 2], pv[:, 3] + pv[:, 4], pv[:, 0]], axis=1)
    w = np.sort(ps, axis=1)[:, ::-1]
    zero = (w[:, 0] + w[:, 3] <= w[:, 1] + w[:, 2])
    poly = (w[:, 0] <= 0.5)
    out("      %-52s  fires %.5f   polygon says %.5f   SPURIOUS %.5f"
        % (lab, zero.mean(), poly.mean(), poly.mean() - zero.mean()))
out("    Both priors are declared; neither is 'the' answer.  The QUALITATIVE break (polygon")
out("    strictly weaker) is prior-free; the numbers 1/4 and 1/2 are not.")
out()
out("DONE.")
open("r2_1_spectator_carriers.OUT.txt", "w").write("\n".join(L) + "\n")
