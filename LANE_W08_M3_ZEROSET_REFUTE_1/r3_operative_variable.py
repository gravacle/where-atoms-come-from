#!/usr/bin/env python3
# LANE W08 / M3 REFUTER 1 — script 3.  IS THE OPERATIVE VARIABLE NAMED CORRECTLY?
#
# The lane names it:  "THE NUMBER OF FREE RELATIVE PHASES AMONG THE UNIT-MODULUS
# COEFFICIENTS -- equivalently, whether at least one class weight vanishes."
# and its HEADLINE restates that as:
#   "that equivalence is an artefact of p00 = 0, it fails on ANY CARRIER WITH A
#    SPECTATOR VERTEX"
# and its self-flag as:
#   "it is a latent trap that will break the first p00 > 0 carrier anyone reasons
#    about (and the corpus already owns THREE such carriers)".
#
# THREE VARIABLES ARE BEING RUN TOGETHER HERE:
#   (V1) how many relative phases among the OCCUPIED characters are free;
#   (V2) whether at least one of the four class weights is zero;
#   (V3) whether p00 = 0 / whether the carrier has a spectator vertex.
# On K1 alone all three coincide.  This script separates them, twice: once on the
# CORPUS'S OWN carriers (S4 sec 4.2), and once on a carrier with three loops.
# EXACT integer arithmetic for every predicate; FLOAT only where labelled.
import numpy as np
from itertools import product

L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 100)
out("R3  THE OPERATIVE VARIABLE: p00 = 0 / 'SPECTATOR VERTEX' IS NOT IT, AND THE CORPUS")
out("    ALREADY CONTAINS THE COUNTEREXAMPLE")
out("=" * 100)
out()

def D_of(p):
    a, b, c, d = p                       # (p00, p10, p01, p11)
    return (a + d - b - c) * (a + c - b - d) * (a + b - c - d)
def torus_zero(p):
    return D_of(p) <= 0
def polygon(p):
    return 2 * max(p) <= sum(p)

# ------------------------------------------------------------------ (a) the corpus's carriers
out("(a) THE TEN S4 CARRIERS, THEIR CLASS COUNTS AS S4 PUBLISHES THEM (S4:575-585), AND")
out("    WHICH READING IS CORRECT ON EACH.  Weights are the SENSE-U (uniform on vertices)")
out("    states S4 itself uses.  All predicates EXACT integer arithmetic on the counts.")
out()
out("    carrier                       class counts (00,10,01,11)  p00>0  #occupied  torus-zero  polygon  AGREE?")
CARRIERS = [
    ("B0a ring torus, loops disjoint", (2, 4, 3, 0)),      # {00:2, 01:3, 10:4}
    ("B0b ring torus, loops meet   ",  (4, 2, 1, 2)),      # {00:4, 01:1, 10:2, 11:2}
    ("B3  horn torus               ",  (0, 2, 2, 1)),      # {01:2, 10:2, 11:1}
    ("B1  K1 (as handed)           ",  (0, 2, 2, 1)),
    ("B4  spindle                  ",  (1, 1, 1, 3)),      # {00:1, 01:1, 10:1, 11:3}
    ("B2  K1 both filled           ",  (0, 2, 2, 1)),
    ("B1p K1-bridged               ",  (0, 3, 3, 0)),      # {01:3, 10:3}
    ("B1q K1-bridged + SPECTATOR   ",  (1, 3, 3, 0)),      # {00:1, 01:3, 10:3}
    ("B1s K1 subdivided            ",  (0, 5, 5, 1)),      # {01:5, 10:5, 11:1}
]
trap_bites = []
for name, c in CARRIERS:
    occ = sum(1 for x in c if x > 0)
    tz, pg = torus_zero(c), polygon(c)
    agree = (tz == pg)
    if not agree:
        trap_bites.append(name.strip())
    out("    %-30s %-26s  %-5s  %d          %-10s  %-7s  %s"
        % (name, str(c), "YES" if c[0] > 0 else "no", occ, tz, pg, "yes" if agree else "*** NO ***"))
out()
out("    CARRIERS WITH p00 > 0 (a spectator vertex, or loops that miss some vertex): %s"
    % ", ".join(n.strip() for n, c in CARRIERS if c[0] > 0))
out("    CARRIERS WHERE THE POLYGON READING IS ACTUALLY WRONG:                       %s"
    % ", ".join(trap_bites))
out()
out("    ==> THE HEADLINE SENTENCE 'it fails on ANY CARRIER WITH A SPECTATOR VERTEX' IS FALSE,")
out("        AND THE CORPUS'S OWN CARRIER NAMED 'K1-bridged + SPECTATOR VERTEX' IS THE")
out("        COUNTEREXAMPLE.  B1q has p00 = 1/7 > 0 and the polygon/triangle reading is")
out("        EXACTLY RIGHT on it, because its FOURTH class is empty (p11 = 0: no vertex of")
out("        B1q lies on both loops).  S4:592 even computes its rate by the three-term")
out("        Cassaigne-Maillot formula m(1/7 + (3/7)x + (3/7)y), which is the three-class")
out("        object.  B0a is a second such carrier: p00 = 2/9 > 0, p11 = 0, reading correct.")
out("    ==> THE SELF-FLAG'S 'the corpus already owns THREE such carriers' IS WRONG TWICE:")
out("        there are FOUR carriers with p00 > 0 (B0a, B0b, B4, B1q -- the lane lists three")
out("        and omits B0a), and only TWO of them (B0b, B4) are carriers on which the")
out("        polygon reading actually fails.  The lane counted the wrong property.")
out()

# ------------------------------------------------------------------ (b) V2 vs V3, exactly
out("(b) SEPARATING (V2) 'some class weight is zero' FROM (V3) 'p00 = 0'.  EXACT lattice sweep:")
out("    for each of the four slots, force that slot to 0 and ask whether the polygon reading")
out("    equals the torus-zero criterion on the whole remaining face.")
for slot in range(4):
    N = 90
    bad = tot = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            k = N - i - j
            w = [i, j, k]
            p = w[:slot] + [0] + w[slot:]
            p = tuple(p)
            tot += 1
            if torus_zero(p) != polygon(p):
                bad += 1
    out("    slot %d forced to 0 : %5d exact lattice states, #{torus-zero != polygon} = %d"
        % (slot, tot, bad))
out("    => (V2) is the right variable INSIDE the two-loop / four-class setting, and (V3) is")
out("       not: slot 3 (p11 = 0) works exactly as well as slot 0 (p00 = 0), and that is the")
out("       case B0a and B1q actually occupy.  The lane's own C5 tested this and got 0")
out("       mismatches -- and then its headline and its self-flag both reverted to (V3).")
out()

# ------------------------------------------------------------------ (c) V1 vs V2: three loops
out("(c) SEPARATING (V1) 'free relative phases' FROM (V2) 'some class weight vanishes'.")
out("    Two loops can never do this: with four classes, 'a weight vanishes' and 'the")
out("    occupied characters have all relative phases free' are the same condition.  THREE")
out("    loops separate them, and the separation is in BOTH directions.")
out()
out("    A THREE-LOOP CARRIER, built the way S1 builds K1.  Root v0; three triangles")
out("      T_i = v0 -> a_i -> b_i -> v0, i = 1,2,3; T_1 FILLED (a face), T_2, T_3 unfilled;")
out("      plus one spectator vertex w joined by an edge to v0.")
out("      V = 1 + 6 + 1 = 8, E = 9 + 1 = 10, F = 1 -> chi = 8 - 10 + 1 = -1, b1 = 2 ... ")
out("      designated loops gamma_1, gamma_2, gamma_3 = the three triangles; their holonomies")
out("      u_1, u_2, u_3 are independently assignable (gamma_1 bounds -> a curvature; the")
out("      other two are flat).  Class of a vertex = which loops contain it.")
V = ["v0", "a1", "b1", "a2", "b2", "a3", "b3", "w"]
LOOPV = {1: {"v0", "a1", "b1"}, 2: {"v0", "a2", "b2"}, 3: {"v0", "a3", "b3"}}
cls = {}
for v in V:
    cls[v] = tuple(1 if v in LOOPV[i] else 0 for i in (1, 2, 3))
out("      vertex classes: %s" % ", ".join("%s:%s" % (v, "".join(map(str, cls[v]))) for v in V))
out()
out("    CASE 1 -- occupy the SINGLETON classes {100},{010},{001} and the EMPTY class {000}.")
out("      Put the weight on a1, a2, a3 and the spectator w; v0 gets none.")
out("      Occupied characters: {1, u_1, u_2, u_3}.  FOUR of them.  FOUR class weights are")
out("      non-zero and FOUR of the eight classes are empty, so (V2) 'a class weight")
out("      vanishes' HOLDS -- indeed four do.  All three relative phases u_1,u_2,u_3 are")
out("      FREE, so (V1) says the POLYGON inequality is the exact criterion.")
out("    CASE 2 -- occupy {000},{100},{010},{110}: weight on w, a1, a2 and a vertex in both")
out("      loop 1 and loop 2.  Occupied characters {1, u_1, u_2, u_1 u_2}: the third ratio is")
out("      the PRODUCT of the first two, so (V1) says the D-criterion, NOT the polygon one.")
out("      (V2) holds here too -- five of the eight classes are empty.")
out("    So (V2) is TRUE in both cases and the criterion DIFFERS between them.  Checked:")
out()
out("    CASE 1, the classical fact, VERIFIED not assumed.  With all of x,y,z free on T the")
out("    four terms w0, w1 x, w2 y, w3 z have three free relative phases, so")
out("        min_{T^3} |P| = max(0, 2 max(w) - sum(w))       (polygon closure)")
out("    and 'a zero exists' IS the polygon inequality.  Full 3-D brute force, FLOAT:")
rng = np.random.default_rng(31415926)
worst_gap = 0.0
for trial in range(12):
    w = rng.dirichlet([1, 1, 1, 1])
    n = 90
    t = 2 * np.pi * np.arange(n) / n
    e = np.exp(1j * t)
    Pv = (w[0] + w[1] * e[:, None, None] + w[2] * e[None, :, None] + w[3] * e[None, None, :])
    gm = float(np.abs(Pv).min())
    cf = max(0.0, 2 * w.max() - w.sum())
    worst_gap = max(worst_gap, abs(gm - cf))
    if trial < 4:
        out("      w = (%.4f,%.4f,%.4f,%.4f)  90^3 grid min = %.6f  closed form = %.6f  polygon %s"
            % (w[0], w[1], w[2], w[3], gm, cf, 2 * w.max() <= w.sum()))
out("      max over 12 weights of |90^3 grid min - max(0,2max-sum)| = %.4e  (grid error O(1/n))"
    % worst_gap)
out()
out("    CASE 2 uses the SAME weights but the constrained character set {1,u1,u2,u1u2}; its")
out("    criterion is D <= 0, established in r2 by two independent groupings.  Comparing the")
out("    TWO PREDICATES EXACTLY (no grid, no tolerance) on 200000 Dirichlet draws:")
rng2 = np.random.default_rng(27182818)
W = rng2.dirichlet([1, 1, 1, 1], size=200000)
pg = 2 * W.max(axis=1) <= W.sum(axis=1)
a, b, c, d = W[:, 0], W[:, 1], W[:, 2], W[:, 3]
Dv = (a + d - b - c) * (a + c - b - d) * (a + b - c - d)
tz = Dv <= 0
out("      P(CASE 1 has a zero) = P(polygon) = %.5f    (exact 1/2)" % pg.mean())
out("      P(CASE 2 has a zero) = P(D <= 0)  = %.5f    (exact 1/4)" % tz.mean())
out("      #{CASE 1 fires, CASE 2 does not} / N = %.5f   (exact 1/4)" % float((pg & ~tz).mean()))
out("      #{CASE 2 fires, CASE 1 does not} / N = %.5f   (exact 0 -- the containment)"
    % float((tz & ~pg).mean()))
out()
out("    ONE EXHIBITED WITNESS, with a CERTIFIED positive lower bound on the CASE 2 side")
out("    (COR-E discipline: a grid min is an upper bound and cannot certify 'no zero'; the")
out("    closed form min_y |P| = | |w0+w1 x| - |w2+w3 x| | minimised exactly over x can):")
wit = (0.40, 0.24, 0.20, 0.16)
n = 200000
t = 2 * np.pi * (np.arange(n) + 0.5) / n
x = np.exp(1j * t)
g = np.abs(np.abs(wit[0] + wit[1] * x) - np.abs(wit[2] + wit[3] * x))
A = wit[0]**2 + wit[1]**2 - wit[2]**2 - wit[3]**2
B = 2 * (wit[0] * wit[1] - wit[2] * wit[3])
out("      w = %s : polygon 2max-sum = %.6f <= 0 so CASE 1 HAS a zero." % (str(wit), 2*max(wit)-sum(wit)))
out("        CASE 2: A = %.6f, B = %.6f, |A| - |B| = %.6f > 0 so A + B cos t NEVER vanishes;"
    % (A, B, abs(A) - abs(B)))
out("        D = %.6f > 0.  min over a %d-point x-circle of the exact y-minimised |P| = %.6f."
    % ((wit[0]+wit[3]-wit[1]-wit[2])*(wit[0]+wit[2]-wit[1]-wit[3])*(wit[0]+wit[1]-wit[2]-wit[3]),
       n, float(g.min())))
out("        The two readings differ on this single weight vector by 'a zero exists' vs")
out("        'no zero exists'.  Same multiset, same four weights, same 'a class weight")
out("        vanishes' status -- only the CHARACTER LATTICE moved.")
out()
out("    ==> THE 'equivalently, whether at least one class weight vanishes' IN THE LANE'S")
out("        operative_variable FIELD IS A TWO-LOOP-ONLY EQUIVALENCE, STATED UNSCOPED, and")
out("        the lane's headline and self-flag then replace it with the WEAKER AND FALSE")
out("        'p00 = 0 / spectator vertex' form.  The FIRST half of the lane's naming -- the")
out("        number of free relative phases -- is the correct invariant and SURVIVES.")
out("        Its exact form: the occupied characters chi_1..chi_n, as elements of the")
out("        character lattice Z^r of the loop holonomies, have free relative phases iff")
out("        chi_2-chi_1, ..., chi_n-chi_1 are Z-linearly independent AND span a direct")
out("        summand of Z^r.  Exact integer check on every case in play:")
def idet(M):
    """exact integer determinant by fraction-free Gaussian elimination"""
    from fractions import Fraction as F
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    det = F(1)
    for i in range(n):
        piv = None
        for r in range(i, n):
            if A[r][i] != 0:
                piv = r; break
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]; det = -det
        det *= A[i][i]
        inv = A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / inv
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return int(det)

def free_phases(chars):
    import itertools
    from math import gcd
    base = chars[0]
    M = [[c[k] - base[k] for k in range(len(base))] for c in chars[1:]]
    rows, cols = len(M), len(base)
    if rows == 0:
        return True
    if rows > cols:
        return False
    g = 0
    for combo in itertools.combinations(range(cols), rows):
        sub = [[row[j] for j in combo] for row in M]
        g = gcd(g, abs(idet(sub)))
    return g == 1

for label, chars in [("K1, 3 occupied  {u,v,uv}      ", [(1,0),(0,1),(1,1)]),
                     ("4 occupied      {1,u,v,uv}    ", [(0,0),(1,0),(0,1),(1,1)]),
                     ("3 loops CASE 1  {1,u1,u2,u3}  ", [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]),
                     ("3 loops CASE 2  {1,u1,u2,u1u2}", [(0,0,0),(1,0,0),(0,1,0),(1,1,0)]),
                     ("p10 = 0         {1,v,uv}      ", [(0,0),(0,1),(1,1)]),
                     ("p01 = 0         {1,u,uv}      ", [(0,0),(1,0),(1,1)]),
                     ("p11 = 0         {1,u,v}       ", [(0,0),(1,0),(0,1)])]:
    out("          %s  free relative phases: %s" % (label, free_phases(chars)))
out("        -- reproduces every case the lane tested (all four single-slot deletions), plus")
out("        the two three-loop cases its 'equivalently' clause gets wrong.")
out()
out("DONE.")
open("r3_operative_variable.OUT.txt", "w").write("\n".join(L) + "\n")
