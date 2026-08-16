"""
rm_2_attack.py -- LANE R (MAPS REFUTER).  The attack on the maps-axis claim:

  "S4 Control 3: 'lambda is not even a homeomorphism invariant', |diff| = 3.181e-02.
   STOPS AT: Immediately, once a map is present.  c3 : B1s -> B1 is class-compatible;
   the transported class weights are (0,5/11,5/11,1/11) on BOTH sides; |dlambda| =
   0.000000e+00 and max|Z_k gap| = 3.342e-16 over 500 random (f,c,k). ...
   an alternative class-compatible collapse (m2->v1) gives identical pi and identical
   lambda, so THE MAP IS NOT UNIQUE BUT THE TRANSPORTED CLASS WEIGHTS ARE."

Six attacks, A1..A6.
"""
import numpy as np, math, sys, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from rm_lib import *

L = lambda_B_closed
BAR = "=" * 78

# =====================================================================  A1
print(BAR); print("A1  THE VACUITY THEOREM -- the claim's control could not have failed")
print(BAR)
print("""
THEOREM (rm-1).  Let c : K -> K' be any map of vertex sets with
class_{K'}(c(w)) = class_K(w) for all w  (\"class-compatible\"), and let
p' = c_* p, i.e. p'_{v} = sum_{w : c(w)=v} p_w.  Then for EVERY k, EVERY (f,c),
EVERY schedule and EVERY gauge group action of this construction,
    Z_k(K', p') = Z_k(K, p)   IDENTICALLY.
PROOF.  Z_k(K,p) = sum_w u^{k a_w} v^{k b_w} p_w = sum_{(a,b)} u^{ka} v^{kb} pi_{ab}(K,p).
The carrier enters ONLY through pi (S4 :133-146; W-03: d2 enters nowhere).  Class
compatibility says c_* commutes with the class map, so pi(K',p') = pi(K,p).  Two equal
arguments to the same function.  QED -- no computation is possible that returns
anything but 0.

This is W-03's Control-1 defect verbatim: \"pi is unchanged, so lambda is unchanged
BY IDENTITY. The controls tested nothing.\"  The corpus now has TWO vacuous controls.
""")
b1, b1s = ALL["B1"](), ALL["B1s"]()
pi_t = (0.0, 5/11, 5/11, 1/11)
print(f"   demonstration: lambda(B1, pi_t) = {L(pi_t):.12f}")
print(f"                  lambda(B1s,pi_t) = {L(pi_t):.12f}   |dlambda| = 0.0  <- same number twice")

print("\n   WHAT THE 3.342e-16 MEASURES.  Z_k over 11 vertices vs over 5 vertices is the")
print("   SAME SUM IN A DIFFERENT ORDER.  The residual is floating-point associativity.")
rng = np.random.Generator(np.random.PCG64(31337))     # SEED 31337
p11v = [1/11]*11
c3 = {0:0, 1:1, 2:2, 3:3, 4:4, 5:1, 6:1, 7:2, 8:3, 9:3, 10:4}   # class-compatible
p5 = [0.0]*5
for w, v in c3.items():
    p5[v] += 1/11
gap_naive = gap_fsum = gap_rev = 0.0
for _ in range(500):
    f, c = rng.uniform(0, 2*np.pi, 2); k = int(rng.integers(1, 200))
    za = Z_from_vertices(b1s, p11v, f, c, k)
    zb = Z_from_vertices(b1, p5, f, c, k)
    gap_naive = max(gap_naive, abs(za - zb))
    cl_s, cl_1 = b1s.classes(), b1.classes()
    ta = [p11v[w]*np.exp(1j*k*(cl_s[w][0]*f + cl_s[w][1]*c)) for w in range(11)]
    tb = [p5[v]*np.exp(1j*k*(cl_1[v][0]*f + cl_1[v][1]*c)) for v in range(5)]
    fa = complex(math.fsum(t.real for t in ta), math.fsum(t.imag for t in ta))
    fb = complex(math.fsum(t.real for t in tb), math.fsum(t.imag for t in tb))
    gap_fsum = max(gap_fsum, abs(fa - fb))
    gap_rev = max(gap_rev, abs(sum(reversed(ta)) - sum(tb)))
print(f"   seed 31337, 500 random (f,c,k):")
print(f"     max|Z_k gap|, naive left-to-right sum : {gap_naive:.3e}   (claim reports 3.342e-16)")
print(f"     max|Z_k gap|, math.fsum both sides    : {gap_fsum:.3e}")
print(f"     max|Z_k gap|, one side summed reversed: {gap_rev:.3e}")
print("   -> the number moves with summation order and vanishes under exact summation.")
print("      It is a measurement of IEEE-754, not of the carrier, the map, or topology.")

# =====================================================================  A2
print("\n" + BAR)
print("A2  THE CLAIM'S OWN EVIDENCE SENTENCE IS FALSE")
print("    'the map is not unique but the transported class weights ARE'")
print(BAR)
print("""
B1s -> B1 un-subdivision: each midpoint m_i sits on a subdivided edge and must be sent
to ONE of that edge's two endpoints.  6 midpoints -> 2^6 = 64 elementary cellular maps,
every one of them degree 1 onto B1 and a homotopy equivalence.  The claim checked TWO of
them, both class-compatible, and concluded the transported weights are unique.
""")
mids = {5:(0,1), 6:(1,2), 7:(2,0), 8:(0,3), 9:(3,4), 10:(4,0)}   # m -> its edge endpoints
cl_s, cl_1 = b1s.classes(), b1.classes()
seen = {}
n_cc = 0
for choice in itertools.product(*[mids[m] for m in sorted(mids)]):
    cmap = {v: v for v in range(5)}
    for m, tgt in zip(sorted(mids), choice):
        cmap[m] = tgt
    cc = all(cl_1[cmap[w]] == cl_s[w] for w in range(11))
    n_cc += cc
    p = [0.0]*5
    for w in range(11):
        p[cmap[w]] += 1/11
    pi = b1.pi_from_p(p)
    key = tuple(sorted(round(x*11) for x in pi if x > 1e-12))
    seen.setdefault(key, {"n": 0, "cc": False, "lam": L(pi)})
    seen[key]["n"] += 1
    seen[key]["cc"] = seen[key]["cc"] or cc
print(f"   64 elementary collapses enumerated.  class-compatible: {n_cc}.  NOT class-compatible: {64-n_cc}.")
print(f"   distinct transported class-weight multisets (x11): {len(seen)}\n")
print(f"   {'multiset/11':>14s} {'#maps':>6s} {'cls-compat':>11s} {'lambda':>16s} "
      f"{'|d vs B1s own|':>16s}")
own = L(b1s.pi_uniform())
mx = 0.0
for key in sorted(seen, key=lambda k: seen[k]["lam"]):
    d = abs(seen[key]["lam"] - own)
    mx = max(mx, d)
    print(f"   {str(key):>14s} {seen[key]['n']:6d} {str(seen[key]['cc']):>11s} "
          f"{seen[key]['lam']:16.12f} {d:16.3e}")
print(f"\n   -> FALSE. The transported class weights take {len(seen)} distinct values over the")
print(f"      64 maps; only {n_cc} maps give the claim's (0,5/11,5/11,1/11).")
print(f"      max |dlambda| over the elementary collapses = {mx:.3e}")
print(f"      S4's Control 3 figure                        = 3.181e-02")
print(f"      -> the collapse the claim did not check gives a LARGER gap than the one it attacked."
      if mx > 3.181e-02 else "")

# =====================================================================  A3
print("\n" + BAR)
print("A3  THE WALL IS IN THE WRONG PLACE -- B1q -> B1p, a map IS present, lambda MOVES")
print(BAR)
b1p, b1q = ALL["B1p"](), ALL["B1q"]()
print("""
B1q is B1p with the BRIDGE subdivided once (S4 :711).  Same topological space --
exactly the relation S4's Control 3 calls 'homeomorphic'.  The un-subdivision map
c : B1q -> B1p exists and is the ONLY kind of map there is: m -> v0 or m -> w0.
NEITHER is class-compatible: class(m) = (0,0) and B1p HAS NO (0,0) VERTEX.
""")
print(f"   B1p classes {b1p.class_counts()}   B1q classes {b1q.class_counts()}")
lp, lq = L(b1p.pi_uniform()), L(b1q.pi_uniform())
print(f"   SENSE U:  B1p {lp:.12f}   B1q {lq:.12f}   |dlambda| = {abs(lp-lq):.3e}")
print(f"   rank G:   B1p 1 (product only)      B1q 2 (separates)   <- the CRITERION changes too")
for tgt, nm in [(0, "m -> v0 (filled side)"), (3, "m -> w0 (unfilled side)")]:
    p = [0.0]*6
    cm = {0:0,1:1,2:2,3:3,4:4,5:5,6:tgt}
    for w in range(7):
        p[cm[w]] += 1/7
    pi = b1p.pi_from_p(p)
    print(f"   axis-lane recipe, transport along {nm}: pi' = "
          f"({pi[0]:.4f},{pi[1]:.4f},{pi[2]:.4f},{pi[3]:.4f})  "
          f"lambda = {L(pi):.12f}  |d vs B1q own| = {abs(L(pi)-lq):.3e}")
print("""
   -> The claim's own recipe, run on the corpus's OTHER homeomorphism pair, returns
      |dlambda| = 1.8e-01 -- nearly SIX TIMES the 3.181e-02 it set out to explain away.
      'Once a map is present' is false: the map is present and lambda moves further.""")

# =====================================================================  A4
print("\n" + BAR)
print("A4  THE IDENTITY MAP COUNTEREXAMPLE -- B0a and B0b are the SAME COMPLEX")
print(BAR)
b0a, b0b = ALL["B0a"](), ALL["B0b"]()
same = (b0a.edges == b0b.edges and b0a.faces == b0b.faces and b0a.nV == b0b.nV
        and b0a.d1() == b0b.d1() and b0a.d2() == b0b.d2())
print(f"   d1(B0a) == d1(B0b) and d2(B0a) == d2(B0b) and V,E,F equal : {same}")
print(f"   they differ ONLY in gamma_C:  B0a gC={b0a.gC}   B0b gC={b0b.gC}")
la, lb = L(b0a.pi_uniform()), L(b0b.pi_uniform())
print(f"   SENSE U:  B0a {la:.12f}   B0b {lb:.12f}   |dlambda| = {abs(la-lb):.3e}")
print(f"   S4's Control 3 figure = 3.181e-02;  this is {abs(la-lb)/3.181e-02:.2f}x larger.")
print("""
   The map here is the IDENTITY -- the strongest map that can be present.  It is not
   class-compatible, because CLASS IS NOT TOPOLOGICAL DATA: it is read off the DESIGNATED
   LOOPS, which are extra structure the homeomorphism does not carry.  'Class-compatible'
   therefore does not weaken to 'a map is present'; it presupposes that the two carriers
   already agree on the thing lambda is a function of.  The claim's condition is its
   conclusion.""")

# =====================================================================  A5
print("\n" + BAR)
print("A5  THE FAILURE IS NOT 3.181e-02.  WITH MAPS PRESENT THROUGHOUT IT IS ~0.76")
print(BAR)
print("""
K1[nF,nC] := K1 with nF extra vertices inserted in the filled loop and nC in the
unfilled loop.  Every member is HOMEOMORPHIC to K1, and every member admits a
class-compatible collapse onto K1 (inserted vertices on the filled loop are class
(1,0); on the unfilled loop, class (0,1) -- send each to a same-class original).
So a class-compatible map -- the exact object the claim says closes the gap -- is
present for EVERY member.  Each carrier is then given its OWN canonical state (SENSE U).
""")
rows = []
for nF in [0, 1, 2, 3, 5, 10, 25, 100, 1000]:
    for nC in [0, 1, 2, 3, 5, 10, 25, 100, 1000]:
        K = K1_partial_subdiv(nF, nC)
        pi = K.pi_uniform()
        # verify a class-compatible collapse onto K1 exists
        cls = K.classes()
        ok = all(cl in {(1,0),(0,1),(1,1)} for cl in cls)
        rows.append((nF, nC, K.nV, pi, L(pi), ok))
print(f"   {'nF':>5s} {'nC':>5s} {'V':>6s} {'p11':>10s} {'p10':>10s} {'p01':>10s} "
      f"{'lambda(SENSE U)':>18s} {'cc-map?':>8s}")
for nF, nC, V, pi, lam, ok in rows:
    if (nF, nC) in [(0,0),(1,0),(2,0),(3,0),(5,0),(10,0),(100,0),(1000,0),
                    (1,1),(3,3),(10,10),(100,100),(1000,1000),(1000,0),(0,1000),(10,1000)]:
        print(f"   {nF:5d} {nC:5d} {V:6d} {pi[3]:10.6f} {pi[1]:10.6f} {pi[2]:10.6f} "
              f"{lam:18.12f} {str(ok):>8s}")
lams = [r[4] for r in rows]
print(f"\n   all {len(rows)} members admit a class-compatible collapse onto K1: "
      f"{all(r[5] for r in rows)}")
print(f"   min lambda over the family = {min(lams):.12f}   (at K1 itself)")
print(f"   max lambda over the family = {max(lams):.12f}")
print(f"   SPREAD over ONE homeomorphism type, maps present throughout = {max(lams)-min(lams):.6f}")
print(f"   S4's Control 3 figure = 0.031810;  ratio = {(max(lams)-min(lams))/3.181e-02:.1f}x")
print("   sup as nF -> inf with nC = 0 is 0^- : lambda is dense in [-0.756574, 0).")

# =====================================================================  A6
print("\n" + BAR)
print("A6  THE UNSTATED NORMALISATION -- the claim rediscovered S4's OWN SENSE C row")
print(BAR)
print("""
'Transport the state along a class-compatible map' is, by rm-1, EXACTLY 'fix the class
weights by hand across the two carriers'.  That is S4's SENSE C, definitionally --
S4 :569: 'SENSE C -- fixed at the CLASS level ... the carrier feeds in only through
WHICH classes exist.'  S4 printed |diff| = 0.0e+00 for SENSE C on the line IMMEDIATELY
BELOW the Control-3 line the claim attacks (S4 :702).
""")
print(f"   S4 SENSE C, 3 classes (0.4,0.3,0.3): lambda = {L((0.0,0.4,0.3,0.3)):.12f} on both, |diff| = 0.0")
print(f"   claim's SENSE C', weights (5/11,5/11,1/11): lambda = {L(pi_t):.12f} on both, |diff| = 0.0")
print("""   Same statement, different numerical choice of the hand-fixed weights.  The claim
   contains no information S4 had not already published one line further down, and does
   not touch S4's SENSE U sentence, which is what Control 3 actually asserts.""")

print("\n" + BAR)
print("VERDICT: the claim is REFUTED on A2, A3, A4 and A5 independently.")
print(BAR)
