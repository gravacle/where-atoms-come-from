"""
rm_4_homeo.py -- LANE R (MAPS REFUTER).
Verify AT THE INCIDENCE that each pair this refutation calls 'homeomorphic' really is
an elementary-subdivision pair, so that the counterexamples are counterexamples.
An elementary subdivision of edge (a,b) inserts a vertex m of degree exactly 2 whose
two neighbours are a and b, and deletes (a,b).  |K| and |K'| are then homeomorphic.
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from rm_lib import *

def deg(K, v):
    return sum((t == v) + (h == v) for (t, h) in K.edges)

def check_subdivision(Kbig, Ksmall, new_vertices, label):
    print(f"\n   {label}")
    print(f"     V {Ksmall.nV} -> {Kbig.nV}   E {Ksmall.nE} -> {Kbig.nE}   "
          f"F {Ksmall.nF} -> {Kbig.nF}")
    ok = (Kbig.nV - Ksmall.nV == len(new_vertices)
          and Kbig.nE - Ksmall.nE == len(new_vertices)
          and Kbig.nF == Ksmall.nF)
    print(f"     +1 vertex and +1 edge per inserted vertex, faces unchanged: {ok}")
    allok = ok
    small_edges = set(frozenset(e) for e in Ksmall.edges)
    for m in new_vertices:
        d = deg(Kbig, m)
        nb = sorted({x for (t, h) in Kbig.edges if m in (t, h) for x in (t, h)} - {m})
        parent_present = frozenset(nb) in small_edges if len(nb) == 2 else False
        print(f"     inserted vertex {m:2d}: degree {d}, neighbours {nb}, "
              f"parent edge {tuple(nb) if len(nb)==2 else '?'} present in the small "
              f"complex: {parent_present}")
        allok = allok and d == 2 and len(nb) == 2 and parent_present
    b_big, b_small = Kbig.betti()[:3], Ksmall.betti()[:3]
    print(f"     Betti (b0,b1,b2): {b_small} -> {b_big}   equal: {b_big == b_small}")
    print(f"     chi: {Ksmall.chi()} -> {Kbig.chi()}   equal: {Ksmall.chi()==Kbig.chi()}")
    print(f"     VERDICT: elementary-subdivision pair, |K| homeomorphic to |K'|: "
          f"{allok and b_big == b_small}")
    return allok and b_big == b_small

print("=" * 78)
print("RM-4  THE HOMEOMORPHISM CHECKS")
print("=" * 78)
r = []
r.append(check_subdivision(ALL["B1s"](), ALL["B1"](), [5, 6, 7, 8, 9, 10],
                           "B1s vs B1  -- S4's Control 3 pair (every edge subdivided)"))
r.append(check_subdivision(ALL["B1q"](), ALL["B1p"](), [6],
                           "B1q vs B1p -- the pair this lane uses as counterexample (A3)"))
r.append(check_subdivision(K1_partial_subdiv(1, 0), ALL["B1"](), None or [3],
                           "K1[1,0] vs B1 -- one insertion in the filled loop (A5)")
         if False else None)
# K1_partial_subdiv relabels vertices, so verify that family by invariants instead.
print("\n   K1[nF,nC] family (A5): verified by construction and by invariants --")
base = ALL["B1"]()
print(f"     {'nF':>4s} {'nC':>4s} {'V':>5s} {'E':>5s} {'F':>3s} {'chi':>4s} "
      f"{'(b0,b1,b2)':>12s} {'all deg-2 inserts':>18s}")
famok = True
for nF, nC in [(0,0),(1,0),(0,1),(2,0),(3,3),(5,2),(10,10)]:
    K = K1_partial_subdiv(nF, nC)
    b = K.betti()[:3]
    inserted = [v for v in range(1, K.nV) if deg(K, v) == 2]
    good = (K.chi() == base.chi() and b == base.betti()[:3]
            and K.nV == 5 + nF + nC and K.nE == 6 + nF + nC and K.nF == 1
            and len(inserted) == K.nV - 1)
    famok = famok and good
    print(f"     {nF:4d} {nC:4d} {K.nV:5d} {K.nE:5d} {K.nF:3d} {K.chi():4d} "
          f"{str(b):>12s} {str(len(inserted)==K.nV-1):>18s}")
print(f"     every member has chi=0, (b0,b1,b2)=(1,1,0), V-1 vertices of degree 2 and one")
print(f"     vertex of degree 4 (the pinch v0): the homeomorphism type of K1.  ok = {famok}")

print("\n" + "=" * 78)
print(f"ALL HOMEOMORPHISM CHECKS PASS: {all(x for x in r if x is not None) and famok}")
print("=" * 78)
print("""
NOTE ON WHAT THE HOMEOMORPHISM DOES AND DOES NOT CARRY.
The homeomorphism |B1s| -> |B1| is the IDENTITY of the underlying space.  It does NOT
send vertices to vertices: the six midpoints go to interior points of B1's edges, which
carry no fibre.  So there is NO transport of a ready state along the homeomorphism at
all.  The claim's c3 is therefore not the homeomorphism; it is a degree-1 cellular
homotopy equivalence chosen from 64, and then chosen class-compatible from 4.  Two
substitutions, neither ledgered.""")
