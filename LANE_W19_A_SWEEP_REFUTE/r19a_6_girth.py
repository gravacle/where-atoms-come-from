# W19-A REFUTER / step 6.  F8 AND THE CLAIM "ONLY GIRTH LENGTHENS THE PLATEAU".
#
# Lane A, F8: "Square-plaquette lattices cannot reach the plateau bar at any volume.  Girth 4 fixes
#              d = 3 and therefore exactly 2 plateau points, independent of lattice size ...
#              Volume does not lengthen the plateau; only girth does."   status: PROVED
#
# Girth 4 does NOT fix d = 3.  d is set by girth_through(l), the shortest cycle THROUGH the system
# link.  A graph can have girth 4 while the system link lies on no short cycle at all.  sq_chain15
# below is exactly that: min degree 3, girth 4, d = 5 -> 4 plateau points at L = 15.  (Its cycle
# basis is 4 squares + 2 hexagons, not all squares; that limit is printed in the output, not hidden.)
#
# The same table shows girth is neither necessary nor sufficient: tri_chain12 has girth 3 -- WORSE
# than every square lattice -- and d = 5, while petersen has girth 5 and only d = 4.
import numpy as np, sys, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP_REFUTE")
from zn_gauge import ZNGauge, S_of, mutual_information, level_cuts, nested_fragments, girth_through
from carriers import heawood, petersen, cube, grid, ladder, LADDER
from refute_carriers import tri_chain12, dbl_chain9, tri_chain11, mg_chain, min_degree

DELTA = 0.10
print("=" * 122)
print("W19-A REFUTER / 6  GIRTH DOES NOT SET THE PLATEAU LENGTH; GIRTH-THROUGH-l DOES")
print("=" * 122)

def sq_chain15():
    """SIMPLE, min degree 3, GIRTH 4, d = 5.  Minimum cycle basis = 4 squares + 2 hexagons.
       Levels from u in G-l: {0} {1,2} {3,4} {5,6} {7,8} {9}=v, joined as K_{2,2} rungs.
       V = 10, L = 15."""
    u, a, b, c, d_, e, f, gg, h, v = range(10)
    E = [(u, v),                                    # l, index 0
         (u, a), (u, b),
         (a, c), (a, d_), (b, c), (b, d_),          # K_{2,2}: makes 4-cycles, girth 4
         (c, e), (d_, f),
         (e, gg), (e, h), (f, gg), (f, h),          # K_{2,2}
         (gg, v), (h, v)]
    return 10, E

def girth(V, E):
    return min(girth_through(V, E, l)[0] for l in range(len(E)))

def run(nm, V, E, gsq=0.5):
    g = ZNGauge(nm, V, E, 2)
    l = max(range(g.L), key=lambda i: girth_through(V, E, i)[1] or -1)
    d = girth_through(V, E, l)[1]
    psi, E0, gap = g.ground(2.0 / gsq, 2.0 * gsq); Psi = g.full_vector(psi)
    HS = S_of(Psi, g.L, 2, [l]); frs, _ = nested_fragments(V, E, l); cuts, _ = level_cuts(V, E, l)
    rr = [mutual_information(Psi, g.L, 2, [l], F) / HS for F in frs]
    pts = sum(abs(x - 1) <= DELTA for x in rr)
    Rd = sum(mutual_information(Psi, g.L, 2, [l], C) / HS >= 1 - DELTA for C in cuts)
    ver = "EXHIBITED" if pts >= 4 else ("MARGINAL" if pts == 3 else "FAIL")
    if HS < 0.10: ver += " (WEIGHTLESS)"
    return dict(nm=nm, V=V, L=g.L, C=g.C, girth=girth(V, E), d=d, HS=HS, pts=int(pts), R=int(Rd),
                ver=ver, plaq=sorted(int(abs(p).sum()) for p in g.plaq),
                sizes=[len(F) for F in frs], ratios=rr)

print("\n[6a] THE SQUARE-PLAQUETTE COUNTEREXAMPLE TO F8.  g^2 = 0.50, Z_2, lane A's criterion.")
print(f"     {'carrier':<14}{'V':>4}{'L':>4}{'C':>4}{'mindeg':>8}{'girth':>7}{'d':>4}{'H(S)':>12}"
      f"{'pts':>5}{'R_d':>5}   plaquette lengths (min cycle basis)")
rows = []
for nm, (V, E) in [("sq_chain15", sq_chain15()), ("grid_3x3_open", grid(3, 3, False)),
                   ("cube_Q3", cube()), ("ladder_5sq", ladder(5)),
                   ("tri_chain12", tri_chain12()), ("petersen", petersen()), ("heawood", heawood())]:
    r = run(nm, V, E); rows.append(r)
    print(f"     {r['nm']:<14}{r['V']:>4}{r['L']:>4}{r['C']:>4}{min_degree(V,E):>8}{r['girth']:>7}{r['d']:>4}"
          f"{r['HS']:>12.9f}{r['pts']:>5}{r['R']:>5}   {r['plaq']}  {r['ver']}")
r = rows[0]
print(f"\n     sq_chain15 detail: |F| = {r['sizes']},  I/H(S) = " + " ".join(f"{x:.9f}" for x in r['ratios']))
nsq = sum(1 for x in r['plaq'] if x == 4)
print(f"     Its minimum cycle basis is {r['plaq']}: {nsq} of {len(r['plaq'])} plaquettes are SQUARES,"
      f" the rest hexagons.")
print("     HONEST LIMIT OF THIS EXHIBIT: it is a girth-4 carrier, NOT a pure square-plaquette one --")
print("     with d = 5 the cycle through l has length 6, so a hexagon is forced into any basis.  I did")
print("     not find, and do not claim, a carrier whose basis is ALL squares with d = 5.")
print(f"     WHAT IS REFUTED: F8's stated mechanism, 'Girth 4 fixes d = 3 and therefore exactly 2")
print(f"     plateau points'.  sq_chain15 has girth {r['girth']} and d = {r['d']}, giving {r['pts']} plateau points at L = 15.")
print("     WHAT SURVIVES: the narrower claim that on a lattice in which EVERY LINK lies on a square")
print("     face, d = 3 for every link.  That is a property of the embedding, not of the girth, and")
print("     F8 states it as the latter.")

print("\n[6b] GIRTH vs GIRTH-THROUGH-l vs d, across every carrier in play.  If girth set the plateau")
print("     length, these two columns would agree.  They do not.")
print(f"     {'carrier':<16}{'L':>4}{'mindeg':>8}{'girth':>7}{'max girth_through(l)':>22}{'d':>4}"
      f"{'pts = d-1':>11}")
ALL = list(LADDER) + [("sq_chain15", sq_chain15()), ("tri_chain12", tri_chain12()),
                      ("tri_chain11", tri_chain11()), ("dbl_chain9", dbl_chain9()),
                      ("mg_chain6", mg_chain(6)), ("mg_chain7", mg_chain(7))]
for nm, (V, E) in ALL:
    gt = max(girth_through(V, E, i)[0] for i in range(len(E)))
    dd = max(girth_through(V, E, i)[1] for i in range(len(E)))
    print(f"     {nm:<16}{len(E):>4}{min_degree(V,E):>8}{girth(V,E):>7}{gt:>22}{dd:>4}{dd-1:>11}")
print("     tri_chain12: girth 3 (worse than any square lattice) yet d = 5 (better than petersen's 4).")
print("     The ordering by girth and the ordering by d are DIFFERENT ORDERINGS.")
json.dump(rows, open("out_6_girth.json", "w"), indent=1, default=float)
print("\nDONE 6.")
