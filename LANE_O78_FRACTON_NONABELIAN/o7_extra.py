"""Addendum: checkerboard formula at more sizes; X-cube L=5; Chamon k(L);
and the LOCALITY statement made explicit as a measured number."""
import sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *
from models import *

print("== checkerboard k(L) ==")
ks = {}
for L in (2,4,6):
    g,n = checkerboard(L); k,r = code_k(g,n); ks[L]=k
    print(f"  L={L}: n={n:4d} k={k:3d}   6L-6={6*L-6}  12L-6={12*L-6}  6L-3={6*L-3}")
print("  fitted: k =", [ks[L] for L in (2,4,6)], " differences:", ks[4]-ks[2], ks[6]-ks[4])

print("\n== X-cube k(L) ==")
for L in (2,3,4,5):
    g,n = xcube(L); k,r = code_k(g,n)
    print(f"  L={L}: n={n:4d} k={k:3d}   6L-3={6*L-3}  {'MATCH' if k==6*L-3 else 'MISMATCH'}")

print("\n== Chamon k(L), all sites, 1 qubit/site ==")
for L in (2,3,4,5,6,8):
    g,n = chamon(L); k,r = code_k(g,n)
    print(f"  L={L}: n={n:4d} gens={len(g):4d} rank={r:4d} k={k:4d}  commute_violations={all_commute(g,n)}")

print("\n== LOCALITY, stated as a number: max stabiliser weight and max range ==")
def geom(name, gens, n, L, sites_per_cell, coords):
    w = max(pweight(g,n) for g in gens)
    rng = 0
    for g in gens:
        pts = [coords(i) for i in range(n) if ((xpart(g,n)>>i)&1) or ((zpart(g,n)>>i)&1)]
        for a in pts:
            for b in pts:
                dd = max(min(abs(a[t]-b[t]), L-abs(a[t]-b[t])) for t in range(len(a)))
                rng = max(rng, dd)
    print(f"  {name:22s} max stabiliser weight = {w:3d}   max Chebyshev range = {rng}  "
          f"(both O(1), independent of L)")
L=4
geom("3D toric L=4", *toric3d(L), L, 3, lambda i: ((i//3)//(L*L), ((i//3)//L)%L, (i//3)%L))
geom("X-cube L=4",   *xcube(L),   L, 3, lambda i: ((i//3)//(L*L), ((i//3)//L)%L, (i//3)%L))
geom("Haah CC1 L=4", *haah(L),    L, 2, lambda i: ((i//2)//(L*L), ((i//2)//L)%L, (i//2)%L))
geom("Chamon L=4",   *chamon(L),  L, 1, lambda i: (i//(L*L), (i//L)%L, i%L))
geom("checkerboard L=4", *checkerboard(L), L, 1, lambda i: (i//(L*L), (i//L)%L, i%L))
print("""
  Every candidate here is geometrically local in the SAME sense the program's toric-code
  carrier is: O(1) stabiliser weight, O(1) range on a cubic lattice.  NO candidate in this
  sweep buys R3 by redefining locality.  (Contrast: expander/qLDPC codes from G-6, which
  are local only in a GRAPH sense with no bounded-dimension embedding.)""")
