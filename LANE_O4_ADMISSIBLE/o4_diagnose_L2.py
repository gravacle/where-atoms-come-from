#!/usr/bin/env python3
"""
DIAGNOSTIC.  o4_toric.py's Section 2 reported 4 admissible R-flippers supported inside
a single PLAQUETTE at L = 2.  If that survives, clause (v) fails even under DEF-A.
Print them and decide whether it is a result or a small-lattice artefact.
"""
import numpy as np

L = 2
N = 2 * L * L


def eidx(x, y, d):
    return 2 * (L * (y % L) + (x % L)) + d


NAME = {eidx(x, y, d): f"({x},{y},{'H' if d == 0 else 'V'})"
        for y in range(L) for x in range(L) for d in range(2)}

STARS = [sorted({eidx(x, y, 0), eidx(x - 1, y, 0), eidx(x, y, 1), eidx(x, y - 1, 1)})
         for y in range(L) for x in range(L)]
PLAQS = [sorted({eidx(x, y, 0), eidx(x, y + 1, 0), eidx(x, y, 1), eidx(x + 1, y, 1)})
         for y in range(L) for x in range(L)]


def vec(edges):
    a = np.zeros(N, dtype=np.int8)
    for e in edges:
        a[e] = 1
    return a


gens = [(vec(s), np.zeros(N, np.int8)) for s in STARS] + \
       [(np.zeros(N, np.int8), vec(p)) for p in PLAQS]
R_sym = (np.zeros(N, np.int8), vec({eidx(0, 0, 0), eidx(1, 0, 0)}))


def sympl(p, q):
    return int((p[0] @ q[1] + p[1] @ q[0]) % 2)


print("R = Z on", [NAME[e] for e in sorted({eidx(0, 0, 0), eidx(1, 0, 0)})])
print()
for pi, p in enumerate(PLAQS):
    print(f"PLAQUETTE {pi}: edges {[NAME[e] for e in p]}")
    for xb in range(16):
        for zb in range(16):
            x = np.zeros(N, np.int8)
            z = np.zeros(N, np.int8)
            for k, e in enumerate(p):
                if xb >> k & 1:
                    x[e] = 1
                if zb >> k & 1:
                    z[e] = 1
            if not x.any() and not z.any():
                continue
            if all(sympl((x, z), g) == 0 for g in gens) and sympl((x, z), R_sym) == 1:
                lab = "".join(f"{'Y' if x[e] and z[e] else 'X' if x[e] else 'Z'}_{NAME[e]} "
                              for e in range(N) if x[e] or z[e])
                print(f"    ADMISSIBLE FLIPPER: {lab}")
    print()

print("""
IS THE OFFENDING SUPPORT A CONTRACTIBLE SET?  On a 2x2 torus a plaquette's four edges
are  bottom (x,y,H), top (x,y+1,H), left (x,y,V), right (x+1,y,V).  With L = 2 the
left and right edges of ANY plaquette are the two DISTINCT vertical edges of that row,
and X on both of them is a closed dual loop that WINDS AROUND THE TORUS.  So a single
plaquette's edge set at L = 2 already carries a non-contractible cycle:
""")
for pi, p in enumerate(PLAQS[:1]):
    print(f"   plaquette {pi} edges: {[NAME[e] for e in p]}")
    print(f"   its two V edges     : {[NAME[e] for e in p if e % 2 == 1]}"
          f"   <- an X-string on these two is a horizontal dual loop, weight 2 = d")
print("""
CONCLUSION: at L = 2 the code distance is d = 2 and a plaquette holds 4 edges, so a
"contractible region" as an edge set is NOT contractible -- it contains a winding
cycle.  This is a SMALL-LATTICE ARTEFACT, not a failure of clause (v).  The clause-(v)
test must be run at L >= 3, where d = 3 > the 4-edge plaquette cannot hold a logical.
o4_L3.py does that with exact F_2 linear algebra rather than enumeration.
""")
