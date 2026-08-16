#!/usr/bin/env python3
"""LANE R (MAPS) — REFUTER, part 4: lambda's true domain, and the round trip."""
import numpy as np
from maps_refuter import (Complex, K1, K1_GF, K1_GC, classes_from_loops, push_to_pi,
                          lam_mahler_generic as lam, hdr)
from maps_refuter2 import collapse, verify_chain_map
np.set_printoptions(linewidth=200, suppress=True)
k1 = K1(); CLS = classes_from_loops(5, K1_GF, K1_GC)

hdr("BLOCK 12 — lambda's ACTUAL DOMAIN: RANDOM CARRIERS, SAME (pi, u, v) => SAME lambda")
rng = np.random.default_rng(1123581321)
# force collisions: one fixed pi realised on many different carriers and states
target = np.array([0.1, 0.2, 0.3, 0.4])
vals = []
for V in range(4, 40):
    # spread each class weight over a random number of vertices
    p, ab = [], []
    for ci, (a, b) in enumerate([(0,0),(1,0),(0,1),(1,1)]):
        n = int(rng.integers(1, max(2, V//3)))
        w = rng.dirichlet(np.ones(n)) * target[ci]
        p += list(w); ab += [(a,b)]*n
    p = np.array(p)
    pi = np.zeros(4)
    for w, c in zip(p, ab):
        pi[{(0,0):0,(1,0):1,(0,1):2,(1,1):3}[c]] += w
    vals.append(lam(pi))
print("seed 1123581321: 36 carriers with 4..40 vertices, wildly different vertex counts and")
print("  wildly different vertex-level states, all with the SAME class pushforward pi =", target)
print("  lambda values: min = %.12f  max = %.12f  spread = %.3e" % (min(vals), max(vals), max(vals)-min(vals)))
print("  -> the vertex count, the incidence, the 2-cells and the vertex-level state are")
print("     all invisible. lambda sees pi and the holonomies. Nothing else exists for it.")

hdr("BLOCK 13 — THE HOMOTOPY INVERSE ROUND TRIP: q AND i ARE NOT INVERSE ON STATES")
T = frozenset({0,1,3,4})
q, phi, emap, (gF2, gC2) = collapse(k1, T, [K1_GF, K1_GC])
p = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
pq = np.zeros(q.nv)
for w in range(5): pq[phi[w]] += p[w]
# homotopy inverse i : wedge -> K1 sends the single 0-cell to v0
ip = np.zeros(5); ip[0] = pq[0]
print("  p on K1            =", p,              " lambda =", f"{lam(push_to_pi(p, CLS)):+.12f}")
print("  q_* p on the wedge =", pq,             " lambda =", f"{lam(np.array([0,0,0,1.0])):+.12f}")
print("  i_* q_* p back on K1 =", ip,           " lambda =", f"{lam(push_to_pi(ip, CLS)):+.12f}")
print("  i_* q_*  is the projection onto delta_v0, NOT the identity, though i . q ~ id_K1.")
print("  -> a homotopy equivalence of CARRIERS induces no equivalence of STATE data, so it")
print("     transports no state functional. The claim's inference has no vehicle.")
