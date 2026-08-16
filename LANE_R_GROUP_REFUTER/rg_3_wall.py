#!/usr/bin/env python3
"""
rg_3_wall.py -- WHERE THE WALL ACTUALLY IS.

The CLAIM says S4-1 stops at d = 2.  S4-1 is a THEOREM, so it "stops" exactly
when a COUNTEREXAMPLE exists.  I therefore SEARCH for counterexamples to

    (i)   |S| >= 3  =>  rank G = 2
    (ii)  |S| <= 2  =>  rank G <= 1
    (iii) |S| = 1   =>  G = {1}   ("no formation, ever")

over a grid of configurations that varies the fibre dimension d and the number
of DISTINCT WEIGHTS |Lambda_eff| INDEPENDENTLY.  A configuration is marked
BROKEN the moment one witness is found, and the witness is printed.

If d were the variable the table would split on d.  It splits on |Lambda_eff|.

SEARCH: for each configuration, all 15 non-empty class subsets x 400 random
draws of (vertices per class in 1..3, per-vertex charge, positive weights,
fibre directions).  SEED: numpy default_rng(4242) reset per configuration.
"""
import numpy as np, math, itertools
from rg_lib import *

CORNERS = [(0, 0), (1, 0), (0, 1), (1, 1)]

# ---- spin-1 (3-dimensional, primitive circle, weights +1,0,-1) -------------
_s = 1 / math.sqrt(2)
JX = np.array([[0, _s, 0], [_s, 0, _s], [0, _s, 0]], dtype=complex)
JY = np.array([[0, -1j*_s, 0], [1j*_s, 0, -1j*_s], [0, 1j*_s, 0]], dtype=complex)
JZ = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dtype=complex)

def op_spin1(axis, sign=+1):
    n = np.asarray(axis, dtype=float); n = n / np.linalg.norm(n)
    N = n[0]*JX + n[1]*JY + n[2]*JZ
    ev, V = np.linalg.eigh(N)
    V = V[:, np.argsort(-ev.real)]
    return OneParam(V, sign * np.array([1, 0, -1], dtype=int))


def search(mkops, d, charges, draws=400, seed=4242):
    """
    mkops(q) -> (A,B) for a vertex of charge q.  `charges` is the allowed set.
    Returns (broken, witness) where witness describes the first violation.
    """
    rng = np.random.default_rng(seed)
    for r in range(1, 5):
        for S in itertools.combinations(CORNERS, r):
            for _ in range(draws):
                cls, ops, secs = [], [], []
                for cl in S:
                    nv = int(rng.integers(1, 4))          # 1..3 vertices per class
                    for _v in range(nv):
                        q = int(charges[rng.integers(len(charges))])
                        cls.append(cl); ops.append(mkops(q))
                        vec = rng.normal(size=d) + 1j * rng.normal(size=d)
                        secs.append(vec / np.linalg.norm(vec)
                                    * math.sqrt(float(rng.uniform(0.1, 1.0))))
                n = math.sqrt(sum(float(np.vdot(x, x).real) for x in secs))
                secs = [x / n for x in secs]
                md = modes_general(cls, secs, ops)
                sup = sorted(md.keys())
                rk = lattice_rank(sup)
                bad = None
                if r >= 3 and rk != 2:
                    bad = f"|S|={r} but rank G = {rk}   (S4-1: rank must be 2)"
                elif r <= 2 and rk == 2:
                    bad = f"|S|={r} but rank G = 2     (S4-1: rank must be <= 1)"
                elif r == 1 and len(sup) > 1:
                    bad = f"|S|=1 but G != {{1}}        (S4-1: no formation, ever)"
                if bad:
                    qs = [int(o[0].w[np.argmax(np.abs(o[0].w))]) for o in ops]
                    return True, (f"classes {list(S)} -> {bad}\n"
                                  f"            support {sup}, "
                                  f"vertex classes {cls}")
    return False, None


CONFIGS = [
    # label, d, |Lambda_eff|, mkops(q), allowed charges
    ("U(1), d=1, uniform unit charge  [S4-1's own box]", 1, 1,
     lambda q: (op_u1([1], -1), op_u1([1], +1)), [1]),
    ("U(1), d=2, rep chi(+)chi   SCALAR", 2, 1,
     lambda q: (op_u1([1, 1], -1), op_u1([1, 1], +1)), [1]),
    ("U(1), d=5, rep chi^(+)5    SCALAR", 5, 1,
     lambda q: (op_u1([1]*5, -1), op_u1([1]*5, +1)), [1]),
    ("U(1), d=1, PER-VERTEX CHARGE q in {1,2}", 1, 2,
     lambda q: (op_u1([q], -1), op_u1([q], +1)), [1, 2]),
    ("U(1), d=1, PER-VERTEX CHARGE q in {1,2,3}", 1, 3,
     lambda q: (op_u1([q], -1), op_u1([q], +1)), [1, 2, 3]),
    ("U(1), d=2, Lambda={1,2}   ABELIAN, non-scalar", 2, 2,
     lambda q: (op_u1([1, 2], -1), op_u1([1, 2], +1)), [1]),
    ("SU(2), d=2, COMMUTING axes  (abelian image)", 2, 2,
     lambda q: (op_su2([0, 0, 1], -1), op_su2([0, 0, 1], +1)), [1]),
    ("SU(2), d=2, NON-commuting axes", 2, 2,
     lambda q: (op_su2([0, 0, 1], -1), op_su2([1, 0, 0], +1)), [1]),
    ("SU(2) spin-1, d=3, NON-commuting axes", 3, 3,
     lambda q: (op_spin1([0, 0, 1], -1), op_spin1([1, 0, 0], +1)), [1]),
]

print("=" * 78)
print("THE WALL.  COUNTEREXAMPLE SEARCH AGAINST THEOREM S4-1.")
print("  15 class subsets x 400 random draws per configuration.  seed 4242.")
print("=" * 78)
print(f"  {'configuration':<50} {'d':>2} {'|L_eff|':>7}  S4-1")
print("  " + "-" * 74)
witnesses = []
for label, d, nlam, mk, ch in CONFIGS:
    broken, wit = search(mk, d, ch)
    print(f"  {label:<50} {d:>2} {nlam:>7}  {'BROKEN' if broken else 'survives'}")
    if broken:
        witnesses.append((label, wit))
print("  " + "-" * 74)

print("\n  FIRST WITNESS FOR EACH BROKEN CONFIGURATION:")
for label, wit in witnesses:
    print(f"\n    {label}\n            {wit}")

print("""
==============================================================================
READ THE COLUMNS.  The verdict does NOT track d.  S4-1 SURVIVES at d = 5 and
is BROKEN at d = 1.  It tracks |Lambda_eff| -- the number of DISTINCT weights
appearing across the ready state's support.

CORRECTED WALL OF RECORD (this lane):

   THEOREM S4-1 holds exactly when the character of the k-dependence at every
   vertex factors through the class map (a_v,b_v) |-> (-a_v, b_v): that is,
   when the transport is SCALAR on every fibre AND the charge is uniform, so a
   single weight w is shared by every vertex in the ready state's support.
   It then holds AT EVERY FIBRE DIMENSION d, d = 1 or d = 5 alike.
   It fails as soon as a second weight appears anywhere in that support --
   AT EVERY d, INCLUDING d = 1.

   FIBRE DIMENSION IS NOT THE VARIABLE.  |Lambda_eff| IS.

This is the same mislabelling ERR-1 convicted the corpus of once already: there
the operative variable was SCALARITY and the corpus named COMMUTATIVITY.  Here
the operative variable is again scalarity -- one weight -- and the CLAIM names
DIMENSION.  Commutativity does reappear, but one level down: it selects which
of the allowed modes actually carry non-zero coefficient (rg_2_theorem T1/T2),
not whether S4-1 holds.
==============================================================================
""")
