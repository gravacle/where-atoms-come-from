"""STEP 8 -- WHAT THE OPERATIONAL RELATION IS A FUNCTION OF, and the headline table.

CORRECTION TO THE READ OF 7D.  7D regressed C_ij on two SCALARS (the anticommutation bit and
the inner product of the two site profiles) and got R^2 falling from 0.90 to 0.22 as n grows.
That does NOT mean the operational relation carries information beyond the symplectic bit and
the support.  SC-12 already PROVED the stronger statement: the sector decomposition computes
C_ij from the triple (anticommutation bit, profile w_i, profile w_j) and from nothing else, so
C_ij is an exact deterministic function of that triple.  The falling R^2 only says the function
is not linear in those two scalars.  This step counts the degrees of freedom instead: how many
DISTINCT (anti, w_i, w_j) triples the matrix actually contains.  The functional dependence is not
asserted -- it is a consequence of SC-12, which showed the sector construction, whose only inputs
are that triple, reproduces the full-space number to 5e-15.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

NQH = 3
P("=" * 112)
P("LANE_SCALE_B_METRIC  STEP 8 -- DEGREES OF FREEDOM IN THE OPERATIONAL RELATION")
P("=" * 112)
P("")
P("  Ordered pairs in the crowding matrix vs the number of DISTINCT (anti, w_i, w_j) keys the")
P("  exact sector decomposition actually evaluates.  Every pair sharing a key has an identical")
P("  C value, so the key count is the true number of independent entries.")
P("")
P("  %-4s %-6s %-12s %-14s %-14s %-12s" %
  ("n", "2k", "ordered prs", "distinct keys", "distinct profs", "keys/pairs"))
P("  " + "-" * 68)
for n in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    profs = []
    for v in vs:
        w = np.zeros(NQH)
        for q in support(v, n): w[q % NQH] += 1
        profs.append(tuple(np.round(w / w.sum(), 9)))
    keys = set()
    for i in range(m):
        for j in range(m):
            if i == j: continue
            keys.add((bool(sp_form(vs[i], vs[j], n)), profs[i], profs[j]))
    P("  %-4d %-6d %-12d %-14d %-14d %-12.3f" %
      (n, m, m * (m - 1), len(keys), len(set(profs)), len(keys) / (m * (m - 1))))
P("")
P("  READ: the operational relation between two records is fixed by the pair's ANTICOMMUTATION")
P("  BIT and the two records' SUPPORT PROFILES.  Those are exactly the ingredients of relation")
P("  (a) and relation (d).  The operational measurement adds no independent relational content;")
P("  it re-weighs the same two ingredients.")

P("")
P("=" * 112)
P("  HEADLINE TABLE -- intrinsic dimension d90 of the record-relation geometry, by definition,")
P("  with both controls in the same table.  Sources: 02, 03, 06, 07.")
P("")
P("  %-4s %-4s | %-5s %-5s %-6s %-6s %-6s %-6s %-5s | %-6s %-6s %-6s" %
  ("n", "2k", "(a)S", "(c)C", "(d)GS", "(d)RS", "(e)GS", "(e)RS", "(b2)", "FREE-d", "FREE-b2", "RAND-b2"))
P("  " + "-" * 92)
P("  (a) symplectic   (c) code overlap   (d) support Jaccard   (e) Pauli letters   (b2) crowding")
P("  GS = the Gram-Schmidt basis symplectic_logicals returns;  RS = a random symplectic basis of")
P("  records on the SAME carrier (worst of 8 draws, step 3A).  FREE = k unentangled qubits, H = 0.")
P("")
DATA = {
    #     (a) (c) (d)GS (d)RS (e)GS (e)RS (b2) FREEd FREEb2 RANDb2
    4:  (1,  3,  1,  1,  2,  2,  2,  1,  2,  2),
    6:  (3,  7,  3,  3,  4,  4,  4,  3,  5,  4),
    8:  (5, 10,  3,  4,  3,  6,  6,  5,  6,  6),
    10: (7, 14,  4,  5,  3,  8,  7,  7,  7,  8),
    12: (9, 18,  4,  6,  3,  9,  9,  9,  8,  9),
    14: (10, 21, 4,  8,  3, 11, 11, 10,  9, 10),
    16: (12, 25, 4,  8,  2, 12, 13, 12, 10, 11),
    18: (14, 28, 4, 10,  2, 14, 14, 14, 11, 13),
    20: (16, 32, 4, 11,  2, 15, 16, 16, 12, 14),
}
for n, row in DATA.items():
    P("  %-4d %-4d | %-5d %-5d %-6d %-6d %-6d %-6d %-5d | %-6d %-6d %-6d" %
      ((n, 2 * (n - 2)) + row))
xs = np.array([2 * (n - 2) for n in DATA], float)
P("")
names = ["(a) symplectic", "(c) code overlap", "(d) support  GS basis", "(d) support  RS basis",
         "(e) letters  GS basis", "(e) letters  RS basis", "(b2) crowding",
         "FREE control under (d)", "FREE control under (b2)", "RAND control under (b2)"]
for idx, name in enumerate(names):
    ys = np.array([DATA[n][idx] for n in DATA], float)
    sl = np.polyfit(xs, ys, 1)[0]
    P("  %-26s d90 vs 2k slope = %+.4f   %s" %
      (name, sl, "SATURATES" if abs(sl) < 0.05 else "GROWS WITH n"))

P("")
P("  RELATION (b1), THE DRIVE-READ MATRIX (step 5): EXACTLY DIAGONAL at every n from 4 to 20 and")
P("  at lam = 0.4 / 0.8 / 1.2 -- max |off-diagonal| between 3e-16 and 7e-15, i.e. machine zero.")
P("  Its distance transform d = sqrt(Mmax - M_sym) is therefore a regular simplex on 2k points,")
P("  d90 = 2k - 1, the largest dimension the point count allows.  The D-15 positive control in")
P("  5C shows the SAME instrument returning 0.026 / 0.108 / 0.257 as the record overlap is")
P("  dialled to 0.25 / 0.50 / 0.75, so the zero is a real zero and not a dead instrument.")
P("")
P("  BASIS-FREE CHECK (step 6): every logical class at minimum weight, support relation.")
P("  %-6s %-8s %-10s %-10s" % ("n", "#classes", "CODE d90", "FREE d90"))
CLS = {4: (15, 2, 1), 6: (255, 5, 3), 8: (1500, 7, 5), 10: (1500, 9, 7),
       12: (1500, 10, 9), 14: (1500, 12, 10)}
for n, (c, a, b) in CLS.items():
    P("  %-6d %-8d %-10d %-10d" % (n, c, a, b))
xs2 = np.array(list(CLS.keys()), float)
ys2 = np.array([CLS[n][1] for n in CLS], float)
P("  CODE class-level d90 vs n slope = %+.4f  -> %s" %
  (np.polyfit(xs2, ys2, 1)[0], "SATURATES" if abs(np.polyfit(xs2, ys2, 1)[0]) < 0.05 else "GROWS WITH n"))
P("")
P("  READ, filled from the numbers above and never in advance:")
P("  Under EVERY relation definition tried, and on the basis-free class object as well, the")
P("  intrinsic dimension GROWS LINEARLY with the number of records.  The only two columns that")
P("  saturate -- (d)GS and (e)GS -- lose their saturation the moment the arbitrary Gram-Schmidt")
P("  basis is replaced by another equally valid basis of records on the SAME carrier.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/08_synthesis.txt", "w").write("\n".join(OUT) + "\n")
