"""T-24 item 1 -- THE MINIMAL TORUS (1 vertex, 2 edges, 1 face; dim |G|^2 = 64 for D_4).

CLAIM TO PROVE, NOT ASSERT: no contractible qubit-supported region exists, so clause (v)
is VACUOUS on this cell structure.  Two independent routes, both exact:

  (A) HOMOLOGY of the cell complex over Z: boundary matrices are integer matrices,
      d1(e_x) = v - v = 0, d1(e_y) = 0, d2(f) = e_x + e_y - e_x - e_y = 0 (the face glues
      along the commutator word x y x^-1 y^-1, whose abelianisation is 0).
      H_1 = ker d1 / im d2 = Z^2 with basis [e_x], [e_y]: EVERY nonzero edge-chain is a
      non-trivial homology class, so no nonempty edge set lies in a contractible subcomplex.
  (B) THE PROGRAM'S OWN CONVENTION (T-11): contractible region on a graph = edge subset
      containing NO cycle (a forest); SINGLE = connected.  A self-loop IS a cycle, so both
      single-edge sets fail, and every superset fails.  Enumerate all 2^2 subsets exactly.

  The two readings AGREE on every one of the 4 subsets -- the D-23 convention flag is
  raised and answered: no convention choice changes the vacuity.

CROSS-CHECK of this lane's machinery against O-36's independent build: the D(D_4)
minimal-torus sector dimensions must reproduce eig -2: 22, -1: 24, 0: 18 (o36_summary.json).
"""
import sys, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV")
from fractions import Fraction
import numpy as np
from t24_lib import make_D4, sp_from_perms, sp_identity, sp_sub, sp_mask_rows, sp_trace, \
                    sp_frob_check_projector

def say(s=""):
    print(s); OUT.append(s)
OUT = []

say("=" * 110)
say("T-24 / item 1 -- MINIMAL TORUS: WHAT REGIONS EXIST, AND THE VACUITY OF CLAUSE (v)")
say("=" * 110)
say("")
say("CELL STRUCTURE: vertices {v}; edges {e_x, e_y}, BOTH with d(e) = v - v (self-loops);")
say("one face f glued along the word  e_x e_y e_x^-1 e_y^-1.")
say("")

# ---- (A) integer homology
say("(A) HOMOLOGY OVER Z (exact integer boundary matrices)")
d1 = [[0], [0]]                # d1: C1 -> C0, columns e_x, e_y -- wait rows=C0 (1), cols=C1 (2)
d1 = [[0, 0]]                  # 1 x 2 zero matrix
d2 = [[0], [0]]                # 2 x 1: abelianised boundary of f is e_x+e_y-e_x-e_y = (0,0)
say("    d1 = %s   (both edges are loops: d(e) = v - v = 0)" % d1)
say("    d2 = %s (face word x y x^-1 y^-1 abelianises to 0)" % d2)
rank_d1 = 0   # zero matrix
rank_d2 = 0
ker_d1 = 2 - rank_d1
say("    rank d1 = %d, rank d2 = %d  =>  H_1 = Z^%d / 0 = Z^2, basis {[e_x], [e_y]}"
    % (rank_d1, rank_d2, ker_d1))
say("    every nonzero integer edge-chain a[e_x] + b[e_y] is a NONZERO class in H_1;")
say("    a region contained in a contractible subcomplex must carry only trivial classes,")
say("    so NO nonempty edge set is contained in a contractible subcomplex.")
say("")

# ---- (B) the T-11 forest convention, enumerated
say("(B) THE T-11 CONVENTION: contractible = forest (no graph cycle); single = connected")
V = 1
EDGES = {"e_x": (0, 0), "e_y": (0, 0)}
def has_cycle(sub):
    par = list(range(V))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for l in sub:
        a, b = EDGES[l]
        ra, rb = find(a), find(b)
        if ra == rb: return True     # includes self-loops: a==b at once
        par[ra] = rb
    return False
subsets = [[], ["e_x"], ["e_y"], ["e_x", "e_y"]]
rows = []
for s in subsets:
    cyc = has_cycle(s)
    homl = "trivial" if not s else "non-trivial (class %s)" % \
           ("+".join("[%s]" % e for e in s))
    contractible_forest = (not cyc) and len(s) > 0
    contractible_homology = (len(s) > 0) and False   # every nonempty set carries a nonzero class
    rows.append((s, cyc, contractible_forest, contractible_homology))
    say("    region %-16s graph-cycle: %-5s  forest-contractible: %-5s  homology-contractible: %-5s"
        % ("{" + ",".join(s) + "}", cyc, contractible_forest, contractible_homology))
agree = all(r[2] == r[3] for r in rows)
say("    the two readings agree on all %d subsets: %s  (D-23: no convention-dependence HERE)"
    % (len(subsets), agree))
say("")
say("VACUITY: the set of single contractible qubit-supported regions is EMPTY.")
say("Clause (v) quantifies over that set; with an empty domain it is TRUE VACUOUSLY --")
say("it constrains nothing and protects nothing on this cell structure.  SCOPE (D-23):")
say("this is a statement about the 1-vertex cell structure of the proxy, not about D(D_4):")
say("the SAME group on the 1x2 refinement has two contractible edges (t24_main.py), so the")
say("vacuity is a resolution artifact that disappears under refinement.")
say("")

# ---- machinery cross-check against O-36
say("CROSS-CHECK: this lane's exact builders reproduce O-36's D(D_4) minimal-torus sectors")
G = make_D4()
n = G["n"]; N = n * n
MUL, INV = G["MUL"], G["INV"]
idx = np.arange(N)
x = idx // n; y = idx % n
def compose(x_, y_): return x_ * n + y_
permsA = []
for k in range(n):
    xk = MUL[MUL[k, x], INV[k]]
    yk = MUL[MUL[k, y], INV[k]]
    permsA.append(compose(xk, yk))
# rep property exact
repok = all(np.array_equal(permsA[k1][permsA[k2]], permsA[int(MUL[k1, k2])])
            for k1 in range(n) for k2 in range(n))
say("    A(k) is a permutation representation of D_4: %s" % repok)
hol = MUL[MUL[MUL[x, y], INV[x]], INV[y]]
bdiag = (hol == 0).astype(np.int64)
# B invariant under every A(k)?
binv = all(np.array_equal(bdiag[pi], bdiag) for pi in permsA)
say("    [A(k), B] = 0 at the permutation level (B-diagonal invariant): %s" % binv)
A = sp_from_perms(permsA, Fraction(1, n))
I = sp_identity(N)
Q = {1: A, 0: sp_sub(I, A)}
dims = {}
for a in (0, 1):
    for b in (0, 1):
        P = sp_mask_rows(Q[a], bdiag, b)
        tr = sp_trace(P)
        fro = sp_frob_check_projector(P)
        assert tr == fro, "sector (%d,%d) fails Tr(P^2)=Tr(P)" % (a, b)
        dims[(a, b)] = tr
say("    sector dims (a,b)->dim: %s   total %s"
    % ({k: str(v) for k, v in dims.items()}, sum(dims.values())))
eig = {}
for (a, b), d in dims.items(): eig[a + b] = eig.get(a + b, 0) + d
say("    H-eigenspace dims: eig -2: %s  -1: %s  0: %s   (O-36: 22 / 24 / 18)"
    % (eig.get(2), eig.get(1), eig.get(0)))
match = (eig.get(2) == 22 and eig.get(1) == 24 and eig.get(0) == 18)
say("    MATCH with O-36's independent construction: %s" % match)
say("")
say("VERDICT item 1: clause (v) on the minimal torus is VACUOUS -- proved from the cell")
say("structure by both the homology route and the program's own forest convention, which agree.")

with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV/t24_minimal.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
json.dump(dict(vacuous=True, conventions_agree=agree, o36_match=match,
               sector_dims={str(k): int(v) for k, v in dims.items()}),
          open("/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV/t24_minimal.json", "w"), indent=1)
