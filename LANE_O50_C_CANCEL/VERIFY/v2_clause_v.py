"""
VERIFY 2.  THE CLAUSE-(v) ATTACK, AND THE 'H1 IS A MEASUREMENT' ATTACK.

(1) The lane's own stated criterion for a CONTRACTIBLE region is, verbatim from s1_carrier.py:
      "a region is the edge set of a d x d block of sites with d < L, so the block
       lifts to the plane and is CONTRACTIBLE."
    The operative property is 'lifts to the plane'.  Any subgraph of the torus lattice that is a
    TREE lifts to the plane (its fundamental group is trivial), so it satisfies the lane's own
    criterion.  This file builds such a tree and asks the LANE'S OWN region test whether it
    contains a non-trivial logical.

(2) The lane reports 'rank of writer map = m at every (L,k)' as a MEASURED property of the torus.
    This file shows the rank is m for ANY code and ANY independent commuting record family --
    it is the definition of 'independently writable', not a fact about the toric code.
"""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

def line(s=""): print(s, flush=True)

# ---------------------------------------------------------------- toric code, rebuilt independently
class Toric:
    def __init__(self, L):
        self.L = L; self.n = 2 * L * L
    def h(self, i, j): L = self.L; return (i % L) * L + (j % L)
    def v(self, i, j): L = self.L; return L * L + (i % L) * L + (j % L)
    def A(self, i, j): return [self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)]
    def B(self, i, j): return [self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)]
    def stabilisers(self):
        n = self.n; rows = []
        for i in range(self.L):
            for j in range(self.L):
                r = [0] * (2 * n)
                for e in self.A(i, j): r[e] ^= 1
                rows.append(r)
        for i in range(self.L):
            for j in range(self.L):
                r = [0] * (2 * n)
                for e in self.B(i, j): r[n + e] ^= 1
                rows.append(r)
        return rows
    def endpoints(self, e):
        L = self.L
        if e < L * L: i, j = divmod(e, L); return ((i, j), (i, (j + 1) % L))
        e -= L * L; i, j = divmod(e, L); return ((i, j), ((i + 1) % L, j))

def sp(a, b, n): return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def rref(rows, width):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(width):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

# ---- the LANE'S OWN region test, copied verbatim from s1_carrier.py (region_logical_dim)
def region_logical_dim(T, S, edges):
    n = T.n
    cols = []
    for e in edges:
        vx = [0] * (2 * n); vx[e] = 1; cols.append(vx)
        vz = [0] * (2 * n); vz[n + e] = 1; cols.append(vz)
    Srows, _ = rref(S, 2 * n)
    M = [[sp(cols[b], Srows[a], n) for b in range(len(cols))] for a in range(len(Srows))]
    Mr, piv = rref(M, len(cols))
    free = [c for c in range(len(cols)) if c not in piv]
    NS = []
    for f in free:
        coef = [0] * len(cols); coef[f] = 1
        for i, c in enumerate(piv): coef[c] = Mr[i][f]
        v = [0] * (2 * n)
        for b, cf in enumerate(coef):
            if cf: v = [(x + y) % 2 for x, y in zip(v, cols[b])]
        if any(v): NS.append(v)
    NSr, _ = rref(NS, 2 * n)
    dimN = len(NSr)
    both = rref(NS + Srows, 2 * n)[0]
    dimS_in_region = dimN + len(Srows) - len(both)
    return dimN, dimS_in_region, dimN - dimS_in_region

def is_tree(T, edges):
    """the subgraph of the torus lattice on `edges` is a forest with one component => a TREE,
       hence simply connected, hence LIFTS TO THE PLANE (the lane's own contractibility test)."""
    parent = {}
    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    verts = set()
    for e in edges:
        u, w = T.endpoints(e); verts.add(u); verts.add(w)
    for e in edges:
        u, w = T.endpoints(e)
        ru, rw = find(u), find(w)
        if ru == rw: return False, "contains a cycle"
        parent[ru] = rw
    comps = len({find(x) for x in verts})
    return (comps == 1), f"{len(verts)} vertices, {len(edges)} edges, {comps} component(s)"

line("=" * 100)
line("V2-A.  A CONTRACTIBLE REGION ON THE TORUS THAT CONTAINS A LOGICAL  --  CLAUSE (v) ATTACK")
line("=" * 100)
line("  Region built: the vertical path v(0,0)..v(L-2,0) joining sites (0,0)->(L-1,0), PLUS the")
line("  horizontal stub h(i,0) at every row i.  This is a TREE, so it is simply connected and")
line("  lifts to the plane -- the LANE'S OWN definition of contractible.  It is NOT a d x d block,")
line("  which is the only shape the lane tested.")
line()
line(f"  {'L':>3} {'#edges':>7} {'tree?':>6} {'shape':>34} {'dimN':>5} {'dimS':>5} {'NON-TRIVIAL LOGICALS':>21} {'clause (v)':>11}")
for L in (3, 4, 5):
    T = Toric(L); S = T.stabilisers()
    edges = [T.v(i, 0) for i in range(L - 1)] + [T.h(i, 0) for i in range(L)]
    edges = sorted(set(edges))
    tree, why = is_tree(T, edges)
    a, b, c = region_logical_dim(T, S, edges)
    line(f"  {L:>3} {len(edges):>7} {str(tree):>6} {why:>34} {a:>5} {b:>5} {c:>21} "
         f"{('holds' if c == 0 else 'FAILS'):>11}")
line()
line("  NEGATIVE CONTROL (D-15): the same tree with ONE horizontal stub removed, so the X-loop is")
line("  incomplete.  If the instrument is sound this must return 0 logicals.")
line(f"  {'L':>3} {'#edges':>7} {'tree?':>6} {'NON-TRIVIAL LOGICALS':>21}")
for L in (3, 4, 5):
    T = Toric(L); S = T.stabilisers()
    edges = sorted(set([T.v(i, 0) for i in range(L - 1)] + [T.h(i, 0) for i in range(L - 1)]))
    tree, why = is_tree(T, edges)
    a, b, c = region_logical_dim(T, S, edges)
    line(f"  {L:>3} {len(edges):>7} {str(tree):>6} {c:>21}")
line()
line("  AND THE LANE'S OWN d<L BLOCK ROWS, REPRODUCED, so the two are read in the same table:")
def region_edges(T, i0, j0, d):
    E = []
    for i in range(i0, i0 + d):
        for j in range(j0, j0 + d):
            if j + 1 < j0 + d: E.append(T.h(i, j))
            if i + 1 < i0 + d: E.append(T.v(i, j))
    return sorted(set(E))
line(f"  {'L':>3} {'d':>3} {'#edges':>7} {'NON-TRIVIAL LOGICALS':>21}")
for L in (3, 4):
    T = Toric(L); S = T.stabilisers()
    for d in range(2, L):
        worst = max(region_logical_dim(T, S, region_edges(T, i0, j0, d))[2]
                    for i0 in range(L) for j0 in range(L))
        line(f"  {L:>3} {d:>3} {len(region_edges(T,0,0,d)):>7} {worst:>21}")
line()
line("  DIAMETER OF THE TREE REGION vs THE BLOCKS (the quantity that actually separates them):")
for L in (3, 4, 5):
    line(f"    L={L}:  tree region spans all {L} rows -> diam = L ;  lane's blocks have diam <= d < L")

# ---------------------------------------------------------------- B. H1 is definitional
line()
line("=" * 100)
line("V2-B.  IS 'rank of writer map = m' A MEASUREMENT ABOUT THE TORUS?")
line("=" * 100)
line("  The lane takes recs = [pair[0] for pair in symplectic_logicals(...)] and gens = both halves")
line("  of every pair.  symplectic_logicals RETURNS CONJUGATE PAIRS, i.e. a symplectic basis, so")
line("  sp(pair[j][1], pair[i][0]) = delta_ij BY CONSTRUCTION.  The rank is then m for any code.")
line()
line(f"  {'code':<40} {'k':>3} {'m':>3} {'rank of w->flip map':>20} {'simply transitive?':>19}")

def rank_of_flip_map(S, n, label):
    pairs = symplectic_logicals(S, n)
    recs = [p[0] for p in pairs]
    gens = [a for p in pairs for a in p]
    img = [[sp(g, R, n) for R in recs] for g in gens]
    r = len(rref(img, len(recs))[0])
    line(f"  {label:<40} {len(pairs):>3} {len(recs):>3} {r:>20} {str(r == len(recs)):>19}")
    return r

# (a) the toric code L=2
T = Toric(2); rank_of_flip_map(T.stabilisers(), T.n, "toric code L=2 (TORUS)")
T = Toric(3); rank_of_flip_map(T.stabilisers(), T.n, "toric code L=3 (TORUS)")
# (b) a code with NO topology at all: 4 bare qubits, no stabilisers -> k = 4 logical qubits
n = 4; rank_of_flip_map([], n, "4 BARE QUBITS, no stabilisers, no code")
# (c) the 5-qubit perfect code [[5,1,3]] x 4 copies -> k = 4
def five_qubit_stabs(offset, ntot):
    gens = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    rows = []
    for g in gens:
        r = [0] * (2 * ntot)
        for i, ch in enumerate(g):
            q = offset + i
            if ch in "XY": r[q] ^= 1
            if ch in "ZY": r[ntot + q] ^= 1
        rows.append(r)
    return rows
ntot = 20; S = []
for c in range(4): S += five_qubit_stabs(5 * c, ntot)
rank_of_flip_map(S, ntot, "four copies of the [[5,1,3]] perfect code")
# (d) a repetition-code-like stabiliser group on 6 qubits, k=5
ntot = 6
S = [[0]*(2*ntot)]
S[0][ntot+0] = 1; S[0][ntot+1] = 1     # Z0 Z1 only
rank_of_flip_map(S, ntot, "one stabiliser Z0Z1 on 6 qubits (k=5)")
line()
line("  READ: the rank is m for every one of these, including a carrier with NO topology, NO")
line("  protection and NO geometry (4 bare qubits).  'rank = m' is the DEFINITION of an")
line("  independent commuting record family, not a measured property of the toric code.")
