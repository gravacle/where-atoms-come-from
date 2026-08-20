"""LANE T48_A -- THE STRUCTURAL DERIVATION (O-58 N2, probe 1), machine-checked.

QUESTION: does the record surface's writer algebra FORCE the induced one-step kernel on
configurations to be measure-conserving (doubly stochastic)?  And does measure
conservation plus the venue's own symmetry force the per-link amplitude to 1/deg -- the
criticality value mu_c located by T-44 -- or does uniformity need MORE?

CANDIDATE ROUTE UNDER TEST (refuted where it fails, kept where it holds):
  (R-a) admissible writers are unitaries (O-4 corner, DEF-A) or energy-conserving
        dilation unitaries (iv'); an ensemble of unitaries is a unital, trace-preserving
        channel; a unital channel's induced transition matrix on the configuration basis
        is doubly stochastic (counting-measure conserving).  [Owner of the standard
        facts: unital/doubly-stochastic correspondence and Birkhoff's theorem, standard
        quantum information / linear algebra.  Machine-checked here on explicit writers.]
  (R-b) on a deg-regular venue with the venue's own automorphism symmetry, double
        stochasticity forces per-link amplitude 1/deg.  [This HALF IS FALSE as stated --
        the counterexample is computed in S4; what more is needed is named and computed.]

CONVENTIONS (fixed before measurement):
  * Kernel convention: K[x][y] = one-step amplitude x -> y.  Row sum = forward
    (probability) normalization; column sum = counting-measure conservation.
    DOUBLY STOCHASTIC := every row sum == 1 AND every column sum == 1 (exact).
  * Exact arithmetic ONLY: python ints and Fractions on every path.  This lane uses NO
    floats anywhere -- there is no comparison gate that needs one.
  * Every verdict is a computed boolean.  Controls that must FAIL a gate are run through
    the SAME gate function and scored on the computed False.
  * deg is COMPUTED from each venue's own row sums (gated constant), never declared;
    1/deg is formed from the computed deg.

DECLARED CONSTANTS (all rational, fixed here):
  N1D        = 8                     -- Ising-ring carrier sites (D=1)
  L2D        = 3, L2D_B = 4          -- toric carrier sides (D=2); full syndrome
                                        channel at L=3 (2^8 = 256 even syndromes)
  L3D        = 3                     -- Z_3^3 grain venue (D=3, venue-level)
  LEAK_S     = 9/10                  -- declared survival bias of the leaky ensemble
  LAZY_C     = 0, 1/5, 1/3, 1/2      -- declared lazy/trivial-writer diagonal rows
  X_ROWS     = 1/3, 1/2, 9/10        -- resolvent dial rows for the resummation identity
  BIAS_2D    = a=1/3 (x-links), b=1/6 (y-links)          -- doubly-stochastic, non-uniform
  BIAS_3D    = 1/4, 1/6, 1/12 per axis                   -- doubly-stochastic, non-uniform
  BATH_P     = 1/2, 3/4, 1           -- bath ground-state weights for the dilation control
  W-ROWS     = three declared non-uniform rational weight lists per writer set (below)
  ZSAMPLE    = z_k = (k * 2654435761) mod 2^18, k = 0..511   -- declared state sample

INSTRUMENT REUSE: toric_stabilizers(L) from model/geometry.py (C-78 lineage; Kitaev
carrier, quant-ph/9707021 -- owner named there).  Plaquette adjacency is rebuilt here
from the plaquette supports alone, the T-44-A convention (one walk step = one shared
carrier edge = one grain-boundary crossing = one unit of writer weight, the Gamma price).

RELEVANCE TEST (borrowed idea -> named program variable, before use): the borrowed
unital-channel / doubly-stochastic correspondence is applied to the named variable
T[x][y] -- the one-step writer kernel on record configurations, whose per-link entry is
the INDUCED mu that O-58 N2 asks for -- and to nothing else.  Gershgorin (borrowed) is
applied only to the named nonsingularity witness for I - LEAK_S*T.  Perron (borrowed,
already the T-44 owner) is applied only to the criticality reading of row sums.

RETURNS ARE DATA: t48a_derivation.RESULT.json.  Output: gate lines, computed booleans.
"""
import json
import sys
from fractions import Fraction as Fr

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from geometry import toric_stabilizers  # noqa: E402  (C-78 lineage instrument)

F0, F1 = Fr(0), Fr(1)

# ---------------------------------------------------------------------------- gate kit
GATES = []


def gate(name, ok):
    ok = bool(ok)
    GATES.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name)
    return ok


def section(title):
    print()
    print("=" * 88)
    print(title)
    print("=" * 88)


# ------------------------------------------------------------------------- exact matrix
def zeros(n):
    return [[F0] * n for _ in range(n)]


def eye(n):
    m = zeros(n)
    for i in range(n):
        m[i][i] = F1
    return m


def madd(A, B):
    n = len(A)
    return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]


def msca(c, A):
    return [[c * x for x in row] for row in A]


def msub(A, B):
    n = len(A)
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def meq(A, B):
    n = len(A)
    return all(A[i][j] == B[i][j] for i in range(n) for j in range(n))


def mvec(A, v):
    n = len(A)
    return [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]


def row_sums(K):
    return [sum(row) for row in K]


def col_sums(K):
    n = len(K)
    return [sum(K[i][j] for i in range(n)) for j in range(n)]


def is_stochastic(K):
    return all(s == F1 for s in row_sums(K))


def is_doubly_stochastic(K):
    return is_stochastic(K) and all(s == F1 for s in col_sums(K))


def det_exact(M):
    """Exact determinant by fraction Gaussian elimination (partial pivot on nonzero)."""
    n = len(M)
    A = [row[:] for row in M]
    det = F1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return F0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = F1 / A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] * inv
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return det


# --------------------------------------------------------------------------- union-find
class UF:
    def __init__(self, items):
        self.p = {x: x for x in items}

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb

    def classes(self, items):
        return len({self.find(x) for x in items})


def orbit_closure(start, gens):
    """BFS closure of `start` (a hashable) under callables gens."""
    seen = {start}
    frontier = [start]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = g(x)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return seen


# ============================================================================ S0
section("S0  THE UNITALITY LEMMA ON EXPLICIT OPERATORS (owner-attributed, machine-checked)")

# S0.1 -- convex mixtures of permutations (non-involutive included) are doubly stochastic.
n0 = 6
perm_r = [(i + 1) % n0 for i in range(n0)]          # 6-cycle (non-involutive)
perm_m = [(n0 - i) % n0 for i in range(n0)]         # reflection
perm_t = [1, 0] + list(range(2, n0))                # transposition
W0_ROWS = [[Fr(1, 2), Fr(1, 3), Fr(1, 6)],
           [Fr(1, 6), Fr(1, 6), Fr(2, 3)],
           [F1, F0, F0]]
ok_all = True
for wi, w in enumerate(W0_ROWS):
    T = zeros(n0)
    for wgt, perm in zip(w, [perm_r, perm_m, perm_t]):
        for x in range(n0):
            T[x][perm[x]] += wgt
    ok_all &= is_doubly_stochastic(T)
gate("S0.1 every declared permutation-ensemble kernel is doubly stochastic (3 weight rows)", ok_all)

# S0.2 -- a rational-orthogonal (non-permutation) unitary induces a doubly stochastic kernel.
O = [[Fr(3, 5), Fr(4, 5)], [Fr(-4, 5), Fr(3, 5)]]
OtO = [[sum(O[k][i] * O[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
gate("S0.2a rational orthogonal writer: O^T O == I exactly", meq(OtO, eye(2)))
T_O = [[O[i][j] ** 2 for j in range(2)] for i in range(2)]
gate("S0.2b induced kernel |O[i][j]|^2 is doubly stochastic (unitality, not permutation-ness)",
     is_doubly_stochastic(T_O))

# S0.3 -- CONTROL: a trace-preserving NON-unitary (reset Kraus) channel is NOT doubly stochastic.
# Kraus K1=|0><0|, K2=|0><1|; convention kept: K[x][y] = x -> y.  Reset sends every x to 0.
T_reset = [[F1, F0], [F1, F0]]
gate("S0.3a reset channel is stochastic (probability-preserving forward): all row sums == 1",
     is_stochastic(T_reset))
gate("S0.3b reset channel FAILS double stochasticity (computed False == the control's job)",
     not is_doubly_stochastic(T_reset))

# ============================================================================ S1
section("S1  CARRIERS AND VENUES -- built from supports, energy sectors computed")

VENUES = {}  # name -> dict(adj=list[set], deg=int, n=int)

# ---- D=1: Ising ring, N1D sites, bonds b=(b,b+1); H = -sum Z_b Z_{b+1}; walls on bonds.
N1D = 8
bond_sites = [(b, (b + 1) % N1D) for b in range(N1D)]
site_bonds = [[] for _ in range(N1D)]
for b, (s1, s2) in enumerate(bond_sites):
    site_bonds[s1].append(b)
    site_bonds[s2].append(b)
gate("S1.1 D=1 carrier: every site lies in exactly 2 bonds (boundary map well-formed)",
     all(len(bs) == 2 for bs in site_bonds))
DEL1 = [frozenset(site_bonds[i]) for i in range(N1D)]  # del(X_i): the two bonds at site i

adj1 = [set() for _ in range(N1D)]
for i in range(N1D):
    a, b = sorted(DEL1[i])
    adj1[a].add(b)
    adj1[b].add(a)
gate("S1.2 D=1 venue (bond graph through shared sites) is 2-regular on 8 vertices",
     all(len(x) == 2 for x in adj1) and len(adj1) == 8)
VENUES["C8"] = dict(adj=adj1, deg=2, n=N1D)


def syndrome1(z):
    return sum(((z >> b) & 1) ^ ((z >> ((b + 1) % N1D)) & 1) and (1 << b) for b in range(N1D))


def maskset1(fs):
    m = 0
    for b in fs:
        m |= 1 << b
    return m


ok = True
for z in range(1 << N1D):
    sz = syndrome1(z)
    for i in range(N1D):
        ok &= syndrome1(z ^ (1 << i)) == (sz ^ maskset1(DEL1[i]))
gate("S1.3 D=1 induced writer action grounded in the carrier: syndrome(X_i z) == "
     "syndrome(z) XOR del(i) for ALL 256 states x 8 writers", ok)

ok = True
moving = 0
for z in range(1 << N1D):
    sz = syndrome1(z)
    for i in range(N1D):
        dE = 2 * (bin(sz ^ maskset1(DEL1[i])).count("1") - bin(sz).count("1"))
        one_adj = bin(sz & maskset1(DEL1[i])).count("1") == 1
        ok &= (dE == 0) == one_adj and dE in (-4, 0, 4)
        moving += 1 if dE == 0 else 0
gate("S1.4 D=1 energy sectors: dE == 0 exactly on the defect-MOVING sector, "
     "dE == +/-4 on create/annihilate (all states x writers)", ok)
print(f"DATA  S1.4 moving (dE=0) cases: {moving} of {(1 << N1D) * N1D} (positive count beside the zero)")

# ---- D=2: toric carrier at L=3 and L=4 (reused instrument), plaquette venue from supports.
TOR = {}
for L in (3, 4):
    n_e, _stab, stars, plaqs = toric_stabilizers(L)
    edge_plaqs = [[p for p in range(len(plaqs)) if (plaqs[p] >> e) & 1] for e in range(n_e)]
    gate(f"S1.5 D=2 L={L}: every carrier edge lies in exactly 2 plaquettes "
         f"({n_e} edges checked)", all(len(ep) == 2 for ep in edge_plaqs))
    npl = len(plaqs)
    adj = [set() for _ in range(npl)]
    shared = {}
    for e, (p, q) in enumerate(edge_plaqs):
        adj[p].add(q)
        adj[q].add(p)
        shared[(p, q)] = shared.get((p, q), 0) + 1
    gate(f"S1.6 D=2 L={L}: adjacent plaquettes share exactly ONE carrier edge "
         "(no doubled links at this side)", all(v == 1 for v in shared.values()))
    gate(f"S1.7 D=2 L={L}: plaquette venue is 4-regular on {npl} vertices",
         all(len(x) == 4 for x in adj))
    TOR[L] = dict(n_e=n_e, stars=stars, plaqs=plaqs, edge_plaqs=edge_plaqs, adj=adj)
VENUES["T3"] = dict(adj=TOR[3]["adj"], deg=4, n=9)
VENUES["T4"] = dict(adj=TOR[4]["adj"], deg=4, n=16)

ok = all(bin(s & p).count("1") % 2 == 0 for s in TOR[3]["stars"] for p in TOR[3]["plaqs"])
gate("S1.8 D=2 L=3: every star shares an EVEN number of edges with every plaquette "
     "(81 pairs) -- star writers act trivially on the plaquette syndrome", ok)

# Induced action of X_e on plaquette syndromes, grounded in the Z-basis parity identity:
# b_p(z) = (-1)^{|z AND supp(p)|};  b_p(z XOR e) = b_p(z) * (-1)^{[e in supp(p)]}.
ZSAMPLE = [(k * 2654435761) % (1 << 18) for k in range(512)]
ok = True
for z in ZSAMPLE:
    for e in range(TOR[3]["n_e"]):
        for p in range(9):
            pm = TOR[3]["plaqs"][p]
            lhs = bin((z ^ (1 << e)) & pm).count("1") % 2
            rhs = (bin(z & pm).count("1") % 2) ^ ((pm >> e) & 1)
            ok &= lhs == rhs
gate("S1.9 D=2 L=3: induced action grounded in the Z-basis -- parity identity verified on "
     "512 declared states x 18 writers x 9 plaquettes", ok)
DEL2 = []
for e in range(TOR[3]["n_e"]):
    p, q = TOR[3]["edge_plaqs"][e]
    DEL2.append((1 << p) | (1 << q))
gate("S1.10 D=2 L=3: the 18 writer boundaries del(e) are pairwise distinct "
     "(distinct writers move distinct plaquette pairs)", len(set(DEL2)) == 18)

ok = True
moving2 = 0
for s in range(1 << 9):
    if bin(s).count("1") % 2:
        continue
    for e in range(18):
        d = DEL2[e]
        dE = 2 * (bin(s ^ d).count("1") - bin(s).count("1"))
        one = bin(s & d).count("1") == 1
        ok &= (dE == 0) == one and dE in (-4, 0, 4)
        moving2 += 1 if dE == 0 else 0
gate("S1.11 D=2 L=3 energy sectors under H = -sum B_p: dE == 0 exactly on the moving "
     "sector, +/-4 on create/annihilate (256 even syndromes x 18 writers)", ok)
print(f"DATA  S1.11 moving (dE=0) cases: {moving2} of {256 * 18}")

# ---- D=3: Z_3^3 grain venue (face adjacency), venue-level (T44-B world lineage).
L3 = 3
g_idx = {(i, j, k): i * 9 + j * 3 + k for i in range(3) for j in range(3) for k in range(3)}
g_lst = sorted(g_idx, key=g_idx.get)
adj3 = [set() for _ in range(27)]
for (i, j, k), gi in g_idx.items():
    for ax, dd in ((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1)):
        c = [i, j, k]
        c[ax] = (c[ax] + dd) % 3
        adj3[gi].add(g_idx[tuple(c)])
gate("S1.12 D=3 venue: grain face-adjacency is 6-regular on 27 grains", all(len(x) == 6 for x in adj3))
VENUES["Z27"] = dict(adj=adj3, deg=6, n=27)

# deg is COMPUTED per venue (gated constant); 1/deg formed from it.
for nm, V in VENUES.items():
    degs = {len(x) for x in V["adj"]}
    gate(f"S1.13 venue {nm}: computed degree is the constant {V['deg']}", degs == {V["deg"]})

# ============================================================================ S2
section("S2  THE VENUE'S OWN SYMMETRY -- automorphisms, edge-transitivity, carrier lifts")

AUTS = {}


def check_aut(adj, gmap):
    n = len(adj)
    return all(sorted(gmap[y] for y in adj[x]) == sorted(adj[gmap[x]]) for x in range(n))


# C8 generators: rotation (from site rotation), reflection.
rot1 = [(b + 1) % 8 for b in range(8)]
ref1 = [(-b - 1) % 8 for b in range(8)]  # site i -> -i maps bond (b,b+1) -> bond -b-1
AUTS["C8"] = [rot1, ref1]

# T3 / T4 generators: two translations, axis swap, reflection (on plaquette cells (i,j)).
for L, nm in ((3, "T3"), (4, "T4")):
    def c2i(i, j, L=L):
        return (i % L) * L + (j % L)
    gens = []
    for f in (lambda i, j: (i + 1, j), lambda i, j: (i, j + 1),
              lambda i, j: (j, i), lambda i, j: (-i, j)):
        gens.append([c2i(*f(i, j)) for i in range(L) for j in range(L)])
    AUTS[nm] = gens

# Z27 generators: three translations, two coordinate transpositions, one reflection.
gens3 = []
for f in (lambda i, j, k: (i + 1, j, k), lambda i, j, k: (i, j + 1, k),
          lambda i, j, k: (i, j, k + 1), lambda i, j, k: (j, i, k),
          lambda i, j, k: (i, k, j), lambda i, j, k: (-i, j, k)):
    gens3.append([g_idx[tuple(x % 3 for x in f(*g))] for g in g_lst])
AUTS["Z27"] = gens3

for nm, V in VENUES.items():
    adj = V["adj"]
    gate(f"S2.1 {nm}: every declared generator preserves adjacency exactly",
         all(check_aut(adj, g) for g in AUTS[nm]))
    vg = [lambda x, g=g: g[x] for g in AUTS[nm]]
    gate(f"S2.2 {nm}: vertex orbit of vertex 0 is ALL {V['n']} vertices (vertex-transitive)",
         len(orbit_closure(0, vg)) == V["n"])
    e0 = (0, sorted(adj[0])[0])
    eg = [lambda e, g=g: (g[e[0]], g[e[1]]) for g in AUTS[nm]]
    ndir = sum(len(x) for x in adj)
    gate(f"S2.3 {nm}: directed-edge orbit of one declared edge is ALL {ndir} directed edges "
         "(edge-transitive: the elementary-writer set is ONE orbit)",
         len(orbit_closure(e0, eg)) == ndir)

# Carrier lifts (D=2, L=3): each venue generator lifts to an edge permutation of the
# carrier that maps plaquette supports onto plaquette supports and star supports onto
# star supports -- the symmetry is the RECORD SURFACE's own, not only the derived graph's.
plaq_set3 = set(TOR[3]["plaqs"])
star_set3 = set(TOR[3]["stars"])
for gi, g in enumerate(AUTS["T3"]):
    lift = {}
    ok_unique = True
    for e in range(18):
        p, q = TOR[3]["edge_plaqs"][e]
        common = TOR[3]["plaqs"][g[p]] & TOR[3]["plaqs"][g[q]]
        ok_unique &= bin(common).count("1") == 1
        lift[e] = common.bit_length() - 1
    ok_bij = sorted(lift.values()) == list(range(18))

    def lift_mask(m):
        out = 0
        for e in range(18):
            if (m >> e) & 1:
                out |= 1 << lift[e]
        return out
    ok_p = all(lift_mask(TOR[3]["plaqs"][p]) == TOR[3]["plaqs"][g[p]] for p in range(9))
    ok_s = {lift_mask(s) for s in star_set3} == star_set3
    gate(f"S2.4 T3 generator {gi}: lifts to a carrier edge bijection mapping plaquette "
         "supports to plaquette supports and star supports onto star supports",
         ok_unique and ok_bij and ok_p and ok_s)

# ============================================================================ S3
section("S3  INDUCED KERNELS OF THE HONEST ENSEMBLE CONSTRUCTIONS (all of them)")

RES = {"venues": {nm: {"n": V["n"], "deg": V["deg"]} for nm, V in VENUES.items()},
       "ensembles": {}, "controls": {}}


def kernel_uniform(V, c):
    """c*I + ((1-c)/deg)*A -- the invariant lazy family; c=0 is E-LOC."""
    n, deg, adj = V["n"], V["deg"], V["adj"]
    t = (F1 - c) / deg
    K = zeros(n)
    for x in range(n):
        K[x][x] = c
        for y in adj[x]:
            K[x][y] = t
    return K


def link_amplitudes(V, K):
    return {K[x][y] for x in range(V["n"]) for y in V["adj"][x]}


def crit_witness(V, K):
    """Exact singularity of I - K: det == 0 (venue-sized, exact)."""
    return det_exact(msub(eye(V["n"]), K))


# ---- E-LOC: local ensemble, uniform over the deg elementary writers at the position.
for nm, V in VENUES.items():
    K = kernel_uniform(V, F0)
    dsK = is_doubly_stochastic(K)
    amps = link_amplitudes(V, K)
    uni = amps == {Fr(1, V["deg"])}
    d0 = crit_witness(V, K)
    gate(f"S3.1 E-LOC on {nm}: doubly stochastic (measure-conserving)", dsK)
    gate(f"S3.2 E-LOC on {nm}: link-uniform with per-link amplitude == 1/deg == 1/{V['deg']}", uni)
    gate(f"S3.3 E-LOC on {nm}: CRITICAL -- det(I - K) == 0 exactly", d0 == 0)
    dleak = crit_witness(V, msca(Fr(9, 10), K))
    gate(f"S3.4 E-LOC on {nm}: positive control beside the zero -- det(I - (9/10)K) != 0", dleak != 0)
    print(f"DATA  S3.4 {nm}: det(I - (9/10)K) = {dleak}")
    RES["ensembles"][f"E-LOC/{nm}"] = dict(ds=dsK, uniform=uni, per_link=str(Fr(1, V["deg"])),
                                           crit_det="0", leak_det=str(dleak))

# ---- E-LAZY-I: identity-augmented local ensemble (uniform over {1} + incident writers).
for nm, V in VENUES.items():
    deg = V["deg"]
    c = Fr(1, deg + 1)
    K = kernel_uniform(V, c)
    naive = {K[x][y] for x in range(V["n"]) for y in V["adj"][x]}
    gate(f"S3.5 E-LAZY-I on {nm}: doubly stochastic", is_doubly_stochastic(K))
    gate(f"S3.5b E-LAZY-I on {nm}: naive PER-STEP link amplitude == 1/(deg+1) != 1/deg "
         "(the honest ambiguity, registered not hidden)",
         naive == {Fr(1, deg + 1)} and Fr(1, deg + 1) != Fr(1, deg))
    mu_eff = Fr(1, deg + 1) / (F1 - c)
    gate(f"S3.6 E-LAZY-I on {nm}: resummed PER-CROSSING amplitude t/(1-c) == 1/deg exactly "
         "(the Gamma-priced object -- weight counts crossings, not dwell steps)",
         mu_eff == Fr(1, deg))
    gate(f"S3.7 E-LAZY-I on {nm}: still CRITICAL -- det(I - K) == 0 exactly",
         crit_witness(V, K) == 0)
    RES["ensembles"][f"E-LAZY-I/{nm}"] = dict(ds=True, per_step=str(Fr(1, deg + 1)),
                                              per_crossing=str(mu_eff))

# ---- E-GLOB: global uniform ensemble on the FULL configuration (syndrome) space.
def glob_kernel(nsynd_bits, deltas, weights, extra_trivial=0):
    """T[s][s XOR delta] += w on the even-parity syndrome space; trivial writers add
    weight `extra_trivial` to the diagonal.  Returns dict-of-dict."""
    synds = [s for s in range(1 << nsynd_bits) if bin(s).count("1") % 2 == 0]
    T = {s: {} for s in synds}
    for s in synds:
        if extra_trivial:
            T[s][s] = T[s].get(s, F0) + extra_trivial
        for d, w in zip(deltas, weights):
            t = s ^ d
            T[s][t] = T[s].get(t, F0) + w
    return synds, T


def dict_ds(synds, T):
    rows_ok = all(sum(T[s].values()) == F1 for s in synds)
    cols = {s: F0 for s in synds}
    for s in synds:
        for t, w in T[s].items():
            cols[t] += w
    return rows_ok and all(v == F1 for v in cols.values())


D1_DELTAS = [maskset1(d) for d in DEL1]
synds1, T1 = glob_kernel(N1D, D1_DELTAS, [Fr(1, 8)] * 8)
gate("S3.8 E-GLOB on D=1 syndrome space (128 even configs, 8 writers, uniform): "
     "doubly stochastic on the FULL configuration space", dict_ds(synds1, T1))

synds2, T2 = glob_kernel(9, DEL2, [Fr(1, 18)] * 18)
gate("S3.9 E-GLOB on D=2 L=3 syndrome space (256 even configs, 18 writers, uniform): "
     "doubly stochastic on the FULL configuration space", dict_ds(synds2, T2))

# E-GLOB-S: stabilizer-augmented (trivial-action writers included from the ALGEBRA:
# the 9 stars on D=2 -- S1.8's even-overlap fact; the global flip on D=1).
synds1s, T1s = glob_kernel(N1D, D1_DELTAS, [Fr(1, 9)] * 8, extra_trivial=Fr(1, 9))
gate("S3.10 E-GLOB-S on D=1 (8 movers + global flip, uniform over 9): doubly stochastic; "
     "diagonal weight == 1/9 from the trivial-action writer", dict_ds(synds1s, T1s)
     and all(T1s[s].get(s, F0) == Fr(1, 9) for s in synds1s))
synds2s, T2s = glob_kernel(9, DEL2, [Fr(1, 27)] * 18, extra_trivial=Fr(9, 27))
gate("S3.11 E-GLOB-S on D=2 L=3 (18 movers + 9 stars, uniform over 27): doubly stochastic; "
     "diagonal weight == 9/27 == 1/3 from the trivial-action stabilizer writers",
     dict_ds(synds2s, T2s) and all(T2s[s].get(s, F0) == Fr(1, 3) for s in synds2s))

# Arbitrary DECLARED non-uniform weights: measure conservation needs NO uniformity
# anywhere on this path (the no-circle gate).
W_D1 = [[Fr(i + 1, 36) for i in range(8)],
        [Fr(2 ** i, 255) for i in range(8)],
        [Fr(1, 6) if i % 2 == 0 else Fr(1, 12) for i in range(8)]]
W_D2 = [[Fr(e + 1, 171) for e in range(18)],
        [Fr(2 ** e, (1 << 18) - 1) for e in range(18)],
        [Fr(1, 12) if e % 2 == 0 else Fr(1, 36) for e in range(18)]]
ok = all(sum(w) == F1 and dict_ds(*glob_kernel(N1D, D1_DELTAS, w)) for w in W_D1)
ok &= all(sum(w) == F1 and dict_ds(*glob_kernel(9, DEL2, w)) for w in W_D2)
gate("S3.12 EVERY declared non-uniform weight row (3 on D=1, 3 on D=2) still induces a "
     "doubly stochastic kernel: measure conservation comes from the writer ALGEBRA "
     "(involutive permutation action), never from uniform weights", ok)

# Pair-sector reading of E-GLOB (the T-44 connecting-string model's per-crossing mu):
# tracked endpoint p = plaquette (1,1), held partner q0 = plaquette (0,0) (non-adjacent).
p_track = 1 * 3 + 1
q_hold = 0
gate("S3.13 declared pair (p=(1,1), q0=(0,0)) is NON-adjacent on T3 (no annihilation moves)",
     p_track not in VENUES["T3"]["adj"][q_hold])
s_pair = (1 << p_track) | (1 << q_hold)
row = T2[s_pair]
tracked, origin_side, creation = {}, {}, F0
for t, w in row.items():
    if bin(t).count("1") == 2 and (t >> q_hold) & 1:
        pt = (t & ~(1 << q_hold)).bit_length() - 1
        if pt != p_track:
            tracked[pt] = w                      # tracked end moved, origin held
    elif bin(t).count("1") == 2 and (t >> p_track) & 1:
        qt = (t & ~(1 << p_track)).bit_length() - 1
        origin_side[qt] = w                      # origin end moved (string growth there)
    else:
        creation += w                            # pair created: 4-defect config
nbrs = VENUES["T3"]["adj"][p_track]
amps_equal = set(tracked) == nbrs and len(set(tracked.values())) == 1
per_link_glob = next(iter(tracked.values()))
act_at_p = sum(tracked.values())
cond = {k: v / act_at_p for k, v in tracked.items()}
gate("S3.14 E-GLOB pair sector: the tracked endpoint's 4 link amplitudes are EQUAL "
     "(== 1/18 each, the global normalization)", amps_equal and per_link_glob == Fr(1, 18))
gate("S3.15 E-GLOB pair sector: CONDITIONED on the step acting at the tracked endpoint, "
     "the per-link amplitude == 1/deg == 1/4 exactly (the connecting-string per-crossing mu)",
     set(cond.values()) == {Fr(1, 4)})
gate("S3.16 E-GLOB pair sector, full step decomposition computed: tracked-end motion "
     "4/18 + origin-end motion 4/18 + PAIR CREATION 10/18 == 1 -- the creation share "
     "leaves the 2-defect sector (the DISCONNECTED pieces the declared H1 string model "
     "does not price; registered, not hidden)",
     sum(tracked.values()) == Fr(4, 18) and sum(origin_side.values()) == Fr(4, 18)
     and creation == Fr(10, 18)
     and sum(tracked.values()) + sum(origin_side.values()) + creation == F1)
RES["ensembles"]["E-GLOB/T3-pair"] = dict(per_link_global=str(per_link_glob),
                                          per_crossing_conditional="1/4",
                                          step_split="4/18 tracked + 4/18 origin + 10/18 creation")

# ============================================================================ S4
section("S4  THE FORCED FORM UNDER SYMMETRY -- and the counterexample that bounds it")

LAZY_C = [F0, Fr(1, 5), Fr(1, 3), Fr(1, 2)]
X_ROWS = [Fr(1, 3), Fr(1, 2), Fr(9, 10)]

for nm, V in VENUES.items():
    n, deg, adj = V["n"], V["deg"], V["adj"]
    positions = [(x, x) for x in range(n)] + [(x, y) for x in range(n) for y in adj[x]]
    uf = UF(positions)
    for g in AUTS[nm]:
        for (x, y) in positions:
            uf.union((x, y), (g[x], g[y]))
    diag_classes = uf.classes([(x, x) for x in range(n)])
    link_classes = uf.classes([(x, y) for x in range(n) for y in adj[x]])
    gate(f"S4.1 {nm}: under the venue's own automorphisms the kernel entries form exactly "
         f"ONE diagonal class and ONE link class (computed: {diag_classes}, {link_classes}) "
         "-- every invariant kernel supported on I+A is c*I + t*A",
         diag_classes == 1 and link_classes == 1)
    ok_ds = ok_crit = ok_mu = ok_res = True
    for c in LAZY_C:
        t = (F1 - c) / deg
        K = kernel_uniform(V, c)
        ok_ds &= is_doubly_stochastic(K)
        ok_crit &= crit_witness(V, K) == 0
        ok_mu &= (t / (F1 - c)) == Fr(1, deg)
        for x in X_ROWS:
            m = (x * t) / (F1 - x * c)
            A = zeros(n)
            for u in range(n):
                for v in adj[u]:
                    A[u][v] = F1
            lhs = msca(F1 - x * c, msub(eye(n), msca(m, A)))
            rhs = msub(eye(n), msca(x, K))
            ok_res &= meq(lhs, rhs)
    gate(f"S4.2 {nm}: for EVERY declared lazy row c, the forced kernel c*I + ((1-c)/deg)*A "
         "is doubly stochastic AND critical (det(I-K) == 0): laziness never opens the gap", ok_ds and ok_crit)
    gate(f"S4.3 {nm}: resolvent identity (1-xc)(I - m(x)A) == I - xK holds exactly at every "
         "declared dial x, and the per-crossing amplitude t/(1-c) == 1/deg IDENTICALLY for "
         "every lazy row: the induced mu does not depend on the trivial-writer share", ok_mu and ok_res)

# THE COUNTEREXAMPLE: measure conservation does NOT force uniformity.
def biased_kernel_T3():
    """a on x-links, b on y-links of the T3 venue, 2a+2b == 1: doubly stochastic, non-uniform."""
    a, b = Fr(1, 3), Fr(1, 6)
    K = zeros(9)
    for i in range(3):
        for j in range(3):
            x = i * 3 + j
            for di, w in (((1, 0), a), ((-1, 0), a), ((0, 1), b), ((0, -1), b)):
                y = ((i + di[0]) % 3) * 3 + ((j + di[1]) % 3)
                K[x][y] += w
    return K


def biased_kernel_Z27():
    ws = [Fr(1, 4), Fr(1, 6), Fr(1, 12)]
    K = zeros(27)
    for (i, j, k), x in g_idx.items():
        for ax in range(3):
            for dd in (1, -1):
                c = [i, j, k]
                c[ax] = (c[ax] + dd) % 3
                K[x][g_idx[tuple(c)]] += ws[ax]
    return K


for nm, K, Vn in (("T3", biased_kernel_T3(), VENUES["T3"]),
                  ("Z27", biased_kernel_Z27(), VENUES["Z27"])):
    amps = link_amplitudes(Vn, K)
    gate(f"S4.4 {nm} CTRL-BIAS-LINK: doubly stochastic (measure IS conserved)", is_doubly_stochastic(K))
    gate(f"S4.5 {nm} CTRL-BIAS-LINK: link-uniformity FAILS (computed: {len(amps)} distinct "
         "link amplitudes) -- measure conservation alone does NOT force 1/deg", len(amps) > 1)
    gate(f"S4.6 {nm} CTRL-BIAS-LINK: still CRITICAL -- det(I - K) == 0 exactly "
         "(conservation pins the gap shut even without uniformity)", crit_witness(Vn, K) == 0)
    swap_gen = AUTS[nm][2] if nm == "T3" else AUTS[nm][3]  # the axis-transposing generator
    viol = any(K[x][y] != K[swap_gen[x]][swap_gen[y]]
               for x in range(Vn["n"]) for y in Vn["adj"][x])
    gate(f"S4.7 {nm} CTRL-BIAS-LINK: invariance FAILS exactly at the axis-transposing "
         "generator (computed witness) -- the MORE that uniformity needs is ensemble "
         "invariance under the venue's edge-transitive automorphisms", viol)
    RES["controls"][f"BIAS/{nm}"] = dict(ds=True, n_link_amps=len(amps), critical=True,
                                         invariance_broken=True)

# ============================================================================ S5
section("S5  NON-MEASURE-CONSERVING CONTROLS (D-15): the leak and the polarized bath")

LEAK_S = Fr(9, 10)
for nm, V in VENUES.items():
    K = msca(LEAK_S, kernel_uniform(V, F0))
    rs = set(row_sums(K))
    gate(f"S5.1 {nm} CTRL-LEAK (declared survival 9/10): double stochasticity FAILS "
         "(computed row sums == 9/10 != 1): measure NOT conserved", not is_doubly_stochastic(K)
         and rs == {LEAK_S})
    dl = crit_witness(V, K)
    gate(f"S5.2 {nm} CTRL-LEAK: OFF criticality -- det(I - K) != 0 exactly "
         "(the gap is open; the conserving partner's det == 0 sits beside it in S3.3)", dl != 0)
    print(f"DATA  S5.2 {nm}: det(I - K_leak) = {dl}")
    amp = {K[x][y] for x in range(V["n"]) for y in V["adj"][x]}
    gate(f"S5.3 {nm} CTRL-LEAK: induced per-link amplitude == 9/(10*deg) != 1/deg "
         "(computes OFF the measure-conserving normalization)",
         amp == {LEAK_S / V["deg"]} and LEAK_S / V["deg"] != Fr(1, V["deg"]))
    v1 = [F1] * V["n"]
    mass_ok = True
    v = v1[:]
    for k in range(1, 9):
        v = mvec(K, v)
        mass_ok &= all(x == LEAK_S ** k for x in v)
    Kc = kernel_uniform(V, F0)
    vc = v1[:]
    for k in range(1, 9):
        vc = mvec(Kc, vc)
    gate(f"S5.4 {nm} CTRL-LEAK: surviving measure after k steps == (9/10)^k exactly "
         "(k = 1..8) while the conserving kernel holds measure 1 at the same depths "
         "(the mass gap IS the declared bias)", mass_ok and all(x == F1 for x in vc))
    RES["controls"][f"LEAK/{nm}"] = dict(row_sum=str(LEAK_S), det_I_minus_K=str(dl),
                                         per_link=str(LEAK_S / V["deg"]))

# ---- CTRL-BATH: the iv' dilation route, machine-checked at its smallest honest size.
# System qubit (excitation = defect) + bath qubit; U = exchange: |1,0> <-> |0,1>, else fix.
U = [[F0] * 4 for _ in range(4)]  # basis order |s b>: 00, 01, 10, 11
U[0][0] = F1
U[3][3] = F1
U[1][2] = F1  # |01> -> |10>
U[2][1] = F1  # |10> -> |01>
UtU = [[sum(U[k][i] * U[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
gate("S5.5 CTRL-BATH: the exchange dilation U is unitary (U^T U == I, permutation)", meq(UtU, eye(4)))
exc = [0, 1, 1, 2]  # total excitation number per basis state
img = [0, 2, 1, 3]  # U's permutation action
gate("S5.6 CTRL-BATH: U conserves total excitation number on every basis state "
     "(the energy-conserving dilation, iv')", all(exc[img[i]] == exc[i] for i in range(4)))

BATH_P = [Fr(1, 2), Fr(3, 4), F1]
bath_rows = {}
for p in BATH_P:
    # induced kernel on the system: T[s][s'] (s -> s'), bath traced with weights (p, 1-p)
    T = [[F0, F0], [F0, F0]]
    for s in (0, 1):
        for b, wb in ((0, p), (1, F1 - p)):
            out = img[2 * s + b]
            s_out = out >> 1
            T[s][s_out] += wb
    tp = all(x == F1 for x in row_sums(T))  # forward probability preserved
    ds = is_doubly_stochastic(T)
    bath_rows[p] = (tp, ds, col_sums(T))
gate("S5.7 CTRL-BATH: the induced system channel is trace-preserving (forward-stochastic) "
     "at EVERY bath bias", all(v[0] for v in bath_rows.values()))
gate("S5.8 CTRL-BATH: measure conservation holds at the UNBIASED bath p == 1/2 and FAILS "
     "at p == 3/4 and p == 1 (computed) -- a polarized bath is exactly a non-measure-"
     "conserving writer ensemble; the iv' route needs the bath condition, and the "
     "condition is the physics (bias == mass), not a gap in the route",
     bath_rows[Fr(1, 2)][1] and not bath_rows[Fr(3, 4)][1] and not bath_rows[F1][1])
for p in BATH_P:
    print(f"DATA  S5.8 bath p = {p}: incoming-measure (column) sums = "
          f"[{bath_rows[p][2][0]}, {bath_rows[p][2][1]}]")
RES["controls"]["BATH"] = {str(p): dict(trace_preserving=bath_rows[p][0], ds=bath_rows[p][1],
                                        col_sums=[str(x) for x in bath_rows[p][2]])
                           for p in BATH_P}

# ============================================================================ S6
section("S6  SYNTHESIS -- the computed verdict booleans")

def allpass(prefixes):
    return all(ok for n, ok in GATES if any(n.startswith(p) for p in prefixes))


gate("S6.1 EVERY unitary-writer ensemble computed in this lane (local, lazy, global, "
     "stabilizer-augmented, all declared non-uniform weight rows, every venue) induced a "
     "doubly stochastic kernel: measure conservation is FORCED by the writer algebra",
     allpass(["S0.1", "S0.2", "S3.1 ", "S3.5 ", "S3.8", "S3.9", "S3.10", "S3.11", "S3.12"]))
gate("S6.2 EVERY measure-conserving kernel computed here is CRITICAL (exact singular "
     "witness), uniform or not, lazy or not: masslessness == measure conservation at "
     "kernel level, with no tuning anywhere",
     allpass(["S3.3", "S3.7", "S4.2", "S4.6"]))
gate("S6.3 EVERY non-measure-conserving control failed double stochasticity AND the leak "
     "controls compute OFF criticality with the declared bias as the gap (DONE_WHEN control)",
     allpass(["S0.3", "S5.1", "S5.2", "S5.3", "S5.4", "S5.8"]))
gate("S6.4 uniformity: FORCED to c*I + t*A with t == (1-c)/deg under ensemble invariance "
     "on every venue (S4.1-S4.3); NOT forced by measure conservation alone (S4.4-S4.7 "
     "counterexample): the needed supplement is exactly invariance under the venue's "
     "edge-transitive automorphisms -- which the surface supplies as the writer-set "
     "ORBIT (S2.3) and as CARRIER symmetry (S2.4), while the ensemble WEIGHTS' symmetry "
     "remains extra data", allpass(["S4.1", "S4.2", "S4.3", "S4.4", "S4.5", "S4.6", "S4.7"]))
gate("S6.5 the induced per-crossing amplitude == 1/deg EXACTLY for every measure-"
     "conserving symmetric construction computed (E-LOC, every lazy row, E-GLOB "
     "conditioned at the moving end): within the declared one-parameter string model, "
     "the surface's own ensemble induces mu == mu_c identically",
     allpass(["S3.2", "S3.6", "S3.15", "S4.3"]))

# ---------------------------------------------------------------------------- summary
print()
print("=" * 88)
npass = sum(1 for _, ok in GATES if ok)
nfail = len(GATES) - npass
print(f"GATES: {npass} PASS, {nfail} FAIL, {len(GATES)} total")
print("=" * 88)

RES["gates"] = {"pass": npass, "fail": nfail, "total": len(GATES)}
RES["gate_list"] = [[n, ok] for n, ok in GATES]
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_T48_A_DERIVATION/"
          "t48a_derivation.RESULT.json", "w") as f:
    json.dump(RES, f, indent=1)
print("RESULT written: t48a_derivation.RESULT.json")
