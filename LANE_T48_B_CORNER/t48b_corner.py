"""T48-B CORNER: THE INDUCED PER-LINK AMPLITUDE OF THE SURFACE'S OWN WRITER ENSEMBLE (O-58 N2, probe 2).

QUESTION (O-58 N2, the labeled hypothesis of the T-44 judge, to be TESTED not assumed):
mu_c = 1/deg on every computed venue, and mu = 1/deg is exactly the stochastic
measure-conserving normalization of one writer step; if Gamma's writer kernel conserves
measure, mu = mu_c IDENTICALLY -- masslessness is measure conservation, not tuning.

THE TASK: compute, from the record surface's OWN writer ensemble -- not a declared mu --
the induced per-link amplitude on the corner venue (the D=2 plaquette adjacency of
LANE_T44_A_CORNER, rebuilt here from carrier supports exactly as its PUBLISHED_CONVENTIONS
prescribes), and determine whether it is the measure-conserving one.

GUARD AGAINST THE CIRCLE (binding): uniformity over links must be EARNED from the writer
algebra and venue symmetry (or shown to fail), never inserted.  If more than one honest
construction of "the writer ensemble" exists, ALL are computed and reported.

HONEST CONSTRUCTIONS SWEPT (each induces its own amplitude; none picked silently):
  K0  BARE COUNTING: the admissible-writer coset's own counting measure (each elementary
      writer once).  Induced per-link amplitude 1.  Conservation is COMPUTED, not assumed.
  K1  ALGEBRA-MEASURABLE: ensemble weight may depend only on the writer algebra's own
      computed invariants of an elementary writer (Gamma price, syndrome, coset data).
      Whether those invariants distinguish any two links is COMPUTED (Section 3); if they
      cannot, this class is link-uniform with ONE unknown scalar, and the measure-conserving
      member (if any) is computed exactly.
  K2  SYMMETRY-MEASURABLE (price-blind): ensemble weight may depend on the venue-graph
      automorphism orbit of the link.  The orbit structure is COMPUTED (Section 2); the
      conserving subfamily and its criticality are computed member by member.
D-15 CONTROLS (demanded by the assignment):
  (a) biased ensembles with DECLARED bias parameters, uniform and anisotropic, that must
      land computably OFF criticality (and off 1/deg);
  (b) the D=1 chain venue, where the same constructions must return its own 1/2;
  (+) a SQUARE (5,5) control venue (same instrument, in-lane control only) to test whether
      any symmetry-orbit freedom found on the rectangle venues is a two-scale artifact.

DISCIPLINE:
  D-24: venue rebuilt from carrier supports; the ensemble constructions are gated in weight
        and count against the O-54-C admissible-writer coset exactly as T-44 did; mu_c is
        re-located IN-LANE by computation (Perron row sums + exact sector sandwich +
        exact-rational resolvent), never imported as a number.
  D-8 : every verdict is a computed boolean; equality targets are COMPUTED references
        (deg from row sums, mu_c from the in-lane location), never literal expected values.
  D-15: every reported zero carries a positive control beside it that registers non-zero.
  D-1 : no gravitational form appears anywhere in this lane (not even in comparison).
  Exact ints and Fractions on every measurement path; NO floats anywhere in this lane.

RELEVANCE TEST (borrowed idea -> named program variable, named FIRST): the borrowed
Perron-Frobenius / Gershgorin / stochastic-kernel machinery is applied to the named
variable W(c -> c') -- the measure the surface's writer ensemble assigns to one elementary
admissible writer crossing the link c -> c' of the venue's own dual lattice -- and to no
other quantity.  Owners: Perron 1907 / Frobenius 1912; Gershgorin 1931; stochastic-matrix
row-sum facts and invariant-measure uniqueness: standard Markov-chain theory (Feller I);
uniqueness of the invariant (Haar) measure on a finite group: standard.  All applied only
through exact row-sum and resolvent computations on the venue.
"""
import json
import sys
from collections import deque
from fractions import Fraction
from math import comb

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T48_B_CORNER"
T44A = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T44_A_CORNER"
for p in (LANE, T44A):
    if p not in sys.path:
        sys.path.insert(0, p)

from t44a_lib import (Torus, sp_pair, rank_f2, independent_subset, generator_graph_dist,
                      coset_min_np, plaquette_adjacency, walk_counts_torus, bfs_dist,
                      cycle_adjacency, resolvent_exact)

GATES = []


def gate(name, ok, extra=""):
    GATES.append((name, bool(ok)))
    print(("PASS  " if ok else "FAIL  ") + name + (("  " + extra) if extra else ""))
    return bool(ok)


def ff(x, nd=6):
    """Display-only decimal rendering of an exact Fraction."""
    x = Fraction(x)
    sign = "-" if x < 0 else ""
    x = abs(x)
    scaled = (x.numerator * 10 ** nd) // x.denominator
    return "%s%d.%0*d" % (sign, scaled // 10 ** nd, nd, scaled % 10 ** nd)


# ---------------- DECLARED PARAMETERS AND THRESHOLDS (the whole declaration surface) ----
# The corner venues and pair lists are LANE_T44_A_CORNER's own (its declared HOLE_CASES).
SAMPLE_PAIRS = {
    (4, 6): [(1, 0), (2, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2), (2, 2), (1, 3), (2, 3)],
    (3, 7): [(1, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3)],
}
SQUARE_CONTROL = (5, 5)          # in-lane symmetry-control venue (square: both scales equal)
SQUARE_PAIRS = [(1, 0), (0, 1)]
CHAIN_L = 24                     # D=1 chain venue (T-44-A's own D-15 discriminator venue)
K_MASS = 40                      # depth of exact one-step-measure propagation gates
RET_K = 40                       # depth of the even-step return-term criticality witness
NEU_K = 60                       # Neumann-vs-resolvent bracket depth on subcritical controls
SECTOR_M = 300                   # direct range of the exact sector-sandwich binomial gates
IND_M = 10000                    # range of the exact induction-step inequality gates
BIAS_UNIF = [Fraction(9, 10), Fraction(11, 10)]     # control (a): uniform bias parameters beta
ALPHAS = [Fraction(1, 3), Fraction(3, 5)]           # K2 conserving-family sweep parameters
ANISO_BIAS_ALPHAS = [Fraction(1, 3), Fraction(-1, 3)]  # control (a): one-orbit-only bias
MU_BESIDE = [Fraction(1, 8), Fraction(23, 100)]     # nonsingular-beside-the-pole checks

print("=" * 100)
print("T48-B  CORNER: the induced per-link amplitude of the surface's OWN writer ensemble (O-58 N2)")
print("  No declared mu anywhere on the measurement path.  The ensemble constructions, their induced")
print("  amplitudes, conservation, and criticality are all computed exactly on the venue's own lattice.")
print("=" * 100)


# ---------------------------------------------------------------- generic exact helpers
def build_W(links, orbit_of, amp_of_orbit, n_cells):
    """Weighted one-step writer kernel: W[c][c'] = amplitude of the elementary writer
    crossing link {c,c'} (per direction; the writer is one involution serving both).
    Exact Fractions.  links: list of (i, j, mult); orbit_of: link index -> orbit id."""
    W = [[Fraction(0)] * n_cells for _ in range(n_cells)]
    for li, (i, j, m) in enumerate(links):
        a = Fraction(amp_of_orbit[orbit_of[li]]) * m
        W[i][j] += a
        W[j][i] += a
    return W

def row_sums(W):
    return [sum(row) for row in W]

def mat_vec_T(W, v):
    """v_{k+1}[c'] = sum_c v_k[c] W[c][c']  (measure propagation, exact)."""
    n = len(W)
    out = [Fraction(0)] * n
    for c in range(n):
        vc = v[c]
        if vc:
            Wc = W[c]
            for cp in range(n):
                if Wc[cp]:
                    out[cp] += vc * Wc[cp]
    return out

def resolvent_weighted(W, src):
    """EXACT rational solve of (I - W) x = e_src.  None if singular (zero pivot column)."""
    n = len(W)
    M = [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] -= W[i][j]
    rhs = [Fraction(0)] * n
    rhs[src] = Fraction(1)
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        rhs[col], rhs[piv] = rhs[piv], rhs[col]
        inv = 1 / M[col][col]
        M[col] = [x * inv for x in M[col]]
        rhs[col] *= inv
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [xr - f * xc for xr, xc in zip(M[r], M[col])]
                rhs[r] -= f * rhs[col]
    return rhs

def all_automorphisms(nbrsets):
    """ALL graph automorphisms of a simple connected graph (exhaustive backtracking with
    BFS-order pruning; adjacency preserved in both directions at every extension)."""
    n = len(nbrsets)
    order, seen = [0], {0}
    dq = deque([0])
    while dq:
        u = dq.popleft()
        for w in sorted(nbrsets[u]):
            if w not in seen:
                seen.add(w)
                order.append(w)
                dq.append(w)
    assert len(order) == n
    deg = [len(s) for s in nbrsets]
    autos = []
    img = [-1] * n
    used = [False] * n

    def bt(i):
        if i == n:
            autos.append(tuple(img))
            return
        v = order[i]
        cands = None
        for w in nbrsets[v]:
            iw = img[w]
            if iw >= 0:
                cands = set(nbrsets[iw]) if cands is None else (cands & nbrsets[iw])
        if cands is None:
            cands = set(range(n))
        for c in sorted(cands):
            if used[c] or deg[c] != deg[v]:
                continue
            ok = True
            for j in range(i):
                u2 = order[j]
                if (u2 in nbrsets[v]) != (img[u2] in nbrsets[c]):
                    ok = False
                    break
            if ok:
                img[v] = c
                used[c] = True
                bt(i + 1)
                img[v] = -1
                used[c] = False

    bt(0)
    return autos

class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb

def orbit_partition(items, autos, act):
    idx = {it: i for i, it in enumerate(items)}
    uf = UF(len(items))
    for a in autos:
        for it in items:
            uf.union(idx[it], idx[act(a, it)])
    roots = {}
    lab = []
    for i, it in enumerate(items):
        r = uf.find(i)
        if r not in roots:
            roots[r] = len(roots)
        lab.append(roots[r])
    return lab, len(roots)


# ================================================================ S0 venues from supports
print("\n-- SECTION 0: venues REBUILT FROM CARRIER SUPPORTS (T-44-A's own construction); the")
print("   elementary writer stratum gated in weight AND count against the O-54-C admissible-")
print("   writer coset, exactly as T-44 did --")

VEN = {}
S0_ROWS = []
for venue in [(4, 6), (3, 7), SQUARE_CONTROL]:
    Lx, Ly = venue
    tagc = " [IN-LANE SQUARE CONTROL VENUE]" if venue == SQUARE_CONTROL else ""
    T = Torus(Lx, Ly)
    n = T.n
    stars = T.all_stars()
    plaqs = T.all_plaqs()
    gate("S0 carrier (%d,%d)%s: stars x plaqs all commute" % (Lx, Ly, tagc),
         all(sp_pair(s, p, n) == 0 for s in stars for p in plaqs))
    cells, idx, adj = plaquette_adjacency(T)
    nc = len(cells)
    gate("S0 venue (%d,%d): all link multiplicities exactly 1 (no L=2 direction)" % venue,
         all(m == 1 for row in adj for _, m in row))
    degs = [sum(m for _, m in row) for row in adj]
    deg_v = degs[0]
    gate("S0 venue (%d,%d): row sums all equal => deg = %d COMPUTED from supports" % (Lx, Ly, deg_v),
         all(d == deg_v for d in degs))
    gate("S0 venue (%d,%d): no self-links (no cell shares a carrier edge with itself) -- the"
         " one-step ensemble has NO stay option; every elementary writer crosses exactly one link"
         % venue, all(all(j != i for j, _ in adj[i]) for i in range(nc)))
    # each carrier edge lies in exactly 2 plaquette supports -> elementary writers <-> links
    edge_in = [0] * n
    for (x, y) in cells:
        for e in T.plaq_edges(x, y):
            edge_in[e] += 1
    gate("S0 venue (%d,%d): every carrier edge lies in exactly 2 plaquette supports => the"
         " single-edge X writers are in bijection with the links of the dual lattice" % venue,
         all(c == 2 for c in edge_in))
    links = []
    edge_of_link = {}
    for i, (x, y) in enumerate(cells):
        si = set(T.plaq_edges(x, y))
        for j in range(i + 1, nc):
            common = si & set(T.plaq_edges(*cells[j]))
            if common:
                links.append((i, j, len(common)))
                edge_of_link[(i, j)] = sorted(common)
    gate("S0 venue (%d,%d): #links == #carrier edges == %d (the bijection counted both ways)"
         % (Lx, Ly, n), len(links) == n and sum(m for _, _, m in links) == n)
    src = idx[(0, 0)]
    counts = walk_counts_torus(adj, src, 12)
    dist = bfs_dist(adj, src)
    local = stars + plaqs
    plaq_index = {(x, y): y * Lx + x for y in range(Ly) for x in range(Lx)}
    gen_u = len(stars) + plaq_index[(0, 0)]
    stars_ind = independent_subset(stars)
    gens_X = stars_ind + [T.xbar1(), T.xbar2()]
    assert rank_f2(gens_X) == len(gens_X)
    xg = [g & ((1 << n) - 1) for g in gens_X]
    VEN[venue] = dict(T=T, cells=cells, idx=idx, adj=adj, links=links, edge_of_link=edge_of_link,
                      counts=counts, stars=stars, plaqs=plaqs, xg=xg, n=n, nc=nc, deg=deg_v,
                      src=src)
    pair_list = SAMPLE_PAIRS.get(venue, SQUARE_PAIRS)
    print("  venue (%d,%d)%s:  v  d_gen  d_dual  w_min  N_min  k0(walk)  N_k0" % (Lx, Ly, tagc))
    ok_dd = ok_w = ok_k0 = ok_nm = True
    for v in pair_list:
        d_gen = generator_graph_dist(local, n, gen_u, len(stars) + plaq_index[v])
        d_dual = dist[idx[v]]
        p = T.dual_path_x((0, 0), v)
        wmin, nmin, tot = coset_min_np(xg, p & ((1 << n) - 1))
        k0 = next(k for k in range(13) if counts[k][idx[v]] > 0)
        nk0 = counts[k0][idx[v]]
        ok_dd &= (d_dual == d_gen)
        ok_w &= (wmin == d_gen)
        ok_k0 &= (k0 == d_gen)
        ok_nm &= (nk0 == nmin)
        S0_ROWS.append(dict(venue=list(venue), v=list(v), d_gen=d_gen, d_dual=d_dual,
                            w_min=wmin, N_min=nmin, k0=k0, N_k0=nk0))
        print("   v=%-7s  %d      %d       %d      %d      %d        %d"
              % (str(v), d_gen, d_dual, wmin, nmin, k0, nk0))
    gate("S0 (%d,%d): dual BFS distance == earned separation d_gen on every sample pair" % venue, ok_dd)
    gate("S0 (%d,%d): coset w_min == d_gen (O-54-C landscape reconfirmed in-lane: w_min = d)" % venue, ok_w)
    gate("S0 (%d,%d): minimal walk length k0 == d on every sample pair" % venue, ok_k0)
    gate("S0 (%d,%d): minimal-walk count == coset N_min EXACTLY (gated in weight AND count,"
         " as T-44 did, before any ensemble is built)" % venue, ok_nm)
    # THE ELEMENTARY STRATUM: adjacent pairs, single-edge writer as coset representative
    ok_e1 = ok_ec = ok_syn = True
    for (i, j, m) in links:
        es = edge_of_link[(i, j)]
        ok_e1 &= (len(es) == m == 1)
        wmask = 1 << es[0]
        # syndrome of the elementary writer: exactly the two plaquettes of its link
        syn = [ci for ci, (x, y) in enumerate(cells) if sp_pair(wmask, T.plaq(x, y), n) == 1]
        ok_syn &= (syn == sorted([i, j]))
    for v in [vv for vv in pair_list if dist[idx[vv]] == 1]:
        e = edge_of_link[tuple(sorted([src, idx[v]]))][0]
        wmin, nmin, tot = coset_min_np(xg, 1 << e)
        ok_ec &= (wmin == 1 and nmin == 1)
    gate("S0 (%d,%d): every link carries exactly ONE elementary writer (multiplicity 1) whose"
         " syndrome is exactly its two cells (computed by symplectic pairing on every link)" % venue,
         ok_e1 and ok_syn)
    gate("S0 (%d,%d): coset of the single-edge writer on adjacent sample pairs has (w_min, N_min)"
         " == (1, 1): the weight-1 stratum of the admissible-writer coset IS the link, in weight"
         " AND count" % venue, ok_ec)
    gate("S0 (%d,%d): number of weight-1 admissible writers incident to each cell == deg == %d"
         " (the one-step continuation set at a cell is its link set)" % (Lx, Ly, deg_v),
         all(sum(1 for (i, j, m) in links if i == c or j == c) == deg_v for c in range(nc)))

# chain venue (D-15 control (b); T-44-A's own D=1 discriminator instrument)
CH = dict(adj=cycle_adjacency(CHAIN_L), nc=CHAIN_L)
CH["links"] = [(i, (i + 1) % CHAIN_L, 1) for i in range(CHAIN_L)]
CH["links"] = [(min(i, j), max(i, j), m) for (i, j, m) in CH["links"]]
CH["deg"] = 2
gate("S0 chain venue C_%d: row sums all exactly 2 => deg = 2 COMPUTED (D=1 control venue)" % CHAIN_L,
     all(sum(m for _, m in row) == 2 for row in CH["adj"]))

# ================================================================ S1 mu_c located in-lane
print("\n-- SECTION 1: mu_c RE-LOCATED IN-LANE (never imported as a number).  Method (T-44-A's,")
print("   re-run): radius of convergence of the positive walk series = 1/spectral radius;")
print("   Perron row-sum sandwich + exact sector sandwich + exact-rational resolvent --")
okL2b = all((2 * m + 1) * (2 * m + 3) <= (2 * m + 2) ** 2 for m in range(1, IND_M + 1))
okL2a = all(comb(2 * m, m) ** 2 * (2 * m + 1) <= 16 ** m for m in range(1, SECTOR_M + 1))
okL3low = all(comb(2 * m, m) * (2 * m + 1) >= 4 ** m for m in range(1, SECTOR_M + 1))
okmono = all((2 * m + 1) ** 2 * (m + 1) >= 4 * m * (m + 1) ** 2 for m in range(1, IND_M + 1))
gate("S1 sector sandwich lemmas: (C(2m,m))^2(2m+1) <= 16^m (m<=%d, exact induction step to %d)"
     " and C(2m,m)(2m+1) >= 4^m (m<=%d, exact ratio induction to %d)"
     % (SECTOR_M, IND_M, SECTOR_M, IND_M), okL2a and okL2b and okL3low and okmono)
gate("S1 => Z^2 sector: 16^m/(2m+1)^2 <= N_2m(0,0) = C(2m,m)^2 <= 16^m: the venue-limit"
     " return series has radius EXACTLY 1/4 (mu_c = 1/4 in the venue limit, computed)",
     okL2a and okL3low)
MUC = {}
for venue in [(4, 6), (3, 7), SQUARE_CONTROL]:
    V = VEN[venue]
    deg_v = V["deg"]
    mu_c = Fraction(1, deg_v)
    MUC[venue] = mu_c
    ok_perron = all(sum(m for _, m in row) == deg_v for row in V["adj"])
    gate("S1 (%d,%d): A.1 == %d.1 exactly (Perron eigenvector with positive entries) AND"
         " Gershgorin row bound %d => spectral radius = %d EXACT => mu_c(venue) = 1/%d"
         % (venue[0], venue[1], deg_v, deg_v, deg_v, deg_v), ok_perron)
    sing = resolvent_exact(V["adj"], mu_c, V["src"])
    gate("S1 (%d,%d): exact resolvent (I - A/%d) SINGULAR (zero pivot: the pole sits AT the"
         " computed mu_c)" % (venue[0], venue[1], deg_v), sing is None)
    ok_beside = all(resolvent_exact(V["adj"], mub, V["src"]) is not None for mub in MU_BESIDE)
    gate("S1 (%d,%d): exact resolvent NONSINGULAR beside it at mu in {1/8, 23/100} (positive"
         " control for the singularity zero)" % venue, ok_beside)
okch = all(sum(m for _, m in row) == 2 for row in CH["adj"])
MUC["chain"] = Fraction(1, 2)
gate("S1 chain: A.1 == 2.1 exactly => mu_c(chain) = 1/2 EXACT; sector side: C(2m,m)(2m+1) >="
     " 4^m re-verified above (the D=1 venue's own number, by the same instrument)", okch and okL3low)
sing_ch = resolvent_exact(CH["adj"], Fraction(1, 2), 0)
ok_ch_beside = resolvent_exact(CH["adj"], Fraction(9, 20), 0) is not None
gate("S1 chain: exact resolvent SINGULAR at 1/2, NONSINGULAR at 9/20", sing_ch is None and ok_ch_beside)

# ================================================================ S2 venue symmetry (computed)
print("\n-- SECTION 2: VENUE SYMMETRY, COMPUTED.  The full automorphism group of the venue's own")
print("   dual lattice (exhaustive backtracking) and its orbits on cells, links, directed links.")
print("   Uniformity by SYMMETRY is earned exactly where the orbit count on links is 1 --")
SYM = {}
for key in [(4, 6), (3, 7), SQUARE_CONTROL, "chain"]:
    if key == "chain":
        nc = CH["nc"]
        nbr = [set(j for j, _ in CH["adj"][i]) for i in range(nc)]
        links = CH["links"]
        name = "chain C_%d" % CHAIN_L
    else:
        V = VEN[key]
        nc = V["nc"]
        nbr = [set(j for j, _ in V["adj"][i]) for i in range(nc)]
        links = V["links"]
        name = "(%d,%d)" % key
    autos = all_automorphisms(nbr)
    litems = [(i, j) for (i, j, m) in links]
    llab, nlorb = orbit_partition(litems, autos, lambda a, it: (min(a[it[0]], a[it[1]]),
                                                               max(a[it[0]], a[it[1]])))
    ditems = [(i, j) for (i, j) in litems] + [(j, i) for (i, j) in litems]
    dlab, ndorb = orbit_partition(ditems, autos, lambda a, it: (a[it[0]], a[it[1]]))
    clab, ncorb = orbit_partition(list(range(nc)), autos, lambda a, it: a[it])
    orbit_sizes = [llab.count(r) for r in range(nlorb)]
    # per-cell orbit-degree vector
    odeg = []
    for c in range(nc):
        vec = [0] * nlorb
        for li, (i, j) in enumerate(litems):
            if i == c or j == c:
                vec[llab[li]] += 1
        odeg.append(tuple(vec))
    odeg_const = all(v == odeg[0] for v in odeg)
    SYM[key] = dict(autos=autos, litems=litems, llab=llab, nlorb=nlorb, ncorb=ncorb,
                    ndorb=ndorb, orbit_sizes=orbit_sizes, odeg=odeg[0])
    print("  venue %s: |Aut| = %d (computed exhaustively); orbits: cells %d, links %d %s,"
          % (name, len(autos), ncorb, nlorb, orbit_sizes)
          + " directed links %d; per-cell orbit-degree vector %s" % (ndorb, str(odeg[0])))
    gate("S2 %s: Aut acts TRANSITIVELY on cells (orbit count == 1: no cell is special;"
         " symmetry-measurable ensembles are cell-independent)" % name, ncorb == 1)
    gate("S2 %s: per-cell link-orbit degree vector CONSTANT across cells (the conservation"
         " condition will be one identical linear condition at every cell)" % name, odeg_const)
    gate("S2 %s: directed-link orbits == undirected-link orbits (each undirected orbit is one"
         " directed orbit: reflections reverse every link; no honest DIRECTED asymmetry exists"
         " in the symmetry-measurable class)" % name, ndorb == nlorb)
# the earned/failed uniformity finding, with D-15 zero/positive pairing
for venue in [(4, 6), (3, 7)]:
    S = SYM[venue]
    rep0 = next(it for li, it in enumerate(S["litems"]) if S["llab"][li] == 0)
    repo = {S["llab"][li]: it for li, it in enumerate(S["litems"])}
    cross = 0
    within = 0
    for a in S["autos"]:
        im = (min(a[rep0[0]], a[rep0[1]]), max(a[rep0[0]], a[rep0[1]]))
        li_im = S["litems"].index(im)
        if S["llab"][li_im] != 0:
            cross += 1
        elif im != rep0:
            within += 1
    gate("S2 (%d,%d): link-orbit count %d >= 2 with ZERO automorphisms mapping the orbit-0"
         " representative across orbits (computed %d) BESIDE the positive control: %d"
         " automorphisms move it non-trivially WITHIN its orbit and the orbit is full size %d"
         " => ON THIS TWO-SCALE VENUE, SYMMETRY ALONE DOES NOT EARN UNIFORMITY ACROSS ALL LINKS"
         " (it earns it only orbit-wise) -- an honest FAILURE, reported, not papered over"
         % (venue[0], venue[1], S["nlorb"], cross, within, S["orbit_sizes"][0]),
         S["nlorb"] == 2 and cross == 0 and within > 0 and
         S["orbit_sizes"][0] == len(S["litems"]) // 2)
Ssq = SYM[SQUARE_CONTROL]
gate("S2 (%d,%d) SQUARE CONTROL: link-orbit count == 1 (the axis swap is an automorphism"
     " here) => symmetry ALONE earns full link-uniformity when the venue's two scales are"
     " equal: the rectangle anisotropy freedom is a two-scale artifact of the finite venue,"
     " COMPUTED, not asserted" % SQUARE_CONTROL, Ssq["nlorb"] == 1)
gate("S2 chain: link-orbit count == 1 => symmetry alone earns uniformity on the D=1 venue",
     SYM["chain"]["nlorb"] == 1)

# ================================================================ S3 algebra invariants
print("\n-- SECTION 3: THE WRITER ALGEBRA'S OWN INVARIANTS of an elementary writer, computed on")
print("   EVERY link: Gamma price; syndrome size; star/plaquette memberships of its edge; and")
print("   (per orbit representative) the full coset data of its pair.  If these cannot separate")
print("   any two links, every ALGEBRA-MEASURABLE ensemble is link-uniform: EARNED, not inserted --")
ALG = {}
for venue in [(4, 6), (3, 7), SQUARE_CONTROL]:
    V = VEN[venue]
    T = V["T"]
    n = V["n"]
    cells = V["cells"]
    invs = []
    for (i, j, m) in V["links"]:
        e = V["edge_of_link"][(i, j)][0]
        wmask = 1 << e
        price = 1 if bin(wmask).count("1") == 1 else None
        syn = sum(1 for (x, y) in cells if sp_pair(wmask, T.plaq(x, y), n) == 1)
        in_stars = sum(1 for s in V["stars"] if (s >> e) & 1)
        in_plaqs = sum(1 for (x, y) in cells if e in T.plaq_edges(x, y))
        invs.append((price, syn, in_stars, in_plaqs, m))
    same = all(v == invs[0] for v in invs)
    ALG[venue] = dict(inv=invs[0], all_equal=same)
    print("  venue (%d,%d): invariant tuple (price, |syndrome|, #stars, #plaqs, mult) = %s on"
          % (venue[0], venue[1], str(invs[0])) + " all %d links" % len(invs))
    gate("S3 (%d,%d): the invariant tuple is IDENTICAL on every link (in particular identical"
         " across the two symmetry orbits): the writer algebra cannot distinguish any two links"
         % venue, same and invs[0] == (1, 2, 2, 2, 1))
    # coset data at one representative per symmetry orbit
    S = SYM[venue]
    repdata = []
    for r in range(S["nlorb"]):
        it = next(itm for li, itm in enumerate(S["litems"]) if S["llab"][li] == r)
        e = V["edge_of_link"][it][0]
        wmin, nmin, tot = coset_min_np(V["xg"], 1 << e)
        repdata.append((wmin, nmin, tot))
    gate("S3 (%d,%d): coset data (w_min, N_min, |coset|) of the elementary writer's pair is"
         " IDENTICAL across the %d symmetry orbit(s): %s (exhaustive scans; the coset cannot"
         " distinguish the orbits either)" % (venue[0], venue[1], S["nlorb"], str(repdata[0])),
         all(rd == repdata[0] for rd in repdata) and repdata[0][0] == 1 and repdata[0][1] == 1)
gate("S3 CONCLUSION (computed above): every ensemble whose weights depend only on the writer"
     " algebra's computed invariants is LINK-UNIFORM on every venue -- uniformity EARNED from"
     " the algebra; the ONLY computed structure that separates links is the graph-global"
     " symmetry-orbit label of Section 2 (price-blind), swept honestly as construction K2",
     all(ALG[v]["all_equal"] for v in ALG))

# ================================================================ S4 constructions swept
print("\n-- SECTION 4: THE HONEST CONSTRUCTIONS AND THEIR INDUCED PER-LINK AMPLITUDES.  For each:")
print("   conservation is a COMPUTED row-sum boolean; criticality is a COMPUTED Perron statement")
print("   (constant row sums = spectral radius sandwich) plus exact resolvent and return-term")
print("   witnesses.  Every zero paired with a non-zero positive control (D-15) --")

CONSTR_ROWS = []

def sweep_venue(key, deep):
    if key == "chain":
        nc, links, deg_v, name = CH["nc"], CH["links"], CH["deg"], "chain C_%d" % CHAIN_L
        src = 0
        mu_c = MUC["chain"]
    else:
        V = VEN[key]
        nc, links, deg_v, name = V["nc"], V["links"], V["deg"], "(%d,%d)" % key
        src = V["src"]
        mu_c = MUC[key]
    S = SYM[key]
    llab, R = S["llab"], S["nlorb"]
    odeg = S["odeg"]
    ret_lb = Fraction(1, nc)
    rows_here = []

    def run_member(tag, amps, expect_cons):
        """amps: per-orbit amplitude vector.  Everything computed; expect_cons only labels
        the printout row -- the booleans are the verdicts."""
        W = build_W(links, llab, amps, nc)
        rs = row_sums(W)
        cons = all(r == 1 for r in rs)
        sym_ok = all(W[i][j] == W[j][i] for i in range(nc) for j in range(nc))
        beta = rs[0]
        rs_const = all(r == beta for r in rs)
        # exact measure propagation from e_src
        v = [Fraction(0)] * nc
        v[src] = Fraction(1)
        mass_ok = True
        ret_ok = True
        ret_seen = 0
        for k in range(1, K_MASS + 1):
            v = mat_vec_T(W, v)
            mass_ok &= (sum(v) == beta ** k)
            if cons and k % 2 == 0:
                ret_seen += 1
                ret_ok &= (v[src] >= ret_lb)
        res = resolvent_weighted(W, src)
        singular = res is None
        row = dict(venue=name, construction=tag,
                   amplitudes=[str(Fraction(a)) for a in amps],
                   row_sum=str(beta), conserving=cons, row_sums_constant=rs_const,
                   symmetric=sym_ok, mass_law="beta^k exact, k<=%d" % K_MASS,
                   resolvent_singular=singular)
        rows_here.append(row)
        return W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res

    # ---- K0: bare counting measure (the coset's own combinatorics, no normalization)
    W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res = \
        run_member("K0 bare counting", [Fraction(1)] * R, False)
    gate("S4 %s K0 BARE COUNTING (amplitude 1 per elementary writer): row sums == %d != 1:"
         " NOT measure-conserving (COMPUTED).  The coset's own counting measure multiplies"
         " measure by deg each step (mass after k steps == %d^k exactly, k<=%d): the venue"
         " combinatorics alone does NOT hand the surface criticality -- conservation is a"
         " substantive property, not automatic" % (name, deg_v, deg_v, K_MASS),
         (not cons) and rs_const and beta == deg_v and mass_ok)

    # ---- K1: algebra-measurable => link-uniform (earned, S3); the conserving member
    t_star = Fraction(1, deg_v)      # computed: unique solution of deg * t == 1, deg from S0
    W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res = \
        run_member("K1 conserving (algebra-uniform)", [t_star] * R, True)
    gate("S4 %s K1: the algebra-measurable class is link-uniform (S3) with ONE unknown scalar t;"
         " measure conservation deg*t == 1 has the UNIQUE solution t* = 1/%d, and t* == mu_c"
         " == %s (the in-lane located critical point): THE INDUCED PER-LINK AMPLITUDE OF THE"
         " CONSERVING ENSEMBLE IS EXACTLY THE CRITICAL ONE" % (name, deg_v, str(mu_c)),
         t_star * deg_v == 1 and t_star == mu_c)
    gate("S4 %s K1 conservation COMPUTED: all row sums of W == 1 exactly AND the propagated"
         " one-step measure from e_src totals 1 exactly at every k <= %d (measure conserved,"
         " never created or lost)" % (name, K_MASS), cons and mass_ok)
    gate("S4 %s K1 criticality COMPUTED: constant row sums 1 (Perron sandwich => spectral"
         " radius == 1: the mass gap -ln(rowsum) is EXACTLY ZERO) AND exact resolvent (I - W)"
         " SINGULAR AND even-step return terms >= 1/%d at all %d computed depths (Cauchy-Schwarz"
         " floor for a symmetric conserving kernel: the coupling series DIVERGES: critical)"
         % (name, nc, RET_K // 2), cons and rs_const and singular and ret_ok and sym_ok)

    # ---- K2: symmetry-measurable (price-blind): the conserving subfamily, member by member
    fam = []
    if R >= 2:
        for alpha in ALPHAS:
            for r0 in range(R):
                amps = []
                for r in range(R):
                    if r == r0:
                        amps.append((1 + alpha) / deg_v)
                    else:
                        amps.append((1 - alpha * odeg[r0] / (odeg[r] * (R - 1))) / deg_v)
                fam.append((alpha, r0, amps))
    fam_ok = True
    fam_pins = True
    for (alpha, r0, amps) in fam:
        W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res = \
            run_member("K2 conserving (alpha=%s on orbit %d)" % (alpha, r0), amps, True)
        fam_ok &= cons and rs_const and singular and ret_ok and mass_ok and sym_ok
        fam_pins &= all(a == t_star for a in amps)
        print("      K2 member alpha=%s orbit %d: amplitudes %s  row sum %s  singular %s"
              % (alpha, r0, [str(a) for a in amps], beta, singular))
    if R >= 2:
        gate("S4 %s K2 SYMMETRY-MEASURABLE FAMILY: every swept conserving member (%d members,"
             " declared alpha recipe) has row sums == 1, resolvent singular, return terms"
             " floored, mass exactly conserved: EVERY member sits EXACTLY at criticality --"
             " masslessness <=> measure conservation holds in the anisotropic honest"
             " construction too" % (name, len(fam)), fam_ok and len(fam) == 2 * len(ALPHAS))
        gate("S4 %s K2 HONEST LIMIT: the conserving family does NOT pin the amplitude to 1/%d"
             " (swept members with per-orbit amplitudes != t* exist: computed).  Symmetry +"
             " conservation leave an (R-1)-parameter freedom on this two-scale venue; the"
             " VALUE 1/deg needs the algebra ground (S3) or link-transitivity (square/chain)"
             % (name, deg_v), (not fam_pins) and fam_ok)
        # the price-measurable member of the family is the uniform one
        gate("S4 %s K2 intersection with the algebra class: the unique conserving member with"
             " equal orbit amplitudes is t* == 1/%d == mu_c (computed: solve %s"
             " with equal amplitudes)" % (name, deg_v,
             " + ".join("%d t" % d for d in odeg) + " == 1"),
             sum(odeg) * t_star == 1 and t_star == mu_c)
    else:
        gate("S4 %s K2: link-orbit count == 1 (S2) => the symmetry-measurable conserving family"
             " COLLAPSES to the single member t = 1/%d == mu_c: on a link-transitive venue,"
             " SYMMETRY ALONE + conservation pins the critical amplitude" % (name, deg_v),
             sum(odeg) * t_star == 1 and t_star == mu_c)

    # ---- CONTROL (a): biased ensembles, declared parameters, must land OFF computably
    for bu in BIAS_UNIF:
        t_b = bu * t_star
        W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res = \
            run_member("CONTROL biased uniform beta=%s" % bu, [t_b] * R, False)
        off_amp = (t_b != mu_c)
        if bu < 1:
            # subcritical: nonsingular resolvent + Neumann bracket + geometric return bound
            v = [Fraction(0)] * nc
            v[src] = Fraction(1)
            S_part = [Fraction(0)] * nc
            S_part[src] += 1
            ret_bound_ok = True
            for k in range(1, NEU_K + 1):
                v = mat_vec_T(W, v)
                for c in range(nc):
                    S_part[c] += v[c]
                ret_bound_ok &= (v[src] <= bu ** k)
            tail = bu ** (NEU_K + 1) / (1 - bu)
            neu_ok = (res is not None) and all(S_part[c] <= res[c] <= S_part[c] + tail
                                              for c in range(nc))
            gate("S4 %s CONTROL(a) beta=%s: induced amplitude %s != mu_c (COMPUTED off 1/%d);"
                 " row sums %s != 1 (conservation FAILS); one-step mass == %s^k exactly (the"
                 " mass gap -ln(%s) is OPEN); resolvent NONSINGULAR with the exact Neumann"
                 " bracket containing it at depth %d; return terms <= %s^k: the coupling series"
                 " CONVERGES: computably OFF criticality, subcritical side"
                 % (name, bu, t_b, deg_v, beta, bu, bu, NEU_K, bu),
                 off_amp and (not cons) and rs_const and beta == bu and mass_ok
                 and neu_ok and ret_bound_ok)
        else:
            gate("S4 %s CONTROL(a) beta=%s: induced amplitude %s != mu_c (COMPUTED off 1/%d);"
                 " row sums %s != 1 (conservation FAILS the other way); one-step mass == %s^k"
                 " exactly, GROWING (measure created): computably OFF criticality,"
                 " mass-creation side" % (name, bu, t_b, deg_v, beta, bu),
                 off_amp and (not cons) and rs_const and beta == bu and mass_ok)
    if R >= 2 and deep:
        for alpha in ANISO_BIAS_ALPHAS:
            amps = [(1 + alpha) / deg_v] + [t_star] * (R - 1)
            W, beta, cons, rs_const, sym_ok, mass_ok, ret_ok, singular, res = \
                run_member("CONTROL biased aniso (alpha=%s on orbit 0 only)" % alpha, amps, False)
            expected_beta = 1 + odeg[0] * alpha / deg_v
            sub = beta < 1
            extra = ""
            okres = True
            if sub:
                okres = (res is not None)
                extra = "; resolvent NONSINGULAR (subcritical side)"
            gate("S4 %s CONTROL(a) anisotropic bias alpha=%s on orbit 0 only: row sums == %s"
                 " != 1 (conservation fails; == 1 + %d*alpha/%d computed)%s: an unbalanced"
                 " orbit ensemble lands computably OFF criticality"
                 % (name, alpha, beta, odeg[0], deg_v, extra),
                 (not cons) and rs_const and beta == expected_beta and mass_ok and okres)
    return rows_here

for key in [(4, 6), (3, 7)]:
    print("\n  ---- venue (%d,%d) ----" % key)
    CONSTR_ROWS += sweep_venue(key, deep=True)
print("\n  ---- SQUARE CONTROL venue (%d,%d) (symmetry control: family must collapse) ----" % SQUARE_CONTROL)
CONSTR_ROWS += sweep_venue(SQUARE_CONTROL, deep=False)
print("\n  ---- CONTROL (b): chain venue C_%d (the same constructions must return 1/2) ----" % CHAIN_L)
CONSTR_ROWS += sweep_venue("chain", deep=False)

# ================================================================ S5 verdict
print("\n" + "=" * 100)
print("SECTION 5 -- VERDICT (assembled from the computed booleans above; nothing new measured)")
print("=" * 100)
npass = sum(1 for _, ok in GATES if ok)
nfail = sum(1 for _, ok in GATES if not ok)
V1 = all(ok for nm, ok in GATES if "K1" in nm)
V2 = all(ok for nm, ok in GATES if "K2" in nm or nm.startswith("S3") or nm.startswith("S2"))
V3 = all(ok for nm, ok in GATES if "CONTROL" in nm)
V0 = all(ok for nm, ok in GATES if "K0" in nm)
verdict = [
    "1. THE INDUCED AMPLITUDE: on the corner venue the surface's own writer ensemble, in every",
    "   construction measurable with respect to the writer algebra's computed invariants, is",
    "   link-uniform (EARNED, Section 3) with one unknown scalar; the measure-conserving member",
    "   is UNIQUE and its induced per-link amplitude is EXACTLY 1/deg = 1/4 = mu_c (in-lane",
    "   located).  The same construction returns 1/2 on the D=1 chain -- each venue its own",
    "   number, nothing imported.  [V1: %s]" % V1,
    "2. THE IDENTITY IS STRUCTURAL, NOT NUMERICAL: conservation (row sums 1) and criticality",
    "   (Perron root at 1, pole of the resolvent, zero mass gap) are THE SAME computed row-sum",
    "   fact read twice.  In the price-blind symmetry-measurable class on the two-scale venues,",
    "   conservation leaves a computed 1-parameter anisotropy family -- and EVERY member of it",
    "   is still EXACTLY critical: masslessness <=> measure conservation survives every honest",
    "   construction; the specific VALUE 1/deg is pinned by the algebra's price-uniformity or",
    "   by link-transitivity (square control, chain).  [V2: %s]" % V2,
    "3. THE CONTROLS LAND OFF: every declared-bias ensemble (uniform 9/10, 11/10; anisotropic",
    "   one-orbit) fails conservation computably and sits computably off criticality (open mass",
    "   gap and convergent series below; growing mass above).  [V3: %s]" % V3,
    "4. THE HONEST NEGATIVE: the coset's own bare counting measure is NOT conserving (row sums",
    "   deg): the surface's combinatorics alone does not put the coupling at criticality.",
    "   WHETHER Gamma's actual writing conserves measure remains the open physical bit -- what",
    "   this lane computed is that the conserving normalization exists, is unique, carries no",
    "   tuning freedom, and IS the critical point identically.  [V0: %s]" % V0,
]
for line in verdict:
    print(line)
print("\nO-58 N2 AT CORNER-VENUE LEVEL: the hypothesis's conditional is now a computed identity --")
print("IF the writer kernel conserves measure THEN mu = mu_c IDENTICALLY (both directions, with")
print("uniqueness, uniformity earned, and off-criticality controls).  Masslessness = measure")
print("conservation = one bit, not a tuned value.  NOT decided here (named, kept open): whether")
print("the record surface's GR3 writing dynamics actually conserves writer measure.")
print("\nNEXT STEP (no route closes without one): formulate the conservation bit as a computable")
print("property of the GR3 writing tier (does one writing event redistribute writer measure or")
print("create/destroy it?), and run this probe's twin on the D=3 world venue (deg 6 -> 1/6).")
print("\nGATES: %d PASS, %d FAIL, total %d" % (npass, nfail, len(GATES)))

result = dict(
    lane="LANE_T48_B_CORNER",
    task="O-58 N2 probe 2 (corner venue): the induced per-link amplitude of the surface's own"
         " writer ensemble -- is it the measure-conserving one, and is that mu_c?",
    date="2026-08-20",
    declared=dict(
        venues="T-44-A corner venues (4,6) and (3,7) rebuilt from carrier supports; in-lane"
               " square control venue (5,5); chain C_24 (T-44-A's D=1 control)",
        constructions="K0 bare counting; K1 algebra-measurable; K2 symmetry-measurable"
                      " (price-blind); declared-bias controls",
        thresholds=dict(K_MASS=K_MASS, RET_K=RET_K, NEU_K=NEU_K, SECTOR_M=SECTOR_M,
                        IND_M=IND_M, BIAS_UNIF=[str(b) for b in BIAS_UNIF],
                        ALPHAS=[str(a) for a in ALPHAS],
                        ANISO_BIAS_ALPHAS=[str(a) for a in ANISO_BIAS_ALPHAS],
                        MU_BESIDE=[str(m) for m in MU_BESIDE]),
        no_floats="no floats anywhere in this lane; exact ints and Fractions only",
    ),
    coset_gates=S0_ROWS,
    mu_c_in_lane={"(4,6)": "1/4", "(3,7)": "1/4", "(5,5)": "1/4", "chain": "1/2",
                  "method": "Perron row-sum sandwich (A.1 = deg.1 exact + Gershgorin) + exact"
                            " sector sandwich 16^m/(2m+1)^2 <= C(2m,m)^2 <= 16^m + exact"
                            " resolvent singular at 1/deg, nonsingular beside"},
    symmetry=dict((str(k), dict(aut=len(SYM[k]["autos"]), link_orbits=SYM[k]["nlorb"],
                                orbit_sizes=SYM[k]["orbit_sizes"],
                                cell_orbits=SYM[k]["ncorb"],
                                per_cell_orbit_degree=list(SYM[k]["odeg"])))
                  for k in SYM),
    algebra_invariants=dict((("(%d,%d)" % k), dict(tuple=list(ALG[k]["inv"]),
                                                   identical_on_all_links=ALG[k]["all_equal"]))
                            for k in ALG),
    construction_rows=CONSTR_ROWS,
    finding=dict(
        induced_amplitude="the measure-conserving member of every algebra-measurable ensemble"
                          " class induces per-link amplitude EXACTLY 1/deg (1/4 corner, 1/2"
                          " chain), and 1/deg == mu_c in-lane: the induced amplitude IS the"
                          " measure-conserving one and IS the critical one, identically",
        uniformity="EARNED from the writer algebra (invariant tuple identical on every link,"
                   " computed); by symmetry alone only on link-transitive venues (chain, square"
                   " control); on the two-scale rectangles symmetry leaves a computed"
                   " 1-parameter conserving anisotropy family, every member exactly critical",
        masslessness="masslessness <=> measure conservation is the same computed row-sum fact"
                     " read twice (Perron root at 1 <=> row sums 1); it survives every honest"
                     " construction including the anisotropic family",
        honest_negative="the bare counting ensemble is NOT conserving (row sums deg): the venue"
                        " combinatorics does not force conservation; whether Gamma's writing"
                        " conserves measure is the remaining physical bit (GR3 tier)",
        controls="biased uniform 9/10 and 11/10 and one-orbit anisotropic biases all land"
                 " computably off conservation and off criticality; chain venue returns its own"
                 " 1/2 end to end",
    ),
    owners=dict(perron_frobenius="Perron 1907 / Frobenius 1912; Gershgorin 1931 (applied only"
                                 " through exact row-sum computations on the venue)",
                stochastic_kernels="standard Markov-chain theory (Feller I): row-stochastic ="
                                   " measure-conserving; invariant-measure facts",
                haar_uniform="uniqueness of the invariant measure on a finite group: standard",
                cauchy_schwarz="return-term floor W^2m(s,s) >= 1/n for symmetric conserving"
                               " kernels: Cauchy-Schwarz, verified computationally at every"
                               " reported depth",
                coset_instrument="LANE_O54_C_ATTEMPT / LANE_T44_A_CORNER machinery imported"
                                 " read-only at runtime (C-78 lineage)"),
    relevance_test="borrowed Perron/stochastic machinery applied to the named variable"
                   " W(c->c'), the measure the surface's writer ensemble assigns to one"
                   " elementary admissible writer crossing link c->c', and to nothing else",
    next_step="formulate the conservation bit as a computable property of the GR3 writing"
              " tier (does one writing event redistribute writer measure or create/destroy"
              " it?); run the twin probe on the D=3 world venue (deg 6 -> 1/6)",
    gates=dict(npass=npass, nfail=nfail, failed=[nm for nm, ok in GATES if not ok]),
)
with open(LANE + "/t48b_corner.RESULT.json", "w") as f:
    json.dump(result, f, indent=1)
print("\nRESULT JSON written: t48b_corner.RESULT.json")
