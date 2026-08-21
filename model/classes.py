"""THE CLASSES LAYER -- the reachable-class machinery (C-87, C-90) folded into the model
(T-54, family: classes).

THE CLAIM THIS LAYER CARRIES (C-87, updated at the C-90 landing): a Gamma-priced coupling
G_mu(d) = sum over admissible strings of mu^weight, strings priced by Gamma at w_min = d,
reaches EXACTLY THREE CLASSES, each located by computation:
  (1) EXPONENTIAL decay for mu < mu_c, leading term N_min mu^d (the confinement cost with
      its exact degeneracy d!/(a!b!c!));
  (2) at mu_c EXACTLY -- the venue's own computed number, located three independent ways
      (Perron row-sum, exact-rational resolvent pole, exact sector sandwich) -- the earned
      dimension's own critical Green's class: LINEAR in D=1, MARGINAL/LOG in D=2, POWER LAW
      exponent 1 (the 1/d member, C-90) in D=3; the critical identity IS the venue's
      discrete Poisson equation;
  (3) DIVERGENT term-by-term above mu_c: no mediated coupling in the venue limit.

Everything here is PORTED from the sealed lanes, with fidelity to sealed behavior over
elegance.  Exact arithmetic throughout the measurement path: python ints and Fractions
only; square roots enter only as certified rational brackets (math.isqrt) and only in
owner-attributed COMPARISON quantities, never in measurements.  Class labels are emitted
by COMPUTED BOOLEANS over exact rationals against the declared threshold constants below
(the sealed lanes' own declaration surface).  Returns are DATA.

WHAT EACH SECTION IMPLEMENTS, WITH ITS CLAIM ROW AND SEALED SOURCE:
  1  VENUE ENTRY (observation-entry; D-25)    C-87   (LANE_T44_B_WORLD/t44b_lib.py
     a venue graph is an adjacency structure          torus3_adjacency, walk_counts_adj,
     with DECLARED PROVENANCE; the world venue        bfs_adj; LANE_T44_A_CORNER/t44a_lib.py
     is the census access geometry, the corner        plaquette_adjacency -- rebuilt here on
     venue is computed from the model's own           the model's own ported carrier
     ported carrier (geometry.Torus), the             geometry.Torus -- and cycle_adjacency)
     chain is the D=1 control
  2  mu_c LOCATED (the resolvent route)       C-87/  (LANE_T44_B_WORLD S0/S2: Perron
     exact-rational Gaussian elimination:     C-90    row-sum; resolvent singular exactly
     (I - mu A) singular exactly AT the               at 1/6 on the D=3 venue, solvable at
     venue's own 1/deg, solvable beside;              19/120 and 7/40; LANE_T44_A_CORNER
     never a hard-coded critical point                S2: singular at 1/4 on the D=2 venue,
                                                      solvable at 1/8 and 23/100)
  3  THE PRICED SERIES G_mu(d)                C-87   (LANE_T44_B_WORLD/t44b_lib.py W1,
     exact partial sums, PROVEN geometric             N3_ref, n3_even_row, series_3d,
     tails; the subcritical class booleans            dp3_counts; S3 class booleans;
     and the divergence witnesses                     S2 witnesses)
  4  THE CRITICAL KERNEL (the 1/d member)     C-90   (LANE_T44_B_WORLD/t44b_lib.py
     a_M(x) = sum (N_2m(0)-N_2m(x))/36^m,             crit_kernel_3d + the certified tail
     certified tails B5 (m-2)^{-3/2} and the          machinery wallis_brackets, q3_constant,
     assembled L3-3D difference bound; the            diff_tail_bound, abs_tail_bound;
     doubling-increment ratio windows                 verifier deepening loop
                                                      VERIFY/adv_verify.py E-section,
                                                      register row C-90)
  5  CROSS-DIMENSION DISCRIMINATOR            C-87   (LANE_T44_A_CORNER/t44a_lib.py
     one instrument, three venues, three              potential_kernel_2d, potential_kernel_1d,
     DISJOINT declared windows: D=1 LIN,              partial_return_sum_2d, series_target_2d,
     D=2 LOG, D=3 INV                                 series_target_1d; LANE_T44_B_WORLD S5)
  6  THE CLASS VERDICT                        C-87   (LANE_T44_B_WORLD taxonomy; the triple
     per (venue, mu): a computed boolean              is computed against the VENUE'S OWN
     triple (exponential, critical,                   mu_c from section 2 -- exact rational
     divergent), with optional evidence               comparison, no literal critical point)
     tiers running the section-3/4/5
     instruments

BORROWED INSTRUMENTS, OWNERS NAMED WHERE THE LANES NAMED THEM: walk generating functions /
lattice Green's functions and radius of convergence of positive series are Spitzer 1964 /
Lawler territory; transience of D=3 is Polya 1921; G(0) at criticality on the simple-cubic
venue is Watson 1939 (1.5163860591...); the 1/d asymptote coefficient 3/(2 pi) is Spitzer
1964 P26.1; the axis-split bijection is standard combinatorics (Feller I); Wallis-product
brackets are classical; Perron-Frobenius/Gershgorin is standard linear algebra; the
Ornstein-Zernike subcritical rate is standard asymptotics; the confinement cost w_min = d
is the standing C-80/O-54 result (Wegner/Wilson).  Owner values appear ONLY in labeled
comparison quantities after each class is computed.  OURS: the Gamma-priced generating
function, criticality-is-the-venue's-Poisson, the located mu_c, the class windows, the
D=1/D=2/D=3 discriminator, per the register rows above.

OBSERVATION ENTRY (the URM directive, 2026-08-21): a NEW venue graph -- a new record
surface's access geometry -- enters as an adjacency structure through venue() with
DECLARED PROVENANCE (D-25: world-tier requires the pinned real-geometry source; an exact
idealisation must self-declare "DEF-A").  Then mu_c_of() locates ITS OWN critical point by
the resolvent route, venue_series() prices couplings on it, and class_verdict() places any
declared mu in the three-class taxonomy against that venue's own number.  Nothing in this
module hard-codes a critical point, a dimension, or a class label.

The sealed lanes remain the source of truth; checks_classes.py gates every number this
module reproduces against its SEALED value."""
import os as _os
import sys as _sys
from collections import deque as _deque
from fractions import Fraction
from math import comb, isqrt

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)


# =====================================================================================
# DECLARED THRESHOLD CONSTANTS (the sealed lanes' own declaration surface -- these are the
# instrument's DECLARED parameters, published in LANE_T44_B_WORLD/t44b_world.RESULT.json
# "declared.thresholds" and LANE_T44_A_CORNER equivalents; never fitted, never tuned here)
# =====================================================================================
MARGIN_LT1 = Fraction(1, 20)    # exponential row: every ratio upper bound <= 1 - margin
POWER_SCOPE = 8                 # power exclusion: d*(1-r) must exceed this
CAUCHY_TOL = Fraction(1, 40)    # ratio stabilization tolerance at the top of the d range
INV_LO, INV_HI = Fraction(2, 5), Fraction(3, 5)    # doubling-increment window: 1/d class
LOG_LO, LOG_HI = Fraction(4, 5), Fraction(5, 4)    # doubling-increment window: log class
LIN_LO, LIN_HI = Fraction(9, 5), Fraction(11, 5)   # doubling-increment window: linear class
ISO_TOL = Fraction(1, 20)       # coefficient agreement across rays
COMP_TOL = Fraction(1, 25)      # owner-comparison tolerance (comparisons only)
OOS_TOL = Fraction(1, 12)       # out-of-sample increment prediction tolerance
M0_WAL = 240                    # Wallis bracket anchor
M0_TRI = 300                    # max-trinomial telescoping anchor (multiple of 3)
EDGE_C = Fraction(40)           # L3-3D edge-region constant (lane-gated on range)
RHO = Fraction(199, 200)        # L3-3D edge decay rational (lane-gated)
WITNESS_M = 150                 # supercritical divergence witness depth
# owner constants (COMPARISON ONLY; brackets, never on the measurement path):
WATSON = (Fraction(1516386059, 10 ** 9), Fraction(1516386060, 10 ** 9))  # Watson 1939 G(0)
C_3D = (Fraction(477464829, 10 ** 9), Fraction(477464830, 10 ** 9))      # 3/(2 pi), Spitzer P26.1


def ff(x, nd=6):
    """Display only: decimal rendering (truncation) of an exact Fraction.
       (t44b_world.py, verbatim -- the sealed OUT/RESULT strings were printed with this.)"""
    x = Fraction(x)
    sign = "-" if x < 0 else ""
    x = abs(x)
    scaled = (x.numerator * 10 ** nd) // x.denominator
    return "%s%d.%0*d" % (sign, scaled // 10 ** nd, nd, scaled % 10 ** nd)


# =====================================================================================
# 1  VENUE ENTRY -- observation-entry gate (D-25 pattern; C-87 venue lineage)
# =====================================================================================
WORLD_VENUE_PROVENANCE = (
    "census access geometry: GR1 census grains, adjacency = shares a face (T42_C/T43_B "
    "world instruments); one walk step = one grain-boundary crossing = one unit of writer "
    "weight, the Gamma price (C-80/O-54); LANE_T44_B_WORLD, sealed")


class VenueGraph:
    """A declared venue graph: the walk lattice a Gamma-priced coupling lives on.
       adj: list over nodes of [(neighbor, multiplicity), ...] -- multiplicity kept
       exactly (a shared double edge counts twice, per the sealed lanes).  sector names
       the declared venue limit ("Z3"/"Z2"/"Z1") for the evidence instruments, or None.
       Constructed ONLY through venue() -- the provenance gate."""

    def __init__(self, name, adj, provenance, tier, sector):
        self.name, self.adj = name, adj
        self.provenance, self.tier, self.sector = provenance, tier, sector
        self.n = len(adj)
        self._mu_c = {}

    def row_sums(self):
        """Exact integer row sums of the adjacency (the Perron/Gershgorin input)."""
        return [sum(m for _, m in row) for row in self.adj]

    def degree(self):
        """The uniform degree, or None where row sums are not uniform (the Perron
           row-sum route then does not apply and mu_c_of DECLINES)."""
        rs = self.row_sums()
        return rs[0] if rs and all(r == rs[0] for r in rs) else None


def venue(name, adj, provenance=None, tier="world", sector=None):
    """THE OBSERVATION-ENTRY GATE (D-25, the URM pattern): a new venue graph -- a new
       record surface's access geometry -- enters HERE, as an adjacency structure with
       declared provenance.  REFUSES a world-tier venue without provenance; an exact
       idealisation (corner tier) must self-declare provenance='DEF-A'.

       CLAIM ROW: C-87 (venue lineage; LANE_T44_B_WORLD S0, LANE_T44_A_CORNER S0)."""
    if tier == "corner":
        if provenance != "DEF-A":
            raise ValueError(
                "classes layer REFUSES: a corner venue must self-declare provenance='DEF-A' "
                "-- the exact idealisation may never silently pose as the world (D-25).")
    else:
        if not provenance or not str(provenance).strip():
            raise ValueError(
                "classes layer REFUSES: a world-tier venue graph requires PROVENANCE -- the "
                "real access geometry it models and its source (D-25; the principal "
                "2026-08-20: the model is grounded in real record data, never the toy "
                "category).  Pass provenance=... naming the geometry's pinned source.")
    return VenueGraph(name, adj, provenance, tier, sector)


def torus3_adjacency(n):
    """Periodic n x n x n grain lattice (the census access geometry closed up for the
       universal-cover wrap gate): cells share a face <-> adjacent, multiplicity kept.
       CLAIM ROW: C-87/C-90 venue (LANE_T44_B_WORLD/t44b_lib.py, verbatim-in-substance)."""
    cells = [(x, y, z) for z in range(n) for y in range(n) for x in range(n)]
    idx = {c: i for i, c in enumerate(cells)}
    adj = [dict() for _ in cells]
    for i, (x, y, z) in enumerate(cells):
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx[((x + dx) % n, (y + dy) % n, (z + dz) % n)]
            adj[i][j] = adj[i].get(j, 0) + 1
    return cells, idx, [sorted(d.items()) for d in adj]


def world_venue(n):
    """The world venue (earned D = 3): the n^3 census-grain torus, entered through the
       gate with the PINNED provenance (D-25).  Returns (VenueGraph, cells, idx).
       CLAIM ROW: C-90 (LANE_T44_B_WORLD S0, sealed: deg-6, BFS == earned L1-with-wraps,
       torus counts == wrap-summed Z^3 counts)."""
    cells, idx, adj = torus3_adjacency(n)
    v = venue("world grain torus n=%d" % n, adj, provenance=WORLD_VENUE_PROVENANCE,
              tier="world", sector="Z3")
    return v, cells, idx


def corner_venue(Lx, Ly):
    """The corner venue (earned D = 2): the dual lattice of the model's own ported toric
       carrier -- plaquettes adjacent iff their supports SHARE A CARRIER EDGE, computed
       from the plaquette masks alone (geometry.Torus, the o54c conventions), multiplicity
       kept.  One walk step = one shared carrier edge = one unit of writer weight.
       DEF-A self-declared (exact idealisation).  Returns (VenueGraph, cells, idx).
       CLAIM ROW: C-87 (LANE_T44_A_CORNER/t44a_lib.py plaquette_adjacency,
       verbatim-in-substance on the model's ported carrier; sealed: row sums exactly 4
       on (4,6) and (3,7))."""
    import geometry as GE
    T = GE.Torus(Lx, Ly)
    cells = [(x, y) for y in range(Ly) for x in range(Lx)]
    idx = {c: i for i, c in enumerate(cells)}
    supp = {c: T.plaq(*c) for c in cells}
    adj = [[] for _ in cells]
    for i, c in enumerate(cells):
        for j in range(i + 1, len(cells)):
            common = GE.pc(supp[c] & supp[cells[j]])
            if common:
                adj[i].append((j, common))
                adj[j].append((i, common))
    v = venue("corner dual lattice (%d,%d)" % (Lx, Ly), adj, provenance="DEF-A",
              tier="corner", sector="Z2")
    return v, cells, idx


def chain_venue(L):
    """The D = 1 control venue: the cycle C_L (the chain's own dual), DEF-A self-declared.
       CLAIM ROW: C-87 discriminator (LANE_T44_A_CORNER/t44a_lib.py cycle_adjacency,
       verbatim; sealed: row sums exactly 2, mu_c = 1/2 the venue's own)."""
    adj = [[(((i + 1) % L), 1), (((i - 1) % L), 1)] for i in range(L)]
    return venue("chain cycle L=%d" % L, adj, provenance="DEF-A", tier="corner",
                 sector="Z1")


def walk_counts(v, src, K):
    """Exact integer walk counts N_k(src -> node), k = 0..K, transfer-matrix action on the
       venue's own lattice.  (t44b_lib.py walk_counts_adj, verbatim.)"""
    vec = [0] * v.n
    vec[src] = 1
    out = [list(vec)]
    for _ in range(K):
        new = [0] * v.n
        for i, vi in enumerate(vec):
            if vi:
                for j, mult in v.adj[i]:
                    new[j] += vi * mult
        vec = new
        out.append(list(vec))
    return out


def bfs_venue(v, src):
    """BFS distance on the venue graph == the earned separation instrument.
       (t44b_lib.py bfs_adj, verbatim.)"""
    dist = {src: 0}
    dq = _deque([src])
    while dq:
        i = dq.popleft()
        for j, _ in v.adj[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return dist


def venue_series(v, mu, src, targets, K):
    """The priced kernel ON THE DECLARED VENUE GRAPH ITSELF: exact partial sums
       S_K(src -> t) = sum_{k<=K} N_k mu^k (Fractions) for each target node, plus the
       exact geometric tail bound (deg mu)^{K+1}/(1 - deg mu) where the venue is
       deg-regular and deg*mu < 1; tail is None otherwise (on a finite venue the sum is
       always finite -- the class statements live in the venue limit; the honest note of
       LANE_T44_A_CORNER S7).  CLAIM ROW: C-87 (t44a_lib.py series_G_torus idiom)."""
    mu = Fraction(mu)
    counts = walk_counts(v, src, K)
    deg = v.degree()
    tail = None
    if deg is not None and deg * mu < 1:
        dm = deg * mu
        tail = (dm ** (K + 1)) / (1 - dm)
    out = {}
    for t in targets:
        s = Fraction(0)
        mp = Fraction(1)
        for k in range(K + 1):
            nk = counts[k][t]
            if nk:
                s += nk * mp
            mp *= mu
        out[t] = (s, tail)
    return out


# =====================================================================================
# 2  mu_c LOCATED -- the resolvent route, exact rationals (C-87/C-90)
# =====================================================================================
def resolvent_exact(adj, mu, src):
    """EXACT rational solve of (I - mu A) x = e_src on the finite venue (Fraction
       Gaussian elimination).  Returns list of Fractions, or None if singular.
       CLAIM ROW: C-87/C-90 mu_c route (LANE_T44_A_CORNER/t44a_lib.py, verbatim; reused
       by LANE_T44_B_WORLD S0)."""
    n = len(adj)
    mu = Fraction(mu)
    M = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = Fraction(1)
        for j, mult in adj[i]:
            M[i][j] -= mu * mult
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


def annihilates_constant(v, mu):
    """Computed check that (I - mu A) annihilates the constant vector on the venue
       (locates the resolvent pole AT mu = 1/deg exactly).  (t44a_lib.py
       kernel_at_quarter, generalized to the venue's own computed degree -- the literal
       1/4 was that lane's venue's own number, never this instrument's.)"""
    mu = Fraction(mu)
    for row in v.adj:
        s = Fraction(1)
        for _, mult in row:
            s -= mu * mult
        if s != 0:
            return False
    return True


def mu_c_of(v, certify="full", src=0):
    """LOCATE the venue's own critical point -- never a literal (D-8).
       Route (LANE_T44_B_WORLD S0/S2, the judge's N1 route; LANE_T44_A_CORNER S2):
         Perron row-sum: a deg-regular venue has A 1 = deg 1 exactly, spectral radius
         deg (Gershgorin ceiling + Perron attainment), so the positive walk series has
         radius of convergence exactly 1/deg: candidate mu_c = 1/deg.
         Resolvent certification (certify='full'): (I - mu_c A) annihilates the constant
         vector AND exact elimination meets a zero pivot AT mu_c (singular -- the D-15
         zero), while the venue is SOLVABLE at mu_c*(19/20) and mu_c*(21/20) (the
         positive controls beside the zero; on the sealed deg-6 venue these probe points
         are exactly the sealed 19/120 and 7/40).
       certify='rowsum' runs only the row-sum + annihilation tier (for large venues; the
       full tier is O(n^3) exact-rational elimination).
       DECLINES (located=False) where the venue is not degree-regular -- the ported route
       does not apply, and the model declines rather than returning a number."""
    if certify in v._mu_c:
        return v._mu_c[certify]
    deg = v.degree()
    if deg is None:
        out = dict(located=False,
                   why="venue is not degree-regular: the Perron row-sum route (the ported "
                       "sealed instrument) does not apply; the model declines")
        v._mu_c[certify] = out
        return out
    cand = Fraction(1, deg)
    pole = annihilates_constant(v, cand)
    out = dict(mu_c=cand, degree=deg, pole_at_candidate=pole, certify=certify)
    if certify == "full":
        singular = resolvent_exact(v.adj, cand, src) is None
        below = resolvent_exact(v.adj, cand * Fraction(19, 20), src) is not None
        above = resolvent_exact(v.adj, cand * Fraction(21, 20), src) is not None
        out.update(resolvent_singular_at_mu_c=singular, solvable_below=below,
                   solvable_above=above,
                   located=bool(pole and singular and below and above))
    else:
        out.update(located=bool(pole))
    v._mu_c[certify] = out
    return out


# =====================================================================================
# 3  THE PRICED SERIES G_mu(d) -- exact partial sums, proven tails (C-87)
# =====================================================================================
def W1(k, a):
    """Exact N_k^{Z}(0 -> a).  (t44b_lib.py, verbatim.)"""
    if (k + a) % 2 or abs(a) > k:
        return 0
    return comb(k, (k + a) // 2)


def N3_ref(k, a, b, c):
    """Reference axis-split evaluation of N_k^{Z^3}((a,b,c)) -- gated against brute-force
       DP before any measurement use.  Owner: standard combinatorics (Feller I).
       (t44b_lib.py, verbatim.)"""
    u, v = b + c, b - c
    tot = 0
    for k1 in range(abs(a), k + 1):
        w1 = W1(k1, a)
        if w1:
            k23 = k - k1
            w2, w3 = W1(k23, u), W1(k23, v)
            if w2 and w3:
                tot += comb(k, k1) * w1 * w2 * w3
    return tot


def dp3_counts(K, inside=None):
    """Brute-force DP walk counts on Z^3 (or restricted to `inside`): the gate reference.
       (t44b_lib.py, verbatim.)"""
    grid = {(0, 0, 0): 1}
    out = {}
    for k in range(K + 1):
        for v, c in grid.items():
            out[(k,) + v] = c
        new = {}
        for (x, y, z), c in grid.items():
            for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                key = (x + d[0], y + d[1], z + d[2])
                if inside is None or inside(key):
                    new[key] = new.get(key, 0) + c
        grid = new
    return out


def n3_even_row(m, a, b, c):
    """Exact N_{2m}^{Z^3}((a,b,c)) for an even-split target by the axis-split identity
       with ONE incremental running product; every update is integer multiply + EXACT
       integer divide (asserted).  CLAIM ROW: C-90 (t44b_lib.py, verbatim)."""
    assert a % 2 == 0 and (b + c) % 2 == 0
    s, t1, t2 = a // 2, (b + c) // 2, (b - c) // 2
    s, t1, t2 = abs(s), abs(t1), abs(t2)
    tmax = max(t1, t2)
    if s > m or tmax > m - s:
        return 0
    j1 = s
    r = m - j1
    T = comb(2 * m, 2 * j1) * comb(2 * r, r + t1) * comb(2 * r, r + t2)
    tot = T
    while j1 + 1 <= m - tmax:
        num = (r * r - t1 * t1) * (r * r - t2 * t2)
        den = (j1 + 1 + s) * (j1 + 1 - s) * (2 * r) * (2 * r - 1)
        q, rem = divmod(T * num, den)
        assert rem == 0
        T = q
        j1 += 1
        r -= 1
        tot += T
    return tot


def series_3d(mu, target, K):
    """Exact partial sum S_K = sum_{k<=K} N_k^{Z^3}(target) mu^k (Fraction) plus the
       exact geometric tail bound (6 mu)^{K+1}/(1-6 mu) (from N_k <= 6^k, the venue's
       row sum).  General-parity target.  CLAIM ROW: C-87 class (1) / C-90 subcritical
       rows (t44b_lib.py, verbatim)."""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    a, b, c = (abs(t) for t in target)
    u, v = b + c, abs(b - c)
    d = a + b + c
    num = 0
    ppow = 1
    for k in range(K + 1):
        if k:
            ppow *= p
        Nk = 0
        if k >= d and (k - d) % 2 == 0:
            k1, k23 = a, k - a
            T = comb(k, k1) * W1(k23, u) * W1(k23, v)
            Nk = T
            while k1 + 2 <= k - max(u, v):
                nmr = (((k23 + u) // 2) * ((k23 - u) // 2)
                       * (((k23 + v) // 2)) * ((k23 - v) // 2))
                dnm = ((k1 + a) // 2 + 1) * ((k1 - a) // 2 + 1) * k23 * (k23 - 1)
                qq, rem = divmod(T * nmr, dnm)
                assert rem == 0
                T = qq
                k1 += 2
                k23 -= 2
                Nk += T
        num = num * q + Nk * ppow
    S = Fraction(num, q ** K)
    six_mu = 6 * mu
    assert six_mu < 1
    tail = (six_mu ** (K + 1)) / (1 - six_mu)
    return S, tail


def ratio_interval(g_num, g_den):
    """[lo,hi] Fractions for G(t2)/G(t1) with G in [S, S+tail] on both sides.
       (t44a_lib.py, verbatim.)"""
    (s2, t2), (s1, t1) = g_num, g_den
    return (s2 / (s1 + t1), (s2 + t2) / s1)


def sqrt_bracket(fr, prec=10 ** 12):
    """Certified rational bracket lo <= sqrt(fr) <= hi via isqrt on scaled integers.
       COMPARISON USE ONLY; never on the measurement path.  (t44a_lib.py, verbatim.)"""
    fr = Fraction(fr)
    scaled = fr.numerator * prec * prec // fr.denominator
    r = isqrt(scaled)
    return Fraction(r, prec), Fraction(r + 2, prec)


def subcritical_row(mu, K, dmax):
    """One subcritical row of the D=3 sector: exact G_mu(d) partial sums d = 1..dmax at
       depth K, the ratio intervals, and the SEALED CLASS BOOLEANS (LANE_T44_B_WORLD S3,
       verbatim logic): every ratio upper bound <= 1 - MARGIN_LT1; Cauchy stabilization
       within CAUCHY_TOL at the top; power exclusion d(1-r) > POWER_SCOPE.  The class
       label is emitted by the booleans, never asserted.  The Ornstein-Zernike owner rate
       (cosh route, X = 1/(2mu) - 2) is returned as a labeled COMPARISON quantity only.
       CLAIM ROW: C-87 class (1) / C-90 subcritical rows."""
    mu = Fraction(mu)
    Gs = {d: series_3d(mu, (d, 0, 0), K) for d in range(1, dmax + 1)}
    rints = {d: ratio_interval(Gs[d + 1], Gs[d]) for d in range(1, dmax)}
    exp_ok = all(hi <= 1 - MARGIN_LT1 for lo, hi in rints.values())
    cau_ok = abs(rints[dmax - 1][0] - rints[dmax - 2][1]) <= CAUCHY_TOL
    qpow = (dmax - 1) * (1 - rints[dmax - 1][1])
    pow_ok = qpow > POWER_SCOPE
    cls = "EXPONENTIAL" if (exp_ok and cau_ok and pow_ok) else "UNRESOLVED"
    X = 1 / (2 * mu) - 2
    slo, shi = sqrt_bracket(X * X - 1)
    r_owner = (X - shi, X - slo)
    rl, rh = rints[dmax - 1]
    oz_within = max(rl - r_owner[1], r_owner[0] - rh) <= COMP_TOL
    return dict(mu=mu, K=K, dmax=dmax, cls=cls, G=Gs, ratios=rints,
                exp_ok=exp_ok, cauchy_ok=cau_ok, power_exclusion=qpow, pow_ok=pow_ok,
                comparison_owner_rate=r_owner, comparison_oz_within_tol=oz_within)


def divergence_witness(mu, m=WITNESS_M):
    """Supercritical witness (LANE_T44_B_WORLD S2, verbatim logic): the single term
       N_2m(0) mu^2m exceeds the half-depth term and the term ratio
       (N_{2(m+1)}(0)/N_2m(0)) mu^2 exceeds 1 + MARGIN_LT1 -- terms GROW, the series
       diverges term-by-term.  Computed booleans over exact rationals.
       CLAIM ROW: C-87 class (3)."""
    mu = Fraction(mu)
    NW1 = n3_even_row(m, 0, 0, 0)
    NW2 = n3_even_row(m + 1, 0, 0, 0)
    NH = n3_even_row(m // 2, 0, 0, 0)
    t_w = NW1 * mu ** (2 * m)
    t_h = NH * mu ** m
    ratio = Fraction(NW2, NW1) * mu * mu
    return dict(mu=mu, witness_m=m, term=t_w, half_term=t_h, term_ratio=ratio,
                grows=bool(ratio > 1 + MARGIN_LT1 and t_w > t_h))


# =====================================================================================
# 4  THE CRITICAL KERNEL -- the 1/d member (C-90) with certified tails
# =====================================================================================
def wallis_brackets(M0):
    """Rational c_l = M0*W_{M0}^2 and c_u = (M0+1/2)*W_{M0}^2, W_m = C(2m,m)/4^m, so that
       (lane-gated monotonicities) c_l/m <= W_m^2 <= c_u/(m+1/2) for all m >= M0.
       Owner: classical Wallis product.  (t44b_lib.py, verbatim.)"""
    W = Fraction(comb(2 * M0, M0), 4 ** M0)
    return M0 * W * W, (Fraction(2 * M0 + 1, 2)) * W * W


def max_trinomial(m):
    """Exact max multinomial coefficient m!/(j1! j2! j3!) (scan near the balanced split;
       window == full scan lane-gated to m = 60, verifier-gated to 80).
       (t44b_lib.py, verbatim.)"""
    best = 0
    lo = max(0, m // 3 - 2)
    for j1 in range(lo, m // 3 + 3):
        for j2 in range(max(0, (m - j1) // 2 - 2), (m - j1) // 2 + 3):
            j3 = m - j1 - j2
            if j3 < 0:
                continue
            v = comb(m, j1) * comb(m - j1, j2)
            best = max(best, v)
    return best


def q3_constant(m0):
    """Q3 rational with q_m := c_max(m)/3^m <= Q3/(m-2) for all m >= m0+2 (telescoping
       steps lane-gated).  (t44b_lib.py, verbatim.)"""
    assert m0 % 3 == 0
    q = Fraction(max_trinomial(m0), 3 ** m0)
    return (m0 * q) / (1 - Fraction(2, 3 * m0))


def sum_m52_bound(M):
    """Certified rational upper bound for sum_{m>M} (m-2)^{-5/2} (telescoping lemma
       4m^5 > (m-1)^3(2m+3)^2, lane-gated).  (t44b_lib.py, verbatim.)"""
    Mm = M - 2
    return Fraction(2, 3 * Mm * isqrt(Mm))


def sum_m32_bound(M):
    """Certified rational upper bound for sum_{m>M} (m-2)^{-3/2} (telescoping lemma
       4m^3 > (m-1)(2m+1)^2, lane-gated).  (t44b_lib.py, verbatim.)"""
    return Fraction(2, isqrt(M - 2))


def tail_constants():
    """Assemble the certified tail constants at the SEALED anchors M0_WAL = 240,
       M0_TRI = 300: c_l, c_u (Wallis), Q3 (max-trinomial), B5 = Q3 * sqrt_hi(c_u).
       Sealed displays: Q3 = 0.826995, B5 = 0.466824 (LANE_T44_B_WORLD S1).
       CLAIM ROW: C-90 tails."""
    c_l, c_u = wallis_brackets(M0_WAL)
    Q3 = q3_constant(M0_TRI)
    B5 = Q3 * sqrt_bracket(c_u)[1]
    return dict(c_l=c_l, c_u=c_u, Q3=Q3, B5=B5)


def diff_tail_bound(target, M, B5, RHO_=RHO, EDGE_C_=EDGE_C):
    """Certified tail of the regularized kernel beyond depth M (the assembled L3-3D
       bound; every ingredient lane-gated).  (t44b_lib.py, verbatim.)"""
    a, b, c = target
    s, t1, t2 = abs(a) // 2, abs(b + c) // 2, abs(b - c) // 2
    d2 = 4 * (s * s + t1 * t1 + t2 * t2)
    main = d2 * B5 * sum_m52_bound(M)
    Mm = M - 2
    edge = EDGE_C_ * B5 * Fraction(1, Mm * isqrt(Mm)) * (RHO_ ** (M + 1)) / (1 - RHO_)
    return main + edge


def abs_tail_bound(M, B5, p2M0):
    """Certified tail of the critical return series beyond depth 2M (even part
       p_2m(x) <= p_2m(0) <= B5 (m-2)^{-3/2}; odd part one-step decomposition).
       (t44b_lib.py, verbatim.)"""
    return 2 * B5 * sum_m32_bound(M) + p2M0


def gate_bound_L3(m, t):
    """The per-m assembled L3-3D difference bound (t44b_world.py gate_bound, verbatim):
       (N_2m(0)-N_2m(x))/N_2m(0) <= 4(s^2+t1^2+t2^2)/m + EDGE_C rho^m."""
    a, b, c = t
    s, t1, t2 = abs(a) // 2, abs(b + c) // 2, abs(b - c) // 2
    return Fraction(4 * (s * s + t1 * t1 + t2 * t2), m) + EDGE_C * RHO ** m


def crit_kernel_3d(targets, M, gate_range=None, gate_bound=None):
    """Regularized critical kernel on the Z^3 sector at the COMPUTED mu_c = 1/6 (located
       by mu_c_of on the venue; 36 = 6^2 is the venue's own degree squared), even-split
       targets:  a_M(x) = sum_{m<=M} (N_2m(0) - N_2m(x)) / 36^m  (EXACT rational).
       Returns (kernels, S_M(0), p_2M(0), gate_ok).  The L1-3D domination N_2m(x) <=
       N_2m(0) is ASSERTED en route.  CLAIM ROW: C-90 (t44b_lib.py, verbatim)."""
    acc = {t: 0 for t in targets}
    acc0 = 0
    gate_ok = True
    p2M0 = None
    for m in range(M + 1):
        N0 = n3_even_row(m, 0, 0, 0)
        acc0 = acc0 * 36 + N0
        for t in targets:
            Nt = n3_even_row(m, *t)
            assert Nt <= N0
            acc[t] = acc[t] * 36 + (N0 - Nt)
            if gate_range and gate_range[0] <= m <= gate_range[1]:
                lhs = Fraction(N0 - Nt, 36 ** m)
                rhs = Fraction(N0, 36 ** m) * gate_bound(m, t)
                if lhs > rhs:
                    gate_ok = False
        if m == M:
            p2M0 = Fraction(N0, 36 ** m)
    den = 36 ** M
    out = {t: Fraction(acc[t], den) for t in targets}
    return out, Fraction(acc0, den), p2M0, gate_ok


def kernel_pass(stops, M, gate_range=None, gate_bound=None, snapshots=()):
    """ONE combined exact pass of the critical kernel with PER-TARGET stop depths and
       snapshot depths -- the sealed verifier's deepening loop
       (LANE_T44_B_WORLD/VERIFY/adv_verify.py section E, verbatim-in-substance; register
       row C-90 carries its M=2800 coefficient bracket).  Identical accumulation to
       crit_kernel_3d (gated equal at spot depth by the check block -- one instrument,
       two drivers, never a shortcut).
         stops: {target: M_t} with M_t <= M; the origin row always runs to M.
       Returns dict(ker={t: a_{M_t}(t)}, ker_at={snap: {t: a_snap(t)} for active t},
                    s0={depth: S_depth(0)}, p2m0={depth: p_2depth(0)}, asm_ok=bool)
       with s0/p2m0 recorded at every snapshot and at M.  CLAIM ROW: C-90."""
    targets = list(stops)
    acc = {t: 0 for t in targets}
    acc0 = 0
    gate_ok = True
    snapset = set(snapshots)
    ker, ker_at, s0, p2m0 = {}, {}, {}, {}
    for m in range(M + 1):
        N0 = n3_even_row(m, 0, 0, 0)
        acc0 = acc0 * 36 + N0
        for t in targets:
            if m > stops[t]:
                continue
            Nt = n3_even_row(m, *t)
            assert Nt <= N0
            acc[t] = acc[t] * 36 + (N0 - Nt)
            if gate_range and gate_range[0] <= m <= gate_range[1]:
                lhs = Fraction(N0 - Nt, 36 ** m)
                rhs = Fraction(N0, 36 ** m) * gate_bound(m, t)
                if lhs > rhs:
                    gate_ok = False
            if m == stops[t]:
                ker[t] = Fraction(acc[t], 36 ** m)
        if m in snapset or m == M:
            den = 36 ** m
            s0[m] = Fraction(acc0, den)
            p2m0[m] = Fraction(N0, den)
            if m in snapset:
                ker_at[m] = {t: Fraction(acc[t], den) for t in targets if stops[t] >= m}
    return dict(ker=ker, ker_at=ker_at, s0=s0, p2m0=p2m0, asm_ok=gate_ok)


def increment_interval(ker, tails, t1, t2):
    """Interval for the kernel increment H = a(t2) - a(t1) = G(t1) - G(t2) from exact
       kernel values plus certified tails.  (t44b_world.py H_int, verbatim.)"""
    return (ker[t2] - ker[t1] - tails[t1], ker[t2] + tails[t2] - ker[t1])


def doubling_ratio(ker, tails, pair_lo, pair_hi):
    """Ratio interval of consecutive doubling increments H(pair_hi)/H(pair_lo), the class
       discriminant: INV window [2/5,3/5] is the 1/d power law (increments halve), LOG
       window [4/5,5/4] the marginal class (constant), LIN window [9/5,11/5] the linear
       class (double).  pair_* are (t1, t2) target pairs.  CLAIM ROW: C-90."""
    H_lo = increment_interval(ker, tails, *pair_lo)
    H_hi = increment_interval(ker, tails, *pair_hi)
    return (H_hi[0] / H_lo[1], H_hi[1] / H_lo[0])


STAB_TOL = Fraction(1, 8)       # evidence tier: kernel deepening M/2 -> M must move each
#                                 value by <= STAB_TOL of the smallest increment used


def critical_evidence_3d(M=350):
    """Cheap self-contained critical-row EVIDENCE on the Z^3 sector (the class_verdict
       evidence tier): exact kernel at depth M for the axis targets (2),(4),(8), the
       doubling-increment ratio H(4->8)/H(2->4) from the POINT kernel values, and the
       window booleans -- in INV [2/5,3/5], outside LOG and LIN (the 1/d-class signature:
       increments halve).  Truncation honesty at this tier is a COMPUTED STABILIZATION
       gate, not a certified tail: the deepening M/2 -> M must move every kernel value by
       <= STAB_TOL of the smallest increment used (the certified-interval statements --
       and the exponent-bracket-contains-1 claim, which needs the deepest pair -- live at
       the sealed depths M = 1400/2800, gated in checks_classes.py; the lane's own
       certified tails are only honest for M >~ 1000, where the edge term RHO^M has
       decayed).  CLAIM ROW: C-90 (evidence tier)."""
    tgts = [(2, 0, 0), (4, 0, 0), (8, 0, 0)]
    kp = kernel_pass({t: M for t in tgts}, M, snapshots=(M // 2,))
    ker = kp["ker"]
    half = kp["ker_at"][M // 2]
    H24 = ker[(4, 0, 0)] - ker[(2, 0, 0)]
    H48 = ker[(8, 0, 0)] - ker[(4, 0, 0)]
    stab = max(abs(ker[t] - half[t]) for t in tgts)
    stable = stab <= STAB_TOL * min(H24, H48)
    r24 = H48 / H24
    in_inv = INV_LO <= r24 <= INV_HI
    out_log = r24 < LOG_LO
    out_lin = r24 < LIN_LO
    return dict(M=M, H_24=H24, H_48=H48, ratio_24=r24, stabilization=stab,
                stabilized=bool(stable), in_inv_window=bool(in_inv),
                outside_log_window=bool(out_log), outside_lin_window=bool(out_lin),
                inv_class_signature=bool(stable and in_inv and out_log and out_lin))


# =====================================================================================
# 5  CROSS-DIMENSION DISCRIMINATOR -- one instrument, three venues (C-87)
# =====================================================================================
class BinStepper:
    """Exact incremental binomial C(k,(k+u)/2) along k = |u|, |u|+2, ...; every update is
       integer multiply + EXACT integer divide (asserted).  (t44a_lib.py, verbatim.)"""

    def __init__(self, u):
        self.u = abs(u)
        self.k = self.u
        self.val = 1

    def step2(self):
        k, u = self.k, self.u
        num = self.val * (k + 1) * (k + 2)
        d1 = (k + u) // 2 + 1
        d2 = (k - u) // 2 + 1
        q, r = divmod(num, d1 * d2)
        assert r == 0
        self.val = q
        self.k = k + 2
        return q


def series_target_2d(mu, a, b, K):
    """Exact S_K = sum_{k<=K} N_k^{Z^2}((a,b)) mu^k plus the exact geometric tail
       (4mu)^{K+1}/(1-4mu).  Rotation bijection owner: Feller I (lane-gated against DP).
       CLAIM ROW: C-87 D=2 rows (t44a_lib.py, verbatim)."""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    u, v = a + b, a - b
    assert (u - v) % 2 == 0
    k0 = max(abs(u), abs(v))
    su, sv = BinStepper(u), BinStepper(v)
    while su.k < k0:
        su.step2()
    while sv.k < k0:
        sv.step2()
    k = k0
    num = su.val * sv.val * (p ** k)
    scale_k = k
    p2 = p * p
    q2 = q * q
    ppow = p ** k
    while k + 2 <= K:
        su.step2()
        sv.step2()
        k += 2
        ppow *= p2
        num = num * q2 + su.val * sv.val * ppow
        scale_k = k
    S = Fraction(num, q ** scale_k)
    four_mu = 4 * mu
    assert four_mu < 1
    tail = (four_mu ** (K + 1)) / (1 - four_mu)
    return S, tail


def series_target_1d(mu, d, K):
    """Same on the 1D chain; tail from N_k <= 2^k.  (t44a_lib.py, verbatim.)"""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    s = BinStepper(d)
    k = s.k
    num = s.val * (p ** k)
    scale_k = k
    p2, q2 = p * p, q * q
    ppow = p ** k
    while k + 2 <= K:
        s.step2()
        k += 2
        ppow *= p2
        num = num * q2 + s.val * ppow
        scale_k = k
    S = Fraction(num, q ** scale_k)
    two_mu = 2 * mu
    assert two_mu < 1
    tail = (two_mu ** (K + 1)) / (1 - two_mu)
    return S, tail


def partial_return_sum_2d(Klist):
    """S_K(0,0) at the D=2 venue's own mu_c = 1/4, exact at each K in Klist (divergence
       witness data for the D=2 critical row).  (t44a_lib.py, verbatim.)"""
    Kmax = max(Klist)
    want = set(Klist)
    out = {}
    acc = 1
    c = BinStepper(0)
    k = 0
    if 0 in want:
        out[0] = Fraction(acc, 1)
    while k + 2 <= Kmax:
        c.step2()
        k += 2
        acc = acc * 16 + c.val * c.val
        if k in want or (k + 1) in want:
            out[k] = Fraction(acc, 4 ** k)
    return out


def potential_kernel_2d(dlist_axis, dlist_diag, K):
    """Regularized critical kernel on Z^2 at the D=2 venue's own computed mu_c = 1/4,
       with the PROVEN tail bound (u^2+v^2)/(4K).  CLAIM ROW: C-87 class (2), D=2 member
       (t44a_lib.py, verbatim; owners Stohr 1950 / Spitzer for the comparison only)."""
    targets = [('axis', d, d, d) for d in dlist_axis] + \
              [('diag', c, 2 * c, 0) for c in dlist_diag]
    acc = {}
    steppers = {}
    for kind, lab, u, v in targets:
        acc[(kind, lab)] = 1
        steppers[(kind, lab)] = None
    c0 = BinStepper(0)
    k = 0
    while k + 2 <= K:
        c0.step2()
        k += 2
        c0sq = c0.val * c0.val
        for kind, lab, u, v in targets:
            key = (kind, lab)
            st = steppers[key]
            if st is None and max(abs(u), abs(v)) <= k:
                su, sv = BinStepper(u), BinStepper(v)
                while su.k < k:
                    su.step2()
                while sv.k < k:
                    sv.step2()
                st = (su, sv)
                steppers[key] = st
            elif st is not None:
                st[0].step2()
                st[1].step2()
            nk = st[0].val * st[1].val if st is not None else 0
            acc[key] = acc[key] * 16 + (c0sq - nk)
    den = 4 ** K
    out = {}
    for kind, lab, u, v in targets:
        tail = Fraction(u * u + v * v, 4 * K)
        out[(kind, lab)] = (Fraction(acc[(kind, lab)], den), tail)
    return out


def potential_kernel_1d(dlist, K):
    """Regularized critical kernel on Z at the chain venue's own computed mu_c = 1/2,
       with the PROVEN tail bound d^2/(2 isqrt(2M)).  CLAIM ROW: C-87 class (2), D=1
       member (t44a_lib.py, verbatim; owner comparison: Spitzer a_Z(d) = |d| exactly)."""
    acc = {d: 1 for d in dlist}
    steppers = {d: None for d in dlist}
    c0 = BinStepper(0)
    k = 0
    while k + 2 <= K:
        c0.step2()
        k += 2
        for d in dlist:
            st = steppers[d]
            if st is None and d <= k:
                st = BinStepper(d)
                while st.k < k:
                    st.step2()
                steppers[d] = st
            elif st is not None:
                st.step2()
            nk = st.val if st is not None else 0
            acc[d] = acc[d] * 4 + (c0.val - nk)
    den = 2 ** K
    M = K // 2
    out = {}
    for d in dlist:
        tail = Fraction(d * d, 2 * isqrt(2 * M))
        out[d] = (Fraction(acc[d], den), tail)
    return out


def discriminator(K2=6000, K1=80000):
    """THE D-15 CROSS-DIMENSION DISCRIMINATOR (LANE_T44_B_WORLD S5, verbatim logic):
       one instrument (regularized critical kernel, doubling-increment ratios), three
       venues, three pairwise-DISJOINT declared windows.  Computes the D=2 kernel at its
       own mu_c = 1/4 (depth K2) and the D=1 kernel at its own mu_c = 1/2 (depth K1) and
       places each ratio in its window: D=2 LOG (constant increments), D=1 LIN (doubling
       increments).  The D=3 INV placement is the critical_evidence_3d / sealed-kernel
       side.  CLAIM ROW: C-87 (the critical class is the venue's own)."""
    K2d = potential_kernel_2d([2, 4, 8, 16], [], K2)
    I2a = (K2d[("axis", 4)][0] - K2d[("axis", 2)][0] - K2d[("axis", 2)][1],
           K2d[("axis", 4)][0] + K2d[("axis", 4)][1] - K2d[("axis", 2)][0])
    I2b = (K2d[("axis", 8)][0] - K2d[("axis", 4)][0] - K2d[("axis", 4)][1],
           K2d[("axis", 8)][0] + K2d[("axis", 8)][1] - K2d[("axis", 4)][0])
    I2c = (K2d[("axis", 16)][0] - K2d[("axis", 8)][0] - K2d[("axis", 8)][1],
           K2d[("axis", 16)][0] + K2d[("axis", 16)][1] - K2d[("axis", 8)][0])
    r2 = (I2b[0] / I2a[1], I2b[1] / I2a[0])
    r2b = (I2c[0] / I2b[1], I2c[1] / I2b[0])
    K1d = potential_kernel_1d([2, 4, 8], K1)
    I1a = (K1d[4][0] - K1d[2][0] - K1d[2][1], K1d[4][0] + K1d[4][1] - K1d[2][0])
    I1b = (K1d[8][0] - K1d[4][0] - K1d[4][1], K1d[8][0] + K1d[8][1] - K1d[4][0])
    r1 = (I1b[0] / I1a[1], I1b[1] / I1a[0])
    d2_log = LOG_LO <= r2[0] and r2[1] <= LOG_HI and LOG_LO <= r2b[0] and r2b[1] <= LOG_HI
    d1_lin = LIN_LO <= r1[0] and r1[1] <= LIN_HI
    disjoint = INV_HI < LOG_LO < LOG_HI < LIN_LO
    return dict(K2=K2, K1=K1, ratio_D2=r2, ratio_D2_deeper=r2b, ratio_D1=r1,
                D2_in_log_window=bool(d2_log), D1_in_lin_window=bool(d1_lin),
                windows_disjoint=bool(disjoint))


# =====================================================================================
# 6  THE CLASS VERDICT -- the computed boolean triple (C-87)
# =====================================================================================
def class_verdict(v, mu, evidence=False, K=160, dmax=12, M_crit=350):
    """THE REACHABLE-CLASS VERDICT for a declared coupling price mu on a declared venue:
       a COMPUTED BOOLEAN TRIPLE (exponential, critical, divergent) by exact rational
       comparison of mu against THE VENUE'S OWN mu_c, located by mu_c_of (the resolvent
       route) -- never a literal critical point (D-8).  Exactly one of the three booleans
       is True whenever mu_c is located.

       evidence=True additionally runs the class's own series-level instrument on the
       venue's DECLARED sector (the venue limit; sector "Z3" fully, "Z2"/"Z1" via the
       discriminator kernels): subcritical_row booleans below mu_c, critical_evidence
       window booleans at mu_c, divergence_witness above.  Evidence depths (K, dmax,
       M_crit) are the caller's declared instrument settings; the sealed full-depth
       numbers are gated in checks_classes.py.  DECLINES (applies=False) where the venue
       is not degree-regular.  CLAIM ROW: C-87 (taxonomy; LANE_T44_B_WORLD)."""
    mu = Fraction(mu)
    loc = mu_c_of(v, certify="full" if v.n <= 100 else "rowsum")
    if not loc.get("located"):
        return dict(applies=False, why=loc.get("why", "mu_c not located"), mu=mu)
    mc = loc["mu_c"]
    out = dict(applies=True, mu=mu, mu_c=mc, mu_c_located=loc,
               exponential=bool(mu < mc), critical=bool(mu == mc),
               divergent=bool(mu > mc), evidence=None)
    if evidence:
        if out["divergent"]:
            out["evidence"] = divergence_witness(mu)
        elif v.sector == "Z3":
            out["evidence"] = (critical_evidence_3d(M_crit) if out["critical"]
                               else subcritical_row(mu, K, dmax))
        elif v.sector == "Z2":
            if out["critical"]:
                out["evidence"] = discriminator(K2=2000, K1=20000)
            else:
                rows = {d: series_target_2d(mu, d, 0, K) for d in (6, 8, 10)}
                rr = ratio_interval(rows[10], rows[8])
                out["evidence"] = dict(ratio=rr,
                                       exp_ok=bool(rr[1] <= 1 - MARGIN_LT1))
        elif v.sector == "Z1":
            if out["critical"]:
                out["evidence"] = discriminator(K2=2000, K1=20000)
            else:
                rows = {d: series_target_1d(mu, d, K) for d in (6, 8, 10)}
                rr = ratio_interval(rows[10], rows[8])
                out["evidence"] = dict(ratio=rr,
                                       exp_ok=bool(rr[1] <= 1 - MARGIN_LT1))
    return out
