"""T48-C shared library -- exact machinery for PROBE 3 (O-58 N2, world venue):
the induced per-link amplitude of one writer step, computed from the world ensemble
the model actually has.

VENUE: rebuilt per LANE_T44_B_WORLD/PUBLISHED_CONVENTIONS.txt -- the census access
geometry's grain lattice (GR1 grains, adjacency = shares a face; one walk step = one
grain-boundary crossing).  The venue is rebuilt fresh in-lane and every identity used
(block-interior DP, torus wrap identity, axis-split binomial identity) is gated against
brute force in the driver BEFORE any kernel is built on it.

EVERYTHING ON THE MEASUREMENT PATH IS EXACT: python ints and Fractions only.  Floats
appear only in the driver's labeled COMPARISON section, run AFTER results are computed.

OWNERS of standard mathematics used (attributed; applied only to the named program
variables, see D24_AUDIT.txt):
  - Perron-Frobenius / row-sum spectral radius of a nonnegative matrix: Perron 1907,
    Frobenius 1912 (applied to the writer-step kernel W and to the directed-edge
    non-backtracking operator B only through exact row/column sums on the venue).
  - Markov kernels, stationary measures, detailed balance: standard (Kolmogorov
    criterion); applied only to the two-state record kernel the model itself declares
    (project_model.py activation convention, corrected form) and to its lattice lifts.
  - Trace preservation of dilation unitaries (iv' writers are CPTP on the record):
    Stinespring 1955 -- attribution for WHY row sums 1 is the structural property of an
    energy-conserving dilation writer; the row sums themselves are computed, never cited.
  - Non-backtracking (directed-edge) operator: Hashimoto 1989 -- used only as the
    criticality REFERENCE for the extension-only ensemble, earned in-lane by the same
    row-sum instrument.
  - Persistent random walks renormalize diffusivity, not class (comparison-only remark
    in the driver): Goldstein 1951, Kac 1974.
  - Axis-split identity / rotation bijection for free Z^3 walk counts: standard
    combinatorics (Feller I); GATED against brute-force DP here before use.
"""
from fractions import Fraction
from math import comb
from collections import deque

DIRS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
OPP = (1, 0, 3, 2, 5, 4)


# ---------------------------------------------------------------- venue rebuild
def torus3(n):
    """Periodic n x n x n grain lattice (census access geometry closed up for the wrap
    gate): cells share a face <-> adjacent.  Returns cells, index, nbr[i][d]."""
    cells = [(x, y, z) for z in range(n) for y in range(n) for x in range(n)]
    idx = {c: i for i, c in enumerate(cells)}
    nbr = []
    for (x, y, z) in cells:
        row = []
        for (dx, dy, dz) in DIRS:
            row.append(idx[((x + dx) % n, (y + dy) % n, (z + dz) % n)])
        nbr.append(row)
    return cells, idx, nbr


def bfs(nbr, src):
    dist = {src: 0}
    dq = deque([src])
    while dq:
        i = dq.popleft()
        for j in nbr[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return dist


def l1_wrap(c, n):
    return sum(min(v, n - v) for v in c)


def torus_counts(nbr, src, K):
    """Exact integer walk counts on the venue, transfer-matrix action, k = 0..K."""
    vec = [0] * len(nbr)
    vec[src] = 1
    out = [list(vec)]
    for _ in range(K):
        new = [0] * len(nbr)
        for i, v in enumerate(vec):
            if v:
                for j in nbr[i]:
                    new[j] += v
        vec = new
        out.append(list(vec))
    return out


def dp3(K, inside=None):
    """Brute-force DP walk counts on Z^3 (or restricted to predicate `inside`):
    dict {(k, x, y, z): N}.  Gate reference; exact ints."""
    grid = {(0, 0, 0): 1}
    out = {}
    for k in range(K + 1):
        for v, c in grid.items():
            out[(k,) + v] = c
        if k == K:
            break
        new = {}
        for (x, y, z), c in grid.items():
            for (dx, dy, dz) in DIRS:
                key = (x + dx, y + dy, z + dz)
                if inside is None or inside(key):
                    new[key] = new.get(key, 0) + c
        grid = new
    return out


def W1(k, a):
    """Exact N_k^{Z}(0 -> a)."""
    a = abs(a)
    if (k + a) % 2 or a > k:
        return 0
    return comb(k, (k + a) // 2)


def N3_free(k, a, b, c):
    """Free Z^3 walk count by the axis-split identity (gated against DP in driver)."""
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


# ---------------------------------------------------------------- writer-step kernels
def kernel_pos(nbr, diag, link):
    """Position-space one-writer-step kernel on the venue: W[i][i] = diag,
    W[i][j] = link on each face-neighbor j.  All Fractions."""
    W = []
    for i, row in enumerate(nbr):
        r = {i: Fraction(diag)}
        for j in row:
            r[j] = r.get(j, 0) + Fraction(link)
        W.append(r)
    return W


def kernel_edge(nbr, stay, fresh, back):
    """Directed-edge (tip-state) kernel: state s = (i, d) = tip at grain i, last step
    along DIRS[d].  Per attempt: continue fresh e != OPP[d] with amplitude `fresh` to
    (nbr[i][e], e); retreat (tip erases, string backtracks) with amplitude `back` to
    (nbr[i][OPP[d]], OPP[d]); stay with `stay`.  back=None DROPS the retreat channel
    (the extension-only ensemble: that measure is lost).  Direction marker after any
    move = the direction just traveled (declared convention)."""
    S = []
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            r = {s: Fraction(stay)}
            for e in range(6):
                t = row[e] * 6 + e
                if e == OPP[d]:
                    if back is not None:
                        r[t] = r.get(t, 0) + Fraction(back)
                else:
                    r[t] = r.get(t, 0) + Fraction(fresh)
            S.append(r)
    return S


def nb_edge_adjacency(nbr):
    """The venue's non-backtracking directed-edge adjacency (0/1): from (i, d) to
    (nbr[i][e], e) for every e != OPP[d].  Criticality reference for the
    extension-only ensemble, earned by row sums (driver gate)."""
    S = []
    for i, row in enumerate(nbr):
        for d in range(6):
            r = {}
            for e in range(6):
                if e != OPP[d]:
                    t = row[e] * 6 + e
                    r[t] = r.get(t, 0) + 1
            S.append(r)
    return S


# ---------------------------------------------------------------- exact linear algebra helpers
def row_sums(W):
    return [sum(r.values()) for r in W]


def col_sums(W, n):
    cs = [Fraction(0)] * n
    for r in W:
        for j, v in r.items():
            cs[j] += v
    return cs


def apply_forward(W, vec):
    """measure evolution: new[j] = sum_i vec[i] W[i][j] (exact)."""
    new = [Fraction(0)] * len(W)
    for i, v in enumerate(vec):
        if v:
            for j, w in W[i].items():
                new[j] += v * w
    return new


def extract_pos(W, nbr):
    """Induced per-link amplitude of the summed writer propagator, position space:
    I - W = c (I - sum_e m_e S_e); returns (set of c across cells,
    set of m across all links/directions).  Uniformity is GATED downstream,
    never inserted."""
    cs, ms = set(), set()
    for i, r in enumerate(W):
        c = 1 - r[i]
        cs.add(c)
        for j in nbr[i]:
            ms.add(r[j] / c)
    return cs, ms


def extract_edge(W, nbr):
    """Induced per-link amplitudes, directed-edge space: for each state, c = 1 - stay;
    returns (set of c, set of m_fresh, set of m_back, set of per-state totals)."""
    cs, mf, mb, tots = set(), set(), set(), set()
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            r = W[s]
            c = 1 - r.get(s, Fraction(0))
            cs.add(c)
            tot = Fraction(0)
            for e in range(6):
                t = row[e] * 6 + e
                a = r.get(t, Fraction(0))
                if e == OPP[d]:
                    if a:
                        mb.add(a / c)
                    tot += a / c
                else:
                    mf.add(a / c)
                    tot += a / c
            tots.add(tot)
    return cs, mf, mb, tots


def op_identity_pos(W, nbr, c, m):
    """Entrywise gate: I - W == c (I - m A) on the venue."""
    for i, r in enumerate(W):
        # diagonal
        if 1 - r[i] != c:
            return False
        # off-diagonal
        row_counts = {}
        for j in nbr[i]:
            row_counts[j] = row_counts.get(j, 0) + 1
        for j, mult in row_counts.items():
            if -r[j] != -c * m * mult:
                return False
        # no other entries
        if set(r.keys()) - {i} - set(row_counts.keys()):
            return False
    return True


def drift_per_state(W, nbr):
    """Exact one-step expected displacement of each edge state (local step vectors,
    not torus coordinates).  Returns dict s -> (Fraction, Fraction, Fraction)."""
    out = {}
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            v = [Fraction(0)] * 3
            for t, a in W[s].items():
                if t != s:
                    e = t % 6
                    for ax in range(3):
                        v[ax] += a * DIRS[e][ax]
            out[s] = tuple(v)
    return out
