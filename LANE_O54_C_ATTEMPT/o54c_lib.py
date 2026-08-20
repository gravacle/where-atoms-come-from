"""O54-C shared library -- exact F_2 machinery for the ADVERSARIAL PURE-GAMMA FALLOFF ATTEMPT.

Everything exact: Paulis as 2n-bit python ints (low n bits = x part, high n bits = z part),
F_2 ranks by bitmask elimination, coset minima by exhaustive numpy span scans.
No floats anywhere on the measurement path.  Every PASS/FAIL in the driver is gated by a
computed boolean.

EARNED-INSTRUMENT REUSE (cited, not recomputed):
  - Pauli mask conventions and weight/pairing functions follow LANE_T42_A_DISTANCE/t42_lib.py
    (C-78 instrument lineage).
  - Stabilizer-state region entropy S(A) = |A| - dim S_A: BORROWED formula,
    owner Fattal-Cubitt-Yamamoto-Bravyi-Chuang, quant-ph/0406168 (2004); perimeter-law
    instance Hamma-Ionicioiu-Zanardi 2005.
"""
import sys
from collections import deque

import numpy as np

MODEL = "/Users/bgm/MB Work/where-atoms-come-from/model"
if MODEL not in sys.path:
    sys.path.insert(0, MODEL)


def pc(x):
    """popcount of a python int (py3.9: no int.bit_count)."""
    return bin(x).count("1")


# ------------------------------------------------------------------ F_2 rank on bitmasks
def rank_f2(rows):
    """Rank over F_2 of the row set (python int masks). Exact."""
    piv = {}
    r = 0
    for m in rows:
        while m:
            t = m.bit_length() - 1
            if t in piv:
                m ^= piv[t]
            else:
                piv[t] = m
                r += 1
                break
    return r


# ------------------------------------------------------------------ Pauli helpers (2n-bit)
def supp_mask(g, n):
    """Qubit-support mask (n bits): union of x and z supports."""
    return (g & ((1 << n) - 1)) | (g >> n)


def weight_xz(g, n):
    return pc(supp_mask(g, n))


def sp_pair(u, v, n):
    """Symplectic F_2 pairing (0 = commute, 1 = anticommute)."""
    mask = (1 << n) - 1
    xu, zu = u & mask, u >> n
    xv, zv = v & mask, v >> n
    return (pc(xu & zv) + pc(zu & xv)) % 2


# ------------------------------------------------------------------ toric carrier
class Torus:
    """L_x x L_y toric code. Qubits on edges, n = 2*Lx*Ly.
    h(x,y) = edge (x,y)-(x+1,y), index y*Lx+x.
    v(x,y) = edge (x,y)-(x,y+1), index Lx*Ly + y*Lx + x.
    Star (X-type) at vertex (x,y); plaquette (Z-type) at cell (x,y)."""

    def __init__(self, Lx, Ly):
        self.Lx, self.Ly = Lx, Ly
        self.n = 2 * Lx * Ly

    def h(self, x, y):
        return (y % self.Ly) * self.Lx + (x % self.Lx)

    def v(self, x, y):
        return self.Lx * self.Ly + (y % self.Ly) * self.Lx + (x % self.Lx)

    def star(self, x, y):
        """X-type: mask in low n bits."""
        m = 0
        for e in (self.h(x, y), self.h(x - 1, y), self.v(x, y), self.v(x, y - 1)):
            m |= 1 << e
        return m

    def plaq(self, x, y):
        """Z-type: mask in high n bits."""
        m = 0
        for e in (self.h(x, y), self.h(x, y + 1), self.v(x, y), self.v(x + 1, y)):
            m |= 1 << (self.n + e)
        return m

    def plaq_edges(self, x, y):
        return [self.h(x, y), self.h(x, y + 1), self.v(x, y), self.v(x + 1, y)]

    def all_stars(self):
        return [self.star(x, y) for y in range(self.Ly) for x in range(self.Lx)]

    def all_plaqs(self):
        return [self.plaq(x, y) for y in range(self.Ly) for x in range(self.Lx)]

    def zbar1(self):
        """Z loop: vertical edges at x=0 (all y)."""
        m = 0
        for y in range(self.Ly):
            m |= 1 << (self.n + self.v(0, y))
        return m

    def zbar2(self):
        """Z loop: horizontal edges at y=0 (all x)."""
        m = 0
        for x in range(self.Lx):
            m |= 1 << (self.n + self.h(x, 0))
        return m

    def xbar1(self):
        """X dual loop: horizontal edges at x=0 (all y). Conjugate to zbar2."""
        m = 0
        for y in range(self.Ly):
            m |= 1 << self.h(0, y)
        return m

    def xbar2(self):
        """X dual loop: vertical edges at y=0 (all x). Conjugate to zbar1."""
        m = 0
        for x in range(self.Lx):
            m |= 1 << self.v(x, 0)
        return m

    def dual_path_x(self, u, vv):
        """X-type connector for removed plaquettes u=(ux,uy) -> vv=(vx,vy): dual path,
        x first then y, each step crossing the shared edge. Returns 2n-bit mask (x part)."""
        (ux, uy), (vx, vy) = u, vv
        m = 0
        x, y = ux, uy
        while x != vx:
            nx = (x + 1) % self.Lx
            # step (x,y)->(x+1,y) crosses v(x+1,y)
            m ^= 1 << self.v(x + 1, y)
            x = nx
        while y != vy:
            ny = (y + 1) % self.Ly
            # step (x,y)->(x,y+1) crosses h(x,y+1)
            m ^= 1 << self.h(x, y + 1)
            y = ny
        return m

    def patch(self, x0, y0, w):
        """Region: all h and v edges with x in [x0,x0+w), y in [y0,y0+w) (mod wrap)."""
        qs = set()
        for dx in range(w):
            for dy in range(w):
                qs.add(self.h(x0 + dx, y0 + dy))
                qs.add(self.v(x0 + dx, y0 + dy))
        return sorted(qs)

    def band(self, x0, b):
        """Non-contractible band: all edges with x in [x0,x0+b), all y."""
        qs = set()
        for dx in range(b):
            for y in range(self.Ly):
                qs.add(self.h(x0 + dx, y))
                qs.add(self.v(x0 + dx, y))
        return sorted(qs)


# ------------------------------------------------------------------ entropy / MI
def region_entropy(gens, G, n, region):
    """S(A) = |A| - dim S_A for a PURE stabilizer state with independent-generator count
    G == n (caller gates this).  dim S_A = G - rank(gens restricted to complement).
    BORROWED formula (Fattal et al. 2004). Exact integers."""
    rs = set(region)
    comp = [q for q in range(n) if q not in rs]
    restricted = []
    for g in gens:
        out = 0
        for i, q in enumerate(comp):
            if (g >> q) & 1:
                out |= 1 << (2 * i)
            if (g >> (n + q)) & 1:
                out |= 1 << (2 * i + 1)
        restricted.append(out)
    dim_SA = G - rank_f2(restricted)
    return len(rs) - dim_SA


def mutual_info(gens, G, n, A, B):
    """I(A:B) = S_A + S_B - S_AB, exact integer (in bits, i.e. log2 units x1)."""
    SA = region_entropy(gens, G, n, A)
    SB = region_entropy(gens, G, n, B)
    SAB = region_entropy(gens, G, n, sorted(set(A) | set(B)))
    return SA, SB, SAB, SA + SB - SAB


def straddle_count(local_gens, n, A, B):
    """Number of generators from the LOCAL list whose support intersects both A and B."""
    mA = 0
    for q in A:
        mA |= 1 << q
    mB = 0
    for q in B:
        mB |= 1 << q
    c = 0
    for g in local_gens:
        s = supp_mask(g, n)
        if (s & mA) and (s & mB):
            c += 1
    return c


# ------------------------------------------------------------------ earned separation
def qubit_graph(local_gens, n):
    """Interaction graph on qubits: q ~ q' iff they co-occur in some LOCAL Hamiltonian
    generator.  Computed from generator supports alone -- no coordinates imported."""
    adj = [set() for _ in range(n)]
    for g in local_gens:
        s = supp_mask(g, n)
        qs = [q for q in range(n) if (s >> q) & 1]
        for i in range(len(qs)):
            for j in range(i + 1, len(qs)):
                adj[qs[i]].add(qs[j])
                adj[qs[j]].add(qs[i])
    return adj

def graph_dist_regions(adj, A, B):
    """Multi-source BFS distance from qubit set A to qubit set B."""
    B = set(B)
    dist = {a: 0 for a in A}
    dq = deque(A)
    while dq:
        q = dq.popleft()
        if q in B:
            return dist[q]
        for r in adj[q]:
            if r not in dist:
                dist[r] = dist[q] + 1
                dq.append(r)
    return None

def generator_graph_dist(local_gens, n, gi, gj):
    """Distance in the GENERATOR-overlap graph (generators adjacent iff supports
    intersect).  BFS from generator gi to gj (indices into local_gens)."""
    supps = [supp_mask(g, n) for g in local_gens]
    N = len(local_gens)
    adj = [[] for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            if supps[i] & supps[j]:
                adj[i].append(j)
                adj[j].append(i)
    dist = {gi: 0}
    dq = deque([gi])
    while dq:
        i = dq.popleft()
        if i == gj:
            return dist[i]
        for j in adj[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return None


# ------------------------------------------------------------------ exhaustive coset scans
def _np_pop(a):
    return np.bitwise_count(a)

def coset_min_np(gens, rep, base_bits=22, weight="pop", n=None):
    """Exhaustive min weight (and count of minima, and total size) over the coset
    rep * span(gens).  gens/rep are python ints that must fit in uint64 under the chosen
    weight convention:
      weight='pop'  : masks are single-part (pure X or pure Z), weight = popcount.
      weight='xz'   : masks are 2n-bit with 2n <= 64, weight = popcount(x|z), needs n.
    Chunked: base span of min(len,base_bits) generators, remaining generators looped."""
    base = gens[:base_bits]
    rest = gens[base_bits:]
    arr = np.zeros(1, dtype=np.uint64)
    for g in base:
        arr = np.concatenate([arr, arr ^ np.uint64(g)])
    best = None
    cnt = 0
    total = (1 << len(gens))
    for i in range(1 << len(rest)):
        off = rep
        k = i
        j = 0
        while k:
            if k & 1:
                off ^= rest[j]
            k >>= 1
            j += 1
        v = arr ^ np.uint64(off)
        if weight == "pop":
            w = _np_pop(v)
        else:
            maskn = np.uint64((1 << n) - 1)
            w = _np_pop((v & maskn) | (v >> np.uint64(n)))
        m = int(w.min())
        c = int((w == m).sum())
        if best is None or m < best:
            best, cnt = m, c
        elif m == best:
            cnt += c
    return best, cnt, total
