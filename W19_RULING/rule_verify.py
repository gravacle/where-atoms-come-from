# rule_verify.py -- W-19 RULING, INDEPENDENT VERIFICATION.
#
# Written from the mathematics, not from either lane's source.  No import of zn_gauge.py, lib_b.py,
# rlib.py, carriers.py or any other lane file.  numpy only (no scipy on this machine).
#
# CONVENTIONS (fixed here; identical in content to the brief's):
#   Z_2 on each link.  Z|j> = (-1)^j |j>,  X|j> = |j+1 mod 2>.
#   G_v = prod_{l incident to v} X_l.   PHYSICAL = +1 eigenspace of every G_v.
#   W_p = prod_{l in p} Z_l  for p in the GF(2) cycle space.
#   H = -(1/g2) sum_p W_p - g2 sum_l X_l      (lane A's H is exactly 2x this; same eigenvectors)
#   PLAQUETTES = a minimum-weight cycle basis obtained by MATROID GREEDY over the WHOLE cycle space
#   (enumerate all 2^C elements, sort by weight, take GF(2)-independent ones).  Independent route
#   from lane A's Horton enumeration.
#
# INDEXING.  A state on L links is a flat vector of length 2^L reshaped to [2]*L.
#   AXIS(l) = L-1-l, i.e. link 0 is the MOST significant bit.  Everything below is written in AXES.
#
# THE THREE CHANNELS, measured on the SAME state with the SAME fragments and the SAME delta.
# Exactly one thing moves between them: the ALGEBRA assigned to the system link.
#   EXT  I(B(H_l) : B(H_F))       full 2x2 matrix algebra on the system link, extended Hilbert
#                                 space.  This is what lane A and lane B READING 1 both compute.
#   CHI  I(alg{X_l} : B(H_F))     system algebra = the only non-trivial GAUGE-INVARIANT subalgebra
#                                 of a single link (the electric flux).  Holevo quantity.
#   CL   I(X_l : X_F)             both sides gauge-invariant and commuting: a classical mutual
#                                 information between electric fluxes.
import numpy as np
from collections import deque

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)

# ------------------------------------------------------------------ graph / GF(2)

def fundamental_cycles(V, edges):
    L = len(edges)
    adj = [[] for _ in range(V)]
    for i, (a, b) in enumerate(edges):
        adj[a].append((b, i)); adj[b].append((a, i))
    par = {0: (None, None)}; tree = set(); seen = {0}; dq = deque([0])
    while dq:
        x = dq.popleft()
        for (y, i) in adj[x]:
            if y not in seen:
                seen.add(y); tree.add(i); par[y] = (x, i); dq.append(y)
    assert len(seen) == V, "disconnected"
    ch = {v: [] for v in range(V)}
    for v in range(V):
        if v: ch[par[v][0]].append(v)
    path = {0: 0}; dq = deque([0])
    while dq:
        x = dq.popleft()
        for y in ch[x]:
            path[y] = path[x] ^ (1 << par[y][1]); dq.append(y)
    return [(1 << i) ^ path[edges[i][0]] ^ path[edges[i][1]] for i in range(L) if i not in tree]

def gf2_reduce(basis, v):
    for b in basis:
        v = min(v, v ^ b)
    return v

def min_weight_cycle_basis(V, edges):
    fund = fundamental_cycles(V, edges); C = len(fund)
    allc = []
    for mask in range(1, 1 << C):
        v = 0
        for i in range(C):
            if mask >> i & 1: v ^= fund[i]
        allc.append(v)
    allc.sort(key=lambda v: bin(v).count("1"))
    red = []; chosen = []
    for v in allc:
        r = gf2_reduce(red, v)
        if r:
            red = sorted(red + [r], reverse=True); chosen.append(v)
            if len(chosen) == C: break
    assert len(chosen) == C
    return chosen

def bfs_dist(V, edges, src, skip):
    adj = [[] for _ in range(V)]
    for i, (a, b) in enumerate(edges):
        if i == skip: continue
        adj[a].append(b); adj[b].append(a)
    dist = {src: 0}; dq = deque([src])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if y not in dist: dist[y] = dist[x] + 1; dq.append(y)
    return dist

def rule_A_fragments(V, edges, l):
    """nested fragments grown by BFS distance from tail(l) inside G-l.  Returns (fragments, d)."""
    u, v = edges[l]
    dist = bfs_dist(V, edges, u, l); d = dist[v]
    lev = {i: max(dist[a], dist[b]) for i, (a, b) in enumerate(edges) if i != l}
    out = []
    for k in range(1, d + 1):
        f = sorted([i for i in lev if lev[i] <= k])
        if not out or f != out[-1]: out.append(f)
    return out, d

def girth(V, edges):
    best = 10 ** 9
    for l, (a, b) in enumerate(edges):
        dd = bfs_dist(V, edges, a, l)
        if b in dd: best = min(best, dd[b] + 1)
    return best

# ------------------------------------------------------------------ physical sector

class Carrier:
    def __init__(self, name, V, edges):
        self.name, self.V, self.edges = name, V, list(edges)
        self.L = len(edges); self.C = self.L - V + 1
        self.plaq = min_weight_cycle_basis(V, edges)
        self.cyc = self.plaq
        self.tvec = []
        for l in range(self.L):
            t = 0
            for i, c in enumerate(self.cyc):
                if c >> l & 1: t |= (1 << i)
            self.tvec.append(t)
        self.dimP = 1 << self.C
        deg = [0] * V
        for (a, b) in edges: deg[a] += 1; deg[b] += 1
        self.mindeg = min(deg)
        self.girth = girth(V, edges)

    def H_phys(self, g2):
        D = self.dimP; H = np.zeros((D, D)); m = np.arange(D)
        for i in range(self.C):
            H[m, m] += -(1.0 / g2) * (1.0 - 2.0 * ((m >> i) & 1))
        for l in range(self.L):
            H[m ^ self.tvec[l], m] += -g2
        return H

    def orbit_label(self):
        z = np.arange(1 << self.L, dtype=np.uint64)
        m = np.zeros(1 << self.L, dtype=np.int64)
        for i, c in enumerate(self.cyc):
            m |= ((np.bitwise_count(z & np.uint64(c)) & 1).astype(np.int64) << i)
        return m

    def lift(self, psi):
        return psi[self.orbit_label()] / (2.0 ** ((self.V - 1) / 2.0))

    def ground(self, g2):
        w, U = np.linalg.eigh(self.H_phys(g2))
        return w, U[:, 0]

# ------------------------------------------------------------------ entropies, in AXES

def ent_rho(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-13]
    return float(-(ev * np.log2(ev)).sum())

def rdm_ax(vec, n, axesA):
    axesA = sorted(axesA); axesB = [a for a in range(n) if a not in axesA]
    T = np.transpose(vec.reshape([2] * n), axesA + axesB)
    T = T.reshape(1 << len(axesA), 1 << len(axesB))
    return T @ T.conj().T

def S_ax(vec, n, axesA):
    """entropy of the axis set for a GLOBAL PURE vector; computed on the smaller side."""
    axesA = sorted(axesA); axesB = [a for a in range(n) if a not in axesA]
    B = axesA if len(axesA) <= len(axesB) else axesB
    return ent_rho(rdm_ax(vec, n, B))

def AX(L, links):
    return [L - 1 - l for l in links]

def hadamard_all(vec, L):
    v = vec.reshape([2] * L).astype(complex)
    for a in range(L):
        v = np.moveaxis(v, a, 0)
        p, q = v[0].copy(), v[1].copy()
        v[0] = (p + q) / np.sqrt(2.0); v[1] = (p - q) / np.sqrt(2.0)
        v = np.moveaxis(v, 0, a)
    return v.reshape(-1)

# ------------------------------------------------------------------ channels

def channel_EXT(vec, L, l, F):
    aS = AX(L, [l]); aF = AX(L, F)
    return S_ax(vec, L, aS) + S_ax(vec, L, aF) - S_ax(vec, L, aS + aF)

def channel_CHI(vecX, L, l, F):
    """I(alg{X_l} : B(H_F)) = Holevo chi of the cq state.  vecX = state in the X product basis."""
    al = L - 1 - l
    T = np.moveaxis(vecX.reshape([2] * L), al, 0)
    rest = [a for a in range(L) if a != al]                    # original axes, in order
    newidx = {a: j for j, a in enumerate(rest)}
    aF = [newidx[L - 1 - f] for f in F]
    ps, Ss = [], []
    for s in (0, 1):
        c = T[s].reshape(-1); p = float(np.vdot(c, c).real)
        ps.append(p)
        Ss.append(0.0 if p < 1e-14 else S_ax(c / np.sqrt(p), L - 1, aF))
    return S_ax(vecX, L, AX(L, F)) - sum(p * s for p, s in zip(ps, Ss))

def channel_CL(probX, L, l, F):
    keepax = sorted(AX(L, [l] + list(F)))
    T = probX.reshape([2] * L)
    dropax = tuple(a for a in range(L) if a not in keepax)
    pj = T.sum(axis=dropax) if dropax else T
    pos = keepax.index(L - 1 - l)
    pj = np.moveaxis(pj, pos, 0).reshape(2, -1)
    def h(p):
        p = np.asarray(p); p = p[p > 1e-15]
        return float(-(p * np.log2(p)).sum())
    return h(pj.sum(axis=1)) + h(pj.sum(axis=0)) - h(pj.reshape(-1))

# ------------------------------------------------------------------ full-space controls

def gauss_masks(car):
    out = []
    for v in range(car.V):
        m = 0
        for i, (a, b) in enumerate(car.edges):
            if a == v or b == v: m ^= (1 << i)
        out.append(m)
    return out

def gauss_residual(car, vec):
    idx = np.arange(len(vec)); worst = 0.0
    for m in gauss_masks(car):
        worst = max(worst, float(np.abs(vec[idx ^ m] - vec).max()))
    return worst

def H_full(car, g2):
    D = 1 << car.L; z = np.arange(D); H = np.zeros((D, D))
    for p in car.plaq:
        H[z, z] += -(1.0 / g2) * (1.0 - 2.0 * (np.bitwise_count(z & np.int64(p)) & 1))
    for l in range(car.L):
        H[z ^ (1 << l), z] += -g2
    return H

def haar_physical(car, seed):
    rng = np.random.default_rng(seed)
    psi = rng.standard_normal(car.dimP) + 1j * rng.standard_normal(car.dimP)
    return psi / np.linalg.norm(psi)

# ------------------------------------------------------------------ measurement

def measure(car, vec, l, tag):
    L = car.L
    frags, d = rule_A_fragments(car.V, car.edges, l)
    vecX = hadamard_all(vec, L); probX = np.abs(vecX) ** 2
    HS = S_ax(vec, L, AX(L, [l]))
    T = probX.reshape([2] * L)
    pm = T.sum(axis=tuple(a for a in range(L) if a != L - 1 - l))
    HX = float(-(pm[pm > 1e-15] * np.log2(pm[pm > 1e-15])).sum())
    rows = []
    for F in frags:
        rows.append((len(F),
                     channel_EXT(vec, L, l, F),
                     channel_CHI(vecX, L, l, F),
                     channel_CL(probX, L, l, F)))
    nz = lambda x, h: (x / h if h > 1e-12 else 0.0)
    P("  %-22s d=%d  H_EXT(S)=%.9f  H_elec(S)=%.9f" % (tag, d, HS, HX))
    P("      |F|  :" + "".join("%11d" % r[0] for r in rows))
    P("      EXT  :" + "".join("%11.6f" % nz(r[1], HS) for r in rows))
    P("      CHI  :" + "".join("%11.6f" % nz(r[2], HX) for r in rows))
    P("      CL   :" + "".join("%11.6f" % nz(r[3], HX) for r in rows))
    pts = lambda vs: sum(1 for v in vs if abs(v - 1.0) <= 0.10)
    nE = pts([nz(r[1], HS) for r in rows])
    nC = pts([nz(r[2], HX) for r in rows])
    nK = pts([nz(r[3], HX) for r in rows])
    P("      plateau points (delta=0.10):  EXT=%d  CHI=%d  CL=%d" % (nE, nC, nK))
    return dict(d=d, HS=HS, HX=HX, rows=rows, nE=nE, nC=nC, nK=nK)

# ------------------------------------------------------------------ carriers

def theta(L):     return 2, [(0, 1)] * L

def mg_chain(d):
    m = [1] * d; m[0] = 2; m[-1] = 2
    for i in range(1, d - 1):
        if m[i - 1] + m[i] < 3: m[i] = 3 - m[i - 1]
    if d >= 2 and m[-2] + m[-1] < 3: m[-2] = 3 - m[-1]
    edges = [(0, d)]
    for i in range(d):
        for _ in range(m[i]): edges.append((i, i + 1))
    return d + 1, edges

def tri_chain12():
    return 8, [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
               (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]

def heawood():
    V = 14; e = [(i, (i + 1) % 14) for i in range(14)]
    lcf = [5, -5] * 7; seen = set()
    for i in range(14):
        k = tuple(sorted((i, (i + lcf[i]) % 14)))
        if k not in seen: seen.add(k); e.append(k)
    return V, e

def petersen():
    e = [(i, (i + 1) % 5) for i in range(5)]
    e += [(i, i + 5) for i in range(5)]
    e += [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return 10, e
