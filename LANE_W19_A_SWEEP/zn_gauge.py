# zn_gauge.py -- LANE W19-A. Z_N pure lattice gauge theory on a finite graph.
#
# EXACT construction. No penalty terms, no projection: the Hamiltonian is built DIRECTLY on the
# gauge-invariant (physical) sector by gauge-fixing to a spanning tree. Every state produced here
# is physical by construction. (Attempt one's T1 used a 1e6 penalty; this is the exact version.)
#
# CONVENTIONS (fixed here once, quoted in every output file):
#   * Graph G = (V vertices, L directed links). Multi-edges allowed. Link e = (tail_e, head_e).
#   * Each link carries Z_N: Z|j> = w^j |j>, X|j> = |j+1 mod N>, w = exp(2 pi i / N).
#   * Gauge operator at vertex v:  G_v = prod_e X_e^{ d[e,v] },  d[e,v] = delta(v,head_e) - delta(v,tail_e).
#     PHYSICAL = the +1 eigenspace of every G_v.  dim_phys = N^(L-V+1) = N^C, C = cyclomatic number.
#   * Wilson loop on a cycle p (signed edge vector, boundary zero): W_p = prod_e Z_e^{p_e}.
#   * Hamiltonian, brief's form:
#         H = -(1/g^2) sum_p (W_p + W_p^dag) - g^2 sum_e (X_e + X_e^dag)
#     implemented as  H = -(mag/2) sum_p (W_p + W_p^dag) - (elec/2) sum_e (X_e + X_e^dag)
#     with the STANDARD choice mag = 2/g^2, elec = 2*g^2.  g^2 IS THE DIMENSIONLESS SLOT.
#   * PLAQUETTE SET P = a MINIMUM-WEIGHT CYCLE BASIS (Horton-style: enumerate all simple cycles up
#     to a length cutoff, sort by length, greedily take GF(2)-independent ones until rank = C).
#     On a planar lattice this returns the faces; on the theta graph it returns two 2-cycles, which
#     is exactly what the sealed T1 run used.
#
# GAUGE FIXING / ORBIT BASIS.
#   Pick a spanning tree T. Every gauge orbit in the Z-basis has a UNIQUE representative with
#   z_e = 0 on all tree links. The orbit label is m = Gamma z mod N, where Gamma is the C x L matrix
#   whose rows are the FUNDAMENTAL CYCLES of the C chords. Consequences, all exact:
#       X_e   translates m by  Gamma[:,e]   (mod N)
#       W_p   is diagonal with phase  w^{a_p . m},  a_p[c] = p[chord_c]
#       Psi(z) = psi(Gamma z mod N) / N^((V-1)/2)     (full-space amplitude, normalised)
#   A bridge link has Gamma[:,e] = 0: X_e acts as the identity on the physical sector, which is
#   correct (it is a product of Gauss operators).

import numpy as np, itertools
from collections import deque

# ---------------------------------------------------------------- graph utilities

def build_adj(V, edges):
    adj = [[] for _ in range(V)]
    for i, (a, b) in enumerate(edges):
        adj[a].append((b, i)); adj[b].append((a, i))
    return adj

def spanning_tree(V, edges):
    adj = build_adj(V, edges)
    seen = [False]*V; seen[0] = True; tree = set(); par = {0: (None, None)}
    dq = deque([0])
    while dq:
        x = dq.popleft()
        for (y, i) in adj[x]:
            if not seen[y]:
                seen[y] = True; tree.add(i); par[y] = (x, i); dq.append(y)
    assert all(seen), "graph not connected"
    return tree, par, adj

def tree_path_vectors(V, edges, tree, par):
    """tp[w] in Z^L : signed sum of tree links along the root->w tree path."""
    L = len(edges); tp = [None]*V; tp[0] = np.zeros(L, dtype=np.int64)
    order = sorted(range(V), key=lambda w: 0)          # BFS order via par chain
    todo = deque([0]); child = {w: [] for w in range(V)}
    for w in range(V):
        if w != 0: child[par[w][0]].append(w)
    while todo:
        x = todo.popleft()
        for y in child[x]:
            i = par[y][1]; a, b = edges[i]
            v = tp[x].copy()
            v[i] += 1 if (a == x and b == y) else -1   # traverse x -> y
            tp[y] = v; todo.append(y)
    return tp

def fundamental_cycles(V, edges):
    """returns Gamma (C x L int) and the chord index list."""
    L = len(edges)
    tree, par, _ = spanning_tree(V, edges)
    tp = tree_path_vectors(V, edges, tree, par)
    chords = [i for i in range(L) if i not in tree]
    G = np.zeros((len(chords), L), dtype=np.int64)
    for r, i in enumerate(chords):
        a, b = edges[i]
        G[r] = tp[a] - tp[b]
        G[r, i] += 1                                    # e_i + tp[a] - tp[b] has zero boundary
    # verify boundary zero
    for r in range(len(chords)):
        bd = np.zeros(V, dtype=np.int64)
        for e, (a, b) in enumerate(edges):
            bd[b] += G[r, e]; bd[a] -= G[r, e]
        assert not bd.any(), "fundamental cycle has nonzero boundary"
    return G, chords

def simple_cycles(V, edges, maxlen):
    """all simple cycles (as signed edge vectors) of length <= maxlen; multigraph-safe."""
    adj = build_adj(V, edges); out = {}
    def dfs(start, cur, vec, used_e, used_v, length):
        if length >= maxlen: return
        for (y, i) in adj[cur]:
            if i in used_e: continue
            a, b = edges[i]; s = 1 if (a == cur and b == y) else -1
            if y == start and length + 1 >= 2:
                key = frozenset(used_e | {i})
                if key not in out:
                    v = vec.copy(); v[i] += s; out[key] = v
                continue
            if y in used_v or y < start: continue
            v = vec.copy(); v[i] += s
            dfs(start, y, v, used_e | {i}, used_v | {y}, length+1)
    for s in range(V):
        dfs(s, s, np.zeros(len(edges), dtype=np.int64), frozenset(), frozenset([s]), 0)
    return sorted(out.values(), key=lambda v: int(np.abs(v).sum()))

def gf2_rank_add(basis, vec):
    """basis: list of int bitmasks in reduced form. returns (added?, new basis)"""
    v = vec
    for b in basis:
        v = min(v, v ^ b)
    if v == 0: return False, basis
    nb = sorted(basis + [v], reverse=True)
    return True, nb

def min_cycle_basis(V, edges, C, maxlen=None):
    """Horton-style minimum weight cycle basis. Raises the cutoff until rank C is reached."""
    L = len(edges)
    cut = maxlen or 4
    while True:
        cyc = simple_cycles(V, edges, cut)
        basis = []; chosen = []
        for v in cyc:
            mask = 0
            for e in range(L):
                if v[e] % 2: mask |= (1 << e)
            ok, basis = gf2_rank_add(basis, mask)
            if ok:
                chosen.append(v)
                if len(chosen) == C: return chosen
        cut += 1
        if cut > L: raise RuntimeError("no cycle basis found")

def girth_through(V, edges, l):
    """length of the shortest cycle through link l, and d = dist_{G-l}(tail,head)."""
    a, b = edges[l]
    adj = build_adj(V, edges)
    dist = {a: 0}; dq = deque([a])
    while dq:
        x = dq.popleft()
        for (y, i) in adj[x]:
            if i == l: continue
            if y not in dist:
                dist[y] = dist[x] + 1; dq.append(y)
    d = dist.get(b, None)
    return (None, None) if d is None else (d + 1, d)

def bfs_levels(V, edges, l):
    """BFS levels from tail(l) in G - l. returns dist dict."""
    a, b = edges[l]; adj = build_adj(V, edges)
    dist = {a: 0}; dq = deque([a])
    while dq:
        x = dq.popleft()
        for (y, i) in adj[x]:
            if i == l: continue
            if y not in dist:
                dist[y] = dist[x] + 1; dq.append(y)
    return dist

# ---------------------------------------------------------------- physical-sector Hamiltonian

class ZNGauge:
    def __init__(self, name, V, edges, N, maxlen=None):
        self.name, self.V, self.edges, self.N = name, V, list(edges), N
        self.L = len(edges); self.C = self.L - V + 1
        self.Gam, self.chords = fundamental_cycles(V, edges)
        assert self.Gam.shape[0] == self.C
        self.plaq = min_cycle_basis(V, edges, self.C, maxlen)
        self.dimP = N ** self.C
        # translation vector per link, in Z_N^C
        self.tvec = [tuple(int(x) % N for x in self.Gam[:, e]) for e in range(self.L)]
        # plaquette coefficient vector a_p in Z_N^C
        self.avec = [tuple(int(p[self.chords[c]]) % N for c in range(self.C)) for p in self.plaq]

    def _digits(self):
        N, C = self.N, self.C
        idx = np.arange(self.dimP)
        return np.stack([(idx // (N**c)) % N for c in range(C)], axis=1)   # m_c = digit c

    def hamiltonian(self, mag, elec):
        N, C, D = self.N, self.C, self.dimP
        m = self._digits()
        H = np.zeros((D, D))
        # magnetic (diagonal)
        diag = np.zeros(D)
        for a in self.avec:
            ph = (m @ np.array(a, dtype=np.int64)) % N
            diag += -(mag/2.0) * 2.0*np.cos(2*np.pi*ph/N)
        np.fill_diagonal(H, diag)
        # electric (translations)
        pw = np.array([N**c for c in range(C)], dtype=np.int64)
        for t in self.tvec:
            t = np.array(t, dtype=np.int64)
            if not t.any():
                H[np.arange(D), np.arange(D)] += -(elec/2.0)*2.0     # bridge: X = identity
                continue
            for s in (+1, -1):
                mm = (m + s*t) % N
                j = mm @ pw
                H[j, np.arange(D)] += -(elec/2.0)
        return H

    def ground(self, mag, elec):
        H = self.hamiltonian(mag, elec)
        w, v = np.linalg.eigh(H)
        psi = v[:, 0].real.copy()
        gap = float(w[1] - w[0])
        return psi/np.linalg.norm(psi), float(w[0]), gap

    def full_vector(self, psi):
        """Psi(z) over the FULL N^L Z-basis, normalised."""
        N, L, C = self.N, self.L, self.C
        tot = N**L
        if tot > 4_500_000:
            raise MemoryError(f"full vector {N}^{L} = {tot} exceeds the declared ceiling 4.5e6")
        acc = np.zeros((C, tot), dtype=np.int8)
        idx = np.arange(tot, dtype=np.int64)
        for e in range(L):
            d = ((idx // (N**(L-1-e))) % N).astype(np.int8)     # link 0 = most significant
            for c in range(C):
                g = int(self.Gam[c, e]) % N
                if g: acc[c] = (acc[c] + g*d) % N
        code = np.zeros(tot, dtype=np.int64)
        for c in range(C):
            code += (N**c) * acc[c].astype(np.int64)
        Psi = psi[code] / np.sqrt(float(N)**(self.V-1))
        return Psi

# ---------------------------------------------------------------- entropies

def _entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-13]
    return float(-(w*np.log2(w)).sum())

def rdm(Psi, L, N, A):
    A = sorted(A); B = [i for i in range(L) if i not in A]
    T = Psi.reshape([N]*L).transpose(A+B).reshape(N**len(A), N**len(B))
    return T @ T.T

def S_of(Psi, L, N, A):
    """von Neumann entropy in bits of the reduced state on link set A (uses the smaller side)."""
    A = sorted(A)
    if len(A) == 0: return 0.0
    if len(A) == L: return 0.0
    B = [i for i in range(L) if i not in A]
    use = A if len(A) <= len(B) else B
    return _entropy(rdm(Psi, L, N, use))

def mutual_information(Psi, L, N, A, B):
    A = sorted(A); B = sorted(B)
    return S_of(Psi, L, N, A) + S_of(Psi, L, N, B) - S_of(Psi, L, N, sorted(set(A)|set(B)))

# ---------------------------------------------------------------- fragment rules

def level_cuts(V, edges, l):
    """RULE C. The d edge-disjoint cuts through l: C_i = links from level i-1 to level >= i,
       in G - l, BFS from tail(l). Each satisfies X_l = X(C_i)^{-1} EXACTLY (Gauss law)."""
    dist = bfs_levels(V, edges, l); b = edges[l][1]; d = dist[b]
    cuts = []
    for i in range(1, d+1):
        Ci = [e for e, (x, y) in enumerate(edges)
              if e != l and ((dist.get(x, 10**9) <= i-1) != (dist.get(y, 10**9) <= i-1))]
        cuts.append(sorted(Ci))
    return cuts, d

def has_uv_path(V, edges, l, F):
    """does fragment F contain a tail(l)->head(l) path in G-l ?  Equivalently: does S u F contain a
       CYCLE through l, i.e. is a Wilson loop through l available inside S u F ?  That is exactly
       the condition under which the conjugate (magnetic) bit becomes readable and I jumps to 2H(S)."""
    a, b = edges[l]; adj = [[] for _ in range(V)]
    for i in F:
        x, y = edges[i]; adj[x].append(y); adj[y].append(x)
    seen = {a}; dq = deque([a])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if y not in seen:
                if y == b: return True
                seen.add(y); dq.append(y)
    return b in seen

def nested_fragments(V, edges, l):
    """RULE A. F_k = links of G-l with at least one endpoint at BFS distance <= k-1 from tail(l).
       F_1 = star(tail)-l ; F_k contains no tail->head path for k <= d-1 ; F_d does."""
    dist = bfs_levels(V, edges, l); b = edges[l][1]; d = dist[b]
    frs = []
    for k in range(1, d+1):
        Fk = [e for e, (x, y) in enumerate(edges)
              if e != l and (dist.get(x, 10**9) <= k-1 or dist.get(y, 10**9) <= k-1)]
        frs.append(sorted(Fk))
    return frs, d
