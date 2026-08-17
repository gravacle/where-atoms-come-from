# refute_carriers.py -- COUNTEREXAMPLE CARRIERS for the W19-A minimality claim.
#
# Lane A, out_04_threshold.txt block 4b, states the minimality argument verbatim:
#   "if every vertex has degree >= 3 and the shortest cycle through the system link l has length >= 6
#    (needed for d >= 5 ...), then L >= 21 ... 1+2+4 from u, 1+2+4 from v (disjoint) = 14 vertices"
# The step |L2| >= 4 says the two level-1 neighbours of u have DISJOINT further neighbourhoods.
# That needs NO 4-CYCLE THROUGH u -- i.e. GLOBAL girth >= 5.  The hypothesis actually assumed is only
# "girth THROUGH l >= 6", which constrains cycles that USE l and says nothing about cycles avoiding l.
# Every carrier below satisfies the STATED hypothesis (min degree >= 3, girth_through(l) = d+1 >= 6)
# and violates the conclusion (L >= 21).

# ---------------------------------------------------------------- named exhibits

def tri_chain12():
    """SIMPLE, cubic, V=8, L=12, girth 3, d = dist_{G-l}(u,v) = 5 -> girth_through(link 0) = 6.
       Two triangles joined by a bridge; the system link closes the long way round.
       BFS levels from u in G-l: {0} {1,2} {3} {4} {5,6} {7} = 1,2,1,1,2,1 -> 8 vertices."""
    return 8, [(0, 7),                          # l, index 0
               (0, 1), (0, 2), (1, 2),          # triangle at u
               (1, 3), (2, 3),                  # -> c
               (3, 4),                          # bridge
               (4, 5), (4, 6),                  # -> e,f
               (5, 6), (5, 7), (6, 7)]          # triangle at v

def tri_chain11():
    """SIMPLE, min degree 3 (one vertex of degree 4), V=7, L=11, d=4 -> 3 plateau points.
       Refutes the lane's secondary floor 'V >= 10, L >= 15, attained by Petersen'."""
    return 7, [(0, 6),
               (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
               (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]

def mg_chain(d):
    """MULTIGRAPH floor family (lane A's own framework allows multi-edges -- its theta carrier is
       three parallel links).  V = d+1, all degrees >= 3, L = ceil(3(d+1)/2), dist_{G-l}(u,v) = d.
       Multiplicities alternate 2,1,2,1,... with a 2 forced at each end."""
    m = [1] * d
    m[0] = 2; m[-1] = 2
    for i in range(1, d - 1):                     # enforce m[i-1]+m[i] >= 3 greedily, cheapest first
        if m[i - 1] + m[i] < 3: m[i] = 3 - m[i - 1]
    if d >= 2 and m[-2] + m[-1] < 3: m[-2] = 3 - m[-1]
    edges = [(0, d)]
    for i in range(d):
        for _ in range(m[i]): edges.append((i, i + 1))
    return d + 1, edges

def dbl_chain9():
    """mg_chain(5): V=6, L=9, all degrees exactly 3, d=5.  THE ABSOLUTE FLOOR for 4 plateau points
       under 'min degree >= 3': d=5 forces V >= 6, min degree 3 forces L >= ceil(3*6/2) = 9."""
    return mg_chain(5)

# ---------------------------------------------------------------- exhaustive search helpers

def min_degree(V, edges):
    deg = [0] * V
    for (a, b) in edges: deg[a] += 1; deg[b] += 1
    return min(deg) if V else 0

def connected(V, edges):
    adj = [[] for _ in range(V)]
    for (a, b) in edges: adj[a].append(b); adj[b].append(a)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == V

def max_d(V, edges):
    """max over links l of dist_{G-l}(tail,head); None if some removal disconnects the pair."""
    from collections import deque
    best = -1
    for l, (a, b) in enumerate(edges):
        adj = [[] for _ in range(V)]
        for i, (x, y) in enumerate(edges):
            if i == l: continue
            adj[x].append(y); adj[y].append(x)
        dist = {a: 0}; dq = deque([a])
        while dq:
            x = dq.popleft()
            for y in adj[x]:
                if y not in dist: dist[y] = dist[x] + 1; dq.append(y)
        if b in dist: best = max(best, dist[b])
    return best

def all_simple_graphs(V, maxL, mindeg=3):
    """every connected labelled simple graph on V vertices with <= maxL edges and min degree >= mindeg."""
    import itertools
    pairs = list(itertools.combinations(range(V), 2))
    P = len(pairs)
    for mask in range(1 << P):
        if bin(mask).count("1") > maxL: continue
        E = [pairs[i] for i in range(P) if mask >> i & 1]
        if len(E) < (V * mindeg + 1) // 2: continue
        if min_degree(V, E) < mindeg: continue
        if not connected(V, E): continue
        yield E

def cubic_graphs(V):
    """all labelled cubic simple graphs on V vertices, by recursive degree saturation."""
    res = []; adj = [set() for _ in range(V)]; E = []
    def rec(v):
        if v == V:
            if all(len(a) == 3 for a in adj): res.append(list(E))
            return
        need = 3 - len(adj[v])
        if need < 0: return
        cand = [u for u in range(v + 1, V) if len(adj[u]) < 3 and u not in adj[v]]
        if need > len(cand): return
        import itertools
        for combo in itertools.combinations(cand, need):
            for u in combo: adj[v].add(u); adj[u].add(v); E.append((v, u))
            rec(v + 1)
            for u in combo: adj[v].discard(u); adj[u].discard(v); E.pop()
    rec(0)
    return res
