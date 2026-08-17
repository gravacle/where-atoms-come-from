# rule_main4.py -- W-19 RULING, part 4.  Independent exhaustive check of the CORRECTED FLOOR.
# Claim under test (refuter A, out_2_minimality.txt blocks 2c/2d):
#   no connected SIMPLE graph with min degree >= 3 on V <= 7 vertices reaches d = 5,
#   and cubic graphs on V = 8 do, so the simple floor for 4 rule-A plateau points is L = 12.
import itertools, time
from collections import deque
from rule_verify import P, LOG

t0 = time.time()
P("=" * 118)
P("W-19 RULING -- PART 4.  INDEPENDENT EXHAUSTIVE SEARCH FOR THE CORRECTED FLOOR.")
P("=" * 118)

def maxd(V, edges):
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

def connected(V, edges):
    adj = [[] for _ in range(V)]
    for (a, b) in edges: adj[a].append(b); adj[b].append(a)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == V

for V in (5, 6, 7):
    pool = list(itertools.combinations(range(V), 2))
    best = -1; hist = {}; firstL = {}
    n = 0
    for mask in range(1 << len(pool)):
        edges = [pool[i] for i in range(len(pool)) if mask >> i & 1]
        if len(edges) * 2 < 3 * V: continue
        deg = [0] * V
        for (a, b) in edges: deg[a] += 1; deg[b] += 1
        if min(deg) < 3: continue
        if not connected(V, edges): continue
        n += 1
        d = maxd(V, edges)
        hist[d] = hist.get(d, 0) + 1
        if d not in firstL or len(edges) < firstL[d]: firstL[d] = len(edges)
        best = max(best, d)
    P("  V=%d : %d connected simple graphs with min degree >= 3;  MAX d = %d;  d-histogram %s;  "
      "smallest L per d %s" % (V, n, best, dict(sorted(hist.items())), dict(sorted(firstL.items()))))

# V = 8 cubic
V = 8; pool = list(itertools.combinations(range(V), 2))
cnt = 0; win = 0; girths = {}
def girth_of(V, edges):
    g = 10 ** 9
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
        if b in dist: g = min(g, dist[b] + 1)
    return g
for comb in itertools.combinations(range(len(pool)), 12):
    edges = [pool[i] for i in comb]
    deg = [0] * V
    for (a, b) in edges: deg[a] += 1; deg[b] += 1
    if any(x != 3 for x in deg): continue
    cnt += 1
    if not connected(V, edges): continue
    if maxd(V, edges) >= 5:
        win += 1; g = girth_of(V, edges); girths[g] = girths.get(g, 0) + 1
P("  V=8 cubic: %d labelled cubic graphs on 8 vertices; %d of them reach d >= 5 at L = 12; "
  "girth histogram of the winners %s" % (cnt, win, dict(sorted(girths.items()))))
P("  min degree 3 forces 2L >= 3V, so V >= 8 forces L >= 12: the SIMPLE floor for four rule-A")
P("  plateau points is therefore exactly L = 12, not 21.")
P("")
P("elapsed %.1f s" % (time.time() - t0))
open("OUT_rule_main4.txt", "w").write("\n".join(LOG) + "\n")
