"""ADVERSARIAL VERIFY for LANE_T43_B_WORLD -- independent instruments, default refuted.

Independent of the lane's code: max-flow here is plain BFS augmenting-path Ford-Fulkerson
(the lane used Dinic); constructions rebuilt from the model statement, not copied.

ATTACKS:
  A1  RIG TEST -- is the certifiability definition rigged toward a verdict?  Recompute
      CERT/epoch under the lane's definition (grain-disjoint) AND two alternative honest
      operationalisations: (a) EDGE-disjoint chains (relay grains shareable, channel
      multiplicity counted), (b) CHANNEL-disjoint only (weakest independence: chains just
      use distinct channels and distinct targets).  If any honest reading flips the
      per-epoch degree from 2 to 3, the lane's headline is rigged.  Also: expose how much
      of the grain-disjoint per-epoch bound is DEFINITIONAL (capacity == #channel-bearing
      grains by a two-line argument) vs geometric.
  A2  VOLUME-CONTROL SENSITIVITY -- recompute the fully-channelled control with the
      independent instrument; add a control the lane did NOT run: a single-channel block
      (the instrument must also see the OTHER extreme, capacity 1, E_min = n^3).
  A3  SPLIT HONESTY -- independently recompute N_opt(2) on the uncapped range n = 8..12
      (own time-expanded construction, own flow) and the E_min table; confirm both arms
      of the split (per-epoch boundary, cumulative volume) are real computed facts.
  A5  TABLE CROSS-CHECK -- finding's table values vs independent recomputation.

All checks print computed booleans; exit 0 iff all pass.
"""
import sys
from collections import deque

def say(*a):
    print(*a); sys.stdout.flush()

ALL = True
def chk(label, ok, detail=""):
    global ALL
    ALL &= bool(ok)
    say(f"  [{'OK' if ok else 'REFUTE'}] {label}" + (f"  {detail}" if detail else ""))

# ---------------- independent max-flow: BFS augmenting paths (Edmonds-Karp style) -------
class Flow:
    def __init__(self):
        self.adj = {}
    def add(self, u, v, c):
        self.adj.setdefault(u, {}); self.adj.setdefault(v, {})
        self.adj[u][v] = self.adj[u].get(v, 0) + c
        self.adj[v].setdefault(u, 0)
    def maxflow(self, s, t):
        total = 0
        while True:
            par = {s: None}
            q = deque([s])
            while q and t not in par:
                u = q.popleft()
                for v, c in self.adj[u].items():
                    if c > 0 and v not in par:
                        par[v] = u
                        q.append(v)
            if t not in par:
                return total
            # bottleneck
            b = float('inf'); v = t
            while par[v] is not None:
                u = par[v]; b = min(b, self.adj[u][v]); v = u
            v = t
            while par[v] is not None:
                u = par[v]
                self.adj[u][v] -= b
                self.adj[v][u] += b
                v = u
            total += b

NB = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

def cells_of(n):
    return [(x,y,z) for x in range(n) for y in range(n) for z in range(n)]

def channel_count(g, S):
    return sum(1 for d in NB if (g[0]+d[0], g[1]+d[1], g[2]+d[2]) not in S)

def grain_disjoint_cap(n, entry=None):
    """Lane definition, independent instrument: vertex-capacity-1 disjoint chains."""
    cells = cells_of(n); S = set(cells)
    if entry is None:
        entry = [g for g in cells if channel_count(g, S)]
    f = Flow()
    for g in cells:
        f.add(('in', g), ('out', g), 1)
        for d in NB:
            h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
            if h in S:
                f.add(('out', g), ('in', h), 1)
        f.add(('out', g), 'T', 1)
    for g in entry:
        f.add('S', ('in', g), 1)
    return f.maxflow('S', 'T')

def edge_disjoint_cap(n):
    """Alternative honest reading (a): chains share relay grains but not edges; each of
    the 6n^2 channels may carry its own chain (multiplicity counted at corner/edge grains).
    Certificates still target distinct grains (unit sink per grain)."""
    cells = cells_of(n); S = set(cells)
    f = Flow()
    seen = set()
    for g in cells:
        for d in NB:
            h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
            if h in S and (h, g) not in seen:
                seen.add((g, h))
                f.add(('v', g), ('v', h), 1)
                f.add(('v', h), ('v', g), 1)
        f.add(('v', g), ('t', g), 1)   # certify target g once
        f.add(('t', g), 'T', 1)
        c = channel_count(g, S)
        if c:
            f.add('S', ('v', g), c)    # channel multiplicity
    return f.maxflow('S', 'T')

def channel_disjoint_cap(n):
    """Alternative honest reading (b), the WEAKEST independence: chains only need
    distinct channels and distinct targets; grains and edges freely shared.
    Every grain is reachable from every channel, so capacity = min(6n^2, n^3)."""
    return min(6*n*n, n**3)

def n_opt_indep(n, E, entry=None):
    """Independent time-expanded optimum, rebuilt from the model statement."""
    cells = cells_of(n); S = set(cells)
    if entry is None:
        entry = [g for g in cells if channel_count(g, S)]
    f = Flow()
    for e in range(E):
        for g in cells:
            f.add(('in', e, g), ('out', e, g), 1)
            for d in NB:
                h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
                if h in S:
                    f.add(('out', e, g), ('in', e, h), 1)
            f.add(('out', e, g), ('cert', g), 1)
        for g in entry:
            f.add('S', ('in', e, g), 1)
    for g in cells:
        f.add(('cert', g), 'T', 1)
    return f.maxflow('S', 'T')

def diff_order(seq):
    s = list(seq)
    for k in range(len(s) - 1):
        if len(set(s)) == 1:
            return k
        s = [b - a for a, b in zip(s, s[1:])]
    return None

say("=" * 96)
say("A1  RIG TEST -- the definition vs two alternative honest operationalisations")
say("=" * 96)

NS = list(range(2, 10))
gd = [grain_disjoint_cap(n) for n in NS]
say("    n:               " + " ".join(f"{n:>5}" for n in NS))
say("    grain-disjoint:  " + " ".join(f"{v:>5}" for v in gd))
chk("independent instrument reproduces CERT/epoch == 6n^2-12n+8 (n=2..9)",
    all(v == 6*n*n - 12*n + 8 for v, n in zip(gd, NS)),
    f"degree={diff_order(gd)}")

# the definitional-content audit: two-line argument, checked computationally
ent = [sum(1 for g in cells_of(n) if channel_count(g, set(cells_of(n)))) for n in NS]
chk("AUDIT (definitional content): capacity == #channel-bearing grains exactly "
    "(upper: grain-disjoint chains use distinct entry grains; lower: each entry "
    "grain certifies itself by a length-1 chain) -- the per-epoch boundary bound "
    "under THIS definition follows from the entry-set size; what is MEASURED is "
    "that the entry set is boundary-sized in this geometry (the controls vary it)",
    all(v == e for v, e in zip(gd, ent)))

NS_ed = list(range(2, 12))
ed = [edge_disjoint_cap(n) for n in NS_ed]
say("    edge-disjoint:   " + " ".join(f"{v:>5}" for v in ed))
say("    (channel mult.)")
chk("alternative (a) EDGE-disjoint: capacity == min(n^3, 6n^2) at every n=2..11 -- the "
    "binding cut moves from the grain layer to the CHANNEL COUNT; below n=6 the whole "
    "volume fits under the channel count, a small-n saturation, not a scaling",
    all(v == min(n**3, 6*n*n) for v, n in zip(ed, NS_ed)))
o_ed = diff_order(ed[4:])   # n >= 6, the regime where 6n^2 <= n^3
chk("alternative (a) on the uncapped regime n=6..11: degree 2 BOUNDARY -- no flip",
    o_ed == 2 and all(v == 6*n*n for v, n in zip(ed[4:], NS_ed[4:])),
    f"degree(n>=6)={o_ed}")

cd = [channel_disjoint_cap(n) for n in NS + [10, 11, 12]]
say("    channel-only:    " + " ".join(f"{v:>5}" for v in cd))
o_cd = diff_order(cd[4:])   # n >= 6, where 6n^2 <= n^3: the uncapped regime
chk("alternative (b) CHANNEL-disjoint (weakest honest independence): capacity = "
    "min(6n^2, n^3); on n >= 6 it is 6n^2, degree 2 BOUNDARY -- even the weakest "
    "reading does not go volume per epoch", o_cd == 2, f"degree(n>=6)={o_cd}")

chk("RIG VERDICT: no honest per-epoch operationalisation tried reaches degree 3 in its "
    "uncapped regime; the boundary-RATE arm is not an artifact of the grain-disjoint "
    "choice", (diff_order(gd) == 2) and o_ed == 2 and o_cd == 2)

say("")
say("cumulative arm under alternatives: MORE permissive independence certifies the "
    "volume at least as fast (any grain-disjoint schedule is also edge-/channel-legal), "
    "so E_min <= 3 holds under every alternative; the volume-TOTAL arm is "
    "definition-robust in the permissive direction by set inclusion (computed check "
    "at n = 4: edge-disjoint 2-epoch optimum >= grain-disjoint's 64):")
ed_n4_e2 = None
# 2-epoch edge-disjoint optimum at n=4 via time expansion with edge capacities
def n_opt_edge(n, E):
    cells = cells_of(n); S = set(cells)
    f = Flow()
    for e in range(E):
        for g in cells:
            for d in NB:
                h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
                if h in S:
                    f.add(('v', e, g), ('v', e, h), 1)
            f.add(('v', e, g), ('cert', g), 1)
            c = channel_count(g, S)
            if c:
                f.add('S', ('v', e, g), c)
    for g in cells:
        f.add(('cert', g), 'T', 1)
    return f.maxflow('S', 'T')
ed_n4_e2 = n_opt_edge(4, 2)
chk("n=4, E=2 edge-disjoint optimum >= 64 (= whole volume, matching grain-disjoint)",
    ed_n4_e2 >= 64, f"value={ed_n4_e2}")

say("")
say("=" * 96)
say("A2  VOLUME-CONTROL SENSITIVITY -- can the instrument actually say volume?  and 1?")
say("=" * 96)

fc = [grain_disjoint_cap(n, entry=cells_of(n)) for n in NS]
say("    fully-channelled: " + " ".join(f"{v:>5}" for v in fc))
chk("fully-channelled CERT/epoch == n^3 (independent instrument), degree 3 VOLUME",
    all(v == n**3 for v, n in zip(fc, NS)) and diff_order(fc) == 3,
    f"degree={diff_order(fc)}")

# control the lane did NOT run: exactly one channel (a corner grain)
sc = [grain_disjoint_cap(n, entry=[(0, 0, 0)]) for n in (2, 3, 4, 5)]
chk("single-channel block: CERT/epoch == 1 at n=2..5 -- the instrument spans the whole "
    "range down to 1; boundary reads are not an instrument floor", all(v == 1 for v in sc),
    f"values={sc}")
sc_curve = [n_opt_indep(3, E, entry=[(0, 0, 0)]) for E in (1, 5, 26, 27)]
chk("single-channel n=3: N_opt(E) == E for E=1,5 and == min(E,27) at E=26,27 "
    "(cumulative volume needs n^3 epochs at rate 1 -- rate bound arithmetic exact)",
    sc_curve == [1, 5, 26, 27], f"curve at E=1,5,26,27: {sc_curve}")

say("")
say("=" * 96)
say("A3  SPLIT HONESTY -- both arms recomputed independently")
say("=" * 96)

# per-epoch arm on the full published range
gd_full = gd + [grain_disjoint_cap(n) for n in (10, 11, 12)]
chk("per-epoch arm, full range n=2..12: CERT/epoch == 6n^2-12n+8, degree 2",
    all(v == 6*n*n - 12*n + 8 for v, n in zip(gd_full, range(2, 13)))
    and diff_order(gd_full) == 2)

n2 = {n: n_opt_indep(n, 2) for n in range(8, 13)}
say("    N_opt(2), n=8..12 (independent): " + str([n2[n] for n in range(8, 13)]))
chk("N_opt(2) == 12n^2 - 36n + 16 at n = 8..12 (independent time-expanded flow)",
    all(n2[n] == 12*n*n - 36*n + 16 for n in range(8, 13)))
chk("N_opt(2) second difference constant 24 (degree 2: even a 2-epoch budget is "
    "boundary-order)", diff_order([n2[n] for n in range(8, 13)]) == 2)

n3 = {n: n_opt_indep(n, 3) for n in (8, 12)}
chk("cumulative arm: N_opt(3) == n^3 at n = 8 and n = 12 (whole stored volume "
    "certified at E_min = 3; independent instrument)",
    n3[8] == 512 and n3[12] == 1728, f"n=8:{n3[8]}, n=12:{n3[12]}")
emin_small = {n: next(E for E in range(1, 5) if n_opt_indep(n, E) == n**3)
              for n in range(2, 8)}
chk("E_min table n=2..7 == 1,2,2,2,2,2 (independent)",
    [emin_small[n] for n in range(2, 8)] == [1, 2, 2, 2, 2, 2],
    str(emin_small))
chk("E_min rate-bound arithmetic: ceil(n^3/(6n^2-12n+8)) = "
    + str([-(-n**3 // (6*n*n-12*n+8)) for n in range(2, 13)])
    + " <= published E_min = [1,2,2,2,2,2,3,3,3,3,3] at every n",
    all(-(-n**3 // (6*n*n-12*n+8)) <= e
        for n, e in zip(range(2, 13), [1,2,2,2,2,2,3,3,3,3,3])))

say("")
say("=" * 96)
say("A5  TABLE CROSS-CHECK -- finding's table vs independent numbers")
say("=" * 96)
finding = {2:(8,24,8,1), 3:(27,54,26,2), 4:(64,96,56,2), 5:(125,150,98,2),
           6:(216,216,152,2), 7:(343,294,218,2), 8:(512,384,296,3), 9:(729,486,386,3),
           10:(1000,600,488,3), 11:(1331,726,602,3), 12:(1728,864,728,3)}
chk("STORED/IFACE/CERT columns match the closed forms at every published n",
    all(finding[n] == (n**3, 6*n*n, 6*n*n-12*n+8, (1 if n == 2 else 2 if n < 8 else 3))
        for n in finding))
chk("finding curve entries n=8 [296,496,512] and n=12 [728,1312,1728] match "
    "independent N_opt", n2[8] == 496 and n3[8] == 512 and n2[12] == 1312
    and n3[12] == 1728 and grain_disjoint_cap(8) == 296)

say("")
say(f"ALL CHECKS PASS: {ALL}")
sys.exit(0 if ALL else 1)
