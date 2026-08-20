"""T-43-B   EXTERNALLY CERTIFIABLE CONTENT -- WORLD TIER, census access model, exact counting.

THE PRE-REGISTERED QUESTION (REGISTER_V001, T-43 rule, registered before this ran): of the
volume-many records a region stores, how many can the outside world independently CERTIFY
(learn the value of) or WRITE, through the interface?  This lane computes the WORLD-tier
numbers.  It reports; the rule decides at the register, across both tiers, not here.

MODEL (INSERTED, declared): the GR1/T42_C access model.  Grains at the integer points of an
n x n x n block; one barrier record per grain (census, LANE_GR1_CENSUS); adjacency = shares a
face; the outside touches the block ONLY through its 6n^2 block-complement face adjacencies
(the interface channels, LANE_T42_C_BOUNDCAP).  An outside operation enters through a channel
and traverses grains in adjacency steps.  Nothing continuous; exact integers throughout; no
fit anywhere (scaling read by constant finite differences, the T42_C instrument).

OPERATIONAL DEFINITIONS ON THE SURFACE (D-24 -- defined in access-graph terms, not imported):

  CHAIN         a path of pairwise-distinct grains g_1, ..., g_d with g_1 channel-bearing
                (adjacent to the complement) and consecutive grains adjacent.  d = the chain's
                traversal cost; the minimal d for a grain is its BFS DEPTH (computed, T42_C).
  CERTIFY       the outside learns the stored value of g_d when a chain ending at g_d relays
                the value out through g_1's channel.  Relay grains are read THROUGH; a chain
                certifies only its endpoint.
  WRITE         the same chain run inward, flipping the record at g_d.  The access model does
                not distinguish read-traversal from write-traversal: WRITABLE-per-epoch is the
                SAME counting problem BY MODEL STRUCTURE (declared, not a finding).
  INDEPENDENT   two certificates are independent when their chains share no grain: a shared
                relay would make both certificates depend on one grain's faithful relaying.
                (Grain-disjointness also retires the channel-multiplicity question: a corner
                grain has 3 channels but can host at most one chain.)
  EPOCH         (model choice, DECLARED): one synchronous round in which each grain
                participates in at most one chain and each channel carries at most one chain;
                a depth-d certification completes within the epoch it is scheduled in.
                The alternative reading (one adjacency step per tick, a depth-d chain costing
                d ticks) is declared beside it; under that reading every per-epoch count below
                stays a per-round count and the time axis stretches by at most the max depth
                (<= ceil(n/2)) -- the per-epoch vs cumulative split reported below is the same
                split in either reading.

INSTRUMENTS (owners named).  Simultaneous independent certification = grain-disjoint chains
from the complement to distinct targets = unit-vertex-capacity max-flow; max-flow / min-cut
and the disjoint-path equivalence are BORROWED (Menger 1927; Ford-Fulkerson 1956; Dinic
1970).  The E-epoch optimum is computed on a TIME-EXPANDED network (E copies of the block
network, one per epoch, coupled only through a shared unit certification sink per grain);
time-expanded flow networks are BORROWED (Ford-Fulkerson flows-over-time territory).  OURS:
the operational certify/write framing, the per-epoch = traversal-free-capacity identity, the
shell-cut epoch bounds and the computed places they are NOT tight, the exact epoch-budget
optimum N_opt(E) and its closed forms, and the two-tier test under the pre-registered rule.
D-1: no classical gravitational form appears anywhere in this file.

D-15 CONTROLS (the brief's two): (1) FULLY-CHANNELLED block -- a channel to every grain; the
certifiable-per-epoch count must go VOLUME, so a boundary bound is a measurement, not
blindness of the instrument.  (2) SHUFFLED channel graph -- the same 6n^2 channels attached
to seeded-random grains; the depth structure must collapse and the per-epoch count must track
the channel-bearing grain count, not the block's geometry.

NO NARRATED VERDICTS: every PASS/FAIL line is gated by a computed boolean.
"""
import sys
from collections import deque
from itertools import groupby, permutations, product

sys.setrecursionlimit(200000)

def say(*a):
    print(*a); sys.stdout.flush()

def verdict(ok):   # gate: computed boolean -> word.  Never called with a literal.
    return "PASS" if bool(ok) else "FAIL"

ALL = True
def gate(label, ok, detail=""):
    global ALL
    ALL &= bool(ok)
    say(f"  [{verdict(ok)}] {label}" + (f"  {detail}" if detail else ""))
    return bool(ok)

# ------------------------------------------------------------------ max-flow (Dinic, exact)
class Dinic:
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]
    def add(self, u, v, c):
        self.g[u].append([v, c, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])
        return (u, len(self.g[u]) - 1)
    def _bfs(self, s, t):
        self.level = [-1] * self.n
        self.level[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v, c, _ in self.g[u]:
                if c > 0 and self.level[v] < 0:
                    self.level[v] = self.level[u] + 1
                    q.append(v)
        return self.level[t] >= 0
    def _dfs(self, u, t, f):
        if u == t:
            return f
        while self.it[u] < len(self.g[u]):
            e = self.g[u][self.it[u]]
            v, c, r = e
            if c > 0 and self.level[v] == self.level[u] + 1:
                d = self._dfs(v, t, min(f, c))
                if d > 0:
                    e[1] -= d
                    self.g[v][r][1] += d
                    return d
            self.it[u] += 1
        return 0
    def maxflow(self, s, t):
        fl = 0
        while self._bfs(s, t):
            self.it = [0] * self.n
            while True:
                f = self._dfs(s, t, 1 << 30)
                if not f:
                    break
                fl += f
        return fl

# ------------------------------------------------------------------ the access model
NB = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))

def block_cells(n):
    return sorted((x, y, z) for x in range(n) for y in range(n) for z in range(n))

def channels_block(cells):
    """GR1/T42_C interface: one channel per block-complement face adjacency."""
    S = set(cells)
    ch = {}
    for g in cells:
        c = sum(1 for d in NB if (g[0]+d[0], g[1]+d[1], g[2]+d[2]) not in S)
        if c:
            ch[g] = c
    return ch

def depths(cells, entry):
    """BFS chain-length depth from the channel-bearing entry set (depth 1 at entry)."""
    S = set(cells)
    dep = {g: 1 for g in entry}
    frontier = sorted(entry)
    while frontier:
        nxt = []
        for g in frontier:
            for d in NB:
                h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
                if h in S and h not in dep:
                    dep[h] = dep[g] + 1
                    nxt.append(h)
        frontier = sorted(nxt)
    assert set(dep) == S
    return dep

def sim_capacity(cells, entry, targets):
    """One epoch: max number of grain-disjoint chains to distinct targets (exact max-flow)."""
    idx = {g: i for i, g in enumerate(cells)}
    dn = Dinic(2 + 2 * len(cells))
    S, T = 0, 1
    Sset = set(cells)
    for g in cells:
        i = idx[g]
        dn.add(2 + 2*i, 3 + 2*i, 1)                    # each grain hosts <= 1 chain
        for d in NB:
            h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
            if h in Sset:
                dn.add(3 + 2*i, 2 + 2*idx[h], 1)       # traversal step
    for g in entry:
        dn.add(S, 2 + 2*idx[g], 1)                     # channel entry
    for g in targets:
        dn.add(3 + 2*idx[g], T, 1)
    return dn.maxflow(S, T)

def n_opt(cells, entry, E):
    """EXACT optimum over ALL E-epoch schedules: max grains certifiable within E epochs.
       Time-expanded network: E copies of the block network (grain-disjointness enforced
       per copy), coupled ONLY through one shared unit certification sink per grain.  A
       flow decomposition IS a schedule; max-flow = the optimum; the min cut is its
       certificate.  Exact integers throughout."""
    idx = {g: i for i, g in enumerate(cells)}
    C = len(cells)
    dn = Dinic(2 + 2 * E * C + C)
    S, T = 0, 1
    cert = lambda i: 2 + 2 * E * C + i
    Sset = set(cells)
    for e in range(E):
        base = e * C
        IN = lambda i: 2 + 2 * (base + i)
        OUT = lambda i: 3 + 2 * (base + i)
        for g in cells:
            i = idx[g]
            dn.add(IN(i), OUT(i), 1)
            dn.add(OUT(i), cert(i), 1)
            for d in NB:
                h = (g[0]+d[0], g[1]+d[1], g[2]+d[2])
                if h in Sset:
                    dn.add(OUT(i), IN(idx[h]), 1)
        for g in entry:
            dn.add(S, IN(idx[g]), 1)
    for i in range(C):
        dn.add(cert(i), T, 1)
    return dn.maxflow(S, T)

def epoch_curve(cells, entry):
    """N_opt(E) for E = 1.. until the whole content is certified; returns the curve."""
    vals = []
    E = 1
    while True:
        v = n_opt(cells, entry, E)
        vals.append(v)
        if v == len(cells):
            return vals
        E += 1

def const_diff_order(seq):
    """T42_C instrument: smallest k with the k-th finite difference constant; None if none."""
    s = list(seq)
    for k in range(len(seq) - 1):
        if len(set(s)) == 1:
            return k
        s = [b - a for a, b in zip(s, s[1:])]
    return None

# =====================================================================================
say("=" * 100)
say("T-43-B   EXTERNALLY CERTIFIABLE CONTENT -- WORLD TIER (exact counts, n = 2..12)")
say("=" * 100)

NS = list(range(2, 13))
rows = {}
for n in NS:
    cells = block_cells(n)
    ch = channels_block(cells)
    entry = sorted(ch)
    dep = depths(cells, entry)
    maxd = max(dep.values())
    L = {d: sum(1 for g in cells if dep[g] == d) for d in range(1, maxd + 1)}
    Vdeep = {k: sum(1 for g in cells if dep[g] >= k) for k in range(1, maxd + 1)}
    sim = sim_capacity(cells, entry, cells)
    sim_shell = {d: sim_capacity(cells, entry, [g for g in cells if dep[g] == d])
                 for d in L}
    sim_deep = {k: sim_capacity(cells, entry, [g for g in cells if dep[g] >= k])
                for k in Vdeep}
    curve = epoch_curve(cells, entry)
    LB = max(-(-Vdeep[k] // sim_deep[k]) for k in Vdeep)     # exact ceil, shell-cut bound
    rows[n] = dict(cells=cells, ch=ch, entry=entry, dep=dep, maxd=maxd, L=L,
                   Vdeep=Vdeep, sim=sim, sim_shell=sim_shell, sim_deep=sim_deep,
                   curve=curve, Emin=len(curve), LB=LB,
                   iface=sum(ch.values()), cap_d1=len(entry),
                   depth_sum=sum(dep.values()))

say("")
say("(a)+(b)+(c)  THE THREE COUNTS SIDE BY SIDE, per n  (all exact):")
say(f"{'n':>3} {'STORED':>7} {'IFACE':>6} {'CERT/epoch':>10} {'WRITE/epoch':>11} "
    f"{'E_min':>6} {'LB_shell':>8} {'DEPTH_SUM':>10}")
for n in NS:
    r = rows[n]
    say(f"{n:>3} {n**3:>7} {r['iface']:>6} {r['sim']:>10} {r['sim']:>11} "
        f"{r['Emin']:>6} {r['LB']:>8} {r['depth_sum']:>10}")

say("")
say("GATES -- model replication and the per-epoch capacity:")
gate("IFACE == 6n^2 at every n (the inserted channel count, replicated)",
     all(rows[n]['iface'] == 6 * n * n for n in NS))
gate("CAP_depth(1) == n^3 - (n-2)^3 at every n (T42_C replicated)",
     all(rows[n]['cap_d1'] == n**3 - (n - 2)**3 for n in NS))
gate("DEPTH_SUM recurrence DEPTH_SUM(n) == DEPTH_SUM(n-2) + n^3, every n >= 4 (T42_C)",
     all(rows[n]['depth_sum'] == rows[n - 2]['depth_sum'] + n**3 for n in NS[2:]))
gate("shell counts: L(d) == (n-2d+2)^3 - (n-2d)^3 (clamped at 0), every n, every d",
     all(rows[n]['L'][d] == max(n - 2*d + 2, 0)**3 - max(n - 2*d, 0)**3
         for n in NS for d in rows[n]['L']))
say("")
say("THE PER-EPOCH CAPACITY (simultaneous independent certification, exact max-flow):")
gate("CERT/epoch == n^3 - (n-2)^3 == CAP_depth(1) at every n -- THE IDENTITY: allowing "
     "chains to any depth buys NOTHING per epoch; grain-disjointness pins the simultaneous "
     "capacity to the traversal-free capacity",
     all(rows[n]['sim'] == rows[n]['cap_d1'] for n in NS))
gate("CERT/epoch is degree 2 in n (BOUNDARY), read by constant finite differences",
     const_diff_order([rows[n]['sim'] for n in NS]) == 2,
     f"order={const_diff_order([rows[n]['sim'] for n in NS])}")
gate("CERT/epoch < IFACE at every n (the channel bound 6n^2 is never tight; the binding "
     "cut is the first grain layer)",
     all(rows[n]['sim'] < rows[n]['iface'] for n in NS))
gate("STORED content is degree 3 (VOLUME) -- the census, replicated beside the capacity",
     const_diff_order([n**3 for n in NS]) == 3)
gate("IFACE is degree 2 (BOUNDARY)",
     const_diff_order([rows[n]['iface'] for n in NS]) == 2)
say("")
say("(b) WRITABLE-FROM-OUTSIDE per epoch: the identical counting problem BY MODEL STRUCTURE")
say("    (the access model does not distinguish read-traversal from write-traversal; declared")
say("    in the header, stated as structure, not a finding).  WRITE/epoch column above.")

# ------------------------------------------------------------------ depth structure
say("")
say("=" * 100)
say("DEPTH STRUCTURE -- certification cost vs depth, and the shell cuts")
say("=" * 100)
say("")
say("Cost of one certification at depth d = d (the chain is the cost; depth is the minimal")
say("chain length, computed by BFS).  Cost to certify everything once = DEPTH_SUM (table).")
say("")
n = 12
r = rows[n]
say(f"witness table, n = {n}:  shell d | L(d) | SIM_shell(d) | V_deep(d) | SIM_deep(d)")
for d in sorted(r['L']):
    say(f"    {d:>2}    {r['L'][d]:>6} {r['sim_shell'][d]:>10} "
        f"{r['Vdeep'][d]:>10} {r['sim_deep'][d]:>9}")
gate("SIM_shell(d) == L(d) at every n and d: EVERY shell is certifiable whole in ONE epoch, "
     "at chain cost d per record",
     all(rows[n]['sim_shell'][d] == rows[n]['L'][d] for n in NS for d in rows[n]['L']))
gate("SIM_deep(k) == L(k) at every n and k: the k-th shell is the exact min-cut for all "
     "content at depth >= k -- per epoch, at most L(k) records deeper than k-1 are "
     "certifiable, however scheduled",
     all(rows[n]['sim_deep'][k] == rows[n]['L'][k] for n in NS for k in rows[n]['Vdeep']))
gate("V_deep(k) == (n-2k+2)^3 (clamped): content at depth >= k is the inner cube -- "
     "deep certifiable content shrinks as a VOLUME two layers thinner per step",
     all(rows[n]['Vdeep'][k] == max(n - 2*k + 2, 0)**3
         for n in NS for k in rows[n]['Vdeep']))

# ------------------------------------------------------------------ epoch budgets
say("")
say("=" * 100)
say("FIXED EPOCH BUDGET -- N_opt(E), the EXACT effective certifiable content, no heuristics")
say("=" * 100)
say("")
say("N_opt(E) = the optimum over ALL E-epoch schedules, computed as ONE exact max-flow on")
say("the time-expanded network (header: instruments).  Beside it, the shell-cut bound")
say("UB_shell(E) = min over k of [content at depth < k  +  E * SIM_deep(k)], capped at n^3.")
say("Scheduling heuristics tried first (shell round-robin, z-slabs, local search) all fell")
say("short of the cuts; the time-expanded flow replaced them with the true optimum.")
say("")

def UB_shell(n, E):
    r = rows[n]
    best = min((n**3 - r['Vdeep'][k]) + E * r['sim_deep'][k] for k in r['Vdeep'])
    return min(best, n**3)

say(f"{'n':>3}  N_opt(E) for E=1..E_min  [beside each: /UB_shell(E); * marks UB not met]")
strict_pairs = []
for n in NS:
    r = rows[n]
    cur = []
    for E, v in enumerate(r['curve'], start=1):
        ub = UB_shell(n, E)
        cur.append(f"{v}/{ub}" + ("*" if v < ub else ""))
        if v < ub:
            strict_pairs.append((n, E))
    say(f"{n:>3}  " + "  ".join(cur))
say("")
gate("N_opt(1) == CERT/epoch at every n (the two instruments agree where they overlap)",
     all(rows[n]['curve'][0] == rows[n]['sim'] for n in NS))
gate("N_opt(E) <= UB_shell(E) always (the shell cuts are valid bounds)",
     all(v <= UB_shell(n, E) for n in NS
         for E, v in enumerate(rows[n]['curve'], start=1)))
gate("the shell cuts are NOT the whole story: N_opt(E) < UB_shell(E) exactly at the "
     "computed pairs (n,E) = " + str(strict_pairs) + " -- at every such pair the "
     "time-expanded min cut is strictly inside the shell-cut family (reported as found)",
     strict_pairs == [(8, 2), (9, 2), (10, 2), (11, 2), (12, 2)])
gate("N_opt(2) on the uncapped range n = 8..12 has CONSTANT SECOND DIFFERENCE 24: "
     "degree 2 (BOUNDARY), exact form 2*(6n^2-12n+8) - 12n == 12n^2 - 36n + 16 -- two "
     "epochs pay an exact 12n congestion toll against twice the one-epoch capacity",
     const_diff_order([rows[n]['curve'][1] for n in range(8, 13)]) == 2 and
     all(rows[n]['curve'][1] == 12*n*n - 36*n + 16 for n in range(8, 13)))
gate("N_opt(E) >= n^3 - (n-2E)^3 (clamped) at every n, E computed: E epochs certify at "
     "least the outer 2E layers' volume (the shell-by-shell schedule is never better)",
     all(v >= n**3 - max(n - 2*E, 0)**3 for n in NS
         for E, v in enumerate(rows[n]['curve'], start=1)))
say("")
say("MINIMAL EPOCHS TO CERTIFY THE WHOLE STORED CONTENT (exact, from the curve):")
say("    E_min: " + ", ".join(f"{n}:{rows[n]['Emin']}" for n in NS))
gate("LB_shell <= E_min at every n (cut validity)",
     all(rows[n]['LB'] <= rows[n]['Emin'] for n in NS))
gate("E_min == LB_shell everywhere EXCEPT the computed exceptions {n=8, n=9}, where "
     "E_min = 3 > LB_shell = 2: the shell cuts alone under-count the epochs needed "
     "(consistent with the N_opt(2) toll: n^3 > 12n^2-36n+16 first at n = 8)",
     all((rows[n]['Emin'] == rows[n]['LB']) == (n not in (8, 9)) for n in NS)
     and all(rows[n]['Emin'] == 3 and rows[n]['LB'] == 2 for n in (8, 9)))
gate("E_min >= ceil(n^3 / CERT_epoch) at every n (arithmetic floor of the rate bound; "
     "over the computed range E_min stays within +1 of it)",
     all(rows[n]['Emin'] >= -(-n**3 // rows[n]['sim']) and
         rows[n]['Emin'] <= -(-n**3 // rows[n]['sim']) + 1 for n in NS))
say("")
say("PER-EPOCH vs CUMULATIVE -- the split the brief requires, kept explicit (both exact):")
say("  * PER EPOCH:   certifiable content is BOUNDARY-bounded: 6n^2 - 12n + 8 (degree 2),")
say("    and no schedule beats it (min-cut).  Per epoch the interface is the whole story.")
say("  * CUMULATIVE:  re-use over epochs is NOT boundary-bounded: N_opt(E) grows with time")
say("    (never faster than E * (6n^2-12n+8), and with the computed 12n toll at E = 2), and")
say("    at E_min <= 3 on the whole computed range the ENTIRE stored volume n^3 is")
say("    certified.  Certifiable CONTENT, given time, is VOLUME; the boundary bounds the")
say("    RATE, not the total.")

# ------------------------------------------------------------------ controls
say("")
say("=" * 100)
say("D-15 CONTROLS -- the instrument must be able to see volume, and to lose the geometry")
say("=" * 100)

say("")
say("CONTROL 1 -- FULLY-CHANNELLED block: a channel to EVERY grain (punctured geometry).")
fc_sim = []
fc_E = []
for n in NS:
    cells = block_cells(n)
    entry = cells
    fc_sim.append(sim_capacity(cells, entry, cells))
    fc_E.append(len(epoch_curve(cells, entry)))
say("    n:          " + " ".join(f"{n:>5}" for n in NS))
say("    CERT/epoch: " + " ".join(f"{s:>5}" for s in fc_sim))
gate("fully-channelled CERT/epoch == n^3 at every n, degree 3 (VOLUME): when the outside "
     "reaches every cell, the instrument REPORTS volume -- the boundary bound above is a "
     "measurement of the access geometry, not blindness",
     all(s == n**3 for s, n in zip(fc_sim, NS)) and const_diff_order(fc_sim) == 3,
     f"order={const_diff_order(fc_sim)}")
gate("fully-channelled E_min == 1 at every n (everything certified in one epoch)",
     all(e == 1 for e in fc_E))

say("")
say("CONTROL 2 -- SHUFFLED channel graph: the same 6n^2 channels attached to seeded-random")
say("grains (deterministic seed 11, drawn with repetition from the sorted cell list).")
import random
sh_rows = []
for n in NS:
    cells = block_cells(n)
    rnd = random.Random(11)
    draws = [cells[rnd.randrange(len(cells))] for _ in range(6 * n * n)]
    entry = sorted(set(draws))
    dep = depths(cells, entry)
    s = sim_capacity(cells, entry, cells)
    curve = epoch_curve(cells, entry)
    sh_rows.append((n, len(draws), len(entry), s, max(dep.values()), len(curve)))
say(f"    {'n':>3} {'channels':>8} {'distinct':>8} {'CERT/epoch':>10} {'max depth':>9} "
    f"{'E_min':>6}")
for t in sh_rows:
    say(f"    {t[0]:>3} {t[1]:>8} {t[2]:>8} {t[3]:>10} {t[4]:>9} {t[5]:>6}")
gate("shuffled: CERT/epoch == number of DISTINCT channel-bearing grains at every n (one "
     "chain per channel grain; grain-disjointness; the count follows the channel graph, "
     "not the block's shape)",
     all(t[3] == t[2] for t in sh_rows))
gate("shuffled: the depth structure COLLAPSES -- max depth <= 4 at every n (block: "
     "ceil(n/2), = 6 at n=12); the boundary-vs-bulk story above lives in the geometry the "
     "shuffle destroys",
     all(t[4] <= 4 for t in sh_rows),
     "block max depths: " + ", ".join(f"{n}:{rows[n]['maxd']}" for n in NS))
gate("shuffled: E_min <= 3 at every n (with channels everywhere in the bulk, the whole "
     "volume is certified in a depth-free handful of epochs)",
     all(t[5] <= 3 for t in sh_rows),
     "E_min: " + ", ".join(f"{t[0]}:{t[5]}" for t in sh_rows))

# ------------------------------------------------------------------ venue check
say("")
say("=" * 100)
say("D-22 VENUE / INSTRUMENT-INVARIANCE CHECK (n = 5 witness)")
say("=" * 100)
n = 5
cells = block_cells(n)
ch = channels_block(cells)
dep = depths(cells, sorted(ch))
sym_ok = True
count = 0
for perm in permutations(range(3)):
    for signs in product((1, -1), repeat=3):
        def m(g, perm=perm, signs=signs):
            w = [g[perm[0]], g[perm[1]], g[perm[2]]]
            return tuple((w[i] if signs[i] == 1 else (n - 1) - w[i]) for i in range(3))
        img = {m(g) for g in cells}
        ok = (img == set(cells)) and all(dep[m(g)] == dep[g] for g in cells)
        sym_ok &= ok
        count += 1
gate(f"all {count} signed-permutation symmetries of the block preserve the cell set and "
     "the computed depth field (the instrument respects the venue's own symmetry)",
     sym_ok and count == 48)

# ------------------------------------------------------------------ reading
say("")
say("=" * 100)
say("READING FOR T-43 (world tier only; the rule decides at the register, across tiers)")
say("=" * 100)
say("""
  STORED content:                 n^3            VOLUME  (degree 3; census, replicated)
  CERTIFIABLE per epoch:          6n^2-12n+8     BOUNDARY (degree 2; exact min-cut; equals
                                                 the traversal-free capacity CAP_depth(1) --
                                                 chains to depth buy nothing per epoch)
  WRITABLE per epoch:             same count     (model structure, declared)
  CERTIFIABLE, budget E:          N_opt(E)       exact (time-expanded flow); at E = 2 it is
                                                 12n^2-36n+16 on the uncapped range --
                                                 still BOUNDARY order, with a 12n toll
  CERTIFIABLE cumulative:         whole n^3      at E_min <= 3 on the computed range;
                                                 E_min >= ceil(n^3/(6n^2-12n+8)) always
  Depth cost:                     d per record at depth d; shell k is the exact per-epoch
                                                 min-cut (L(k)) for everything deeper.

  The world tier's honest answer under the pre-registered rule is SPLIT, and the split is
  the physics: the boundary bounds the certification RATE exactly (per-epoch and even
  per-epoch-pair capacity are boundary-scaling, and the fully-channelled control shows the
  instrument would have said volume if the access geometry allowed it); the boundary does
  NOT bound the certifiable TOTAL (re-use over epochs certifies the whole stored volume,
  in at most 3 epochs on the computed range, at a rate that keeps E_min >= ~n/6 as n
  grows -- the cut arithmetic, stated as a bound, not an extrapolation).  Whether
  'externally certifiable content' in the rule's sense means the per-epoch count
  (boundary-bounded) or the anytime count (volume) is exactly the per-epoch vs cumulative
  declaration of D-24 -- reported here as found, both exact, NOT collapsed.
""")
say(f"ALL GATES PASS: {ALL}")
sys.exit(0 if ALL else 1)
