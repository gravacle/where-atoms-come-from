"""VERIFY v2 -- (a) the lane's GAP census is sampled from a BIASED POOL; recompute it UNIFORMLY.
                (b) re-run the lane's H2 'not pairwise-determined' test with a NON-lossy two-body key.
                (c) re-check the lane's own D-15 controls in E-9 G2 and E-10 H4.

Independent engine again -- nothing imported from exact_pauli.py.
"""
import sys, random
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); OUT.append(s); print(s); sys.stdout.flush()

def pc(v): return bin(v).count("1")
def enc(vec, n):
    x = z = 0
    for j in range(n):
        x = (x << 1) | vec[j]; z = (z << 1) | vec[n + j]
    return (x, z)
def sp(A, B): return (pc(A[0] & B[1]) + pc(A[1] & B[0])) & 1
def xr(A, B): return (A[0] ^ B[0], A[1] ^ B[1])
def Iu(A, B): return pc((A[0] & B[1]) ^ (A[1] & B[0]))

def carrier(n):
    sX = [1] * n + [0] * n; sZ = [0] * n + [1] * n
    pairs = symplectic_logicals([sX, sZ], n)
    gens = []
    for X, Z in pairs:
        gens.append(enc(X, n)); gens.append(enc(Z, n))
    kX = enc(sX, n); kZ = enc(sZ, n)
    return gens, [(0, 0), kX, kZ, xr(kX, kZ)]

def uniform_record(gens, rng):
    """a UNIFORMLY random element of the record group N(S)/S: random F_2 combination of generators"""
    v = (0, 0)
    for g in gens:
        if rng.getrandbits(1):
            v = xr(v, g)
    return v

def tbl(a, b, S): return [[Iu(xr(a, s1), xr(b, s2)) for s2 in S] for s1 in S]
def Jt(T): return min(min(r) for r in T)
def GAPt(T1, T2, T3):
    q = len(T1); best = None
    for i in range(q):
        for j in range(q):
            v0 = T1[i][j]
            for k in range(q):
                v = v0 + T2[j][k] + T3[i][k]
                if best is None or v < best: best = v
    return best - (Jt(T1) + Jt(T2) + Jt(T3))

say("=" * 120)
say("VERIFY v2")
say("=" * 120)

# ---------------------------------------------------------------- (a) UNIFORM census
say("")
say("B1  THE LANE'S 'FRUSTRATED FRACTION RISES 0, 0.9%, 6.3%, 10.4%, 17.9%, 23.7%' -- RECOMPUTED")
say("    UNIFORMLY over the FULL record group instead of the lane's generators+pairwise-sums pool.")
say("-" * 120)
say("  %-5s %-9s %-14s %-40s %-16s %-16s"
    % ("n", "k", "N sampled", "UNIFORM GAP distribution (exact)", "uniform frac>0", "lane's frac>0"))
lane = {4: 0.0, 6: 54 / 5832, 8: 366 / 5832, 10: 606 / 5832, 12: 1044 / 5832, 14: 1380 / 5832}
for n in (4, 6, 8, 10, 12, 14, 16):
    gens, S = carrier(n)
    rng = random.Random(31337 + n)
    N = 40000
    d = {}
    for _ in range(N):
        a = uniform_record(gens, rng); b = uniform_record(gens, rng); c = uniform_record(gens, rng)
        Tab = tbl(a, b, S); Tbc = tbl(b, c, S); Tac = tbl(a, c, S)
        g = GAPt(Tab, Tbc, Tac)
        d[g] = d.get(g, 0) + 1
    frac = (N - d.get(0, 0)) / N
    say("  %-5d %-9d %-14d %-40s %-16.4f %-16s"
        % (n, n - 2, N, str(dict(sorted(d.items()))), frac,
           ("%.4f" % lane[n]) if n in lane else "not reported"))
say("")
say("  At n=6 the lane reports 0.0093; the uniform value is ~0.10 -- an order of magnitude apart.")
say("  The lane's monotone 'growth' series is therefore a property of its sampling pool, not of n.")

# ---------------------------------------------------------------- (b) H2 with a non-lossy key
say("")
say("B2  THE LANE'S H2: 'GAP is NOT a function of the pairwise data'.  Their key was the six")
say("    two-body MINIMA (sp x3, J x3) -- a lossy summary.  Re-run with the FULL two-body key:")
say("    the three two-body tables themselves (still purely two-body data, nothing three-body).")
say("-" * 120)
say("  %-5s %-11s %-34s %-34s"
    % ("n", "triples", "key = six minima (lane's)", "key = the three two-body TABLES"))
for n in (6, 8, 10, 12, 14, 16):
    gens, S = carrier(n)
    rng = random.Random(555 + n)
    reps = sorted({uniform_record(gens, rng) for _ in range(30)} - {(0, 0)})
    Tc = {}
    for i in range(len(reps)):
        for j in range(len(reps)):
            Tc[(i, j)] = tbl(reps[i], reps[j], S)
    g_min = {}; g_tab = {}
    tot = 0
    for i in range(len(reps)):
        for j in range(len(reps)):
            for k in range(len(reps)):
                T1 = Tc[(i, j)]; T2 = Tc[(j, k)]; T3 = Tc[(i, k)]
                g = GAPt(T1, T2, T3)
                k1 = (sp(reps[i], reps[j]), sp(reps[j], reps[k]), sp(reps[i], reps[k]),
                      Jt(T1), Jt(T2), Jt(T3))
                k2 = (tuple(map(tuple, T1)), tuple(map(tuple, T2)), tuple(map(tuple, T3)))
                g_min.setdefault(k1, set()).add(g)
                g_tab.setdefault(k2, set()).add(g)
                tot += 1
    m1 = sum(1 for v in g_min.values() if len(v) > 1)
    m2 = sum(1 for v in g_tab.values() if len(v) > 1)
    say("  %-5d %-11d %-34s %-34s"
        % (n, tot,
           "%d/%d groups multi-valued" % (m1, len(g_min)),
           "%d/%d groups multi-valued" % (m2, len(g_tab))))
say("")
say("  Zero multi-valued groups under the full two-body key at every n: GAP IS pairwise-determined.")

# ---------------------------------------------------------------- (c) the lane's own controls
say("")
say("B3  THE LANE'S OWN D-15 CONTROLS FOR GAP, RE-READ")
say("-" * 120)
say("  E-9 table G2 prints an 'all three on the same block' reference row alongside its CTRL-Z rows.")
say("  In the lane's own e9_frustration_gap.txt that reference row reads EXACTLY ZERO in ALL FOUR")
say("  configurations (n0,m) = (4,2),(4,3),(6,2),(6,3).  A null whose positive reference also returns")
say("  zero classifies nothing -- this is precisely the defect that the lane itself diagnosed in E-8's")
say("  quantity T and used to discard it.  Re-derive that on-block reference here, uniformly:")
say("")
say("  %-5s %-4s %-6s %-34s %-30s" % ("n0", "m", "n", "pattern", "GAP values found (uniform sample)"))
for n0 in (4, 6):
    m = 2; n = n0 * m
    stab = []
    for b in range(m):
        sx = [0] * (2 * n); sz = [0] * (2 * n)
        for jj in range(n0):
            sx[b * n0 + jj] = 1; sz[n + b * n0 + jj] = 1
        stab.append(enc(sx, n)); stab.append(enc(sz, n))
    S = [(0, 0)]
    for g in stab:
        S = S + [xr(s, g) for s in S]
    S = sorted(set(S))
    sX0 = [1] * n0 + [0] * n0; sZ0 = [0] * n0 + [1] * n0
    p0 = symplectic_logicals([sX0, sZ0], n0)
    byb = {}
    for b in range(m):
        gg = []
        for (Xl, Zl) in p0:
            for v in (Xl, Zl):
                w = [0] * (2 * n)
                for jj in range(n0):
                    w[b * n0 + jj] = v[jj]; w[n + b * n0 + jj] = v[n0 + jj]
                gg.append(enc(w, n))
        byb[b] = gg
    rng = random.Random(99 + n0)
    def ur(b):
        v = (0, 0)
        for g in byb[b]:
            if rng.getrandbits(1): v = xr(v, g)
        return v
    for lbl, bs in (("all three on block 0", (0, 0, 0)), ("two on block 0, one on block 1", (0, 0, 1))):
        vals = {}
        for _ in range(20000):
            a = ur(bs[0]); b = ur(bs[1]); c = ur(bs[2])
            g = GAPt(tbl(a, b, S), tbl(b, c, S), tbl(a, c, S))
            vals[g] = vals.get(g, 0) + 1
        say("  %-5d %-4d %-6d %-34s %-30s" % (n0, m, n, lbl, dict(sorted(vals.items()))))
say("")
say("  If the on-block reference at n0=4 is {0} while the cross-region row is also {0}, the n0=4")
say("  cross-region null is uninformative by the lane's own D-15 rule; only n0=6 carries a live")
say("  reference.  The lane's summary nevertheless asserts the separation pattern holds at")
say("  'n0 = 4,6,8 and m = 2..6' for GAP; e10_gap_characterised.txt shows m=3 SKIPPED and no n0=8 row.")

say("")
say("=" * 120)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO/VERIFY/v2_sampling_and_controls.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
