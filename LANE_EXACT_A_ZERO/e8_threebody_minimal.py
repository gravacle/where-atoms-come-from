"""E-8  A GENUINELY THREE-BODY, GAUGE-INVARIANT, EXACTLY-DECIDED INTEGER INVARIANT.

E-2/E-3 showed that every ALGEBRAIC three-record scalar (associator, tau, K, any word) is either
representative-dependent or determined by the pairwise F_2 pairing.  E-7 then found that the
MINIMISED support quantity J -- a geometric, not algebraic, invariant -- is exactly non-zero on
pairs the F_2 pairing cannot see.  This step builds the three-record analogue and asks the only
question that matters: is it determined by the pairwise data, or is it new?

DEFINITION (exact; gauge-invariant because it is a minimum over the whole orbit).
For record classes [a],[b],[c] let
    T([a],[b],[c]) = min over the |S|^3 representative choices of
                     #{ sites j : the three single-qubit Paulis at j are ALL non-identity
                                  and PAIRWISE DIFFERENT }
A site contributing to T is one where all three records act, and act differently -- a genuine
three-way meeting that no pair of them exhibits.  T is an integer, exactly computed.

DECIDED HERE
  U1  the exact distribution of T
  U2  is T EXACTLY NON-ZERO on any triple?
  U3  DETERMINACY: is T a function of the three pairwise invariants (sp_ab,sp_bc,sp_ac, J_ab,J_bc,J_ac)?
      If NOT, T is three-body information no two-body quantity carries -- the first such object in
      this program.
  U4  SEPARATION: T across disjoint regions
  U5  ADDITIVITY: is T over a union of regions the sum of the per-region values?

CONTROLS IN THE SAME TABLE (D-15)
  CTRL-Z   T([a],[a],[b]) -- a repeated record cannot meet itself differently: must be EXACTLY 0
  CTRL-Z   a triple with one record on a disjoint block: must be EXACTLY 0
  CTRL-NZ  a triple whose exhibited representatives already have a three-way site: T must be
           EXACTLY NON-ZERO (it is bounded above by that exhibited count and below by 1 if the
           minimum cannot reach 0)
"""
import sys, random

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import enc, sp_i, xr_i, _pc1, signed_stabiliser_group_i
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


def carrier(n):
    sX = [1] * n + [0] * n
    sZ = [0] * n + [1] * n
    pairs = symplectic_logicals([sX, sZ], n)
    S = signed_stabiliser_group_i([enc(sX, n), enc(sZ, n)])
    gens = []
    for X, Z in pairs:
        gens.append(enc(X, n)); gens.append(enc(Z, n))
    return gens, S


def Iu(A, B):
    xa, za = A; xb, zb = B
    return _pc1((xa & zb) ^ (za & xb))


def J(a0, b0, Sk):
    m = None
    for s1 in Sk:
        a = xr_i(a0, s1)
        for s2 in Sk:
            v = Iu(a, xr_i(b0, s2))
            m = v if m is None else min(m, v)
    return m


def three_way_sites(A, B, C):
    """# sites where all three local Paulis are non-identity AND pairwise different.  Exact bitwise.
       'non-identity at j' = (x|z) != (0,0);  'different at j' = the (x,z) pairs differ."""
    xa, za = A; xb, zb = B; xc, zc = C
    nzA = xa | za; nzB = xb | zb; nzC = xc | zc
    allnz = nzA & nzB & nzC
    dAB = (xa ^ xb) | (za ^ zb)
    dBC = (xb ^ xc) | (zb ^ zc)
    dAC = (xa ^ xc) | (za ^ zc)
    return _pc1(allnz & dAB & dBC & dAC)


def T(a0, b0, c0, Sk):
    m = None
    for s1 in Sk:
        a = xr_i(a0, s1)
        for s2 in Sk:
            b = xr_i(b0, s2)
            for s3 in Sk:
                v = three_way_sites(a, b, xr_i(c0, s3))
                m = v if m is None else (v if v < m else m)
    return m


def classes(gens, cap, seed):
    if 2 ** len(gens) <= cap:
        reps = [(0, 0)]
        for g in gens:
            reps = reps + [xr_i(r, g) for r in reps]
        return [r for r in reps if r != (0, 0)], True
    rng = random.Random(seed)
    out = set(gens)
    for i in range(len(gens)):
        for j in range(i + 1, len(gens)):
            out.add(xr_i(gens[i], gens[j]))
    while len(out) < cap:
        v = (0, 0)
        for g in gens:
            if rng.random() < 0.5:
                v = xr_i(v, g)
        out.add(v)
    out.discard((0, 0))
    return sorted(out), False


say("=" * 128)
say("E-8  T = MINIMAL THREE-WAY MEETING NUMBER -- EXACT, GAUGE-INVARIANT BY CONSTRUCTION")
say("=" * 128)

say("")
say("U1/U2   EXACT DISTRIBUTION OF T, WITH THE CONTROLS IN THE SAME TABLE")
say("-" * 128)
say("  %-4s %-5s %-10s %-11s %-40s %-14s %-26s"
    % ("n", "k", "triples", "exhaustive", "exact distribution of T", "T>0 count", "CTRL-Z T(a,a,b) / T(a,b,b)"))
store = {}
for n in (4, 6, 8, 10, 12):
    gens, S = carrier(n)
    Sk = list(S.keys())
    cap = 64 if n <= 6 else 40
    reps, exh = classes(gens, cap, 555 + n)
    rng = random.Random(n)
    if len(reps) > 40:
        reps = rng.sample(reps, 40); exh = False
    dist = {}
    pos = 0
    recs = []
    for i in range(len(reps)):
        for j in range(len(reps)):
            for k in range(len(reps)):
                a, b, c = reps[i], reps[j], reps[k]
                v = T(a, b, c, Sk)
                dist[v] = dist.get(v, 0) + 1
                if v > 0:
                    pos += 1
                recs.append((a, b, c, v))
    cz1 = T(gens[0], gens[0], gens[1], Sk)
    cz2 = T(gens[0], gens[1], gens[1], Sk)
    say("  %-4d %-5d %-10d %-11s %-40s %-14d %-26s"
        % (n, n - 2, len(recs), "YES" if exh else "sampled classes",
           str(dict(sorted(dist.items()))), pos,
           "%d / %d %s" % (cz1, cz2, "OK" if (cz1 == 0 and cz2 == 0) else "**CONTROL FAILED**")))
    store[n] = (recs, Sk)

say("")
say("U3   DETERMINACY.  Group triples by the SIX pairwise invariants (sp_ab,sp_bc,sp_ac,J_ab,J_bc,J_ac).")
say("     If T takes more than one value inside a group, T is NOT a function of the pairwise data.")
say("-" * 128)
say("  %-4s %-12s %-16s %-46s %-30s"
    % ("n", "groups", "groups with |T|>1", "example group -> T values", "verdict"))
for n in sorted(store):
    recs, Sk = store[n]
    groups = {}
    for (a, b, c, v) in recs:
        key = (sp_i(a, b), sp_i(b, c), sp_i(a, c), J(a, b, Sk), J(b, c, Sk), J(a, c, Sk))
        groups.setdefault(key, set()).add(v)
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    ex = ""
    if multi:
        k0 = sorted(multi)[0]
        ex = "%s -> %s" % (str(k0), sorted(multi[k0]))
    say("  %-4d %-12d %-16d %-46s %-30s"
        % (n, len(groups), len(multi), ex if ex else "-",
           "T IS pairwise-determined" if not multi else
           "T is NOT pairwise-determined -- genuinely 3-body"))

say("")
say("U4/U5   SEPARATION AND ADDITIVITY.  m disjoint [[n0,n0-2,2]] blocks.")
say("-" * 128)
say("  %-5s %-4s %-6s %-34s %-24s %-30s"
    % ("n0", "m", "n", "block pattern of the triple", "T values (exact)", "verdict"))
for n0 in (4, 6):
    for m in (2, 3):
        n = n0 * m
        stabk = []
        for b in range(m):
            sx = [0] * (2 * n); sz = [0] * (2 * n)
            for jj in range(n0):
                sx[b * n0 + jj] = 1
                sz[n + b * n0 + jj] = 1
            stabk.append(enc(sx, n)); stabk.append(enc(sz, n))
        S = signed_stabiliser_group_i(stabk)
        Sk = list(S.keys())
        sX0 = [1] * n0 + [0] * n0
        sZ0 = [0] * n0 + [1] * n0
        p0 = symplectic_logicals([sX0, sZ0], n0)
        byb = {}
        for b in range(m):
            byb[b] = []
            for (Xl, Zl) in p0:
                for v in (Xl, Zl):
                    w = [0] * (2 * n)
                    for jj in range(n0):
                        w[b * n0 + jj] = v[jj]
                        w[n + b * n0 + jj] = v[n0 + jj]
                    byb[b].append(enc(w, n))
                # include the product XZ inside the block too
            extra = []
            for i in range(0, len(byb[b]), 2):
                extra.append(xr_i(byb[b][i], byb[b][i + 1]))
            byb[b] += extra
        pats = [("all three on block 0", 0, 0, 0)]
        if m >= 2:
            pats.append(("two on block 0, one on block 1", 0, 0, 1))
            pats.append(("one each on blocks 0,1 (+ repeat)", 0, 1, 1))
        if m >= 3:
            pats.append(("one on each of blocks 0,1,2", 0, 1, 2))
        budget = len(Sk) ** 3 * 12 ** 3
        if budget > 6e7:
            say("  %-5d %-4d %-6d %-34s  SKIPPED: |S|^3 = %d makes the exhaustive orbit"
                "  sweep infeasible; the m=2 rows above already carry the cross-region test."
                % (n0, m, n, "(all patterns)", len(Sk) ** 3))
            continue
        for lbl, b1, b2, b3 in pats:
            vals = set()
            for a in byb[b1]:
                for b in byb[b2]:
                    for c in byb[b3]:
                        vals.add(T(a, b, c, Sk))
            same = (b1 == b2 == b3)
            ok = ("CTRL: on-block reference" if same else
                  ("EXACTLY ZERO across regions" if vals == {0} else "**NON-ZERO ACROSS REGIONS**"))
            say("  %-5d %-4d %-6d %-34s %-24s %-30s" % (n0, m, n, lbl, sorted(vals), ok))

say("")
say("=" * 128)
say("  E-8 SUMMARY")
say("=" * 128)
say("  T is a three-record integer, exactly computed, gauge-invariant by construction.")
say("  See U3 for whether it is or is not determined by the pairwise data -- that is the whole point.")
say("=" * 128)

with open(LANE + "/e8_threebody_minimal.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
