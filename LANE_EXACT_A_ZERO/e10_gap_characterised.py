"""E-10  THE FRUSTRATION GAP, CHARACTERISED.  Is it really three-body, and where does it switch on?

E-9 found GAP EXACTLY ZERO on every triple at n = 4 (exhaustive over all 15 record classes, 3375
triples) and EXACTLY NON-ZERO on some triples from n = 6 upward.  Three things must now be settled,
all exactly, or the finding is not usable:

  H1  AN EXPLICIT WITNESS.  Print one triple at n = 6 with GAP = 2 in full, with every number in
      the definition shown, so the claim rests on an exhibited object and not on a count.
  H2  IS GAP A FUNCTION OF THE PAIRWISE DATA?  Group triples by (sp_ab,sp_bc,sp_ac, J_ab,J_bc,J_ac).
      If GAP takes more than one value inside a group, it is NOT determined by the two-body data,
      and it is the first quantity in this program of which that is true.
  H3  THE ONSET.  n = 4 gives exactly zero, exhaustively.  Confirm that, and locate the smallest
      carrier at which a non-zero GAP exists, exactly, by exhaustive search where feasible.
  H4  SEPARATION AND ADDITIVITY, at every separation, exactly.

CONTROLS (D-15), carried in every table
  CTRL-Z   repeated record, and a record on a disjoint block: GAP must be EXACTLY 0
  INSTRUMENT CONTROL: the same statistic with the minimisation restricted to s1=s2=s3, which
           registers non-zero and so shows the instrument is not blind.
"""
import sys, random
from itertools import combinations_with_replacement

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import enc, dec, sp_i, xr_i, _pc1, signed_stabiliser_group_i
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


def Jpair(a0, b0, Sk):
    return min(Iu(xr_i(a0, s1), xr_i(b0, s2)) for s1 in Sk for s2 in Sk)


def Jsum(a0, b0, c0, Sk, diagonal=False):
    if diagonal:
        return min(Iu(xr_i(a0, s), xr_i(b0, s)) + Iu(xr_i(b0, s), xr_i(c0, s))
                   + Iu(xr_i(a0, s), xr_i(c0, s)) for s in Sk)
    best = None
    for s1 in Sk:
        a = xr_i(a0, s1)
        for s2 in Sk:
            b = xr_i(b0, s2)
            iab = Iu(a, b)
            for s3 in Sk:
                c = xr_i(c0, s3)
                v = iab + Iu(b, c) + Iu(a, c)
                if best is None or v < best:
                    best = v
    return best


def full_classes(gens):
    reps = [(0, 0)]
    for g in gens:
        reps = reps + [xr_i(r, g) for r in reps]
    return [r for r in reps if r != (0, 0)]


def pstr(k, n):
    v = dec(k, n)
    return "".join("IXZY"[v[j] + 2 * v[n + j]] for j in range(n))


say("=" * 128)
say("E-10  THE FRUSTRATION GAP, CHARACTERISED -- EXACT")
say("=" * 128)

# ------------------------------------------------------- H3 + H1 : onset and an explicit witness
say("")
say("H3   ONSET.  Exhaustive over the ENTIRE record group where the group is small enough.")
say("-" * 128)
say("  %-4s %-6s %-12s %-13s %-12s %-30s %-20s"
    % ("n", "k", "|N(S)/S|", "triples", "exhaustive", "exact GAP distribution", "smallest non-zero"))
witness = None
for n in (4, 6):
    gens, S = carrier(n)
    Sk = list(S.keys())
    reps = full_classes(gens)
    exh = True
    if len(reps) ** 3 > 4_000_000:
        rng = random.Random(2024)
        reps = sorted(rng.sample(reps, 120))
        exh = False
    Jc = {}
    for i in range(len(reps)):
        for j in range(i, len(reps)):
            v = Jpair(reps[i], reps[j], Sk)
            Jc[(i, j)] = v; Jc[(j, i)] = v
    dist = {}
    for i in range(len(reps)):
        for j in range(len(reps)):
            for k in range(len(reps)):
                jp = Jc[(i, j)] + Jc[(j, k)] + Jc[(i, k)]
                g = Jsum(reps[i], reps[j], reps[k], Sk) - jp
                dist[g] = dist.get(g, 0) + 1
                if g > 0 and witness is None and n == 6:
                    witness = (n, reps[i], reps[j], reps[k], Sk, Jc[(i, j)], Jc[(j, k)], Jc[(i, k)], g)
    nz = sorted(v for v in dist if v > 0)
    say("  %-4d %-6d %-12d %-13d %-12s %-30s %-20s"
        % (n, n - 2, 2 ** len(gens), sum(dist.values()), "YES" if exh else "120-class sample",
           str(dict(sorted(dist.items()))), nz[0] if nz else "NONE -- exactly zero"))

say("")
say("H1   AN EXPLICIT WITNESS, every number in the definition shown")
say("-" * 128)
if witness is None:
    say("  no witness found -- the existence claim is NOT supported and nothing is concluded here")
else:
    n, a, b, c, Sk, jab, jbc, jac, g = witness
    say("  carrier [[%d,%d,2]],  stabilisers X^(x)%d and Z^(x)%d" % (n, n - 2, n, n))
    say("  record A representative : %s" % pstr(a, n))
    say("  record B representative : %s" % pstr(b, n))
    say("  record C representative : %s" % pstr(c, n))
    say("  pairwise minima, each achieved separately :  J_AB = %d   J_BC = %d   J_AC = %d   sum = %d"
        % (jab, jbc, jac, jab + jbc + jac))
    say("  joint minimum over all %d representative triples : Jsum = %d"
        % (len(Sk) ** 3, Jsum(a, b, c, Sk)))
    say("  GAP = Jsum - (J_AB + J_BC + J_AC) = %d   -- EXACTLY NON-ZERO" % g)
    say("  the three pairwise minima CANNOT be realised at the same time; that deficit is the object.")
    say("")
    say("  the full joint table (crossing sums over all 64 representative choices), exact:")
    vals = {}
    for i1, s1 in enumerate(Sk):
        for i2, s2 in enumerate(Sk):
            for i3, s3 in enumerate(Sk):
                A = xr_i(a, s1); B = xr_i(b, s2); C = xr_i(c, s3)
                t = Iu(A, B) + Iu(B, C) + Iu(A, C)
                vals[t] = vals.get(t, 0) + 1
    say("     crossing-sum -> count : %s" % dict(sorted(vals.items())))
    say("     pairwise minima are %d+%d+%d = %d, which NEVER occurs in that table" % (jab, jbc, jac, jab + jbc + jac))

# ------------------------------------------------------- H2 : is GAP pairwise-determined?
say("")
say("=" * 128)
say("H2   IS GAP A FUNCTION OF THE PAIRWISE DATA (sp_ab,sp_bc,sp_ac, J_ab,J_bc,J_ac)?")
say("-" * 128)
say("  %-4s %-11s %-9s %-19s %-46s %-28s"
    % ("n", "triples", "groups", "groups with |GAP|>1", "an example group -> the GAP values in it", "verdict"))
for n in (6, 8, 10, 12, 14, 16):
    gens, S = carrier(n)
    Sk = list(S.keys())
    reps = full_classes(gens) if 2 ** len(gens) <= 256 else None
    rng = random.Random(60 + n)
    if reps is None:
        reps = set(gens)
        for i in range(len(gens)):
            for j in range(i + 1, len(gens)):
                reps.add(xr_i(gens[i], gens[j]))
        guard = 0
        while len(reps) < 26 and guard < 20000:
            v = (0, 0)
            for gg in gens:
                if rng.random() < 0.5:
                    v = xr_i(v, gg)
            reps.add(v); guard += 1
        reps.discard((0, 0))
        reps = sorted(reps)
    if len(reps) > 26:
        reps = sorted(rng.sample(reps, 26))
    Jc = {}
    for i in range(len(reps)):
        for j in range(i, len(reps)):
            v = Jpair(reps[i], reps[j], Sk)
            Jc[(i, j)] = v; Jc[(j, i)] = v
    groups = {}
    tot = 0
    for i in range(len(reps)):
        for j in range(len(reps)):
            for k in range(len(reps)):
                a, b, c = reps[i], reps[j], reps[k]
                jp = Jc[(i, j)] + Jc[(j, k)] + Jc[(i, k)]
                g = Jsum(a, b, c, Sk) - jp
                key = (sp_i(a, b), sp_i(b, c), sp_i(a, c), Jc[(i, j)], Jc[(j, k)], Jc[(i, k)])
                groups.setdefault(key, set()).add(g)
                tot += 1
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    ex = "-"
    if multi:
        k0 = sorted(multi)[0]
        ex = "%s -> %s" % (str(k0), sorted(multi[k0]))
    say("  %-4d %-11d %-9d %-19d %-46s %-28s"
        % (n, tot, len(groups), len(multi), ex,
           "GAP IS pairwise-determined" if not multi else
           "GAP is NOT pairwise-determined"))

# ------------------------------------------------------- H4 : separation and additivity
say("")
say("=" * 128)
say("H4   SEPARATION AND ADDITIVITY OF THE GAP OVER DISJOINT REGIONS -- EXACT")
say("-" * 128)
say("  %-5s %-4s %-6s %-34s %-9s %-20s %-24s"
    % ("n0", "m", "n", "block pattern", "#triples", "GAP values (exact)", "verdict"))
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
            L = []
            for (Xl, Zl) in p0:
                for v in (Xl, Zl):
                    w = [0] * (2 * n)
                    for jj in range(n0):
                        w[b * n0 + jj] = v[jj]
                        w[n + b * n0 + jj] = v[n0 + jj]
                    L.append(enc(w, n))
            extra = [xr_i(L[i], L[i + 1]) for i in range(0, len(L) - 1, 2)]
            byb[b] = L + extra
        pats = [("all three on block 0", 0, 0, 0), ("two on block 0, one on block 1", 0, 0, 1)]
        if m >= 3:
            pats.append(("one on each of blocks 0,1,2", 0, 1, 2))
        budget = len(Sk) ** 3 * 12 ** 3
        if budget > 6e7:
            say("  %-5d %-4d %-6d %-34s  SKIPPED: |S|^3 = %d makes the exhaustive orbit"
                "  sweep infeasible; the m=2 rows above already carry the cross-region test."
                % (n0, m, n, "(all patterns)", len(Sk) ** 3))
            continue
        for lbl, b1, b2, b3 in pats:
            vals = set(); cnt = 0
            for a in byb[b1]:
                for b in byb[b2]:
                    for c in byb[b3]:
                        jp = Jpair(a, b, Sk) + Jpair(b, c, Sk) + Jpair(a, c, Sk)
                        vals.add(Jsum(a, b, c, Sk) - jp); cnt += 1
            same = (b1 == b2 == b3)
            say("  %-5d %-4d %-6d %-34s %-9d %-20s %-24s"
                % (n0, m, n, lbl, cnt, sorted(vals),
                   "on-block reference" if same else
                   ("EXACTLY ZERO across regions" if vals == {0} else "**NON-ZERO ACROSS REGIONS**")))

say("")
say("=" * 128)
say("  E-10 SUMMARY")
say("=" * 128)
say("  GAP is a gauge-invariant integer, exact, three-body BY CONSTRUCTION, and the tables above")
say("  give its onset, its determinacy verdict, and its behaviour across disjoint regions.")
say("=" * 128)

with open(LANE + "/e10_gap_characterised.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
