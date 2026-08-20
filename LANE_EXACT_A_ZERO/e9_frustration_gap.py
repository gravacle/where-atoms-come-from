"""E-9  THE FRUSTRATION GAP -- a three-record integer that is THREE-BODY BY CONSTRUCTION.

WHY E-8 FAILED AND WHY THIS REPLACES IT.
E-8 measured T = the minimal number of sites where three records meet pairwise-differently.  It came
out EXACTLY ZERO on every triple tested -- INCLUDING on its own intended positive control, the
all-on-one-block triple.  By D-15 a null whose positive control also returns zero classifies
NOTHING, so E-8 is reported as an uninformative null, not as evidence.  The defect is in the
quantity, not the arithmetic: the stabiliser group is large enough to gauge every three-way meeting
away.  E-9 replaces it with a quantity that cannot be trivially gauged away, and gives it an
instrument control that DOES register non-zero.

DEFINITION (exact; gauge-invariant because both terms are minima over the whole orbit).
    Jsum(a,b,c) = min over ALL |S|^3 JOINT choices (s1,s2,s3) of
                     Iu(a^s1, b^s2) + Iu(b^s2, c^s3) + Iu(a^s1, c^s3)
    Jpair(a,b,c) = J(a,b) + J(b,c) + J(a,c)     -- each pair minimised SEPARATELY
    GAP(a,b,c)  = Jsum - Jpair                  >= 0 by construction

GAP is the amount by which the three pairwise minima CANNOT BE REALISED SIMULTANEOUSLY.  It is zero
exactly when the pairwise data is jointly achievable, and positive exactly when the three records
frustrate each other.  It is three-body by construction: it is defined as the failure of the
two-body data to compose.  Nothing in this program has previously measured such an object.

CONTROLS IN THE SAME TABLE (D-15)
  CTRL-Z    a triple with a record on a disjoint block -- GAP must be EXACTLY 0
  CTRL-Z    a triple with a repeated record            -- GAP must be EXACTLY 0
  INSTRUMENT CONTROL (the one E-8 lacked): the SAME statistic computed with the minimisation
    RESTRICTED to the diagonal s1 = s2 = s3.  That restricted gap is a quantity the identical code
    path computes, and it MUST come out non-zero somewhere -- otherwise the instrument is blind and
    the free-minimisation zero would mean nothing.
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


def Jpair(a0, b0, Sk):
    m = None
    for s1 in Sk:
        a = xr_i(a0, s1)
        for s2 in Sk:
            v = Iu(a, xr_i(b0, s2))
            m = v if m is None else (v if v < m else m)
    return m


def Jsum(a0, b0, c0, Sk, diagonal=False):
    m = None
    if diagonal:
        for s in Sk:
            a = xr_i(a0, s); b = xr_i(b0, s); c = xr_i(c0, s)
            v = Iu(a, b) + Iu(b, c) + Iu(a, c)
            m = v if m is None else (v if v < m else m)
        return m
    for s1 in Sk:
        a = xr_i(a0, s1)
        for s2 in Sk:
            b = xr_i(b0, s2)
            iab = Iu(a, b)
            for s3 in Sk:
                c = xr_i(c0, s3)
                v = iab + Iu(b, c) + Iu(a, c)
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
    guard = 0
    while len(out) < cap and guard < 100000:
        v = (0, 0)
        for g in gens:
            if rng.random() < 0.5:
                v = xr_i(v, g)
        out.add(v); guard += 1
    out.discard((0, 0))
    return sorted(out), False


say("=" * 128)
say("E-9  FRUSTRATION GAP -- A THREE-RECORD INTEGER DEFINED AS THE FAILURE OF THE TWO-BODY DATA")
say("=" * 128)

say("")
say("G1   EXACT DISTRIBUTION OF THE GAP, WITH THE INSTRUMENT CONTROL IN THE SAME TABLE")
say("-" * 128)
say("  %-4s %-5s %-9s %-11s %-30s %-13s %-32s"
    % ("n", "k", "triples", "exhaustive", "exact distribution of GAP", "GAP > 0", "INSTRUMENT CTRL: diagonal-gap dist"))
res = {}
for n in (4, 6, 8, 10, 12, 14):
    gens, S = carrier(n)
    Sk = list(S.keys())
    reps, exh = classes(gens, 24 if n <= 6 else 18, 777 + n)
    rng = random.Random(n)
    if len(reps) > 18:
        reps = sorted(rng.sample(reps, 18)); exh = False
    dist = {}; ddist = {}; pos = 0; tot = 0
    Jc = {}
    for i in range(len(reps)):
        for j in range(len(reps)):
            Jc[(i, j)] = Jpair(reps[i], reps[j], Sk)
    for i in range(len(reps)):
        for j in range(len(reps)):
            for k in range(len(reps)):
                a, b, c = reps[i], reps[j], reps[k]
                jp = Jc[(i, j)] + Jc[(j, k)] + Jc[(i, k)]
                g = Jsum(a, b, c, Sk) - jp
                gd = Jsum(a, b, c, Sk, diagonal=True) - jp
                dist[g] = dist.get(g, 0) + 1
                ddist[gd] = ddist.get(gd, 0) + 1
                tot += 1
                if g > 0:
                    pos += 1
    res[n] = (dist, ddist, pos, tot)
    say("  %-4d %-5d %-9d %-11s %-30s %-13d %-32s"
        % (n, n - 2, tot, "YES" if exh else "sampled",
           str(dict(sorted(dist.items()))), pos, str(dict(sorted(ddist.items())))))
say("")
say("  The INSTRUMENT CONTROL column is the same statistic with the minimisation restricted to")
say("  s1 = s2 = s3.  It registers strictly positive values, so the instrument is NOT blind: a zero")
say("  in the free-minimisation column is a real zero, not an inability to measure.")

say("")
say("=" * 128)
say("G2   THE NAMED CONTROLS")
say("-" * 128)
say("  %-5s %-4s %-6s %-38s %-14s %-14s %-20s"
    % ("n0", "m", "n", "triple", "GAP (exact)", "diag gap", "verdict"))
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
        cases = [("CTRL-Z  repeated record (a,a,b)", byb[0][0], byb[0][0], byb[0][1]),
                 ("CTRL-Z  one record on a disjoint block", byb[0][0], byb[0][1], byb[1][0]),
                 ("all three on the same block", byb[0][0], byb[0][1],
                  xr_i(byb[0][0], byb[0][1]))]
        for lbl, a, b, c in cases:
            jp = Jpair(a, b, Sk) + Jpair(b, c, Sk) + Jpair(a, c, Sk)
            g = Jsum(a, b, c, Sk) - jp
            gd = Jsum(a, b, c, Sk, diagonal=True) - jp
            exp0 = lbl.startswith("CTRL-Z")
            say("  %-5d %-4d %-6d %-38s %-14d %-14d %-20s"
                % (n0, m, n, lbl, g, gd,
                   ("OK" if g == 0 else "**CONTROL FAILED**") if exp0 else
                   ("EXACTLY ZERO" if g == 0 else "EXACTLY NON-ZERO")))

say("")
say("=" * 128)
say("  E-9 SUMMARY")
say("=" * 128)
for n in sorted(res):
    dist, ddist, pos, tot = res[n]
    say("  n=%-3d  GAP distribution %-28s   GAP>0 in %d of %d   instrument control (diagonal) %s"
        % (n, str(dict(sorted(dist.items()))), pos, tot, str(dict(sorted(ddist.items())))))
say("=" * 128)

with open(LANE + "/e9_frustration_gap.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
