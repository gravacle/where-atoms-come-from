"""ADVERSARIAL RECHECK #2 (independent verifier, default refuted).

Fresh code throughout: own torus dual-lattice construction, own shortest-path DP, own
Fraction Gaussian elimination, own binomial kernel stepper, own coset span enumeration
(numpy XOR-fold + unpackbits popcount).  Only the carrier definition (Torus, xbar) is
taken from the O-54-C instrument lineage, as the lane itself declares.

Attack axes:
 (1) mu_c measured not imported  -- own elimination pole at 1/4; exact binomial sandwich
     to m=3000; sharpness +-3/10000 (terms grow at mu=2503/10000); shifted-mu fake test
     (kernel construction at subcritical mu=249/1000 leaves the log window).
 (2) truncation honesty          -- own K=24000 kernel reproduces sealed brackets
     digit-for-digit; short-K bracket must contain long-K value; owner anchors 4/pi,
     4-8/pi, (2/pi)ln2 from 18-digit rational constants.
 (3) earned dimension = computed -- own BFS/DP on the dual lattice reproduces every sealed
     (d, N_min) row; own 2^22 coset span reproduces (w_min, N_min) on two (3,7) pairs;
     dual lattice built from plaquette SUPPORTS equals the plain torus grid (computed, not
     assumed).
 (4) READ vs table               -- sealed OUT strings parsed and matched against my own
     exact recomputation (subcritical row, critical kernel, 1D kernel).
 (5) controls                    -- 1D discriminator, supercritical witness, S_16384
     divergence witness recomputed.
"""
import re
import sys
from fractions import Fraction
from math import comb, isqrt

import numpy as np

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T44_A_CORNER"
O54C = "/Users/bgm/MB Work/where-atoms-come-from/LANE_O54_C_ATTEMPT"
sys.path.insert(0, O54C)
from o54c_lib import Torus, rank_f2  # carrier lineage only

OUT = open(LANE + "/t44a_corner.OUT.txt").read()
CHECKS = []


def chk(name, ok, extra=""):
    CHECKS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("  " + extra) if extra else ""))


def ffrac(x, nd=6):
    x = Fraction(x)
    sign = "-" if x < 0 else ""
    x = abs(x)
    scaled = (x.numerator * 10 ** nd) // x.denominator
    return "%s%d.%0*d" % (sign, scaled // 10 ** nd, nd, scaled % 10 ** nd)


# 18-digit rational brackets for owner constants (my own, from decimal expansions)
PI_LO, PI_HI = Fraction(314159265358979323, 10**17), Fraction(314159265358979324, 10**17)
LN2_LO, LN2_HI = Fraction(693147180559945309, 10**18), Fraction(693147180559945310, 10**18)
FOUR_OVER_PI = (4 / PI_HI, 4 / PI_LO)
FOUR_MINUS_8_PI = (4 - 8 / PI_LO, 4 - 8 / PI_HI)
TWO_PI_LN2 = (2 * LN2_LO / PI_HI, 2 * LN2_HI / PI_LO)

# ---------------------------------------------------------------- (3) venue + earned d
SEALED_S0 = {
    (4, 6): {(1, 0): (1, 1), (2, 0): (2, 2), (0, 1): (1, 1), (0, 2): (2, 1), (0, 3): (3, 2),
             (1, 1): (2, 2), (2, 1): (3, 6), (1, 2): (3, 3), (2, 2): (4, 12), (1, 3): (4, 8),
             (2, 3): (5, 40)},
    (3, 7): {(1, 0): (1, 1), (0, 1): (1, 1), (0, 2): (2, 1), (0, 3): (3, 1), (1, 1): (2, 2),
             (1, 2): (3, 3), (1, 3): (4, 4)},
}
for (Lx, Ly), rows in SEALED_S0.items():
    T = Torus(Lx, Ly)
    cells = [(x, y) for y in range(Ly) for x in range(Lx)]
    idx = {c: i for i, c in enumerate(cells)}
    supp = {c: frozenset(T.plaq_edges(*c)) for c in cells}
    # my own dual lattice from supports
    adj = {i: {} for i in range(len(cells))}
    for i, ci in enumerate(cells):
        for j in range(i + 1, len(cells)):
            m = len(supp[ci] & supp[cells[j]])
            if m:
                adj[i][j] = m
                adj[j][i] = m
    # plain torus grid, built independently of supports
    grid = {i: {} for i in range(len(cells))}
    for (x, y) in cells:
        i = idx[(x, y)]
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            j = idx[((x + dx) % Lx, (y + dy) % Ly)]
            grid[i][j] = grid[i].get(j, 0) + 1
    chk("V1 (%d,%d): dual lattice FROM SUPPORTS == plain torus grid graph (computed, not assumed; 4-regular)"
        % (Lx, Ly), adj == grid and all(sum(r.values()) == 4 for r in adj.values()))
    # my own BFS + shortest-path-count DP
    import collections
    dist = {idx[(0, 0)]: 0}
    cnt = {idx[(0, 0)]: 1}
    dq = collections.deque([idx[(0, 0)]])
    order = [idx[(0, 0)]]
    while dq:
        i = dq.popleft()
        for j, m in adj[i].items():
            if j not in dist:
                dist[j] = dist[i] + 1
                cnt[j] = 0
                dq.append(j)
    for i in sorted(dist, key=lambda t: dist[t]):
        if dist[i] == 0:
            continue
        cnt[i] = sum(cnt[j] * m for j, m in adj[i].items() if dist[j] == dist[i] - 1)
    ok = all(dist[idx[v]] == d and cnt[idx[v]] == n for v, (d, n) in rows.items())
    chk("V2 (%d,%d): my own BFS distance + shortest-path count reproduce EVERY sealed S0 row (d, N_min)"
        % (Lx, Ly), ok)

# own 2^22 coset span on (3,7): two pairs, own popcount route
T37 = Torus(3, 7)
n37 = T37.n
stars = T37.all_stars()
# greedy F2-independent subset (rewritten)
piv, ind = {}, []
for g in stars:
    m = g
    while m:
        t = m.bit_length() - 1
        if t in piv:
            m ^= piv[t]
        else:
            piv[t] = m
            ind.append(g)
            break
gens = ind + [T37.xbar1(), T37.xbar2()]
assert rank_f2(gens) == len(gens) == 22
mask = (1 << n37) - 1
gens = [g & mask for g in gens]


def own_hist(rep):
    arr = np.zeros(1, dtype=np.uint64)
    for g in gens:
        arr = np.concatenate([arr, arr ^ np.uint64(g)])
    arr ^= np.uint64(rep)
    b = np.unpackbits(arr.view(np.uint8))
    w = b.reshape(-1, 8 * arr.itemsize).sum(axis=1)
    return np.bincount(w.astype(np.int64), minlength=n37 + 1)


for v, (d, nmin) in [((0, 3), (3, 1)), ((1, 3), (4, 4))]:
    rep = T37.dual_path_x((0, 0), v) & mask
    h = own_hist(rep)
    wm = int(np.nonzero(h)[0][0])
    chk("V3 (3,7) v=%s: OWN 2^22 coset span (unpackbits popcount): w_min=%d==d and N_min=%d (total 2^22)"
        % (v, wm, int(h[wm])), wm == d and int(h[wm]) == nmin and int(h.sum()) == 1 << 22)

# ---------------------------------------------------------------- (1) mu_c
# own Fraction Gaussian elimination on the (4,6) dual lattice
T46 = Torus(4, 6)
cells46 = [(x, y) for y in range(6) for x in range(4)]
idx46 = {c: i for i, c in enumerate(cells46)}
N46 = len(cells46)
A46 = [[0] * N46 for _ in range(N46)]
for (x, y) in cells46:
    i = idx46[(x, y)]
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        A46[i][idx46[((x + dx) % 4, (y + dy) % 6)]] += 1


def own_solve(mu):
    mu = Fraction(mu)
    M = [[Fraction(int(i == j)) - mu * A46[i][j] for j in range(N46)] for i in range(N46)]
    rhs = [Fraction(0)] * N46
    rhs[idx46[(0, 0)]] = Fraction(1)
    for c in range(N46):
        p = next((r for r in range(c, N46) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        rhs[c], rhs[p] = rhs[p], rhs[c]
        iv = 1 / M[c][c]
        M[c] = [x * iv for x in M[c]]
        rhs[c] *= iv
        for r in range(N46):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
                rhs[r] -= f * rhs[c]
    return rhs


chk("M1 (4,6): OWN exact elimination: (I-muA) SINGULAR at mu=1/4, nonsingular at 1/8 and 249/1000",
    own_solve(Fraction(1, 4)) is None and own_solve(Fraction(1, 8)) is not None
    and own_solve(Fraction(249, 1000)) is not None)
r18 = own_solve(Fraction(1, 8))

# exact binomial sandwich direct to m=3000
ok = all(16 ** m <= comb(2 * m, m) ** 2 * (2 * m + 1) ** 2 and comb(2 * m, m) ** 2 <= 16 ** m
         for m in range(1, 3001))
chk("M2 sandwich 16^m/(2m+1)^2 <= N_2m(0,0)=C(2m,m)^2 <= 16^m EXACT, m<=3000 (radius exactly 1/4)", ok)
# induction steps for all m (algebraic identities, checked exactly on a long range)
ok = all((2 * m + 1) ** 2 * (m + 1) - 4 * m * (m + 1) ** 2 == m + 1 for m in range(1, 200001))
chk("M3 induction step (2m+1)^2(m+1) - 4m(m+1)^2 == m+1 > 0 EXACT for m<=2*10^5 (lower bound all m)", ok)

# sharpness +-3/10000: at mu=2503/10000 the return terms exceed 10^6 and grow
mu_hi = Fraction(2503, 10000)
m = 20000
c2m = comb(2 * m, m)
t = Fraction(c2m) ** 2 * mu_hi ** (2 * m)
rat = Fraction(comb(2 * m + 2, m + 1)) ** 2 * mu_hi ** 2 / Fraction(c2m) ** 2
chk("M4 sharpness: at mu=1/4 + 3/10000 the return term at m=20000 is > 10^6 (%s) and term ratio > 1 (%s)"
    % (ffrac(t, 1)[:12], ffrac(rat, 6)), t > 10 ** 6 and rat > 1)
mu_lo = Fraction(2497, 10000)
chk("M5 sharpness: at mu=1/4 - 3/10000 the geometric envelope (4mu)^k contracts (4mu=%s < 1)"
    % ffrac(4 * mu_lo, 4), 4 * mu_lo < 1)

# ---------------------------------------------------------------- own kernel machinery
class Step:
    def __init__(self, u):
        self.u = abs(u)
        self.k = self.u
        self.v = 1

    def adv(self):
        k, u = self.k, self.u
        num = self.v * (k + 1) * (k + 2)
        den = ((k + u) // 2 + 1) * ((k - u) // 2 + 1)
        q, r = divmod(num, den)
        assert r == 0
        self.v = q
        self.k += 2
        return q


def own_kernel_2d(targets, K):
    """a_K(u,v) = sum_{k<=K}(N_k(0,0)-N_k(target))/4^k exact; tail (u^2+v^2)/(4K)."""
    acc = {t: 1 for t in targets}
    st = {t: None for t in targets}
    c0 = Step(0)
    k = 0
    while k + 2 <= K:
        c0.adv()
        k += 2
        c2 = c0.v * c0.v
        for (u, v) in targets:
            s = st[(u, v)]
            if s is None and max(abs(u), abs(v)) <= k:
                su, sv = Step(u), Step(v)
                while su.k < k:
                    su.adv()
                while sv.k < k:
                    sv.adv()
                s = (su, sv)
                st[(u, v)] = s
            elif s is not None:
                s[0].adv()
                s[1].adv()
            nk = s[0].v * s[1].v if s else 0
            acc[(u, v)] = acc[(u, v)] * 16 + (c2 - nk)
    den = 4 ** K
    return {t: (Fraction(acc[t], den), Fraction(t[0] ** 2 + t[1] ** 2, 4 * K)) for t in targets}


# axis d -> (u,v)=(d,d); diag c -> (u,v)=(2c,0)
TARGETS = [(2, 2), (4, 4), (8, 8), (16, 16), (2, 0), (4, 0), (8, 0), (16, 0)]
KER24 = own_kernel_2d(TARGETS, 24000)
KER6 = own_kernel_2d(TARGETS, 6000)

# (2) truncation honesty: short-K bracket must contain long-K value (tails are real bounds)
ok = all(KER6[t][0] <= KER24[t][0] + KER24[t][1] and KER24[t][0] <= KER6[t][0] + KER6[t][1]
         for t in TARGETS)
chk("T1 truncation honesty: K=6000 bracket and K=24000 bracket INTERSECT for all 8 targets "
    "(proven tails are real bounds, not decoration)", ok)

# READ vs table: sealed kernel brackets digit-for-digit
SEALED_KER = {}
for line in OUT.splitlines():
    m6 = re.match(r"\s+(axis|diag) (\d+)\s+\(earned d =\s+\d+\):\s+a in \[([0-9.]+), ([0-9.]+)\]", line)
    if m6:
        SEALED_KER[(m6.group(1), int(m6.group(2)))] = (m6.group(3), m6.group(4))
ok = len(SEALED_KER) == 8
for (kind, lab), (slo, shi) in SEALED_KER.items():
    t = (lab, lab) if kind == "axis" else (2 * lab, 0)
    s, tl = KER24[t]
    ok &= (ffrac(s, 6) == slo and ffrac(s + tl, 6) == shi)
chk("T2 READ==table: all 8 sealed critical-kernel brackets reproduced DIGIT-FOR-DIGIT by my own "
    "stepper at K=24000", ok)

# owner anchors from my own rational constants
a11 = KER24[(2, 0)]   # diag 1 target (1,1) -> (u,v)=(2,0): sealed a(1,1)
a20 = KER24[(2, 2)]   # axis 2 target (2,0) -> (u,v)=(2,2): sealed a(2,0)
chk("T3 anchors: a(1,1) bracket contains 4/pi and a(2,0) bracket contains 4-8/pi "
    "(constants from 18-digit rational brackets, my own)",
    a11[0] <= FOUR_OVER_PI[0] and FOUR_OVER_PI[1] <= a11[0] + a11[1]
    and a20[0] <= FOUR_MINUS_8_PI[0] and FOUR_MINUS_8_PI[1] <= a20[0] + a20[1])


def inc(kd2, kd1):
    (s2, t2), (s1, t1) = kd2, kd1
    return (s2 - s1 - t1, s2 + t2 - s1)


IAX = [inc(KER24[(4, 4)], KER24[(2, 2)]), inc(KER24[(8, 8)], KER24[(4, 4)]),
       inc(KER24[(16, 16)], KER24[(8, 8)])]
IDG = [inc(KER24[(4, 0)], KER24[(2, 0)]), inc(KER24[(8, 0)], KER24[(4, 0)]),
       inc(KER24[(16, 0)], KER24[(8, 0)])]
rats = [(b[0] / a[1], b[1] / a[0]) for a, b in zip(IAX[:-1] + IDG[:-1], IAX[1:] + IDG[1:])]
oklog = all(Fraction(4, 5) <= lo and hi <= Fraction(5, 4) for lo, hi in rats)
oklin_excl = all(hi < Fraction(9, 5) for lo, hi in rats)
chk("T4 class recomputed: all doubling-increment ratios inside [4/5,5/4] and outside the linear "
    "window -> MARGINAL_LOG is what the numbers say", oklog and oklin_excl)
chk("T5 largest axis doubling increment contains (2/pi)ln2 (owner value, my own rational bracket)",
    IAX[2][0] <= TWO_PI_LN2[0] and TWO_PI_LN2[1] <= IAX[2][1])

# shifted-mu fake test: same construction at subcritical mu=249/1000 leaves the log window
def sub_kernel(mu, targets, K):
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    acc = {t: 1 for t in targets}
    st = {t: None for t in targets}
    c0 = Step(0)
    k = 0
    ppow = 1
    while k + 2 <= K:
        c0.adv()
        k += 2
        ppow = ppow * p * p
        c2 = c0.v * c0.v
        for t in targets:
            s = st[t]
            if s is None and max(abs(t[0]), abs(t[1])) <= k:
                su, sv = Step(t[0]), Step(t[1])
                while su.k < k:
                    su.adv()
                while sv.k < k:
                    sv.adv()
                s = (su, sv)
                st[t] = s
            elif s is not None:
                s[0].adv()
                s[1].adv()
            nk = s[0].v * s[1].v if s else 0
            acc[t] = acc[t] * q * q + (c2 - nk) * ppow
    den = q ** K
    return {t: Fraction(acc[t], den) for t in targets}


SK = sub_kernel(Fraction(249, 1000), [(4, 4), (8, 8), (16, 16)], 12000)
r_shift = (SK[(16, 16)] - SK[(8, 8)]) / (SK[(8, 8)] - SK[(4, 4)])
chk("M6 shifted-mu fake test: kernel doubling ratio at mu=249/1000 (0.4%% below mu_c) is %s -- "
    "far OUTSIDE the log window [4/5,5/4]: the log class exists only AT the computed mu_c"
    % ffrac(r_shift, 4), r_shift < Fraction(4, 5))

# ---------------------------------------------------------------- subcritical row + resolvent tie
def own_series_2d(mu, u, v, K):
    mu = Fraction(mu)
    su, sv = Step(u), Step(v)
    k0 = max(abs(u), abs(v))
    while su.k < k0:
        su.adv()
    while sv.k < k0:
        sv.adv()
    S = Fraction(su.v * sv.v) * mu ** k0
    k = k0
    while k + 2 <= K:
        su.adv()
        sv.adv()
        k += 2
        S += Fraction(su.v * sv.v) * mu ** k
    tail = (4 * mu) ** (K + 1) / (1 - 4 * mu)
    return S, tail


G18 = {d: own_series_2d(Fraction(1, 8), d, d, 260) for d in range(1, 10)}
sealed_G1 = re.search(r"\n   1    \[(0\.\d{8}), 0\.\d{8}\]   \[(0\.\d{6})",
                      OUT.split("\n  mu = 1/8  (4mu")[1])
r8 = (G18[9][0] / (G18[8][0] + G18[8][1]), (G18[9][0] + G18[9][1]) / G18[8][0])
chk("R1 READ==table (mu=1/8): my G(1)=%s matches sealed %s; r(8)=[%s] matches sealed 0.160658; "
    "ratio < 1 - 1/20" % (ffrac(G18[1][0], 8), sealed_G1.group(1), ffrac(r8[0], 6)),
    ffrac(G18[1][0], 8) == sealed_G1.group(1) == "0.14636401"
    and ffrac(r8[0], 6) == "0.160658" and r8[1] <= Fraction(19, 20))
# owner rate 3 - 2 sqrt2 within 1/25 (comparison)
s2lo = Fraction(isqrt(2 * 10 ** 24), 10 ** 12)
rinf = 3 - 2 * s2lo
chk("R2 (mu=1/8) r(8) within 1/25 of owner OZ rate 3-2sqrt2 = %s (comparison reproduced)"
    % ffrac(rinf, 6), abs(r8[0] - rinf) < Fraction(1, 25))
# my own series+tail brackets my own exact resolvent? (different objects: torus vs Z^2 --
# tie them on the torus itself via my own torus walk counts)
counts = [[0] * N46 for _ in range(41)]
counts[0][idx46[(0, 0)]] = 1
for k in range(1, 41):
    for i in range(N46):
        counts[k][i] = sum(counts[k - 1][j] * A46[j][i] for j in range(N46) if A46[j][i])
mu = Fraction(1, 8)
ok = True
for v in SEALED_S0[(4, 6)]:
    i = idx46[v]
    S = sum(Fraction(counts[k][i]) * mu ** k for k in range(41))
    tl = (4 * mu) ** 41 / (1 - 4 * mu)
    ok &= (S <= r18[i] <= S + tl)
chk("R3 (4,6) mu=1/8: OWN torus series + geometric tail brackets OWN exact resolvent on every "
    "sealed pair (series machinery vs exact linear algebra, both mine)", ok)

# ---------------------------------------------------------------- critical divergence witness
acc = 1
c0 = Step(0)
k = 0
want = {256: None, 1024: None, 4096: None, 16384: None}
while k + 2 <= 16384:
    c0.adv()
    k += 2
    acc = acc * 16 + c0.v * c0.v
    if k in want:
        want[k] = Fraction(acc, 4 ** k)
chk("C1 critical S_K(0,0): my own values %s match sealed 2.6125/3.0524/3.4933/3.9345 and "
    "quadrupling increments >= 1/4"
    % [ffrac(want[K], 4) for K in sorted(want)],
    [ffrac(want[K], 4) for K in sorted(want)] == ["2.6125", "3.0524", "3.4933", "3.9345"]
    and all(want[b] - want[a] >= Fraction(1, 4) for a, b in [(256, 1024), (1024, 4096), (4096, 16384)]))

# ---------------------------------------------------------------- 1D discriminator + super
def own_series_1d(mu, d, K):
    mu = Fraction(mu)
    s = Step(d)
    S = Fraction(s.v) * mu ** s.k
    while s.k + 2 <= K:
        s.adv()
        S += Fraction(s.v) * mu ** s.k
    return S, (2 * mu) ** (K + 1) / (1 - 2 * mu)


G14 = {d: own_series_1d(Fraction(1, 4), d, 260) for d in range(1, 10)}
ratios = [(G14[d + 1][0] / (G14[d][0] + G14[d][1]), (G14[d + 1][0] + G14[d + 1][1]) / G14[d][0])
          for d in range(1, 9)]
spread = max(h for _, h in ratios) - min(l for l, _ in ratios)
s3lo = Fraction(isqrt(3 * 10 ** 24), 10 ** 12)
r_owner = 2 - s3lo   # 2 - sqrt(3): root of r^2/4 - r + 1/4
chk("D1 1D mu=1/4 (2D's critical value): ratio constant (spread %s < 1/100) and equals 2-sqrt3=%s "
    "-> SUBCRITICAL pure exponential here: mu_c is the venue's own number"
    % (ffrac(spread, 8), ffrac(r_owner, 6)),
    spread < Fraction(1, 100) and abs(ratios[-1][0] - r_owner) < Fraction(1, 1000))

# 1D critical kernel, own code, K=20000: a(d)=|d| anchors
def own_kernel_1d(dl, K):
    acc = {d: 1 for d in dl}
    st = {d: None for d in dl}
    c0 = Step(0)
    k = 0
    while k + 2 <= K:
        c0.adv()
        k += 2
        for d in dl:
            s = st[d]
            if s is None and d <= k:
                s = Step(d)
                while s.k < k:
                    s.adv()
                st[d] = s
            elif s is not None:
                s.adv()
            nk = s.v if s else 0
            acc[d] = acc[d] * 4 + (c0.v - nk)
    M = K // 2
    return {d: (Fraction(acc[d], 2 ** K), Fraction(d * d, 2 * isqrt(2 * M))) for d in dl}


K1D = own_kernel_1d([2, 4, 8], 20000)
ok = all(K1D[d][0] <= d <= K1D[d][0] + K1D[d][1] for d in (2, 4, 8))
i1 = inc(K1D[4], K1D[2])
i2 = inc(K1D[8], K1D[4])
r1d = (i2[0] / i1[1], i2[1] / i1[0])
chk("D2 1D critical kernel (own, K=20000): brackets contain 2,4,8 EXACTLY (Spitzer a_Z=|d|); "
    "doubling ratio [%s,%s] in LINEAR window [9/5,11/5], DISJOINT from 2D log window"
    % (ffrac(r1d[0], 3), ffrac(r1d[1], 3)),
    ok and Fraction(9, 5) <= r1d[0] and r1d[1] <= Fraction(11, 5))

# supercritical witness mu=26/100 at m=400 (sealed row)
mu = Fraction(26, 100)
t400 = Fraction(comb(800, 400)) ** 2 * mu ** 800
rat = Fraction(comb(802, 401)) ** 2 * mu ** 2 / Fraction(comb(800, 400)) ** 2
chk("S1 supercritical mu=26/100: my term at m=400 = %s (sealed 33665893174.1) > 10^6, ratio %s > 1"
    % (ffrac(t400, 1), ffrac(rat, 4)),
    ffrac(t400, 1) == "33665893174.1" and rat > 1)

n_ok = sum(CHECKS)
print("\nADVERSARIAL RECHECK #2: %d checks, %d pass, %d fail" % (len(CHECKS), n_ok, len(CHECKS) - n_ok))
