"""ADVERSARIAL RERUN (independent verifier, default refuted) for LANE_T44_A_CORNER.

Every check below is REIMPLEMENTED from scratch (math.comb + Fraction only; the lane's
t44a_lib is NOT imported on any measurement path).  The only shared code is the carrier
definition (Torus) from the O-54-C instrument lineage, used solely to rebuild the (3,7)
coset histogram check R — the venue itself is common ground, not the lane's claim.

Attack surface:
  (1) mu_c measured-vs-imported: independent radius sandwich; independent resolvent
      singularity; sharpness probe (mu = 255/1024 converges-type, mu = 257/1024 has
      growing terms) — a shifted mu_c is excluded within +-1/500;
      class-at-shifted-mu probe: the log kernel signature DISAPPEARS at subcritical mu.
  (2) truncation honesty: tail-bound lemma chain re-checked beyond the lane's range
      (k <= 3000, d = 16); independent kernel brackets at K = 24000 must reproduce the
      lane's intervals; independent K = 4000 bracket must contain the owner values.
  (3) earned separation: independent BFS + shortest-path DP on my own torus grid must
      reproduce every (d, N_min) row of Section 0.
  (4) attribution anchors: computed brackets must contain 4/pi, 4 - 8/pi, (2/pi)ln2,
      3 - 2sqrt(2), 2 - sqrt(3) — all via certified rational brackets built here.
  (5) READ-vs-table: independent numbers must match the OUT table digits.
"""
import sys
from fractions import Fraction
from math import comb, isqrt

F = Fraction
PASS = []


def gate(name, ok):
    PASS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)
    return ok


def ff(x, nd=6):
    x = F(x); s = "-" if x < 0 else ""; x = abs(x)
    sc = (x.numerator * 10 ** nd) // x.denominator
    return "%s%d.%0*d" % (s, sc // 10 ** nd, nd, sc % 10 ** nd)


# certified constant brackets (verifier's own; display-independent)
PI_LO, PI_HI = F(3141592653589793, 10**15), F(3141592653589794, 10**15)
LN2_LO, LN2_HI = F(693147180559945, 10**15), F(693147180559946, 10**15)
def sqrt_br(fr, prec=10**12):
    fr = F(fr); s = fr.numerator * prec * prec // fr.denominator
    r = isqrt(s)
    return F(r, prec), F(r + 2, prec)


# my own Z^2 walk count (independent formula use, gated against my own DP below)
def n2(k, a, b):
    u, v = a + b, a - b
    if (k + u) % 2 or abs(u) > k or abs(v) > k:
        return 0
    return comb(k, (k + u) // 2) * comb(k, (k + v) // 2)


# ---- A: formula vs my own DP -------------------------------------------------------
grid = {(0, 0): 1}
okA = True
for k in range(13):
    for (a, b), c in grid.items():
        okA &= (n2(k, a, b) == c)
    new = {}
    for (a, b), c in grid.items():
        for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            new[(a + da, b + db)] = new.get((a + da, b + db), 0) + c
    grid = new
gate("A  Z^2 rotation-bijection formula == my own DP, k<=12, every endpoint", okA)

# ---- B: mu_c = 1/4 by independent sandwich + sharpness probe -----------------------
okBlo = all(comb(2*m, m) * (2*m + 1) >= 4**m for m in range(1, 401))
okBhi = all(comb(2*m, m)**2 * (2*m + 1) <= 16**m for m in range(1, 401))
okBind1 = all((2*m+1)*(2*m+3) <= (2*m+2)**2 for m in range(1, 100001))
okBind2 = all((2*m+1)**2 * (m+1) >= 4*m*(m+1)**2 for m in range(1, 100001))
gate("B1 radius sandwich 16^m/(2m+1)^2 <= N_2m(0,0) <= 16^m: direct m<=400 AND both "
     "induction steps to m<=10^5 => radius of the return series is EXACTLY 1/4", okBlo and okBhi and okBind1 and okBind2)
# sharpness: mu = 257/1024 > 1/4 has growing terms; mu = 255/1024 < 1/4 has geometric tail
mu_hi = F(257, 1024)
m = 4096
t_m  = F(n2(2*m, 0, 0)) * mu_hi**(2*m)
t_m1 = F(n2(2*m+2, 0, 0)) * mu_hi**(2*m+2)
mu_lo = F(255, 1024)
okB2 = (t_m1 > t_m) and (t_m > 10**3) and (4*mu_lo < 1)
gate("B2 sharpness within +-1/500 of 1/4: at mu=257/1024 the return terms GROW "
     "(t_%d=%s < t_%d=%s...) and at mu=255/1024 the geometric tail (4mu)^k contracts "
     "=> a shifted mu_c is excluded; criticality sits AT the computed 1/4" % (2*m, ff(t_m, 1), 2*m+2, ff(t_m1, 1)), okB2)

# ---- C: independent (4,6) dual venue: adjacency, resolvent pole, d + N_min ---------
Lx, Ly = 4, 6
cells = [(x, y) for y in range(Ly) for x in range(Lx)]
idx = {c: i for i, c in enumerate(cells)}
def nbrs(x, y):
    return [((x+1) % Lx, y), ((x-1) % Lx, y), (x, (y+1) % Ly), (x, (y-1) % Ly)]
okrow = all(len(nbrs(*c)) == 4 for c in cells)
# BFS distances + shortest-path counts by my own DP
import collections
dist = {(0, 0): 0}
dq = collections.deque([(0, 0)])
while dq:
    c = dq.popleft()
    for w in nbrs(*c):
        if w not in dist:
            dist[w] = dist[c] + 1
            dq.append(w)
cnt = {c: 0 for c in cells}; cnt[(0, 0)] = 1
order = sorted(cells, key=lambda c: dist[c])
for c in order:
    if c == (0, 0):
        continue
    cnt[c] = sum(cnt[w] for w in nbrs(*c) if dist[w] == dist[c] - 1)
SEALED_46 = {(1,0):(1,1),(2,0):(2,2),(0,1):(1,1),(0,2):(2,1),(0,3):(3,2),(1,1):(2,2),
             (2,1):(3,6),(1,2):(3,3),(2,2):(4,12),(1,3):(4,8),(2,3):(5,40)}
okDN = all(dist[v] == d and cnt[v] == N for v, (d, N) in SEALED_46.items())
gate("C1 (4,6) my own venue: 4-regular AND independent BFS d + shortest-path count "
     "reproduce EVERY sealed Section-0 row (d, N_min) incl. (2,3)->(5,40)", okrow and okDN)
# exact Fraction elimination: singular at 1/4, nonsingular at 1/8; constant vector killed
def resolvent(mu):
    n = len(cells)
    M = [[F(0)] * n for _ in range(n)]
    for i, c in enumerate(cells):
        M[i][i] = F(1)
        for w in nbrs(*c):
            M[i][idx[w]] -= F(mu)
    rhs = [F(0)] * n; rhs[0] = F(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]; rhs[col], rhs[piv] = rhs[piv], rhs[col]
        inv = 1 / M[col][col]
        M[col] = [x * inv for x in M[col]]; rhs[col] *= inv
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [xr - f * xc for xr, xc in zip(M[r], M[col])]
                rhs[r] -= f * rhs[col]
    return rhs
okC2 = (resolvent(F(1, 4)) is None) and (resolvent(F(1, 8)) is not None) \
    and all(1 - sum(F(1, 4) for _ in nbrs(*c)) == 0 for c in cells)
gate("C2 (4,6) exact resolvent: SINGULAR at mu=1/4, nonsingular at 1/8, constant vector "
     "annihilated by (I - A/4) — the pole sits AT the computed mu_c (my own elimination)", okC2)

# ---- D: subcritical row mu=1/8 independent; READ-vs-table --------------------------
def G2(mu, a, b, K):
    mu = F(mu)
    S = sum(F(n2(k, a, b)) * mu**k for k in range(K + 1))
    tail = (4*mu)**(K+1) / (1 - 4*mu)
    return S, tail
mu = F(1, 8); K = 260
Gd = {d: G2(mu, d, 0, K) for d in (1, 2, 8, 9)}
r8 = (Gd[9][0] / (Gd[8][0] + Gd[8][1]), (Gd[9][0] + Gd[9][1]) / Gd[8][0])
sl, sh = sqrt_br(2)
owner = (3 - 2*sh, 3 - 2*sl)          # 3 - 2 sqrt2
okD = (ff(Gd[1][0], 8) == "0.14636401") and (ff(Gd[2][0], 8) == "0.02062932") \
    and (ff(r8[0]) == "0.160658") and r8[1] <= 1 - F(1, 20) \
    and abs((r8[0]+r8[1])/2 - (owner[0]+owner[1])/2) <= F(1, 25)
gate("D  mu=1/8 independent series: G(1)=0.14636401, G(2)=0.02062932, r(8)=0.160658 "
     "MATCH the OUT table digit-for-digit; r < 1 - 1/20; within 1/25 of owner 3-2sqrt2 "
     "(READ==table verified on the exponential row)", okD)

# ---- E: tail lemma pushed BEYOND the lane's range ----------------------------------
okE = True
for d in (8, 16):
    for k in range(d, 3001, 2):
        diff = n2(k, 0, 0) - n2(k, d, 0)
        okE &= (diff * 2 * k * k <= (d*d + d*d) * 4**k)
gate("E  critical tail-term bound c_k(d,0) <= (u^2+v^2)/(2k^2) holds to k<=3000 "
     "(beyond the lane's gated k<=2000): the proven tail is not an edge artifact", okE)

# ---- F: critical row — divergence + kernel, INDEPENDENT stepper, K=24000 -----------
class Step:
    """my own incremental C(k,(k+u)/2), exact divisions asserted"""
    def __init__(s, u):
        s.u = abs(u); s.k = s.u; s.v = 1
    def two(s):
        num = s.v * (s.k + 1) * (s.k + 2)
        q, r = divmod(num, ((s.k + s.u)//2 + 1) * ((s.k - s.u)//2 + 1))
        assert r == 0
        s.v = q; s.k += 2
print("  [F] running independent K=24000 kernel (exact ints)...")
targets = [("axis", 2, 2, 2), ("axis", 4, 4, 4), ("axis", 8, 8, 8), ("axis", 16, 16, 16),
           ("diag", 2, 2, 0), ("diag", 4, 4, 0), ("diag", 8, 8, 0), ("diag", 16, 16, 0)]
KC = 24000
acc = {t[:2]: 1 for t in targets}
st = {t[:2]: None for t in targets}
c0 = Step(0)
k = 0
S_at = {}
Swant = {256, 1024, 4096, 16384}
Sacc = 1
while k + 2 <= KC:
    c0.two(); k += 2
    c0sq = c0.v * c0.v
    Sacc = Sacc * 16 + c0sq
    if k in Swant or k + 1 in Swant:
        S_at[k] = F(Sacc, 4**k)
    for kind, lab, u, v in targets:
        key = (kind, lab)
        t = st[key]
        if t is None and max(u, v) <= k:
            su, sv = Step(u), Step(v)
            while su.k < k: su.two()
            while sv.k < k: sv.two()
            st[key] = t = (su, sv)
        elif t is not None:
            t[0].two(); t[1].two()
        nk = t[0].v * t[1].v if t is not None else 0
        acc[key] = acc[key] * 16 + (c0sq - nk)
A = {}
for kind, lab, u, v in targets:
    A[(kind, lab)] = (F(acc[(kind, lab)], 4**KC), F(u*u + v*v, 4*KC))
# divergence: quadrupling increments + harmonic bound (my own floor-scaled harmonic)
hlo = F(sum((1 << 40) // (4*m) for m in range(1, 8193)), 1 << 40)
okF1 = (S_at[16384] >= hlo) and all(S_at[c*4] - S_at[c] >= F(1, 4) for c in (256, 1024, 4096)) \
    and ff(S_at[16384], 4) == "3.9345"
gate("F1 critical G(0,0): independent S_16384 = 3.9345 (table match), >= harmonic bound "
     "%s, quadrupling increments non-shrinking => marginal divergence reproduced" % ff(hlo, 3), okF1)
# owner anchors: 4/pi and 4 - 8/pi inside independent brackets
a11 = A[("diag", 2)]; a20 = A[("axis", 2)]
four_over_pi = (4/PI_HI, 4/PI_LO)
four_m_8pi = (4 - 8/PI_LO, 4 - 8/PI_HI)
okF2 = (a11[0] <= four_over_pi[0] and four_over_pi[1] <= a11[0] + a11[1]) \
    and (a20[0] <= four_m_8pi[0] and four_m_8pi[1] <= a20[0] + a20[1])
gate("F2 owner anchors: independent a(1,1) bracket [%s,%s] contains 4/pi and a(2,0) "
     "bracket [%s,%s] contains 4-8/pi (attribution anchors reproduced)"
     % (ff(a11[0]), ff(a11[0]+a11[1]), ff(a20[0]), ff(a20[0]+a20[1])), okF2)
# lane's sealed intervals must contain/intersect mine (same K, so equality expected)
SEALED_A = {("axis",2): ("1.453467","1.453551"), ("axis",4): ("1.907762","1.908095"),
            ("axis",8): ("2.351479","2.352812"), ("axis",16): ("2.790863","2.796196"),
            ("diag",2): ("1.273213","1.273254"), ("diag",4): ("1.697546","1.697713"),
            ("diag",8): ("2.133767","2.134434"), ("diag",16): ("2.572540","2.575207")}
okF3 = all(ff(A[key][0]) == lo and ff(A[key][0] + A[key][1]) == hi
           for key, (lo, hi) in SEALED_A.items())
gate("F3 all 8 sealed kernel intervals reproduced DIGIT-FOR-DIGIT at K=24000 by the "
     "independent stepper (READ==table on the critical row)", okF3)
# class booleans recomputed
def inc(kind, d2, d1):
    (s2, t2), (s1, t1) = A[(kind, d2)], A[(kind, d1)]
    return (s2 - s1 - t1, s2 + t2 - s1)
IAX = [inc("axis", 4, 2), inc("axis", 8, 4), inc("axis", 16, 8)]
IDG = [inc("diag", 4, 2), inc("diag", 8, 4), inc("diag", 16, 8)]
rats = [(i2[0]/i1[1], i2[1]/i1[0]) for i1, i2 in [(IAX[0],IAX[1]),(IAX[1],IAX[2]),(IDG[0],IDG[1]),(IDG[1],IDG[2])]]
oklog = all(F(4,5) <= lo and hi <= F(5,4) for lo, hi in rats)
oknotlin = all(hi < F(9,5) for lo, hi in rats)
tol = (2*LN2_LO/PI_HI, 2*LN2_HI/PI_LO)
okF4 = oklog and oknotlin and (IAX[2][0] <= tol[0] and tol[1] <= IAX[2][1])
gate("F4 class booleans recomputed: all doubling-increment ratios inside [4/5,5/4] "
     "(LOG), outside linear window; largest increment contains (2/pi)ln2 "
     "=> MARGINAL_LOG is the computed class, not a label", okF4)
# K=4000 independent shorter run: owner value must still be inside (truncation honesty)
acc4 = 1; s4 = None; c4 = Step(0); k = 0
while k + 2 <= 4000:
    c4.two(); k += 2
    if s4 is None and k >= 2:
        su, sv = Step(2), Step(2)
        while su.k < k: su.two()
        while sv.k < k: sv.two()
        s4 = (su, sv)
    elif s4 is not None:
        s4[0].two(); s4[1].two()
    nk = s4[0].v * s4[1].v if s4 is not None else 0
    acc4 = acc4 * 16 + (c4.v * c4.v - nk)
a20_s = F(acc4, 4**4000); a20_t = F(8, 16000)
okF5 = a20_s <= four_m_8pi[0] and four_m_8pi[1] <= a20_s + a20_t
gate("F5 truncation honesty: at the SHORTER K=4000 the proven tail (u^2+v^2)/(4K) still "
     "brackets the owner value 4-8/pi ([%s,%s]): tails are real bounds, not decoration"
     % (ff(a20_s), ff(a20_s + a20_t)), okF5)

# ---- G: shifted-mu probe — the log signature DISAPPEARS off criticality ------------
mu = F(23, 100); K = 2400
Gg = {d: G2(mu, d, 0, K) for d in (2, 4, 8, 16)}
b = {d: Gg[2][0] - Gg[d][0] for d in (4, 8, 16)}   # saturating differences
r_off = (b[16] - b[8]) / (b[8] - b[4])
okG = r_off < F(4, 5) / 4
gate("G  class-at-shifted-mu: at subcritical mu=23/100 the doubling-increment ratio of "
     "the same kernel construction is %s — far BELOW the log window [4/5,5/4]: the log "
     "class exists only AT the computed mu_c; a mislocated mu_c could not fake it" % ff(r_off, 4), okG)

# ---- H: 1D discriminator independent ----------------------------------------------
def G1(mu, d, K):
    mu = F(mu)
    S = sum(F(comb(k, (k+d)//2)) * mu**k for k in range(abs(d), K+1) if (k+d) % 2 == 0)
    return S, (2*mu)**(K+1) / (1 - 2*mu)
mu = F(1, 4)
g = {d: G1(mu, d, 260) for d in range(1, 14)}
R = {d: (g[d+1][0]/(g[d][0]+g[d][1]), (g[d+1][0]+g[d+1][1])/g[d][0]) for d in range(1, 13)}
los = [lo for lo, hi in R.values()]; his = [hi for lo, hi in R.values()]
s3l, s3h = sqrt_br(3)
okH1 = (max(his) - min(los) <= F(1, 100)) and (2 - s3h <= min(los)) and (max(his) <= 2 - s3l + F(1, 10**6))
gate("H1 1D mu=1/4 (2D's critical value): ratio CONSTANT (spread <= 1/100) and equals "
     "2-sqrt3 (root of mu r^2 - r + mu) — pure exponential, SUBCRITICAL here: mu_c is "
     "the venue's own number, independently confirmed", okH1)
# 1D critical kernel at K=20000 my own; a(d) must track d
K1, M1 = 20000, 10000
acc1 = {d: 1 for d in (2, 4, 8)}
st1 = {d: None for d in (2, 4, 8)}
c1 = Step(0); k = 0
while k + 2 <= K1:
    c1.two(); k += 2
    for d in (2, 4, 8):
        t = st1[d]
        if t is None and d <= k:
            t = Step(d)
            while t.k < k: t.two()
            st1[d] = t
        elif t is not None:
            t.two()
        acc1[d] = acc1[d] * 4 + (c1.v - (t.v if t else 0))
A1 = {d: (F(acc1[d], 2**K1), F(d*d, 2*isqrt(2*M1))) for d in (2, 4, 8)}
okH2 = all(A1[d][0] <= d <= A1[d][0] + A1[d][1] for d in (2, 4, 8))
I1 = (A1[4][0] - A1[2][0] - A1[2][1], A1[4][0] + A1[4][1] - A1[2][0])
I2 = (A1[8][0] - A1[4][0] - A1[4][1], A1[8][0] + A1[8][1] - A1[4][0])
r1d = (I2[0]/I1[1], I2[1]/I1[0])
okH3 = F(9,5) <= r1d[0] and r1d[1] <= F(11,5) and r1d[0] > F(5,4)
gate("H2 1D critical kernel independent (K=20000): brackets contain 2, 4, 8 EXACTLY "
     "(owner a_Z(d)=|d|) and doubling ratio [%s,%s] sits in the LINEAR window, disjoint "
     "from the 2D log window => D-15 discriminator reproduced" % (ff(r1d[0],3), ff(r1d[1],3)), okH2 and okH3)

# ---- I: supercritical witness independent ------------------------------------------
mu = F(1, 3); m = 200
t0 = F(n2(2*m, 0, 0)) * mu**(2*m)
rat = F(n2(2*m+2, 0, 0), n2(2*m, 0, 0)) * mu * mu
okI = t0 > 10**6 and rat > 1
gate("I  supercritical mu=1/3: term at m=200 = %s > 10^6, term ratio %s > 1 "
     "(divergence witness reproduced independently)" % (ff(t0, 1), ff(rat, 4)), okI)

# ---- R: (3,7) true-coset histogram, INDEPENDENT enumeration ------------------------
O54C = "/Users/bgm/MB Work/where-atoms-come-from/LANE_O54_C_ATTEMPT"
if O54C not in sys.path:
    sys.path.insert(0, O54C)
try:
    from o54c_lib import Torus  # carrier definition only (shared venue, not lane code)
    T = Torus(3, 7)
    n = T.n
    stars = T.all_stars()
    # my own greedy F2-independent subset
    out, piv = [], {}
    for gph in stars:
        mm = gph
        while mm:
            tt = mm.bit_length() - 1
            if tt in piv:
                mm ^= piv[tt]
            else:
                piv[tt] = mm; out.append(gph); break
    gens = [g & ((1 << n) - 1) for g in out + [T.xbar1(), T.xbar2()]]
    assert len(gens) == 22
    import numpy as np
    arr = np.zeros(1, dtype=np.uint64)
    for gph in gens:
        arr = np.concatenate([arr, arr ^ np.uint64(gph)])
    def hist_for(v):
        rep = T.dual_path_x((0, 0), v) & ((1 << n) - 1)
        x = arr ^ np.uint64(rep)
        # my own popcount (numpy bit tricks, exact)
        w = np.zeros_like(x)
        for shift in range(0, 64, 8):
            byte = (x >> np.uint64(shift)) & np.uint64(0xFF)
            lut = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint64)
            w += lut[byte.astype(np.int64)]
        return np.bincount(w.astype(np.int64), minlength=n + 1)
    # expected values read from the SEALED RESULT.json (6-digit excess, 8-digit G_true)
    checks = {(0, 3): (3, 3, 1, "0.477035", 6, "0.00288483"),
              (1, 3): (4, 4, 4, "0.731146", 7, "0.00169057")}
    okR = True
    for v, (d, wm, Nm, exc_txt, amax13, g18_txt) in checks.items():
        h = hist_for(v)
        wmin = next(w for w in range(len(h)) if h[w])
        mu8 = F(1, 8)
        g18 = sum(F(int(c)) * mu8**w for w, c in enumerate(h) if c)
        excess = g18 / (int(h[wmin]) * mu8**d) - 1
        mu3 = F(1, 3)
        am3 = max(range(len(h)), key=lambda w: F(int(h[w])) * mu3**w if h[w] else F(-1))
        am8 = max(range(len(h)), key=lambda w: F(int(h[w])) * mu8**w if h[w] else F(-1))
        okR &= (wmin == wm == d and int(h[wmin]) == Nm and int(h.sum()) == 1 << 22
                and ff(excess) == exc_txt and ff(g18, 8) == g18_txt
                and am8 == d and am3 == amax13 and h[d + 1] > 0)
    gate("R  (3,7) true-coset histograms REBUILT (my own span walk + byte-LUT popcount): "
         "v=(0,3) and v=(1,3) reproduce (w_min=d, N_min), total 2^22, the sealed G_true "
         "and excess values DIGIT-FOR-DIGIT (0.477035 / 0.731146), argmax w=d at mu=1/8, "
         "argmax 6/7 at mu=1/3, and the wrap stratum at w=d+1 — the corrected S7 claims "
         "hold on the true ensemble", okR)
except Exception as e:
    gate("R  (3,7) true-coset histogram rebuild (EXCEPTION: %r)" % e, False)

print("\nADVERSARIAL RERUN: %d checks, %d pass, %d fail"
      % (len(PASS), sum(PASS), len(PASS) - sum(PASS)))
sys.exit(0 if all(PASS) else 1)
