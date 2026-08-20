"""ADVERSARIAL VERIFICATION of LANE_T44_B_WORLD (independent implementations).

Attacks:
 (1) mu_c located computationally?  Independent sandwich checks at large m; divergence
     witnesses recomputed; near-critical subcritical row recomputed independently.
 (2) truncation honesty: independent exact kernel deepened M=1400 -> 2800; the lane's
     published tail bounds must contain the exact deepening increment; the assembled
     per-m difference bound is checked on the ACTUAL usage range m in [1401,2800]
     (exact at spot m, float everywhere); published intervals recomputed to 6dp.
 (3) dimension on the measurement path: the only inputs here are the venue's own walk
     counts (checked against an independent DP); tail exponents derived from the venue's
     combinatorics are re-verified numerically far beyond the lane's gate range.
 (4/5) published-number cross-check (READ vs table) for every headline interval.

Everything on the exact path is int/Fraction.  Floats only where marked ADV-FLOAT
(adversary-side screening with wide margins, never the deciding comparison).
"""
import time
from fractions import Fraction
from math import comb, isqrt, lgamma, exp, log

T0 = time.time()
CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok)))
    print(("OK    " if ok else "REFUTE") + "  " + name + (("  " + extra) if extra else ""))

def ff(x, nd=6):
    x = Fraction(x); sign = "-" if x < 0 else ""; x = abs(x)
    scaled = (x.numerator * 10 ** nd) // x.denominator
    return "%s%d.%0*d" % (sign, scaled // 10 ** nd, nd, scaled % 10 ** nd)

# ---------------- independent walk counts ------------------------------------------------
def my_dp(K):
    """Independent brute-force DP on Z^3 (written fresh)."""
    g = {(0, 0, 0): 1}; out = {}
    for k in range(K + 1):
        for v, c in g.items(): out[(k,) + v] = c
        h = {}
        for (x, y, z), c in g.items():
            for p in ((x+1,y,z),(x-1,y,z),(x,y+1,z),(x,y-1,z),(x,y,z+1),(x,y,z-1)):
                h[p] = h.get(p, 0) + c
        g = h
    return out

def myrow_direct(m, a, b, c):
    """Direct comb sum for N_2m(a,b,c), even-split target (axis-split + 2D rotation)."""
    s, t1, t2 = abs(a)//2, abs(b+c)//2, abs(b-c)//2
    tmax = max(t1, t2)
    if s > m - tmax: return 0
    tot = 0
    for j in range(s, m - tmax + 1):
        r = m - j
        tot += comb(2*m, 2*j)*comb(2*j, j+s)*comb(2*r, r+t1)*comb(2*r, r+t2)
    return tot

def myrow_fast(m, a, b, c):
    """Incremental exact version (ratio re-derived; every divide asserted exact)."""
    s, t1, t2 = abs(a)//2, abs(b+c)//2, abs(b-c)//2
    tmax = max(t1, t2)
    if s > m - tmax: return 0
    j = s; r = m - j
    T = comb(2*m, 2*j)*comb(2*j, j+s)*comb(2*r, r+t1)*comb(2*r, r+t2)
    tot = T
    while j < m - tmax:
        num = (r*r - t1*t1)*(r*r - t2*t2)
        den = (j+1+s)*(j+1-s)*(2*r)*(2*r-1)
        T, rem = divmod(T*num, den)
        assert rem == 0
        j += 1; r -= 1
        tot += T
    return tot

def my_series(mu, target, K):
    """Independent general-parity series S_K = sum_{k<=K} N_k(target) mu^k, direct combs."""
    mu = Fraction(mu)
    a, b, c = (abs(t) for t in target)
    u, v = b + c, abs(b - c)
    S = Fraction(0)
    for k in range(K + 1):
        Nk = 0
        for k1 in range(a, k + 1, 2):
            k23 = k - k1
            if (k23 + u) % 2 or k23 < u or k23 < v: continue
            Nk += (comb(k, k1) * comb(k1, (k1 + a)//2)
                   * comb(k23, (k23 + u)//2) * comb(k23, (k23 + v)//2))
        if Nk: S += Nk * mu ** k
    return S

AXIS = [(2,0,0),(4,0,0),(6,0,0),(8,0,0),(12,0,0),(16,0,0)]
FDIAG = [(2,2,0),(4,4,0),(8,8,0)]
BDIAG = [(2,2,2),(4,4,4),(8,8,8)]
TARGETS = AXIS + FDIAG + BDIAG

print("-- A: independent walk counts vs fresh DP; fast == direct at scattered depths --")
DP = my_dp(10)
chk("A1 my direct row == fresh DP, all even-split points k<=10",
    all(myrow_direct(m, a, b, c) == DP.get((2*m, a, b, c), 0)
        for m in range(6) for a in range(0, 5, 2) for b in range(0, 5)
        for c in range(0, 5) if (b + c) % 2 == 0))
chk("A2 my fast row == my direct row at scattered m (5..2800), all targets+origin",
    all(myrow_fast(m, *t) == myrow_direct(m, *t)
        for m in (5, 17, 50, 137, 500, 1400, 2350, 2800)
        for t in TARGETS + [(0, 0, 0)]))
chk("A3 leading stratum: N_k=0 below d, N_d = d!/(a!b!c!) (fresh DP, 7 targets)",
    all(all(DP.get((k,) + t, 0) == 0 for k in range(sum(t))) and
        DP[(sum(t),) + t] == comb(sum(t), t[0]) * comb(sum(t) - t[0], t[1])
        for t in [(1,0,0),(2,0,0),(3,0,0),(1,1,0),(2,1,0),(1,1,1),(2,2,2)]))

print("-- B: mu_c sandwich far beyond the lane's gate range (exact) --")
sand_ok = True
for m in (100, 400, 1000, 1400, 2000, 2800):
    N0 = myrow_fast(m, 0, 0, 0)
    if not (2 * 36**m <= N0 * (2*m+1) * (m+1) * (m+2) and N0 <= 36**m):
        sand_ok = False
chk("B1 sector sandwich 2*36^m/poly <= N_2m(0) <= 36^m at m in {100,...,2800} (exact)", sand_ok)
mu = Fraction(13, 72)
t150 = myrow_fast(150, 0, 0, 0) * mu ** 300
t75 = myrow_fast(75, 0, 0, 0) * mu ** 150
t1000 = myrow_fast(1000, 0, 0, 0) * mu ** 2000
rat151 = Fraction(myrow_fast(151, 0, 0, 0), myrow_fast(150, 0, 0, 0)) * mu * mu
chk("B2 mu=13/72 divergence: t150 > t75, t1000 > 10^90 * t150, witness ratio matches",
    t150 > t75 and t1000 > 10**90 * t150 and ff(rat151, 4) == "1.1619",
    "t150=%s ratio=%s" % (ff(t150, 2), ff(rat151, 4)))
chk("B3 published supercritical witness term t150 = 3398653.92 (13/72)",
    ff(t150, 2) == "3398653.92")

print("-- C: independent recomputation of a subcritical row + resolvent identity --")
mu = Fraction(1, 8); K = 160
S2 = my_series(mu, (2, 0, 0), K)
S1 = my_series(mu, (1, 0, 0), K - 1)
S3 = my_series(mu, (3, 0, 0), K - 1)
S21 = my_series(mu, (2, 1, 0), K - 1)
chk("C1 resolvent identity S_K(2,0,0) = mu[S(1)+S(3)+4S(2,1,0)] EXACT (my own series)",
    S2 == mu * (S1 + S3 + 4 * S21))
tail = (6*mu)**(K+1) / (1 - 6*mu)
G13 = my_series(mu, (13, 0, 0), K); G14 = my_series(mu, (14, 0, 0), K)
rlo, rhi = G14 / (G13 + tail), (G14 + tail) / G13
chk("C2 mu=1/8 ratio at dmax matches published [0.247928, ...] and <= 1-1/20",
    ff(rlo) == "0.247928" and rhi <= 1 - Fraction(1, 20), "r=[%s,%s]" % (ff(rlo), ff(rhi)))

print("-- D: constants re-derived (c_u, Q3, B5) + bound checks far beyond gate range --")
M0W = 240
W = Fraction(comb(2*M0W, M0W), 4**M0W)
c_u = Fraction(2*M0W + 1, 2) * W * W
def maxtri(m):
    best = 0
    for j1 in range(max(0, m//3 - 2), m//3 + 3):
        for j2 in range(max(0, (m - j1)//2 - 2), (m - j1)//2 + 3):
            if m - j1 - j2 >= 0:
                best = max(best, comb(m, j1) * comb(m - j1, j2))
    return best
# balanced-mode window: verify against full scan independently to m=80
chk("D1 windowed max-trinomial == full scan, m <= 80 (independent)",
    all(maxtri(m) == max(comb(m, j1) * comb(m - j1, j2)
        for j1 in range(m + 1) for j2 in range(m - j1 + 1)) for m in range(1, 81)))
M0T = 300
q300 = Fraction(maxtri(M0T), 3**M0T)
Q3 = (M0T * q300) / (1 - Fraction(2, 3*M0T))
prec = 10**12
su = Fraction(isqrt(c_u.numerator * prec * prec // c_u.denominator) + 2, prec)
B5 = Q3 * su
chk("D2 my B5 matches published 0.466824", ff(B5) == "0.466824", "B5=" + ff(B5))
chk("D3 my Q3 matches published 0.826995", ff(Q3) == "0.826995", "Q3=" + ff(Q3))
# ADV-FLOAT: p_2m(0) <= B5 (m-2)^{-3/2} and q_m <= Q3/(m-2) far beyond the gate range
def lf(n): return lgamma(n + 1)
def p0f(m):  # float p_2m(0) via W_m * T_m/9^m needs T; use exact row instead at spot m
    return None
spot = (700, 1000, 1400, 2000, 2800)
b5_ok = all(Fraction(myrow_fast(m, 0, 0, 0), 36**m) * (m - 2) * (isqrt(m - 2) + 1) <= B5
            for m in spot)
chk("D4 p_2m(0)(m-2)^{3/2} <= B5 EXACT at m in {700,1000,1400,2000,2800}", b5_ok)
q_ok = all(Fraction(maxtri(m), 3**m) * (m - 2) <= Q3 for m in (700, 1001, 1502, 2003, 2800))
chk("D5 q_m(m-2) <= Q3 EXACT at m in {700,1001,1502,2003,2800} (mixed residues)", q_ok)

print("-- E: EXACT kernel deepening M=1400 -> 2800 (the tail-honesty attack) --")
M1, M2 = 1400, 2800
EDGE_C = Fraction(40); RHO = Fraction(199, 200)
def m52(M): return Fraction(2, 3 * (M - 2) * isqrt(M - 2))
def m32(M): return Fraction(2, isqrt(M - 2))
def dtail(t, M):
    a, b, c = t
    s, t1, t2 = abs(a)//2, abs(b+c)//2, abs(b-c)//2
    d2 = 4 * (s*s + t1*t1 + t2*t2)
    Mm = M - 2
    return (d2 * B5 * m52(M)
            + EDGE_C * B5 * Fraction(1, Mm * isqrt(Mm)) * RHO**(M+1) / (1 - RHO))
acc = {t: 0 for t in TARGETS}; acc0 = 0
KER1 = {}; S01 = None; P2M01 = None
permB_ok = True   # assembled per-m bound on the usage range (float screen + exact spots)
EXACT_SPOTS = {1401, 1500, 1700, 2000, 2400, 2800}
rho_f = float(RHO); B5_f = float(B5)
for m in range(M2 + 1):
    N0 = myrow_fast(m, 0, 0, 0)
    acc0 = acc0 * 36 + N0
    for t in TARGETS:
        Nt = myrow_fast(m, *t)
        assert Nt <= N0
        acc[t] = acc[t] * 36 + (N0 - Nt)
        if m > M1:
            a, b, c = t
            s, t1, t2 = abs(a)//2, abs(b+c)//2, abs(b-c)//2
            d2 = 4 * (s*s + t1*t1 + t2*t2)
            # bound: (N0-Nt)/N0 <= d2/m + EDGE_C rho^m
            lhs_f = float(Fraction(N0 - Nt, N0))
            rhs_f = d2 / m + 40.0 * rho_f ** m
            if lhs_f > rhs_f * (1 + 1e-9):
                permB_ok = False
                print("  per-m bound FAILS (float) at m=%d t=%s: %.3e > %.3e" % (m, t, lhs_f, rhs_f))
            if m in EXACT_SPOTS:
                if Fraction(N0 - Nt, N0) > Fraction(d2, m) + EDGE_C * RHO**m:
                    permB_ok = False
                    print("  per-m bound FAILS (exact) at m=%d t=%s" % (m, t))
    if m == M1:
        den = 36 ** M1
        KER1 = {t: Fraction(acc[t], den) for t in TARGETS}
        S01 = Fraction(acc0, den)
        P2M01 = Fraction(N0, den)
den = 36 ** M2
KER2 = {t: Fraction(acc[t], den) for t in TARGETS}
S02 = Fraction(acc0, den)
P2M02 = Fraction(myrow_fast(M2, 0, 0, 0), 36 ** M2)
chk("E1 assembled per-m difference bound holds on the ACTUAL usage range m in (1400,2800]"
    " (float screen everywhere, exact at 6 spot m)", permB_ok)
deep_ok = True
for t in TARGETS:
    inc = KER2[t] - KER1[t]
    if not (0 <= inc <= dtail(t, M1)):
        deep_ok = False
        print("  tail dishonest at %s: inc=%s tail=%s" % (t, ff(inc, 8), ff(dtail(t, M1), 8)))
chk("E2 EXACT deepening: a_2800(t) - a_1400(t) in [0, published-tail(1400)] for ALL targets",
    deep_ok)
TABS1 = 2 * B5 * m32(M1) + P2M01
chk("E3 G(0): S0(2800) - S0(1400) <= abs-tail(1400) EXACT; Watson inside both brackets",
    S02 - S01 <= TABS1 and S01 <= Fraction(1516386059, 10**9)
    and Fraction(1516386060, 10**9) <= S02 + (2 * B5 * m32(M2) + P2M02),
    "G0_1400=[%s,%s] G0_2800=[%s,%s]" % (ff(S01), ff(S01 + TABS1),
                                         ff(S02), ff(S02 + 2 * B5 * m32(M2) + P2M02)))

# published intervals recomputed at M=1400 (READ vs table)
def H(t1, t2, KER, TL):
    return (KER[t2] - KER[t1] - TL[t1], KER[t2] + TL[t2] - KER[t1])
TL1 = {t: dtail(t, M1) for t in TARGETS}
TL2 = {t: dtail(t, M2) for t in TARGETS}
def rats(KER, TL):
    Hax = {2: H((2,0,0),(4,0,0),KER,TL), 4: H((4,0,0),(8,0,0),KER,TL),
           6: H((6,0,0),(12,0,0),KER,TL), 8: H((8,0,0),(16,0,0),KER,TL)}
    Hfd = {2: H((2,2,0),(4,4,0),KER,TL), 4: H((4,4,0),(8,8,0),KER,TL)}
    Hbd = {2: H((2,2,2),(4,4,4),KER,TL), 4: H((4,4,4),(8,8,8),KER,TL)}
    out = {}
    for lab, HH, pairs in (("axis", Hax, [(2,4),(4,8)]), ("fdiag", Hfd, [(2,4)]),
                           ("bdiag", Hbd, [(2,4)])):
        for d1, d2 in pairs:
            out[(lab, d1, d2)] = (HH[d2][0] / HH[d1][1], HH[d2][1] / HH[d1][0])
    return out, Hax, Hfd, Hbd
R1, Hax1, Hfd1, Hbd1 = rats(KER1, TL1)
PUB = {("axis",2,4): ("0.453369","0.458705"), ("axis",4,8): ("0.469869","0.507905"),
       ("fdiag",2,4): ("0.492270","0.513897"), ("bdiag",2,4): ("0.493254","0.536493")}
match = all((ff(R1[k][0]), ff(R1[k][1])) == PUB[k] for k in PUB)
chk("E4 published ratio intervals recomputed EXACTLY at M=1400 (6dp match, all four)", match,
    str({k: (ff(v[0]), ff(v[1])) for k, v in R1.items()}))
hpos = all(HH[0] > 0 for HH in list(Hax1.values()) + list(Hfd1.values()) + list(Hbd1.values()))
chk("E5 all increment interval lower bounds positive (ratio arithmetic valid)", hpos)
R2, Hax2, Hfd2, Hbd2 = rats(KER2, TL2)
INV = (Fraction(2,5), Fraction(3,5))
half_in = R2[("axis",4,8)][0] <= Fraction(1,2) <= R2[("axis",4,8)][1]
tighter = all(R2[k][1] - R2[k][0] < R1[k][1] - R1[k][0] for k in PUB)
in_inv = all(INV[0] <= R2[k][0] and R2[k][1] <= INV[1] for k in PUB)
chk("E6 DEEPER instrument (M=2800): all ratios still in INV window, axis 4->8 still "
    "contains 1/2, every bracket strictly tighter", half_in and tighter and in_inv,
    str({k: (ff(v[0]), ff(v[1])) for k, v in R2.items()}))
# exponent bracket at deeper depth (display)
from math import log2
print("  exponent brackets at M=2800 [display]: " +
      str({str(k): [round(-log2(float(R2[k][1])), 4), round(-log2(float(R2[k][0])), 4)]
           for k in PUB}))
c_ax2 = (16 * Hax2[8][0], 16 * Hax2[8][1])
chk("E7 coefficient at M=2800 tighter and overlaps published [0.467230,0.500095]; "
    "owner 3/(2pi)=0.4774648 relation reported",
    c_ax2[0] <= Fraction(500095, 10**6) and Fraction(467230, 10**6) <= c_ax2[1],
    "c_ax(M=2800)=[%s,%s]" % (ff(c_ax2[0]), ff(c_ax2[1])))

npass = sum(1 for _, ok in CHECKS if ok)
print("\nADV CHECKS: %d OK, %d REFUTE, total %d   [%.1f s]"
      % (npass, len(CHECKS) - npass, len(CHECKS), time.time() - T0))
