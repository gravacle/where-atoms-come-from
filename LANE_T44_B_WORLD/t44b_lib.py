"""T44-B shared library -- exact machinery for the WORLD (earned D = 3) exponent-class lane.

INSTRUMENT LINEAGE (reuse, cited, not re-derived):
  - The class-discriminator windows, ratio-interval arithmetic, certified sqrt brackets,
    exact-rational resolvent, and the D=2/D=1 control kernels are IMPORTED at runtime from
    LANE_T44_A_CORNER/t44a_lib.py (the judge's N1 instruction: reuse the corner instrument;
    read-only import, nothing written outside this lane).
  - This file adds ONLY: the world venue's own 3D lattices (census access geometry: block
    grains with face adjacency, GR1/T42_C/T43_B lineage; periodic torus for the wrap gate),
    exact integer 3D walk counts by brute-force DP (gate reference) and by the AXIS-SPLIT
    identity (gated against DP before any measurement use), exact-rational partial sums
    with proven tails, the exact regularized critical kernel a_M(x), and the certified
    3D tail lemmas (Wallis-monotonicity bracket, max-trinomial telescoping bound,
    per-factor L3 difference bound with edge-region geometric control).

EVERYTHING ON THE MEASUREMENT PATH IS EXACT: python ints and Fractions only.  Square
roots enter only as certified rational brackets via math.isqrt, and only in COMPARISON
blocks (owner-attributed), never in measurements.

OWNERS of standard mathematics used (ATTRIBUTED; every load-bearing inequality is
re-verified computationally in the driver before use):
  - Walk generating functions / lattice Green's functions, radius of convergence of
    positive series: standard potential theory (Spitzer 1964; Lawler).
  - Transience of D=3 (finiteness of the critical walk sum): Polya 1921.
  - G(0) at criticality on the simple-cubic venue: Watson 1939 (the Watson integral),
    1.5163860591519780...; asymptote G(x) -> (3/(2 pi))/|x|: Spitzer 1964 (P26.1),
    Lawler.  Used ONLY as named comparisons AFTER the class is computed.
  - Axis-split identity N_k^{Z^3}(x) = sum_{k1} C(k,k1) N_{k1}^{Z}(a) N_{k-k1}^{Z^2}(b,c)
    and the rotation bijection for N^{Z^2}: standard combinatorics (Feller I); GATED
    against brute-force DP here before use.
  - Wallis-product bracket (m W_m^2 increasing, (m+1/2) W_m^2 decreasing,
    W_m = C(2m,m)/4^m): classical; both induction steps verified exactly in the driver.
  - Max-trinomial bound q_m <= Q3/(m-2) by exact telescoping product: elementary; every
    step (trinomial Pascal, balanced-mode ratio, telescoped sum) gated in the driver.
  - Perron-Frobenius / Gershgorin row-sum bound: standard linear algebra.
  - Ornstein-Zernike subcritical axis rate cosh(kappa) = 1/(2 mu) - 2 for the D=3
    lattice Green's function: standard asymptotics; COMPARISON ONLY.
"""
import sys
from fractions import Fraction
from math import comb, isqrt

T44A = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T44_A_CORNER"
if T44A not in sys.path:
    sys.path.insert(0, T44A)

# corner-lane instrument reuse (read-only import; judge's N1)
from t44a_lib import (ratio_interval, sqrt_bracket, resolvent_exact,          # noqa: F401
                      potential_kernel_2d, potential_kernel_1d,
                      partial_return_sum_2d, series_target_2d, series_target_1d,
                      cycle_adjacency)


# ---------------------------------------------------------------- venue's own 3D lattices
def torus3_adjacency(n):
    """Periodic n x n x n grain lattice (the census access geometry closed up for the
    universal-cover wrap gate): cells share a face <-> adjacent.  Multiplicity kept."""
    cells = [(x, y, z) for z in range(n) for y in range(n) for x in range(n)]
    idx = {c: i for i, c in enumerate(cells)}
    adj = [dict() for _ in cells]
    for i, (x, y, z) in enumerate(cells):
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            j = idx[((x + dx) % n, (y + dy) % n, (z + dz) % n)]
            adj[i][j] = adj[i].get(j, 0) + 1
    return cells, idx, [sorted(d.items()) for d in adj]


def walk_counts_adj(adj, src, K):
    """Exact integer walk counts N_k(src -> v), k = 0..K, transfer-matrix action."""
    vec = [0] * len(adj)
    vec[src] = 1
    out = [list(vec)]
    for _ in range(K):
        new = [0] * len(adj)
        for i, vi in enumerate(vec):
            if vi:
                for j, mult in adj[i]:
                    new[j] += vi * mult
        vec = new
        out.append(list(vec))
    return out


def bfs_adj(adj, src):
    from collections import deque
    dist = {src: 0}
    dq = deque([src])
    while dq:
        i = dq.popleft()
        for j, _ in adj[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return dist


def dp3_counts(K, inside=None):
    """Brute-force DP walk counts on Z^3 (or restricted to the predicate `inside`):
    returns dict {(k, x, y, z): N}.  Gate reference; exact ints."""
    grid = {(0, 0, 0): 1}
    out = {}
    for k in range(K + 1):
        for v, c in grid.items():
            out[(k,) + v] = c
        new = {}
        for (x, y, z), c in grid.items():
            for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                key = (x + d[0], y + d[1], z + d[2])
                if inside is None or inside(key):
                    new[key] = new.get(key, 0) + c
        grid = new
    return out


def W1(k, a):
    """Exact N_k^{Z}(0 -> a)."""
    if (k + a) % 2 or abs(a) > k:
        return 0
    return comb(k, (k + a) // 2)


def N3_ref(k, a, b, c):
    """Reference axis-split evaluation of N_k^{Z^3}((a,b,c)) (gated against DP)."""
    u, v = b + c, b - c
    tot = 0
    for k1 in range(abs(a), k + 1):
        w1 = W1(k1, a)
        if w1:
            k23 = k - k1
            w2, w3 = W1(k23, u), W1(k23, v)
            if w2 and w3:
                tot += comb(k, k1) * w1 * w2 * w3
    return tot


# ---------------------------------------------------------------- fast exact target counts
def n3_even_row(m, a, b, c):
    """Exact N_{2m}^{Z^3}((a,b,c)) for an even-split target (a even, b+c even) by the
    axis-split identity with ONE incremental running product over j1 = k1/2:
      T(j1) = C(2m,2j1) C(2j1,j1+s) C(2r,r+t1) C(2r,r+t2),  r = m-j1,
      s = a/2, t1 = (b+c)/2, t2 = (b-c)/2,
      T(j1+1)/T(j1) = (r^2-t1^2)(r^2-t2^2) / ((j1+1+s)(j1+1-s)(2r)(2r-1)).
    Every update is integer multiply + EXACT integer divide (asserted)."""
    assert a % 2 == 0 and (b + c) % 2 == 0
    s, t1, t2 = a // 2, (b + c) // 2, (b - c) // 2
    s, t1, t2 = abs(s), abs(t1), abs(t2)
    tmax = max(t1, t2)
    if s > m or tmax > m - s:
        return 0
    j1 = s
    r = m - j1
    T = comb(2 * m, 2 * j1) * comb(2 * r, r + t1) * comb(2 * r, r + t2)  # C(2s,2s)=1
    tot = T
    while j1 + 1 <= m - tmax:
        num = (r * r - t1 * t1) * (r * r - t2 * t2)
        den = (j1 + 1 + s) * (j1 + 1 - s) * (2 * r) * (2 * r - 1)
        q, rem = divmod(T * num, den)
        assert rem == 0
        T = q
        j1 += 1
        r -= 1
        tot += T
    return tot


def series_3d(mu, target, K):
    """Exact partial sum S_K = sum_{k<=K} N_k^{Z^3}(target) mu^k (Fraction) plus the exact
    geometric tail bound (6 mu)^{K+1}/(1-6 mu) (from N_k <= 6^k, gated in driver).
    General-parity target; per-k inner sum over k1 (steps on axis 1) by one incremental
    running product with exact integer divides (asserted):
      Term(k1) = C(k,k1) C(k1,(k1+a)/2) C(k23,(k23+u)/2) C(k23,(k23+v)/2), k23 = k-k1,
      Term(k1+2)/Term(k1) = ((k23^2-u^2)/4)((k23^2-v^2)/4)
                            / ( ((k1+a)/2+1) ((k1-a)/2+1) k23 (k23-1) )."""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    a, b, c = (abs(t) for t in target)
    u, v = b + c, abs(b - c)
    d = a + b + c
    num = 0
    ppow = 1  # p^k
    for k in range(K + 1):
        if k:
            ppow *= p
        Nk = 0
        if k >= d and (k - d) % 2 == 0:
            k1, k23 = a, k - a
            T = comb(k, k1) * W1(k23, u) * W1(k23, v)  # C(2a?,..)=C(a,a)=1 at k1=a
            Nk = T
            while k1 + 2 <= k - max(u, v):
                nmr = (((k23 + u) // 2) * ((k23 - u) // 2)
                       * ((k23 + v) // 2) * ((k23 - v) // 2))
                dnm = ((k1 + a) // 2 + 1) * ((k1 - a) // 2 + 1) * k23 * (k23 - 1)
                qq, rem = divmod(T * nmr, dnm)
                assert rem == 0
                T = qq
                k1 += 2
                k23 -= 2
                Nk += T
        num = num * q + Nk * ppow
    S = Fraction(num, q ** K)
    six_mu = 6 * mu
    assert six_mu < 1
    tail = (six_mu ** (K + 1)) / (1 - six_mu)
    return S, tail


def crit_kernel_3d(targets, M, gate_range=None, gate_bound=None):
    """Regularized critical kernel on Z^3 at the COMPUTED mu_c = 1/6, even-split targets:
        a_M(x) = sum_{m<=M} (N_{2m}(0) - N_{2m}(x)) / 36^m   (EXACT rational; odd-k walk
                                                              counts vanish for these
                                                              parities -- gated)
    Also returns S_M(0) = sum_{m<=M} N_{2m}(0)/36^m (exact) and p_{2M}(0) (exact).
    If gate_range=(mlo,mhi) and gate_bound(m, target)->Fraction are given, verifies the
    assembled difference bound c_{2m}(x) <= p_{2m}(0)*gate_bound EXACTLY on that range and
    returns the boolean conjunction as third output."""
    acc = {t: 0 for t in targets}
    acc0 = 0
    gate_ok = True
    p2M0 = None
    for m in range(M + 1):
        N0 = n3_even_row(m, 0, 0, 0)
        acc0 = acc0 * 36 + N0
        for t in targets:
            Nt = n3_even_row(m, *t)
            assert Nt <= N0  # L1-3D (per-factor central max), enforced en route
            acc[t] = acc[t] * 36 + (N0 - Nt)
            if gate_range and gate_range[0] <= m <= gate_range[1]:
                lhs = Fraction(N0 - Nt, 36 ** m)
                rhs = Fraction(N0, 36 ** m) * gate_bound(m, t)
                if lhs > rhs:
                    gate_ok = False
        if m == M:
            p2M0 = Fraction(N0, 36 ** m)
    den = 36 ** M
    out = {t: Fraction(acc[t], den) for t in targets}
    return out, Fraction(acc0, den), p2M0, gate_ok


# ---------------------------------------------------------------- certified tail machinery
def wallis_brackets(M0):
    """Rational c_l = M0*W_{M0}^2 and c_u = (M0+1/2)*W_{M0}^2, W_m = C(2m,m)/4^m, so that
    (gated monotonicities)  c_l/m <= W_m^2 <= c_u/(m+1/2) for all m >= M0."""
    W = Fraction(comb(2 * M0, M0), 4 ** M0)
    return M0 * W * W, (Fraction(2 * M0 + 1, 2)) * W * W


def max_trinomial(m):
    """Exact max multinomial coefficient m!/(j1! j2! j3!) (scan near the balanced split)."""
    best = 0
    lo = max(0, m // 3 - 2)
    for j1 in range(lo, m // 3 + 3):
        for j2 in range(max(0, (m - j1) // 2 - 2), (m - j1) // 2 + 3):
            j3 = m - j1 - j2
            if j3 < 0:
                continue
            v = comb(m, j1) * comb(m - j1, j2)
            best = max(best, v)
    return best


def q3_constant(m0):
    """Q3 rational with q_m := c_max(m)/3^m <= Q3/(m-2) for all m >= m0+2 (driver gates the
    telescoping steps).  Q3 = m0 * q_{m0} / (1 - 2/(3 m0)), m0 a multiple of 3."""
    assert m0 % 3 == 0
    q = Fraction(max_trinomial(m0), 3 ** m0)
    return (m0 * q) / (1 - Fraction(2, 3 * m0))


def sum_m52_bound(M):
    """Certified rational upper bound for sum_{m>M} (m-2)^{-5/2}: <= (2/3)(M-2)^{-3/2}
    <= 2/(3 (M-2) isqrt(M-2)).  (Telescoping lemma gated in driver: 4m^5 > (m-1)^3(2m+3)^2.)"""
    Mm = M - 2
    return Fraction(2, 3 * Mm * isqrt(Mm))


def sum_m32_bound(M):
    """Certified rational upper bound for sum_{m>M} (m-2)^{-3/2}: <= 2 (M-2)^{-1/2}
    <= 2/isqrt(M-2).  (Telescoping lemma gated in driver: 4m^3 > (m-1)(2m+1)^2.)"""
    return Fraction(2, isqrt(M - 2))


def diff_tail_bound(target, M, B5, RHO, EDGE_C):
    """Certified tail sum_{m>M} (N_{2m}(0)-N_{2m}(x))/36^m
       <= (4 s^2 + 4 t1^2 + 4 t2^2) * B5 * sum_{m>M}(m-2)^{-5/2}
          + EDGE_C * B5 * (M-2)^{-3/2} * RHO^{M+1}/(1-RHO)
    with B5 = Q3 * sqrt_hi(c_u) (driver-assembled), RHO the gated rational with
    RHO^12 >= 2/3.  All rational."""
    a, b, c = target
    s, t1, t2 = abs(a) // 2, abs(b + c) // 2, abs(b - c) // 2
    d2 = 4 * (s * s + t1 * t1 + t2 * t2)
    main = d2 * B5 * sum_m52_bound(M)
    Mm = M - 2
    edge = EDGE_C * B5 * Fraction(1, Mm * isqrt(Mm)) * (RHO ** (M + 1)) / (1 - RHO)
    return main + edge


def abs_tail_bound(M, B5, p2M0):
    """Certified tail sum_{k>2M} p_k(x) for ANY x:
       <= 2 * B5 * sum_{m>M}(m-2)^{-3/2} + p_{2M}(0)
    (even part: p_{2m}(x) <= p_{2m}(0) <= B5 (m-2)^{-3/2}; odd: p_{2m+1}(x) <= p_{2m}(0))."""
    return 2 * B5 * sum_m32_bound(M) + p2M0
