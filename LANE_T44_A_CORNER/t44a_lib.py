"""T44-A shared library -- exact machinery for the CORNER (earned D = 2) exponent-class lane.

INSTRUMENT LINEAGE (reuse, cited, not re-derived):
  - Pauli mask conventions, F_2 rank, Torus carrier, generator-graph earned-separation
    instrument, and the exhaustive coset scan are IMPORTED at runtime from
    LANE_O54_C_ATTEMPT/o54c_lib.py (C-78 instrument lineage; O-54-C landscape w_min = d).
  - This file adds ONLY: exact integer walk counting (transfer matrix on the venue's own
    plaquette-adjacency lattice), the Z^2 binomial-product walk-count formula (GATED against
    brute-force DP and against torus counts through the universal-cover wrap identity before
    any measurement use), exact-rational partial sums with PROVEN tail bounds, and a coset
    WEIGHT HISTOGRAM variant of the O-54-C scan.

EVERYTHING ON THE MEASUREMENT PATH IS EXACT: python ints and Fractions only.
No floats anywhere on the measurement path.  Square roots enter only as certified rational
brackets via math.isqrt, and only in COMPARISON blocks (owner-attributed), never in
measurements.

OWNERS of standard mathematics used (ATTRIBUTED; every load-bearing inequality is
re-verified computationally in the driver before use):
  - Walk generating functions / lattice Green's functions, radius of convergence of
    positive series: standard potential theory (Spitzer 1964, "Principles of Random Walk";
    Lawler, "Intersections of Random Walks").
  - Recurrence/marginality of the D<=2 critical walk sum: Polya 1921.
  - Potential kernel of the 2D walk, a(d) = (2/pi) ln d + kappa + O(d^-2):
    Stohr 1950; Spitzer 1964 (P12.3, P28.4); exact-value algorithm McCrea-Whipple 1940.
    Used ONLY as a named comparison AFTER the class is computed.
  - Perron-Frobenius / Gershgorin row-sum bound for the transfer-matrix spectral radius:
    standard linear algebra.
  - Rotation bijection N_k^{Z^2}(a,b) = C(k,(k+a+b)/2) C(k,(k+a-b)/2): standard
    combinatorics (Feller vol. I); GATED against DP here before use.
  - Wallis-product bound (C(2m,m)/4^m)^2 <= 1/(2m+1): classical; induction step verified
    exactly in the driver.
  - Subcritical Ornstein-Zernike axis decay rate cosh(kappa) = 1/(2 mu) - 1 for the D=2
    lattice Green's function: standard asymptotics; COMPARISON ONLY.
"""
import sys
from fractions import Fraction
from math import comb, isqrt

O54C = "/Users/bgm/MB Work/where-atoms-come-from/LANE_O54_C_ATTEMPT"
if O54C not in sys.path:
    sys.path.insert(0, O54C)

import numpy as np
from o54c_lib import (pc, rank_f2, supp_mask, weight_xz, sp_pair, Torus,  # noqa: F401
                      qubit_graph, graph_dist_regions, generator_graph_dist,
                      coset_min_np, _np_pop)


def independent_subset(gens):
    """Greedy F_2-independent subset (same routine as the O-54-C driver)."""
    out, piv = [], {}
    for g in gens:
        m = g
        while m:
            t = m.bit_length() - 1
            if t in piv:
                m ^= piv[t]
            else:
                piv[t] = m
                out.append(g)
                break
    return out


# ---------------------------------------------------------------- venue's own walk lattice
def plaquette_adjacency(T):
    """Adjacency lists of the venue's own dual lattice: plaquettes adjacent iff they SHARE
    A CARRIER EDGE (computed from the plaquette supports alone -- no geometry imported
    beyond the carrier's own labels).  One walk step = one shared carrier edge = one unit
    of writer weight (the per-link Gamma price).  Multiplicity kept exactly (an L=2
    direction would share two edges)."""
    Lx, Ly = T.Lx, T.Ly
    cells = [(x, y) for y in range(Ly) for x in range(Lx)]
    idx = {c: i for i, c in enumerate(cells)}
    supp = {c: set(T.plaq_edges(*c)) for c in cells}
    adj = [[] for _ in cells]
    for i, c in enumerate(cells):
        for j in range(i + 1, len(cells)):
            common = supp[c] & supp[cells[j]]
            if common:
                adj[i].append((j, len(common)))
                adj[j].append((i, len(common)))
    return cells, idx, adj


def walk_counts_torus(adj, src, K):
    """Exact integer walk counts N_k(src -> v), k = 0..K, by transfer-matrix action."""
    n = len(adj)
    vec = [0] * n
    vec[src] = 1
    out = [list(vec)]
    for _ in range(K):
        new = [0] * n
        for i in range(n):
            vi = vec[i]
            if vi:
                for j, mult in adj[i]:
                    new[j] += vi * mult
        vec = new
        out.append(list(vec))
    return out


def bfs_dist(adj, src):
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


def cycle_adjacency(L):
    """1D chain venue: cycle C_L adjacency (the D=1 analog of the dual lattice)."""
    return [[(((i + 1) % L), 1), (((i - 1) % L), 1)] for i in range(L)]


# ---------------------------------------------------------------- exact walk-count formulas
def walk_count_z2(k, a, b):
    """Exact N_k^{Z^2}((0,0)->(a,b)) via the rotation bijection (GATED in driver)."""
    u, v = a + b, a - b
    if (k + u) % 2 or abs(u) > k or abs(v) > k:
        return 0
    return comb(k, (k + u) // 2) * comb(k, (k + v) // 2)


def walk_count_z1(k, d):
    """Exact N_k^{Z}(0->d)."""
    if (k + d) % 2 or abs(d) > k:
        return 0
    return comb(k, (k + d) // 2)


def walk_counts_z2_dp(K):
    """Brute-force DP walk counts on Z^2 for the formula gate: {(k,a,b): N}, k <= K."""
    grid = {(0, 0): 1}
    out = {}
    for k in range(K + 1):
        for ab, c in grid.items():
            out[(k, ab[0], ab[1])] = c
        new = {}
        for (a, b), c in grid.items():
            for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                key = (a + da, b + db)
                new[key] = new.get(key, 0) + c
        grid = new
    return out


class BinStepper:
    """Exact incremental binomial C(k,(k+u)/2) along k = |u|, |u|+2, |u|+4, ...
    Every update is integer multiply + EXACT integer divide (asserted)."""

    def __init__(self, u):
        self.u = abs(u)
        self.k = self.u
        self.val = 1          # C(|u|, |u|) = 1

    def step2(self):
        k, u = self.k, self.u
        num = self.val * (k + 1) * (k + 2)
        d1 = (k + u) // 2 + 1
        d2 = (k - u) // 2 + 1
        q, r = divmod(num, d1 * d2)
        assert r == 0
        self.val = q
        self.k = k + 2
        return q


# ---------------------------------------------------------------- exact partial sums
def series_target_2d(mu, a, b, K):
    """S_K = sum_{k<=K} N_k^{Z^2}((a,b)) mu^k, EXACT Fraction, plus the exact geometric
    tail bound sum_{k>K} N_k mu^k <= (4mu)^{K+1}/(1-4mu) (from N_k(v) <= 4^k = total
    walk count, gated in driver).  u = a+b and v = a-b must share parity (asserted)."""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    u, v = a + b, a - b
    assert (u - v) % 2 == 0
    k0 = max(abs(u), abs(v))
    if (abs(u) - abs(v)) % 2:          # impossible by parity assert, kept for clarity
        raise ValueError
    su, sv = BinStepper(u), BinStepper(v)
    # advance the smaller-|.| stepper to k0
    while su.k < k0:
        su.step2()
    while sv.k < k0:
        sv.step2()
    k = k0
    num = su.val * sv.val * (p ** k)   # scaled by q^0 at exponent-k reference
    scale_k = k
    p2 = p * p
    q2 = q * q
    ppow = p ** k
    while k + 2 <= K:
        su.step2()
        sv.step2()
        k += 2
        ppow *= p2
        num = num * q2 + su.val * sv.val * ppow
        scale_k = k
    S = Fraction(num, q ** scale_k)
    four_mu = 4 * mu
    assert four_mu < 1
    tail = (four_mu ** (K + 1)) / (1 - four_mu)
    return S, tail


def series_target_1d(mu, d, K):
    """Same on the 1D chain; tail from N_k <= 2^k (gated in driver)."""
    mu = Fraction(mu)
    p, q = mu.numerator, mu.denominator
    s = BinStepper(d)
    k = s.k
    num = s.val * (p ** k)
    scale_k = k
    p2, q2 = p * p, q * q
    ppow = p ** k
    while k + 2 <= K:
        s.step2()
        k += 2
        ppow *= p2
        num = num * q2 + s.val * ppow
        scale_k = k
    S = Fraction(num, q ** scale_k)
    two_mu = 2 * mu
    assert two_mu < 1
    tail = (two_mu ** (K + 1)) / (1 - two_mu)
    return S, tail


def series_G_torus(counts, mu, targets, K):
    """Partial sums on the FINITE torus venue from precomputed counts[k][v]; exact
    Fractions + geometric tail (row sums = 4, gated in driver)."""
    mu = Fraction(mu)
    four_mu = 4 * mu
    assert four_mu < 1
    tail = (four_mu ** (K + 1)) / (1 - four_mu)
    out = {}
    for t in targets:
        s = Fraction(0)
        mp = Fraction(1)
        for k in range(K + 1):
            nk = counts[k][t]
            if nk:
                s += nk * mp
            mp *= mu
        out[t] = (s, tail)
    return out


# ---------------------------------------------------------------- critical-row kernels
def partial_return_sum_2d(Klist):
    """S_K(0,0) = sum_{k<=K} N_k(0,0)/4^k at mu_c = 1/4, EXACT Fractions at each K in
    Klist (ascending).  Divergence witness data for the critical row."""
    Kmax = max(Klist)
    want = set(Klist)
    out = {}
    acc = 1                       # k = 0 term, scale 4^0
    c = BinStepper(0)
    k = 0
    if 0 in want:
        out[0] = Fraction(acc, 1)
    while k + 2 <= Kmax:
        c.step2()
        k += 2
        acc = acc * 16 + c.val * c.val
        if k in want or (k + 1) in want:
            out[k] = Fraction(acc, 4 ** k)
    return out


def potential_kernel_2d(dlist_axis, dlist_diag, K):
    """Regularized critical kernel on Z^2 at the COMPUTED mu_c = 1/4, even-parity targets:
        a_K(t) = sum_{k<=K} (N_k(0,0) - N_k(t)) / 4^k      (EXACT rational; k even only,
                                                            odd-k counts vanish -- gated)
    with the PROVEN tail bound sum_{k>K} c_k(t) <= (u^2+v^2)/(4K), u=a+b, v=a-b, derived
    from three elementary lemmas, each GATED in the driver before use:
      L1 nonnegativity: N_k(0,0) >= N_k(t)                    [central-binomial max]
      L2 Wallis:  (C(2m,m)/4^m)^2 <= 1/(2m+1)                 [induction step exact]
      L3 ratio:   C(2m,m+s)/C(2m,m) >= 1 - s^2/m              [telescoping union bound]
    chain: c_{2m}(t) <= p_{2m}(0,0) (s^2+t^2)/m <= (u^2+v^2)/(2(2m)^2);
           sum_{m>K/2} 1/m^2 < 2/K  =>  tail <= (u^2+v^2)/(4K)... (exact algebra:
           sum_{m>M} (u^2+v^2)/(8 m^2) <= (u^2+v^2)/(8M) = (u^2+v^2)/(4K)).
    Returns {('axis',d): (a_K, tail)} and {('diag',c): ...} merged (diag earned d = 2c)."""
    targets = [('axis', d, d, d) for d in dlist_axis] + \
              [('diag', c, 2 * c, 0) for c in dlist_diag]   # (kind, label, u, v)
    acc = {}
    steppers = {}
    for kind, lab, u, v in targets:
        acc[(kind, lab)] = 1            # k = 0 term: c_0 = 1 - 0 = 1 (target != origin)
        steppers[(kind, lab)] = None
    c0 = BinStepper(0)
    k = 0
    while k + 2 <= K:
        c0.step2()
        k += 2
        c0sq = c0.val * c0.val
        for kind, lab, u, v in targets:
            key = (kind, lab)
            st = steppers[key]
            if st is None and max(abs(u), abs(v)) <= k:
                su, sv = BinStepper(u), BinStepper(v)
                while su.k < k:
                    su.step2()
                while sv.k < k:
                    sv.step2()
                st = (su, sv)
                steppers[key] = st
            elif st is not None:
                st[0].step2()
                st[1].step2()
            nk = st[0].val * st[1].val if st is not None else 0
            acc[key] = acc[key] * 16 + (c0sq - nk)
    den = 4 ** K
    out = {}
    for kind, lab, u, v in targets:
        tail = Fraction(u * u + v * v, 4 * K)
        out[(kind, lab)] = (Fraction(acc[(kind, lab)], den), tail)
    return out


def potential_kernel_1d(dlist, K):
    """Regularized critical kernel on Z at the COMPUTED mu_c = 1/2, even d:
        a_K(d) = sum_{k<=K,even} (C(k,k/2) - C(k,(k+d)/2)) / 2^k     (EXACT rational)
    with PROVEN tail bound sum_{k>K} c_k <= d^2/(2 sqrt(2M)), M = K/2, certified as the
    rational d^2/(2 isqrt(2M)) (floor sqrt only ENLARGES the bound):
      c_{2m}(d) <= p_{2m}(0) s^2/m, s=d/2;  p_{2m}(0) <= 1/sqrt(2m+1)   [from L2]
      => c_{2m} <= d^2/(4 m sqrt(2m));  sum_{m>M} m^{-3/2} <= 2/sqrt(M)."""
    acc = {d: 1 for d in dlist}     # k = 0 term: c_0 = 1
    steppers = {d: None for d in dlist}
    c0 = BinStepper(0)
    k = 0
    while k + 2 <= K:
        c0.step2()
        k += 2
        for d in dlist:
            st = steppers[d]
            if st is None and d <= k:
                st = BinStepper(d)
                while st.k < k:
                    st.step2()
                steppers[d] = st
            elif st is not None:
                st.step2()
            nk = st.val if st is not None else 0
            acc[d] = acc[d] * 4 + (c0.val - nk)
    den = 2 ** K
    M = K // 2
    out = {}
    for d in dlist:
        tail = Fraction(d * d, 2 * isqrt(2 * M))
        out[d] = (Fraction(acc[d], den), tail)
    return out


# ---------------------------------------------------------------- exact torus resolvent
def resolvent_exact(adj, mu, src):
    """EXACT rational solve of (I - mu A) x = e_src on the finite venue (Fraction
    Gaussian elimination).  Returns list of Fractions, or None if singular."""
    n = len(adj)
    mu = Fraction(mu)
    M = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = Fraction(1)
        for j, mult in adj[i]:
            M[i][j] -= mu * mult
    rhs = [Fraction(0)] * n
    rhs[src] = Fraction(1)
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        rhs[col], rhs[piv] = rhs[piv], rhs[col]
        inv = 1 / M[col][col]
        M[col] = [x * inv for x in M[col]]
        rhs[col] *= inv
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [xr - f * xc for xr, xc in zip(M[r], M[col])]
                rhs[r] -= f * rhs[col]
    return rhs


def kernel_at_quarter(adj):
    """Computed check that (I - A/4) annihilates the constant vector on the venue
    (locates the resolvent pole AT mu_c = 1/4 exactly)."""
    for i in range(len(adj)):
        s = Fraction(1)
        for j, mult in adj[i]:
            s -= Fraction(mult, 4)
        if s != 0:
            return False
    return True


# ---------------------------------------------------------------- interval helpers (exact)
def ratio_interval(g_num, g_den):
    """[lo,hi] Fractions for G(t2)/G(t1) with G in [S, S+tail] on both sides."""
    (s2, t2), (s1, t1) = g_num, g_den
    return (s2 / (s1 + t1), (s2 + t2) / s1)


def sqrt_bracket(fr, prec=10 ** 12):
    """Certified rational bracket lo <= sqrt(fr) <= hi via isqrt on scaled integers.
    COMPARISON USE ONLY; never on the measurement path."""
    fr = Fraction(fr)
    scaled = fr.numerator * prec * prec // fr.denominator
    r = isqrt(scaled)
    return Fraction(r, prec), Fraction(r + 2, prec)


# ---------------------------------------------------------------- coset weight histogram
def coset_weight_histogram(gens, rep, nbits):
    """FULL weight histogram over the coset rep ^ span(gens) (numpy uint64, exact ints).
    Extends the O-54-C exhaustive scan from (min, count) to the whole enumerator."""
    arr = np.zeros(1, dtype=np.uint64)
    for g in gens:
        arr = np.concatenate([arr, arr ^ np.uint64(g)])
    w = _np_pop(arr ^ np.uint64(rep))
    hist = np.bincount(w.astype(np.int64), minlength=nbits + 1)
    return [int(x) for x in hist]
