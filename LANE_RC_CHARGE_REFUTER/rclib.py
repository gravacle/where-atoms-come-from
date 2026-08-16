"""
LANE_RC_CHARGE_REFUTER -- rclib.py
Refuter on the CHARGE axis.  Written from scratch; no code reused from any other lane.
Conventions published here in full, per the W-03 defect-of-record.

CARRIER K1 (S1_CARRIER_K1_V001.md, sha256 verified this session)
  vertices  v0..v4          (index 0..4)
  edges     e1: v0->v1   e2: v1->v2   e3: v2->v0
            e4: v0->v3   e5: v3->v4   e6: v4->v0
  face      F attached along e1.e2.e3   (filled triangle v0 v1 v2)
  gamma_F = e1.e2.e3   (boundary of the filled 2-cell)
  gamma_C = e4.e5.e6   (the unfilled triangle, the 1-cycle)

INCIDENCE MATRICES (published, per the W-03 requirement)
  d1 : C_1 -> C_0 , 5 x 6, d1[v,e] = +1 if e targets v, -1 if e sources v
  d2 : C_2 -> C_1 , 6 x 1, d2[e,F] = coefficient of e in boundary of F

TRANSPORT AND CHARGE
  U_e = exp(i a_e) in U(1).  W_F = exp(i(a1+a2+a3)), W_C = exp(i(a4+a5+a6)).
  u := conj(W_F) = exp(-i f) ,  v := W_C = exp(i c)   with f = a1+a2+a3, c = a4+a5+a6.

  CHARGE.  S2:174 defines the charge knob as U_e |-> exp(i q a_e), still rank one.
  Carried to the vertex fibres (the only place a per-vertex charge can live on a
  rank-one bundle): the fibre at v carries charge q_v in Z, so W-01's loop operator is

      (M_gamma s)(v) = W(gamma)^{q_v} s(v)   if v lies on gamma ;  s(v) otherwise.

  q_v = 1 for all v recovers S4 exactly.  This is the ONLY charge convention used here
  and it is declared, not assumed silently.

  Then with p_v = |s_v|^2,
      Z_k = <M_F^k s, M_C^k s> = sum_v u^{k q_v a_v} v^{k q_v b_v} p_v
  where a_v = [v on gamma_F], b_v = [v on gamma_C].

  EXPONENT VECTOR of vertex v:   E(v) = q_v * (a_v, b_v)  in Z^2.
  EXPONENT POINT SET S = { E(v) : p_v > 0 } with pushforward weights.
  DIFFERENCE LATTICE  Delta(S) = < E(a) - E(b) : a,b in S >  <= Z^2.
  RELATION LATTICE    L = { (m,n) in Z^2 : u^m v^n = 1 }.
  FORMATION CRITERION under test:  formation  <=>  Delta(S) NOT contained in L.

SCHEDULES
  A (uniform)        k_n = 1      lambda_A = log|Z_1|
  B (canonical clock) k_n = n     lambda_B = lim (1/N) sum_{n=1..N} log|Z_n|
"""

import numpy as np
from math import pi, log, cos, sin, sqrt, gcd

# ----------------------------------------------------------------- the carrier
V = 5
E = 6
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]   # (source, target)


def d1_matrix():
    d1 = np.zeros((V, E), dtype=int)
    for j, (s, t) in enumerate(EDGES):
        d1[s, j] -= 1
        d1[t, j] += 1
    return d1


def d2_matrix():
    # face F = e1 + e2 + e3  (indices 0,1,2)
    d2 = np.zeros((E, 1), dtype=int)
    d2[0, 0] = 1
    d2[1, 0] = 1
    d2[2, 0] = 1
    return d2


# vertices on gamma_F = e1.e2.e3 -> {v0,v1,v2};  on gamma_C = e4.e5.e6 -> {v0,v3,v4}
A_INC = np.array([1, 1, 1, 0, 0], dtype=int)   # a_v
B_INC = np.array([1, 0, 0, 1, 1], dtype=int)   # b_v


def exponents(q):
    """E(v) = q_v * (a_v, b_v), shape (5,2)."""
    q = np.asarray(q, dtype=int)
    return np.stack([q * A_INC, q * B_INC], axis=1)


def support_points(q, p, tol=0.0):
    """Collapse vertices onto distinct exponent points; return (points, weights)."""
    Ev = exponents(q)
    acc = {}
    for v in range(V):
        if p[v] > tol:
            key = (int(Ev[v, 0]), int(Ev[v, 1]))
            acc[key] = acc.get(key, 0.0) + float(p[v])
    pts = sorted(acc.keys())
    w = np.array([acc[k] for k in pts], dtype=float)
    return np.array(pts, dtype=int), w


def delta_rank(pts):
    """Rank over Z (== rank over Q) of the lattice of differences of the point set."""
    if len(pts) <= 1:
        return 0
    D = pts[1:] - pts[0]
    return int(np.linalg.matrix_rank(D.astype(float)))


def delta_basis(pts):
    """Hermite-ish basis of Delta(S) as a sublattice of Z^2, returned as a 2x2 int
    matrix whose rows span it (zero rows padded), plus the rank."""
    if len(pts) <= 1:
        return np.zeros((2, 2), dtype=int), 0
    D = [tuple(map(int, x)) for x in (pts[1:] - pts[0])]
    # integer row reduction (2 columns only -- do it by hand, exactly)
    rows = [list(r) for r in D]
    basis = []
    # first: reduce on column 0 using gcd
    while rows:
        rows = [r for r in rows if r != [0, 0]]
        if not rows:
            break
        nz = [r for r in rows if r[0] != 0]
        if not nz:
            g = 0
            for r in rows:
                g = gcd(g, abs(r[1]))
            basis.append([0, g])
            break
        # gcd of column 0
        g = 0
        for r in nz:
            g = gcd(g, abs(r[0]))
        # build a vector with first entry g by integer combination (extended gcd chain)
        cur = None
        for r in nz:
            if cur is None:
                cur = r[:]
            else:
                cur = _bezout_combine(cur, r)
        if cur[0] < 0:
            cur = [-cur[0], -cur[1]]
        basis.append(cur)
        newrows = []
        for r in rows:
            k = r[0] // cur[0] if cur[0] != 0 else 0
            rr = [r[0] - k * cur[0], r[1] - k * cur[1]]
            newrows.append(rr)
        rows = [r for r in newrows if r != [0, 0]]
        if not rows:
            break
        # remaining all have first entry 0
        g = 0
        for r in rows:
            g = gcd(g, abs(r[1]))
        if g:
            basis.append([0, g])
        break
    B = np.zeros((2, 2), dtype=int)
    for i, b in enumerate(basis[:2]):
        B[i] = b
    return B, len(basis)


def _bezout_combine(r1, r2):
    """Return an integer combination of r1,r2 whose first entry is gcd(r1[0],r2[0])."""
    a, b = r1[0], r2[0]
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        qq = old_r // r
        old_r, r = r, old_r - qq * r
        old_s, s = s, old_s - qq * s
        old_t, t = t, old_t - qq * t
    return [old_s * r1[0] + old_t * r2[0], old_s * r1[1] + old_t * r2[1]]


# ------------------------------------------------------------------ transports
def Z_closed(k, f, c, q, p):
    """Closed form: Z_k = sum_v u^{k q_v a_v} v^{k q_v b_v} p_v ,  u=e^{-if}, v=e^{ic}."""
    Ev = exponents(q)
    ph = k * (-f * Ev[:, 0] + c * Ev[:, 1])
    return complex(np.sum(np.asarray(p, dtype=float) * np.exp(1j * ph)))


def Z_direct(k, f, c, q, p, seed=None):
    """INDEPENDENT re-derivation: build M_F, M_C as 5x5 matrices acting on C^5, take
    a section s with |s_v|^2 = p_v and ARBITRARY phases, and compute <M_F^k s, M_C^k s>.
    Uses the actual edge holonomies, summed around the loops -- not the closed form."""
    rng = np.random.default_rng(seed if seed is not None else 20260816)
    # arbitrary edge angles a_1..a_6 realising the given f and c
    a = np.zeros(6)
    a[0], a[1] = rng.uniform(0, 2 * pi, 2)
    a[2] = f - a[0] - a[1]
    a[3], a[4] = rng.uniform(0, 2 * pi, 2)
    a[5] = c - a[3] - a[4]
    WF = np.exp(1j * (a[0] + a[1] + a[2]))
    WC = np.exp(1j * (a[3] + a[4] + a[5]))
    q = np.asarray(q, dtype=int)
    MF = np.diag([WF ** (q[v] * A_INC[v]) for v in range(V)])
    MC = np.diag([WC ** (q[v] * B_INC[v]) for v in range(V)])
    phases = rng.uniform(0, 2 * pi, V)
    s = np.sqrt(np.asarray(p, dtype=float)) * np.exp(1j * phases)
    xF = np.linalg.matrix_power(MF, k) @ s
    xC = np.linalg.matrix_power(MC, k) @ s
    return complex(np.vdot(xF, xC))          # vdot conjugates the FIRST argument


# ------------------------------------------------------------------- lambda_B
def lambda_B_direct(f, c, q, p, N=2_000_000, chunk=200_000):
    """Schedule B by direct simulation: (1/N) sum_{n=1..N} log|Z_n|."""
    Ev = exponents(q)
    pw = np.asarray(p, dtype=float)
    keep = pw > 0
    Ek = Ev[keep].astype(float)
    pk = pw[keep]
    tot = 0.0
    done = 0
    while done < N:
        m = min(chunk, N - done)
        n = np.arange(done + 1, done + m + 1, dtype=np.float64)
        ph = np.outer(n, (-f * Ek[:, 0] + c * Ek[:, 1]))
        Zn = np.exp(1j * ph) @ pk
        tot += float(np.sum(np.log(np.abs(Zn))))
        done += m
    return tot / N


def lambda_B_torus(q, p, n=4000, seed=0):
    """Generic connection (rank L = 0, H = T^2): lambda_B = Mahler measure
    m(sum_e p_e x^{e1} y^{e2}) computed on an n x n midpoint grid."""
    pts, w = support_points(q, p)
    th = (np.arange(n) + 0.5) * 2 * pi / n
    X, Y = np.meshgrid(th, th, indexing='ij')
    Zs = np.zeros_like(X, dtype=complex)
    for (e1, e2), wt in zip(pts, w):
        Zs += wt * np.exp(1j * (e1 * X + e2 * Y))
    return float(np.mean(np.log(np.abs(Zs) + 1e-300)))


def mahler_3term_jensen(p0, p1, p2, n=2_000_000):
    """m(p0 + p1 x + p2 y) via Jensen in y:
       int log|A + p2 y| dy/2pi = log max(|A|, p2).  Exact 1-D reduction."""
    th = (np.arange(n) + 0.5) * 2 * pi / n
    A = np.abs(p0 + p1 * np.exp(1j * th))
    return float(np.mean(np.log(np.maximum(A, p2))))


def mahler_1var_exact(coeffs):
    """m(sum coeffs[j] w^j) = log|lead| + sum_{|root|>1} log|root|  (Jensen)."""
    c = np.array(coeffs, dtype=float)
    while len(c) > 1 and c[-1] == 0:
        c = c[:-1]
    lead = c[-1]
    r = np.roots(c[::-1])
    return float(np.log(abs(lead)) + np.sum(np.log(np.abs(r[np.abs(r) > 1]))))


def lambda_B_finite_orbit(f, c, q, p, L_basis):
    """rank L = 2 : the orbit {(u^n, v^n)} is finite; average log|Z| over it exactly."""
    # period = smallest T with u^T = v^T = 1
    T = None
    for t in range(1, 20001):
        if abs((t * f) % (2 * pi)) < 1e-9 or abs((t * f) % (2 * pi) - 2 * pi) < 1e-9:
            if abs((t * c) % (2 * pi)) < 1e-9 or abs((t * c) % (2 * pi) - 2 * pi) < 1e-9:
                T = t
                break
    if T is None:
        return None, None
    vals = [np.log(abs(Z_closed(n, f, c, q, p))) for n in range(1, T + 1)]
    return float(np.mean(vals)), T


def relation_lattice_rank(f, c, mmax=200, tol=1e-9):
    """rank of L = {(m,n): -m f + n c = 0 mod 2pi}, searched over |m|,|n| <= mmax."""
    found = []
    for m in range(-mmax, mmax + 1):
        for nn in range(-mmax, mmax + 1):
            if m == 0 and nn == 0:
                continue
            x = (-m * f + nn * c) % (2 * pi)
            if min(x, 2 * pi - x) < tol:
                found.append((m, nn))
    if not found:
        return 0, []
    Fm = np.array(found, dtype=float)
    return int(np.linalg.matrix_rank(Fm)), found[:8]


def in_lattice(vec, gens, tol=1e-9):
    """Is (m,n) in L?  Test directly: u^m v^n == 1."""
    raise NotImplementedError
