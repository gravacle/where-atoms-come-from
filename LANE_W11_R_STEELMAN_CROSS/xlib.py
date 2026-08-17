# LANE W11-R-CROSS  --  shared library.  CROSS-REFUTER of LANE_W11_R_STEELMAN.
# Rebuilt from the SEALED BYTES, not from either W-11 lane's code.  Where I reproduce a number of
# theirs I say so explicitly at the call site.
#   K1 : S1_CARRIER_K1_V001.md sec1  (5 vertices, 6 edges, faces, root v0).
#   B0b: 3x3 torus grid.  I build the incidence MYSELF from an explicit vertex/edge list and CHECK
#        the class multiset against S4_THE_MEASUREMENT_V001.md:575  {00:4, 01:1, 10:2, 11:2}.
# DOUBLE PRECISION IS THE DEFAULT.  Exact claims are marked and done in fractions.Fraction.
import numpy as np

# ------------------------------------------------------------------ K1, from S1 sec1 :16-22
# e1: v0->v1  e2: v1->v2  e3: v2->v0   e4: v0->v3  e5: v3->v4  e6: v4->v0
K1_NV = 5
K1_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
# a loop is a list of (src, dst, edge_index, sign); sign +1 means traversed with the orientation
K1_LOOP_F = [(0, 1, 0, +1), (1, 2, 1, +1), (2, 0, 2, +1)]      # gamma_F = e1.e2.e3, filled
K1_LOOP_C = [(0, 3, 3, +1), (3, 4, 4, +1), (4, 0, 5, +1)]      # gamma_C = e4.e5.e6, unfilled

# ------------------------------------------------------------------ B0b, 3x3 torus, built here
def _build_b0b():
    """3x3 torus grid.  V(i,j) = 3j+i, i,j in Z_3.  Horizontal edge H(i,j): V(i,j)->V(i+1,j),
    index 3j+i.  Vertical edge W(i,j): V(i,j)->V(i,j+1), index 9+3j+i.
    gamma_F = boundary of face (0,0), length 4 :  V(0,0)->V(1,0)->V(1,1)->V(0,1)->V(0,0)
    gamma_C = the horizontal row j=0, length 3 :  V(0,0)->V(1,0)->V(2,0)->V(0,0)"""
    V = lambda i, j: 3 * (j % 3) + (i % 3)
    H = lambda i, j: 3 * (j % 3) + (i % 3)
    Wd = lambda i, j: 9 + 3 * (j % 3) + (i % 3)
    E = [None] * 18
    for j in range(3):
        for i in range(3):
            E[H(i, j)] = (V(i, j), V(i + 1, j))
            E[Wd(i, j)] = (V(i, j), V(i, j + 1))
    # face (0,0) boundary, walked as a closed path
    gF = [(V(0, 0), V(1, 0), H(0, 0), +1),
          (V(1, 0), V(1, 1), Wd(1, 0), +1),
          (V(1, 1), V(0, 1), H(0, 1), -1),      # traversed against orientation
          (V(0, 1), V(0, 0), Wd(0, 0), -1)]
    gC = [(V(0, 0), V(1, 0), H(0, 0), +1),
          (V(1, 0), V(2, 0), H(1, 0), +1),
          (V(2, 0), V(0, 0), H(2, 0), +1)]
    return E, gF, gC, 9

B0B_E, B0B_LOOP_F, B0B_LOOP_C, B0B_NV = _build_b0b()

# ------------------------------------------------------------------ operators
def loop_vertices(loop):
    return {t[0] for t in loop}

def holonomy(loop, a):
    z = 1.0 + 0j
    for (src, dst, e, sg) in loop:
        z *= np.exp(1j * a[e]) if sg > 0 else np.exp(-1j * a[e])
    return z

def T_edge(loop, a, NV):
    """COR-F's EDGE TICK (S3_THE_CROSSING_AUDIT_V001.md:160-209, :794): move each fibre value one
       edge along the loop, identity off the loop."""
    on = loop_vertices(loop)
    T = np.zeros((NV, NV), dtype=complex)
    for v in range(NV):
        if v not in on:
            T[v, v] = 1.0
    for (src, dst, e, sg) in loop:
        T[dst, src] = np.exp(1j * a[e]) if sg > 0 else np.exp(-1j * a[e])
    return T

def M_circuit(loop, a, NV):
    """THE CORPUS'S OPERATOR (W-01, REGISTER:31-35; S3 sec2.3): multiply s_v by the WHOLE-CIRCUIT
       holonomy for every v on the loop, identity elsewhere."""
    W = holonomy(loop, a)
    M = np.eye(NV, dtype=complex)
    for v in loop_vertices(loop):
        M[v, v] = W
    return M

def D_uniform(loop, a, NV):
    """the uniform (principal) L-th root: diag(W^{1/L} on the loop, 1 off).  D^L = M exactly."""
    Lg = len(loop)
    W = holonomy(loop, a)
    r = np.exp(1j * np.angle(W) / Lg)
    D = np.eye(NV, dtype=complex)
    for v in loop_vertices(loop):
        D[v, v] = r
    return D

# ------------------------------------------------------------------ classes / pushforward
def classes(loopF, loopC, NV):
    F, C = loop_vertices(loopF), loop_vertices(loopC)
    return [(int(v in F), int(v in C)) for v in range(NV)]

def pi_of(s, loopF, loopC, NV):
    cl = classes(loopF, loopC, NV)
    w = np.abs(s) ** 2
    p = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
    for v in range(NV):
        p[cl[v]] += w[v]
    return np.array([p[(0, 0)], p[(1, 0)], p[(0, 1)], p[(1, 1)]])

def m_jensen(p, n=1 << 20):
    a, b, c, d = p
    t = 2 * np.pi * np.arange(n) / n
    ct = np.cos(t)
    A = np.sqrt(np.maximum(a * a + b * b + 2 * a * b * ct, 0))
    B = np.sqrt(np.maximum(c * c + d * d + 2 * c * d * ct, 0))
    return np.log(np.maximum(A, B) + 1e-300).mean()

def Z(opF, opC, s, nF, nC):
    return np.vdot(np.linalg.matrix_power(opF, nF) @ s, np.linalg.matrix_power(opC, nC) @ s)

def rate(opF, opC, s, N, stepF=1, stepC=1):
    AF = np.linalg.matrix_power(opF, stepF)
    AC = np.linalg.matrix_power(opC, stepC)
    xF, xC, tot = s.copy(), s.copy(), 0.0
    for _ in range(N):
        xF = AF @ xF
        xC = AC @ xC
        z = abs(np.vdot(xF, xC))
        tot += np.log(z) if z > 0 else -700.0
    return tot / N

def arms_differ(*arrays):
    """DIFF YOUR ARMS.  True only if every pair differs in BYTES."""
    bs = [np.ascontiguousarray(x).tobytes() for x in arrays]
    return all(bs[i] != bs[j] for i in range(len(bs)) for j in range(i + 1, len(bs)))

# ------------------------------------------------------------------ ready states with equal pi
def pi_identical_states_K1(rng):
    """three states with IDENTICAL pi, differing only in within-class distribution or in phase.
       Same construction as both W-11 lanes so the comparison is like-for-like."""
    sA = np.sqrt(np.array([0.40, 0.15, 0.15, 0.15, 0.15])) + 0j
    sB = np.sqrt(np.array([0.40, 0.30, 0.00, 0.05, 0.25])) + 0j
    sC = sA * np.exp(1j * np.array([0.0, 1.3, -0.7, 2.2, 0.4]))
    return sA, sB, sC

def random_pi_identical(rng, loopF, loopC, NV, base_w, k=8):
    """k states, all with the SAME pi as base_w, redistributed freely WITHIN each class and given
       random phases.  Used to test pi-only-ness against many states, not three."""
    cl = classes(loopF, loopC, NV)
    p = {}
    for v in range(NV):
        p.setdefault(cl[v], []).append(v)
    out = []
    for _ in range(k):
        w = np.zeros(NV)
        for c, vs in p.items():
            tot = sum(base_w[v] for v in vs)
            r = rng.random(len(vs))
            r = r / r.sum() if r.sum() > 0 else np.ones(len(vs)) / len(vs)
            for i, v in enumerate(vs):
                w[v] = tot * r[i]
        out.append(np.sqrt(w) * np.exp(1j * rng.uniform(0, 2 * np.pi, NV)))
    return out

# ------------------------------------------------------------------ root variety of M_gamma
def _rand_unitary(d, rng):
    Q, R = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))

def random_root(loop, a, NV, rng, kind="generic"):
    """sample a unitary U with U^L = M_gamma(a) EXACTLY (to machine precision).
       Any such U commutes with M, so it is block diagonal on (loop set, complement); on the loop
       block U = w.V with w^L = W and V^L = I, off it U = V' with V'^L = I."""
    Lg = len(loop)
    on = sorted(loop_vertices(loop))
    off = [v for v in range(NV) if v not in on]
    W = holonomy(loop, a)
    w0 = np.exp(1j * (np.angle(W) + 2 * np.pi * int(rng.integers(0, Lg))) / Lg)
    z = np.exp(2j * np.pi / Lg)

    def blk(d):
        ks = rng.integers(0, Lg, size=d)
        if kind == "generic":
            Q = _rand_unitary(d, rng)
            return Q @ np.diag(z ** ks.astype(float)) @ Q.conj().T
        if kind == "diag":
            return np.diag(z ** ks.astype(float))
        raise ValueError(kind)

    U = np.zeros((NV, NV), dtype=complex)
    if on:
        B = w0 * blk(len(on))
        for i, v in enumerate(on):
            for j, u in enumerate(on):
                U[v, u] = B[i, j]
    if off:
        B = blk(len(off))
        for i, v in enumerate(off):
            for j, u in enumerate(off):
                U[v, u] = B[i, j]
    return U

# ------------------------------------------------------------------ gauge
def gauge_apply(a, th, EDGES):
    """S1 :59-63.  a_e -> a_e + theta_dst - theta_src on e: src->dst."""
    a2 = a.copy()
    for e, (u, v) in enumerate(EDGES):
        a2[e] = a[e] + th[v] - th[u]
    return a2
