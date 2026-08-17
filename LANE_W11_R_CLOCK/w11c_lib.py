# LANE W-11 R/C — THE CLOCK CORRESPONDENCE — shared library.
# Written from the sealed bytes, NOT copied from LANE_W11_CONVENTION_TEST.
#   K1 incidence:  S1_CARRIER_K1_V001.md sec1 (:16-22 edge list, :10-14 vertices/root)
#   B0b incidence: LANE_W10_A_CARRIERS/w10a_lib.py:182-213 (sealed), checked against
#                  S4_THE_MEASUREMENT_V001.md:575 class multiset {00:4,01:1,10:2,11:2}
#   COR-F edge tick T: S3_THE_CROSSING_AUDIT_V001.md:160-209, :794
#   CIRCUIT operator M_gamma: S3 sec2.3 / REGISTER:31-35 (W-01)
# Double precision is the default throughout. Exact/closed-form checks are marked where used.
import numpy as np

# ---------------------------------------------------------------- carriers
class Carrier:
    """A carrier is: nv vertices, an edge list, and two directed closed walks."""
    def __init__(self, name, nv, edges, walkF, walkC):
        self.name, self.nv, self.edges = name, nv, edges
        self.walkF, self.walkC = walkF, walkC          # [(src,dst,edge_index,sign)]
        self.LF, self.LC = len(walkF), len(walkC)
        self.VF = tuple(sorted({u for u, _, _, _ in walkF}))
        self.VC = tuple(sorted({u for u, _, _, _ in walkC}))
        # closure check: a walk must be a closed cycle visiting each of its vertices once
        for w in (walkF, walkC):
            for i in range(len(w)):
                assert w[i][1] == w[(i + 1) % len(w)][0], "walk is not closed"
            assert len({u for u, _, _, _ in w}) == len(w), "walk revisits a vertex"

    def cls(self, v):
        return (int(v in self.VF), int(v in self.VC))

    def multiset(self):
        from collections import Counter
        return dict(Counter("%d%d" % self.cls(v) for v in range(self.nv)))

def K1():
    # S1:16-22   e1: v0->v1  e2: v1->v2  e3: v2->v0   e4: v0->v3  e5: v3->v4  e6: v4->v0
    E = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    wF = [(0, 1, 0, +1), (1, 2, 1, +1), (2, 0, 2, +1)]      # filled triangle boundary
    wC = [(0, 3, 3, +1), (3, 4, 4, +1), (4, 0, 5, +1)]      # unfilled cycle
    return Carrier("K1", 5, E, wF, wC)

def B0b():
    # 3x3 square-grid torus. v(i,j)=3*(j%3)+(i%3); h(i,j) index 3j+i; w(i,j) index 9+3j+i.
    V = lambda i, j: 3 * (j % 3) + (i % 3)
    H = lambda i, j: 3 * (j % 3) + (i % 3)
    Wv = lambda i, j: 9 + 3 * (j % 3) + (i % 3)
    E = [None] * 18
    for j in range(3):
        for i in range(3):
            E[H(i, j)] = (V(i, j), V(i + 1, j))
            E[Wv(i, j)] = (V(i, j), V(i, j + 1))
    gF = [(H(0, 0), +1), (Wv(1, 0), +1), (H(0, 1), -1), (Wv(0, 0), -1)]   # face (0,0) boundary
    gC = [(H(0, 0), +1), (H(1, 0), +1), (H(2, 0), +1)]                    # row j=0, flat
    def to_walk(g):
        out = []
        for (e, s) in g:
            u, v = E[e]
            out.append((u, v, e, +1) if s > 0 else (v, u, e, -1))
        return out
    return Carrier("B0b", 9, E, to_walk(gF), to_walk(gC))

# ---------------------------------------------------------------- operators
def Tedge(K, walk, a):
    """COR-F's EDGE TICK: move each fibre value one edge along the loop, identity off it.
       (Ts)(dst) = U_e^{+-1} s(src).  Unitary; T^L = M_gamma."""
    U = np.exp(1j * np.asarray(a, dtype=float))
    T = np.zeros((K.nv, K.nv), dtype=complex)
    on = {u for u, _, _, _ in walk}
    for v in range(K.nv):
        if v not in on:
            T[v, v] = 1.0
    for (u, v, e, sgn) in walk:
        T[v, u] = U[e] if sgn > 0 else np.conj(U[e])
    return T

def hol(walk, a):
    z = 1.0 + 0j
    for (_, _, e, sgn) in walk:
        z *= np.exp(1j * a[e]) if sgn > 0 else np.exp(-1j * a[e])
    return z

def Mcirc(K, vs, W):
    """CIRCUIT operator: multiply every fibre ON the loop by the whole-circuit holonomy."""
    M = np.eye(K.nv, dtype=complex)
    for v in vs:
        M[v, v] = W
    return M

def ops(K, a):
    TF, TC = Tedge(K, K.walkF, a), Tedge(K, K.walkC, a)
    WF, WC = hol(K.walkF, a), hol(K.walkC, a)
    return TF, TC, Mcirc(K, K.VF, WF), Mcirc(K, K.VC, WC), WF, WC

# ---------------------------------------------------------------- states
def pi_of(K, s):
    w = np.abs(s) ** 2
    p = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
    for v in range(K.nv):
        p[K.cls(v)] += w[v]
    return np.array([p[(0, 0)], p[(1, 0)], p[(0, 1)], p[(1, 1)]])

def classes(K):
    out = {}
    for v in range(K.nv):
        out.setdefault(K.cls(v), []).append(v)
    return out

def states_same_pi(K, pi, n, rng, phases=True):
    """n ready states with EXACTLY the class sums pi, differing only within classes / in phase."""
    cl, out = classes(K), []
    order = [(0, 0), (1, 0), (0, 1), (1, 1)]
    for _ in range(n):
        w = np.zeros(K.nv)
        for idx, c in enumerate(order):
            vs = cl.get(c, [])
            if not vs:
                assert abs(pi[idx]) < 1e-14, "pi puts weight on an unoccupied class"
                continue
            x = rng.dirichlet(np.ones(len(vs)))
            for v, xv in zip(vs, x):
                w[v] = pi[idx] * xv
        ph = rng.uniform(0, 2 * np.pi, K.nv) if phases else np.zeros(K.nv)
        out.append(np.sqrt(w) * np.exp(1j * ph))
    return out

# ---------------------------------------------------------------- observable
def Zlat(s, TF, TC, mF, mC):
    """THE GENERAL OBJECT.  Z(mF,mC) = <T_F^{mF} s , T_C^{mC} s>.
       Both conventions are RAYS in this lattice:
         EDGE clock    (n, n)          CIRCUIT clock  (L_F k, L_C k)"""
    return np.vdot(np.linalg.matrix_power(TF, mF) @ s,
                   np.linalg.matrix_power(TC, mC) @ s)

def m_jensen(p, n=1 << 22):
    """log Mahler measure of p00 + p10 x + p01 y + p11 xy, by Jensen's formula (N1)."""
    a, b, c, d = p
    t = 2 * np.pi * np.arange(n) / n
    ct = np.cos(t)
    A = np.sqrt(np.maximum(a * a + b * b + 2 * a * b * ct, 0.0))
    B = np.sqrt(np.maximum(c * c + d * d + 2 * c * d * ct, 0.0))
    return float(np.log(np.maximum(np.maximum(A, B), 1e-300)).mean())

def arms_differ(name, *arrays):
    """ANTI-ZERO-VARIABLE GUARD (W-08's finding; W-10 N-6).  Hash the ARMS THEMSELVES."""
    import hashlib
    hs = [hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()[:16] for x in arrays]
    ok = len(set(hs)) == len(hs)
    print(f"    ARMS-DIFF [{name}]: {hs}  ->  {'DISTINCT' if ok else '*** BYTE-IDENTICAL ***'}")
    assert ok, "two arms are byte-identical — the control is void"

def generic_conn(K, rng, f=1.0, c=2**0.5):
    """A connection with EVERY edge phase non-zero and holonomies (W_F,W_C) = (e^{if}, e^{ic}).
    (f,c) = (1.0, sqrt(2)) is THE ONLY GENERIC CONNECTION THE CORPUS PUBLISHES (W-10 N-4; S4:603).
    Genericity in one line: a relation needs m f + n c = 2 pi j with (m,n) != 0; n = 0 forces
    m = j = 0 since pi is irrational, and n != 0 makes c rational (j=0) or transcendental (j!=0),
    while sqrt(2) is neither."""
    a = rng.uniform(0.2, 2.0, len(K.edges))
    setC = {ee for (_, _, ee, _) in K.walkC}
    setF = {ee for (_, _, ee, _) in K.walkF}
    eF = [e for (_, _, e, _) in K.walkF if e not in setC]
    eC = [e for (_, _, e, _) in K.walkC if e not in setF]
    assert eF and eC, "need one private edge in each loop"
    def phase(walk):
        return sum(a[e] if s > 0 else -a[e] for (_, _, e, s) in walk)
    sF = [s for (_, _, e, s) in K.walkF if e == eF[0]][0]
    a[eF[0]] += sF * (f - phase(K.walkF))
    sC = [s for (_, _, e, s) in K.walkC if e == eC[0]][0]
    a[eC[0]] += sC * (c - phase(K.walkC))
    return a
