"""
r2_lib.py -- LANE W10-A REFUTE-2.  Written fresh.  Opposite d1 sign convention to the lane's,
transposed B0b indexing, independent B4 list.  See CONVENTIONS.txt.
"""
from fractions import Fraction
from itertools import combinations
import numpy as np

CLASSES = [(0, 0), (1, 0), (0, 1), (1, 1)]          # order (00,10,01,11) -- same as the lane's,
CNAME = {(0, 0): "00", (1, 0): "10", (0, 1): "01", (1, 1): "11"}   # because pi's order is fixed
                                                                   # by the register's own P(x,y)

# ------------------------------------------------------------------ exact rank over Q (mine)
def rank_Q(M):
    if not M or not M[0]:
        return 0
    A = [[Fraction(x) for x in r] for r in M]
    m, n = len(A), len(A[0])
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == m:
            break
    return r

def in_span(cols, vec):
    """EXACT: is vec in the Q-span of the list of column vectors `cols`?"""
    if not cols:
        return all(x == 0 for x in vec)
    M = [list(t) for t in zip(*cols)]
    return rank_Q(M) == rank_Q([row + [vec[i]] for i, row in enumerate(M)])


class Cx:
    """Oriented regular CW complex, dim <= 2.  d1[v,e] = +1 at TAIL, -1 at HEAD  (OPPOSITE of
    the lane's and of S1:57 -- deliberately, to prove no number here depends on the sign)."""
    def __init__(self, name, nV, edges, faces, gF, gC):
        self.name, self.nV, self.edges, self.faces = name, nV, list(edges), list(faces)
        self.gF, self.gC = list(gF), list(gC)
    @property
    def nE(self): return len(self.edges)
    @property
    def nF(self): return len(self.faces)
    def d1(self):
        M = [[0] * self.nE for _ in range(self.nV)]
        for e, (t, h) in enumerate(self.edges):
            M[t][e] += 1          # OPPOSITE SIGN CONVENTION
            M[h][e] -= 1
        return M
    def d2(self):
        M = [[0] * self.nF for _ in range(self.nE)]
        for f, cyc in enumerate(self.faces):
            for (e, s) in cyc:
                M[e][f] += s
        return M
    def chainvec(self, signed):
        v = [0] * self.nE
        for (e, s) in signed:
            v[e] += s
        return v
    def is_cycle(self, signed):
        d1 = self.d1(); v = self.chainvec(signed)
        return all(sum(d1[i][e] * v[e] for e in range(self.nE)) == 0 for i in range(self.nV))
    def bounds(self, signed):
        d2 = self.d2()
        cols = [[d2[e][f] for e in range(self.nE)] for f in range(self.nF)]
        return in_span(cols, self.chainvec(signed))
    def d1d2_zero(self):
        A = np.array(self.d1(), dtype=np.int64); B = np.array(self.d2(), dtype=np.int64)
        return 0 if self.nF == 0 else int(np.max(np.abs(A @ B)))
    def betti(self):
        r1 = rank_Q(self.d1()) if self.nE else 0
        r2 = rank_Q(self.d2()) if self.nF else 0
        return self.nV - r1, (self.nE - r1) - r2, self.nF - r2, r1, r2
    def support(self, signed):
        S = set()
        for (e, s) in signed:
            if s != 0:
                t, h = self.edges[e]; S.add(t); S.add(h)
        return S
    def multiset(self):
        VF, VC = self.support(self.gF), self.support(self.gC)
        cl = [(1 if v in VF else 0, 1 if v in VC else 0) for v in range(self.nV)]
        return {ab: sum(1 for c in cl if c == ab) for ab in CLASSES}
    def pi_uniform(self):
        ms = self.multiset()
        return [Fraction(ms[c], self.nV) for c in CLASSES]
    def independent(self):
        """gF, gC independent in the cycle space over Q."""
        return rank_Q([[self.chainvec(self.gF)[e], self.chainvec(self.gC)[e]]
                       for e in range(self.nE)]) == 2


# ------------------------------------------------------------------ carriers, built fresh
def my_K1():
    E = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    return Cx("B1 K1", 5, E, [[(0, 1), (1, 1), (2, 1)]],
              [(0, 1), (1, 1), (2, 1)], [(3, 1), (4, 1), (5, 1)])

def my_B1q():
    E = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 6), (6, 3)]
    return Cx("B1q", 7, E, [[(0, 1), (1, 1), (2, 1)]],
              [(0, 1), (1, 1), (2, 1)], [(3, 1), (4, 1), (5, 1)])

def my_B0b():
    """3x3 grid torus.  MY indexing: vertex (i,j) -> 3*i + j  (the lane uses 3*j+i)."""
    def V(i, j): return 3 * (i % 3) + (j % 3)
    H, W = {}, {}
    E = []
    for i in range(3):
        for j in range(3):
            H[(i, j)] = len(E); E.append((V(i, j), V(i + 1, j)))     # "horizontal": i -> i+1
    for i in range(3):
        for j in range(3):
            W[(i, j)] = len(E); E.append((V(i, j), V(i, j + 1)))     # "vertical":   j -> j+1
    F = []
    for i in range(3):
        for j in range(3):
            F.append([(H[(i, j)], 1), (W[((i + 1) % 3, j)], 1),
                      (H[(i, (j + 1) % 3)], -1), (W[(i, j)], -1)])
    gF = [(H[(0, 0)], 1), (W[(1, 0)], 1), (H[(0, 1)], -1), (W[(0, 0)], -1)]
    gC = [(H[(0, 0)], 1), (H[(1, 0)], 1), (H[(2, 0)], 1)]
    return Cx("B0b ring torus 3x3, loops MEET", 9, E, F, gF, gC)

def my_B4(gC_choice="a1b1"):
    """Spindle: two 2-spheres glued at two points p,q.  p=0 q=1 a1=2 a2=3 b1=4 b2=5.
    Each sphere = a square 1-cycle with TWO 2-cells attached (the two hemispheres)."""
    E = [(0, 2), (2, 1), (1, 3), (3, 0),      # sphere A square p a1 q a2
         (0, 4), (4, 1), (1, 5), (5, 0)]      # sphere B square p b1 q b2
    A = [(0, 1), (1, 1), (2, 1), (3, 1)]
    B = [(4, 1), (5, 1), (6, 1), (7, 1)]
    F = [list(A), list(A), list(B), list(B)]
    CH = {"a1b1": [(0, 1), (1, 1), (5, -1), (4, -1)],       # p a1 q b1
          "a1b2": [(0, 1), (1, 1), (6, 1), (7, 1)],         # p a1 q b2
          "a2b1": [(2, -1), (3, -1), (4, 1), (5, 1)],       # p a2 q b1  (reverse a2 leg)
          "a2b2": [(2, -1), (3, -1), (6, -1), (7, -1)],     # p a2 q b2
          "sqB":  list(B)}
    return Cx(f"B4 spindle (gC={gC_choice})", 6, E, F, A, CH[gC_choice])


# ------------------------------------------------------------------ lattices
def hnf(vecs):
    """EXACT canonical basis of the sublattice of Z^2 generated by vecs.  My own routine."""
    from math import gcd
    V = [(int(a), int(b)) for (a, b) in vecs if (a, b) != (0, 0)]
    if not V:
        return ()
    # column-style HNF: find g = gcd of first coords
    rows = [list(v) for v in V]
    while True:
        nz = [r for r in rows if r[0] != 0]
        if not nz:
            break
        piv = min(nz, key=lambda r: abs(r[0]))
        changed = False
        for r in rows:
            if r is not piv and r[0] != 0:
                q = r[0] // piv[0]
                r[0] -= q * piv[0]; r[1] -= q * piv[1]
                changed = True
        if not changed:
            break
    piv = next((r for r in rows if r[0] != 0), None)
    g = 0
    for r in rows:
        if r[0] == 0:
            g = gcd(g, abs(r[1]))
    if piv is None:
        return ((0, g),) if g else ()
    if piv[0] < 0:
        piv = [-piv[0], -piv[1]]
    if g == 0:
        return ((piv[0], piv[1]),)
    return ((piv[0], piv[1] % g), (0, g))

def L_supp(S):
    return hnf([(a - a2, b - b2) for (a, b), (a2, b2) in combinations(S, 2)])

def L_conn(A, B, q):
    """EXACT relation lattice for alpha = 2pi A/q, beta = 2pi B/q."""
    if q == 1:
        return ((1, 0), (0, 1))
    return hnf([(m, n) for m in range(-q, q + 1) for n in range(-q, q + 1)
                if (m, n) != (0, 0) and (m * A + n * B) % q == 0])

def contained(a, b):
    """EXACT: lattice(a) subset lattice(b)?"""
    if not a:
        return True
    if not b:
        return False
    if len(b) == 1:
        bx, by = b[0]
        for (x, y) in a:
            if bx != 0:
                if x % bx: return False
                t = x // bx
            else:
                if x != 0 or by == 0 or y % by: return False
                t = y // by
            if (t * bx, t * by) != (x, y): return False
        return True
    det = b[0][0] * b[1][1] - b[0][1] * b[1][0]
    if det == 0:
        return False
    for (x, y) in a:
        n1 = x * b[1][1] - y * b[1][0]
        n2 = b[0][0] * y - b[0][1] * x
        if n1 % det or n2 % det:
            return False
    return True


# ------------------------------------------------------------------ the observable
def Zabs(pi4, a, b, K, k0=1):
    """|Z_k| for k = k0..k0+K-1, alpha/2pi = a, beta/2pi = b.  float64."""
    k = np.arange(k0, k0 + K, dtype=np.float64)
    u = np.exp(2j * np.pi * ((k * a) % 1.0))
    v = np.exp(2j * np.pi * ((k * b) % 1.0))
    p00, p10, p01, p11 = [float(x) for x in pi4]
    return np.abs(p00 + p10 * u + p01 * v + p11 * u * v)

def involute(pi4):
    """W-03's involution 00<->11, 10<->01 on the class weights."""
    p00, p10, p01, p11 = pi4
    return [p11, p01, p10, p00]
