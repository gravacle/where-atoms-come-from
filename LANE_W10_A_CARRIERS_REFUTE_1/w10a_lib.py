"""
w10a_lib.py -- LANE W10-A, "Where Atoms Come From".
Carriers, exact ranks, class pushforwards, and exact arithmetic in Q[sqrt3].

Written fresh for this lane.  Declared in PUBLISHED_CONVENTIONS.txt: I read
LANE_R_MAPS_REFUTER/rm_lib.py:386-419 (an earlier round's 3x3 torus) before writing B0b,
so B0b's LOOP DESIGNATIONS are a reproduction, not an independent choice.  B4 has no prior
code anywhere in the corpus.

Conventions: see PUBLISHED_CONVENTIONS.txt.  IEEE double unless a function name says EXACT.
"""

from fractions import Fraction
from itertools import combinations
import numpy as np


# ----------------------------------------------------------------- exact rank over Q

def exact_rank(M):
    """Rank over Q of an integer matrix given as list of lists. Fraction elimination."""
    if not M or not M[0]:
        return 0
    A = [[Fraction(x) for x in row] for row in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                fac = A[i][c]
                A[i] = [a - fac * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r


def in_column_span(cols, vec):
    """EXACT: is integer vector `vec` in the Q-span of the integer columns `cols`?"""
    if not cols:
        return all(x == 0 for x in vec)
    M = [list(row) for row in zip(*cols)]          # rows = ambient coords, cols = generators
    r0 = exact_rank(M)
    M2 = [row + [vec[i]] for i, row in enumerate(M)]
    r1 = exact_rank(M2)
    return r0 == r1


# ----------------------------------------------------------------- carrier

class Carrier:
    def __init__(self, name, vnames, edges, faces, gF, gC):
        self.name = name
        self.vnames = list(vnames)
        self.edges = list(edges)        # (tail, head)
        self.faces = list(faces)        # list of [(edge, sign), ...]
        self.gF = list(gF)              # [(edge, sign), ...]
        self.gC = list(gC)

    nV = property(lambda s: len(s.vnames))
    nE = property(lambda s: len(s.edges))
    nF = property(lambda s: len(s.faces))

    def d1(self):
        M = [[0] * self.nE for _ in range(self.nV)]
        for e, (t, h) in enumerate(self.edges):
            M[t][e] -= 1
            M[h][e] += 1
        return M

    def d2(self):
        M = [[0] * self.nF for _ in range(self.nE)]
        for f, cyc in enumerate(self.faces):
            for (e, s) in cyc:
                M[e][f] += s
        return M

    def chain(self, signed):
        v = [0] * self.nE
        for (e, s) in signed:
            v[e] += s
        return v

    def d1_times(self, chain):
        out = [0] * self.nV
        for e, x in enumerate(chain):
            if x:
                t, h = self.edges[e]
                out[t] -= x
                out[h] += x
        return out

    def d1d2_max(self):
        A = np.array(self.d1(), dtype=np.int64)
        B = np.array(self.d2(), dtype=np.int64)
        if self.nF == 0:
            return 0
        return int(np.max(np.abs(A @ B)))

    def betti(self):
        r1 = exact_rank(self.d1()) if self.nE else 0
        r2 = exact_rank(self.d2()) if self.nF else 0
        return (self.nV - r1, (self.nE - r1) - r2, self.nF - r2, r1, r2)

    def d2_columns(self):
        M = self.d2()
        return [[M[e][f] for e in range(self.nE)] for f in range(self.nF)]

    def loop_vertices(self, signed):
        S = set()
        for (e, s) in signed:
            t, h = self.edges[e]
            S.add(t)
            S.add(h)
        return S

    def classes(self):
        """EXACT: per-vertex class (a,b) in {0,1}^2."""
        VF = self.loop_vertices(self.gF)
        VC = self.loop_vertices(self.gC)
        return [(1 if v in VF else 0, 1 if v in VC else 0) for v in range(self.nV)]

    def class_multiset(self):
        cl = self.classes()
        return {ab: sum(1 for c in cl if c == ab) for ab in [(0, 0), (0, 1), (1, 0), (1, 1)]}

    def pushforward_uniform(self):
        """EXACT: pi under SENSE U (p_v = 1/V), as Fractions in order (00,10,01,11)."""
        cl = self.classes()
        n = Fraction(1, self.nV)
        pi = {ab: Fraction(0) for ab in [(0, 0), (1, 0), (0, 1), (1, 1)]}
        for c in cl:
            pi[c] += n
        return [pi[(0, 0)], pi[(1, 0)], pi[(0, 1)], pi[(1, 1)]]

    def regular(self):
        """Each attaching cycle injective on its boundary: all its vertices distinct and it is a
        single closed edge path visiting each of its vertices once."""
        ok = True
        for cyc in self.faces:
            verts = []
            for (e, s) in cyc:
                t, h = self.edges[e]
                verts.append(t if s > 0 else h)
            ok = ok and (len(set(verts)) == len(verts))
        return ok


# ----------------------------------------------------------------- the carriers

def K1():
    """S1's carrier, from S1_CARRIER_K1_V001.md sec1 verbatim."""
    E = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    F = [[(0, 1), (1, 1), (2, 1)]]
    return Carrier("B1  K1 (the pinch, as handed)", ["v0", "v1", "v2", "v3", "v4"],
                   E, F, gF=[(0, 1), (1, 1), (2, 1)], gC=[(3, 1), (4, 1), (5, 1)])


def B1q():
    """K1-bridged + spectator, S4:519 row B1q.  Built to match S4's PUBLISHED multiset
    {00:1, 01:3, 10:3} and row V=7 E=8 F=1 chi=0 b1=1 b2=0.  Three classes, no 11."""
    #  filled triangle v0 v1 v2 ; unfilled triangle w0 w1 w2 ; bridge v0 - m - w0
    #  v0=0 v1=1 v2=2  w0=3 w1=4 w2=5  m=6
    E = [(0, 1), (1, 2), (2, 0),          # 0,1,2 filled triangle
         (3, 4), (4, 5), (5, 3),          # 3,4,5 unfilled triangle
         (0, 6), (6, 3)]                  # 6,7   bridge
    F = [[(0, 1), (1, 1), (2, 1)]]
    return Carrier("B1q K1-bridged + spectator", ["v0", "v1", "v2", "w0", "w1", "w2", "m"],
                   E, F, gF=[(0, 1), (1, 1), (2, 1)], gC=[(3, 1), (4, 1), (5, 1)])


def B0b():
    """
    B0b -- ring torus, 3x3 square-grid torus, LOOPS MEET.  S4:512, :539, :575.
    vertex (i,j) -> 3*j + i, i,j in Z_3.
    horizontal edge h(i,j): (i,j) -> (i+1,j)      index      3j+i
    vertical   edge w(i,j): (i,j) -> (i,j+1)      index  9 + 3j+i
    face (i,j) attached along  h(i,j) + w(i+1,j) - h(i,j+1) - w(i,j)   (anticlockwise square)
    gamma_F = the boundary of face (0,0)      -- BOUNDS, so W_F is a curvature
    gamma_C = the horizontal row j = 0        -- does NOT bound, so W_C is flat; MEETS gamma_F
    """
    def V(i, j):
        return 3 * (j % 3) + (i % 3)

    def H(i, j):
        return 3 * (j % 3) + (i % 3)

    def W(i, j):
        return 9 + 3 * (j % 3) + (i % 3)

    E = [None] * 18
    for j in range(3):
        for i in range(3):
            E[H(i, j)] = (V(i, j), V(i + 1, j))
            E[W(i, j)] = (V(i, j), V(i, j + 1))
    F = []
    for j in range(3):
        for i in range(3):
            F.append([(H(i, j), 1), (W(i + 1, j), 1), (H(i, j + 1), -1), (W(i, j), -1)])
    gF = [(H(0, 0), 1), (W(1, 0), 1), (H(0, 1), -1), (W(0, 0), -1)]
    gC = [(H(0, 0), 1), (H(1, 0), 1), (H(2, 0), 1)]
    names = [f"({i},{j})" for j in range(3) for i in range(3)]
    return Carrier("B0b ring torus 3x3 grid, loops MEET", names, E, F, gF, gC)


def B4():
    """
    B4 -- spindle: two 2-spheres glued at TWO points.  S4:515, :542, :578, and S4's ledger C8
    ("two spheres glued at two points (= my B4)").  S4 PUBLISHED NO INCIDENCE FOR B4; this is
    the first construction of it anywhere in the corpus.

    Each sphere is a square 1-cycle with TWO 2-cells attached along it (the two hemispheres):
        V=4 E=4 F=2, chi=2, b1=0, b2=1.  Regular: each attaching map is injective.
    Glue the two squares at their two OPPOSITE corners p and q.

        p=0  q=1  a1=2  a2=3  b1=4  b2=5
        sphere A square:  p -> a1 -> q -> a2 -> p    (edges 0,1,2,3)
        sphere B square:  p -> b1 -> q -> b2 -> p    (edges 4,5,6,7)

    gamma_F = sphere A's square       -- BOUNDS (it is the attaching cycle of face A+)
    gamma_C = p -> a1 -> q -> b1 -> p -- does NOT bound (it is the spindle's one free cycle)
    """
    E = [(0, 2), (2, 1), (1, 3), (3, 0),      # 0: p->a1  1: a1->q  2: q->a2  3: a2->p
         (0, 4), (4, 1), (1, 5), (5, 0)]      # 4: p->b1  5: b1->q  6: q->b2  7: b2->p
    sqA = [(0, 1), (1, 1), (2, 1), (3, 1)]
    sqB = [(4, 1), (5, 1), (6, 1), (7, 1)]
    F = [list(sqA), list(sqA), list(sqB), list(sqB)]      # two hemispheres per sphere
    gF = list(sqA)
    gC = [(0, 1), (1, 1), (5, -1), (4, -1)]   # p->a1->q then q->b1->p
    return Carrier("B4  spindle (two spheres, 2 glue pts)",
                   ["p", "q", "a1", "a2", "b1", "b2"], E, F, gF, gC)


# ----------------------------------------------------------------- lattices / the group G

CLASSES = [(0, 0), (1, 0), (0, 1), (1, 1)]
CLASS_NAME = {(0, 0): "00", (1, 0): "10", (0, 1): "01", (1, 1): "11"}
CHAR_NAME = {(0, 0): "1", (1, 0): "u", (0, 1): "v", (1, 1): "uv"}


def hnf2(vecs):
    """EXACT: Hermite-style canonical basis of the sublattice of Z^2 generated by `vecs`.
    Returns a tuple of 0, 1 or 2 basis vectors in a canonical form, so lattices compare by ==."""
    V = [tuple(int(a) for a in w) for w in vecs if (w[0], w[1]) != (0, 0)]
    if not V:
        return ()
    # integer row reduction to Hermite normal form (2 columns)
    rows = [list(w) for w in V]
    # gcd-reduce first column
    piv = None
    for i, r in enumerate(rows):
        if r[0] != 0:
            piv = i
            break
    if piv is None:
        g = 0
        for r in rows:
            g = np.gcd(g, abs(r[1]))
        return ((0, int(g)),)
    rows[0], rows[piv] = rows[piv], rows[0]
    changed = True
    while changed:
        changed = False
        for i in range(1, len(rows)):
            if rows[i][0] != 0:
                if abs(rows[i][0]) < abs(rows[0][0]):
                    rows[0], rows[i] = rows[i], rows[0]
                q = rows[i][0] // rows[0][0]
                rows[i] = [rows[i][0] - q * rows[0][0], rows[i][1] - q * rows[0][1]]
                if rows[i][0] != 0:
                    changed = True
    if rows[0][0] < 0:
        rows[0] = [-rows[0][0], -rows[0][1]]
    g = 0
    for r in rows[1:]:
        g = int(np.gcd(g, abs(r[1])))
    if g == 0:
        return ((rows[0][0], rows[0][1]),)
    r0 = [rows[0][0], rows[0][1] % g]
    return ((r0[0], r0[1]), (0, g))


def L_S(S):
    """EXACT: canonical basis of the support difference lattice of S (a list of classes)."""
    diffs = [(a - a2, b - b2) for (a, b), (a2, b2) in combinations(S, 2)]
    return hnf2(diffs)


def rank_L(basis):
    return len(basis)


def L_conn_rational(pa, qa, pb, qb):
    """EXACT: the relation lattice L = {(m,n): m*alpha + n*beta = 0 mod 2pi} for
    alpha = 2 pi pa/qa, beta = 2 pi pb/qb.  Returns canonical basis of L (rank 2 always here)."""
    from math import gcd
    q = qa * qb // gcd(qa, qb)
    A = pa * (q // qa)      # alpha = 2 pi A / q
    B = pb * (q // qb)      # beta  = 2 pi B / q
    # L = {(m,n) : m A + n B = 0 mod q}.  Generate by brute force over a fundamental box.
    gens = []
    for m in range(-q, q + 1):
        for n in range(-q, q + 1):
            if (m * A + n * B) % q == 0 and (m, n) != (0, 0):
                gens.append((m, n))
    return hnf2(gens)


def sublattice(basisA, basisB):
    """EXACT: is lattice(basisA) contained in lattice(basisB)?"""
    if not basisA:
        return True
    if not basisB:
        return False
    M = [[b[0] for b in basisB], [b[1] for b in basisB]]
    # solve M z = a over Q, then require z integral
    import itertools
    for a in basisA:
        if len(basisB) == 1:
            b = basisB[0]
            if b[0] != 0:
                if a[0] % b[0] != 0:
                    return False
                t = a[0] // b[0]
            else:
                if a[0] != 0:
                    return False
                if b[1] == 0 or a[1] % b[1] != 0:
                    return False
                t = a[1] // b[1]
            if (t * b[0], t * b[1]) != (a[0], a[1]):
                return False
        else:
            det = basisB[0][0] * basisB[1][1] - basisB[0][1] * basisB[1][0]
            if det == 0:
                return False
            x = (a[0] * basisB[1][1] - a[1] * basisB[1][0])
            y = (basisB[0][0] * a[1] - basisB[0][1] * a[0])
            if x % det != 0 or y % det != 0:
                return False
    return True


# ----------------------------------------------------------------- EXACT Q[sqrt3]

class Q3:
    """EXACT arithmetic in Q + Q*sqrt(3).  Used for q-th roots of unity, q in {1,2,3,4,6,12}."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a * o.a + 3 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__
    __radd__ = __add__

    def __eq__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return self.a == o.a and self.b == o.b

    def __repr__(self):
        return f"({self.a}+{self.b}*sqrt3)"

    def is_rational(self):
        return self.b == 0

    def to_frac(self):
        assert self.b == 0, "not rational"
        return self.a

    def float(self):
        return float(self.a) + float(self.b) * (3 ** 0.5)


_COS12 = {0: Q3(1, 0), 1: Q3(0, Fraction(1, 2)), 2: Q3(Fraction(1, 2), 0), 3: Q3(0, 0),
          4: Q3(Fraction(-1, 2), 0), 5: Q3(0, Fraction(-1, 2)), 6: Q3(-1, 0)}


def cos12(m):
    """EXACT cos(2 pi m / 12) in Q[sqrt3]."""
    m %= 12
    return _COS12[m] if m <= 6 else _COS12[12 - m]


def sin12(m):
    """EXACT sin(2 pi m / 12) in Q[sqrt3]."""
    m %= 12
    return cos12(m - 3)


def root12(m):
    """EXACT (cos, sin) of the 12th root of unity omega^m."""
    return (cos12(m), sin12(m))
