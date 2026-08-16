"""
rm_lib.py  --  LANE R (MAPS REFUTER), "Where Atoms Come From".
Written from scratch for this lane.  Reuses NO code from the maps axis lane
(run3.py), from S4, or from any other lane in this project.

PUBLISHED CONVENTIONS (S4 published none of this; this lane does).

  A carrier is a triple (V, d1, d2) plus two designated 1-cycles gamma_F, gamma_C
  given as EDGE SUPPORT SETS (index sets into the edge list).

  d1 : R^E -> R^V  is the boundary map of 1-cells.  Column e of d1 for an edge
       e = (i -> j) has  d1[j,e] = +1, d1[i,e] = -1.   (head minus tail.)
  d2 : R^F -> R^E  is the boundary map of 2-cells.  Column f holds the signed
       edge multiplicities of the attaching cycle of face f, sign +1 if the face
       traverses the edge along its stored orientation, -1 against.

  b0 = V - rank(d1);  b1 = (E - rank(d1)) - rank(d2);  b2 = F - rank(d2).
  Ranks over Q, computed by exact integer row reduction (Fraction), NOT by SVD.

  VERTEX CLASS.  a_v = [v is an endpoint of some edge of gamma_F],
                 b_v = [v is an endpoint of some edge of gamma_C].
                 class(v) = (a_v, b_v) in {0,1}^2.

  CLASS PUSHFORWARD.  pi = (p00, p10, p01, p11), p_ab = sum of p_v over class(v)=(a,b).

  TRANSPORT.  u = e^{i f}, v = e^{i c}   (see CONVENTION NOTE below)
              Z_k = p00 + p10 u^k + p01 v^k + p11 (uv)^k
              P(x,y) = p00 + p10 x + p01 y + p11 x y            (the Z-polynomial)

  CONVENTION NOTE (declared, not hidden).  S4 :155 writes the expansion of |Z_k|^2
  with a cross term 2 q r cos(k(f+c)); the convention above yields 2 q r cos(k(f-c)).
  The two differ exactly by c -> -c.  Every quantity this lane computes (Mahler
  measure, lambda_B, |Z_k| moduli, the relation-lattice tier) is invariant under
  y -> y^{-1}, so the difference is immaterial; it is verified numerically in
  rm_1_validate.py rather than assumed.

  SCHEDULES.  Schedule A: k_n = 1, lambda_A = log|Z_1|.
              Schedule B (canonical clock): k_n = n,
              lambda_B = lim_N (1/N) sum_{n=1..N} log|Z_n|.
              For (f,c) with trivial relation lattice L = {(m,n): u^m v^n = 1} = 0,
              Weyl equidistribution on T^2 gives lambda_B = m(P), the 2-variable
              logarithmic Mahler measure.

  SENSE U: p_v = 1/V (uniform on vertices).      SENSE C: class weights fixed by hand.

  SEEDS.  Every random draw in this lane uses numpy Generator(PCG64(seed)) with the
  seed printed at the point of use.  No global np.random state is used anywhere.
"""

import numpy as np
from fractions import Fraction
import math

# ---------------------------------------------------------------- exact rank

def exact_rank(M):
    """Rank over Q by fraction-free Gaussian elimination. M is a list of lists of int."""
    A = [[Fraction(x) for x in row] for row in M]
    rows, cols = len(A), (len(A[0]) if A else 0)
    r = 0
    for cidx in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i][cidx] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][cidx]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][cidx] != 0:
                f = A[i][cidx]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r


# ---------------------------------------------------------------- carrier

class Carrier:
    def __init__(self, name, nV, edges, faces, gF, gC, vnames=None):
        """
        edges : list of (tail, head) vertex indices
        faces : list of list of (edge_index, sign)
        gF, gC: list of edge indices (the designated loops' edge supports)
        """
        self.name = name
        self.nV = nV
        self.edges = list(edges)
        self.faces = list(faces)
        self.gF = sorted(gF)
        self.gC = sorted(gC)
        self.vnames = vnames if vnames else [f"v{i}" for i in range(nV)]

    @property
    def nE(self):
        return len(self.edges)

    @property
    def nF(self):
        return len(self.faces)

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

    def d1d2(self):
        A = np.array(self.d1(), dtype=np.int64)
        B = np.array(self.d2(), dtype=np.int64)
        if self.nF == 0:
            return np.zeros((self.nV, 0), dtype=np.int64)
        return A @ B

    def betti(self):
        r1 = exact_rank(self.d1()) if self.nE else 0
        r2 = exact_rank(self.d2()) if self.nF else 0
        b0 = self.nV - r1
        b1 = (self.nE - r1) - r2
        b2 = self.nF - r2
        return b0, b1, b2, r1, r2

    def chi(self):
        return self.nV - self.nE + self.nF

    def loop_vertices(self, edge_set):
        S = set()
        for e in edge_set:
            t, h = self.edges[e]
            S.add(t)
            S.add(h)
        return S

    def classes(self):
        """Return list of (a_v,b_v) per vertex."""
        VF = self.loop_vertices(self.gF)
        VC = self.loop_vertices(self.gC)
        return [(1 if v in VF else 0, 1 if v in VC else 0) for v in range(self.nV)]

    def is_cycle(self, edge_set):
        """Check the edge set is a 1-cycle: d1 applied to its indicator (with some
        orientation choice) is zero.  We test the unsigned closed-walk condition:
        every vertex has even degree in the sub-multigraph."""
        deg = {}
        for e in edge_set:
            t, h = self.edges[e]
            deg[t] = deg.get(t, 0) + 1
            deg[h] = deg.get(h, 0) + 1
        return all(d % 2 == 0 for d in deg.values())

    def pi_uniform(self):
        """SENSE U class pushforward."""
        return self.pi_from_p([1.0 / self.nV] * self.nV)

    def pi_from_p(self, p):
        cl = self.classes()
        out = {(0, 0): 0.0, (1, 0): 0.0, (0, 1): 0.0, (1, 1): 0.0}
        for v, w in enumerate(p):
            out[cl[v]] += w
        return (out[(0, 0)], out[(1, 0)], out[(0, 1)], out[(1, 1)])

    def class_counts(self):
        cl = self.classes()
        out = {}
        for c in cl:
            out[c] = out.get(c, 0) + 1
        return out


# ---------------------------------------------------------------- transport

def Z_from_vertices(carrier, p, f, c, k):
    """Z_k computed vertex-by-vertex on the carrier (the definition, not the class form)."""
    cl = carrier.classes()
    tot = 0j
    for v in range(carrier.nV):
        a, b = cl[v]
        tot += p[v] * np.exp(1j * (k * (a * f + b * c)))
    return tot


def Z_from_pi(pi, f, c, k):
    """Z_k computed from the class pushforward alone."""
    p00, p10, p01, p11 = pi
    u = np.exp(1j * k * f)
    v = np.exp(1j * k * c)
    return p00 + p10 * u + p01 * v + p11 * u * v


# ---------------------------------------------------------------- Mahler

def mahler_jensen(pi, n=200000):
    """
    m(P), P = p00 + p10 x + p01 y + p11 x y, by Jensen reduction in x:
        P = (p00 + p01 y) + (p10 + p11 y) x
        m(P) = (1/2pi) int_0^{2pi} log max(|p00+p01 e^{it}|, |p10+p11 e^{it}|) dt
    Trapezoid on a uniform grid; the integrand is periodic so trapezoid is
    spectrally accurate away from the (measure-zero) branch points.
    """
    p00, p10, p01, p11 = pi
    t = 2.0 * np.pi * np.arange(n) / n          # uniform, periodic -> rectangle == trapezoid
    e = np.exp(1j * t)
    A = np.abs(p00 + p01 * e)
    B = np.abs(p10 + p11 * e)
    M = np.maximum(A, B)
    with np.errstate(divide='ignore'):
        L = np.log(M)
    return float(np.mean(L))


def mahler_2d(pi, n=3000):
    """Independent re-derivation: brute 2-D quadrature of (1/(2pi)^2) int int log|P|."""
    p00, p10, p01, p11 = pi
    t = 2.0 * np.pi * np.arange(n) / n
    X = np.exp(1j * t)[:, None]
    Y = np.exp(1j * t)[None, :]
    P = p00 + p10 * X + p01 * Y + p11 * X * Y
    with np.errstate(divide='ignore'):
        L = np.log(np.abs(P))
    return float(np.mean(L))


def bloch_wigner(z):
    """D(z) = Im Li2(z) + arg(1-z) log|z|, via mpmath's polylog."""
    import mpmath as mp
    mp.mp.dps = 30
    zz = mp.mpc(z)
    return float(mp.im(mp.polylog(2, zz)) + mp.arg(1 - zz) * mp.log(abs(zz)))


def mahler_cassaigne_maillot(a, b, c):
    """
    Exact closed form for m(a + b x + c y), a,b,c > 0 (Cassaigne-Maillot 2000).
    IMPORT -- see IMPORT AUDIT.
    If a,b,c satisfy the triangle inequality:
        pi*m = D(|a/b| e^{i gamma}) + alpha log a + beta log b + gamma log c
    with alpha,beta,gamma the angles opposite sides a,b,c.
    Else m = log max(a,b,c).
    """
    a, b, c = float(a), float(b), float(c)
    if not (a < b + c and b < a + c and c < a + b):
        return math.log(max(a, b, c))
    def ang(o, p, q):   # angle opposite o
        return math.acos(max(-1.0, min(1.0, (p * p + q * q - o * o) / (2 * p * q))))
    alpha = ang(a, b, c)
    beta = ang(b, a, c)
    gamma = ang(c, a, b)
    D = bloch_wigner((a / b) * complex(math.cos(gamma), math.sin(gamma)))
    return (D + alpha * math.log(a) + beta * math.log(b) + gamma * math.log(c)) / math.pi


def lambda_B_closed(pi):
    """
    lambda_B at generic (f,c) (relation lattice L = 0) = m(P).
    Uses the exact Cassaigne-Maillot form when only three class weights are
    positive (after clearing a monomial), else the Jensen quadrature.
    """
    nz = [w for w in pi if w > 0]
    if len(nz) == 0:
        return float('-inf')
    if len(nz) == 1:
        return math.log(nz[0])
    if len(nz) == 2:
        return math.log(max(nz))          # Jensen; a binomial after monomial clearing
    if len(nz) == 3:
        # any 3 of the 4 corners of the unit square are affinely independent, so a
        # GL_2(Z) monomial substitution carries them to {1, x, y}: m depends only on
        # the multiset of the three coefficients.
        return mahler_cassaigne_maillot(*sorted(nz))
    return mahler_jensen(pi)


def lambda_B_direct(pi, f, c, N=2000000):
    """Schedule B by direct simulation: (1/N) sum_{n=1..N} log|Z_n|.  No closed form."""
    n = np.arange(1, N + 1, dtype=np.float64)
    p00, p10, p01, p11 = pi
    uf = n * f
    vc = n * c
    Z = p00 + p10 * np.exp(1j * uf) + p01 * np.exp(1j * vc) + p11 * np.exp(1j * (uf + vc))
    return float(np.mean(np.log(np.abs(Z))))


def relation_lattice_rank(f, c, tol=1e-9, K=60):
    """rank of L = {(m,n): m f + n c = 0 mod 2pi}, probed on |m|,|n| <= K."""
    gens = []
    for m in range(-K, K + 1):
        for nn in range(-K, K + 1):
            if m == 0 and nn == 0:
                continue
            r = (m * f + nn * c) / (2 * np.pi)
            if abs(r - round(r)) < tol:
                gens.append((m, nn))
    if not gens:
        return 0, []
    return exact_rank([[g[0], g[1]] for g in gens]), gens[:6]


# ---------------------------------------------------------------- carriers

def K1():
    """B1 = K1 as handed.  v0 root; filled triangle v0,v1,v2; unfilled v0,v3,v4."""
    E = [(0, 1), (1, 2), (2, 0),        # e1,e2,e3  filled triangle  (face F)
         (0, 3), (3, 4), (4, 0)]        # e4,e5,e6  unfilled triangle
    F = [[(0, 1), (1, 1), (2, 1)]]      # the one 2-cell, boundary e1+e2+e3
    return Carrier("B1  K1", 5, E, F, gF=[0, 1, 2], gC=[3, 4, 5],
                   vnames=["v0", "v1", "v2", "v3", "v4"])


def K1_subdivided():
    """B1s = K1 with EVERY edge subdivided once.  11 vertices, 12 edges."""
    # 0..4 original; 5..10 midpoints of e1..e6
    E = [(0, 5), (5, 1),      # e1 -> 0-m1-1
         (1, 6), (6, 2),      # e2 -> 1-m2-2
         (2, 7), (7, 0),      # e3 -> 2-m3-0
         (0, 8), (8, 3),      # e4 -> 0-m4-3
         (3, 9), (9, 4),      # e5 -> 3-m5-4
         (4, 10), (10, 0)]    # e6 -> 4-m6-0
    F = [[(i, 1) for i in range(6)]]                 # subdivided filled triangle
    return Carrier("B1s K1 subdivided", 11, E, F, gF=list(range(6)), gC=list(range(6, 12)),
                   vnames=["v0", "v1", "v2", "v3", "v4",
                           "m1", "m2", "m3", "m4", "m5", "m6"])


def K1_partial_subdiv(nF_sub, nC_sub):
    """
    K1 with nF_sub extra vertices inserted into the FILLED triangle's boundary and
    nC_sub extra vertices inserted into the UNFILLED triangle's boundary.
    Homeomorphic to K1 for every (nF_sub, nC_sub).
    Built as: filled loop = closed walk v0 - (nF_sub+2 others) - v0 ; ditto unfilled.
    """
    # filled loop vertices: v0, then a chain of (2 + nF_sub) vertices back to v0
    nV = 1
    Fchain = []
    for _ in range(2 + nF_sub):
        Fchain.append(nV); nV += 1
    Cchain = []
    for _ in range(2 + nC_sub):
        Cchain.append(nV); nV += 1
    E = []
    fl = [0] + Fchain + [0]
    gF = []
    for i in range(len(fl) - 1):
        gF.append(len(E)); E.append((fl[i], fl[i + 1]))
    cl = [0] + Cchain + [0]
    gC = []
    for i in range(len(cl) - 1):
        gC.append(len(E)); E.append((cl[i], cl[i + 1]))
    F = [[(e, 1) for e in gF]]
    return Carrier(f"K1[nF={nF_sub},nC={nC_sub}]", nV, E, F, gF=gF, gC=gC)


def K1_bridged():
    """B1p: two triangles joined by a BRIDGE edge (not wedged).  6 vertices, 7 edges."""
    # triangle A: 0,1,2 (filled).  triangle B: 3,4,5 (unfilled).  bridge 0-3.
    E = [(0, 1), (1, 2), (2, 0),      # 0,1,2 filled triangle
         (3, 4), (4, 5), (5, 3),      # 3,4,5 unfilled triangle
         (0, 3)]                      # 6 bridge
    F = [[(0, 1), (1, 1), (2, 1)]]
    return Carrier("B1p K1-bridged", 6, E, F, gF=[0, 1, 2], gC=[3, 4, 5])


def K1_bridged_subdiv():
    """B1q: B1p with the BRIDGE subdivided once.  7 vertices, 8 edges.
    Homeomorphic to B1p (a subdivision of a 1-cell)."""
    E = [(0, 1), (1, 2), (2, 0),
         (3, 4), (4, 5), (5, 3),
         (0, 6), (6, 3)]              # bridge 0-m-3, m = vertex 6
    F = [[(0, 1), (1, 1), (2, 1)]]
    return Carrier("B1q K1-bridged+spectator", 7, E, F, gF=[0, 1, 2], gC=[3, 4, 5],
                   vnames=["v0", "v1", "v2", "w0", "w1", "w2", "m"])


def torus33(loops_meet):
    """
    B0a / B0b: the 3x3 grid ring torus.  ONE AND THE SAME COMPLEX for both;
    only the DESIGNATION of gamma_C differs.  9 vertices, 18 edges, 9 faces.
    vertex (i,j) -> index 3*j+i, i,j in Z_3.
    horizontal edge h(i,j): (i,j) -> (i+1,j)          index      3*j+i
    vertical   edge v(i,j): (i,j) -> (i,j+1)          index  9 + 3*j+i
    face (i,j) boundary: h(i,j) + v(i+1,j) - h(i,j+1) - v(i,j)
    gamma_F = boundary of face (0,0)   (4 vertices: (0,0),(1,0),(1,1),(0,1))
    gamma_C = row j=2  (disjoint from gamma_F)   if loops_meet is False
    gamma_C = row j=0  (meets gamma_F in 2 vertices) if loops_meet is True
    """
    def V(i, j):
        return 3 * (j % 3) + (i % 3)
    def H(i, j):
        return 3 * (j % 3) + (i % 3)
    def Vt(i, j):
        return 9 + 3 * (j % 3) + (i % 3)
    E = [None] * 18
    for j in range(3):
        for i in range(3):
            E[H(i, j)] = (V(i, j), V(i + 1, j))
            E[Vt(i, j)] = (V(i, j), V(i, j + 1))
    Fc = []
    for j in range(3):
        for i in range(3):
            Fc.append([(H(i, j), 1), (Vt(i + 1, j), 1), (H(i, j + 1), -1), (Vt(i, j), -1)])
    gF = [H(0, 0), Vt(1, 0), H(0, 1), Vt(0, 0)]        # the square's four edges
    row = 0 if loops_meet else 2
    gC = [H(0, row), H(1, row), H(2, row)]             # a horizontal cycle
    nm = "B0b torus, loops MEET" if loops_meet else "B0a torus, loops DISJOINT"
    return Carrier(nm, 9, E, Fc, gF=gF, gC=gC,
                   vnames=[f"({i},{j})" for j in range(3) for i in range(3)])


ALL = {
    "B1": K1, "B1s": K1_subdivided, "B1p": K1_bridged,
    "B1q": K1_bridged_subdiv,
    "B0a": (lambda: torus33(False)), "B0b": (lambda: torus33(True)),
}
