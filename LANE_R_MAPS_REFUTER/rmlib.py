"""
LANE R-MAPS REFUTER -- library, written from scratch, numpy only (+ mpmath for one
independent cross-check of the Mahler measure).  No corpus code reused.

CONVENTIONS PUBLISHED HERE (S4's failure to publish these is a defect of record):

  d1 : |V| x |E| integer matrix.  d1[v,e] = +1 if v is the TARGET of e,
                                          -1 if v is the SOURCE of e,
                                           0 otherwise.
       (boundary of an edge = target - source; a loop edge gives a zero column)
  d2 : |E| x |F| integer matrix.  d2[e,F] = signed multiplicity of e in the
       attaching word of F (+1 traversed forward, -1 backward).
  Chain complex over Q:  C2 --d2--> C1 --d1--> C0.
  b0 = |V| - rank d1 ;  b1 = |E| - rank d1 - rank d2 ;  b2 = |F| - rank d2.
  Ranks computed by EXACT Gaussian elimination over the rationals (fractions.Fraction),
  never by floating-point SVD.

  Fibre: L_v = C, rank one, <z,w> = conj(z) w   (S1 sec.3).
  Connection: a_e in R/2piZ, U_e = exp(i a_e).
  Loop gamma given as a list of (edge_index, sign); holonomy W = exp(i * sum sign*a_e).
  Vertex incidence a_v = [v lies on gamma_F], b_v = [v lies on gamma_C].
  M_gamma = diagonal matrix, entry W(gamma) on vertices of gamma, 1 elsewhere (W-01/S4 sec.2).
  Ready section s with s_v = sqrt(p_v) >= 0, sum p_v = 1.
  Z_k = < M_F^k s , M_C^k s >  computed BY EXPLICIT MATRIX POWERS on C^V.
  u = conj(W_F), v = W_C.  chi_(a,b) = u^a v^b.
  pi = (p00,p10,p01,p11) = pushforward of p under v |-> (a_v,b_v).
  S = { classes with pi > 0 }.  G = < chi_a/chi_b : a,b in S >.
  rank G = rank over Z of the difference lattice D = { a-b : a,b in S } < Z^2
           (this equals the rank of G for u,v multiplicatively independent).
  Schedule A: k_n = 1, lambda_A = log|Z_1|.
  Schedule B: k_n = n, lambda_B = lim (1/N) sum_{n=1..N} log|Z_n|.
"""

import numpy as np
from fractions import Fraction

# ---------------------------------------------------------------- exact rank

def rank_exact(M):
    """Exact rank over Q of an integer matrix, by fraction-free-ish Gaussian elim."""
    A = [[Fraction(int(x)) for x in row] for row in np.asarray(M, dtype=object)]
    if not A or not A[0]:
        return 0
    nr, nc = len(A), len(A[0])
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(nc)]
        r += 1
        if r == nr:
            break
    return r


# ---------------------------------------------------------------- carrier

class Carrier:
    """(CW complex of dim <=2, designated loops gamma_F and gamma_C, vertex weights p)."""

    def __init__(self, name, verts, edges, faces, loopF, loopC, p):
        """
        verts : list of vertex names
        edges : list of (edge_name, src_name, tgt_name)
        faces : list of (face_name, [(edge_name, sign), ...])   attaching word
        loopF, loopC : [(edge_name, sign), ...]  closed edge-paths
        p     : dict vertex_name -> weight (need not be normalised; normalised here)
        """
        self.name = name
        self.verts = list(verts)
        self.edges = list(edges)
        self.faces = list(faces)
        self.vi = {v: i for i, v in enumerate(self.verts)}
        self.ei = {e[0]: i for i, e in enumerate(self.edges)}
        self.loopF = list(loopF)
        self.loopC = list(loopC)
        tot = sum(p[v] for v in self.verts)
        self.p = np.array([float(p[v]) / tot for v in self.verts])

        V, E, F = len(self.verts), len(self.edges), len(self.faces)
        d1 = np.zeros((V, E), dtype=int)
        for j, (en, s, t) in enumerate(self.edges):
            d1[self.vi[t], j] += 1
            d1[self.vi[s], j] -= 1
        d2 = np.zeros((E, F), dtype=int)
        for j, (fn, word) in enumerate(self.faces):
            for (en, sg) in word:
                d2[self.ei[en], j] += sg
        self.d1, self.d2 = d1, d2

    # -- topology
    def betti(self):
        V, E, F = len(self.verts), len(self.edges), len(self.faces)
        r1, r2 = rank_exact(self.d1), rank_exact(self.d2)
        return (V - r1, E - r1 - r2, F - r2), (r1, r2)

    def chi(self):
        return len(self.verts) - len(self.edges) + len(self.faces)

    def check_d2_zero(self):
        """d1 . d2 must be the zero matrix."""
        return int(np.abs(self.d1 @ self.d2).max()) if self.d2.size else 0

    def loop_is_closed(self, loop):
        """A closed edge-path: its 1-chain must be a cycle (d1 . chain = 0)."""
        ch = np.zeros(len(self.edges), dtype=int)
        for (en, sg) in loop:
            ch[self.ei[en]] += sg
        return int(np.abs(self.d1 @ ch).max())

    def loop_verts(self, loop):
        s = set()
        for (en, sg) in loop:
            _, a, b = self.edges[self.ei[en]]
            s.add(a); s.add(b)
        return s

    # -- incidence / classes
    def incidence(self):
        VF, VC = self.loop_verts(self.loopF), self.loop_verts(self.loopC)
        a = np.array([1 if v in VF else 0 for v in self.verts])
        b = np.array([1 if v in VC else 0 for v in self.verts])
        return a, b

    def pi(self):
        a, b = self.incidence()
        out = np.zeros(4)          # order: 00, 10, 01, 11
        idx = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}
        for i in range(len(self.verts)):
            out[idx[(int(a[i]), int(b[i]))]] += self.p[i]
        return out

    def holon(self, aedge):
        """aedge : dict edge_name -> a_e ; returns (W_F, W_C)."""
        fF = sum(sg * aedge[en] for (en, sg) in self.loopF)
        fC = sum(sg * aedge[en] for (en, sg) in self.loopC)
        return np.exp(1j * fF), np.exp(1j * fC)

    # -- transport, by EXPLICIT matrices
    def Z_matrix(self, WF, WC, k):
        a, b = self.incidence()
        MF = np.diag(np.where(a == 1, WF, 1.0 + 0j))
        MC = np.diag(np.where(b == 1, WC, 1.0 + 0j))
        s = np.sqrt(self.p).astype(complex)
        x = np.linalg.matrix_power(MF, k) @ s
        y = np.linalg.matrix_power(MC, k) @ s
        return np.vdot(x, y)          # np.vdot conjugates the FIRST argument


# ------------------------------------------------- class-level (vectorised) Z

def Z_from_pi(pi, f, c, k):
    """Z_k from the class weights only.  u = e^{-if}, v = e^{ic}."""
    u = np.exp(-1j * f); v = np.exp(1j * c)
    return pi[0] + pi[1] * u**k + pi[2] * v**k + pi[3] * (u * v)**k


def lambda_A(pi, f, c):
    return float(np.log(np.abs(Z_from_pi(pi, f, c, 1))))


def lambda_B_direct(pi, f, c, N, chunk=200000):
    """(1/N) sum_{n=1..N} log|Z_n| with k_n = n.  Vectorised, chunked."""
    tot = 0.0
    n0 = 1
    while n0 <= N:
        n1 = min(N, n0 + chunk - 1)
        n = np.arange(n0, n1 + 1, dtype=np.float64)
        z = (pi[0] + pi[1] * np.exp(-1j * f * n) + pi[2] * np.exp(1j * c * n)
             + pi[3] * np.exp(1j * (c - f) * n))
        tot += float(np.sum(np.log(np.abs(z))))
        n0 = n1 + 1
    return tot / N


def mahler2(pi, n=3000, seed=None):
    """m(p00 + p10 x + p01 y + p11 xy) by a 2D grid average on T^2.

    GRID CONVENTION: theta_j = 2 pi (j + 1/2)/n, phi_l = 2 pi (l + 1/2)/n,
    j,l = 0..n-1  (MIDPOINT rule -- deliberately avoids the lattice point
    (0,0) where the polynomial can be large, and avoids all rational
    resonances of the form theta = phi exactly).  Chunked over j.
    """
    th = 2 * np.pi * (np.arange(n) + 0.5) / n
    ph = 2 * np.pi * (np.arange(n) + 0.5) / n
    ey = np.exp(1j * ph)
    tot = 0.0
    for j0 in range(0, n, 200):
        j1 = min(n, j0 + 200)
        ex = np.exp(1j * th[j0:j1])[:, None]
        z = pi[0] + pi[1] * ex + pi[2] * ey[None, :] + pi[3] * ex * ey[None, :]
        tot += float(np.sum(np.log(np.abs(z))))
    return tot / (n * n)


def mahler2_mc(pi, N=20_000_000, seed=20260816):
    """Monte-Carlo m(P) on T^2.  SEED PUBLISHED: 20260816, numpy default_rng."""
    rng = np.random.default_rng(seed)
    tot = 0.0; done = 0
    while done < N:
        b = min(2_000_000, N - done)
        th = rng.uniform(0, 2 * np.pi, b); ph = rng.uniform(0, 2 * np.pi, b)
        ex = np.exp(1j * th); ey = np.exp(1j * ph)
        z = pi[0] + pi[1] * ex + pi[2] * ey + pi[3] * ex * ey
        tot += float(np.sum(np.log(np.abs(z))))
        done += b
    return tot / N


def cassaigne_maillot(a, b, c, dps=40):
    """IMPORT (see IMPORT AUDIT): Cassaigne-Maillot closed form for m(a+bx+cy),
    a,b,c >= 0.  Used ONLY as a third independent cross-check."""
    import mpmath as mp
    mp.mp.dps = dps
    a, b, c = [mp.mpf(t) for t in (a, b, c)]
    if a + b <= c or b + c <= a or c + a <= b:
        return float(mp.log(max(a, b, c)))
    # angles opposite the sides a, b, c
    al = mp.acos((b**2 + c**2 - a**2) / (2 * b * c))
    be = mp.acos((c**2 + a**2 - b**2) / (2 * c * a))
    ga = mp.acos((a**2 + b**2 - c**2) / (2 * a * b))
    def D(z):
        return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))
    val = (al * mp.log(a) + be * mp.log(b) + ga * mp.log(c)
           + D((a / b) * mp.exp(1j * ga))) / mp.pi
    return float(val)


def rank_G(pi, tol=0.0):
    """|S| and rank over Z of the difference lattice of the occupied classes."""
    cls = [(0, 0), (1, 0), (0, 1), (1, 1)]
    S = [cls[i] for i in range(4) if pi[i] > tol]
    if len(S) <= 1:
        return len(S), 0, S
    D = []
    for A in S:
        for B in S:
            d = (A[0] - B[0], A[1] - B[1])
            if d != (0, 0):
                D.append(d)
    return len(S), rank_exact(np.array(D, dtype=int)), S


def forms(pi, tol=0.0):
    """FORMATION OCCURS <=> G != {1}.  For u,v multiplicatively independent this
    is exactly rank(difference lattice) >= 1, i.e. |S| >= 2."""
    _, r, _ = rank_G(pi, tol)
    return r >= 1


# ---------------------------------------------------------------- K1 itself

def K1(p=None, name="K1"):
    verts = ["v0", "v1", "v2", "v3", "v4"]
    edges = [("e1", "v0", "v1"), ("e2", "v1", "v2"), ("e3", "v2", "v0"),
             ("e4", "v0", "v3"), ("e5", "v3", "v4"), ("e6", "v4", "v0")]
    faces = [("F", [("e1", +1), ("e2", +1), ("e3", +1)])]
    loopF = [("e1", +1), ("e2", +1), ("e3", +1)]
    loopC = [("e4", +1), ("e5", +1), ("e6", +1)]
    if p is None:
        p = {v: 1.0 for v in verts}          # SENSE U: uniform on vertices
    return Carrier(name, verts, edges, faces, loopF, loopC, p)


def contract(K, edge_name, newname=None):
    """Contract one edge (a contractible subcomplex): a homotopy equivalence.
    Weights push forward by SUMMATION (the only pushforward of a measure)."""
    _, s, t = K.edges[K.ei[edge_name]]
    if s == t:
        raise ValueError("loop edge cannot be contracted")
    merged = s + t
    rel = {s: merged, t: merged}
    def m(v): return rel.get(v, v)
    verts = []
    for v in K.verts:
        if m(v) not in verts:
            verts.append(m(v))
    edges = [(en, m(a), m(b)) for (en, a, b) in K.edges if en != edge_name]
    def strip(word): return [(en, sg) for (en, sg) in word if en != edge_name]
    faces = [(fn, strip(w)) for (fn, w) in K.faces]
    p = {v: 0.0 for v in verts}
    for i, v in enumerate(K.verts):
        p[m(v)] += K.p[i]
    return Carrier(newname or (K.name + "/" + edge_name), verts, edges, faces,
                   strip(K.loopF), strip(K.loopC), p)


def collapse_tree(K, tree_edges, newname=None):
    """Quotient by a subtree of the 1-skeleton: a homotopy equivalence."""
    out = K
    for e in tree_edges:
        out = contract(out, e)
    out.name = newname or (K.name + "/T")
    return out
