"""
r1lib.py -- LANE W10-A REFUTER 1 (LENS: MATHEMATICS).  Written from scratch for this refutation.

I did NOT import w10a_lib for any construction or any rank computation.  Where I re-run the
lane's own code it is called out explicitly and by module name.

CONVENTIONS (mine; deliberately DIFFERENT indexing from the lane's, so that agreement is
evidence and not a shared bug):
  a complex is (nV, edges, faces, gF, gC) with edges[e] = (tail, head) and faces/loops given
  as signed edge chains {edge: coeff}.
  d1[v,e] = +1 if v = head, -1 if v = tail.
  ranks over Q by Fraction elimination (independent implementation).
  b0 = V - rk d1 ; b1 = (E - rk d1) - rk d2 ; b2 = F - rk d2.
  class(v) = (v in V(gF), v in V(gC)) where V(g) = endpoints of edges in supp(g).
"""
from fractions import Fraction
from itertools import combinations


# ------------------------------------------------------------------ exact linear algebra (mine)
def rk(M):
    """Rank over Q.  M = list of rows (lists of ints/Fractions).  Independent implementation:
    column-major pivot search with explicit row scaling."""
    if not M:
        return 0
    A = [[Fraction(x) for x in r] for r in M]
    nr, nc = len(A), len(A[0])
    row = 0
    for col in range(nc):
        p = -1
        for i in range(row, nr):
            if A[i][col] != 0:
                p = i
                break
        if p < 0:
            continue
        A[row], A[p] = A[p], A[row]
        inv = Fraction(1) / A[row][col]
        A[row] = [inv * x for x in A[row]]
        for i in range(nr):
            if i != row and A[i][col] != 0:
                f = A[i][col]
                A[i] = [x - f * y for x, y in zip(A[i], A[row])]
        row += 1
        if row == nr:
            break
    return row


def in_span(gens, vec):
    """Is vec in the Q-span of gens (all same length)?"""
    if not gens:
        return all(x == 0 for x in vec)
    return rk([list(g) for g in gens]) == rk([list(g) for g in gens] + [list(vec)])


# ------------------------------------------------------------------ complex
class CW:
    def __init__(self, name, nV, edges, faces, gF, gC, vnames=None):
        self.name = name
        self.nV = nV
        self.edges = list(edges)
        self.faces = [dict(f) for f in faces]
        self.gF = dict(gF)
        self.gC = dict(gC)
        self.vnames = vnames or [f"v{i}" for i in range(nV)]

    @property
    def nE(self):
        return len(self.edges)

    @property
    def nF(self):
        return len(self.faces)

    def d1_rows(self):
        M = [[0] * self.nE for _ in range(self.nV)]
        for e, (t, h) in enumerate(self.edges):
            M[t][e] -= 1
            M[h][e] += 1
        return M

    def face_chains(self):
        out = []
        for f in self.faces:
            v = [0] * self.nE
            for e, s in f.items():
                v[e] += s
            out.append(v)
        return out

    def chain(self, g):
        v = [0] * self.nE
        for e, s in g.items():
            v[e] += s
        return v

    def boundary_of(self, chain):
        out = [0] * self.nV
        for e, x in enumerate(chain):
            if x:
                t, h = self.edges[e]
                out[t] -= x
                out[h] += x
        return out

    def report(self):
        d1 = self.d1_rows()
        fc = self.face_chains()
        r1 = rk(d1)
        # rank d2 = rank of the face-chain matrix (rows = faces, cols = edges)
        r2 = rk(fc) if fc else 0
        b0 = self.nV - r1
        b1 = (self.nE - r1) - r2
        b2 = self.nF - r2
        chi = self.nV - self.nE + self.nF
        gauge = self.nV - b0
        inv = self.nE - gauge
        # d1 . d2 == 0
        d1d2 = 0
        for f in fc:
            bd = self.boundary_of(f)
            d1d2 = max(d1d2, max(abs(x) for x in bd) if bd else 0)
        cF, cC = self.chain(self.gF), self.chain(self.gC)
        Fcyc = all(x == 0 for x in self.boundary_of(cF))
        Ccyc = all(x == 0 for x in self.boundary_of(cC))
        Fb = in_span(fc, cF)
        Cb = in_span(fc, cC)
        indep = rk([[cF[e], cC[e]] for e in range(self.nE)]) == 2
        return dict(V=self.nV, E=self.nE, F=self.nF, chi=chi, b0=b0, b1=b1, b2=b2,
                    gauge=gauge, inv=inv, curv=r2, flat=b1, d1d2=d1d2,
                    gF_cycle=Fcyc, gC_cycle=Ccyc, gF_bounds=Fb, gC_bounds=Cb,
                    independent=indep)

    def loop_vertices(self, g):
        S = set()
        for e, s in g.items():
            if s:
                t, h = self.edges[e]
                S.add(t)
                S.add(h)
        return S

    def multiset(self):
        VF, VC = self.loop_vertices(self.gF), self.loop_vertices(self.gC)
        ms = {(0, 0): 0, (1, 0): 0, (0, 1): 0, (1, 1): 0}
        for v in range(self.nV):
            ms[(1 if v in VF else 0, 1 if v in VC else 0)] += 1
        return ms

    def pi_U(self):
        ms = self.multiset()
        n = Fraction(1, self.nV)
        return [ms[(0, 0)] * n, ms[(1, 0)] * n, ms[(0, 1)] * n, ms[(1, 1)] * n]


# ------------------------------------------------------------------ my carriers
def my_K1():
    E = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    return CW("K1 (S1 sec1 verbatim)", 5, E, [{0: 1, 1: 1, 2: 1}],
              {0: 1, 1: 1, 2: 1}, {3: 1, 4: 1, 5: 1})


def my_B0b():
    """3x3 grid torus, MY indexing: vertex (i,j) -> i*3 + j  (the lane uses 3*j+i).
    Edge order: all VERTICAL first, then all HORIZONTAL (the lane does the opposite).
    gamma_F = boundary of the square with corners (0,0),(1,0),(1,1),(0,1).
    gamma_C = the row j = 0 traversed in i  ->  meets gamma_F in 2 vertices."""
    def V(i, j):
        return (i % 3) * 3 + (j % 3)
    E = []
    _idx = {}
    for i in range(3):                       # vertical edges (i,j)->(i,j+1)  FIRST
        for j in range(3):
            _idx[('v', i, j)] = len(E)
            E.append((V(i, j), V(i, j + 1)))
    for i in range(3):                       # horizontal edges (i,j)->(i+1,j)
        for j in range(3):
            _idx[('h', i, j)] = len(E)
            E.append((V(i, j), V(i + 1, j)))

    def idx(kind, i, j):
        return _idx[(kind, i % 3, j % 3)]

    faces = []
    for i in range(3):
        for j in range(3):
            faces.append({idx('h', i, j): 1, idx('v', i + 1, j): 1,
                          idx('h', i, j + 1): -1, idx('v', i, j): -1})
    gF = {idx('h', 0, 0): 1, idx('v', 1, 0): 1, idx('h', 0, 1): -1, idx('v', 0, 0): -1}
    gC = {idx('h', 0, 0): 1, idx('h', 1, 0): 1, idx('h', 2, 0): 1}
    names = [f"({i},{j})" for i in range(3) for j in range(3)]
    return CW("B0b 3x3 torus, loops meet (MY indexing)", 9, E, faces, gF, gC, names), _idx


def my_B4_square():
    """The lane's B4 shape, rebuilt with MY labels and a DIFFERENT edge order:
    two spheres, each a SQUARE with two 2-cells, glued at two OPPOSITE corners."""
    #  0=p 1=q 2=a1 3=a2 4=b1 5=b2 ; edge order interleaved A/B, not blocked
    E = [(0, 2), (0, 4), (2, 1), (4, 1), (1, 3), (1, 5), (3, 0), (5, 0)]
    sqA = {0: 1, 2: 1, 4: 1, 6: 1}      # p->a1->q->a2->p
    sqB = {1: 1, 3: 1, 5: 1, 7: 1}      # p->b1->q->b2->p
    faces = [dict(sqA), dict(sqA), dict(sqB), dict(sqB)]
    gF = dict(sqA)
    gC = {0: 1, 2: 1, 3: -1, 1: -1}     # p->a1->q->b1->p
    return CW("B4 spindle, SQUARE|SQUARE at opposite corners (MY labels)", 6, E, faces, gF, gC,
              ["p", "q", "a1", "a2", "b1", "b2"])


def my_B4_tri_pent():
    """AN ALTERNATIVE SPINDLE, ALSO 'two 2-spheres glued at two points', ALSO V=6 E=8 F=4.
    Sphere A = a TRIANGLE p,q,r with two 2-cells (v=3,e=3,f=2, chi=2).
    Sphere B = a PENTAGON p,s1,q,s2,s3 with two 2-cells (v=5,e=5,f=2, chi=2).
    Glued at the two points p and q.  V = 3+5-2 = 6, E = 3+5 = 8, F = 2+2 = 4."""
    #  0=p 1=q 2=r 3=s1 4=s2 5=s3
    E = [(0, 1), (1, 2), (2, 0),                 # 0,1,2 triangle p->q->r->p
         (0, 3), (3, 1), (1, 4), (4, 5), (5, 0)]  # 3..7 pentagon p->s1->q->s2->s3->p
    tri = {0: 1, 1: 1, 2: 1}
    pen = {3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
    faces = [dict(tri), dict(tri), dict(pen), dict(pen)]
    gF = dict(tri)
    gC = {0: 1, 5: 1, 6: 1, 7: 1}     # p->q (triangle edge) then q->s2->s3->p (pentagon)
    return CW("B4' spindle, TRIANGLE|PENTAGON at 2 points", 6, E, faces, gF, gC,
              ["p", "q", "r", "s1", "s2", "s3"])


def my_B1q():
    """B1q: two triangles joined by a 2-edge bridge; the bridge midpoint is the spectator."""
    E = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 6), (6, 3)]
    return CW("B1q K1-bridged + spectator", 7, E, [{0: 1, 1: 1, 2: 1}],
              {0: 1, 1: 1, 2: 1}, {3: 1, 4: 1, 5: 1},
              ["v0", "v1", "v2", "w0", "w1", "w2", "m"])


# ------------------------------------------------------------------ simple-cycle enumeration
def simple_cycles(cw):
    """All simple cycles of the underlying multigraph, as signed edge chains.
    Returns list of (chain dict, frozenset of vertices)."""
    adj = {}
    for e, (t, h) in enumerate(cw.edges):
        adj.setdefault(t, []).append((e, h, +1))
        adj.setdefault(h, []).append((e, t, -1))
    found = {}

    def walk(start, cur, path_edges, path_signs, visited):
        for (e, nxt, s) in adj.get(cur, []):
            if e in path_edges:
                continue
            if nxt == start and len(path_edges) >= 1:
                if len(path_edges) < 2:      # need at least 2 edges for a cycle (no self loops)
                    continue
                ch = {}
                for ee, ss in zip(path_edges + [e], path_signs + [s]):
                    ch[ee] = ch.get(ee, 0) + ss
                key = tuple(sorted(path_edges + [e]))
                if key not in found:
                    vs = set()
                    for ee in ch:
                        vs.add(cw.edges[ee][0])
                        vs.add(cw.edges[ee][1])
                    found[key] = (ch, frozenset(vs))
                continue
            if nxt in visited:
                continue
            walk(start, nxt, path_edges + [e], path_signs + [s], visited | {nxt})

    for v0 in range(cw.nV):
        walk(v0, v0, [], [], {v0})
    return list(found.values())
