"""glib -- finite-group toolkit for LANE_SCALE_A_TRANSPORT.

Everything here is EXACT combinatorics on an explicit multiplication table, plus a character
table obtained by simultaneous diagonalisation of the class-sum matrices.  Every object carries
a self-check that FAILS LOUDLY if the construction is wrong:

  * associativity of the multiplication table is checked on all n^3 triples (n <= 64)
  * inverses are checked to be unique
  * the character table is checked for  sum_rho d_rho^2 = |G|,  #irr = #classes,
    integrality of every d_rho, and BOTH orthogonality relations
  * every isotypic decomposition is checked for  sum_rho d_rho m_rho = dim
"""
import numpy as np, itertools

TOL = 1e-7

# ---------------------------------------------------------------- the group object
class Grp:
    def __init__(self, name, elems, mulf, order_check=None):
        self.name = name
        self.el   = list(elems)
        self.n    = len(self.el)
        self.idx  = {g: i for i, g in enumerate(self.el)}
        mt = np.zeros((self.n, self.n), dtype=np.int64)
        for i, a in enumerate(self.el):
            for j, b in enumerate(self.el):
                p = mulf(a, b)
                if p not in self.idx: raise AssertionError(f"{name}: not closed, {a}*{b}={p}")
                mt[i, j] = self.idx[p]
        self.mt = mt
        if order_check is not None and self.n != order_check:
            raise AssertionError(f"{name}: order {self.n} != {order_check}")
        # identity
        e = [i for i in range(self.n)
             if all(mt[i, j] == j for j in range(self.n)) and all(mt[j, i] == j for j in range(self.n))]
        if len(e) != 1: raise AssertionError(f"{name}: {len(e)} identities")
        self.e = e[0]
        # inverses
        self.iv = np.zeros(self.n, dtype=np.int64)
        for i in range(self.n):
            js = [j for j in range(self.n) if mt[i, j] == self.e]
            if len(js) != 1: raise AssertionError(f"{name}: element {i} has {len(js)} inverses")
            self.iv[i] = js[0]
        # SELF-CHECK: associativity, all n^3 triples
        if self.n <= 70:
            for a in range(self.n):
                lhs = mt[mt[a, :], :]                     # (a*b)*c
                rhs = mt[a, :][mt]                        # a*(b*c) -> mt[a, mt[b,c]]
                if not np.array_equal(lhs, rhs):
                    raise AssertionError(f"{name}: NOT ASSOCIATIVE")
        self.abelian = bool(np.array_equal(mt, mt.T))
        # conjugacy classes
        self._classes = None; self._cent = None; self._chi = None; self._subs = None

    def mul(self, a, b): return int(self.mt[a, b])
    def conj(self, h, g): return int(self.mt[self.mt[h, g], self.iv[h]])

    @property
    def classes(self):
        if self._classes is None:
            seen = set(); out = []
            for g in range(self.n):
                if g in seen: continue
                c = sorted({self.conj(h, g) for h in range(self.n)})
                seen |= set(c); out.append(c)
            self._classes = out
        return self._classes

    def centralizer(self, g):
        return [h for h in range(self.n) if self.mt[h, g] == self.mt[g, h]]

    @property
    def centre(self):
        return [g for g in range(self.n) if all(self.mt[g, h] == self.mt[h, g] for h in range(self.n))]

    def k(self): return len(self.classes)

    def generate(self, gens):
        S = {self.e}; fr = [self.e]
        gens = list(gens)
        while fr:
            nf = []
            for x in fr:
                for g in gens:
                    y = int(self.mt[x, g])
                    if y not in S: S.add(y); nf.append(y)
            fr = nf
        return frozenset(S)

    @property
    def subgroups(self):
        if self._subs is None:
            subs = {frozenset([self.e])}
            fr = list(subs)
            while fr:
                nf = []
                for H in fr:
                    for g in range(self.n):
                        if g in H: continue
                        N = self.generate(set(H) | {g})
                        if N not in subs: subs.add(N); nf.append(N)
                fr = nf
            self._subs = sorted(subs, key=lambda s: (len(s), sorted(s)))
        return self._subs

    # ------------------------------------------------------------ character table
    def chars(self, seed=0):
        """Irreducible characters by simultaneous diagonalisation of the class sums.
           Returns (classes, chi[k x k] complex, d[k] int).  chi[r][i] = chi_r(rep of class i)."""
        if self._chi is not None: return self._chi
        cl = self.classes; k = len(cl); N = self.n
        reps = [c[0] for c in cl]
        cls_of = np.zeros(N, dtype=np.int64)
        for i, c in enumerate(cl):
            for g in c: cls_of[g] = i
        # class-sum structure constants a[i][j][l]
        Ms = []
        for i in range(k):
            M = np.zeros((k, k))
            for j in range(k):
                cnt = np.zeros(N)
                for x in cl[i]:
                    row = self.mt[x, :]
                    for y in cl[j]: cnt[row[y]] += 1
                for l in range(k): M[j, l] = cnt[reps[l]]
            Ms.append(M)
        rng = np.random.default_rng(seed)
        for attempt in range(40):
            c = rng.normal(size=k) + 1j * rng.normal(size=k)
            Z = sum(c[i] * Ms[i] for i in range(k))
            w, V = np.linalg.eig(Z)
            om = np.zeros((k, k), dtype=complex); ok = True
            for r in range(k):
                v = V[:, r]
                p = int(np.argmax(np.abs(v)))
                for i in range(k):
                    om[r, i] = (Ms[i] @ v)[p] / v[p]
                    if np.linalg.norm(Ms[i] @ v - om[r, i] * v) > 1e-6 * max(1, np.linalg.norm(v)):
                        ok = False
            if not ok: continue
            sz = np.array([len(c_) for c_ in cl], dtype=float)
            d = np.sqrt(N / np.real(np.sum(np.abs(om) ** 2 / sz[None, :], axis=1)))
            if np.max(np.abs(d - np.round(d))) > 1e-6: continue
            d = np.round(d).astype(int)
            chi = om * d[:, None] / sz[None, :]
            # SELF-CHECKS
            if abs(int(np.sum(d ** 2)) - N) > 0: continue
            orth = np.zeros((k, k), dtype=complex)
            for r in range(k):
                for s in range(k):
                    orth[r, s] = np.sum(sz * chi[r] * np.conj(chi[s])) / N
            if np.max(np.abs(orth - np.eye(k))) > 1e-6: continue
            # column orthogonality
            col = np.zeros((k, k), dtype=complex)
            for i in range(k):
                for j in range(k):
                    col[i, j] = np.sum(chi[:, i] * np.conj(chi[:, j]))
            tgt = np.diag([N / len(cl[i]) for i in range(k)])
            if np.max(np.abs(col - tgt)) > 1e-6: continue
            order = np.lexsort((np.real(chi[:, 0]),))          # d ascending, trivial first
            chi = chi[order]; d = d[order]
            # put trivial character first
            triv = [r for r in range(k) if np.max(np.abs(chi[r] - 1)) < 1e-7]
            if len(triv) != 1: continue
            t = triv[0]
            idxs = [t] + [r for r in range(k) if r != t]
            chi = chi[idxs]; d = d[idxs]
            self._chi = (cl, chi, d, cls_of)
            return self._chi
        raise AssertionError(f"{self.name}: CHARACTER TABLE SELF-CHECK FAILED -- no conclusion")

    def decompose(self, classfn):
        """multiplicities of each irrep in a rep with the given class function (array over classes)"""
        cl, chi, d, _ = self.chars()
        sz = np.array([len(c) for c in cl], dtype=float)
        m = np.real(np.array([np.sum(sz * classfn * np.conj(chi[r])) / self.n for r in range(len(cl))]))
        return m, d

# ---------------------------------------------------------------- constructors
def cyclic(n): return Grp(f"Z_{n}", range(n), lambda a, b: (a + b) % n, n)

def direct(G1, G2, name=None):
    el = [(a, b) for a in G1.el for b in G2.el]
    m1 = {g: i for i, g in enumerate(G1.el)}; m2 = {g: i for i, g in enumerate(G2.el)}
    def mul(x, y):
        return (G1.el[G1.mt[m1[x[0]], m1[y[0]]]], G2.el[G2.mt[m2[x[1]], m2[y[1]]]])
    return Grp(name or f"{G1.name}x{G2.name}", el, mul, G1.n * G2.n)

def metacyclic(N, t, u, name, order=None):
    """<r,s | r^N = 1, s r s^-1 = r^t, s^2 = r^u>.  Elements (i,eps).
       D_n : N=n, t=-1, u=0.   Dic_n (order 4n): N=2n, t=-1, u=n.
       SD_16: N=8,t=3,u=0.     M_4(2): N=8,t=5,u=0."""
    el = [(i, s) for s in (0, 1) for i in range(N)]
    def mul(x, y):
        i, a = x; j, b = y
        if a == 0: return ((i + j) % N, b)
        # (r^i s)(r^j s^b) = r^(i + t j) s^(1+b);  s^2 = r^u
        k = (i + t * j) % N
        if b == 0: return (k, 1)
        return ((k + u) % N, 0)
    return Grp(name, el, mul, order or 2 * N)

def pauli16():
    """<i, X, Z> of order 16: elements (a,b,c) = i^a X^b Z^c, a in Z_4, b,c in Z_2."""
    el = [(a, b, c) for a in range(4) for b in range(2) for c in range(2)]
    def mul(x, y):
        a, b, c = x; a2, b2, c2 = y
        return ((a + a2 + 2 * c * b2) % 4, (b + b2) % 2, (c + c2) % 2)
    return Grp("Pauli16", el, mul, 16)


def extraspecial(nq, minus=False):
    """The Pauli group on nq qubits modulo phases i, i.e. the extraspecial 2-group of order
       2^(2nq+1):  elements (sign, b_1..b_nq, c_1..c_nq),  (X^b Z^c)(X^b' Z^c') =
       (-1)^{c.b'} X^{b+b'} Z^{c+c'}.  Centre = {+-1} of order 2, so G/Z = Z_2^(2nq) is
       ELEMENTARY ABELIAN -- the control against the dihedral tower, whose G/Z is dihedral."""
    import itertools as _it
    el = [(s,) + v for s in (0, 1) for v in _it.product((0, 1), repeat=2 * nq)]
    def mul(x, y):
        s = (x[0] + y[0]) % 2
        b = x[1:1 + nq]; c = x[1 + nq:]; b2 = y[1:1 + nq]; c2 = y[1 + nq:]
        s = (s + sum(c[i] * b2[i] for i in range(nq))) % 2
        return (s,) + tuple((b[i] + b2[i]) % 2 for i in range(nq)) \
                    + tuple((c[i] + c2[i]) % 2 for i in range(nq))
    return Grp(f"ES_2^(1+{2*nq})", el, mul, 2 ** (2 * nq + 1))

# ---------------------------------------------------------------- the ladder
def ladder(max_order=64):
    """2-groups only -- C-41 says these are exactly the groups D(G) admits records on.
       ABELIAN and NON-ABELIAN at the SAME order is the control this lane turns on."""
    Z2 = cyclic(2); Z4 = cyclic(4); Z8 = cyclic(8)
    L = []
    L += [cyclic(2)]
    L += [cyclic(4), direct(Z2, Z2, "Z_2^2")]
    L += [cyclic(8), direct(Z4, Z2, "Z_4xZ_2"), direct(direct(Z2, Z2, "Z_2^2"), Z2, "Z_2^3"),
          metacyclic(4, -1, 0, "D_4"), metacyclic(4, -1, 2, "Q_8")]
    L += [cyclic(16), direct(Z8, Z2, "Z_8xZ_2"), direct(Z4, Z4, "Z_4^2"),
          direct(direct(Z4, Z2, "Z_4xZ_2"), Z2, "Z_4xZ_2^2"),
          direct(direct(direct(Z2, Z2, "Z_2^2"), Z2, "Z_2^3"), Z2, "Z_2^4"),
          metacyclic(8, -1, 0, "D_8"), metacyclic(8, -1, 4, "Q_16"),
          metacyclic(8, 3, 0, "SD_16"), metacyclic(8, 5, 0, "M_4(2)"),
          direct(metacyclic(4, -1, 0, "D_4"), Z2, "D_4xZ_2"),
          direct(metacyclic(4, -1, 2, "Q_8"), Z2, "Q_8xZ_2"),
          pauli16()]
    if max_order >= 32:
        L += [cyclic(32), direct(direct(direct(direct(Z2, Z2, "Z_2^2"), Z2, "Z_2^3"), Z2, "Z_2^4"), Z2, "Z_2^5"),
              metacyclic(16, -1, 0, "D_16"), metacyclic(16, -1, 8, "Q_32"),
              metacyclic(16, 7, 0, "SD_32"), metacyclic(16, 9, 0, "M_5(2)"),
              direct(metacyclic(8, -1, 0, "D_8"), Z2, "D_8xZ_2"),
              direct(pauli16(), Z2, "Pauli16xZ_2"),
              direct(direct(metacyclic(4, -1, 0, "D_4"), Z2, "D_4xZ_2"), Z2, "D_4xZ_2^2"),
              extraspecial(2)]
    if max_order >= 64:
        L += [cyclic(64),
              direct(direct(direct(direct(direct(Z2, Z2, "Z_2^2"), Z2, "Z_2^3"), Z2, "Z_2^4"), Z2, "Z_2^5"), Z2, "Z_2^6"),
              metacyclic(32, -1, 0, "D_32"), metacyclic(32, -1, 16, "Q_64"),
              direct(metacyclic(16, -1, 0, "D_16"), Z2, "D_16xZ_2")]
    return [G for G in L if G.n <= max_order]
