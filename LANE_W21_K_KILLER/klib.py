"""
klib.py -- LANE W21-K (THE KILLER).  Instrument only; no claims live here.

Everything is built from scratch: Pauli strings as symplectic vectors over F_2 with a
Hermitian representative P(x,z) = i^(x.z) X^x Z^z, Z_2 lattice gauge theory on a graph,
the Gauss projector, exact diagonalisation, and -- the load-bearing piece -- the entropy
of a state RESTRICTED TO AN ARBITRARY PAULI-SPANNED SUBALGEBRA, with its centre and its
type-I block decomposition computed explicitly rather than assumed.

No import from any other lane.  Cross-checked against brute force in k1_typecheck.py.
"""
import numpy as np
import itertools

# ----------------------------------------------------------------------------------
# Pauli strings as symplectic vectors over F_2.  v = (x | z) in F_2^{2n}.
# Hermitian representative  P(x,z) = i^{x.z} X^x Z^z   (so P=P^dag, P^2 = I).
# ----------------------------------------------------------------------------------

_PC = {}


def _popcnt(arr):
    return np.bitwise_count(arr) if hasattr(np, "bitwise_count") else np.vectorize(
        lambda c: bin(int(c)).count("1"))(arr)


def _signs(z, n):
    key = (z, n)
    s = _PC.get(key)
    if s is None:
        cols = np.arange(1 << n, dtype=np.int64)
        s = (1.0 - 2.0 * (_popcnt(cols & int(z)) & 1)).astype(np.complex128)
        _PC[key] = s
    return s


_MATCACHE = {}


def pauli_matrix(x, z, n):
    """Dense d x d Hermitian Pauli.  One nonzero per column: row = col ^ x."""
    key = (int(x), int(z), n)
    M = _MATCACHE.get(key)
    if M is not None:
        return M
    d = 1 << n
    xi = int(x); zi = int(z)
    ph = (1j) ** (bin(xi & zi).count("1"))
    cols = np.arange(d, dtype=np.int64)
    M = np.zeros((d, d), dtype=complex)
    M[cols ^ xi, cols] = ph * _signs(zi, n)
    if len(_MATCACHE) < 4096:
        _MATCACHE[key] = M
    return M


def pauli_expect(psi, x, z, n):
    """<psi| P(x,z) |psi>, O(d), no matrix built."""
    d = 1 << n
    xi = int(x); zi = int(z)
    ph = (1j) ** (bin(xi & zi).count("1"))
    cols = np.arange(d, dtype=np.int64)
    return complex(np.vdot(psi[cols ^ xi], ph * _signs(zi, n) * psi))


def sympl(v1, v2):
    """Symplectic form over F_2 on packed (x,z) integer pairs. 0 = commute, 1 = anticommute."""
    x1, z1 = v1
    x2, z2 = v2
    return (bin(x1 & z2).count("1") + bin(z1 & x2).count("1")) & 1


def vadd(v1, v2):
    return (v1[0] ^ v2[0], v1[1] ^ v2[1])


# ----------------------------------------------------------------------------------
# F_2 linear algebra on packed 2n-bit vectors (encode v=(x,z) as x | z<<n).
# ----------------------------------------------------------------------------------

def pack(v, n):
    return v[0] | (v[1] << n)


def unpack(p, n):
    m = (1 << n) - 1
    return (p & m, (p >> n) & m)


def rref_basis(vecs):
    """Row-reduce a list of ints over F_2; return an independent basis (list of ints)."""
    basis = []
    for v in vecs:
        cur = v
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    # re-reduce to a canonical echelon basis
    out = []
    for v in sorted(basis, reverse=True):
        cur = v
        for b in out:
            cur = min(cur, cur ^ b)
        if cur:
            out.append(cur)
            out.sort(reverse=True)
    return out


def span_elements(basis):
    """All 2^k elements of the span (as ints)."""
    out = [0]
    for b in basis:
        out += [o ^ b for o in out]
    return out


def in_span(v, basis):
    cur = v
    for b in basis:
        cur = min(cur, cur ^ b)
    return cur == 0


def intersect(b1, b2, n):
    """Basis of span(b1) cap span(b2) over F_2 (Zassenhaus, done crudely but exactly)."""
    b1 = rref_basis(b1); b2 = rref_basis(b2)
    # solve for combinations of b1 that lie in span(b2)
    piv = {}
    for v in b2:
        cur = v
        for p, rr in piv.items():
            if (cur >> p) & 1:
                cur ^= rr
        if cur:
            piv[cur.bit_length() - 1] = cur
    def red(v):
        cur = v
        for p in sorted(piv.keys(), reverse=True):
            if (cur >> p) & 1:
                cur ^= piv[p]
        return cur
    # gaussian elimination on (red(b1_i) | e_i)
    rows = [(red(v), 1 << i) for i, v in enumerate(b1)]
    out = []
    piv2 = {}
    for r, tag in rows:
        cur, ct = r, tag
        while cur:
            p = cur.bit_length() - 1
            if p in piv2:
                cur ^= piv2[p][0]; ct ^= piv2[p][1]
            else:
                piv2[p] = (cur, ct); cur = 0; ct = 0; break
        if cur == 0 and ct != 0:
            v = 0
            for i in range(len(b1)):
                if (ct >> i) & 1:
                    v ^= b1[i]
            if v:
                out.append(v)
    return rref_basis(out)


def sympl_perp(basis, n):
    """Symplectic orthogonal complement of span(basis) inside F_2^{2n}, as a basis."""
    # solve: for each generator g, sympl(v,g)=0.  Linear system over F_2 in 2n unknowns.
    rows = []
    for g in basis:
        gx, gz = unpack(g, n)
        # sympl((x,z),(gx,gz)) = |x & gz| + |z & gx| mod 2 -> row over (x,z) coords
        r = gz | (gx << n)
        rows.append(r)
    rows = rref_basis(rows)
    # nullspace of the linear map v -> (<r_i, v>) ; build by Gaussian elimination
    N = 2 * n
    piv = {}
    for r in rows:
        cur = r
        for p, rr in piv.items():
            if (cur >> p) & 1:
                cur ^= rr
        if cur:
            p = cur.bit_length() - 1
            piv[p] = cur
    free = [i for i in range(N) if i not in piv]
    nullb = []
    for f in free:
        v = 1 << f
        # back-substitute
        for p in sorted(piv.keys()):
            rr = piv[p]
            # inner product of v with rr must be 0; adjust bit p
            ip = bin(v & rr).count("1") & 1
            if ip:
                v ^= (1 << p)
        nullb.append(v)
    # verify
    for v in nullb:
        for r in rows:
            assert bin(v & r).count("1") % 2 == 0, "sympl_perp failed"
    return rref_basis(nullb)


# ----------------------------------------------------------------------------------
# The algebra of a Pauli subspace: centre, block structure, restricted entropy.
# ----------------------------------------------------------------------------------

class PauliAlgebra:
    """span{ P(v) : v in span(basis) } as a *-subalgebra of B(C^{2^n})."""

    def __init__(self, basis, n, name=""):
        self.basis = rref_basis(list(basis))
        self.n = n
        self.name = name
        self.dimV = len(self.basis)
        self.dim = 1 << self.dimV                      # complex dim of the algebra
        # radical of the symplectic form restricted to V = centre generators
        elems = span_elements(self.basis)
        rad = []
        for e in self.basis:
            pass
        # centre = { v in V : sympl(v,w)=0 for all w in V }
        cand = []
        for e in elems:
            ve = unpack(e, n)
            if all(sympl(ve, unpack(w, n)) == 0 for w in self.basis):
                cand.append(e)
        self.centre_basis = rref_basis(cand)
        self.r = len(self.centre_basis)                # dim of centre subspace
        assert (self.dimV - self.r) % 2 == 0
        self.k = (self.dimV - self.r) // 2             # blocks are M_{2^k}
        self.nblocks = 1 << self.r
        self.mult = 1 << (n - self.r - self.k)         # multiplicity per block
        assert self.mult >= 1, (n, self.r, self.k)

    def describe(self):
        return (f"{self.name}: dim_F2(V)={self.dimV}  dim_C(A)={self.dim}  "
                f"centre dim={self.nblocks} (r={self.r})  blocks={self.nblocks} x M_{2**self.k}  "
                f"multiplicity m={self.mult}")

    def omega(self, psi, exp_cache=None):
        """Hilbert-Schmidt projection of |psi><psi| onto A.  omega in A, and
        Tr(omega a) = <psi|a|psi> for all a in A."""
        n = self.n
        d = 1 << n
        om = np.zeros((d, d), dtype=complex)
        for e in span_elements(self.basis):
            v = unpack(e, n)
            if exp_cache is not None:
                c = exp_cache.get(e)
                if c is None:
                    c = pauli_expect(psi, v[0], v[1], n); exp_cache[e] = c
            else:
                c = pauli_expect(psi, v[0], v[1], n)
            if abs(c) > 1e-13:
                om += c * pauli_matrix(v[0], v[1], n)
        return om / d

    def entropy(self, psi, exp_cache=None):
        """S(rho|_A) in BITS.  Derivation in NOTES; validated against brute force."""
        om = self.omega(psi, exp_cache)
        w = np.linalg.eigvalsh((om + om.conj().T) / 2)
        w = np.clip(w.real, 0.0, None)
        m = self.mult
        s = 0.0
        for nu in w:
            if nu > 1e-14:
                s -= nu * np.log2(nu * m)
        return s

    def central_distribution(self, psi):
        """p_chi over the 2^r central sectors, and the sector labels."""
        n = self.n
        cb = self.centre_basis
        out = []
        for signs in itertools.product([1, -1], repeat=len(cb)):
            # p = <prod (I + s_i P_i)/2>  -- expand
            p = 0.0
            for sub in range(1 << len(cb)):
                coef = 1.0
                acc = 0
                for i in range(len(cb)):
                    if (sub >> i) & 1:
                        coef *= signs[i]
                        acc ^= cb[i]
                v = unpack(acc, n)
                # NOTE: products of commuting Hermitian Paulis may pick a sign; we take
                # the Hermitian representative of the product and correct its sign.
                sgn = _prod_sign(cb, sub, n)
                p += coef * sgn * pauli_expect(psi, v[0], v[1], n).real
            out.append((signs, p / (1 << len(cb))))
        return out


def _prod_sign(basis, sub, n):
    """Sign s such that prod_{i in sub} P(b_i) = s * P(sum b_i), for mutually
    commuting Hermitian Paulis (so s = +-1).  Computed numerically-free via cocycle."""
    # P(a)P(b) = eps * P(a+b).  Track by explicit small-matrix check on 1 qubit? No --
    # use the standard cocycle: for Hermitian reps P=i^{x.z}X^xZ^z,
    # P(a)P(b) = i^{ax.az + bx.bz - cx.cz} (-1)^{az.bx} P(c), c=a+b.
    cur = 0
    ph = 1 + 0j
    for i in range(len(basis)):
        if not ((sub >> i) & 1):
            continue
        b = basis[i]
        ax, az = unpack(cur, n)
        bx, bz = unpack(b, n)
        cx, cz = ax ^ bx, az ^ bz
        e = (bin(ax & az).count("1") + bin(bx & bz).count("1") - bin(cx & cz).count("1"))
        ph *= (1j) ** e * ((-1) ** (bin(az & bx).count("1")))
        cur = cur ^ b
    assert abs(ph.imag) < 1e-9, "non-commuting centre generators"
    return 1 if ph.real > 0 else -1


def mutual_information(algA, algB, algAB, psi):
    """I = S(A) + S(B) - S(AB) with AB the algebra generated by A and B (join)."""
    return algA.entropy(psi) + algB.entropy(psi) - algAB.entropy(psi)


def join(algA, algB, name=""):
    return PauliAlgebra(algA.basis + algB.basis, algA.n, name)


# ----------------------------------------------------------------------------------
# Z_2 lattice gauge theory on a graph.
# ----------------------------------------------------------------------------------

class Z2Gauge:
    def __init__(self, nverts, edges, name=""):
        self.V = nverts
        self.edges = list(edges)
        self.L = len(self.edges)
        self.n = self.L
        self.name = name
        # Gauss generators G_v = prod_{l ~ v} X_l
        self.gauss = []
        for v in range(nverts):
            x = 0
            for i, (a, b) in enumerate(self.edges):
                if a == v or b == v:
                    if a == v and b == v:
                        continue
                    x |= (1 << i)
            self.gauss.append((x, 0))
        # cycle space basis (independent Wilson loops), from a spanning tree
        self.cycles = self._cycle_basis()
        self.cycle_dim = len(self.cycles)
        assert self.cycle_dim == self.L - self.V + self._ncomp()

    def _ncomp(self):
        par = list(range(self.V))
        def f(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        for a, b in self.edges:
            ra, rb = f(a), f(b)
            if ra != rb:
                par[ra] = rb
        return len({f(v) for v in range(self.V)})

    def _cycle_basis(self):
        par = list(range(self.V))
        def f(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        tree, extra = [], []
        for i, (a, b) in enumerate(self.edges):
            ra, rb = f(a), f(b)
            if ra != rb and a != b:
                par[ra] = rb
                tree.append(i)
            else:
                extra.append(i)
        # adjacency of the tree
        adj = {v: [] for v in range(self.V)}
        for i in tree:
            a, b = self.edges[i]
            adj[a].append((b, i)); adj[b].append((a, i))
        cycles = []
        for i in extra:
            a, b = self.edges[i]
            if a == b:
                cycles.append(1 << i); continue
            # BFS path a->b in tree
            prev = {a: (None, None)}
            q = [a]
            while q:
                u = q.pop(0)
                if u == b:
                    break
                for (w, li) in adj[u]:
                    if w not in prev:
                        prev[w] = (u, li); q.append(w)
            mask = 1 << i
            u = b
            while prev[u][0] is not None:
                mask |= (1 << prev[u][1]); u = prev[u][0]
            cycles.append(mask)
        return cycles

    # --- operators as symplectic vectors ---
    def X(self, l):
        return (1 << l, 0)

    def Zop(self, mask):
        return (0, mask)

    def W(self, cyc_mask):
        return (0, cyc_mask)

    # --- physical projector & ground state ---
    def phys_projector(self):
        d = 1 << self.n
        P = np.eye(d, dtype=complex)
        for (x, z) in self.gauss:
            P = P @ (np.eye(d, dtype=complex) + pauli_matrix(x, z, self.n)) / 2
        return P

    def hamiltonian(self, g2):
        d = 1 << self.n
        H = np.zeros((d, d), dtype=complex)
        for c in self.cycles:
            H -= (1.0 / g2) * pauli_matrix(0, c, self.n)
        for l in range(self.L):
            H -= g2 * pauli_matrix(1 << l, 0, self.n)
        return H

    def ground_state(self, g2):
        """Ground state inside the Gauss-physical sector, returned in the ambient space."""
        P = self.phys_projector()
        w, V = np.linalg.eigh((P + P.conj().T) / 2)
        cols = V[:, w > 0.5]                       # orthonormal basis of the physical sector
        H = self.hamiltonian(g2)
        Hp = cols.conj().T @ H @ cols
        wv, Vv = np.linalg.eigh((Hp + Hp.conj().T) / 2)
        psi = cols @ Vv[:, 0]
        psi = psi / np.linalg.norm(psi)
        return psi, cols.shape[1], float(wv[0])

    def gauss_residual(self, psi):
        r = 0.0
        for (x, z) in self.gauss:
            r = max(r, abs(1.0 - pauli_expect(psi, x, z, self.n).real))
        return r

    # --- region bookkeeping ---
    def region(self, verts):
        vs = set(verts)
        internal, boundary, external = [], [], []
        for i, (a, b) in enumerate(self.edges):
            ina, inb = a in vs, b in vs
            if ina and inb:
                internal.append(i)
            elif ina or inb:
                boundary.append(i)
            else:
                external.append(i)
        return internal, boundary, external


def gauge_invariant_subspace(gauge):
    """Basis of the F_2 subspace of Pauli strings commuting with every Gauss generator."""
    n = gauge.n
    gb = rref_basis([pack(g, n) for g in gauge.gauss])
    return sympl_perp(gb, n)
