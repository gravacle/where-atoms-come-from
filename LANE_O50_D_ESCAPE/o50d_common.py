"""O-50 D  ESCAPE LANE -- common machinery.

CANONICAL CARRIER (O-49, D-23): the TORIC CODE on an L x L torus.  Qubits on edges,
n = 2L^2, H = -sum_v A_v - sum_p B_p.  Clause (v) is realised by manifold homology, so
nothing here rests on the 1D proxy convention.

Everything below is EXACT F_2 linear algebra unless a routine says otherwise.  Dense
Hilbert-space objects are built only at L = 2 (dim 256).
"""
import sys, itertools
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv, build_writer)

def say(*a):
    print(*a); sys.stdout.flush()

# ------------------------------------------------------------------ F_2 linear algebra
def rref(rows, ncols):
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    return M[:r], piv

def rank2(rows, ncols):
    return len(rref(rows, ncols)[0])

def in_span(v, rows, ncols):
    R, piv = rref(rows, ncols)
    v = list(v)
    for i, c in enumerate(piv):
        if v[c]: v = [(x + y) % 2 for x, y in zip(v, R[i])]
    return not any(v)

def nullspace2(M, ncols):
    """basis of {v : M v = 0} over F_2, M given as list of rows of length ncols"""
    R, piv = rref(M, ncols)
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        v = [0] * ncols; v[f] = 1
        for i, c in enumerate(piv): v[c] = R[i][f]
        out.append(v)
    return out

# ------------------------------------------------------------------ the torus
class Torus:
    """L x L torus.  edge h(i,j) = i*L+j  (from vertex (i,j) to (i,j+1))
                     edge v(i,j) = L*L + i*L + j  (from vertex (i,j) to (i+1,j))"""
    def __init__(self, L):
        self.L = L; self.nq = 2 * L * L
        self.h = lambda i, j: (i % L) * L + (j % L)
        self.v = lambda i, j: L * L + (i % L) * L + (j % L)
        # X-type vertex stars and Z-type plaquettes, as edge sets
        self.star = [[self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)]
                     for i in range(L) for j in range(L)]
        self.plaq = [[self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)]
                     for i in range(L) for j in range(L)]
        self.stab = ([self.xz(s, []) for s in self.star] + [self.xz([], p) for p in self.plaq])
        # edge midpoint coordinates, for geometric distance on the torus
        self.pos = {}
        for i in range(L):
            for j in range(L):
                self.pos[self.h(i, j)] = (i, j + 0.5)
                self.pos[self.v(i, j)] = (i + 0.5, j)

    def xz(self, xs, zs):
        v = [0] * (2 * self.nq)
        for e in xs: v[e] ^= 1
        for e in zs: v[self.nq + e] ^= 1
        return v

    def sp(self, a, b):
        n = self.nq
        return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

    def weight(self, a):
        n = self.nq
        return sum(1 for i in range(n) if a[i] or a[n + i])

    def support(self, a):
        n = self.nq
        return [i for i in range(n) if a[i] or a[n + i]]

    def dist(self, e, f):
        L = self.L
        (a, b), (c, d) = self.pos[e], self.pos[f]
        dx = abs(a - c); dy = abs(b - d)
        return min(dx, L - dx) + min(dy, L - dy)

# ------------------------------------------------------------------ dense operators (L=2 only)
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Ym = 1j * Xm @ Zm

def dense(vec, nq):
    """(x|z) F_2 vector -> dense Pauli operator (Hermitian phase convention of xz_to_matrix)"""
    return xz_to_matrix(vec, nq)

def dense_from_sets(xs, zs, nq):
    M = np.array([[1]], dtype=complex)
    for q in range(nq):
        P = I2
        if q in xs and q in zs: P = Xm @ Zm
        elif q in xs: P = Xm
        elif q in zs: P = Zm
        M = np.kron(M, P)
    return M

# ------------------------------------------------------------------ record-blind control (C-61)
def control_carrier(n, seed=0):
    """C-61's control: H' = sum J_i Z_i Z_{i+1} + sum h_i Z_i with distinct h_i.
       Fully NON-DEGENERATE spectrum => by P-1 it holds ZERO records.
       Returned dense, plus its 2^n basis-state 'configuration' labels."""
    rng = np.random.default_rng(seed)
    J = rng.normal(size=n - 1) + 1.5
    h = np.array([0.31 * (k + 1) + 0.017 * k * k for k in range(n)])
    dim = 2 ** n
    diag = np.zeros(dim)
    for b in range(dim):
        s = [1 - 2 * ((b >> (n - 1 - q)) & 1) for q in range(n)]
        diag[b] = sum(J[q] * s[q] * s[q + 1] for q in range(n - 1)) + sum(h[q] * s[q] for q in range(n))
    return np.diag(diag).astype(complex), J, h, diag
