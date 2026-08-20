"""Exact F_2 linear algebra + toric-code construction on an L x L torus.
   No dense Hilbert-space matrices anywhere in this file.
   Convention: a Pauli is a vector (x | z) in F_2^{2n}; sp(a,b) = sum_i a_i b_{n+i} + a_{n+i} b_i.
   This matches record_model.symplectic_logicals and record_model.xz_to_matrix exactly."""

def sp(a, b, n):
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def rref(rows, width):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(width):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

def rank(rows, width):
    return len(rref(rows, width)[1])

def in_span(v, basis, width):
    """basis must be in rref form (as returned by rref)."""
    B, piv = rref(basis, width)
    v = v[:]
    for row, c in zip(B, piv):
        if v[c]: v = [(x + y) % 2 for x, y in zip(v, row)]
    return not any(v)

def nullspace(rows, width):
    """basis of {v : M v = 0} for M with the given rows, over F_2."""
    R, piv = rref(rows, width)
    free = [c for c in range(width) if c not in piv]
    out = []
    for f in free:
        v = [0] * width; v[f] = 1
        for i, c in enumerate(piv): v[c] = R[i][f]
        out.append(v)
    return out

def span(basis, width):
    """every vector in the F_2 span (small bases only)."""
    out = [[0] * width]
    for b in basis:
        out += [[(x + y) % 2 for x, y in zip(v, b)] for v in out]
    return out

# ------------------------------------------------------------------ toric code
class Toric:
    """qubits on edges of an L x L periodic square lattice.
       h(i,j) = edge from vertex (i,j) to (i,j+1)   index  i*L + j
       v(i,j) = edge from vertex (i,j) to (i+1,j)   index  L*L + i*L + j
       A_v(i,j) = X on {h(i,j), h(i,j-1), v(i,j), v(i-1,j)}
       B_p(i,j) = Z on {h(i,j), h(i+1,j), v(i,j), v(i,j+1)}"""
    def __init__(self, L):
        self.L = L; self.n = 2 * L * L
        self.h = lambda i, j: (i % L) * L + (j % L)
        self.v = lambda i, j: L * L + (i % L) * L + (j % L)
        self.A, self.B = [], []
        for i in range(L):
            for j in range(L):
                s = [0] * (2 * self.n)
                for e in (self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)): s[e] = 1
                self.A.append(s)                                   # X-type: x-block
                t = [0] * (2 * self.n)
                for e in (self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)):
                    t[self.n + e] = 1                              # Z-type: z-block
                self.B.append(t)
        self.stab = self.A + self.B

    # --- vertex boundary map (for Z-type: cycles) and plaquette map (for X-type: cocycles)
    def vertex_rows(self):
        """rows over the n edge coordinates: row_v[e] = 1 iff edge e touches vertex v."""
        L = self.L; rows = []
        for i in range(L):
            for j in range(L):
                r = [0] * self.n
                for e in (self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)): r[e] ^= 1
                rows.append(r)
        return rows

    def plaquette_rows(self):
        L = self.L; rows = []
        for i in range(L):
            for j in range(L):
                r = [0] * self.n
                for e in (self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)): r[e] ^= 1
                rows.append(r)
        return rows

    # --- homology invariants (exact, convention-free: crossing parities with the two dual cuts)
    def z_class(self, zs, j0=0, i0=0):
        """zs: length-n support of a Z-type cycle.  (w_x, w_y) in H_1(T^2, F_2)."""
        L = self.L
        wx = sum(zs[self.h(i, j0)] for i in range(L)) % 2
        wy = sum(zs[self.v(i0, j)] for j in range(L)) % 2
        return (wx, wy)

    def x_class(self, xs, j0=0, i0=0):
        """xs: length-n support of an X-type cocycle (dual cycle)."""
        L = self.L
        w1 = sum(xs[self.v(i, j0)] for i in range(L)) % 2
        w2 = sum(xs[self.h(i0, j)] for j in range(L)) % 2
        return (w1, w2)
