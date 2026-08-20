"""O-48-D SHARED MACHINERY.  Nothing here nominates anything; every writer is SEARCHED.

TWO REPRESENTATIONS, KEPT IN LOCKSTEP SO EACH CHECKS THE OTHER:
  (1) EXACT SYMBOLIC.  A Pauli word is (a|b) in F_2^{2n}.  Two Paulis anticommute iff the
      symplectic form a.b' + b.a' = 1.  For a Pauli R and Pauli W:  W-dag R W = -R  IFF  they
      ANTICOMMUTE.  For H = sum_k J_k P_k with the P_k DISTINCT Pauli words and J_k != 0, the
      P_k are linearly independent, so  W H W-dag = sum_k eps_k J_k P_k  equals H IFF eps_k = +1
      for every k.  Hence  [W,H] = 0  IFF  W commutes with EVERY TERM.  That is EXACT -- it is
      not a tolerance and not a fit.
  (2) DENSE NUMERIC.  Build the matrices and measure ||[W,H]||, ||W-dag R W + R||, Tr(P_E R).
      Used at every n where 2^n fits, as a CONTROL on (1).

DIAGONAL SHORTCUT.  When every term is Z-type the Hamiltonian is DIAGONAL in the configuration
basis, energies are exact INTEGER sums when the couplings are integers (D-19: no floats where the
conclusion turns on cancellation), and eigenspaces are exact integer-equality classes.  No dense
2^n x 2^n matrix is ever needed for clauses (i)-(iv).
"""
import numpy as np
from itertools import product

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Ym = 1j * Xm @ Zm


# ------------------------------------------------------------------ Pauli words as (a|b)
def pauli_matrix(a, b):
    """(a|b) -> the Hermitian Pauli word  i^{a.b} X^a Z^b  (real phase, squares to I)."""
    n = len(a)
    M = np.array([[1]], dtype=complex)
    for i in range(n):
        x, z = a[i], b[i]
        P = I2 if (x, z) == (0, 0) else (Xm if (x, z) == (1, 0) else (Zm if (x, z) == (0, 1) else Ym))
        M = np.kron(M, P)
    return M


def symp(a1, b1, a2, b2):
    """1 if the two Pauli words ANTICOMMUTE, 0 if they commute."""
    return (sum(a1[i] * b2[i] + b1[i] * a2[i] for i in range(len(a1)))) % 2


def weight(a, b):
    return sum(1 for i in range(len(a)) if (a[i], b[i]) != (0, 0))


def all_paulis(n):
    """Every Pauli word on n qubits, as (a,b) tuples.  4^n of them -- small n only."""
    for av in product((0, 1), repeat=n):
        for bv in product((0, 1), repeat=n):
            yield av, bv


# ------------------------------------------------------------------ F_2 linear algebra
def f2_rank(rows, ncols):
    rows = [r[:] for r in rows]
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


def f2_nullspace(rows, ncols):
    """Basis of {v : M v = 0} over F_2, M given by `rows` (each a length-ncols 0/1 list)."""
    R = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        p = next((i for i in range(r, len(R)) if R[i][c]), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        for i in range(len(R)):
            if i != r and R[i][c]:
                R[i] = [(x + y) % 2 for x, y in zip(R[i], R[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = [0] * ncols
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = R[i][f]
        basis.append(v)
    return basis


def f2_span(basis, ncols):
    """Every vector in the span (2^dim of them) -- small dims only."""
    out = []
    d = len(basis)
    for m in range(1 << d):
        v = [0] * ncols
        for i in range(d):
            if (m >> i) & 1:
                v = [(x + y) % 2 for x, y in zip(v, basis[i])]
        out.append(v)
    return out


def f2_in_span(v, basis, ncols):
    """Is v in the F_2 span of basis?"""
    return f2_rank(list(basis) + [list(v)], ncols) == f2_rank(list(basis), ncols)


# ------------------------------------------------------------------ commuting-Pauli Hamiltonian
class PauliH:
    """H = sum_k J_k P_k with P_k distinct Pauli words.  Terms need NOT commute; if they do
       not, we say so and refuse the diagonal shortcut."""

    def __init__(self, n, terms):
        self.n = n
        self.terms = [(list(a), list(b), J) for a, b, J in terms]
        self.commuting = all(symp(a1, b1, a2, b2) == 0
                             for i, (a1, b1, _) in enumerate(self.terms)
                             for (a2, b2, _) in self.terms[i + 1:])
        self.z_type = all(all(x == 0 for x in a) for a, b, J in self.terms)

    def matrix(self):
        n = self.n
        M = np.zeros((2 ** n, 2 ** n), dtype=complex)
        for a, b, J in self.terms:
            M = M + J * pauli_matrix(a, b)
        return M

    def admissible(self, a, b):
        """[W,H] = 0.  EXACT: W must commute with every term (terms are distinct Pauli words,
           hence linearly independent, so no cancellation between them is possible)."""
        return all(symp(a, b, ta, tb) == 0 for ta, tb, _ in self.terms)

    def commutation_rows(self):
        """Rows of the F_2 matrix whose nullspace is the admissible Pauli group, in the
           variable order (a_0..a_{n-1}, b_0..b_{n-1})."""
        n = self.n
        rows = []
        for ta, tb, _ in self.terms:
            rows.append([tb[i] for i in range(n)] + [ta[i] for i in range(n)])
        return rows

    def stabiliser_rows(self):
        """The terms themselves as F_2 vectors (a|b) -- the group S they generate."""
        return [list(a) + list(b) for a, b, _ in self.terms]


# ------------------------------------------------------------------ diagonal (Z-type) machinery
def spin_table(n):
    """All 2^n configurations as +-1 spins, shape (2^n, n).  int8."""
    idx = np.arange(1 << n, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int8)
    return (1 - 2 * bits).astype(np.int8)


def diag_energies(n, zterms):
    """zterms = list of (support_tuple, J_int).  Exact INTEGER energies, shape (2^n,)."""
    sig = spin_table(n)
    E = np.zeros(1 << n, dtype=np.int64)
    for sup, J in zterms:
        p = np.ones(1 << n, dtype=np.int64)
        for i in sup:
            p = p * sig[:, i].astype(np.int64)
        E = E + int(J) * p
    return sig, E


def eig_classes(E):
    """Exact eigenspace labels from integer energies."""
    vals, inv = np.unique(E, return_inverse=True)
    sizes = np.bincount(inv)
    return vals, inv, sizes
