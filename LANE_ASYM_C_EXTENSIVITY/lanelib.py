"""LANE_ASYM_C_EXTENSIVITY -- shared helpers.

REPRESENTATIONS USED IN THIS LANE, stated once so every table can name its own:
  [PHYS]  the full dense Hilbert space of 4m physical qubits, dim 16^m.  Used for m<=2
          (dim<=256) as the ground truth against which the reduced pictures are checked.
  [F2]    the symplectic (x|z) representation over F_2^{2n}, n=4m.  Exact, combinatorial,
          no Hilbert space built at all.  Used for weights, supports, commutation.
  [CODE]  the restriction of the full dynamics to  codespace (x) bath.  Exact whenever every
          coupling and the Hamiltonian preserve the code space -- which they do here, because
          H is the stabiliser Hamiltonian and every coupling is a LOGICAL operator, so
          codespace (x) bath is an invariant subspace of H_tot and the initial state lies in
          it.  System dim drops 16^m -> 4^m.  VALIDATED against [PHYS] at m=1,2.

NOTHING outside this directory is written or modified.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Ym = 1j * Xm @ Zm

# ------------------------------------------------------------------ [F2] machinery
def sp(a, b, n):
    """symplectic form; 1 iff the two Paulis ANTICOMMUTE"""
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def weight(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

def support(v, n):
    return frozenset(i for i in range(n) if v[i] or v[n + i])

def block_stab_422():
    """the two stabilisers of [[4,2,2]] as (x|z) over F_2^8"""
    return [[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]

def embed(v, b, m, nb=4):
    """embed a length-2*nb (x|z) vector as block b of an m-block carrier (n = m*nb qubits)"""
    n = m * nb
    out = [0] * (2 * n)
    for i in range(nb):
        out[b * nb + i] = v[i]
        out[n + b * nb + i] = v[nb + i]
    return out

def composite_stab(m, nb=4):
    return [embed(s, b, m, nb) for b in range(m) for s in block_stab_422()]

def composite_records_writers(m, nb=4):
    """RECORDS and their conjugate WRITERS for m disjoint [[4,2,2]] blocks, in [F2].

       The single-block logicals are COMPUTED by symplectic_logicals -- never nominated --
       and then embedded blockwise.  Returns (records, writers, n)."""
    pairs = symplectic_logicals(block_stab_422(), 4)
    recs, wrts = [], []
    for b in range(m):
        for (a, c) in pairs:
            recs.append(embed(a, b, m, nb))
            wrts.append(embed(c, b, m, nb))
    return recs, wrts, m * nb

def in_span_f2(v, basis, L):
    v = v[:]
    for bb in basis:
        h = next((i for i in range(L) if bb[i]), None)
        if h is not None and v[h]:
            v = [(x + y) % 2 for x, y in zip(v, bb)]
    return not any(v)

def rref_f2(rows, L):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(L):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

# ------------------------------------------------------------------ [PHYS] machinery
def stab_hamiltonian(m, nb=4):
    """H = -sum_b (X^(x)nb + Z^(x)nb) on block b.  Dense, dim (2^nb)^m."""
    n = m * nb
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for s in composite_stab(m, nb):
        H = H - xz_to_matrix(s, n)
    return H

# ------------------------------------------------------------------ [CODE] machinery
def pauli_on(k, ops):
    """ops: dict qubit->2x2 matrix, on k qubits"""
    M = np.array([[1]], dtype=complex)
    for i in range(k):
        M = np.kron(M, ops.get(i, I2))
    return M

def code_records_couplings(m):
    """[CODE] picture: m blocks -> 2m LOGICAL qubits.  In the code space the chosen record of
       logical qubit i acts exactly as Z_i and its writer as X_i (verified in s1 against
       [PHYS]).  Returns (list of 2m record matrices, dim)."""
    k = 2 * m
    return [pauli_on(k, {i: Zm}) for i in range(k)], 2 ** k

# ------------------------------------------------------------------ cached propagator
class Propagator:
    """One eigendecomposition of H_tot, many times.  Mirrors RecordModel.evolve exactly
       (verified in s1) but caches the eigh so a 25-point time average costs one solve."""
    def __init__(self, HS, env, coupling_list, lam=0.8, state0=None):
        nS, nB = HS.shape[0], env.dim
        HINT = sum(np.kron(A, env.site[j]) for A, j in coupling_list) \
               if coupling_list else np.zeros((nS * nB, nS * nB), dtype=complex)
        Ht = np.kron(HS, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
        self.w, self.U = np.linalg.eigh(Ht)
        if state0 is None:
            state0 = np.eye(nS, dtype=complex) / nS
        self.r0c = self.U.conj().T @ np.kron(state0, env.thermal()) @ self.U
        self.nS, self.env = nS, env

    def state(self, t):
        ph = np.exp(-1j * self.w * t)
        return self.U @ (ph[:, None] * self.r0c * ph.conj()[None, :]) @ self.U.conj().T

TIMES = np.linspace(1.0, 13.0, 25)      # D: time-average chi over ~25 times in [1,13]

def chi_timeavg(prop, records, times=TIMES):
    """time-averaged chi for EVERY record, from one eigendecomposition"""
    acc = np.zeros(len(records))
    for t in times:
        r = prop.state(t)
        for i, R in enumerate(records):
            acc[i] += prop.env.holevo(r, R, prop.nS)
    return acc / len(times)
