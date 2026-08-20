"""LANE_EXACT_C_FORM shared utilities.

Nothing here is nominated: every logical operator used downstream is DERIVED from
symplectic_logicals() + the stabiliser group by F_2 linear algebra, and then has all
five record clauses VERIFIED before it is used.

Exact objects use integers / F_2 / fractions.  Floats appear only for Holevo chi, and
every float claim downstream carries a noise floor.
"""
import sys, itertools, numpy as np
from fractions import Fraction
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import (RecordModel, Environment, symplectic_logicals,
                          xz_to_matrix, commutant, eigenspaces)

# ----------------------------------------------------------------- F_2 machinery
def sp(a, b, n):
    """symplectic form on (x|z) in F_2^{2n}: 1 iff the two Paulis ANTICOMMUTE"""
    return sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n)) % 2

def rref(rows, width):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(width):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x+y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

def in_span(v, basis, width):
    R, piv = rref(basis, width)
    v = v[:]
    for i, c in enumerate(piv):
        if v[c]: v = [(x+y) % 2 for x, y in zip(v, R[i])]
    return not any(v)

def f2_rank(rows, width):
    return len(rref(rows, width)[1])

# ----------------------------------------------------------------- carriers (as F_2 data)
def stab_nn2(n):
    """[[n,n-2,2]]: X^(x)n and Z^(x)n as (x|z)."""
    return [[1]*n + [0]*n, [0]*n + [1]*n]

def stab_blocks(m, bs=4):
    """m independent [[4,2,2]] blocks on n = 4m qubits."""
    n = bs*m; out = []
    for b in range(m):
        x = [0]*(2*n); z = [0]*(2*n)
        for q in range(b*bs, (b+1)*bs):
            x[q] = 1; z[n+q] = 1
        out += [x, z]
    return out

def pauli_vec(n, xs=(), zs=()):
    v = [0]*(2*n)
    for i in xs: v[i] ^= 1
    for i in zs: v[n+i] ^= 1
    return v

# ----------------------------------------------------------------- derived logicals
def derived_logical_span(stab, n):
    """N(S) as an F_2 basis, DERIVED: stabilisers + the conjugate pairs returned by
       symplectic_logicals.  Returns (S_basis, L_basis, pairs)."""
    pairs = symplectic_logicals(stab, n)
    L = []
    for a, b in pairs: L += [a, b]
    return list(stab), L, pairs

def is_nontrivial_logical(v, stab, Lb, n):
    """v is in N(S) and NOT in S  -- decided exactly in F_2 against the DERIVED basis."""
    if any(sp(v, s, n) for s in stab): return False          # not in N(S)
    if in_span(v, stab, 2*n): return False                   # in S: trivial
    return in_span(v, stab + Lb, 2*n)                        # in N(S) = <S, L>

# ----------------------------------------------------------------- code space of [[n,n-2,2]]
def code_reps(n):
    """Representatives of the +1,+1 stabiliser sector of [[n,n-2,2]].
       BIT CONVENTION, matching xz_to_matrix: qubit i is bit (n-1-i) of the integer index,
       i.e. qubit 0 is the MOST significant bit.  (A first pass used the opposite convention
       and T0(c) caught it -- logged, not hidden.)
       Basis vector for rep v is (|v> + |v-bar>)/sqrt2 with |v| even and qubit0(v) = 0.
       Exactly 2^{n-2} of them."""
    def q0(v): return (v >> (n-1)) & 1
    reps = [v for v in range(1 << n) if bin(v).count('1') % 2 == 0 and not q0(v)]
    return reps, {v: i for i, v in enumerate(reps)}

def qbit(v, i, n):
    return (v >> (n-1-i)) & 1

def _canon(v, n):
    return v if not ((v >> (n-1)) & 1) else (~v) & ((1 << n) - 1)

def compress_XX(i, j, n, reps, idx):
    """X_i X_j compressed to the code space -- an exact 0/1 PERMUTATION matrix."""
    d = len(reps); M = np.zeros((d, d), dtype=np.int64)
    mask = (1 << (n-1-i)) | (1 << (n-1-j))
    for a, v in enumerate(reps):
        M[idx[_canon(v ^ mask, n)], a] = 1
    return M

def compress_ZZ(i, j, n, reps, idx):
    """Z_i Z_j compressed to the code space -- an exact diagonal +-1 matrix."""
    d = len(reps); M = np.zeros((d, d), dtype=np.int64)
    for a, v in enumerate(reps):
        s = qbit(v, i, n) ^ qbit(v, j, n)
        M[a, a] = -1 if s else 1
    return M

def compress_pauli(v, n, reps, idx):
    """General logical (even X-weight, even Z-weight) compressed exactly.
       Hermitian phase i^{x.z} applied, matching xz_to_matrix's per-site factor i*X*Z."""
    xs = [i for i in range(n) if v[i]]
    zs = [i for i in range(n) if v[n+i]]
    xm = 0
    for i in xs: xm |= (1 << (n-1-i))
    d = len(reps); M = np.zeros((d, d), dtype=complex)
    nY = sum(1 for i in range(n) if v[i] and v[n+i])
    phase = (1j) ** (nY % 4)
    for a, u in enumerate(reps):
        w = u ^ xm
        s = sum(qbit(u, i, n) for i in zs) % 2          # X^x Z^z acting on |u>
        M[idx[_canon(w, n)], a] = phase * (-1 if s else 1)
    return M

# ----------------------------------------------------------------- evolution helper
def evolve_times(Hs, env, terms, lam, times, state0):
    """ONE eigendecomposition, many times, STREAMED (a first version stored all 25 density
       matrices and drove the machine into swap at dim 4096 -- logged, not hidden).
       Mirrors RecordModel.evolve exactly; verified against it in t0_clauses.py."""
    nS = Hs.shape[0]; nB = env.dim
    HINT = sum(np.kron(np.asarray(A, dtype=complex), env.site[j % env.nq]) for A, j in terms)
    Ht = np.kron(Hs, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam*HINT
    w, U = np.linalg.eigh(Ht)
    r0 = np.kron(state0, env.thermal())
    Uc = U.conj().T @ r0 @ U
    for t in times:
        ph = np.exp(-1j*w*t)
        yield U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T

TIMES = np.linspace(1.0, 13.0, 25)

def chi_avg(Hs, env, terms, lam, readouts, state0, times=TIMES):
    """Time-averaged Holevo chi for each readout record.  Returns a list of floats."""
    nS = Hs.shape[0]
    acc = np.zeros(len(readouts))
    ro = [np.asarray(R, dtype=complex) for R in readouts]
    for r in evolve_times(Hs, env, terms, lam, times, state0):
        for k, R in enumerate(ro):
            acc[k] += env.holevo(r, R, nS)
        del r
    return list(acc/len(times))
