"""LANE_O48_B_SEPARATION -- shared helpers.  numpy only (no scipy in this environment).

Nothing here nominates an operator.  Writers are SEARCHED (D-18); clause (iv) is also
checked by the carrier-free criterion Tr(P_E R) = 0 on every eigenspace (C-11 / O-4).
"""
import numpy as np
from itertools import product

I2 = np.eye(2, dtype=complex)
X  = np.array([[0, 1], [1, 0]], dtype=complex)
Y  = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z  = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]

def kron_list(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out

def op_on(n, sites_ops):
    """Operator on n qubits: sites_ops = {site: 2x2 matrix}."""
    return kron_list([sites_ops.get(k, I2) for k in range(n)])

# ---------------------------------------------------------------- F2 Pauli algebra
# A Pauli is (x|z) in F_2^{2n}.  P,Q commute iff  x.z' + z.x' = 0 (mod 2).
def sym(p, q, n):
    s = 0
    for k in range(n):
        s ^= (p[k] & q[n + k]) ^ (p[n + k] & q[k])
    return s

def all_paulis(n):
    """Iterate every Pauli on n qubits up to phase, as a tuple of 2n bits."""
    for v in range(4 ** n):
        bits = tuple((v >> b) & 1 for b in range(2 * n))
        yield bits

def commutes_with_H(p, terms, n):
    """P commutes with H = sum_b coef_b Q_b (Q_b DISTINCT Paulis, coef_b != 0)
       iff P commutes with every Q_b: P Q_b P^dag = (-1)^{sym} Q_b and the Q_b are
       linearly independent, so P H P^dag = H forces every sign to be +1.  EXACT."""
    return all(sym(p, q, n) == 0 for q in terms)

def search_admissible_flippers(n, terms, target, limit=8):
    """SEARCH the full Pauli group on n qubits for U with [U,H]=0 and {U,target}=0.
       Returns up to `limit` of them, plus the total count."""
    found, count = [], 0
    for p in all_paulis(n):
        if sym(p, target, n) != 1:
            continue
        if commutes_with_H(p, terms, n):
            count += 1
            if len(found) < limit:
                found.append(p)
    return found, count

def _single(x, z):
    if x == 0 and z == 0: return I2, "I"
    if x == 1 and z == 0: return X, "X"
    if x == 0 and z == 1: return Z, "Z"
    return Y, "Y"

def pauli_matrix(p, n):
    """Matrix of the F2 Pauli (x|z) with the Hermitian convention (XZ on a site -> Y)."""
    return kron_list([_single(p[k], p[n + k])[0] for k in range(n)])

def pauli_label(p, n):
    return "".join(_single(p[k], p[n + k])[1] for k in range(n))

# ---------------------------------------------------------------- clause machinery
def eigenspaces(H, tol=1e-8):
    w, V = np.linalg.eigh(H)
    out, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol:
            j += 1
        Q = V[:, i:j + 1]
        out.append((w[i], Q @ Q.conj().T, j - i + 1))
        i = j + 1
    return out

def clause_i(R, tol=1e-9):
    n = R.shape[0]
    return (np.linalg.norm(R - R.conj().T) < tol) and (np.linalg.norm(R @ R - np.eye(n)) < tol)

def clause_ii(R, H, Ls=(), tol=1e-9):
    ok = np.linalg.norm(H @ R - R @ H) < tol
    for L in Ls:
        ok = ok and np.linalg.norm(L @ R - R @ L) < tol
    return ok

def clause_iii(R, es, tol=1e-9):
    n = R.shape[0]
    for _, P, m in es:
        M = P @ R @ P
        if np.linalg.norm(M - (np.trace(M) / m) * P) > tol:
            return True
    return False

def clause_iv_trace(R, es, tol=1e-8):
    """C-11: an admissible flipper exists IFF Tr(P_E R) = 0 on EVERY eigenspace."""
    worst = max(abs(np.trace(P @ R)) for _, P, _ in es)
    return worst < tol, float(worst)

# ---------------------------------------------------------------- Walsh / Hadamard
def fwht(a):
    """In-place fast Walsh-Hadamard transform.  a has length 2^n, indexed by the bitmask x.
       Result[S] = sum_x (-1)^{popcount(S & x)} a[x].   With z_k = 1 - 2*bit_k(x) this is
       exactly  sum_z (prod_{k in S} z_k) a(z).   EXACT up to float rounding."""
    a = np.array(a, dtype=np.float64, copy=True)
    h = 1
    n = a.shape[0]
    while h < n:
        for i in range(0, n, h * 2):
            x = a[i:i + h].copy()
            y = a[i + h:i + 2 * h].copy()
            a[i:i + h] = x + y
            a[i + h:i + 2 * h] = x - y
        h *= 2
    return a

def walsh_coeffs(E, m):
    """All Walsh coefficients  c_S = 2^{-m} sum_z chi_S(z) E(z)  from the energy table E[x]."""
    return fwht(E) / (2.0 ** m)

def pair_index(i, j, m):
    return (1 << i) | (1 << j)

# ---------------------------------------------------------------- vectorised Pauli search
def search_admissible_vec(n, terms, target):
    """Vectorised FULL search of the 4^n Pauli group for U with [U,H]=0 and {U,target}=0.
       Paulis are packed as v = x + (z << n).  Commutation of (x|z) with (qx|qz) is the parity
       of popcount(x & qz) ^ popcount(z & qx).  Returns (count, first_example_bits or None)."""
    N = 4 ** n
    v = np.arange(N, dtype=np.int64)
    mask = (1 << n) - 1
    x = v & mask
    z = (v >> n) & mask
    def anti(q):
        qx = sum(q[k] << k for k in range(n))
        qz = sum(q[n + k] << k for k in range(n))
        return (np.bitwise_count(x & qz) ^ np.bitwise_count(z & qx)) & 1
    ok = anti(target) == 1
    for q in terms:
        ok &= (anti(q) == 0)
    cnt = int(ok.sum())
    if cnt == 0:
        return 0, None
    v0 = int(v[ok][0])
    bits = tuple(((v0 >> b) & 1) for b in range(n)) + tuple(((v0 >> (n + b)) & 1) for b in range(n))
    return cnt, bits

def fwht_fast(a):
    """Vectorised FWHT, identical output to fwht() but O(m) numpy passes."""
    a = np.array(a, dtype=np.float64, copy=True)
    n = a.shape[0]; h = 1
    while h < n:
        b = a.reshape(-1, 2, h)
        a = np.concatenate([b[:, 0, :] + b[:, 1, :], b[:, 0, :] - b[:, 1, :]], axis=1).reshape(-1)
        h *= 2
    return a

# ---------------------------------------------------------------- model selection (D-20)
def fit_power_vs_exp(rs, vals, split, floor=None):
    """Fit log10|v| against log10(r) (POWER) and against r (EXPONENTIAL) on r <= split, then
       PREDICT r > split.  Returns in-sample and OUT-OF-SAMPLE rms error in log10 units.
       Points whose |v| falls below `floor` are DROPPED and counted -- the noise floor is the
       model-selection knob (D-20) and must be declared, not hidden."""
    rs = np.asarray(rs, dtype=float); vals = np.asarray(vals, dtype=float)
    a = np.abs(vals)
    keep = np.ones(len(a), bool)
    if floor is not None:
        keep = a > floor
    dropped = int((~keep).sum())
    rs, a = rs[keep], a[keep]
    if len(rs) < 4:
        return dict(n=len(rs), dropped=dropped, ok=False)
    y = np.log10(a)
    ins = rs <= split; oos = ~ins
    if ins.sum() < 3 or oos.sum() < 2:
        return dict(n=len(rs), dropped=dropped, ok=False)
    def lin(xs, ys):
        A = np.vstack([xs, np.ones_like(xs)]).T
        c, *_ = np.linalg.lstsq(A, ys, rcond=None)
        return c
    xp = np.log10(rs)
    cp = lin(xp[ins], y[ins]); ce = lin(rs[ins], y[ins])
    rp_in = y[ins] - (cp[0]*xp[ins] + cp[1]); re_in = y[ins] - (ce[0]*rs[ins] + ce[1])
    rp_out = y[oos] - (cp[0]*xp[oos] + cp[1]); re_out = y[oos] - (ce[0]*rs[oos] + ce[1])
    rms = lambda v: float(np.sqrt(np.mean(v**2)))
    return dict(ok=True, n=len(rs), dropped=dropped,
                p_exponent=float(-cp[0]), xi=float(-1.0/(ce[0]*np.log(10)) if ce[0] < 0 else np.inf),
                pow_in=rms(rp_in), exp_in=rms(re_in),
                pow_out=rms(rp_out), exp_out=rms(re_out),
                n_in=int(ins.sum()), n_out=int(oos.sum()),
                r_in=(float(rs[ins].min()), float(rs[ins].max())),
                r_out=(float(rs[oos].min()), float(rs[oos].max())))
