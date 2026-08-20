"""Code-space reduction for the operational relation, with the self-checks that make it usable.

WHY A REDUCTION.  The couplings used below are the records themselves, and every record
preserves the code space, so the whole evolution from a code-space initial state stays inside
code space (x) bath.  Working in the 2^k-dimensional code space instead of the 2^n-dimensional
physical space is therefore EXACT, not an approximation -- and it is what makes n = 10 reachable.

SELF-CHECKS (all must pass or nothing is reported):
  SC-8   the isometry W satisfies W-dagger W = I_{2^k} and W W-dagger = P_g
  SC-9   W-dagger Zbar_i W = Z_i and W-dagger Xbar_i W = X_i  exactly (standard k-qubit Paulis)
  SC-10  W-dagger H W = -2 I
  SC-11  chi computed in the reduced space equals chi computed by RecordModel.evolve +
         Environment.holevo in the FULL 2^n space, to 1e-10
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np

def kron_list(ms):
    M = np.array([[1]], dtype=complex)
    for m in ms: M = np.kron(M, m)
    return M

def std_pauli(k, q, letter):
    return kron_list([({'X': Xm, 'Z': Zm, 'I': I2}[letter] if i == q else I2) for i in range(k)])

def code_isometry(n, pairs):
    """W : C^{2^k} -> C^{2^n}, an isometry onto the code space, in which the logicals become
       the STANDARD k-qubit Paulis."""
    k = len(pairs)
    H = hamiltonian(n)
    w, V = np.linalg.eigh(H)
    kd = int((np.abs(w - w[0]) < 1e-9).sum())
    assert kd == 2 ** k, (kd, 2 ** k)
    Q = V[:, :kd]
    Xb = [xz_to_matrix(list(x), n) for x, z in pairs]
    Zb = [xz_to_matrix(list(z), n) for x, z in pairs]
    Zr = [Q.conj().T @ Z @ Q for Z in Zb]
    Xr = [Q.conj().T @ X @ Q for X in Xb]
    # simultaneous eigenvector of every Zr with all +1: project
    P = np.eye(kd, dtype=complex)
    for Z in Zr: P = P @ ((np.eye(kd) + Z) / 2)
    wv, Vv = np.linalg.eigh((P + P.conj().T) / 2)
    ref = Vv[:, -1]
    ref = ref / np.linalg.norm(ref)
    cols = np.zeros((kd, kd), dtype=complex)
    for b in range(kd):
        v = ref.copy()
        for i in range(k):
            if (b >> (k - 1 - i)) & 1:
                v = Xr[i] @ v
        cols[:, b] = v
    W = Q @ cols
    return W, Xb, Zb

def isometry_checks(W, n, pairs, tol=1e-8):
    k = len(pairs); kd = 2 ** k
    H = hamiltonian(n)
    out = {}
    out["SC-8a"] = float(np.linalg.norm(W.conj().T @ W - np.eye(kd)))
    wv, Vv = np.linalg.eigh(H)
    Q = Vv[:, :kd]; Pg = Q @ Q.conj().T
    out["SC-8b"] = float(np.linalg.norm(W @ W.conj().T - Pg))
    e = 0.0
    for i, (x, z) in enumerate(pairs):
        Xr = W.conj().T @ xz_to_matrix(list(x), n) @ W
        Zr = W.conj().T @ xz_to_matrix(list(z), n) @ W
        e = max(e, float(np.linalg.norm(Xr - std_pauli(k, i, 'X'))),
                float(np.linalg.norm(Zr - std_pauli(k, i, 'Z'))))
    out["SC-9"] = e
    out["SC-10"] = float(np.linalg.norm(W.conj().T @ H @ W + 2 * np.eye(kd)))
    return out

# ---------------------------------------------------------------- fast chi on effective qubits
def chi_fast(rho, nS, nB, k, q, letter):
    """Holevo chi(record : whole bath) for a record that is a STANDARD single-effective-qubit
       Pauli (X or Z on effective qubit q).  Exact; equals Environment.holevo (SC-11)."""
    if letter == 'X':
        Hd = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        t = rho.reshape([2] * k + [nB] + [2] * k + [nB])
        t = np.tensordot(Hd, t, axes=([1], [q]))
        t = np.moveaxis(t, 0, q)
        t = np.tensordot(t, Hd.conj().T, axes=([k + 1 + q], [0]))
        t = np.moveaxis(t, -1, k + 1 + q)
        rho = t.reshape(nS * nB, nS * nB)
    a = 2 ** q; b = 2 ** (k - 1 - q)
    t = rho.reshape(a, 2, b, nB, a, 2, b, nB)
    outs = []
    for s in (0, 1):
        sub = t[:, s, :, :, :, s, :, :]
        rB = np.einsum('abiabj->ij', sub)
        p = float(np.real(np.trace(rB)))
        if p < 1e-12: continue
        outs.append((p, rB / p))
    if len(outs) < 2: return 0.0
    av = sum(p * r for p, r in outs)
    def vn(r):
        e = np.linalg.eigvalsh((r + r.conj().T) / 2); e = e[e > 1e-13]
        return float(-(e * np.log2(e)).sum())
    return max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)

def evolve_cached(Heff, env, HINT, lam, times, state0):
    """One eigendecomposition of the total Hamiltonian, reused across all times."""
    nS = Heff.shape[0]; nB = env.dim
    Ht = np.kron(Heff, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
    w, U = np.linalg.eigh(Ht)
    r0 = np.kron(state0, env.thermal())
    Uc = U.conj().T @ r0 @ U
    for t in times:
        ph = np.exp(-1j * w * t)
        yield U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
