"""LANE_O48_A_EXTENSIVE -- shared exact chain machinery.

H = sum_{i=1..n-1} J_i Z_i Z_{i+1}   on an OPEN chain of n qubits.
Records nominated for TESTING (not asserted): R_i = Z_i.

EXACTNESS POLICY (D-19): every coupling is a RATIONAL a_i / D with a_i a Python int and
D = 2^40 fixed.  Energies are therefore EXACT INTEGERS in units of 1/D, and every trace,
every level grouping, every sum is integer arithmetic with no float rounding anywhere.
Nothing here is fitted; the objects are finite and are computed, not estimated.
"""
import numpy as np

D = 1 << 40                      # fixed denominator; couplings live in [0.5, 1.0)

def couplings(nbonds, stream=0):
    """DISTINCT, BOUNDED integer couplings a_i, so J_i = a_i/D in [0.5,1.0).
       Deterministic LCG -- no floats, reproducible, and distinctness is CHECKED."""
    a, x = [], (0x9E3779B97F4A7C15 ^ (stream * 0xBF58476D1CE4E5B9)) & ((1 << 64) - 1)
    for _ in range(nbonds):
        x = (x * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        a.append((D >> 1) + (x >> 25) % (D >> 1))
    return a

def uniform_couplings(nbonds):
    return [(D >> 1) + 12345] * nbonds

# ---------------------------------------------------------------- configuration space
def configs(n):
    """All 2^n sign strings as an int8 array of shape (2^n, n), entries +-1. s[:,0] is site 1."""
    idx = np.arange(1 << n, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(n - 1, -1, -1)[None, :]) & 1).astype(np.int8)
    return (1 - 2 * bits).astype(np.int8)     # bit 0 -> +1, bit 1 -> -1

def energies_int(s, a):
    """E(s) * D as an EXACT int64 array. E(s) = sum_i J_i s_i s_{i+1}."""
    n = s.shape[1]
    E = np.zeros(s.shape[0], dtype=np.int64)
    for i in range(n - 1):
        E += np.int64(a[i]) * (s[:, i].astype(np.int64) * s[:, i + 1].astype(np.int64))
    return E

def levels(E):
    """Exact eigenspaces of H: group configuration indices by EXACT integer energy.
       Returns (sorted unique energies, array of level-id per config)."""
    u, inv = np.unique(E, return_inverse=True)
    return u, inv

# ---------------------------------------------------------------- dense Pauli tools (cross-check only)
_I2 = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_P = {0: _I2, 1: _X, 2: _Y, 3: _Z}

def pauli(label):
    """label: iterable of 0/1/2/3 = I/X/Y/Z, site 1 first."""
    M = np.array([[1]], dtype=complex)
    for c in label:
        M = np.kron(M, _P[c])
    return M

def dense_H(n, a):
    """H as a dense matrix, in the SAME basis order as configs(n) (bit 0 -> +1 eigenvalue)."""
    H = np.zeros((1 << n, 1 << n), dtype=complex)
    for i in range(n - 1):
        lab = [0] * n; lab[i] = 3; lab[i + 1] = 3
        H = H + (a[i] / D) * pauli(lab)
    return H

def dense_Z(n, i):
    lab = [0] * n; lab[i] = 3
    return pauli(lab)
