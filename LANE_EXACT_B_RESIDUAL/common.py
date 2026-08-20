"""LANE_EXACT_B_RESIDUAL -- shared machinery.

CARRIER: [[n, n-2, 2]], n even.  Stabilisers X^(x)n and Z^(x)n.  H = -(X^(x)n + Z^(x)n).
Logicals are COMPUTED by symplectic_logicals -- NEVER nominated (five failures in this program).

THE EXACT REDUCTION used throughout (proved and then VERIFIED numerically in b1):
  * every logical operator commutes with both stabilisers, hence with H, hence the joint
    evolution H (x) I + I (x) H_B + lam*sum_i Rbar_i (x) X_{site(i)} preserves each syndrome
    sector.  The initial state Pg/k lives entirely in the code (= ground) sector, so the whole
    computation may be carried out on the 2^(n-2)-dimensional code space WITHOUT approximation.
  * on the code space H = -2*I, a c-number, so H_S drops out of the dynamics entirely.
"""
import numpy as np, itertools, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
from record_model import Environment, symplectic_logicals, xz_to_matrix

I2 = np.eye(2)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)

# ------------------------------------------------------------------ F_2 symplectic bookkeeping
def sp(a, b, n):
    """symplectic form on (x|z) in F_2^{2n}: 1 iff the two Paulis ANTICOMMUTE."""
    return (sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n))) % 2

def xor(a, b):
    return [(u + v) % 2 for u, v in zip(a, b)]

def stabilisers(n):
    return [[1] * n + [0] * n, [0] * n + [1] * n]

def logical_pairs(n):
    """COMPUTED conjugate pairs [(X_i, Z_i), ...]; returns a LIST OF PAIRS, not two lists."""
    return symplectic_logicals(stabilisers(n), n)

def check_symplectic(pairs, n):
    """SELF-CHECK: <X_i,Z_i>=1, every other pairing 0.  Returns (ok, detail)."""
    K = len(pairs)
    bad = []
    for i in range(K):
        if sp(pairs[i][0], pairs[i][1], n) != 1:
            bad.append(f"<X{i},Z{i}>!=1")
        for j in range(K):
            if i == j: continue
            for (a, la), (b, lb) in itertools.product(
                    [(pairs[i][0], f"X{i}"), (pairs[i][1], f"Z{i}")],
                    [(pairs[j][0], f"X{j}"), (pairs[j][1], f"Z{j}")]):
                if sp(a, b, n) != 0:
                    bad.append(f"<{la},{lb}>!=0")
    return (len(bad) == 0, bad)

# ------------------------------------------------------------------ the code space, EXACTLY
def code_basis(n):
    """Orthonormal basis of the [[n,n-2,2]] code space, built exactly (no eigensolver):
       Z^(x)n = +1 selects EVEN-weight bitstrings; X^(x)n = +1 pairs x with its complement.
       Basis vectors (|x> + |xbar>)/sqrt2 for even-weight x with x < xbar.  dim = 2^(n-2)."""
    full = 1 << n
    mask = full - 1
    cols = []
    for x in range(full):
        if bin(x).count('1') % 2: continue
        xb = x ^ mask
        if x >= xb: continue
        cols.append((x, xb))
    Q = np.zeros((full, len(cols)), dtype=complex)
    for c, (x, xb) in enumerate(cols):
        Q[x, c] = Q[xb, c] = 1 / np.sqrt(2.0)
    return Q

def reduce_op(P, Q):
    """restrict an operator that PRESERVES the code space to it: Q^dag P Q."""
    return Q.conj().T @ P @ Q

# ------------------------------------------------------------------ chi, the fast exact path
def chi_times(sysops_sites, Rread, env, lam, times, nS=None, HS=None, state0=None):
    """Time-averaged Holevo chi(Rread : whole bath).

       sysops_sites: list of (system operator, bath-site index).  Coupling is
           lam * sum_i A_i (x) X_{site_i}  -- the DISTRIBUTED form the lanes use.
       HS: system Hamiltonian; None means the c-number 0 (correct on the code space).
       state0: system state; None means maximally mixed on the supplied system space.

       Uses ONE eigendecomposition for all times and never forms the full joint state:
           r(t) = W Y W^dag,  W = U diag(e^{-i w t}),  Y = U^dag r0 U
           Tr_S[(Pi_s (x) I) r] = sum_{a,g} (Pi_s A)_{(a,b),g} conj(W)_{(a,b'),g},  A = W Y
       which is what env.holevo consumes.  b1 checks this against env.holevo on the full r."""
    nS = Rread.shape[0] if nS is None else nS
    nB = env.dim
    D = nS * nB
    HINT = sum(np.kron(A, env.site[j % env.nq]) for A, j in sysops_sites)
    Ht = np.kron(np.eye(nS), env.HB) + lam * HINT
    if HS is not None:
        Ht = Ht + np.kron(HS, np.eye(nB))
    w, U = np.linalg.eigh(Ht)
    r0 = np.kron(np.eye(nS) / nS if state0 is None else state0, env.thermal())
    Y = U.conj().T @ r0 @ U
    Uc = U.conj()
    out = []
    for t in times:
        ph = np.exp(-1j * w * t)
        W = U * ph[None, :]
        A = W @ Y
        Wc = Uc * ph.conj()[None, :]
        chi = _holevo_from_AW(A, Wc, Rread, nS, nB)
        out.append(chi)
    return np.array(out)

def _vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def _holevo_from_AW(A, Wc, R, nS, nB):
    D = nS * nB
    A3 = A.reshape(nS, nB, D)
    W3 = Wc.reshape(nS, nB, D)
    res = []
    for s in (+1, -1):
        Pi = (np.eye(nS) + s * R) / 2
        Ap = np.einsum('ac,cbg->abg', Pi, A3, optimize=True)
        rB = np.einsum('abg,acg->bc', Ap, W3, optimize=True)
        p = float(np.real(np.trace(rB)))
        if p < 1e-12: continue
        res.append((p, rB / p))
    if len(res) < 2: return 0.0
    av = sum(p * rb for p, rb in res)
    return max(_vN(av) - sum(p * _vN(rb) for p, rb in res), 0.0)

# ------------------------------------------------------------------ venues
class Venue:
    def __init__(self, name, energies, times, beta=2.0):
        self.name, self.energies, self.times, self.beta = name, tuple(energies), np.asarray(times), beta
    def env(self, nq):
        return Environment(nq, energies=self.energies[:nq], beta=self.beta)

BASE = Venue("V0-baseline", (1.0, 1.4, 0.7, 1.2), np.linspace(1.0, 13.0, 25))

def say(*a):
    print(*a); sys.stdout.flush()
