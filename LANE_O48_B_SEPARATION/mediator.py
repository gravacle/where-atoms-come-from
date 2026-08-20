"""LANE_O48_B_SEPARATION -- the mediated carrier and its exact solution.

THE CARRIER.  2m qubits: m RECORD qubits (0..m-1) and m MEDIATOR qubits (m..2m-1).
The mediator is NOT a record.  Records couple to it ONLY ON-SITE.

  H = - sum_{i<m-1} [ (t_i/2)(Xa_i Xa_{i+1} + Ya_i Ya_{i+1})
                    + (w_i/2)(Xa_i Xa_{i+1} - Ya_i Ya_{i+1}) ]      <- mediator's OWN dynamics
      - g sum_i  Zr_i Za_i                                          <- LOCAL, ON-SITE coupling
      - mu sum_i Za_i                                               <- optional mediator field

Every record operator Zr_i appears in H only as Zr_i, so [H, Zr_i] = 0 identically: clause (ii)
is exact for any parameters.  H is therefore BLOCK DIAGONAL over the record configuration
z in {+-1}^m, and E0(z) := the lowest energy in block z is an exact function of z.

EXACT SOLUTION.  Jordan-Wigner on the OPEN mediator chain (no string obstruction for
nearest-neighbour terms) turns each block into a quadratic fermion Hamiltonian
   H(z) = sum A_ij c^dag_i c_j + (1/2) sum (B_ij c^dag_i c^dag_j + h.c.) + const
with  A_{i,i+1}=A_{i+1,i} = -t_i,  A_ii = -2 g z_i - 2 mu,  B_{i,i+1} = -w_i = -B_{i+1,i},
const = g sum z_i + mu m.  Since A is real symmetric and B real antisymmetric,
(A+B)^T = A-B, so the Bogoliubov energies are the SINGULAR VALUES of M = A+B, and
   E0(z) = (1/2)(tr A - sum_k sigma_k(M)) + const = -(1/2) sum_k sigma_k(M(z)).
The constants cancel EXACTLY -- verified against dense ED below.  Cost: one m x m SVD per z.
"""
import numpy as np

def Mmat(m, z, t, w, g, mu=0.0):
    """The m x m matrix M = A + B for one record configuration z."""
    z = np.asarray(z, dtype=np.float64)
    M = np.zeros((m, m))
    M[np.arange(m), np.arange(m)] = -2.0 * g * z - 2.0 * mu
    for i in range(m - 1):
        M[i, i + 1] = -t[i] - w[i]
        M[i + 1, i] = -t[i] + w[i]
    return M

def Mbatch(m, Zs, t, w, g, mu=0.0):
    """Stacked M for a batch of configurations Zs of shape (N, m)."""
    N = Zs.shape[0]
    M = np.zeros((N, m, m))
    idx = np.arange(m)
    M[:, idx, idx] = -2.0 * g * Zs - 2.0 * mu
    for i in range(m - 1):
        M[:, i, i + 1] = -t[i] - w[i]
        M[:, i + 1, i] = -t[i] + w[i]
    return M

def E0_batch(m, Zs, t, w, g, mu=0.0, chunk=8192):
    """E0(z) for every configuration in Zs.  EXACT (up to float SVD)."""
    N = Zs.shape[0]
    out = np.empty(N)
    for a in range(0, N, chunk):
        b = min(a + chunk, N)
        s = np.linalg.svd(Mbatch(m, Zs[a:b], t, w, g, mu), compute_uv=False)
        out[a:b] = -0.5 * s.sum(axis=1)
    return out

def all_configs(m):
    x = np.arange(2 ** m, dtype=np.int64)
    return (1 - 2 * ((x[:, None] >> np.arange(m)[None, :]) & 1)).astype(np.float64)

# ------------------------------------------------------------------ dense spin operators
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

def _kron(ops):
    o = np.array([[1.0 + 0j]])
    for a in ops: o = np.kron(o, a)
    return o

def spin_op(nq, d):
    return _kron([d.get(k, I2) for k in range(nq)])

def H_med_dense(m, z, t, w, g, mu=0.0):
    """The mediator block Hamiltonian for a FIXED record configuration z, as a 2^m matrix."""
    H = np.zeros((2 ** m, 2 ** m), dtype=complex)
    for i in range(m - 1):
        XX = spin_op(m, {i: SX, i + 1: SX}); YY = spin_op(m, {i: SY, i + 1: SY})
        H += -(t[i] / 2.0) * (XX + YY) - (w[i] / 2.0) * (XX - YY)
    for i in range(m):
        H += -(g * z[i] + mu) * spin_op(m, {i: SZ})
    return H

def H_full_dense(m, t, w, g, mu=0.0, lr=None):
    """The FULL 2m-qubit Hamiltonian.  Record qubits 0..m-1, mediator qubits m..2m-1.
       lr, if given, is a dict {(i,j): coupling} of an INSERTED record-record term."""
    nq = 2 * m
    H = np.zeros((2 ** nq, 2 ** nq), dtype=complex)
    for i in range(m - 1):
        XX = spin_op(nq, {m + i: SX, m + i + 1: SX}); YY = spin_op(nq, {m + i: SY, m + i + 1: SY})
        H += -(t[i] / 2.0) * (XX + YY) - (w[i] / 2.0) * (XX - YY)
    for i in range(m):
        H += -g * spin_op(nq, {i: SZ, m + i: SZ}) - mu * spin_op(nq, {m + i: SZ})
    if lr:
        for (i, j), c in lr.items():
            H += c * spin_op(nq, {i: SZ, j: SZ})
    return H

def H_full_terms(m, t, w, g, mu=0.0, lr=None):
    """The same H as a list of (F2 Pauli, coefficient) on nq=2m qubits, for the exact Pauli
       commutation test.  Ordering: record qubits 0..m-1, mediator qubits m..2m-1."""
    nq = 2 * m
    def pauli(xs, zs):
        x = [0] * nq; zz = [0] * nq
        for k in xs: x[k] = 1
        for k in zs: zz[k] = 1
        return tuple(x + zz)
    out = []
    for i in range(m - 1):
        cxx = -(t[i] / 2.0) - (w[i] / 2.0)
        cyy = -(t[i] / 2.0) + (w[i] / 2.0)
        if abs(cxx) > 1e-14: out.append((pauli([m + i, m + i + 1], []), cxx))
        if abs(cyy) > 1e-14: out.append((pauli([m + i, m + i + 1], [m + i, m + i + 1]), cyy))
    for i in range(m):
        if abs(g) > 1e-14: out.append((pauli([], [i, m + i]), -g))
        if abs(mu) > 1e-14: out.append((pauli([], [m + i]), -mu))
    if lr:
        for (i, j), c in lr.items():
            if abs(c) > 1e-14: out.append((pauli([], [i, j]), c))
    return out

# ------------------------------------------------------------------ exact O(g^2) coefficient
def chi_free(m, t, fermi=0.0):
    """EXACT second-order (O(g^2)) two-body Walsh coefficient for a hopping-only mediator.

    The perturbation is V = -2g sum_i z_i n_i.  Standard second-order perturbation theory over
    the BARE mediator's particle-hole excitations gives
        dE^(2) = -4 g^2 sum_{ij} z_i z_j T_ij,
        T_ij   = sum_{p occupied, q empty} phi_p(i) phi_q(i) phi_p(j) phi_q(j) / (eps_q - eps_p),
    so the two-body Walsh coefficient is  J_eff(i,j) = -8 g^2 T_ij  for i<j.
    NO BACKGROUND ENTERS: at this order the record configuration is what is being expanded in,
    so T is a property of the mediator alone.  Validated against the exact Walsh transform below.

    Returns T (m x m).  Half filling (the state clause (iv) forces) = all eps < 0 occupied.
    """
    A0 = np.zeros((m, m))
    for i in range(m - 1):
        A0[i, i + 1] = A0[i + 1, i] = -t[i]
    eps, phi = np.linalg.eigh(A0)
    occ = np.where(eps < fermi - 1e-12)[0]
    emp = np.where(eps > fermi + 1e-12)[0]
    if len(occ) + len(emp) != m:
        raise RuntimeError("zero mode at the Fermi level: %d occ %d emp of %d" % (len(occ), len(emp), m))
    # W[(p,q), i] = phi_p(i) phi_q(i);  T = W^T diag(1/(eps_q-eps_p)) W
    P = phi[:, occ]                      # (m, no)
    Q = phi[:, emp]                      # (m, ne)
    no, ne = P.shape[1], Q.shape[1]
    d = 1.0 / (eps[emp][None, :] - eps[occ][:, None])       # (no, ne)
    W = (P[:, :, None] * Q[:, None, :]).reshape(m, no * ne)  # (m, no*ne)
    T = (W * d.reshape(-1)[None, :]) @ W.T
    return T

def mediator_gap(m, t, w):
    """The BARE mediator's own spectral gap: the smallest Bogoliubov energy at g=0, mu=0."""
    M = Mmat(m, np.zeros(m), t, w, 0.0, 0.0)
    return float(np.linalg.svd(M, compute_uv=False).min())

def chi_row(m, t, i0, fermi=0.0):
    """One ROW of T (see chi_free), at O(m^3) flops and O(m^2) memory instead of O(m^4)/O(m^3).
       T[i0, j] = sum_{p occ, q emp} phi_p(i0) phi_q(i0) phi_p(j) phi_q(j) / (eps_q - eps_p).
       Identical mathematics to chi_free -- cross-checked against it in the scripts."""
    A0 = np.zeros((m, m))
    for i in range(m - 1):
        A0[i, i + 1] = A0[i + 1, i] = -t[i]
    eps, phi = np.linalg.eigh(A0)
    occ = np.where(eps < fermi - 1e-12)[0]
    emp = np.where(eps > fermi + 1e-12)[0]
    if len(occ) + len(emp) != m:
        raise RuntimeError("zero mode at the Fermi level")
    P_ = phi[:, occ]; Q_ = phi[:, emp]
    u = np.outer(P_[i0, :], Q_[i0, :]) / (eps[emp][None, :] - eps[occ][:, None])
    return ((P_ @ u) * Q_).sum(axis=1)

def chi_row_general(A0, i0, fermi=0.0):
    """chi_row for an ARBITRARY hopping matrix A0 (any graph, any dimension)."""
    m = A0.shape[0]
    eps, phi = np.linalg.eigh(A0)
    occ = np.where(eps < fermi - 1e-10)[0]
    emp = np.where(eps > fermi + 1e-10)[0]
    if len(occ) + len(emp) != m:
        raise RuntimeError("zero mode at the Fermi level: %d + %d != %d"
                           % (len(occ), len(emp), m))
    P_ = phi[:, occ]; Q_ = phi[:, emp]
    u = np.outer(P_[i0, :], Q_[i0, :]) / (eps[emp][None, :] - eps[occ][:, None])
    return ((P_ @ u) * Q_).sum(axis=1)

def E0_fields(A0, B, h):
    """E0 of the quadratic mediator with on-site fields h (h_i = g z_i), EXACT.
       M = A + B with A = A0 + diag(-2h);  E0 = -(1/2) sum sigma(M) (constants cancel, see above)."""
    A = A0 + np.diag(-2.0 * np.asarray(h, dtype=float))
    return -0.5 * float(np.linalg.svd(A + B, compute_uv=False).sum())

def j_eff_fd(A0, B, i0, js, delta=1e-2):
    """The exact O(g^2) two-body coefficient by a CENTRAL SECOND DIFFERENCE of E0 in the fields,
       taken at ZERO background so the mediator stays at half filling and is never doped.
       Returns d^2 E0 / dh_i dh_j, which equals J_eff(i,j)/g^2.  Works for ANY quadratic mediator,
       pairing included, where the closed-form orbital sum does not apply."""
    m = A0.shape[0]
    out = []
    for j in js:
        vals = []
        for si in (+1, -1):
            for sj in (+1, -1):
                h = np.zeros(m); h[i0] = si * delta; h[j] = sj * delta
                vals.append(E0_fields(A0, B, h))
        out.append((vals[0] - vals[1] - vals[2] + vals[3]) / (4 * delta * delta))
    return np.array(out)

def pair_B(m, w):
    B = np.zeros((m, m))
    for i in range(m - 1):
        B[i, i + 1] = -w[i]; B[i + 1, i] = +w[i]
    return B

def hop_A(m, t):
    A = np.zeros((m, m))
    for i in range(m - 1):
        A[i, i + 1] = A[i + 1, i] = -t[i]
    return A

def square_lattice_A(Lx, Ly, t=1.0):
    """Open Lx x Ly square lattice, nearest-neighbour hopping.  Bipartite, so the half-filled
       state is particle-hole symmetric and clause (iv)'s global flip survives."""
    m = Lx * Ly
    A = np.zeros((m, m))
    idx = lambda x, y: x * Ly + y
    for x in range(Lx):
        for y in range(Ly):
            if x + 1 < Lx: A[idx(x, y), idx(x + 1, y)] = A[idx(x + 1, y), idx(x, y)] = -t
            if y + 1 < Ly: A[idx(x, y), idx(x, y + 1)] = A[idx(x, y + 1), idx(x, y)] = -t
    return A, idx
