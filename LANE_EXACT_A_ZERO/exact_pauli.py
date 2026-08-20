"""EXACT PAULI ENGINE -- integers only, no floats anywhere.

A signed Pauli is (m, a) meaning  i^m * W(a),  m in Z_4,  a = (x|z) in F_2^{2n}.

CONVENTION, chosen to agree EXACTLY with record_model.xz_to_matrix:
    xz_to_matrix maps  (x_j,z_j) -> I, X, Z, i*X*Z = Y   per qubit.
    Hence  W(a) = tensor_j sigma(x_j,z_j) = i^{x.z} X^x Z^z   with x.z the INTEGER dot product.

Multiplication.  X^x Z^z X^x' Z^z' = (-1)^{z.x'} X^{x xor x'} Z^{z xor z'}   (integer dot, mod 2)
so
    W(a) W(b) = i^{phi(a,b)} W(a xor b),
    phi(a,b) = x.z + x'.z' - (x xor x').(z xor z') + 2 (z.x')   (mod 4).

Everything below is exact integer arithmetic.  No tolerance, no float, no sampling.
"""
from itertools import product as iproduct


# ------------------------------------------------------------------ core algebra
def dot(u, v):
    return sum(ui * vi for ui, vi in zip(u, v))


def xr(u, v):
    return [(ui + vi) % 2 for ui, vi in zip(u, v)]


def split(a, n):
    return a[:n], a[n:]


def phi(a, b, n):
    """Exact Z_4 exponent with W(a)W(b) = i^phi W(a xor b)."""
    x, z = split(a, n)
    xp, zp = split(b, n)
    s = xr(x, xp)
    t = xr(z, zp)
    return (dot(x, z) + dot(xp, zp) - dot(s, t) + 2 * dot(z, xp)) % 4


def sp(a, b, n):
    """Symplectic form over F_2.  0 = commute, 1 = anticommute."""
    x, z = split(a, n)
    xp, zp = split(b, n)
    return (dot(x, zp) + dot(z, xp)) % 2


def pmul(p, q, n):
    """(m,a)*(m',b) -> (m'',c) exactly."""
    m, a = p
    mp, b = q
    return ((m + mp + phi(a, b, n)) % 4, xr(a, b))


def pdag(p, n):
    """Adjoint of i^m W(a) is i^{-m} W(a) since W(a) is Hermitian."""
    m, a = p
    return ((-m) % 4, a[:])


def is_hermitian(p, n):
    return p[0] % 2 == 0 and (p[0] % 4 in (0, 2))  # i^m real  <=> m even


def pidentity(n):
    return (0, [0] * (2 * n))


def is_identity(p):
    return p[0] % 4 == 0 and not any(p[1])


def is_zero_free(p):
    return not any(p[1])


# ------------------------------------------------------------------ integer intersection lifts
def local_anticommute_sites(a, b, n):
    """Sites where the two single-qubit Paulis anticommute.  Its PARITY is sp(a,b)."""
    x, z = split(a, n)
    xp, zp = split(b, n)
    return [j for j in range(n) if (x[j] * zp[j] + z[j] * xp[j]) % 2 == 1]


def I_unsigned(a, b, n):
    """Integer lift 1: NUMBER of sites at which the Paulis locally anticommute."""
    return len(local_anticommute_sites(a, b, n))


def I_signed(a, b, n):
    """Integer lift 2: signed intersection sum_j (x_j z'_j - z_j x'_j) over Z."""
    x, z = split(a, n)
    xp, zp = split(b, n)
    return sum(x[j] * zp[j] - z[j] * xp[j] for j in range(n))


def I_overlap(a, b, n):
    """Integer lift 3: number of sites where BOTH act non-trivially (support overlap)."""
    x, z = split(a, n)
    xp, zp = split(b, n)
    return sum(1 for j in range(n) if (x[j] or z[j]) and (xp[j] or zp[j]))


# ------------------------------------------------------------------ stabiliser / code data
def signed_stabiliser_group(stabs, n):
    """stabs: list of (x|z) vectors whose W() we DECLARE to be +1-stabilisers.
       Returns the full signed group as a dict a_tuple -> m, with i^m W(a) the group element
       that acts as +1 on the code space.  Exact."""
    grp = {tuple([0] * (2 * n)): 0}
    frontier = [(0, s[:]) for s in stabs]
    changed = True
    gens = [(0, s[:]) for s in stabs]
    while changed:
        changed = False
        for g in list(grp.items()):
            gp = (g[1], list(g[0]))
            for h in gens:
                r = pmul(gp, h, n)
                key = tuple(r[1])
                if key not in grp:
                    grp[key] = r[0]
                    changed = True
                else:
                    assert grp[key] == r[0], "stabiliser group is not consistently signed"
    return grp


def ground_projector_trace_ratio(p, stab_grp, n):
    """Tr(Pi * i^m W(a)) / Tr(Pi) for Pi the projector onto the joint +1 eigenspace of the
       signed stabiliser group.  EXACT: returns a Z_4 exponent, or None when the trace is 0.

       Pi = (1/|S|) sum_{s in S} s.  Tr(Pi P) / Tr(Pi) = (1/|S|) sum_s Tr(s P)/Tr(Pi).
       Tr(W(c)) = 0 unless c = 0.  So only the single s with a_s = a contributes, and then
       s * P = i^{...} I whose trace is 2^n * i^{...}, while Tr(Pi) = 2^n/|S|.
       Result: i^{(m_s + m + phi)} where s = (m_s, a)."""
    m, a = p
    key = tuple(a)
    if key not in stab_grp:
        return None            # EXACTLY ZERO
    ms = stab_grp[key]
    # s * p = i^{ms + m + phi(a,a)} W(0) = i^{ms+m} I  (phi(a,a)=0 always)
    e = (ms + m + phi(a, a, n)) % 4
    return e                    # value is i^e, exactly non-zero


def zint_str(e):
    """Render i^e exactly."""
    return {0: "1", 1: "i", 2: "-1", 3: "-i"}[e % 4]


# ------------------------------------------------------------------ float-free cross check
def to_gaussian_matrix(p, n):
    """EXACT matrix over Z[i] as nested tuples of (re, im) integer pairs.  Only for small n."""
    m, a = p
    x, z = split(a, n)
    # single-qubit factors as 2x2 Gaussian-integer matrices
    I2 = (((1, 0), (0, 0)), ((0, 0), (1, 0)))
    X2 = (((0, 0), (1, 0)), ((1, 0), (0, 0)))
    Z2 = (((1, 0), (0, 0)), ((0, 0), (-1, 0)))
    Y2 = (((0, 0), (0, -1)), ((0, 1), (0, 0)))          # i*X*Z
    M = (((1, 0),),)
    for j in range(n):
        f = I2 if (x[j], z[j]) == (0, 0) else (X2 if (x[j], z[j]) == (1, 0)
                                               else (Z2 if (x[j], z[j]) == (0, 1) else Y2))
        M = gkron(M, f)
    ph = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}[m % 4]
    return gscale(M, ph)


def gmul_s(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gadd_s(u, v):
    return (u[0] + v[0], u[1] + v[1])


def gkron(A, B):
    ra, ca = len(A), len(A[0])
    rb, cb = len(B), len(B[0])
    return tuple(tuple(gmul_s(A[i // rb][j // cb], B[i % rb][j % cb])
                       for j in range(ca * cb)) for i in range(ra * rb))


def gscale(A, s):
    return tuple(tuple(gmul_s(e, s) for e in row) for row in A)


def gmatmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return tuple(tuple(_gsum(gmul_s(A[i][t], B[t][j]) for t in range(k)) for j in range(m))
                 for i in range(n))


def _gsum(it):
    r = (0, 0)
    for v in it:
        r = gadd_s(r, v)
    return r


def gsub(A, B):
    return tuple(tuple((A[i][j][0] - B[i][j][0], A[i][j][1] - B[i][j][1])
                       for j in range(len(A[0]))) for i in range(len(A)))


def gadd(A, B):
    return tuple(tuple(gadd_s(A[i][j], B[i][j]) for j in range(len(A[0]))) for i in range(len(A)))


def gis_zero(A):
    return all(e == (0, 0) for row in A for e in row)


def gfrob2(A):
    """Exact integer squared Frobenius norm."""
    return sum(e[0] * e[0] + e[1] * e[1] for row in A for e in row)


def gtrace(A):
    return _gsum(A[i][i] for i in range(len(A)))


# ------------------------------------------------------------------ EXACT numpy int64 backend
# numpy int64 arrays are EXACT INTEGER arithmetic -- no floating point is involved anywhere below.
# Pauli matrices have entries in {0,+-1,+-i}; products and short sums stay far inside int64.
import numpy as _np

_I2 = (_np.array([[1, 0], [0, 1]], dtype=_np.int64), _np.zeros((2, 2), dtype=_np.int64))
_X2 = (_np.array([[0, 1], [1, 0]], dtype=_np.int64), _np.zeros((2, 2), dtype=_np.int64))
_Z2 = (_np.array([[1, 0], [0, -1]], dtype=_np.int64), _np.zeros((2, 2), dtype=_np.int64))
_Y2 = (_np.zeros((2, 2), dtype=_np.int64), _np.array([[0, -1], [1, 0]], dtype=_np.int64))


def _nkron(A, B):
    ar, ai = A; br, bi = B
    return (_np.kron(ar, br) - _np.kron(ai, bi), _np.kron(ar, bi) + _np.kron(ai, br))


def np_matrix(p, n):
    """EXACT int64 Gaussian-integer matrix for i^m W(a)."""
    m, a = p
    x, z = split(a, n)
    M = (_np.array([[1]], dtype=_np.int64), _np.array([[0]], dtype=_np.int64))
    for j in range(n):
        f = _I2 if (x[j], z[j]) == (0, 0) else (_X2 if (x[j], z[j]) == (1, 0)
                                                else (_Z2 if (x[j], z[j]) == (0, 1) else _Y2))
        M = _nkron(M, f)
    pr, pi = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}[m % 4]
    return (pr * M[0] - pi * M[1], pr * M[1] + pi * M[0])


def np_mul(A, B):
    ar, ai = A; br, bi = B
    return (ar @ br - ai @ bi, ar @ bi + ai @ br)


def np_sub(A, B):
    return (A[0] - B[0], A[1] - B[1])


def np_add(A, B):
    return (A[0] + B[0], A[1] + B[1])


def np_is_zero(A):
    return not A[0].any() and not A[1].any()


def np_frob2(A):
    return int((A[0].astype(object) ** 2).sum() + (A[1].astype(object) ** 2).sum())


def np_trace(A):
    return (int(_np.trace(A[0].astype(object))), int(_np.trace(A[1].astype(object))))


def np_overflow_safe(A, bound=2 ** 40):
    return int(_np.abs(A[0]).max()) < bound and int(_np.abs(A[1]).max()) < bound


# ------------------------------------------------------------------ EXACT MONOMIAL BACKEND
# A Pauli is a GENERALISED PERMUTATION MATRIX: W(a)|v> = i^{x.z} (-1)^{z.v} |v xor x>.
# So the full 2^n x 2^n matrix is stored EXACTLY in 2^n int8 slots -- one Z_4 exponent per column.
# This is exact integer arithmetic and it reaches n = 20+ where a dense int64 matmul cannot.
#   mono(p, n) -> (x_int, e)  with entry (row j xor x_int, col j) equal to i^{e[j]}.

def _parity_table(zint, d):
    """e_z[j] = parity of popcount(zint & j), as an int8 array of length d.  Exact."""
    out = _np.zeros(d, dtype=_np.int8)
    j = _np.arange(d, dtype=_np.int64)
    v = _np.bitwise_and(j, _np.int64(zint))
    # popcount parity by folding
    for s in (1, 2, 4, 8, 16, 32):
        v = _np.bitwise_xor(v, _np.right_shift(v, s))
    out[:] = _np.bitwise_and(v, 1).astype(_np.int8)
    return out


def mono(p, n):
    """EXACT monomial representation of i^m W(a).  Qubit 0 is the MOST significant bit, matching
       the kron order used by xz_to_matrix / np_matrix."""
    m, a = p
    x, z = split(a, n)
    xint = 0
    zint = 0
    for j in range(n):
        xint = (xint << 1) | x[j]
        zint = (zint << 1) | z[j]
    d = 1 << n
    base = (m + dot(x, z)) % 4
    e = (_np.int8(base) + 2 * _parity_table(zint, d)) % 4
    return (xint, e.astype(_np.int8))


def mono_mul(A, B, n):
    """(A B) as a monomial: perm x_A xor x_B, exponent e[j] = e_A[j xor x_B] + e_B[j]."""
    xa, ea = A
    xb, eb = B
    d = ea.shape[0]
    j = _np.arange(d, dtype=_np.int64)
    return (xa ^ xb, ((ea[_np.bitwise_xor(j, _np.int64(xb))].astype(_np.int16) + eb) % 4).astype(_np.int8))


def mono_equal(A, B):
    return A[0] == B[0] and bool(_np.array_equal(A[1], B[1]))


def mono_diff_frob2(A, B):
    """EXACT ||A - B||_F^2 for two monomial matrices.  Integer valued."""
    if A[0] != B[0]:
        return 2 * A[1].shape[0]                      # disjoint supports, each entry modulus 1
    dd = (A[1].astype(_np.int16) - B[1].astype(_np.int16)) % 4
    # |i^a - i^b|^2 = 0, 2, 4, 2  for (a-b) mod 4 = 0, 1, 2, 3
    tab = _np.array([0, 2, 4, 2], dtype=_np.int64)
    return int(tab[dd].astype(object).sum())


def mono_trace(A):
    """EXACT trace as a Gaussian integer (re, im)."""
    x, e = A
    if x != 0:
        return (0, 0)
    c = _np.bincount(e.astype(_np.int64), minlength=4)
    return (int(c[0]) - int(c[2]), int(c[1]) - int(c[3]))


# ------------------------------------------------------------------ EXACT PAULI COMBINATIONS
# Any polynomial in Pauli operators is a SPARSE combination sum_a c_a W(a) with c_a in Z[i].
# The Paulis are orthogonal:  Tr(W(a) W(b)) = 2^n delta_{ab}.  So
#     ||X||_F^2 = 2^n sum_a |c_a|^2      and      Tr(X) = 2^n c_0
# are EXACT INTEGERS with no 2^n-sized object ever built.  This reaches ANY n.

def pc_from(p, n):
    m, a = p
    return {tuple(a): {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}[m % 4]}


def pc_add(A, B):
    out = dict(A)
    for k, v in B.items():
        w = out.get(k, (0, 0))
        w = (w[0] + v[0], w[1] + v[1])
        if w == (0, 0):
            out.pop(k, None)
        else:
            out[k] = w
    return out


def pc_neg(A):
    return {k: (-v[0], -v[1]) for k, v in A.items()}


def pc_sub(A, B):
    return pc_add(A, pc_neg(B))


def pc_mul(A, B, n):
    out = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            m = phi(list(ka), list(kb), n)
            ph = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}[m]
            c = gmul_s(gmul_s(va, vb), ph)
            kc = tuple((x + y) % 2 for x, y in zip(ka, kb))
            w = out.get(kc, (0, 0))
            w = (w[0] + c[0], w[1] + c[1])
            if w == (0, 0):
                out.pop(kc, None)
            else:
                out[kc] = w
    return out


def pc_comm(A, B, n):
    return pc_sub(pc_mul(A, B, n), pc_mul(B, A, n))


def pc_is_zero(A):
    return not A


def pc_frob2_over_dim(A):
    """||X||_F^2 / 2^n, an EXACT non-negative integer."""
    return sum(v[0] * v[0] + v[1] * v[1] for v in A.values())


def pc_trace_over_dim(A):
    """Tr(X)/2^n as an exact Gaussian rational (here always a Gaussian integer)."""
    z = tuple([0] * len(next(iter(A)))) if A else None
    if z is None:
        return (0, 0)
    return A.get(z, (0, 0))


def pc_ground_trace_ratio(A, stab_grp, n):
    """Tr(Pi X)/Tr(Pi) EXACTLY, Pi the joint +1 projector of the signed stabiliser group.
       Pi = (1/|S|) sum_s s.  Tr(Pi W(a)) is non-zero only for a in S, and then it contributes
       the stabiliser's own sign.  Returns an exact Gaussian integer (re, im)."""
    tot = (0, 0)
    for k, v in A.items():
        if k in stab_grp:
            ms = stab_grp[k]
            # s * W(a) = i^{ms} W(a)W(a) = i^{ms} I ; Tr(Pi W(a))/Tr(Pi) = i^{-ms}... see below
            ph = {0: (1, 0), 1: (0, -1), 2: (-1, 0), 3: (0, 1)}[ms % 4]
            tot = gadd_s(tot, gmul_s(v, ph))
    return tot


# ------------------------------------------------------------------ FAST INTEGER-ENCODED LAYER
# Identical mathematics, with (x|z) packed into two Python ints.  Still exact integer arithmetic;
# this exists only so the exhaustive three-body censuses finish.  Validated against the list layer.

def _pc1(v):
    return bin(v).count("1")


def enc(a, n):
    x = 0; z = 0
    for j in range(n):
        x = (x << 1) | a[j]
        z = (z << 1) | a[n + j]
    return (x, z)


def dec(k, n):
    x, z = k
    return [(x >> (n - 1 - j)) & 1 for j in range(n)] + [(z >> (n - 1 - j)) & 1 for j in range(n)]


def phi_i(A, B):
    xa, za = A; xb, zb = B
    return (_pc1(xa & za) + _pc1(xb & zb) - _pc1((xa ^ xb) & (za ^ zb)) + 2 * _pc1(za & xb)) % 4


def sp_i(A, B):
    xa, za = A; xb, zb = B
    return (_pc1(xa & zb) + _pc1(za & xb)) & 1


def xr_i(A, B):
    return (A[0] ^ B[0], A[1] ^ B[1])


_PH = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}


def qc_from(k, m=0):
    return {k: _PH[m % 4]}


def qc_mul(A, B):
    out = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            ph = _PH[phi_i(ka, kb)]
            c = gmul_s(gmul_s(va, vb), ph)
            kc = (ka[0] ^ kb[0], ka[1] ^ kb[1])
            w = out.get(kc, (0, 0))
            w = (w[0] + c[0], w[1] + c[1])
            if w == (0, 0):
                out.pop(kc, None)
            else:
                out[kc] = w
    return out


def qc_add(A, B):
    out = dict(A)
    for k, v in B.items():
        w = out.get(k, (0, 0))
        w = (w[0] + v[0], w[1] + v[1])
        if w == (0, 0):
            out.pop(k, None)
        else:
            out[k] = w
    return out


def qc_comm(A, B):
    out = dict(qc_mul(A, B))
    for k, v in qc_mul(B, A).items():
        w = out.get(k, (0, 0))
        w = (w[0] - v[0], w[1] - v[1])
        if w == (0, 0):
            out.pop(k, None)
        else:
            out[k] = w
    return out


def qc_frob2_over_dim(A):
    return sum(v[0] * v[0] + v[1] * v[1] for v in A.values())


def qc_ground_trace_ratio(A, sgrp):
    """sgrp: dict (xint,zint) -> m for the SIGNED stabiliser group."""
    tot = (0, 0)
    for k, v in A.items():
        ms = sgrp.get(k)
        if ms is not None:
            ph = {0: (1, 0), 1: (0, -1), 2: (-1, 0), 3: (0, 1)}[ms % 4]
            tot = gadd_s(tot, gmul_s(v, ph))
    return tot


def signed_stabiliser_group_i(stab_keys):
    grp = {(0, 0): 0}
    gens = [(0, k) for k in stab_keys]
    changed = True
    while changed:
        changed = False
        for k, m in list(grp.items()):
            for gm, gk in gens:
                nm = (m + gm + phi_i(k, gk)) % 4
                nk = (k[0] ^ gk[0], k[1] ^ gk[1])
                if nk not in grp:
                    grp[nk] = nm; changed = True
                else:
                    assert grp[nk] == nm, "inconsistent stabiliser sign"
    return grp
