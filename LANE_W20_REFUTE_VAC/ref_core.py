# ref_core.py -- LANE W20_REFUTE_VAC.  AN INDEPENDENT INSTRUMENT.  Prints nothing.
#
# THIS FILE SHARES NO CODE WITH LANE_W20_C_CHARGE OR LANE_W20_R_LEDGER.
# It deliberately shares no ALGORITHM either, at the two places where a shared bug would be
# invisible:
#
#   (1) THE PHYSICAL SECTOR.  Both lanes build it combinatorially, by enumerating the 128-element
#       cocycle group B^1 and forming cosets (LANE_R: PhysSector.orb/sgn; LANE_C: Sector.U over
#       CYC).  If that enumeration were wrong in the same way twice, nothing in either lane would
#       notice.  Here the sector is built by BRUTE LINEAR ALGEBRA: the 4096x4096 matrices G_v are
#       assembled by a 12-fold Kronecker product of 2x2 Paulis, the projector
#       P = prod_v (I + eta_v G_v)/2 is applied to a random 4096x64 block, and an orthonormal
#       basis of its range is taken by SVD.  No bit arithmetic, no group enumeration, no GF(2).
#       Call this the GOLD route.  A faster coset route (FAST) is provided for bulk work and is
#       checked against GOLD before it is used for anything.
#
#   (2) THE ALGEBRAIC ENTROPY.  Both lanes compute it from a GF(2) SYMPLECTIC DECOMPOSITION:
#       find a radical and hyperbolic pairs over GF(2), lift them to Pauli strings with an
#       i^{|a&c|} phase, reconstruct a block density matrix from expectation values, diagonalise.
#       That route has exactly one place where a missing complex conjugation deletes the whole
#       "Y" Bloch direction -- LANE_C found that bug live in its own instrument, and a shared
#       convention error there would be invisible to a cross-check that used the same route.
#       Here the entropy is computed with NO GF(2) and NO PHASE CONVENTION AT ALL:
#
#         * the algebra A is generated as a LINEAR SPAN of explicit 32x32 matrices,
#         * its commutant A' is the nullspace of the explicit commutator map,
#         * its centre Z = A n A' is an explicit subspace intersection,
#         * the minimal central projections come from the eigenspaces of a generic Hermitian
#           element of Z,
#         * sigma = the Hilbert-Schmidt projection of rho = |psi><psi| onto A,
#         * and then the IDENTITY
#
#               S(rho|A)  =  S_vN(sigma)  -  sum_k p_k log2 m_k
#
#           where m_k is the multiplicity of block k (A|_k = B(C^{a_k}) (x) 1_{m_k}).
#
#       PROOF OF THE IDENTITY, since it is the load-bearing part of this instrument.
#       Write H = (+)_k H_k^A (x) H_k^B with A = (+)_k B(H_k^A) (x) 1_{m_k}, m_k = dim H_k^B.
#       The HS projection onto B(H^A) (x) 1 inside a block is X |-> Tr_B(X) (x) 1/m, because
#       <Y(x)1, X> = Tr(Y^dag Tr_B X) = <Y(x)1, Tr_B(X)(x)1/m>.  Off-block parts are killed.
#       So sigma|_k = (p_k rho_k) (x) 1_{m_k}/m_k, whose eigenvalues are p_k lam_j / m_k with
#       multiplicity m_k.  Then
#         sum_{eigs in k} -nu log2 nu = -p_k log2(p_k/m_k) + p_k S(rho_k),
#       and summing over k,
#         S_vN(sigma) = sum_k [-p_k log2 p_k] + sum_k p_k log2 m_k + sum_k p_k S(rho_k)
#                     = S(rho|A) + sum_k p_k log2 m_k.                                    QED
#       Sanity anchors, all asserted in ref_verify.py: A = B(H) gives 0 on a pure state;
#       A = C.I gives 0; A maximal abelian gives the Shannon entropy of the diagonal.
#
# CONVENTIONS (fixed here, matched to the program's, and checked against both lanes numerically):
#   Z_2 on each of L = 12 links of tri_chain12 (V = 8, cubic, girth 3).
#   G_v = prod_{l incident to v} X_l ;  physical sector G_v = eta_v = (-1)^{q_v}.
#   H = -(1/g2) sum_{p in PLAQ} W_p  -  g2 sum_l X_l ,  W_p = prod_{l in p} Z_l.
#   PLAQ = [1,2,3],[3,4,5],[7,8,9],[9,10,11],[0,1,4,6,7,10].
#   Basis index j in 0..4095, bit l of j is link l, link 0 = LEAST significant.

import numpy as np
import itertools, math

# ------------------------------------------------------------------ carrier
V_N = 8
EDGES = [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
         (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
L = len(EDGES)
DIM_FREE = 1 << L

PLAQ_LINKS = [[1, 2, 3], [3, 4, 5], [7, 8, 9], [9, 10, 11], [0, 1, 4, 6, 7, 10]]
S_LINKS = [1, 2, 3]
SIGMA_LINKS = [0, 4, 5]
ENV_LINKS = [0, 4, 5, 6, 7, 8, 9, 10, 11]
FRAG_P = [("F1", [0, 4, 5]), ("F2", [7, 8, 9]), ("F3", [11]), ("F4", [6, 10])]
FRAG_S = [("G1", [0, 4, 5]), ("G2", [6, 9]), ("G3", [7, 11]), ("G4", [8, 10])]
GRID = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20, 3.00, 5.00]

STAR = [[l for l, (a, b) in enumerate(EDGES) if a == v or b == v] for v in range(V_N)]
STAR_MASK = [sum(1 << l for l in s) for s in STAR]


def mask_of(links):
    return sum(1 << l for l in links)


def links_of(m):
    return [l for l in range(L) if (m >> l) & 1]


def popc(m):
    return bin(m).count("1")


def boundary(mask):
    """d(mask): the vertex set (as a bitmask over V) that the link set has odd degree at."""
    r = 0
    for l in links_of(mask):
        a, b = EDGES[l]
        r ^= (1 << a) ^ (1 << b)
    return r


# ------------------------------------------------------------------ GOLD: brute Kronecker route
_I2 = np.eye(2)
_X2 = np.array([[0.0, 1.0], [1.0, 0.0]])
_Z2 = np.array([[1.0, 0.0], [0.0, -1.0]])


def kron_full(ops):
    """12-fold Kronecker product; ops[l] is the 2x2 factor on link l; link 0 least significant."""
    M = np.array([[1.0]])
    for l in range(L - 1, -1, -1):
        M = np.kron(M, ops[l])
    return M


def _pauli_string_matrix(a_mask, c_mask):
    """Returns the two dense factors (X^a, Z^c) separately, to avoid a 4096^3 matmul."""
    ox = [_I2] * L
    for l in links_of(a_mask):
        ox[l] = _X2
    oz = [_I2] * L
    for l in links_of(c_mask):
        oz[l] = _Z2
    return ox, oz


def gold_sector_basis(eta, seed=20260817):
    """Orthonormal 4096 x 32 basis Q of the physical sector, by explicit projector.
    eta is a length-8 list of +-1 with product +1."""
    assert len(eta) == V_N and int(np.prod(eta)) == 1
    rng = np.random.default_rng(seed)
    R = rng.normal(size=(DIM_FREE, 64))
    for v in range(V_N):
        ops = [_I2] * L
        for l in STAR[v]:
            ops[l] = _X2
        G = kron_full(ops)
        R = 0.5 * (R + eta[v] * (G @ R))
        del G
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-8 * s[0]
    Q = U[:, keep]
    return Q, s


def gold_restrict(Q, a_mask, c_mask):
    """32x32 matrix of X^a Z^c restricted to the sector spanned by Q (real, no phase)."""
    ox, oz = _pauli_string_matrix(a_mask, c_mask)
    W = Q
    if c_mask:
        MZ = kron_full(oz)
        W = MZ @ W
        del MZ
    if a_mask:
        MX = kron_full(ox)
        W = MX @ W
        del MX
    return Q.T @ W


def gold_H_parts(Q):
    """(mag32, elec32) with H(g2) = -(1/g2) mag32 - g2 elec32."""
    mag = np.zeros((Q.shape[1], Q.shape[1]))
    for p in PLAQ_LINKS:
        mag += gold_restrict(Q, 0, mask_of(p))
    ele = np.zeros((Q.shape[1], Q.shape[1]))
    for l in range(L):
        ele += gold_restrict(Q, 1 << l, 0)
    return mag, ele


# ------------------------------------------------------------------ FAST: coset route (mine)
class Sector:
    """32-dim physical sector for a charge set, by canonical coset representatives.

    Representation: basis vector k is  |k> = 128^{-1/2} sum_{m in B} s(m) |r_k xor m>,
    B = {d(A) : A subset V} (128 masks), s(d(A)) = prod_{v in A} eta_v.
    In this basis Z^c is DIAGONAL for a cycle c (since <c, d(A)> = 0), with entry
    (-1)^{|c and r_k|}; and X^a is a signed permutation.
    """

    _bgrp = None

    def __init__(self, charges):
        self.charges = tuple(sorted(charges))
        self.qmask = sum(1 << v for v in self.charges)
        assert popc(self.qmask) % 2 == 0, "odd total charge is inadmissible"
        self.eta = [1 - 2 * ((self.qmask >> v) & 1) for v in range(V_N)]
        if Sector._bgrp is None:
            g = {}
            for A in range(1 << V_N):
                m = 0
                for v in range(V_N):
                    if (A >> v) & 1:
                        m ^= STAR_MASK[v]
                g.setdefault(m, []).append(A)
            assert len(g) == 128
            Sector._bgrp = g
        self.sign = {}
        for m, As in Sector._bgrp.items():
            vals = set()
            for A in As:
                s = 1
                for v in range(V_N):
                    if (A >> v) & 1:
                        s *= self.eta[v]
                vals.add(s)
            assert len(vals) == 1, "eta inconsistent"
            self.sign[m] = vals.pop()
        # canonical representatives
        seen = {}
        reps = []
        for c in range(DIM_FREE):
            if c in seen:
                continue
            orb = [c ^ m for m in self.sign]
            r = min(orb)
            k = len(reps)
            reps.append(r)
            for x in orb:
                seen[x] = k
        self.reps = reps
        self.index = seen
        self.dim = len(reps)
        assert self.dim == 32
        self._xc = {}
        self._zc = {}

    def xmat(self, a_mask):
        """32x32 signed permutation matrix of X^a."""
        if a_mask not in self._xc:
            M = np.zeros((32, 32))
            for k, r in enumerate(self.reps):
                t = r ^ a_mask
                k2 = self.index[t]
                m = t ^ self.reps[k2]
                M[k2, k] = self.sign[m]
            self._xc[a_mask] = M
        return self._xc[a_mask]

    def zdiag(self, c_mask):
        """diagonal of Z^c (c must be a cycle)."""
        if c_mask not in self._zc:
            assert boundary(c_mask) == 0, "Z^c leaves the sector unless c is a cycle"
            self._zc[c_mask] = np.array(
                [1.0 - 2.0 * (popc(c_mask & r) & 1) for r in self.reps])
        return self._zc[c_mask]

    def zmat(self, c_mask):
        return np.diag(self.zdiag(c_mask))

    def op(self, a_mask, c_mask):
        """X^a Z^c as a 32x32 real matrix (no i^{|a&c|} phase: we never need one)."""
        return self.xmat(a_mask) * self.zdiag(c_mask)[None, :]

    def H_parts(self):
        mag = np.zeros((32, 32))
        for p in PLAQ_LINKS:
            mag += self.zmat(mask_of(p))
        ele = np.zeros((32, 32))
        for l in range(L):
            ele += self.xmat(1 << l)
        return mag, ele

    def H(self, g2, mag=None, ele=None, elec_flip=0):
        """elec_flip: a link mask whose electric term carries the opposite sign (used to test
        the claim that a charge pair IS a sign-flipped electric coupling along a string)."""
        if mag is None or ele is None:
            mag, ele = self.H_parts()
        if elec_flip:
            ele = np.zeros((32, 32))
            for l in range(L):
                s = -1.0 if (elec_flip >> l) & 1 else 1.0
                ele = ele + s * self.xmat(1 << l)
        return -(1.0 / g2) * mag - g2 * ele

    def ground(self, g2, **kw):
        Hm = self.H(g2, **kw)
        w, v = np.linalg.eigh(Hm)
        return w, v[:, 0].astype(complex)


# ------------------------------------------------------------------ matrix algebra machinery
def _orth(cols, tol=1e-9):
    """orthonormal basis of the column space of `cols` (n x m)."""
    if cols.shape[1] == 0:
        return cols
    U, s, _ = np.linalg.svd(cols, full_matrices=False)
    k = int((s > tol * max(1.0, s[0])).sum())
    return U[:, :k]


def _vec(M):
    return M.reshape(-1)


class MAlg:
    """A *-subalgebra of B(C^32) given by explicit Hermitian generator MATRICES.
    Nothing about GF(2), cocycles, Pauli phases or the graph enters here."""

    def __init__(self, gens, label=""):
        self.label = label
        d = gens[0].shape[0] if gens else 32
        self.d = d
        base = [np.eye(d, dtype=complex)]
        cols = np.array([_vec(base[0])]).T
        cols = _orth(cols)
        frontier = [np.eye(d, dtype=complex)]
        allm = [np.eye(d, dtype=complex)]
        for _ in range(64):
            newm = []
            for M in frontier:
                for g in gens:
                    P = M @ g
                    v = _vec(P)[:, None]
                    test = np.hstack([cols, v])
                    if _orth(test).shape[1] > cols.shape[1]:
                        cols = _orth(test)
                        newm.append(P)
                        allm.append(P)
            if not newm:
                break
            frontier = newm
        self.dim = cols.shape[1]
        # HS-orthonormal basis of A: columns of `cols` are already orthonormal in the vec inner
        # product <u,v> = sum conj(u) v = Tr(U^dag V).  Reshape back to matrices.
        self.basis = [cols[:, i].reshape(d, d) for i in range(self.dim)]
        self.Bstack = np.array(self.basis)                       # (dim, d, d)
        self._commutant = None
        self._centre = None
        self._blocks = None

    # ---- commutant, centre, blocks
    def commutant_basis(self, gens=None):
        if self._commutant is None:
            d = self.d
            rows = []
            for B in self.basis[1:]:
                Lm = np.kron(np.eye(d), B) - np.kron(B.T, np.eye(d))
                rows.append(Lm)
            Amat = np.vstack(rows) if rows else np.zeros((1, d * d))
            U, s, Vh = np.linalg.svd(Amat)
            tol = 1e-9 * max(1.0, s[0] if s.size else 1.0)
            ns = Vh[(np.concatenate([s, np.zeros(Vh.shape[0] - s.size)]) <= tol)].conj().T
            self._commutant = [ns[:, i].reshape(d, d) for i in range(ns.shape[1])]
        return self._commutant

    def centre_basis(self):
        if self._centre is None:
            Ba = np.array([_vec(M) for M in self.basis]).T          # d^2 x dimA
            Bc = np.array([_vec(M) for M in self.commutant_basis()]).T
            Mm = np.hstack([Ba, -Bc])
            U, s, Vh = np.linalg.svd(Mm)
            tol = 1e-9 * max(1.0, s[0] if s.size else 1.0)
            sfull = np.concatenate([s, np.zeros(Vh.shape[0] - s.size)])
            coeff = Vh[sfull <= tol].conj().T
            Zc = Ba @ coeff[:Ba.shape[1], :]
            Zo = _orth(Zc)
            self._centre = [Zo[:, i].reshape(self.d, self.d) for i in range(Zo.shape[1])]
        return self._centre

    def blocks(self, seed=99):
        """minimal central projections; returns list of (P_k, a_k, m_k)."""
        if self._blocks is None:
            Zb = self.centre_basis()
            rng = np.random.default_rng(seed)
            Zh = np.zeros((self.d, self.d), dtype=complex)
            for M in Zb:
                Zh = Zh + rng.normal() * (M + M.conj().T) / 2 + 1j * rng.normal() * (
                    M - M.conj().T) / 2
            Zh = (Zh + Zh.conj().T) / 2
            w, U = np.linalg.eigh(Zh)
            groups = []
            cur = [0]
            for i in range(1, len(w)):
                if abs(w[i] - w[i - 1]) > 1e-7 * max(1.0, abs(w).max()):
                    groups.append(cur)
                    cur = []
                cur.append(i)
            groups.append(cur)
            out = []
            for gidx in groups:
                Uk = U[:, gidx]
                Pk = Uk @ Uk.conj().T
                # dim of A restricted to block k
                proj = np.array([_vec(Pk @ B @ Pk) for B in self.basis]).T
                ak2 = _orth(proj).shape[1]
                ak = int(round(math.sqrt(ak2)))
                assert abs(ak * ak - ak2) < 1e-9, (ak2, self.label)
                rk = len(gidx)
                assert rk % ak == 0, (rk, ak, self.label)
                out.append((Pk, ak, rk // ak))
            assert sum(ak * ak for (_, ak, _) in out) == self.dim, (
                [x[1] for x in out], self.dim, self.label)
            self._blocks = out
        return self._blocks

    # ---- the entropy
    def entropy(self, psi):
        """S(rho|A) for rho = |psi><psi|, psi a normalised complex vector."""
        # sigma = HS projection of rho onto A.  NOTE THE CONJUGATION: the coefficient is
        # Tr(B^dag rho) = <psi|B^dag|psi> = conj(<psi|B|psi>).  Getting this wrong is exactly the
        # bug LANE_C found in its own instrument, so it is written out rather than inlined.
        amp = np.einsum('i,kij,j->k', psi.conj(), self.Bstack, psi)     # <psi|B_k|psi>
        coef = np.conj(amp)                                            # Tr(B_k^dag rho)
        sigma = np.einsum('k,kij->ij', coef, self.Bstack)
        sigma = (sigma + sigma.conj().T) / 2
        ev = np.linalg.eigvalsh(sigma).real
        S = 0.0
        for x in ev:
            if x > 1e-13:
                S -= x * math.log2(x)
            elif x < -1e-8:
                raise AssertionError("negative sigma eigenvalue %r in %s" % (x, self.label))
        corr = 0.0
        for (Pk, ak, mk) in self.blocks():
            pk = float(np.real(psi.conj() @ (Pk @ psi)))
            if pk > 1e-14:
                corr += pk * math.log2(mk)
        return S - corr

    def block_state(self, psi):
        """(p_k, rho_k) per block, for trace-distance work.  rho_k is a_k x a_k."""
        amp = np.einsum('i,kij,j->k', psi.conj(), self.Bstack, psi)
        coef = np.conj(amp)
        sigma = np.einsum('k,kij->ij', coef, self.Bstack)
        sigma = (sigma + sigma.conj().T) / 2
        out = []
        for (Pk, ak, mk) in self.blocks():
            pk = float(np.real(psi.conj() @ (Pk @ psi)))
            w, U = np.linalg.eigh(Pk)
            cols = U[:, w > 0.5]
            sk = cols.conj().T @ sigma @ cols                    # (ak*mk) x (ak*mk)
            ev = np.linalg.eigvalsh(sk).real
            ev = np.sort(ev)[::-1]
            lam = ev[::mk][:ak] * mk / pk if pk > 1e-14 else np.zeros(ak)
            out.append((pk, np.clip(lam, 0.0, None)))
        return out

    def maxent(self):
        return sum(math.log2(ak) for (_, ak, _) in self.blocks()) / len(self.blocks()) + \
               math.log2(len(self.blocks()))


_ALG_CACHE = {}


def alg_from_links(sec, links, label=""):
    """the gauge-invariant algebra supported on a link set: X_l for l in the set, plus W_z for
    every cycle z inside it."""
    key = (sec.charges, tuple(sorted(links)), "links")
    if key in _ALG_CACHE:
        return _ALG_CACHE[key]
    m = mask_of(links)
    gens = [sec.op(1 << l, 0).astype(complex) for l in links]
    for z in cycles_inside(m):
        gens.append(sec.op(0, z).astype(complex))
    A = MAlg(gens, label or str(sorted(links)))
    _ALG_CACHE[key] = A
    return A


def alg_from_ops(sec, ops, label=""):
    key = (sec.charges, tuple(ops), "ops")
    if key in _ALG_CACHE:
        return _ALG_CACHE[key]
    A = MAlg([sec.op(a, c).astype(complex) for (a, c) in ops], label)
    _ALG_CACHE[key] = A
    return A


_CYC = None


def all_cycles():
    global _CYC
    if _CYC is None:
        _CYC = [m for m in range(DIM_FREE) if boundary(m) == 0]
    return _CYC


def cycles_inside(m):
    """a GF(2) basis of the cycles supported inside link mask m (independent generators only)."""
    out = []
    basis = []
    for z in all_cycles():
        if z == 0 or (z & ~m):
            continue
        x = z
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis = sorted(basis + [x], reverse=True)
            out.append(z)
    return out


def joinA(sec, A, B, label=""):
    gens = [x for x in A.basis[1:]] + [x for x in B.basis[1:]]
    gens = [(g + g.conj().T) / 2 for g in gens] + [1j * (g - g.conj().T) / 2 for g in gens]
    return MAlg(gens, label or ("(%s v %s)" % (A.label, B.label)))


_JOIN_CACHE = {}


def MI(sec, A, B, psi):
    key = (id(A), id(B))
    if key not in _JOIN_CACHE:
        _JOIN_CACHE[key] = joinA(sec, A, B)
    AB = _JOIN_CACHE[key]
    return A.entropy(psi) + B.entropy(psi) - AB.entropy(psi)


# ------------------------------------------------------------------ named algebras
def named_algebras(sec):
    S = mask_of(S_LINKS)
    A_S = alg_from_ops(sec, [(1 << 1, 0), (1 << 2, 0), (1 << 3, 0), (0, S)], "A_S")
    A_CEN = alg_from_ops(sec, [((1 << 1) | (1 << 2), 0), ((1 << 1) | (1 << 3), 0)], "Z(A_S)")
    A_BLK = alg_from_ops(sec, [(1 << 1, 0), (0, S)], "BLOCK")
    A_MAG = alg_from_ops(sec, [(0, S)], "MAG")
    A_SIG = alg_from_links(sec, SIGMA_LINKS, "A_Sigma")
    A_ENV = alg_from_links(sec, ENV_LINKS, "A_env")
    return dict(S=A_S, CEN=A_CEN, BLK=A_BLK, MAG=A_MAG, SIG=A_SIG, ENV=A_ENV)


# ------------------------------------------------------------------ random states
def haar_state(seed, d=32):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=d) + 1j * rng.normal(size=d)
    return v / np.linalg.norm(v)


def random_unitary_in(A, seed):
    """a random unitary lying in the *-algebra A (used to twirl without disturbing A's commutant)."""
    rng = np.random.default_rng(seed)
    Hh = np.zeros((A.d, A.d), dtype=complex)
    for B in A.basis:
        Hh = Hh + (rng.normal() + 1j * rng.normal()) * B
    Hh = (Hh + Hh.conj().T) / 2
    w, U = np.linalg.eigh(Hh)
    return U @ np.diag(np.exp(1j * w * 3.0)) @ U.conj().T


def trace_distance(p, q):
    """between two lists of (p_k, lam_k) block states with matching block count."""
    tot = 0.0
    for (pa, la), (pb, lb) in zip(p, q):
        a = np.sort(pa * la)[::-1]
        b = np.sort(pb * lb)[::-1]
        n = max(len(a), len(b))
        a = np.pad(a, (0, n - len(a)))
        b = np.pad(b, (0, n - len(b)))
        tot += 0.5 * np.abs(a - b).sum()
    return tot
