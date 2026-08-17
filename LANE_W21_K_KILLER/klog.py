"""
klog.py -- the PHYSICAL SECTOR as it actually is.

Diagnosis that forced this file: in k2 v1 all 374 assignable algebras returned the SAME
entropy on carrier ladder3.  That is not a finding, it is a degeneracy: the Gauss operators
are X-type, they lie in the gauge-invariant algebra AND in its symplectic radical, and on that
carrier V_max(A) = V_min(A) + span(Gauss) exactly, so the 374 algebras are 374 distinct
subalgebras of B(C^128) that all RESTRICT TO THE SAME ALGEBRA on the 4-dimensional physical
sector.  The ambiguity has to be counted MODULO GAUSS or it counts nothing.

So: the gauge-invariant algebra modulo Gauss is the full Pauli algebra on c = L - V + 1
"logical" qubits acting on the physical sector, and every region algebra is a subalgebra of
that.  Physical basis states are the Gauss orbits of Z-basis states.  A gauge-invariant Pauli
(x,z) with z orthogonal to every Gauss mask acts on orbits as a monomial matrix.  Everything
below is exact; no approximation, no gauge fixing.
"""
import numpy as np
import itertools
from klib import rref_basis, span_elements, in_span, pack, unpack, sympl


def _reduce(s, ech):
    for b in ech:
        s = min(s, s ^ b)
    return s


class Sector:
    def __init__(self, G):
        self.G = G
        self.L = G.L
        self.gmasks = rref_basis([g[0] for g in G.gauss])
        self.gdim = len(self.gmasks)
        self.c = self.L - self.gdim
        self.D = 1 << self.c
        reps = sorted({_reduce(s, self.gmasks) for s in range(1 << self.L)})
        assert len(reps) == self.D, (len(reps), self.D)
        self.reps = reps
        self.idx = {r: i for i, r in enumerate(reps)}
        self.orbsize = 1 << self.gdim

    def mat(self, x, z):
        """D x D compression of the ambient Hermitian Pauli P(x,z) to the physical sector.
        Requires z orthogonal to every Gauss mask (i.e. P gauge-invariant)."""
        for g in self.gmasks:
            assert bin(g & z).count("1") % 2 == 0, "operator is not gauge invariant"
        ph = (1j) ** (bin(x & z).count("1"))
        M = np.zeros((self.D, self.D), dtype=complex)
        for o, r in enumerate(self.reps):
            tgt = self.idx[_reduce(r ^ x, self.gmasks)]
            M[tgt, o] = ph * ((-1) ** (bin(r & z).count("1")))
        return M

    def lift(self, psi_log):
        """Physical-sector vector -> ambient 2^L vector."""
        v = np.zeros(1 << self.L, dtype=complex)
        for s in range(1 << self.L):
            v[s] = psi_log[self.idx[_reduce(s, self.gmasks)]]
        return v / np.sqrt(self.orbsize)

    def hamiltonian(self, g2):
        H = np.zeros((self.D, self.D), dtype=complex)
        for cyc in self.G.cycles:
            H -= (1.0 / g2) * self.mat(0, cyc)
        for l in range(self.L):
            H -= g2 * self.mat(1 << l, 0)
        return H

    def ground_state(self, g2):
        H = self.hamiltonian(g2)
        w, V = np.linalg.eigh((H + H.conj().T) / 2)
        return V[:, 0].copy(), float(w[0]), float(w[1] - w[0])

    # ---- logical classes: GI / Gauss ----
    def logical_class(self, v, gauss_pack):
        """Canonical representative of the coset v + span(Gauss) inside F_2^{2L}."""
        return _reduce(pack(v, self.L) if isinstance(v, tuple) else v, gauss_pack)


class LogAlgebra:
    """Subalgebra of B(physical sector) spanned by a set of gauge-invariant Paulis whose
    ambient symplectic vectors form a subspace modulo Gauss.  Same block theory as klib,
    with the ambient dimension replaced by D and the ambient trace by Tr on the sector."""

    def __init__(self, basis_amb, sector, gauss_pack, name=""):
        self.sec = sector
        self.n = sector.L
        self.name = name
        # reduce modulo Gauss and re-basis
        red = rref_basis([_reduce(b, gauss_pack) for b in basis_amb])
        red = [b for b in red if b != 0]
        self.basis = red
        self.dimV = len(red)
        self.dim = 1 << self.dimV
        elems = span_elements(self.basis)
        cand = [e for e in elems
                if all(sympl(unpack(e, self.n), unpack(w, self.n)) == 0 for w in self.basis)]
        self.centre_basis = rref_basis(cand)
        self.r = len(self.centre_basis)
        assert (self.dimV - self.r) % 2 == 0
        self.k = (self.dimV - self.r) // 2
        self.nblocks = 1 << self.r
        m = sector.D / (2 ** (self.r + self.k))
        assert abs(m - round(m)) < 1e-9 and m >= 1, (sector.D, self.r, self.k)
        self.mult = int(round(m))

    def describe(self):
        return (f"{self.name}: dim_F2={self.dimV} dim_C(A)={self.dim} "
                f"blocks={self.nblocks} x M_{2**self.k} mult m={self.mult}")

    def omega(self, rho, cache):
        D = self.sec.D
        om = np.zeros((D, D), dtype=complex)
        for e in span_elements(self.basis):
            M = cache.get(e)
            if M is None:
                x, z = unpack(e, self.n)
                M = self.sec.mat(x, z)
                cache[e] = M
            om += np.trace(rho @ M) * M
        return om / D

    def entropy(self, rho, cache):
        om = self.omega(rho, cache)
        w = np.clip(np.linalg.eigvalsh((om + om.conj().T) / 2).real, 0.0, None)
        m = self.mult
        return float(-sum(nu * np.log2(nu * m) for nu in w if nu > 1e-13))
