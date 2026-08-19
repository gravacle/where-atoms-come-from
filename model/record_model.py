"""THE RECORD MODEL — constructs records from first principles, importing no values.

INPUT: a Hamiltonian H and a set of Lindblad operators {L_k}.  NOTHING ELSE.
No lattice, no gauge group, no temperature, no coupling constant, no code, no geometry.

OUTPUT: every record the pair admits, its writer, and the multi-record structure.

WHAT IS FIRST-PRINCIPLES AND WHAT IS NOT -- the model makes the boundary explicit:
  clauses (i) BIT, (ii) DURABLE, (iii) NON-TRIVIAL, (iv) WRITABLE   -- computable from (H,{L_k}) ALONE
  clause  (v) PROTECTED                                            -- REQUIRES a locality structure,
                                                                      which is CARRIER DATA, not
                                                                      derivable from (H,{L_k})
That asymmetry is a finding, not a limitation of the code: four of the five clauses are
carrier-free; the fifth is not, and any claim resting on (v) inherits a carrier.

THE CONSTRUCTION, and every step is a theorem in the register:
  C-9   clause (ii) puts R in the commutant of the *-ALGEBRA A = alg{I,H,L_k,L_k-dagger},
        because [L,R]=0 <=> [L-dagger,R]=0 for Hermitian R.
  C-10  a record non-trivial on eigenspace E exists <=> P_E A P_E is a PROPER subalgebra.
  C-11  an admissible U with U-dagger R U = -R exists <=> Tr(P_E R) = 0 on every eigenspace.
  C-12  a record satisfying (i)-(iv) exists <=> A' contains a projection that is non-trivial
        on some eigenspace of H AND trace-balanced.
  O-4   ADMISSIBLE := unitary with [U,H] = 0.
"""
import numpy as np
from itertools import product

TOL = 1e-9

# ---------------------------------------------------------------- the *-algebra
def star_algebra(H, Ls, tol=TOL):
    """Basis of A = alg{I, H, L_k, L_k-dagger}, closed under product and adjoint."""
    n = H.shape[0]
    gens = [np.eye(n, dtype=complex), H] + [np.asarray(L, dtype=complex) for L in Ls] \
                                         + [np.asarray(L, dtype=complex).conj().T for L in Ls]
    basis, M = [], np.zeros((0, n * n), dtype=complex)
    def add(X):
        nonlocal M
        v = X.reshape(-1)
        if M.shape[0] == 0:
            if np.linalg.norm(v) < tol: return False
            M = v[None, :]; basis.append(X); return True
        r0 = np.linalg.matrix_rank(M, tol=tol)
        r1 = np.linalg.matrix_rank(np.vstack([M, v[None, :]]), tol=tol)
        if r1 > r0: M = np.vstack([M, v[None, :]]); basis.append(X); return True
        return False
    for g in gens: add(g)
    frontier = list(basis)
    while frontier:
        nxt = []
        for a in frontier:
            for b in list(basis):
                for p in (a @ b, b @ a):
                    if add(p): nxt.append(p)
        frontier = nxt
        if len(basis) >= n * n: break
    return basis

# ---------------------------------------------------------------- the commutant
def _herm_blocks(H, tol=1e-8):
    """H's eigenbasis and the index ranges of its degenerate blocks."""
    w, V = np.linalg.eigh(H)
    b, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol: j += 1
        b.append((i, j + 1)); i = j + 1
    return V, b

def _project_H(X, V, blocks):
    """Projection onto the commutant of a HERMITIAN generator: in its own eigenbasis, keep
       only the diagonal blocks. O(n^2) and exact -- no solve. H is Hermitian but almost never
       unitary, so it must be handled this way; treating it as a generic generator is what sent
       the whole computation into the n^2 x n^2 fallback."""
    Y = V.conj().T @ X @ V
    Z = np.zeros_like(Y)
    for a, z in blocks: Z[a:z, a:z] = Y[a:z, a:z]
    return V @ Z @ V.conj().T

def _commutant_project(X, unitaries, hprojs, sweeps=200, tol=1e-12):
    """Project X onto the joint commutant. Unitary generators by TWISTED AVERAGING,
       X -> (X + g X g-dagger)/2; Hermitian ones by their exact block projection. The fixed
       points are exactly the operators commuting with every generator."""
    for _ in range(sweeps):
        prev = X
        for V, blocks in hprojs: X = _project_H(X, V, blocks)
        for g in unitaries: X = (X + g @ X @ g.conj().T) / 2
        if np.linalg.norm(X - prev) < tol * max(1.0, np.linalg.norm(X)): break
    return X

def commutant_element(gens, seed=0):
    """ONE generic Hermitian element of the commutant. This is all minimal_projections needs,
       and it costs a single projection.

       The basis is not required and must not be computed on the critical path: with no noise
       the commutant of the toric 2x2 has dimension 27744, so no sampling scheme can span it --
       and nothing in the model ever asks it to."""
    n = gens[0].shape[0]; I = np.eye(n)
    uni, hpr = [], []
    for g in gens:
        if np.linalg.norm(g - I) < 1e-12: continue
        if np.linalg.norm(g.conj().T @ g - I) < 1e-9: uni.append(g)
        elif np.linalg.norm(g - g.conj().T) < 1e-9:
            V, b = _herm_blocks(g); hpr.append((V, b))
        else: uni.append(None)
    if any(u is None for u in uni):
        cb = commutant(gens)
        rng = np.random.default_rng(seed)
        X = sum(rng.normal() * ((A + A.conj().T) / 2) for A in cb) if cb else np.eye(n, dtype=complex)
        return X
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return _commutant_project((A + A.conj().T) / 2, uni, hpr)

def commutant(gens, tol=TOL, samples=None, seed=0):
    """Basis of the commutant of a SET of generators.

    O-19: the commutant of a SET equals the commutant of the algebra it generates, so the
    *-algebra is never built.

    THE COST PROBLEM AND WHY THE OBVIOUS FIXES FAIL. Solving [X,g]=0 directly is n^2 x n^2.
    Restricting to H's eigenbasis first makes X block-diagonal and cuts the unknowns to
    D = sum_E m_E^2 -- but D is large PRECISELY WHEN H IS DEGENERATE, which is exactly what
    a record requires (P-1). Measured on the toric 2x2: multiplicities [4,48,152,48,4],
    D = 27744, a 12.3 GB Gram matrix. The reduction helps the case we do not care about.

    WHAT WORKS: when the generators are UNITARY -- Paulis, stabilisers, and every carrier in
    this program -- the projection onto the commutant is a twisted group average, reachable by
    iterating over generators. Sample random Hermitian matrices, project, keep what is
    independent. Falls back to the direct solve for non-unitary generators."""
    n = gens[0].shape[0]; I = np.eye(n)
    uni, hpr, rest = [], [], []
    for g in gens:
        if np.linalg.norm(g - I) < 1e-12: continue                       # identity: no constraint
        if np.linalg.norm(g.conj().T @ g - I) < 1e-9: uni.append(g)      # unitary: average
        elif np.linalg.norm(g - g.conj().T) < 1e-9:                      # Hermitian: block projection
            V, b = _herm_blocks(g); hpr.append((V, b))
        else: rest.append(g)
    if rest:                                        # neither unitary nor Hermitian: direct solve
        M = np.vstack([np.kron(I, g.T) - np.kron(g, I) for g in gens])
        G = M.conj().T @ M
        ev, U = np.linalg.eigh(G)
        return [U[:, k].reshape(n, n) for k in range(int((ev < tol * max(1.0, abs(ev).max())).sum()))]
    rng = np.random.default_rng(seed)
    if samples is None: samples = min(8 * n, 2000)
    # INCREMENTAL ORTHOGONALISATION. A full matrix_rank per sample over n^2-length vectors is
    # what made this unusable at dim 256; Gram-Schmidt against an orthonormal basis is O(dim A')
    # per sample. Stop once STALL consecutive samples add nothing -- the commutant is spanned.
    STALL = 25
    basis, onb, miss = [], [], 0
    for _ in range(samples):
        A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        X = _commutant_project((A + A.conj().T) / 2, uni, hpr)
        nx = np.linalg.norm(X)
        if nx < 1e-9: miss += 1
        else:
            v = (X / nx).reshape(-1)
            for u in onb: v = v - np.vdot(u, v) * u
            r = np.linalg.norm(v)
            if r > 1e-6:
                onb.append(v / r); basis.append(X / nx); miss = 0
            else: miss += 1
        if miss >= STALL: break
    return basis

def hermitian_span(cb, tol=TOL):
    """Real basis of the Hermitian part of the span of cb."""
    n = cb[0].shape[0]; cand = []
    for X in cb: cand += [(X + X.conj().T) / 2, 1j * (X - X.conj().T) / 2]
    out, M = [], np.zeros((0, n * n), dtype=complex)
    for X in cand:
        v = X.reshape(-1)
        r0 = np.linalg.matrix_rank(M, tol=tol) if M.shape[0] else 0
        r1 = np.linalg.matrix_rank(np.vstack([M, v[None, :]]) if M.shape[0] else v[None, :], tol=tol)
        if r1 > r0:
            M = np.vstack([M, v[None, :]]) if M.shape[0] else v[None, :]; out.append(X)
    return out

def minimal_projections_from(X, n, tol=TOL):
    """Eigenprojections of a GENERIC Hermitian element of A' -- a maximal abelian subalgebra
       of A'. Generic so the split is as fine as A' allows. Takes the element directly, so no
       basis of A' is ever needed."""
    if X is None or np.linalg.norm(X) < 1e-12: return [np.eye(n, dtype=complex)]
    w, V = np.linalg.eigh(X)
    projs, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-7: j += 1
        Q = V[:, i:j + 1]; projs.append(Q @ Q.conj().T); i = j + 1
    return projs

# ---------------------------------------------------------------- the clauses
def eigenspaces(H, tol=1e-8):
    w, V = np.linalg.eigh(H); out, i = [], 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < tol: j += 1
        Q = V[:, i:j + 1]; out.append((w[i], Q @ Q.conj().T, j - i + 1)); i = j + 1
    return out

def clause_iii(R, es, tol=TOL):
    """R is not constant on some eigenspace of H."""
    for _, P, m in es:
        M = P @ R @ P
        if np.linalg.norm(M - (np.trace(M) / m) * P) > tol: return True
    return False

def clause_iv(R, es, tol=TOL):
    """C-11 / O-4: an ADMISSIBLE flipper exists iff Tr(P_E R) = 0 on every eigenspace."""
    return all(abs(np.trace(P @ R)) < tol for _, P, _ in es)

def build_writer(R, es, tol=TOL):
    """Explicit admissible U with [U,H]=0 and U-dagger R U = -R, by a block swap."""
    n = R.shape[0]; U = np.zeros((n, n), dtype=complex)
    for _, P, m in es:
        w, V = np.linalg.eigh(P @ R @ P + (np.eye(n) - P) * 1e3)
        idx = [i for i in range(n) if w[i] < 5e2]
        sub = V[:, idx]; Rs = sub.conj().T @ R @ sub
        ws, Vs = np.linalg.eigh(Rs)
        plus = [i for i in range(len(ws)) if ws[i] > 0]; minus = [i for i in range(len(ws)) if ws[i] <= 0]
        if len(plus) != len(minus): return None
        S = np.zeros((len(ws), len(ws)), dtype=complex)
        for a, b in zip(plus, minus): S[a, b] = 1; S[b, a] = 1
        U += sub @ Vs @ S @ Vs.conj().T @ sub.conj().T
    return U

# ---------------------------------------------------------------- the model
def _vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

# ---------------------------------------------------------------- stabiliser logicals (T-20)
def symplectic_logicals(stab_xz, n):
    """COMPUTE a conjugate basis of logical operators for a stabiliser code. Never nominate them.

       Paulis are (x|z) in F_2^{2n}; P and Q anticommute iff x.z' + z.x' = 1. The logicals are
       N(S)/S, and a symplectic Gram-Schmidt gives pairs (X_i, Z_i) with {X_i,Z_i} = 0 and every
       other pair commuting.

       Nominating logicals has failed FOUR times in this program -- W-62's first build, F-7's,
       [[4,2,2]] in O-16, and [[8,3,2]] in T-20 -- and every time a self-check caught it."""
    def sp(a, b):                                   # symplectic form
        return (sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n))) % 2
    def rref(rows):
        rows = [r[:] for r in rows]; piv = []; r = 0
        for c in range(2 * n):
            p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
            if p is None: continue
            rows[r], rows[p] = rows[p], rows[r]
            for i in range(len(rows)):
                if i != r and rows[i][c]:
                    rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
            piv.append(c); r += 1
        return rows[:r], piv
    S, _ = rref(list(stab_xz))
    # N(S): the symplectic complement, as the nullspace of the form against every stabiliser
    M = [[sp(e, s) for e in [[1 if k == j else 0 for k in range(2 * n)] for j in range(2 * n)]] for s in S]
    M = [list(row) for row in zip(*[[sp([1 if k == j else 0 for k in range(2 * n)], s)
                                     for j in range(2 * n)] for s in S])]
    A = [[sp([1 if k == j else 0 for k in range(2 * n)], s) for j in range(2 * n)] for s in S]
    Ar, piv = rref(A)
    free = [c for c in range(2 * n) if c not in piv]
    N = []
    for f in free:
        v = [0] * (2 * n); v[f] = 1
        for i, c in enumerate(piv): v[c] = Ar[i][f]
        N.append(v)
    # strip anything already in S, then symplectic Gram-Schmidt
    def in_span(v, basis):
        v = v[:]
        for b in basis:
            h = next((i for i in range(2 * n) if b[i]), None)
            if h is not None and v[h]: v = [(x + y) % 2 for x, y in zip(v, b)]
        return not any(v)
    cand = [v for v in N if not in_span(v, S)]
    pairs = []
    pool = cand[:]
    while pool:
        a = pool.pop(0)
        partner = next((b for b in pool if sp(a, b) == 1), None)
        if partner is None: continue
        pool.remove(partner)
        pairs.append((a, partner))
        pool = [[(c[i] + sp(c, partner) * a[i] + sp(c, a) * partner[i]) % 2 for i in range(2 * n)]
                for c in pool]
        pool = [c for c in pool if any(c)]
    return pairs

def xz_to_matrix(v, n):
    I2 = np.eye(2); Xm = np.array([[0, 1], [1, 0]], dtype=complex); Zm = np.array([[1, 0], [0, -1]], dtype=complex)
    M = np.array([[1]], dtype=complex)
    for i in range(n):
        x, z = v[i], v[n + i]
        P = I2 if (x, z) == (0, 0) else (Xm if (x, z) == (1, 0) else (Zm if (x, z) == (0, 1) else 1j * Xm @ Zm))
        M = np.kron(M, P)
    return M

class Environment:
    """A bath of qubits at inverse temperature beta. `probe` is the bath operator the system
       couples through."""
    def __init__(self, nq=3, energies=(1.0, 1.4, 0.7), beta=2.0):
        I2 = np.eye(2)
        Xb = np.array([[0, 1], [1, 0]], dtype=complex)
        Zb = np.array([[1, 0], [0, -1]], dtype=complex)
        def bop(j, P):
            M = np.array([[1]], dtype=complex)
            for k in range(nq): M = np.kron(M, P if k == j else I2)
            return M
        self.nq, self.beta, self.dim = nq, beta, 2 ** nq
        self.HB = sum(energies[j] * bop(j, Zb) for j in range(nq))
        self.probe = sum(bop(j, Xb) for j in range(nq))
        self.site = [bop(j, Xb) for j in range(nq)]     # per-qubit probes, for DISTRIBUTED couplings

    def thermal(self):
        w, V = np.linalg.eigh(self.HB)
        p = np.exp(-self.beta * w); p /= p.sum()
        return (V * p) @ V.conj().T

    def _trace_system(self, r, nS):
        return r.reshape(nS, self.dim, nS, self.dim).trace(axis1=0, axis2=2)

    def _fragment(self, rB, keep):
        t = rB.reshape([2] * self.nq + [2] * self.nq)
        for j in reversed([j for j in range(self.nq) if j not in keep]):
            t = np.trace(t, axis1=j, axis2=j + t.ndim // 2)
        d = 2 ** len(keep)
        return t.reshape(d, d)

    def holevo(self, r, R, nS, fragment=None):
        """chi(R : bath) -- how much the bath holds about the +-1 observable R."""
        out = []
        for s in (+1, -1):
            P = np.kron((np.eye(nS) + s * R) / 2, np.eye(self.dim))
            blk = P @ r @ P; p = np.real(np.trace(blk))
            if p < 1e-12: continue
            rB = self._trace_system(blk / p, nS)
            out.append((p, self._fragment(rB, fragment) if fragment is not None else rB))
        if len(out) < 2: return 0.0
        av = sum(p * rb for p, rb in out)
        return max(_vN(av) - sum(p * _vN(rb) for p, rb in out), 0.0)

class RecordModel:
    def __init__(self, H, Ls=(), seed=0):
        self.H = np.asarray(H, dtype=complex); self.Ls = [np.asarray(L, dtype=complex) for L in Ls]
        self.n = self.H.shape[0]; self.es = eigenspaces(self.H)
        # O-19: generators only -- the *-algebra is never built.
        I = np.eye(self.n, dtype=complex)
        self.gens = [I, self.H] + self.Ls + [L.conj().T for L in self.Ls]
        self.Aelt = commutant_element(self.gens, seed=seed)
        self._A = None; self._Ac = None
        self.projs = minimal_projections_from(self.Aelt, self.n)

    def records(self):
        """Every R satisfying (i)-(iv), constructed -- not searched for."""
        out = []
        k = len(self.projs)
        if k > 20: raise RuntimeError(f"{k} minimal projections: enumeration too large")
        for signs in product((1, -1), repeat=k):
            if signs[0] == -1: continue                      # R and -R are the same record
            R = sum(s * P for s, P in zip(signs, self.projs))
            if np.linalg.norm(R @ R - np.eye(self.n)) > TOL: continue
            if not clause_iii(R, self.es): continue
            if not clause_iv(R, self.es): continue
            out.append(R)
        return out

    @property
    def Ac(self):
        """a BASIS of the commutant -- expensive, and only built if something asks (O-19)"""
        if self._Ac is None: self._Ac = commutant(self.gens)
        return self._Ac

    @property
    def herm(self):
        if self._Ac is None and self.n > 64:
            raise RuntimeError("herm needs a commutant BASIS; too large at dim %d -- use .Aelt" % self.n)
        return hermitian_span(self.Ac)

    @property
    def A(self):
        """the *-algebra, built ONLY if something asks for it (O-19)"""
        if self._A is None: self._A = star_algebra(self.H, self.Ls)
        return self._A

    def report(self, with_algebra=False):
        recs = self.records()
        return dict(dim=self.n, dim_algebra=(len(self.A) if with_algebra else None),
                    dim_commutant=len(self.Ac),
                    n_minimal_projections=len(self.projs),
                    eigenvalue_multiplicities=[m for _, _, m in self.es],
                    n_records=len(recs), records=recs)

    def signature(self, R):
        """A record's sign pattern over the minimal projections -- its F2 coordinates."""
        v = 0
        for P in self.projs:
            v = (v << 1) | (1 if np.real(np.trace(P @ R) / np.trace(P)) > 0 else 0)
        return v

    def commuting_family(self, recs, tol=1e-7):
        """MULTI-RECORD. A maximal family of records that are GENUINELY INDEPENDENT BITS.

           F2-independence of sign vectors is the WRONG notion: k independent record-bits
           must produce 2^k joint eigenspaces, so k <= log2(dim). The right criterion is
           that a new record SPLIT EVERY EXISTING JOINT EIGENSPACE EVENLY -- only then does
           it carry a bit the family did not already hold, and only then can it be flipped
           without disturbing the others."""
        fam = []
        blocks = [C for (_, PE, _) in self.es
                  for C in [np.linalg.eigh(PE)[1][:, np.linalg.eigh(PE)[0] > 0.5]] if C.shape[1]]
        for R in recs:
            if not all(np.linalg.norm(C @ C.conj().T @ R @ C @ C.conj().T
                                      - C @ C.conj().T @ R @ C @ C.conj().T) < 1 for C in blocks):
                continue
            nb, ok = [], True
            for C in blocks:
                Rs = C.conj().T @ R @ C
                if np.linalg.norm(Rs @ Rs - np.eye(Rs.shape[0])) > 1e-6: ok = False; break
                w, V = np.linalg.eigh(Rs)
                pl = [i for i in range(len(w)) if w[i] > 0]; mi = [i for i in range(len(w)) if w[i] <= 0]
                if len(pl) != len(mi) or not pl: ok = False; break      # must SPLIT, and EVENLY
                nb += [C @ V[:, pl], C @ V[:, mi]]
            if ok: fam.append(R); blocks = nb
        return fam

    def joint_basis(self, family, tol=1e-7):
        """Simultaneous eigenbasis of H and every member of the family (all commute).
           Returns {(energy_index, sign_tuple): column-block}."""
        blocks = {}
        for ei, (_, PE, m) in enumerate(self.es):
            w, V = np.linalg.eigh(PE)
            cols = V[:, w > 0.5]
            groups = {(): cols}
            for R in family:
                ng = {}
                for lab, C in groups.items():
                    Rs = C.conj().T @ R @ C
                    ws, Vs = np.linalg.eigh(Rs)
                    for s in (+1, -1):
                        idx = [i for i in range(len(ws)) if (ws[i] > 0) == (s > 0)]
                        if idx: ng[lab + (s,)] = C @ Vs[:, idx]
                groups = ng
            for lab, C in groups.items(): blocks[(ei, lab)] = C
        return blocks

    def independently_writable(self, family, tol=1e-7):
        """Can each record be flipped WITHOUT disturbing the others? Pair joint eigenspaces
           differing ONLY in the j-th sign and swap them -- a permutation, so manifestly
           unitary and automatically [U,H]=0 since it never leaves an energy shell.
           It exists iff the paired blocks have EQUAL dimension."""
        n = self.n; blocks = self.joint_basis(family); out = []
        for j in range(len(family)):
            U = np.zeros((n, n), dtype=complex); ok = True; seen = set()
            for (ei, lab), C in blocks.items():
                if (ei, lab) in seen: continue
                flip = lab[:j] + (-lab[j],) + lab[j+1:]
                D = blocks.get((ei, flip))
                if D is None or D.shape[1] != C.shape[1]: ok = False; break
                U += C @ D.conj().T + D @ C.conj().T
                seen.add((ei, lab)); seen.add((ei, flip))
            if not ok: continue
            Rj = family[j]; others = [R for i, R in enumerate(family) if i != j]
            if np.linalg.norm(U.conj().T @ U - np.eye(n)) > 1e-6: continue
            if np.linalg.norm(U @ self.H - self.H @ U) > tol: continue
            if np.linalg.norm(U.conj().T @ Rj @ U + Rj) > tol: continue
            if all(np.linalg.norm(U.conj().T @ R @ U - R) < tol for R in others): out.append(j)
        return out

    def independence(self, recs):
        fam = self.commuting_family(recs); m = len(fam)
        comm = np.zeros((m, m), dtype=bool)
        for i in range(m):
            for j in range(m):
                comm[i, j] = np.linalg.norm(fam[i] @ fam[j] - fam[j] @ fam[i]) < TOL
        return fam, comm, self.independently_writable(fam)

    # ------------------------------------------------------------ FORMATION (T-4)
    def ground_space(self, tol=1e-9):
        """projector onto H's lowest eigenspace, and its dimension"""
        w, V = np.linalg.eigh(self.H)
        k = int(np.sum(np.abs(w - w[0]) < tol))
        Q = V[:, :k]
        return Q @ Q.conj().T, k

    def formation(self, record, coupling, env, lam=0.8, t=4.0, state0=None, fragment=None):
        """chi(record : environment) after evolving a PRODUCT state under
               H_tot = H (x) I  +  I (x) H_B  +  lam * coupling (x) env.probe
           Formation is the environment coming to hold information about the record; the
           record's VALUE need not change, and under a commuting coupling it does not.
           `fragment` restricts the readout to a subset of bath qubits (redundancy, F-21)."""
        nS, nB = self.n, env.dim
        if state0 is None:
            Pg, k = self.ground_space(); state0 = Pg / k
        # A coupling may be given three ways, and the third matters: the lanes distribute each
        # system term to a SPECIFIC bath qubit rather than to the sum of them, and a product
        # ansatz A (x) probe silently gives a different number.
        if isinstance(coupling, np.ndarray) and coupling.shape[0] == nS * nB:
            HINT = coupling                                          # full interaction operator
        elif isinstance(coupling, (list, tuple)):
            HINT = sum(np.kron(A, env.site[j % env.nq]) for A, j in coupling)   # distributed
        else:
            HINT = np.kron(coupling, env.probe)                      # product, A (x) probe
        Ht = np.kron(self.H, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
        w, U = np.linalg.eigh(Ht)
        r0 = np.kron(state0, env.thermal())
        Uc = U.conj().T @ r0 @ U
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        return env.holevo(r, record, nS, fragment=fragment)

    def evolve(self, coupling, env, lam=0.8, t=4.0, state0=None):
        """The joint state after unitary evolution from a product state. Returned so that many
           readouts -- whole bath, every fragment, several records -- share ONE eigendecomposition.
           formation() redoing eigh per fragment is what made a 7-qubit bath unreachable."""
        nS, nB = self.n, env.dim
        if state0 is None:
            Pg, k = self.ground_space(); state0 = Pg / k
        if isinstance(coupling, np.ndarray) and coupling.shape[0] == nS * nB: HINT = coupling
        elif isinstance(coupling, (list, tuple)):
            HINT = sum(np.kron(A, env.site[j % env.nq]) for A, j in coupling)
        else: HINT = np.kron(coupling, env.probe)
        Ht = np.kron(self.H, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
        w, U = np.linalg.eigh(Ht)
        Uc = U.conj().T @ np.kron(state0, env.thermal()) @ U
        ph = np.exp(-1j * w * t)
        return U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T

    def redundancy(self, record, coupling, env, lam=0.8, t=4.0):
        """chi held by the whole bath and by EACH single-qubit fragment, from one evolution."""
        r = self.evolve(coupling, env, lam=lam, t=t)
        whole = env.holevo(r, record, self.n)
        parts = [env.holevo(r, record, self.n, fragment=[j]) for j in range(env.nq)]
        return whole, np.array(parts)

    def channel(self, record, coupling, tol=1e-9):
        """DOES THIS COUPLING OPEN A CHANNEL TO THIS RECORD?

           G-16 as first registered said: iff the coupling fails to commute with the record's
           CONJUGATE. THAT IS NECESSARY AND NOT SUFFICIENT, and it was validated only on
           weight-2 couplings where the distinction cannot arise. A coupling can anticommute
           with the writer and still be a DIFFERENT logical -- Zbar*Zbar2 anticommutes with
           Xbar but is not Zbar -- and the bath then learns that other logical, which in a
           maximally mixed code space says nothing about this record.

           THE CRITERION IS WHAT THE COUPLING DOES ON THE CODE SPACE: it opens a channel iff
           its compression has a non-zero component along the record itself."""
        Pg, k = self.ground_space()
        M = Pg @ coupling @ Pg
        Rc = Pg @ record @ Pg
        nrm = np.real(np.trace(Rc.conj().T @ Rc))
        if nrm < tol: return dict(component=0.0, opens_channel=False, identity_part=0.0, residual=0.0)
        a = np.trace(Rc.conj().T @ M) / nrm
        b = np.trace(M) / k
        resid = float(np.linalg.norm(M - a * Rc - b * Pg))
        return dict(component=float(abs(a)), opens_channel=bool(abs(a) > tol),
                    identity_part=float(abs(b)), residual=resid)

    # ------------------------------------------------------------ MULTI-RECORD (T-6)
    def channel_map(self, family, couplings, tol=1e-9):
        """WHICH COUPLING OPENS A CHANNEL TO WHICH RECORD. Rows are couplings, columns are
           records. This is the dependency structure of a multi-record environment: a coupling
           serves exactly the records its compression has a component along."""
        return np.array([[self.channel(R, A, tol=tol)['opens_channel'] for R in family]
                         for A in couplings], dtype=bool)

    def formation_independence(self, family, couplings, env, lam=0.8, t=4.0, tol=1e-8):
        """CAN EACH RECORD BE FORMED WITHOUT DISTURBING THE OTHERS? Measured, not asserted.

           For each (coupling, target) pair that opens a channel, evolve and report for EVERY
           record in the family: how much the bath learned about it, and how far its value
           moved. A record is formed INDEPENDENTLY when the bath learns about it and about no
           other member, and no other member's value moves."""
        Pg, k = self.ground_space(); rho0 = Pg / k
        out = []
        for ci, A in enumerate(couplings):
            targets = [j for j, R in enumerate(family) if self.channel(R, A)['opens_channel']]
            if not targets: continue
            nS, nB = self.n, env.dim
            if isinstance(A, (list, tuple)):
                HINT = sum(np.kron(a, env.site[j % env.nq]) for a, j in A)
            else:
                HINT = np.kron(A, env.probe)
            Ht = np.kron(self.H, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * HINT
            w, U = np.linalg.eigh(Ht)
            r0 = np.kron(rho0, env.thermal())
            ph = np.exp(-1j * w * t)
            Uc = U.conj().T @ r0 @ U
            r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
            rS = r.reshape(nS, nB, nS, nB).trace(axis1=1, axis2=3)
            learned = [env.holevo(r, R, nS) for R in family]
            moved = [abs(float(np.real(np.trace(rS @ R)) - np.real(np.trace(rho0 @ R)))) for R in family]
            out.append(dict(coupling=ci, targets=targets, learned=learned, moved=moved,
                            independent=bool(len(targets) == 1
                                             and all(learned[j] < tol for j in range(len(family)) if j not in targets)
                                             and all(m < tol for j, m in enumerate(moved) if j not in targets))))
        return out

    def protection(self, regions):
        """CLAUSE (v). regions = list of projector-lists defining what 'contractible' means.
           NOT first-principles: this is CARRIER DATA and the model requires it to be supplied."""
        raise NotImplementedError("clause (v) needs a locality structure -- supply it explicitly")
