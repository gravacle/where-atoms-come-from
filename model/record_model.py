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
def commutant(basis, tol=TOL):
    """Basis of A' = {X : [X,a] = 0 for all a in A}, by nullspace of the commutator system."""
    n = basis[0].shape[0]; I = np.eye(n)
    rows = [np.kron(I, a.T) - np.kron(a, I) for a in basis]   # vec([X,a]) with row-major vec
    Mx = np.vstack(rows)
    _, s, Vh = np.linalg.svd(Mx)
    null = Vh[np.sum(s > tol * max(Mx.shape) * (s[0] if s.size else 1)):]
    return [v.reshape(n, n) for v in null]

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

def minimal_projections(herm, n, tol=TOL, seed=0):
    """Eigenprojections of a GENERIC Hermitian element of A' -- a maximal abelian
       subalgebra of A'. Generic so the split is as fine as A' allows."""
    if not herm: return [np.eye(n, dtype=complex)]
    rng = np.random.default_rng(seed)
    X = sum(rng.normal() * A for A in herm)
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
class RecordModel:
    def __init__(self, H, Ls=(), seed=0):
        self.H = np.asarray(H, dtype=complex); self.Ls = [np.asarray(L, dtype=complex) for L in Ls]
        self.n = self.H.shape[0]; self.es = eigenspaces(self.H)
        self.A = star_algebra(self.H, self.Ls)
        self.Ac = commutant(self.A)
        self.herm = hermitian_span(self.Ac)
        self.projs = minimal_projections(self.herm, self.n, seed=seed)

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

    def report(self):
        recs = self.records()
        return dict(dim=self.n, dim_algebra=len(self.A), dim_commutant=len(self.Ac),
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

    def protection(self, regions):
        """CLAUSE (v). regions = list of projector-lists defining what 'contractible' means.
           NOT first-principles: this is CARRIER DATA and the model requires it to be supplied."""
        raise NotImplementedError("clause (v) needs a locality structure -- supply it explicitly")
