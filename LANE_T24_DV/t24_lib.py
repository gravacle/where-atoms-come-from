"""T-24 shared library -- EXACT machinery for clause (v) on D(G) quantum doubles, 1x2 torus.

Everything that decides a verdict is integer or Fraction arithmetic:
  * groups as explicit multiplication tables (integers)
  * vertex operators A_v(k) and global-average A_v, plaquette projectors B_p, as
    permutations / 0-1 diagonals on the configuration basis of C[G]^(x)4
  * sector projectors P_s = Q_v0 Q_v1 D_p0 D_p1 as sparse Fraction matrices (exact)
  * the edge-local admissible algebra C_e = { M on C[G]_e : [M (x) I, H] = 0 } via the
    EXACT integer Gram matrix of the commutator constraint (no sampling -- D-21):
      for rational x,  x^T G x = ||sum_c x_c [E_c (x) I, 8H]||_F^2,  so ker_Q(G) = C_e exactly.
  * per-eigenspace traces Tr(P_E S) as exact rationals (integers when S is a unitary
    involution restricted to an invariant subspace).

THE FLIP REDUCTION (proved in the note, verified here):
  U admissible unitary supported on edge e flips record R  <=>  {U,R} = 0 (U unitary)
  and one may take U to be a self-adjoint involution S in C_e (spectral symmetrisation).
  A record anticommuting with involution S exists  <=>  Tr(P_E S) = 0 on EVERY H-eigenspace E
  (S-grading of E must split evenly; then R = off-diagonal unitary pairing, which is
  automatically trace-balanced and non-trivial).  So the single-edge clause-(v) question
  is decided by the INTEGERS Tr(P_E S) over the involutions S of C_e.
"""
from fractions import Fraction
import numpy as np

# ------------------------------------------------------------------ groups
def make_D4():
    """D_4 = <r,s>, elements g = (a,b) = r^a s^b encoded as index a + 4*b, order 8."""
    n = 8
    def enc(a, b): return a % 4 + 4 * (b % 2)
    def dec(g): return (g % 4, g // 4)
    MUL = np.zeros((n, n), dtype=np.int64)
    for g1 in range(n):
        a1, b1 = dec(g1)
        for g2 in range(n):
            a2, b2 = dec(g2)
            MUL[g1, g2] = enc(a1 + (a2 if b1 == 0 else -a2), b1 + b2)
    INV = np.zeros(n, dtype=np.int64)
    for g in range(n):
        for h in range(n):
            if MUL[g, h] == 0: INV[g] = h
    names = {enc(0,0):"e", enc(1,0):"r", enc(2,0):"r2", enc(3,0):"r3",
             enc(0,1):"s", enc(1,1):"rs", enc(2,1):"r2s", enc(3,1):"r3s"}
    # conjugacy classes
    classes = conj_classes(MUL, INV)
    # 1-dim characters: chi_{eps,del}(r^a s^b) = eps^a del^b ; 2-dim: (2,-2,0,0,0) on classes
    chars1 = {}
    for eps in (1,-1):
        for dl in (1,-1):
            v = np.array([eps**dec(g)[0] * dl**dec(g)[1] for g in range(n)], dtype=np.int64)
            chars1[(eps,dl)] = v
    chi2 = np.zeros(n, dtype=np.int64); chi2[enc(0,0)] = 2; chi2[enc(2,0)] = -2
    irreps = [("triv", chars1[(1,1)], 1), ("chi_a(s->-1)", chars1[(1,-1)], 1),
              ("chi_b(r->-1)", chars1[(-1,1)], 1), ("chi_c(both->-1)", chars1[(-1,-1)], 1),
              ("2dim", chi2, 2)]
    return dict(n=n, MUL=MUL, INV=INV, names=names, classes=classes, irreps=irreps,
                center=[g for g in range(n) if all(MUL[g,h]==MUL[h,g] for h in range(n))],
                label="D_4")

def make_Z2():
    n = 2
    MUL = np.array([[0,1],[1,0]], dtype=np.int64)
    INV = np.array([0,1], dtype=np.int64)
    classes = conj_classes(MUL, INV)
    irreps = [("triv", np.array([1,1], dtype=np.int64), 1),
              ("sign", np.array([1,-1], dtype=np.int64), 1)]
    return dict(n=n, MUL=MUL, INV=INV, names={0:"0",1:"1"}, classes=classes,
                irreps=irreps, center=[0,1], label="Z_2")

def conj_classes(MUL, INV):
    n = MUL.shape[0]
    seen = [False]*n; classes = []
    for g in range(n):
        if seen[g]: continue
        cl = sorted({int(MUL[k, MUL[g, INV[k]]]) for k in range(n)})
        for x in cl: seen[x] = True
        classes.append(cl)
    return classes

# ------------------------------------------------------------------ 1x2 torus carrier
# configuration c = (h0, h1, u0, u1); index = ((h0*n + h1)*n + u0)*n + u1
# edges: h0: v0->v1, h1: v1->v0 (contractible arcs), u0: loop at v0, u1: loop at v1
# plaquette p0 holonomy: g_h0 g_u1 g_h0^-1 g_u0^-1 ; p1: g_h1 g_u0 g_h1^-1 g_u1^-1
class Carrier:
    def __init__(self, G):
        self.G = G; n = G["n"]; self.n = n; self.N = n**4
        idx = np.arange(self.N)
        self.u1 = idx % n; r = idx // n
        self.u0 = r % n;  r = r // n
        self.h1 = r % n;  self.h0 = r // n
        self.MUL = G["MUL"]; self.INV = G["INV"]
    def compose(self, h0, h1, u0, u1):
        return ((h0 * self.n + h1) * self.n + u0) * self.n + u1
    # vertex gauge permutations (as index maps pi with  A|c> = |pi(c)>)
    def permA0(self, k):
        M, I = self.MUL, self.INV
        return self.compose(M[k, self.h0], M[self.h1, I[k]],
                            M[M[k, self.u0], I[k]], self.u1)
    def permA1(self, k):
        M, I = self.MUL, self.INV
        return self.compose(M[self.h0, I[k]], M[k, self.h1],
                            self.u0, M[M[k, self.u1], I[k]])
    # plaquette 0-1 diagonals
    def diagB0(self):
        M, I = self.MUL, self.INV
        return (self.u0 == M[M[self.h0, self.u1], I[self.h0]]).astype(np.int64)
    def diagB1(self):
        M, I = self.MUL, self.INV
        return (self.u1 == M[M[self.h1, self.u0], I[self.h1]]).astype(np.int64)
    # single-edge index maps: replace component of edge e by group action
    def edge_comp(self, e):
        return {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}[e]
    def replace(self, e, newcomp):
        c = {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}
        c[e] = newcomp
        return self.compose(c["h0"], c["h1"], c["u0"], c["u1"])
    def edge_perm_right(self, e, z):
        """g_e -> g_e * z  on edge e (index map)."""
        return self.replace(e, self.MUL[self.edge_comp(e), z])
    def edge_perm_left(self, e, z):
        return self.replace(e, self.MUL[z, self.edge_comp(e)])

# ------------------------------------------------------------------ exact sparse (dict-of-cols)
def sp_from_perms(perms, weight):
    """sum over perms of weight * permutation-matrix, cols[j] = {i: Fraction}."""
    N = len(perms[0]); cols = [dict() for _ in range(N)]
    for pi in perms:
        for j in range(N):
            i = int(pi[j])
            cols[j][i] = cols[j].get(i, Fraction(0)) + weight
    return cols

def sp_identity(N):
    return [ {j: Fraction(1)} for j in range(N) ]

def sp_sub(colsA, colsB):
    N = len(colsA); out = [dict(colsA[j]) for j in range(N)]
    for j in range(N):
        for i, v in colsB[j].items():
            w = out[j].get(i, Fraction(0)) - v
            if w: out[j][i] = w
            else: out[j].pop(i, None)
    return out

def sp_mul(colsA, colsB):
    """(A B): col j = sum_k B[k,j] * (col k of A)."""
    N = len(colsA); out = [dict() for _ in range(N)]
    for j in range(N):
        acc = out[j]
        for k, v2 in colsB[j].items():
            for i, v1 in colsA[k].items():
                w = acc.get(i, Fraction(0)) + v1 * v2
                if w: acc[i] = w
                else: acc.pop(i, None)
    return out

def sp_mask_rows(cols, mask01, keep):
    """multiply by diag(mask) (keep=1) or diag(1-mask) (keep=0) -- rows filtered."""
    N = len(cols); out = [dict() for _ in range(N)]
    for j in range(N):
        for i, v in cols[j].items():
            if int(mask01[i]) == keep: out[j][i] = v
    return out

def sp_trace(cols):
    return sum(cols[j].get(j, Fraction(0)) for j in range(len(cols)))

def sp_trace_perm(cols, pi):
    """Tr(P * M(pi)) where M(pi)|c> = |pi(c)> :  sum_c P[c, pi(c)]."""
    t = Fraction(0)
    for c in range(len(cols)):
        t += cols[int(pi[c])].get(c, Fraction(0))
    return t

def sp_trace_diag(cols, d):
    return sum(cols[c].get(c, Fraction(0)) * int(d[c]) for c in range(len(cols)))

def sp_trace_diag_perm(cols, d, pi):
    """Tr(P * diag(d) * M(pi)) = sum_c d(pi(c)) P[c, pi(c)]."""
    t = Fraction(0)
    for c in range(len(cols)):
        p = int(pi[c])
        t += cols[p].get(c, Fraction(0)) * int(d[p])
    return t

def sp_frob_check_projector(cols):
    """Tr(P^2) computed entrywise (P symmetric is checked separately) -- equals Tr(P) iff
       the symmetric contraction is projector-consistent."""
    t = Fraction(0)
    for j in range(len(cols)):
        for i, v in cols[j].items():
            w = cols[i].get(j, Fraction(0))
            t += v * w
    return t

# ------------------------------------------------------------------ Hamiltonian as term list
def hamiltonian_terms(car):
    """8H_int-style integer term list for [.,H]=0 constraints: H = -(A0+A1+B0+B1),
       n*H_int = -(sum_k A0(k) + sum_k A1(k) + n B0 + n B1).  Sign irrelevant for commutants.
       Returns (perm_terms, diag_terms): perms weight 1 each, diags weight n each."""
    n = car.n
    perms = [car.permA0(k) for k in range(n)] + [car.permA1(k) for k in range(n)]
    diags = [(car.diagB0(), n), (car.diagB1(), n)]
    return perms, diags

def h_entries_exact(car):
    """Full H = -(A0+A1+B0+B1) as dict {(i,j): Fraction} (for record-commutator checks)."""
    n = car.n
    ent = {}
    for k in range(n):
        for pi in (car.permA0(k), car.permA1(k)):
            for j in range(car.N):
                key = (int(pi[j]), j)
                ent[key] = ent.get(key, Fraction(0)) - Fraction(1, n)
    for d, _w in [(car.diagB0(), None), (car.diagB1(), None)]:
        for j in np.nonzero(d)[0]:
            key = (int(j), int(j))
            ent[key] = ent.get(key, Fraction(0)) - 1
    return {k: v for k, v in ent.items() if v}

# ------------------------------------------------------------------ edge-local admissible algebra
def edge_admissible_algebra(car, e):
    """C_e = { M in M_n : [M (x) I_rest, H] = 0 } via exact integer Gram of constraints.
       Constraint operator for basis unit E_ab on edge e against  n*H_int (integers):
       returns (kernel_basis as list of n x n Fraction matrices, gram_rank, checks)."""
    n, N = car.n, car.N
    perms, diags = hamiltonian_terms(car)
    comp = car.edge_comp(e)
    # constraint columns: dict pos -> {colidx: int}
    pos_map = {}
    def add(pos, col, val):
        if val == 0: return
        d = pos_map.get(pos)
        if d is None: pos_map[pos] = d = {}
        d[col] = d.get(col, 0) + val
        if d[col] == 0: del d[col]
    idx = np.arange(N)
    for a in range(n):
        for b in range(n):
            col = a * n + b
            Eab_target = car.replace(e, np.full(N, a, dtype=np.int64))  # c(e->a)
            # (E (x) I) P  : for c with pi(c)_e == b : entry (pi(c)(e->a), c) += 1
            for pi in perms:
                pc = pi
                sel = (comp[pc] == b)
                rows = Eab_target[pc[sel]]
                for r, c in zip(rows, idx[sel]):
                    add((int(r), int(c)), col, 1)
                # P (E (x) I): for c with c_e == b: entry (pi(c(e->a)), c) -= 1
                sel2 = (comp == b)
                src = Eab_target[sel2]           # c(e->a)
                rows2 = pc[src]
                for r, c in zip(rows2, idx[sel2]):
                    add((int(r), int(c)), col, -1)
            for d, w in diags:
                # [E, diag]: entries (c(e->a), c) with c_e == b, val w*(d(c) - d(c(e->a)))
                sel = (comp == b)
                tgt = Eab_target[sel]
                vals = w * (d[sel] - d[tgt])
                for r, c, v in zip(tgt, idx[sel], vals):
                    if v: add((int(r), int(c)), col, int(v))
    m = n * n
    G = [[0] * m for _ in range(m)]
    for d in pos_map.values():
        items = list(d.items())
        for i1, (c1, v1) in enumerate(items):
            for c2, v2 in items[i1:]:
                G[c1][c2] += v1 * v2
                if c1 != c2: G[c2][c1] += v1 * v2
    ker = int_kernel(G, m)         # list of Fraction vectors
    basis = []
    for v in ker:
        M = [[v[a * n + b] for b in range(n)] for a in range(n)]
        basis.append(M)
    return basis, m - len(ker), pos_map

def int_kernel(G, m):
    """Exact kernel over Q of integer matrix G (m x m) by Fraction RREF."""
    A = [[Fraction(G[i][j]) for j in range(m)] for i in range(m)]
    pivots = []
    r = 0
    for c in range(m):
        p = None
        for i in range(r, m):
            if A[i][c] != 0: p = i; break
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        pivots.append(c); r += 1
        if r == m: break
    free = [c for c in range(m) if c not in pivots]
    ker = []
    for fc in free:
        v = [Fraction(0)] * m
        v[fc] = Fraction(1)
        for ri, pc in enumerate(pivots):
            v[pc] = -A[ri][fc]
        ker.append(v)
    return ker

def verify_in_commutant(car, e, M):
    """EXACT check [M (x) I, H] = 0 for an n x n Fraction matrix M on edge e."""
    n, N = car.n, car.N
    perms, diags = hamiltonian_terms(car)
    comp = car.edge_comp(e)
    acc = {}
    def add(pos, val):
        if val == 0: return
        w = acc.get(pos, Fraction(0)) + val
        if w: acc[pos] = w
        else: acc.pop(pos, None)
    idx = np.arange(N)
    targets = [car.replace(e, np.full(N, a, dtype=np.int64)) for a in range(n)]
    for a in range(n):
        for b in range(n):
            v = M[a][b]
            if v == 0: continue
            Eab_target = targets[a]
            for pi in perms:
                sel = (comp[pi] == b)
                rows = Eab_target[pi[sel]]
                for r, c in zip(rows, idx[sel]): add((int(r), int(c)), v)
                sel2 = (comp == b)
                rows2 = pi[Eab_target[sel2]]
                for r, c in zip(rows2, idx[sel2]): add((int(r), int(c)), -v)
            for d, w in diags:
                sel = (comp == b)
                tgt = Eab_target[sel]
                vals = (d[sel] - d[tgt]) * w
                for r, c, dv in zip(tgt, idx[sel], vals):
                    if dv: add((int(r), int(c)), v * int(dv))
    return len(acc) == 0, len(acc)

# ------------------------------------------------------------------ sectors
def build_sectors(car):
    """All 16 sector projectors P_s, s=(a0,a1,b0,b1), as exact sparse cols; returns dict."""
    n, N = car.n, car.N
    A0 = sp_from_perms([car.permA0(k) for k in range(n)], Fraction(1, n))
    A1 = sp_from_perms([car.permA1(k) for k in range(n)], Fraction(1, n))
    I = sp_identity(N)
    Q0 = {1: A0, 0: sp_sub(I, A0)}
    Q1 = {1: A1, 0: sp_sub(I, A1)}
    b0, b1 = car.diagB0(), car.diagB1()
    sectors = {}
    for a0 in (0, 1):
        for a1 in (0, 1):
            M = sp_mul(Q0[a0], Q1[a1])
            for k0 in (0, 1):
                Mk0 = sp_mask_rows(M, b0, k0)
                for k1 in (0, 1):
                    sectors[(a0, a1, k0, k1)] = sp_mask_rows(Mk0, b1, k1)
    return sectors

def eigen_projectors(sectors):
    """Group sectors by H-eigenvalue -(a0+a1+b0+b1): returns {k: list-of-sector-keys}."""
    eig = {}
    for s in sectors:
        eig.setdefault(sum(s), []).append(s)
    return eig
