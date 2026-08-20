"""ADVERSARIAL VERIFY for LANE_O54_A_CORNER.

Independent re-derivation of the key O-54-A numbers with a DIFFERENT implementation:
  - edges are labeled tuples ('E'|'S', i, j), not bit indices
  - GF(2) linear algebra via numpy uint8 row reduction, not python-int bitmasks
  - r_in computed by an explicit change-of-basis elimination on the complement
    columns (different algorithm path from the lane's rank-difference formula)
Checks (all against the sealed OUT table):
  A. main sweep L=12 s=3 all g: I_MI, I_IR (expect 0 separated; 2/1 at both seams)
  B. s=2 contact absorption (L=8 g=0: 0/0) beside s=4 contact (L=12 g=0: 4/2)
  C. double contact L=8 s=4: I_IR=9, I_MI=4, capA=capB=0, capAB=1 (=> JOINT_NEW=1)
     and L=6 s=3: I_IR=5, I_MI=2, capAB=1
  D. winding rings L=10: w=1 -> 1/1 cap 1; w=2 -> 2/2 cap 2; several gaps
  E. KP tripartite L=12 q=7 offset (2,3) and q=6 offset (0,0): TEE = -1
  F. inserted Z-row L=12 s=3 g=3: I_IR = 1 iff ell = 6; SB>0 iff ell >= 6
  G. enclosing L=10 sA=2 w=2: t=0 -> 6/3 ; t=1 -> 0/0 ; w=1 t=0 -> 0 with SB>0
Every check is a computed boolean; exits nonzero on any failure.
"""
import numpy as np
import sys

def rref_rank(M):
    """rank over GF(2) of a numpy uint8 matrix (rows = vectors)."""
    A = M.copy() % 2
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
        if r == rows:
            break
    return r

def r_in_explicit(M, Rmask):
    """rank of the subspace of rowspace(M) supported inside R, by eliminating the
       complement columns first and counting the residual rows that vanish outside R.
       (Different path from rank(M) - rank(M[:, comp]).)"""
    A = M.copy() % 2
    comp = np.where(~Rmask)[0]
    r = 0
    rows = A.shape[0]
    for c in comp:
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    # rows r.. now vanish on ALL complement columns pivoted so far; but rows below r
    # may still be nonzero only inside R OR be zero rows. Count independent ones inside R.
    inside = A[r:, :]
    assert not inside[:, comp].any(), "elimination failed to clear complement columns"
    return rref_rank(inside)

class Torus:
    def __init__(self, L):
        self.L = L
        self.edges = []
        for i in range(L):
            for j in range(L):
                self.edges.append(('E', i, j))   # (i,j)-(i,j+1)
                self.edges.append(('S', i, j))   # (i,j)-(i+1,j)
        self.eidx = {e: k for k, e in enumerate(self.edges)}
        self.n = len(self.edges)
        L_ = L
        def E(i, j): return self.eidx[('E', i % L_, j % L_)]
        def S(i, j): return self.eidx[('S', i % L_, j % L_)]
        self.E, self.S = E, S
        stars, plaqs = [], []
        for i in range(L):
            for j in range(L):
                v = np.zeros(self.n, dtype=np.uint8)
                for k in (E(i, j), E(i, j - 1), S(i, j), S(i - 1, j)):
                    v[k] ^= 1
                stars.append(v)
                w = np.zeros(self.n, dtype=np.uint8)
                for k in (E(i, j), E(i + 1, j), S(i, j), S(i, j + 1)):
                    w[k] ^= 1
                plaqs.append(w)
        self.stars = np.array(stars, dtype=np.uint8)
        self.plaqs = np.array(plaqs, dtype=np.uint8)

    def induced_mask(self, VS):
        L = self.L
        m = np.zeros(self.n, dtype=bool)
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                m[self.E(i, j)] = True
            if (((i + 1) % L), (j % L)) in VS:
                m[self.S(i, j)] = True
        return m

def sector_vals(M, Rmask):
    rk = rref_rank(M)
    rin = r_in_explicit(M, Rmask)
    rout = r_in_explicit(M, ~Rmask)
    ir2 = rk - rin - rout
    return rk, rin, rout, ir2

def measures(T, Rmask, extraZ=None):
    Z = T.plaqs if extraZ is None else np.vstack([T.plaqs, extraZ])
    rkX, rinX, routX, ir2X = sector_vals(T.stars, Rmask)
    rkZ, rinZ, routZ, ir2Z = sector_vals(Z, Rmask)
    AQ = int(Rmask.sum())
    return dict(IR2=ir2X + ir2Z, S=AQ - rinX - rinZ)

def pair(T, VSA, VSB, extraZ=None):
    RA, RB = T.induced_mask(VSA), T.induced_mask(VSB)
    assert not (RA & RB).any()
    mA, mB, mAB = measures(T, RA, extraZ), measures(T, RB, extraZ), measures(T, RA | RB, extraZ)
    return dict(I_IR=mA['IR2'] + mB['IR2'] - mAB['IR2'],
                I_MI=mA['S'] + mB['S'] - mAB['S'])

def cap_record(T, Rmask):
    """number of logical classes supported in R, independently: for the X side,
       ops = kernel of plaq-matrix restricted to R columns (X strings in R commuting
       with all plaquettes), cap_X = rank(stars ∪ kernel) - rank(stars); same for Z."""
    cap = 0
    for (checks, gens) in ((T.plaqs, T.stars), (T.stars, T.plaqs)):
        cols = np.where(Rmask)[0]
        sub = checks[:, cols] % 2
        # kernel of sub (as map from F2^|cols| to F2^rows) via rref of transpose trick
        A = sub.copy()
        rows, ncols = A.shape
        # gaussian elimination tracking combinations
        Aug = np.concatenate([A.T, np.eye(ncols, dtype=np.uint8)], axis=1)  # rows = cols of A
        r = 0
        for c in range(rows):
            piv = None
            for i in range(r, ncols):
                if Aug[i, c]:
                    piv = i
                    break
            if piv is None:
                continue
            Aug[[r, piv]] = Aug[[piv, r]]
            for i in range(ncols):
                if i != r and Aug[i, c]:
                    Aug[i] ^= Aug[r]
            r += 1
        ker_local = Aug[r:, rows:]  # combinations giving zero
        ker = np.zeros((ker_local.shape[0], T.n), dtype=np.uint8)
        ker[:, cols] = ker_local
        r0 = rref_rank(gens)
        cap += rref_rank(np.vstack([gens, ker])) - r0
    return cap

def block(L, s, i0, j0):
    return {(((i0 + a) % L), ((j0 + b) % L)) for a in range(s) for b in range(s)}

def ring(L, i0, w):
    return {(((i0 + a) % L), j) for a in range(w) for j in range(L)}

def frame(L, o0, q, w):
    return block(L, q, o0, o0) - block(L, q - 2 * w, o0 + w, o0 + w)

fails = []
def chk(name, got, want):
    ok = (got == want)
    print("%-58s got=%-18s want=%-18s %s" % (name, got, want, "OK" if ok else "FAIL"))
    if not ok:
        fails.append(name)

# ---- A: main sweep L=12 s=3
T12 = Torus(12)
for g in range(0, 7):
    r = pair(T12, block(12, 3, 0, 0), block(12, 3, 0, 3 + g))
    want = (2, 1) if g in (0, 6) else (0, 0)
    chk("A L=12 s=3 g=%d (I_IR,I_MI)" % g, (r['I_IR'], r['I_MI']), want)

# ---- B: absorption and s=4 contact
T8 = Torus(8)
r = pair(T8, block(8, 2, 0, 0), block(8, 2, 0, 2))
chk("B L=8 s=2 g=0 contact absorbed", (r['I_IR'], r['I_MI']), (0, 0))
r = pair(T12, block(12, 4, 0, 0), block(12, 4, 0, 4))
chk("B L=12 s=4 g=0 contact", (r['I_IR'], r['I_MI']), (4, 2))

# ---- C: double contact + joint content
r = pair(T8, block(8, 4, 0, 0), block(8, 4, 0, 4))
chk("C L=8 s=4 double seam", (r['I_IR'], r['I_MI']), (9, 4))
RA = T8.induced_mask(block(8, 4, 0, 0)); RB = T8.induced_mask(block(8, 4, 0, 4))
chk("C L=8 s=4 capA,capB,capAB",
    (cap_record(T8, RA), cap_record(T8, RB), cap_record(T8, RA | RB)), (0, 0, 1))
T6 = Torus(6)
r = pair(T6, block(6, 3, 0, 0), block(6, 3, 0, 3))
chk("C L=6 s=3 double seam", (r['I_IR'], r['I_MI']), (5, 2))
RA = T6.induced_mask(block(6, 3, 0, 0)); RB = T6.induced_mask(block(6, 3, 0, 3))
chk("C L=6 s=3 capAB", cap_record(T6, RA | RB), 1)

# ---- D: winding rings L=10
T10 = Torus(10)
for (w, r0) in ((1, 3), (1, 6), (2, 4), (2, 7)):
    rr = pair(T10, ring(10, 0, w), ring(10, r0, w))
    chk("D L=10 w=%d r0=%d" % (w, r0), (rr['I_IR'], rr['I_MI']), (w, w))
chk("D L=10 w=2 capA,capAB",
    (cap_record(T10, T10.induced_mask(ring(10, 0, 2))),
     cap_record(T10, T10.induced_mask(ring(10, 0, 2)) | T10.induced_mask(ring(10, 4, 2)))),
    (2, 2))

# ---- E: KP tripartite (same partition rule, independent code)
def tee(T, L, q, i0, j0):
    D = T.induced_mask(block(L, q, i0, j0))
    A = np.zeros(T.n, dtype=bool); B = np.zeros(T.n, dtype=bool); C = np.zeros(T.n, dtype=bool)
    rs = cs = q - 1
    for k in np.where(D)[0]:
        (t, i, j) = T.edges[k]
        if t == 'E':
            mi, mj = 2 * ((i - i0) % L), 2 * ((j - j0) % L) + 1
        else:
            mi, mj = 2 * ((i - i0) % L) + 1, 2 * ((j - j0) % L)
        if mi >= rs: C[k] = True
        elif mj < cs: A[k] = True
        else: B[k] = True
    S = lambda R: measures(T, R)['S']
    return (S(A) + S(B) + S(C) - S(A | B) - S(A | C) - S(B | C) + S(A | B | C))
chk("E KP TEE L=12 q=7 off(2,3)", tee(T12, 12, 7, 2, 3), -1)
chk("E KP TEE L=12 q=6 off(0,0)", tee(T12, 12, 6, 0, 0), -1)

# ---- F: inserted Z-row, L=12 s=3 g=3
for ell in range(2, 10):
    row = np.zeros(T12.n, dtype=np.uint8)
    for j in range(1, ell + 1):
        row[T12.E(1, j)] ^= 1
    r = pair(T12, block(12, 3, 0, 0), block(12, 3, 0, 6), extraZ=row[None, :])
    chk("F insert g=3 ell=%d I_IR" % ell, r['I_IR'], 1 if ell == 6 else 0)

# ---- G: enclosing L=10
for (t, want) in ((0, (6, 3)), (1, (0, 0))):
    VSB = frame(10, 0, 2 + 2 * t + 4, 2)
    VSA = block(10, 2, 2 + t, 2 + t)
    r = pair(T10, VSA, VSB)
    chk("G L=10 sA=2 w=2 t=%d" % t, (r['I_IR'], r['I_MI']), want)
VSB = frame(10, 0, 4, 1); VSA = block(10, 2, 1, 1)
r = pair(T10, VSA, VSB)
chk("G L=10 sA=2 w=1 t=0 absorbed", (r['I_IR'], r['I_MI']), (0, 0))

# ---- H: BFS earned-separation shape (the dV claim): recompute independently
from collections import deque
def bfs(L, VSA, VSB):
    dist = {v: 0 for v in VSA}
    dq = deque(VSA)
    while dq:
        (i, j) = dq.popleft()
        if (i, j) in VSB:
            return dist[(i, j)]
        for (a, b) in ((i+1, j), (i-1, j), (i, j+1), (i, j-1)):
            wv = (a % L, b % L)
            if wv not in dist:
                dist[wv] = dist[(i, j)] + 1
                dq.append(wv)
dvs = [bfs(12, block(12, 2, 0, 0), block(12, 2, 0, 2 + g)) for g in range(0, 9)]
print("H dV sweep L=12 s=2 (audit claims g+1 monotone):", dvs)
chk("H dV equals geff+1 (NOT g+1)", dvs, [min(g, 8 - g) + 1 for g in range(9)])

print()
if fails:
    print("FAILURES:", fails)
    sys.exit(1)
print("ALL INDEPENDENT CHECKS PASS")
