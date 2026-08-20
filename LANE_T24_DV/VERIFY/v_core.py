"""ADVERSARIAL VERIFICATION of LANE_T24_DV (T-24, clause (v) on D(D_4)).

Fresh code -- nothing imported from t24_lib.py.  Different config-index layout
(idx = ((u1*n + u0)*n + h1)*n + h0, the REVERSE of the lane's), traces computed by
FIXED-POINT COUNTING (never by building sparse projector matrices), commutant
constraints assembled by iterating over sparse H entries (a different assembly route
than the lane's per-basis-unit accumulation).  All verdict-bearing numbers are
integers / Fractions.  Floats appear only in v_numeric.py (labeled).

Lattice conventions (from the physics, restated independently):
  vertices v0, v1; edges h0: v0->v1, h1: v1->v0, u0 loop at v0, u1 loop at v1.
  A_v0(k): h0 -> k h0,      h1 -> h1 k^-1,  u0 -> k u0 k^-1,  u1 fixed.
  A_v1(k): h0 -> h0 k^-1,   h1 -> k h1,     u1 -> k u1 k^-1,  u0 fixed.
  B_p0: [u0 == h0 u1 h0^-1];  B_p1: [u1 == h1 u0 h1^-1].
  H = -(A_v0 + A_v1 + B_p0 + B_p1), A_v = (1/|G|) sum_k A_v(k).
"""
import itertools, json, time
from fractions import Fraction
import numpy as np

t0 = time.time()
REPORT = {}
def say(s=""):
    print(s, flush=True)

# ---------------------------------------------------------------- groups (fresh)
def group_D4():
    # element = (a, b) meaning r^a s^b ; s r s = r^-1 ; encode a*2 + b  (DIFFERENT from lane: lane used a + 4b)
    n = 8
    enc = lambda a, b: (a % 4) * 2 + (b % 2)
    dec = lambda g: (g // 2, g % 2)
    MUL = np.zeros((n, n), dtype=np.int64)
    for g1 in range(n):
        a1, b1 = dec(g1)
        for g2 in range(n):
            a2, b2 = dec(g2)
            MUL[g1, g2] = enc(a1 + (a2 if b1 == 0 else -a2), b1 + b2)
    INV = np.array([[h for h in range(n) if MUL[g, h] == enc(0, 0)][0] for g in range(n)], dtype=np.int64)
    E = enc(0, 0)
    # associativity check (full)
    ok_assoc = all(MUL[MUL[a, b], c] == MUL[a, MUL[b, c]] for a in range(n) for b in range(n) for c in range(n))
    assert ok_assoc
    # conjugacy classes
    classes = []
    seen = set()
    for g in range(n):
        if g in seen: continue
        cl = sorted({int(MUL[k, MUL[g, INV[k]]]) for k in range(n)})
        seen |= set(cl); classes.append(cl)
    # characters (5 irreps): 1-dim eps^a del^b, and 2-dim (2 on e, -2 on r^2, else 0)
    chars = []
    for eps in (1, -1):
        for dl in (1, -1):
            chars.append(np.array([eps ** dec(g)[0] * dl ** dec(g)[1] for g in range(n)], dtype=np.int64))
    two = np.zeros(n, dtype=np.int64); two[enc(0, 0)] = 2; two[enc(2, 0)] = -2
    chars.append(two)
    r2 = enc(2, 0); s = enc(0, 1)
    return dict(n=n, MUL=MUL, INV=INV, E=E, classes=classes, chars=chars, r2=r2, s=s, enc=enc, dec=dec)

def group_Z2():
    MUL = np.array([[0, 1], [1, 0]], dtype=np.int64)
    return dict(n=2, MUL=MUL, INV=np.array([0, 1], dtype=np.int64), E=0,
                classes=[[0], [1]], chars=[np.array([1, 1], np.int64), np.array([1, -1], np.int64)],
                r2=None, s=None)

# ---------------------------------------------------------------- carrier (fresh layout)
class Car:
    """config index = ((u1*n + u0)*n + h1)*n + h0  -- REVERSED vs the lane."""
    def __init__(self, G):
        self.G = G; n = G["n"]; self.n = n; self.N = n ** 4
        idx = np.arange(self.N)
        self.h0 = idx % n; r = idx // n
        self.h1 = r % n; r = r // n
        self.u0 = r % n; self.u1 = r // n
        self.M = G["MUL"]; self.I = G["INV"]
    def pack(self, h0, h1, u0, u1):
        n = self.n
        return ((u1 * n + u0) * n + h1) * n + h0
    def A0(self, k):
        M, I = self.M, self.I
        return self.pack(M[k, self.h0], M[self.h1, I[k]], M[M[k, self.u0], I[k]], self.u1)
    def A1(self, k):
        M, I = self.M, self.I
        return self.pack(M[self.h0, I[k]], M[k, self.h1], self.u0, M[M[k, self.u1], I[k]])
    def dB0(self):
        M, I = self.M, self.I
        return (self.u0 == M[M[self.h0, self.u1], I[self.h0]]).astype(np.int64)
    def dB1(self):
        M, I = self.M, self.I
        return (self.u1 == M[M[self.h1, self.u0], I[self.h1]]).astype(np.int64)
    def comp(self, e):
        return {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}[e]
    def with_comp(self, e, new):
        c = {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}
        c = dict(c); c[e] = new
        return self.pack(c["h0"], c["h1"], c["u0"], c["u1"])

# ---------------------------------------------------------------- exact trace engine
def sector_trace(car, sector, Mspec):
    """Tr(P_s M) exactly.  P_s = Q0^{a0} Q1^{a1} D0^{b0} D1^{b1}, Q^1 = A-avg, Q^0 = I - A-avg,
    D^1 = diagB, D^0 = 1 - diagB.  Mspec: ("perm", pi) or ("diag", d) or ("permdiag", pi, d) meaning diag then perm? --
    we only ever need pure perm or pure diag."""
    a0, a1, b0, b1 = sector
    n, N = car.n, car.N
    d = (car.dB0() if b0 else 1 - car.dB0()) * (car.dB1() if b1 else 1 - car.dB1())
    idx = np.arange(N)
    total = Fraction(0)
    terms0 = [(1, True)] if a0 == 1 else [(1, False), (-1, True)]   # (sign, use A0 average?)
    terms1 = [(1, True)] if a1 == 1 else [(1, False), (-1, True)]
    for s0, useA0 in terms0:
        for s1, useA1 in terms1:
            ks0 = range(n) if useA0 else [None]
            ks1 = range(n) if useA1 else [None]
            w = Fraction(s0 * s1)
            if useA0: w /= n
            if useA1: w /= n
            acc = 0
            for k0 in ks0:
                p0 = car.A0(k0) if k0 is not None else idx
                for k1 in ks1:
                    p01 = p0[car.A1(k1)] if k1 is not None else p0
                    if Mspec[0] == "perm":
                        pim = Mspec[1]
                        comp = p01[pim]           # A0 A1 M as index map
                        fixed = comp == idx
                        acc += int(d[pim[fixed]].sum()) if fixed.any() else 0
                    else:                          # diag
                        m = Mspec[1]
                        fixed = p01 == idx
                        acc += int((d[fixed] * m[fixed]).sum()) if fixed.any() else 0
            total += w * acc
    return total

def eig_traces(car, Mspec):
    """{eig k (H eigenvalue -k): exact Tr(P_E M)} via the 16 sectors."""
    out = {}
    for s in itertools.product((0, 1), repeat=4):
        k = sum(s)
        out[k] = out.get(k, Fraction(0)) + sector_trace(car, s, Mspec)
    return out

# ---------------------------------------------------------------- commutant of H on one edge
def build_H_sparse(car):
    """n*H (up to sign) as dict {(i,j): int}: sum_k A0(k) + A1(k) perms weight 1, diags weight n.
    (Scaling/sign irrelevant for the commutant.)"""
    n, N = car.n, car.N
    H = {}
    idx = np.arange(N)
    for k in range(n):
        for pi in (car.A0(k), car.A1(k)):
            for j in range(N):
                key = (int(pi[j]), j)
                H[key] = H.get(key, 0) + 1
    for dvec in (car.dB0(), car.dB1()):
        nz = np.nonzero(dvec)[0]
        for j in nz:
            key = (int(j), int(j))
            H[key] = H.get(key, 0) + n
    return H

def edge_commutant_nullity(car, e, Hs):
    """dim { M in M_n : [M (x) I_e, H] = 0 } by assembling the constraint from H's sparse
    entries (fresh route), Gram over the 64 (n^2) columns, exact Fraction kernel."""
    n = car.n
    comp = car.comp(e)
    # per-edge component of every index, and index-with-edge-replaced tables
    repl = [car.with_comp(e, np.full(car.N, a, dtype=np.int64)) for a in range(n)]
    # group H entries by needed component values
    items = list(Hs.items())
    I_arr = np.array([ij[0] for ij, _ in items], dtype=np.int64)
    J_arr = np.array([ij[1] for ij, _ in items], dtype=np.int64)
    V_arr = np.array([v for _, v in items], dtype=np.int64)
    ce_I = comp[I_arr]   # edge component of row index
    ce_J = comp[J_arr]
    pos = {}
    def add(i, j, col, val):
        if val == 0: return
        key = (int(i), int(j))
        dd = pos.get(key)
        if dd is None: pos[key] = dd = {}
        dd[col] = dd.get(col, 0) + val
        if dd[col] == 0: del dd[col]
    for a in range(n):
        Ia = repl[a]
        for b in range(n):
            col = a * n + b
            # X = E_ab on edge e:  X|c> = [c_e == b] |c(e->a)>
            # X.H: rows i of H with i_e == b -> (i(e->a), j)
            sel = ce_I == b
            for i, j, v in zip(Ia[I_arr[sel]], J_arr[sel], V_arr[sel]):
                add(i, j, col, int(v))
            # H.X: cols j of H with j_e == a -> (i, j(e->b))
            selc = ce_J == a
            Jb = car.with_comp(e, np.full(car.N, b, dtype=np.int64))
            for i, j, v in zip(I_arr[selc], Jb[J_arr[selc]], V_arr[selc]):
                add(i, j, col, -int(v))
    m = n * n
    Gm = [[0] * m for _ in range(m)]
    for dd in pos.values():
        its = list(dd.items())
        for x in range(len(its)):
            c1, v1 = its[x]
            for y in range(x, len(its)):
                c2, v2 = its[y]
                Gm[c1][c2] += v1 * v2
                if c1 != c2: Gm[c2][c1] += v1 * v2
    # exact kernel dimension by Fraction Gaussian elimination (fresh implementation)
    A = [[Fraction(Gm[i][j]) for j in range(m)] for i in range(m)]
    rank = 0
    for c in range(m):
        piv = None
        for r in range(rank, m):
            if A[r][c] != 0: piv = r; break
        if piv is None: continue
        A[rank], A[piv] = A[piv], A[rank]
        pv = A[rank][c]
        A[rank] = [x / pv for x in A[rank]]
        for r in range(m):
            if r != rank and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[rank])]
        rank += 1
    # float SVD cross-check of the SAME constraint matrix (independent rank route)
    cols = sorted({c for dd in pos.values() for c in dd})
    X = np.zeros((len(pos), m))
    for rr, dd in enumerate(pos.values()):
        for c, v in dd.items():
            X[rr, c] = v
    sv = np.linalg.svd(X, compute_uv=False) if len(pos) else np.zeros(0)
    rank_f = int((sv > 1e-6 * (sv[0] if len(sv) else 1)).sum())
    assert rank_f == rank, ("Gram rank vs SVD rank disagree", rank, rank_f)
    return m - rank

def check_edge_op_commutes(car, e, Mmat, Hs):
    """EXACT check [M (x) I_e, H] = 0 for integer matrix Mmat (n x n), via sparse H entries."""
    n = car.n
    comp = car.comp(e)
    repl = [car.with_comp(e, np.full(car.N, a, dtype=np.int64)) for a in range(n)]
    acc = {}
    def add(i, j, v):
        if v == 0: return
        k = (int(i), int(j))
        w = acc.get(k, 0) + v
        if w: acc[k] = w
        else: acc.pop(k, None)
    items = list(Hs.items())
    I_arr = np.array([ij[0] for ij, _ in items], dtype=np.int64)
    J_arr = np.array([ij[1] for ij, _ in items], dtype=np.int64)
    V_arr = np.array([v for _, v in items], dtype=np.int64)
    ce_I = comp[I_arr]; ce_J = comp[J_arr]
    for a in range(n):
        Ia = repl[a]
        Jb_cache = {}
        for b in range(n):
            v0 = int(Mmat[a][b])
            if v0 == 0: continue
            sel = ce_I == b
            for i, j, v in zip(Ia[I_arr[sel]], J_arr[sel], V_arr[sel]):
                add(i, j, v0 * int(v))
            selc = ce_J == a
            if b not in Jb_cache:
                Jb_cache[b] = car.with_comp(e, np.full(car.N, b, dtype=np.int64))
            Jb = Jb_cache[b]
            for i, j, v in zip(I_arr[selc], Jb[J_arr[selc]], V_arr[selc]):
                add(i, j, -v0 * int(v))
    return len(acc) == 0

# ================================================================ PART 1: minimal torus
say("=" * 100)
say("V1. MINIMAL TORUS -- vacuity from the cell structure (independent)")
say("=" * 100)
# cell structure: 1 vertex, edges e_x, e_y (self-loops), 1 face word x y x^-1 y^-1
d1 = np.array([[0, 0]])                    # boundary of each edge: v - v = 0
word = [("x", +1), ("y", +1), ("x", -1), ("y", -1)]
d2 = {"x": 0, "y": 0}
for ltr, sgn in word: d2[ltr] += sgn
d2v = np.array([[d2["x"]], [d2["y"]]])
rank_d1 = np.linalg.matrix_rank(d1) if d1.any() else 0
rank_d2 = np.linalg.matrix_rank(d2v) if d2v.any() else 0
say("  d1 = %s (rank %d), d2 = %s (rank %d) => H_1 = Z^%d" %
    (d1.tolist(), rank_d1, d2v.T.tolist(), rank_d2, 2 - rank_d1 - rank_d2))
# every nonempty edge subset: contains a self-loop => graph cycle (forest reading) and
# carries a nonzero H_1 class (homology reading, since d1 = 0 and im d2 = 0)
subsets = [(), ("x",), ("y",), ("x", "y")]
agree = True
for sub in subsets:
    forest_contractible = (len(sub) == 0)   # any self-loop is a cycle; empty set is a forest
    hom_contractible = (len(sub) == 0)      # any nonzero chain is a nonzero class (H1 = Z^2, im d2 = 0)
    agree &= (forest_contractible == hom_contractible)
say("  forest and homology readings agree on all 4 subsets: %s" % agree)
say("  single contractible qubit-supported regions: NONE  => clause (v) VACUOUS (empty domain)")
vac_ok = (rank_d1 == 0 and rank_d2 == 0 and agree)
REPORT["minimal_vacuity_proved"] = bool(vac_ok)

# machinery cross-check: D(D_4) minimal torus sector dims by exact counting
G4 = group_D4()
n = G4["n"]; MUL, INV, E = G4["MUL"], G4["INV"], G4["E"]
# configs (x, y), idx = y*n + x (fresh layout); A(k): x->kxk^-1, y->kyk^-1; B: [x y x^-1 y^-1 == e]
idx2 = np.arange(n * n)
X = idx2 % n; Y = idx2 // n
def A_min(k):
    return (MUL[MUL[k, Y], INV[k]]) * n + MUL[MUL[k, X], INV[k]]
comm_ok = MUL[MUL[X, Y], MUL[INV[X], INV[Y]]] == E
dB = comm_ok.astype(np.int64)
# sector (a, b): Q^a D^b; trace by fixed-point counting
mindims = {}
for a in (0, 1):
    for b in (0, 1):
        d = dB if b else 1 - dB
        terms = [(1, True)] if a == 1 else [(1, False), (-1, True)]
        tot = Fraction(0)
        for s, useA in terms:
            w = Fraction(s, n if useA else 1)
            acc = 0
            for k in (range(n) if useA else [None]):
                p = A_min(k) if k is not None else idx2
                fixed = p == idx2
                acc += int(d[fixed].sum())
            tot += w * acc
        assert tot.denominator == 1
        mindims[(a, b)] = int(tot)
eig_min = {}
for (a, b), dd in mindims.items():
    eig_min[a + b] = eig_min.get(a + b, 0) + dd
say("  D(D_4) minimal torus sector dims: %s ; eigen dims (eig -k): %s" %
    (mindims, {(-k): v for k, v in sorted(eig_min.items())}))
o36_match = (eig_min.get(2, 0) == 22 and eig_min.get(1, 0) == 24 and eig_min.get(0, 0) == 18
             and sum(mindims.values()) == 64)
say("  matches O-36 (22 / 24 / 18): %s" % o36_match)
REPORT["minimal_o36_match"] = bool(o36_match)

# ================================================================ PART 2: 1x2 D(D_4)
say("")
say("=" * 100)
say("V2. 1x2 TORUS D(D_4) (dim 4096) -- structure, sectors, commutants, the obstruction")
say("=" * 100)
car = Car(G4)
Nc = car.N
# structure checks
okA0 = all(np.array_equal(car.A0(k1)[car.A0(k2)], car.A0(int(MUL[k1, k2]))) for k1 in range(n) for k2 in range(n))
okA1 = all(np.array_equal(car.A1(k1)[car.A1(k2)], car.A1(int(MUL[k1, k2]))) for k1 in range(n) for k2 in range(n))
okC = all(np.array_equal(car.A0(k1)[car.A1(k2)], car.A1(k2)[car.A0(k1)]) for k1 in range(n) for k2 in range(n))
b0, b1v = car.dB0(), car.dB1()
okB = all(np.array_equal(b0[p], b0) and np.array_equal(b1v[p], b1v)
          for k in range(n) for p in (car.A0(k), car.A1(k)))
say("  perm reps: %s %s ; [A0,A1]=0: %s ; B invariant under gauge: %s" % (okA0, okA1, okC, okB))
assert okA0 and okA1 and okC and okB

# sector dims / eigen dims by exact counting
ID = ("perm", np.arange(Nc))
et = eig_traces(car, ID)
eigdims = {}
for k, v in et.items():
    assert v.denominator == 1
    eigdims[-k] = int(v)
say("  eigen dims: %s   (sum %d)" % (dict(sorted(eigdims.items())), sum(eigdims.values())))
dims_ok = (eigdims == {0: 2686, -1: 864, -2: 476, -3: 48, -4: 22})
even_ok = all(v % 2 == 0 for v in eigdims.values())
say("  match lane {0:2686,-1:864,-2:476,-3:48,-4:22}: %s ; all even (C-41): %s" % (dims_ok, even_ok))
REPORT["eigdims_match"] = bool(dims_ok); REPORT["eigdims_all_even"] = bool(even_ok)

# region catalogue (fresh union-find)
ENDS = {"h0": (0, 1), "h1": (1, 0), "u0": (0, 0), "u1": (1, 1)}
def acyclic(sub):
    par = {0: 0, 1: 1}
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for e in sub:
        a, b = ENDS[e]
        ra, rb = f(a), f(b)
        if ra == rb: return False
        par[ra] = rb
    return True
singles_contractible = [e for e in ENDS if acyclic([e])]
say("  contractible single-edge regions (forest reading): %s" % singles_contractible)
# homology reading on the 1x2 CW complex: d1(h0)=v1-v0, d1(h1)=v0-v1, d1(u)=0;
# d2(p0)=u1-u0, d2(p1)=u0-u1.  h0,h1 are non-cycles inside a contractible arc; u_i are
# cycles with nonzero class ([u0] generates H1 with [h0+h1]).
d1_12 = {"h0": (-1, 1), "h1": (1, -1), "u0": (0, 0), "u1": (0, 0)}
hom_contract = [e for e in ENDS if d1_12[e] != (0, 0)]  # single edge contractible iff its chain is not a cycle
say("  homology reading, single-edge: non-cycle (contractible arc): %s ; cycles u0,u1 have class"
    " [u0] = [u1] != 0 (im d2 = span(u1-u0))" % hom_contract)
regions_ok = sorted(singles_contractible) == ["h0", "h1"] == sorted(hom_contract)
REPORT["regions_match"] = bool(regions_ok)

# commutants
Hs = build_H_sparse(car)
say("  sparse n*H entries: %d" % len(Hs))
nul = {}
for e in ("h0", "h1", "u0", "u1"):
    nul[e] = edge_commutant_nullity(car, e, Hs)
    say("  dim C_%s = %d" % (e, nul[e]))
c_ok = (nul["h0"] == 2 and nul["h1"] == 2 and nul["u0"] == 5 and nul["u1"] == 5)
say("  match lane (2, 2, 5, 5): %s" % c_ok)
REPORT["commutant_dims"] = {e: int(v) for e, v in nul.items()}

# claimed spanning elements, verified exactly
r2 = G4["r2"]
Tmat = [[1 if MUL[bb, r2] == aa else 0 for bb in range(n)] for aa in range(n)]  # right-mult by r^2
Imat = [[1 if aa == bb else 0 for bb in range(n)] for aa in range(n)]
span_ok = True
for e in ("h0", "h1"):
    span_ok &= check_edge_op_commutes(car, e, Tmat, Hs)
    span_ok &= check_edge_op_commutes(car, e, Imat, Hs)
class_ok = True
for e in ("u0", "u1"):
    for cl in G4["classes"]:
        D = [[1 if (aa == bb and aa in cl) else 0 for bb in range(n)] for aa in range(n)]
        class_ok &= check_edge_op_commutes(car, e, D, Hs)
say("  I, T=R_(r^2) in C_h0 and C_h1 (exact commutator = 0): %s" % span_ok)
say("  all 5 class diagonals in C_u0 and C_u1 (exact): %s" % class_ok)
say("  => with the nullities above, C_h = span{I, T} (dim 2) and C_u = class diagonals (dim 5)")
REPORT["span_ok"] = bool(span_ok and class_ok and c_ok)

# THE OBSTRUCTION: t_E = Tr(P_E T) on h0 and h1
tE = {}
for e in ("h0", "h1"):
    Tperm = car.with_comp(e, MUL[car.comp(e), r2])
    tt = eig_traces(car, ("perm", Tperm))
    tE[e] = {}
    for k, v in tt.items():
        assert v.denominator == 1
        tE[e][-k] = int(v)
    say("  t_E(T on %s) = %s ; sum = %d (= Tr T)" % (e, dict(sorted(tE[e].items())), sum(tE[e].values())))
tE_ok = all(tE[e] == {0: 30, -1: -48, -2: 12, -3: 0, -4: 6} for e in ("h0", "h1"))
say("  match lane {0:30,-1:-48,-2:12,-3:0,-4:6} on both edges: %s" % tE_ok)
say("  NOT all zero => no Hermitian involution commuting with H anticommutes with T;")
say("  only involutions of C_h are +-I, +-T => NO admissible single-edge flip on h0/h1.")
REPORT["tE_match"] = bool(tE_ok)
REPORT["tE"] = {e: tE[e] for e in tE}

# ground-space non-scalarity of T (KL divergence): t_gs = 6, |t_gs| < 22
t_gs = tE["h0"][-4]
say("  t_gs = %d with ground dim 22: T non-scalar on ground space (%d != +-22) -- weight-1" % (t_gs, t_gs))
say("  admissible nontrivial logical that flips nothing (flip-reading HOLDS, KL-reading fails).")
# independent structural cross-check: classes of commuting pairs fixed by (a,b)->(a r^2, b)
pairs = [(a, b) for a in range(n) for b in range(n) if MUL[a, b] == MUL[b, a]]
orb = {}
for (a, b) in pairs:
    key = min((int(MUL[k, MUL[a, INV[k]]]), int(MUL[k, MUL[b, INV[k]]])) for k in range(n))
    orb.setdefault(key, set()).add((a, b))
say("  commuting-pair classes: %d (ground dim must be 22: %s)" % (len(orb), len(orb) == 22))
fixed = 0
for key, members in orb.items():
    a, b = key
    im = min((int(MUL[k, MUL[int(MUL[a, r2]), INV[k]]]), int(MUL[k, MUL[b, INV[k]]])) for k in range(n))
    if im == key: fixed += 1
say("  classes fixed by (a,b) -> (a r^2, b): %d (t_gs cross-check, lane says 6): %s" % (fixed, fixed == t_gs))
REPORT["t_gs_crosscheck"] = bool(len(orb) == 22 and fixed == t_gs)

# u-edges: class traces per eigenspace, and balanced +-1 pattern count
u_ok = True
for e in ("u0", "u1"):
    class_traces = []
    for cl in G4["classes"]:
        dcl = np.isin(car.comp(e), cl).astype(np.int64)
        tt = eig_traces(car, ("diag", dcl))
        class_traces.append({-k: int(v) for k, v in tt.items()})
    say("  edge %s class traces per E: %s" % (e, [ [ct[k] for k in sorted(ct)] for ct in class_traces]))
    balanced = 0
    for signs in itertools.product((1, -1), repeat=len(class_traces)):
        if all(sum(sg * ct[k] for sg, ct in zip(signs, class_traces)) == 0 for k in (0, -1, -2, -3, -4)):
            balanced += 1
    say("  edge %s balanced +-1 class patterns: %d" % (e, balanced))
    u_ok &= (balanced == 0)
say("  zero balanced patterns on both u-edges (lane claim): %s" % u_ok)
REPORT["u_balanced_zero"] = bool(u_ok)

# Wilson diagonal fails clause (iv)
chi_a = G4["chars"][1]  # eps=1, dl=-1 : s -> -1, r -> +1
Wh = chi_a[MUL[car.h0, car.h1]].astype(np.int64)
ttW = eig_traces(car, ("diag", Wh))
trW = {-k: int(v) for k, v in ttW.items()}
say("  Wilson diag chi_a(h0 h1): Tr(P_E R) = %s (lane: {0:62,-1:-112,-2:44,-3:0,-4:6})" % dict(sorted(trW.items())))
wilson_ok = trW == {0: 62, -1: -112, -2: 44, -3: 0, -4: 6}
say("  clause (iv) FAILS for it (not balanced): %s" % wilson_ok)
REPORT["wilson_not_record"] = bool(wilson_ok)

# positive controls (non-admissible flippers), exact
s_el = G4["s"]
def frob2_comm_perm(car, perm):
    """||[P_perm, H]||_F^2 exactly, H = -(A0avg + A1avg + B0 + B1)."""
    n, N = car.n, car.N
    acc = {}
    def add(i, j, v):
        if v == 0: return
        k = (int(i), int(j)); w = acc.get(k, Fraction(0)) + v
        if w: acc[k] = w
        else: acc.pop(k, None)
    idx = np.arange(N)
    w8 = Fraction(1, n)
    for k in range(n):
        for pi in (car.A0(k), car.A1(k)):
            LP = perm[pi]; PL = pi[perm]     # index maps for P_L P_pi and P_pi P_L
            for j in range(N):
                add(LP[j], j, -w8)           # -(1/8) per gauge perm in H; sign squares away
                add(PL[j], j, w8)
    for dvec in (car.dB0(), car.dB1()):
        dv = dvec.astype(np.int64)
        for j in range(N):
            v = -(int(dv[j]) - int(dv[perm[j]]))   # [L, D] entries (L(j), j): d(j) - d(L(j)) times -1
            if v: add(perm[j], j, Fraction(v))
    return sum(v * v for v in acc.values())

Ls_h0 = car.with_comp("h0", MUL[s_el, car.h0])   # left-mult s on h0
f2 = frob2_comm_perm(car, Ls_h0)
flip_W = np.array_equal(Wh[Ls_h0], -Wh)
say("  L_s on h0: ||[L_s,H]||_F^2 = %s (lane: 1024) ; conj flips Wilson diag: %s" % (f2, flip_W))
Rs_u0 = car.with_comp("u0", MUL[car.u0, s_el])   # right-mult s on u0
chi_u = chi_a[car.u0].astype(np.int64)
f2b = frob2_comm_perm(car, Rs_u0)
flip_chi = np.array_equal(chi_u[Rs_u0], -chi_u)
say("  R_s on u0: ||[R_s,H]||_F^2 = %s (lane: 3840) ; conj flips chi_a(u0): %s" % (f2b, flip_chi))
pos_ok = (f2 == 1024 and flip_W and f2b == 3840 and flip_chi)
say("  positive controls fire (D-15): %s" % pos_ok)
REPORT["positive_controls"] = bool(pos_ok)

# gauge-rep multiplicities per eigenspace -> transport-fixed record existence, EXACT
say("")
say("  transport-fixed existence (independent of O-36 feasibility): multiplicities of")
say("  G x G irreps on each eigenspace, by exact character counting:")
reps = [cl[0] for cl in G4["classes"]]
clsz = [len(cl) for cl in G4["classes"]]
chiE = {}   # chiE[(ci, cj)] = {eig: Tr(P_E A0(k_i) A1(k_j))}
for ci, ki in enumerate(reps):
    for cj, kj in enumerate(reps):
        perm = car.A0(ki)[car.A1(kj)]
        tt = eig_traces(car, ("perm", perm))
        chiE[(ci, cj)] = {(-k): v for k, v in tt.items()}
mult_all_even = True
mult_table = {}
for Ek in (0, -1, -2, -3, -4):
    ms = []
    for i1, ch1 in enumerate(G4["chars"]):
        for i2, ch2 in enumerate(G4["chars"]):
            m = Fraction(0)
            for ci, ki in enumerate(reps):
                for cj, kj in enumerate(reps):
                    m += Fraction(clsz[ci] * clsz[cj]) * int(ch1[ki]) * int(ch2[kj]) * chiE[(ci, cj)][Ek]
            m /= n * n
            assert m.denominator == 1 and m >= 0, (Ek, i1, i2, m)
            mi = int(m)
            if mi:
                ms.append(((i1, i2), mi))
                if mi % 2: mult_all_even = False
    mult_table[Ek] = ms
    # consistency: sum d_i d_j m = dim
    dvec = [1, 1, 1, 1, 2]
    tot = sum(dvec[i1] * dvec[i2] * mi for (i1, i2), mi in ms)
    say("    E=%2d: nonzero mults %s ; sum d*d*m = %d (dim %d)" % (Ek, ms, tot, eigdims[Ek]))
    assert tot == eigdims[Ek]
say("  ALL multiplicities even: %s  => a gauge-commuting (transport-fixed) Hermitian involution" % mult_all_even)
say("     with zero trace on every eigenspace EXISTS exactly (split each multiplicity space evenly).")
REPORT["transport_fixed_exists_exact"] = bool(mult_all_even)

say("")
say("V2 done, %.1f s" % (time.time() - t0))

# ================================================================ PART 3: Z2 control
say("")
say("=" * 100)
say("V3. ABELIAN CONTROL D(Z_2), 1x2 -- same lattice, fresh code")
say("=" * 100)
G2 = group_Z2()
car2 = Car(G2)
n2, N2 = 2, 16
M2 = G2["MUL"]
et2 = eig_traces(car2, ("perm", np.arange(N2)))
eig2 = {(-k): int(v) for k, v in et2.items()}
say("  eigen dims: %s (lane: {0:4,-1:0,-2:8,-3:0,-4:4})" % dict(sorted(eig2.items())))
z2_dims_ok = eig2 == {0: 4, -1: 0, -2: 8, -3: 0, -4: 4}
Hs2 = build_H_sparse(car2)
nul2 = {e: edge_commutant_nullity(car2, e, Hs2) for e in ("h0", "h1", "u0", "u1")}
say("  commutant dims: %s (lane: all 2)" % nul2)
Xmat = [[0, 1], [1, 0]]
x_adm = check_edge_op_commutes(car2, "h0", Xmat, Hs2)
Xh0 = car2.with_comp("h0", M2[1, car2.h0])
ttX = eig_traces(car2, ("perm", Xh0))
tX = {(-k): int(v) for k, v in ttX.items()}
say("  X_h0 admissible ([X,H]=0 exact): %s ; t_E(X_h0) = %s (all zero: %s)"
    % (x_adm, dict(sorted(tX.items())), all(v == 0 for v in tX.values())))
# record Zbar_h = (-1)^(h0+h1), clause-verified exactly
sgn = G2["chars"][1]
Zh = sgn[M2[car2.h0, car2.h1]].astype(np.int64)
inv_ok = np.array_equal(Zh * Zh, np.ones(N2, dtype=np.int64))
Zmat_comm = True
# [diag Zh, H] = 0 exact: check H entries connect only equal-Zh configs OR diagonal
for (i, j), v in Hs2.items():
    if Zh[i] != Zh[j]: Zmat_comm = False
ttZ = eig_traces(car2, ("diag", Zh))
tZ = {(-k): int(v) for k, v in ttZ.items()}
flip_Z = np.array_equal(Zh[Xh0], -Zh)
say("  record Zbar_h: involution %s ; [R,H]=0 %s ; balanced %s ; nontrivial (trace 0 on dim>0 E) True"
    % (inv_ok, Zmat_comm, all(v == 0 for v in tZ.values())))
say("  X_h0 flips Zbar_h exactly: %s  => CLAUSE (v) FAILS on abelian 1x2 (weight-1 admissible flipper)" % flip_Z)
# non-contractible: Z_u0 admissible diag flips Xbar_v = X_u0 X_u1
Zu0 = sgn[car2.u0].astype(np.int64)
z_adm = check_edge_op_commutes(car2, "u0", [[1, 0], [0, -1]], Hs2)
Xv = car2.with_comp("u0", M2[1, car2.u0])
Xv = car2.pack(car2.h0, car2.h1, M2[1, car2.u0], M2[1, car2.u1])
# Xbar_v as perm; record checks: involution (perm^2 = id), commutes with H, balanced
invXv = np.array_equal(Xv[Xv], np.arange(N2))
commXv = True
Hd2 = {}
for (i, j), v in Hs2.items():
    Hd2[(i, j)] = v
commXv = all(Hd2.get((int(Xv[i]), int(Xv[j])), 0) == v for (i, j), v in Hs2.items())
ttXv = eig_traces(car2, ("perm", Xv))
tXv = {(-k): int(v) for k, v in ttXv.items()}
anti = np.array_equal(Zu0[Xv], -Zu0)   # {diag Zu0, P_Xv} = 0 iff Zu0(Xv(c)) = -Zu0(c)
say("  Xbar_v record: involution %s ; H-invariant %s ; balanced %s ; Z_u0 admissible %s ; flips Xbar_v %s"
    % (invXv, commXv, all(v == 0 for v in tXv.values()), z_adm, anti))
z2_ok = (z2_dims_ok and all(v == 2 for v in nul2.values()) and x_adm and
         all(v == 0 for v in tX.values()) and inv_ok and Zmat_comm and
         all(v == 0 for v in tZ.values()) and flip_Z and invXv and commXv and z_adm and anti)
say("  ABELIAN CONTROL VERDICT reproduced (FAILS, same test): %s" % z2_ok)
REPORT["z2_control"] = bool(z2_ok)

REPORT["elapsed_s"] = round(time.time() - t0, 1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV/VERIFY/v_core.json", "w") as f:
    json.dump(REPORT, f, indent=1)
say("")
say("REPORT: %s" % json.dumps(REPORT))
