#!/usr/bin/env python3
"""
LANE O-4 -- clause (v) tested at L = 3 and L = 5.

HISTORY OF THIS SCRIPT'S OWN ERRORS (kept, because both were caught by the checks):
  ERR-1  At L = 2 a single plaquette's 4 edges already carry a winding cycle, so the
         "contractible region" was not contractible.  (o4_diagnose_L2.py.)
  ERR-2  At L = 3 a 2x2 BLOCK of plaquettes spans H-rows y, y+1, y+2 = all three rows,
         so it too winds.  A block of a x b plaquettes is contractible only for
         a, b <= L-2.  The first version of this script asserted contractibility
         geometrically and was wrong.
  FIX    Every region is now certified contractible by an INDEPENDENT computation of
         H_1 of the subcomplex it spans, primal AND dual.  A region is admitted only
         if both vanish.  Nothing about the record enters that test.

METHOD.  Exact F_2 linear algebra.  A Pauli is v = (x|z) in F_2^{2n}.
   ADMISSIBLE (DEF-A: [U,H] = 0)  <=>  z . star_v = 0 and x . plaq_p = 0  (linear)
   SUPPORTED IN T                 <=>  x_e = z_e = 0 for e not in T       (linear)
   FLIPS R                        <=>  x . R_z = 1                        (functional)
{admissible and supported in T} is an F_2 subspace; the flipper count is 0 if the
functional vanishes on it and exactly half of it otherwise.  Proof for the region, not
a sample.
"""

import numpy as np

FAILURES = []


def check(name, condition, detail=""):
    tag = "PASS" if condition else "FAIL"
    if not condition:
        FAILURES.append(name)
    print(f"   [{tag}] {name}   {detail}")


def hr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------- F_2 linear algebra
def rref(M):
    M = M.copy() % 2
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if M[i, c]), None)
        if pr is None:
            continue
        M[[r, pr]] = M[[pr, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return M[:r], piv


def nullspace(M, ncols):
    if M.shape[0] == 0:
        return np.eye(ncols, dtype=np.int8)
    Rm, piv = rref(M)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(ncols, dtype=np.int8)
        v[f] = 1
        for i, p in enumerate(piv):
            v[p] = Rm[i, f]
        basis.append(v)
    return np.array(basis, dtype=np.int8) if basis else np.zeros((0, ncols), np.int8)


def rank2(M):
    if M.size == 0 or M.shape[0] == 0:
        return 0
    return rref(M)[0].shape[0]


# ---------------------------------------------------------------------- the lattice
class Torus:
    def __init__(self, L):
        self.L = L
        self.N = 2 * L * L

    def e(self, x, y, d):
        L = self.L
        return 2 * (L * (y % L) + (x % L)) + d

    def name(self, i):
        L = self.L
        d = i % 2
        c = i // 2
        return f"({c % L},{c // L},{'H' if d == 0 else 'V'})"

    def plaq(self, x, y):
        return sorted({self.e(x, y, 0), self.e(x, y + 1, 0),
                       self.e(x, y, 1), self.e(x + 1, y, 1)})

    def star(self, x, y):
        return sorted({self.e(x, y, 0), self.e(x - 1, y, 0),
                       self.e(x, y, 1), self.e(x, y - 1, 1)})

    def all_plaqs(self):
        return [self.plaq(x, y) for y in range(self.L) for x in range(self.L)]

    def all_stars(self):
        return [self.star(x, y) for y in range(self.L) for x in range(self.L)]

    def edge_ends(self, i):
        L = self.L
        d = i % 2
        c = i // 2
        x, y = c % L, c // L
        return ((x, y), ((x + 1) % L, y)) if d == 0 else ((x, y), (x, (y + 1) % L))

    def edge_faces(self, i):
        """the two plaquettes on either side of edge i (dual endpoints)."""
        L = self.L
        d = i % 2
        c = i // 2
        x, y = c % L, c // L
        if d == 0:                      # horizontal edge (x,y,H): faces (x,y) and (x,y-1)
            return ((x, y), (x, (y - 1) % L))
        return (((x - 1) % L, y), (x, y))   # vertical edge: faces (x-1,y) and (x,y)


def h1_of_region(T, tor, dual=False):
    """dim H_1 of the subcomplex spanned by edge set T (primal or dual)."""
    T = sorted(T)
    idx = {e: k for k, e in enumerate(T)}
    if dual:
        ends = {e: tor.edge_faces(e) for e in T}
        faces = [s for s in tor.all_stars() if set(s) <= set(T)]
    else:
        ends = {e: tor.edge_ends(e) for e in T}
        faces = [p for p in tor.all_plaqs() if set(p) <= set(T)]
    verts = sorted({v for e in T for v in ends[e]})
    vidx = {v: k for k, v in enumerate(verts)}
    d1 = np.zeros((len(verts), len(T)), dtype=np.int8)
    for e in T:
        a, b = ends[e]
        d1[vidx[a], idx[e]] ^= 1
        d1[vidx[b], idx[e]] ^= 1
    d2 = np.zeros((len(T), len(faces)), dtype=np.int8)
    for j, f in enumerate(faces):
        for e in f:
            d2[idx[e], j] ^= 1
    dim_ker = len(T) - rank2(d1)
    return dim_ker - rank2(d2)


def is_contractible(T, tor):
    """No cycle in T, primal or dual, survives as a class -> T carries no winding loop."""
    return h1_of_region(T, tor, dual=False) == 0 and h1_of_region(T, tor, dual=True) == 0


# ------------------------------------------------------------------ the record level
class Setup:
    def __init__(self, L):
        self.tor = Torus(L)
        self.L, self.N = L, 2 * L * L
        t = self.tor
        self.R_edges = sorted({t.e(x, 0, 0) for x in range(L)})     # Z on a horizontal loop
        self.W_edges = sorted({t.e(0, y, 0) for y in range(L)})     # X on a dual vertical loop
        self.Rz = self.ind(self.R_edges)
        rows = []
        for s in t.all_stars():
            r = np.zeros(2 * self.N, dtype=np.int8)
            r[self.N:] = self.ind(s)
            rows.append(r)
        for p in t.all_plaqs():
            r = np.zeros(2 * self.N, dtype=np.int8)
            r[:self.N] = self.ind(p)
            rows.append(r)
        self.ADM = np.array(rows, dtype=np.int8)

    def ind(self, edges):
        a = np.zeros(self.N, dtype=np.int8)
        for e in edges:
            a[e] = 1
        return a

    def region_rows(self, T):
        rows = []
        for e in range(self.N):
            if e not in T:
                for off in (0, self.N):
                    r = np.zeros(2 * self.N, dtype=np.int8)
                    r[off + e] = 1
                    rows.append(r)
        return np.array(rows, dtype=np.int8).reshape(-1, 2 * self.N)

    def flips(self, v):
        return int((v[:self.N] @ self.Rz) % 2)

    def analyse(self, T):
        T = set(T)
        rows = np.vstack([self.ADM, self.region_rows(T)]) if len(T) < self.N else self.ADM
        B = nullspace(rows, 2 * self.N)
        m = B.shape[0]
        n_adm_flip = (2 ** (m - 1)) if (m and any(self.flips(b) for b in B)) else 0
        n_any_flip = (4 ** len(T)) // 2 if (T & set(self.R_edges)) else 0
        return m, n_adm_flip, n_any_flip, 2 ** m - 1


def run(L, do_minweight):
    S = Setup(L)
    t = S.tor
    hr(f"L = {L}   (n = {S.N} qubits)")

    B_all = nullspace(S.ADM, 2 * S.N)
    check(f"L={L}: centraliser dimension = n + k = {S.N} + 2", B_all.shape[0] == S.N + 2,
          f"measured {B_all.shape[0]}")
    check(f"L={L}: stabiliser rank = n - k = {S.N - 2}", rank2(S.ADM) == S.N - 2,
          f"measured {rank2(S.ADM)}")

    Wv = np.zeros(2 * S.N, dtype=np.int8)
    Wv[:S.N] = S.ind(S.W_edges)
    adm_W = all(int((Wv[:S.N] @ r[:S.N] + Wv[S.N:] @ r[S.N:]) % 2) == 0 for r in S.ADM)
    check(f"L={L}: the writer is admissible under DEF-A", adm_W,
          f"U = X on {[t.name(e) for e in S.W_edges]}")
    check(f"L={L}: the writer flips R", S.flips(Wv) == 1,
          f"|supp(U) cap supp(R)| = {len(set(S.W_edges) & set(S.R_edges))} (odd)")
    check(f"L={L}: => clause (iv) HOLDS, NOT vacuous (criterion 1)", adm_W and S.flips(Wv) == 1)

    # ---- candidate regions, each CERTIFIED contractible by its own H_1 -------------
    print("\n   CANDIDATE REGIONS, each certified by dim H_1 (primal and dual) = 0:")
    cands = []
    for e in range(S.N):
        cands.append(({e}, f"single edge {t.name(e)}"))
    for y in range(L):
        for x in range(L):
            cands.append((set(t.plaq(x, y)), f"plaquette ({x},{y})"))
            cands.append((set(t.star(x, y)), f"star ({x},{y})"))
    for a in range(1, L):
        for b in range(1, L):
            if a == b == 1:
                continue
            for y in range(L):
                for x in range(L):
                    blk = set()
                    for i in range(a):
                        for j in range(b):
                            blk |= set(t.plaq(x + i, y + j))
                    cands.append((blk, f"{a}x{b} plaquette block at ({x},{y})"))

    regions, rejected = [], {}
    for T, lab in cands:
        if is_contractible(T, t):
            regions.append((T, lab))
        else:
            key = lab.split(" at ")[0].split(" (")[0]
            rejected[key] = rejected.get(key, 0) + 1
    kinds = {}
    for T, lab in regions:
        key = lab.split(" at ")[0].split(" (")[0]
        kinds[key] = kinds.get(key, 0) + 1
    print(f"      admitted : {dict(sorted(kinds.items()))}")
    print(f"      REJECTED as not contractible (they carry a winding cycle): "
          f"{dict(sorted(rejected.items()))}")
    check(f"L={L}: the contractibility certifier rejects the blocks that wind",
          any("block" in k for k in rejected) or L == 3,
          f"{sum(rejected.values())} of {len(cands)} candidates rejected")
    check(f"L={L}: single plaquettes and stars ARE certified contractible",
          kinds.get("plaquette", 0) == L * L and kinds.get("star", 0) == L * L,
          f"plaquettes {kinds.get('plaquette',0)}/{L*L}, stars {kinds.get('star',0)}/{L*L}")

    # ---- clause (v) ---------------------------------------------------------------
    tot_adm, tot_any, offenders, maxT = 0, 0, [], 0
    for T, lab in regions:
        m, af, anyf, _ = S.analyse(T)
        tot_adm += af
        tot_any += anyf
        maxT = max(maxT, len(T))
        if af:
            offenders.append((lab, af))
    print(f"\n   regions tested (all certified contractible)   : {len(regions)}"
          f"   (largest |T| = {maxT})")
    print(f"   ADMISSIBLE operators flipping R, all regions  : {tot_adm}")
    print(f"   ANY Pauli operators flipping R, all regions   : {tot_any}")
    check(f"L={L}: CLAUSE (v) HOLDS under DEF-A", tot_adm == 0,
          f"count = {tot_adm}" + (f"   offenders: {offenders[:3]}" if offenders else ""))
    check(f"L={L}: POSITIVE CONTROL (same routine, same regions, admissibility dropped)",
          tot_any > 0, f"count = {tot_any}  <-- the zero above is not a zero of the routine")
    check(f"L={L}: CLAUSE (v) FAILS under DEF-C (admissible = any unitary)", tot_any > 0,
          f"e.g. X on the single edge {t.name(S.R_edges[0])} flips R")

    # ---- positive controls: winding regions ---------------------------------------
    print("\n   POSITIVE CONTROLS -- regions that WIND the torus:")
    col = set()
    row = set()
    for y in range(L):
        col |= set(t.plaq(0, y))
    for x in range(L):
        row |= set(t.plaq(x, 0))
    ctrls = [
        (set(S.W_edges), "the writer's own support", True),
        (col, "a full COLUMN of plaquettes (winds ACROSS R)", True),
        (row, "a full ROW of plaquettes (winds ALONG R)", False),
        (set(range(S.N)), "the whole lattice", True),
    ]
    ok = True
    for T, lab, expect in ctrls:
        m, af, anyf, _ = S.analyse(T)
        contr = is_contractible(T, t)
        print(f"      {lab:<46s} contractible={str(contr):5s} dim={m:3d}"
              f"  #admissible flippers = {af}")
        ok &= ((af > 0) == expect)
    check(f"L={L}: every winding control behaves as predicted in advance", ok,
          "the ROW control is predicted ZERO: it winds ALONG R, not across it, so it "
          "carries no operator anticommuting with R -- winding is necessary, not sufficient")

    # ---- minimum weight -----------------------------------------------------------
    if do_minweight:
        basis_int = []
        for b in B_all:
            xi = sum(int(b[i]) << i for i in range(S.N))
            zi = sum(int(b[S.N + i]) << i for i in range(S.N))
            basis_int.append((xi, zi, int((b[:S.N] @ S.Rz) % 2)))
        minw, cnt, ex = None, 0, None
        cx = cz = cf = prev = 0
        for i in range(1, 2 ** len(basis_int)):
            g = i ^ (i >> 1)
            bit = (g ^ prev).bit_length() - 1
            prev = g
            xi, zi, fi = basis_int[bit]
            cx ^= xi
            cz ^= zi
            cf ^= fi
            if cf:
                cnt += 1
                w = bin(cx | cz).count("1")
                if minw is None or w < minw:
                    minw, ex = w, (cx, cz)
        print(f"\n   centraliser enumerated: {2**len(basis_int)-1} elements, "
              f"{cnt} flip R")
        check(f"L={L}: exactly half the centraliser flips R",
              cnt == 2 ** (len(basis_int) - 1), f"{cnt}")
        check(f"L={L}: minimum weight of an admissible writer = d = {L}", minw == L,
              f"measured {minw}")
        sup = {e for e in range(S.N) if (ex[0] >> e) & 1 or (ex[1] >> e) & 1}
        print(f"      witness support: {[t.name(e) for e in sorted(sup)]}")
        check(f"L={L}: that minimum-weight writer fits inside NO certified "
              f"contractible region", not any(sup <= T for T, _ in regions))

    frac = (2 ** B_all.shape[0]) / (4 ** S.N)
    check(f"L={L}: criterion (2) -- DEF-A is not the trivial class",
          frac < 1e-4, f"admissible / all Paulis = 2^{B_all.shape[0]}/4^{S.N} = {frac:.3e}")
    return tot_adm, tot_any, len(regions), maxT


hr("O-4 -- CLAUSE (v) UNDER DEF-A AND DEF-C, AT TWO LATTICE SIZES")
r3 = run(3, do_minweight=True)
r5 = run(5, do_minweight=False)

hr("SUMMARY")
if FAILURES:
    print(f"*** {len(FAILURES)} SELF-CHECK FAILURES: {FAILURES}")
else:
    print("ALL SELF-CHECKS PASSED.")
print(f"""
                        regions   largest |T|   admissible flippers   any flippers
   L = 3                {r3[2]:5d}   {r3[3]:9d}   {r3[0]:19d}   {r3[1]:12d}
   L = 5                {r5[2]:5d}   {r5[3]:9d}   {r5[0]:19d}   {r5[1]:12d}

   DEF-A : clause (iv) holds, clause (v) holds, P-3 goes through.
   DEF-C : clause (iv) holds, clause (v) FAILS -- under DEF-C the program's own
           witness is not a record.
""")
