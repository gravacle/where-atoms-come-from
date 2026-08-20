"""T42-A shared library -- exact F_2 / Z_4 machinery for EARNING DISTANCE IN THE CORNER.

Everything here is exact: bitmask integers over F_2, small integers over Z_4.
No floats anywhere.  Every PASS/FAIL the callers print is gated by a computed boolean.
"""
import sys

MODEL = "/Users/bgm/MB Work/where-atoms-come-from/model"
if MODEL not in sys.path:
    sys.path.insert(0, MODEL)


def pc(x):
    """popcount of a python int (py3.9: no int.bit_count)."""
    return bin(x).count("1")


# ---------------------------------------------------------------- F_2 linear algebra on bitmasks
def solve_affine_f2(rows, rhs, nvars):
    """Solve M x = b over F_2.  rows: list of int masks (each a row of M, bit j = var j),
       rhs: list of 0/1.  Returns (x0, nullbasis, rank) with x0 a particular solution mask and
       nullbasis a list of masks spanning the solution space of Mx=0; or (None, None, rank) if
       inconsistent."""
    aug = []  # (rowmask, rhsbit)
    for r, b in zip(rows, rhs):
        aug.append([r, b & 1])
    pivots = []  # (col, rowindex-in-reduced)
    red = []
    for r, b in aug:
        for pcol, prow, pb in red:
            if (r >> pcol) & 1:
                r ^= prow
                b ^= pb
        if r == 0:
            if b:
                return None, None, len(red)
            continue
        col = r.bit_length() - 1  # leading (highest) set bit as pivot
        # eliminate this col from existing reduced rows
        nred = []
        for pcol2, prow2, pb2 in red:
            if (prow2 >> col) & 1:
                prow2 ^= r
                pb2 ^= b
            nred.append((pcol2, prow2, pb2))
        red = nred
        red.append((col, r, b))
    rank = len(red)
    # particular solution: set free vars 0, pivot vars from rhs
    x0 = 0
    for col, r, b in red:
        if b:
            x0 |= (1 << col)
    # check (defensive, exact)
    for r, b in zip(rows, rhs):
        assert pc(r & x0) % 2 == (b & 1)
    pivcols = set(c for c, _, _ in red)
    nullbasis = []
    for f in range(nvars):
        if f in pivcols:
            continue
        v = 1 << f
        for col, r, b in red:
            # value of pivot var 'col' so that row is satisfied with free var f = 1, others 0
            if (r >> f) & 1:
                v |= (1 << col)
        nullbasis.append(v)
    for v in nullbasis:
        for r in rows:
            assert pc(r & v) % 2 == 0
    return x0, nullbasis, rank


def span_all(basis):
    """All 2^k elements of the span (python ints).  Use only for k <= ~22."""
    out = [0]
    for g in basis:
        out += [x ^ g for x in out]
    return out


def span_numpy(basis):
    """All elements of the span as a numpy uint64 array (needs all masks < 2^64)."""
    import numpy as np
    arr = np.zeros(1, dtype=np.uint64)
    for g in basis:
        arr = np.concatenate([arr, arr ^ np.uint64(g)])
    return arr


def np_weights(arr):
    import numpy as np
    return np.bitwise_count(arr)


# ---------------------------------------------------------------- symplectic form on (x|z) pairs
def sp_pair(u, v, n):
    """Symplectic F_2 pairing of u=(xu|zu), v=(xv|zv) given as single 2n-bit masks
       (low n bits = x part, high n bits = z part)."""
    mask = (1 << n) - 1
    xu, zu = u & mask, u >> n
    xv, zv = v & mask, v >> n
    return (pc(xu & zv) + pc(zu & xv)) % 2


def weight_xz(u, n):
    """Number of qubits acted on: |supp(x) union supp(z)|."""
    mask = (1 << n) - 1
    return pc((u & mask) | (u >> n))


def crossings(u, v, n):
    """Number of sites where the two Paulis locally anticommute."""
    mask = (1 << n) - 1
    xu, zu = u & mask, u >> n
    xv, zv = v & mask, v >> n
    return pc((xu & zv) ^ (zu & xv))


def vec_to_mask(vec):
    """list of 0/1 -> int mask (bit i = vec[i])."""
    m = 0
    for i, b in enumerate(vec):
        if b:
            m |= (1 << i)
    return m


# ---------------------------------------------------------------- metric axioms, exact, ALL tuples
def metric_axioms(points, d):
    """d: dict[(p,q)] -> exact nonneg int.  Checks over ALL pairs and ALL triples.
       Returns dict of booleans + witnesses.  NOTHING is asserted; caller gates prints on these."""
    rep = {}
    idw = [p for p in points if d[(p, p)] != 0]
    rep["identity_ok"] = (len(idw) == 0)
    rep["identity_witness"] = idw[:3]
    symw = [(p, q) for p in points for q in points if d[(p, q)] != d[(q, p)]]
    rep["symmetry_ok"] = (len(symw) == 0)
    rep["symmetry_witness"] = symw[:3]
    indw = [(p, q) for p in points for q in points if p != q and d[(p, q)] == 0]
    rep["indiscernible_ok"] = (len(indw) == 0)
    rep["indiscernible_witness"] = indw[:3]
    triw = []
    for p in points:
        for q in points:
            dpq = d[(p, q)]
            for r in points:
                if d[(p, r)] > dpq + d[(q, r)]:
                    triw.append((p, q, r, d[(p, r)], dpq, d[(q, r)]))
                    if len(triw) >= 3:
                        rep["triangle_ok"] = False
                        rep["triangle_witness"] = triw
                        rep["n_pairs"] = len(points) ** 2
                        rep["n_triples"] = len(points) ** 3
                        rep["all"] = False
                        return rep
    rep["triangle_ok"] = (len(triw) == 0)
    rep["triangle_witness"] = triw[:3]
    rep["n_pairs"] = len(points) ** 2
    rep["n_triples"] = len(points) ** 3
    rep["all"] = (rep["identity_ok"] and rep["symmetry_ok"]
                  and rep["indiscernible_ok"] and rep["triangle_ok"])
    return rep


# ---------------------------------------------------------------- automorphisms of a finite metric
def aut_count(points, d, cap=10**7):
    """EXACT number of distance-preserving permutations, by backtracking.
       If the metric is uniform (all off-diagonal distances equal) the group is the full
       symmetric group; return ('uniform', n!) computed, not enumerated.
       Returns (kind, count) with kind in {'uniform','enumerated','capped'}."""
    pts = list(points)
    N = len(pts)
    vals = set(d[(p, q)] for p in pts for q in pts if p != q)
    if len(vals) <= 1:
        f = 1
        for i in range(2, N + 1):
            f *= i
        return "uniform", f
    idx = {p: i for i, p in enumerate(pts)}
    D = [[d[(p, q)] for q in pts] for p in pts]
    # distance profile per point (sorted multiset) for pruning
    prof = [tuple(sorted(D[i])) for i in range(N)]
    count = [0]
    capped = [False]

    def bt(k, perm):
        if count[0] > cap:
            capped[0] = True
            return
        if k == N:
            count[0] += 1
            return
        for cand in range(N):
            if cand in perm[:k]:
                continue
            if prof[cand] != prof[k]:
                continue
            ok = True
            for j in range(k):
                if D[k][j] != D[cand][perm[j]]:
                    ok = False
                    break
            if ok:
                perm[k] = cand
                bt(k + 1, perm)
        return

    bt(0, [0] * N)
    return ("capped" if capped[0] else "enumerated"), count[0]


# ---------------------------------------------------------------- carriers
def toric_stabilizers(L):
    """Toric code on the L x L torus.  Qubits on edges, n = 2 L^2.
       h(i,j) = edge from vertex (i,j) to (i,j+1), index i*L+j.
       v(i,j) = edge from vertex (i,j) to (i+1,j), index L^2 + i*L + j.
       Returns (n, stab_vecs list-of-2n-lists, star_masks, plaq_masks) with stars X-type
       and plaquettes Z-type.  One star and one plaquette are dependent; ALL are returned,
       rank handled by the F_2 solver."""
    n = 2 * L * L

    def h(i, j):
        return (i % L) * L + (j % L)

    def v(i, j):
        return L * L + (i % L) * L + (j % L)

    stars, plaqs = [], []
    for i in range(L):
        for j in range(L):
            stars.append(vec_to_mask([1 if e in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j))
                                      else 0 for e in range(n)]))
            plaqs.append(vec_to_mask([1 if e in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1))
                                      else 0 for e in range(n)]))
    stab = []
    for s in stars:
        stab.append([(s >> e) & 1 for e in range(n)] + [0] * n)      # X type
    for p in plaqs:
        stab.append([0] * n + [(p >> e) & 1 for e in range(n)])      # Z type
    return n, stab, stars, plaqs


def toric_edge_masks(L):
    """(horizontal-edge mask, vertical-edge mask, list of L vertical-cut h-edge masks,
        list of L horizontal-cut v-edge masks)."""
    n = 2 * L * L
    hmask = (1 << (L * L)) - 1
    vmask = ((1 << n) - 1) ^ hmask
    # The class-invariant transverse objects are CLOSED LOOPS (elements of ker d), because
    # <x + coboundary, S> = <x, S> iff dS = 0.  The L parallel horizontal loops r_i (all
    # h-edges of row i) are pairwise disjoint, jointly cover every h-edge, and each is a
    # cycle; likewise the L vertical loops c_j.  A writer crossing every r_i oddly carries
    # >= L horizontal edges -- integer arithmetic, no sampling.
    vcuts = []   # r_i: h-edges h(i,j) for all j  (fixed row i)
    for i in range(L):
        m = 0
        for j in range(L):
            m |= 1 << (i * L + j)
        vcuts.append(m)
    hcuts = []   # c_j: v-edges v(i,j) for all i  (fixed column j)
    for j in range(L):
        m = 0
        for i in range(L):
            m |= 1 << (L * L + i * L + j)
        hcuts.append(m)
    return hmask, vmask, vcuts, hcuts
