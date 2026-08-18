"""GF(2) symplectic stabiliser machinery.  Pauli on n qubits packed as one python int:
   bits 0..n-1   = X-part
   bits n..2n-1  = Z-part
No external deps."""

def popcount(v): return bin(v).count("1")

def rank_gf2(rows):
    piv = {}
    r = 0
    for row in rows:
        cur = row
        while cur:
            p = cur.bit_length() - 1
            if p in piv:
                cur ^= piv[p]
            else:
                piv[p] = cur; r += 1; break
    return r, piv

def reduce_by(piv, v):
    cur = v
    while cur:
        p = cur.bit_length() - 1
        if p in piv: cur ^= piv[p]
        else: return cur
    return 0

def in_span(piv, v): return reduce_by(piv, v) == 0

def xpart(p, n): return p & ((1 << n) - 1)
def zpart(p, n): return p >> n
def mk(x, z, n): return x | (z << n)

def symp(p, q, n):
    """symplectic product mod 2: x_p.z_q + x_q.z_p"""
    return (popcount(xpart(p,n) & zpart(q,n)) + popcount(xpart(q,n) & zpart(p,n))) & 1

def pweight(p, n):
    """number of qubits carrying a non-identity Pauli"""
    return popcount(xpart(p,n) | zpart(p,n))

def all_commute(gens, n):
    bad = 0
    for i in range(len(gens)):
        for j in range(i+1, len(gens)):
            bad += symp(gens[i], gens[j], n)
    return bad

def code_k(gens, n):
    r,_ = rank_gf2(gens)
    return n - r, r

def nullspace_basis(rows, ncols):
    """basis of {v in F2^ncols : v.row = 0 for all rows}, rows given as ints"""
    # gaussian elim with column tracking
    M = list(rows)
    piv_col = []
    pr = 0
    R = [m for m in M]
    lead = {}
    for c in range(ncols):
        pr_row = None
        for i in range(pr, len(R)):
            if (R[i] >> c) & 1: pr_row = i; break
        if pr_row is None: continue
        R[pr], R[pr_row] = R[pr_row], R[pr]
        for i in range(len(R)):
            if i != pr and ((R[i] >> c) & 1): R[i] ^= R[pr]
        lead[c] = pr; piv_col.append(c); pr += 1
    free = [c for c in range(ncols) if c not in lead]
    basis = []
    for f in free:
        v = 1 << f
        for c in piv_col:
            if (R[lead[c]] >> f) & 1: v |= (1 << c)
        basis.append(v)
    return basis

def normaliser_basis(gens, n):
    """basis of S-perp inside F2^{2n} under the symplectic form"""
    swapped = [mk(zpart(g,n), xpart(g,n), n) for g in gens]   # v.symp g == v . swap(g) as plain dot
    return nullspace_basis(swapped, 2*n)

def min_logical_weight(gens, n, trials=200000, seed=0, exhaustive_limit=1<<22):
    """upper bound on d: min pweight over S-perp \ S.  Exhaustive when the coset space is small."""
    _, piv_S = rank_gf2(gens)
    N = normaliser_basis(gens, n)
    if not N: return None, 0
    m = len(N)
    best = None
    if (1 << m) <= exhaustive_limit:
        for mask in range(1, 1 << m):
            v = 0; mm = mask; i = 0
            while mm:
                if mm & 1: v ^= N[i]
                mm >>= 1; i += 1
            if in_span(piv_S, v): continue
            w = pweight(v, n)
            if best is None or w < best: best = w
        return best, m
    import random
    rng = random.Random(seed)
    for _ in range(trials):
        v = 0
        for b in N:
            if rng.getrandbits(1): v ^= b
        if v == 0 or in_span(piv_S, v): continue
        # greedy local reduction
        improved = True
        w = pweight(v, n)
        while improved:
            improved = False
            for b in N:
                v2 = v ^ b
                if v2 == 0 or in_span(piv_S, v2): continue
                w2 = pweight(v2, n)
                if w2 < w: v, w = v2, w2; improved = True
        if best is None or w < best: best = w
    return best, m

def css_split(gens, n):
    """rank of the pure-X part, pure-Z part, and total.  CSS in this basis iff rx+rz == r."""
    r, _ = rank_gf2(gens)
    # pure-X elements of S: v in span(gens) with z-part 0
    # solve over coefficient space
    G = list(gens); g = len(G)
    # build matrix: coefficient c -> z-part of sum
    rowsZ = []
    for bit in range(n):
        row = 0
        for i, gg in enumerate(G):
            if (zpart(gg, n) >> bit) & 1: row |= (1 << i)
        rowsZ.append(row)
    coefX = nullspace_basis(rowsZ, g)
    rowsX = []
    for bit in range(n):
        row = 0
        for i, gg in enumerate(G):
            if (xpart(gg, n) >> bit) & 1: row |= (1 << i)
        rowsX.append(row)
    coefZ = nullspace_basis(rowsX, g)
    def realise(coefs):
        out = []
        for c in coefs:
            v = 0; cc = c; i = 0
            while cc:
                if cc & 1: v ^= G[i]
                cc >>= 1; i += 1
            out.append(v)
        return out
    SX = realise(coefX); SZ = realise(coefZ)
    rx, _ = rank_gf2(SX); rz, _ = rank_gf2(SZ)
    rtot, _ = rank_gf2(SX + SZ)
    return rx, rz, r, (rtot == r and rx + rz == r)
