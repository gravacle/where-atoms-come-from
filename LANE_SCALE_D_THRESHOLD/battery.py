"""The scalar order-parameter battery.  Every quantity is EXACT (F_2 or integer arithmetic).

BASIS-DEPENDENT vs BASIS-INVARIANT is marked on every function and is not decoration: the
family of records is fixed only up to a symplectic change of basis of the logical group, so a
"threshold" visible only in a basis-dependent quantity is a property of the CHOSEN GENERATORS,
not of the carrier.  s2 reports both and re-runs the basis-dependent ones over a randomised
symplectic-basis ensemble.

THE CORE OBJECT.  For a region Reg, the ADMISSIBLE operations supported on it are the Paulis
supported in Reg that lie in N(S) (those and only those commute with H).  Their action on the
record family is the linear map  P -> ( sp(P,R_1), ..., sp(P,R_k) )  in F_2^k.  The IMAGE V of
that map is a subspace of dimension at most 2|Reg|, so it can be ENUMERATED EXACTLY -- which
makes "which records can one region flip together" an exact computation, not a search.
"""
import numpy as np
import carriers as C
from f2 import rank, kernel, span

def _cols(car, region):
    n = car["n"]; cols = []
    for q in sorted(region):
        e = [0] * (2 * n); e[q] = 1; cols.append(e)
        e = [0] * (2 * n); e[n + q] = 1; cols.append(e)
    return cols

def region_image(car, region, fam):
    """V = image in F_2^k of the admissible operations supported on `region`.  Returns the
       list of its elements (exact enumeration; dim V <= 2|region|)."""
    n = car["n"]; k = len(fam)
    cols = _cols(car, region)
    M = [[C.sp(c, s, n) for c in cols] for s in car["stabs"]]
    ker = kernel(M, len(cols))
    vs = []
    for c in ker:
        P = [0] * (2 * n)
        for t, bit in enumerate(c):
            if bit: P = [(a + b) % 2 for a, b in zip(P, cols[t])]
        vs.append([C.sp(P, R, n) for R in fam])
    return span(vs, k)

def local_logical_dim(car, region):
    """BASIS-INVARIANT.  dim_F2 of the group of logical CLASSES a single region can write."""
    n = car["n"]; cols = _cols(car, region)
    M = [[C.sp(c, s, n) for c in cols] for s in car["stabs"]]
    nul = len(cols) - rank(M, len(cols))
    # stabilisers supported INSIDE the region, restricted to the region's own coordinates --
    # same rank, but on a 2|region|-column matrix instead of a 2n-column one.
    reg = set(region); idx = sorted(reg)
    inside = [[s[q] for q in idx] + [s[n + q] for q in idx]
              for s in car["stabs"] if C.support(s, n) <= reg]
    ns = rank(inside, 2 * len(idx))
    return nul - ns

def regions(car, r, contiguous=True):
    n = car["n"]
    if contiguous:
        return [set(range(s, s + r)) for s in range(n - r + 1)]
    from itertools import combinations
    return [set(c) for c in combinations(range(n), r)]

# ------------------------------------------------------- graph scalars
def graph_scalars(adj, k, tag):
    A = np.array(adj, dtype=float)
    if k: np.fill_diagonal(A, 0.0)
    deg = A.sum(1) if k else np.zeros(0)
    E = int(A.sum() // 2)
    npairs = k * (k - 1) // 2
    seen = [False] * k; comps = []
    for s in range(k):
        if seen[s]: continue
        st = [s]; seen[s] = True; c = 0
        while st:
            u = st.pop(); c += 1
            for v in range(k):
                if A[u, v] and not seen[v]:
                    seen[v] = True; st.append(v)
        comps.append(c)
    A3 = A @ A @ A
    tri = np.trace(A3) / 6.0
    trip = sum(d * (d - 1) / 2.0 for d in deg)
    clus = (3 * tri / trip) if trip > 0 else 0.0
    L = np.diag(deg) - A
    ev = np.sort(np.linalg.eigvalsh(L)) if k > 1 else np.zeros(2)
    gap = float(ev[1]) if k > 1 else 0.0
    aev = np.linalg.eigvalsh(A) if k else np.zeros(1)
    lam = float(aev.max()) if k else 0.0
    tr = float(deg.sum())
    return {tag + "_edgefrac": (E / npairs) if npairs else 0.0,
            tag + "_ncomp": len(comps),
            tag + "_giant": (max(comps) / k) if k else 0.0,
            tag + "_clus": float(clus),
            tag + "_lapgap": gap,
            tag + "_lammax_over_tr": (lam / tr) if tr > 0 else 0.0,
            tag + "_percolates": bool(max(comps) > k / 2.0) if k else False}

# ------------------------------------------------------- the relations
def overlap_adj(car, fam):
    """BASIS-DEPENDENT.  Records whose SUPPORTS intersect."""
    n = car["n"]; k = len(fam)
    sup = [C.support(v, n) for v in fam]
    return [[1 if (i != j and sup[i] & sup[j]) else 0 for j in range(k)] for i in range(k)]

def coflip(car, fam, r):
    """BASIS-DEPENDENT.  From the exact region images: the adjacency of records that ONE
       admissible region-r operation flips TOGETHER, the largest such set, and the largest
       number of records a single region can flip individually."""
    k = len(fam)
    A = [[0] * k for _ in range(k)]
    maxset = 0; maxsingle = 0
    for reg in regions(car, r):
        V = region_image(car, reg, fam)
        touched = set()
        for v in V:
            w = sum(v)
            if w > maxset: maxset = w
            if w:
                idx = [i for i in range(k) if v[i]]
                touched.update(idx)
                for a in range(len(idx)):
                    for b in range(a + 1, len(idx)):
                        A[idx[a]][idx[b]] = A[idx[b]][idx[a]] = 1
        maxsingle = max(maxsingle, len(touched))
    return A, maxset, maxsingle
