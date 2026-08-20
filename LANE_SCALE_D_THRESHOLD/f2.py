"""Exact F_2 linear algebra used by the battery.  No sampling anywhere."""
def rref(rows, ncol):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

def rank(rows, ncol):
    if not rows: return 0
    return len(rref(rows, ncol)[0])

def kernel(rows, ncol):
    """Basis of {c : rows . c = 0} over F_2."""
    if not rows:
        return [[1 if j == i else 0 for j in range(ncol)] for i in range(ncol)]
    R, piv = rref(rows, ncol)
    free = [c for c in range(ncol) if c not in piv]
    out = []
    for f in free:
        v = [0] * ncol; v[f] = 1
        for i, c in enumerate(piv): v[c] = R[i][f]
        out.append(v)
    return out

def span(vs, k):
    """Every element of the span of vs in F_2^k, as tuples.  Only called with dim <= ~12."""
    R, _ = rref(vs, k) if vs else ([], [])
    out = [tuple([0] * k)]
    for b in R:
        out += [tuple((x + y) % 2 for x, y in zip(o, b)) for o in out]
    return out

def solvable(M, b):
    if not M: return all(x == 0 for x in b)
    ncol = len(M[0])
    return rank(M, ncol) == rank([row + [bb] for row, bb in zip(M, b)], ncol + 1)
