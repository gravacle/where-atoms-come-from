"""O-50-A  step 4.  CLAUSE (v) ON THE TORUS, through homology, with NO convention (D-23).

   A region is CONTRACTIBLE iff it is the induced edge set of an a x b block of vertices with
   a,b <= L-1, i.e. a genuine disk in the manifold -- this is fixed by the topology of T^2, not
   by a diameter convention.  For every such disk we ask, by exact F_2 linear algebra, whether
   ANY Pauli supported inside it is a non-trivial logical.  The answer must be NO for clause (v).

   POSITIVE CONTROL IN THE SAME TABLE (D-15): the same test on a NON-contractible region -- an
   annular band that wraps the torus -- which MUST come out YES, or the instrument is blind."""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from f2lib import Toric, sp, rank, in_span, nullspace, span

def logicals_in_region(T, edges):
    """basis of {Paulis supported in `edges`} cap N(S); returns (dim, n_outside_S, witness)."""
    n = T.n; E = sorted(edges)
    coords = [e for e in E] + [n + e for e in E]          # x_e then z_e for e in E
    rows = []
    for s in T.stab:
        rows.append([sp_basis(T, c, s) for c in coords])
    ns = nullspace(rows, len(coords)) if rows else []
    out, witness = 0, None
    for b in ns:
        v = [0] * (2 * n)
        for c, bit in zip(coords, b):
            if bit: v[c] = 1
        if not in_span(v, T.stab, 2 * n):
            out += 1
            if witness is None: witness = v
    # count the whole subspace, not just basis vectors
    full = []
    for b in ns:
        v = [0] * (2 * n)
        for c, bit in zip(coords, b):
            if bit: v[c] = 1
        full.append(v)
    nlog = 0; wit = None
    if len(full) <= 22:
        for v in span(full, 2 * n):
            if any(v) and not in_span(v, T.stab, 2 * n):
                nlog += 1
                if wit is None: wit = v
    else:
        nlog = -1
        for v in full:
            if not in_span(v, T.stab, 2 * n):
                wit = v; nlog = -1; break
    return len(ns), nlog, wit

def sp_basis(T, coord, s):
    """symplectic product of the unit vector e_coord with stabiliser s."""
    n = T.n
    return (s[n + coord] if coord < n else s[coord - n]) % 2

def disks(T):
    L = T.L
    for a in range(1, L):
        for b in range(1, L):
            for i0 in range(L):
                for j0 in range(L):
                    rows = [(i0 + t) % L for t in range(a)]
                    cols = [(j0 + t) % L for t in range(b)]
                    E = set()
                    for i in rows:
                        for t in range(b - 1): E.add(T.h(i, cols[t]))
                    for j in cols:
                        for t in range(a - 1): E.add(T.v(rows[t], j))
                    yield (a, b, i0, j0), E

def bands(T):
    """NON-contractible control regions: a x L bands that wrap the torus."""
    L = T.L
    for a in range(1, L):
        for i0 in range(L):
            rows = [(i0 + t) % L for t in range(a)]
            E = set()
            for i in rows:
                for j in range(L): E.add(T.h(i, j))
            for j in range(L):
                for t in range(a - 1): E.add(T.v(rows[t], j))
            yield (a, L, i0, 0), E

for L in (2, 3, 4):
    T = Toric(L)
    print("=" * 78)
    print(f"L = {L}   n = {T.n}   distance d = L (exact, s1)")
    worst = 0; nd = 0; bad = []
    for tag, E in disks(T):
        nd += 1
        dim, nlog, wit = logicals_in_region(T, E)
        worst = max(worst, len(E))
        if nlog != 0: bad.append((tag, len(E), nlog, wit))
    print(f"  CONTRACTIBLE regions (disks) tested: {nd}   largest disk: {worst} qubits")
    print(f"  disks containing a NON-TRIVIAL logical: {len(bad)}   -> clause (v) holds: {len(bad)==0}")
    print("  CONTROL, same test, NON-contractible bands:")
    for tag, E in bands(T):
        dim, nlog, wit = logicals_in_region(T, E)
        supp = [i for i in range(T.n) if wit[i] or wit[T.n + i]] if wit else None
        cls = None
        if wit:
            cls = dict(x_class=T.x_class(wit[:T.n]) if any(wit[:T.n]) else None,
                       z_class=T.z_class(wit[T.n:]) if any(wit[T.n:]) else None)
        print(f"    band {tag}: {len(E):3d} qubits, N(S)-dim inside = {dim:2d}, "
              f"non-trivial logicals inside = {nlog}, witness support = {supp}, class = {cls}")
