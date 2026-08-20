"""O-50-A  step 1.  The toric-code carrier in the F_2 symplectic representation.
   Logicals COMPUTED by record_model.symplectic_logicals -- never nominated (D-18).
   Symplectic pairing matrix checked non-degenerate BEFORE anything else.
   Code distance, homology class of the logical support, and the carrier's automorphism
   group (D-22).  Every number here is EXACT integer / F_2 arithmetic."""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import symplectic_logicals
from f2lib import Toric, sp, rref, rank, in_span, nullspace, span

def det_f2(M):
    M = [r[:] for r in M]; N = len(M); d = 1
    for c in range(N):
        p = next((i for i in range(c, N) if M[i][c]), None)
        if p is None: return 0
        M[c], M[p] = M[p], M[c]
        for i in range(N):
            if i != c and M[i][c]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[c])]
    return 1

def wt(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

# ---------------------------------------------------------------- automorphisms (D-22)
def grid_auts(L):
    """Aut of the L x L torus grid graph on VERTICES, by exact backtracking. L>=3."""
    V = L * L
    idx = lambda i, j: (i % L) * L + (j % L)
    adj = [set() for _ in range(V)]
    for i in range(L):
        for j in range(L):
            for (a, b) in ((i, j + 1), (i, j - 1), (i + 1, j), (i - 1, j)):
                adj[idx(i, j)].add(idx(a, b))
    order = list(range(V))
    out = []
    def rec(k, sigma, used):
        if k == V:
            out.append(sigma[:]); return
        u = order[k]
        for c in range(V):
            if c in used: continue
            ok = True
            for w in order[:k]:
                if (w in adj[u]) != (sigma[w] in adj[c]): ok = False; break
            if ok:
                sigma[u] = c; used.add(c); rec(k + 1, sigma, used); used.discard(c); sigma[u] = -1
    rec(0, [-1] * V, set())
    return out, adj

def carrier_auts(T):
    """Exact order of the group of EDGE permutations preserving the stabiliser structure.
       L=2: brute force over all n! edge permutations.
       L>=3: every such permutation induces a vertex-graph automorphism and is determined by
             it (adjacent vertices share exactly one edge), so enumerate those and test."""
    L, n = T.L, T.n
    Aset = set(tuple(a[:n]) for a in T.A)          # x-supports of vertex ops
    Bset = set(tuple(b[n:]) for b in T.B)          # z-supports of plaquette ops
    if L == 2:
        cnt = 0
        for perm in itertools.permutations(range(n)):
            ok = True
            for s in Aset:
                t = [0] * n
                for e in range(n):
                    if s[e]: t[perm[e]] = 1
                if tuple(t) not in Aset: ok = False; break
            if not ok: continue
            for s in Bset:
                t = [0] * n
                for e in range(n):
                    if s[e]: t[perm[e]] = 1
                if tuple(t) not in Bset: ok = False; break
            if ok: cnt += 1
        return cnt, "exact brute force over %d! edge permutations" % n
    sig, adj = grid_auts(L)
    idx = lambda i, j: (i % L) * L + (j % L)
    # edge -> unordered vertex pair
    epair, pair2e = {}, {}
    for i in range(L):
        for j in range(L):
            e = T.h(i, j); p = frozenset((idx(i, j), idx(i, j + 1))); epair[e] = p; pair2e[p] = e
            e = T.v(i, j); p = frozenset((idx(i, j), idx(i + 1, j))); epair[e] = p; pair2e[p] = e
    cnt = 0
    for s in sig:
        perm = [None] * n; ok = True
        for e in range(n):
            p = epair[e]; q = frozenset(s[x] for x in p)
            if q not in pair2e: ok = False; break
            perm[e] = pair2e[q]
        if not ok or len(set(perm)) != n: continue
        good = True
        for b in Bset:
            t = [0] * n
            for e in range(n):
                if b[e]: t[perm[e]] = 1
            if tuple(t) not in Bset: good = False; break
        if good: cnt += 1
    return cnt, "exact: induced by the %d vertex-graph automorphisms, plaquette set checked" % len(sig)

# ---------------------------------------------------------------- main
def analyse(L, do_distance=True, do_aut=True):
    T = Toric(L); n = T.n
    rec = {}
    rec['L'] = L; rec['n_qubits'] = n; rec['hilbert_dim'] = 2 ** n
    r = rank(T.stab, 2 * n)
    rec['n_stabiliser_generators_supplied'] = len(T.stab)
    rec['stabiliser_rank'] = r
    rec['k_logical_qubits'] = n - r
    pairs = symplectic_logicals(T.stab, n)
    rec['n_conjugate_pairs_returned'] = len(pairs)
    flat = [v for pr in pairs for v in pr]
    G = [[sp(a, b, n) for b in flat] for a in flat]
    rec['symplectic_gram'] = G
    rec['gram_det_f2'] = det_f2(G)
    rec['gram_nondegenerate'] = (det_f2(G) == 1)
    # each returned vector must be in N(S) and not in S
    Srr, _ = rref(T.stab, 2 * n)
    rec['all_logicals_commute_with_stabilisers'] = all(sp(v, s, n) == 0 for v in flat for s in T.stab)
    rec['no_logical_lies_in_S'] = all(not in_span(v, T.stab, 2 * n) for v in flat)
    rec['logical_weights'] = [wt(v, n) for v in flat]
    rec['logical_types'] = ['X-type' if not any(v[n:]) else ('Z-type' if not any(v[:n]) else 'mixed')
                            for v in flat]
    rec['logical_homology'] = []
    for v in flat:
        xs, zs = v[:n], v[n:]
        rec['logical_homology'].append(dict(z_class=T.z_class(zs) if any(zs) else None,
                                            x_class=T.x_class(xs) if any(xs) else None))
    if do_distance:
        Zc = nullspace(T.vertex_rows(), n)          # Z-type: cycles
        Xc = nullspace(T.plaquette_rows(), n)       # X-type: cocycles
        rec['cycle_space_dim'] = len(Zc); rec['cocycle_space_dim'] = len(Xc)
        Bz = [b[n:] for b in T.B]; Ax = [a[:n] for a in T.A]
        dz = min((sum(v) for v in span(Zc, n) if any(v) and not in_span(v, Bz, n)), default=None)
        dx = min((sum(v) for v in span(Xc, n) if any(v) and not in_span(v, Ax, n)), default=None)
        rec['d_Z'] = dz; rec['d_X'] = dx; rec['distance'] = min(dz, dx)
        # a minimum-weight Z logical, and its homology class
        best = min((v for v in span(Zc, n) if any(v) and not in_span(v, Bz, n)), key=sum)
        rec['min_weight_Z_logical_support'] = [i for i in range(n) if best[i]]
        rec['min_weight_Z_logical_class'] = T.z_class(best)
        bestx = min((v for v in span(Xc, n) if any(v) and not in_span(v, Ax, n)), key=sum)
        rec['min_weight_X_logical_class'] = T.x_class(bestx)
        if L == 2:                                   # cross-check over the WHOLE of N(S)
            NS = nullspace([[sp([1 if k == j else 0 for k in range(2 * n)], s, n)
                             for j in range(2 * n)] for s in T.stab], 2 * n)
            rec['N(S)_dim'] = len(NS)
            rec['distance_full_NS_check'] = min(wt(v, n) for v in span(NS, 2 * n)
                                                if any(v) and not in_span(v, T.stab, 2 * n))
    if do_aut:
        a, how = carrier_auts(T)
        rec['aut_group_order_edge_permutations'] = a; rec['aut_method'] = how
        rec['sym_group_order_n_factorial'] = None
        import math; rec['n_factorial'] = math.factorial(n)
        rec['permutation_symmetric'] = (a == math.factorial(n))
    return rec, T, pairs

if __name__ == "__main__":
    out = []
    for L in (2, 3, 4):
        r, T, pairs = analyse(L, do_distance=True, do_aut=True)
        out.append(r)
    for r in out:
        print("=" * 78)
        for k, v in r.items():
            print(f"  {k:42s} {v}")
