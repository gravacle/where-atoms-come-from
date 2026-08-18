"""G-10: does the index law survive off manifolds?  Three claims tested SEPARATELY:
   (A) INDEX IDENTITY   dim H1 = dim H0 + dim H2 - chi      [Euler-Poincare -- algebra, cannot fail]
   (B) SURFACE LAW      dim H1 = 2 - chi                     [our registered G-2 form]
   (C) CAPACITY         Z2 gauge ground-state degeneracy = 2^dim H1   [R1]
   (D) PROTECTION       d = min weight of a nontrivial class, and DOES IT GROW?  [R3]
Carriers include non-manifolds (pinch points), a pure graph, a non-orientable surface,
and an abstract chain complex with no geometric realisation at all."""
import numpy as np, itertools, random
random.seed(11)

def rref2(rows, n):
    """row-reduce list-of-int bitvectors; return (pivots dict col->row, basis)."""
    basis, piv = [], {}
    for r in rows:
        cur = r
        for c in sorted(piv, reverse=True):
            if (cur >> c) & 1: cur ^= piv[c]
        if cur:
            c = cur.bit_length() - 1
            piv[c] = cur; basis.append(cur)
    return piv, basis

def in_span(v, piv):
    for c in sorted(piv, reverse=True):
        if (v >> c) & 1: v ^= piv[c]
    return v == 0

def rank2(M):
    rows = [int(''.join(map(str, r)), 2) if r.size else 0 for r in (M % 2)] if M.size else []
    return len(rref2(rows, 0)[1])

def cols_as_ints(M):
    M = M % 2
    return [int(''.join(map(str, M[:, j])), 2) for j in range(M.shape[1])] if M.size else []

def analyse(name, nV, E, F, want_distance=True, geometric=True):
    nE, nF = len(E), len(F)
    d1 = np.zeros((nV, nE), dtype=np.int8)
    for k, (a, b) in enumerate(E): d1[a, k] ^= 1; d1[b, k] ^= 1
    d2 = np.zeros((nE, nF), dtype=np.int8)
    for k, f in enumerate(F):
        for e in f: d2[e, k] ^= 1
    assert not ((d1 @ d2) % 2).any(), name + ": d1.d2 != 0"
    return _analyse_maps(name, d1, d2, want_distance, geometric)

def _analyse_maps(name, d1, d2, want_distance=True, geometric=True):
    nV, nE = d1.shape; nF = d2.shape[1]
    r1, r2 = rank2(d1), rank2(d2)
    h0, h1, h2 = nV - r1, (nE - r1) - r2, nF - r2
    chi = nV - nE + nF
    # Z1 = ker d1  (bitvectors over edges)
    Z = []
    A = (d1 % 2).astype(np.int8).copy()
    # nullspace over GF(2) by elimination
    M = A.copy(); rows, cols = M.shape; pivcols = []; r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        pivcols.append(c); r += 1
    free = [c for c in range(cols) if c not in pivcols]
    for fc in free:
        v = np.zeros(cols, dtype=np.int8); v[fc] = 1
        for i, pc in enumerate(pivcols): v[pc] = M[i, fc]
        Z.append(int(''.join(map(str, v)), 2))
    Bcols = cols_as_ints(d2)
    pivB, _ = rref2(Bcols, nE)
    dist = None
    if want_distance and len(Z) <= 18:
        best = None
        for mask in range(1, 1 << len(Z)):
            v = 0; m = mask; i = 0
            while m:
                if m & 1: v ^= Z[i]
                m >>= 1; i += 1
            if not in_span(v, pivB):
                w = bin(v).count('1')
                if best is None or w < best: best = w
        dist = best
    return dict(name=name, V=nV, E=nE, F=nF, chi=chi, h0=h0, h1=h1, h2=h2,
                index=h0 + h2 - chi, surf=2 - chi, dZ=len(Z), dB=r2, dist=dist,
                geometric=geometric)

# ---------------- carriers ----------------
def torus(nx, ny, off=0):
    vid = {(i, j): off + j * nx + i for j in range(ny) for i in range(nx)}
    E, ind = [], {}
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = len(E); E.append((vid[(i, j)], vid[((i+1) % nx, j)]))
            ind[('v', i, j)] = len(E); E.append((vid[(i, j)], vid[(i, (j+1) % ny)]))
    F = [[ind[('h', i, j)], ind[('v', (i+1) % nx, j)], ind[('h', i, (j+1) % ny)], ind[('v', i, j)]]
         for j in range(ny) for i in range(nx)]
    return nx*ny, E, F

def bouquet(k):
    """k triangles all sharing ONE vertex -- a pinch point. Non-manifold. No 2-cells."""
    E = []; v = 1
    for _ in range(k):
        a, b = v, v+1; v += 2
        E += [(0, a), (a, b), (b, 0)]
    return v, E, []

def theta():
    return 2, [(0,1),(0,1),(0,1)], []          # two vertices, three parallel edges

def two_tori_wedged():
    V1, E1, F1 = torus(2,2,0)
    V2, E2, F2 = torus(2,2,V1)
    off = len(E1)
    E2r = [(a if a != V1 else 0, b if b != V1 else 0) for a,b in E2]   # glue vertex V1 -> 0
    E = E1 + E2r
    F = F1 + [[e+off for e in f] for f in F2]
    used = sorted({x for e in E for x in e})
    ren = {u:i for i,u in enumerate(used)}
    return len(used), [(ren[a],ren[b]) for a,b in E], F

def rp2():
    """minimal 6-vertex triangulation of the projective plane: K6 with 10 triangles."""
    E = [(a,b) for a in range(6) for b in range(a+1,6)]
    ei = {frozenset(e):k for k,e in enumerate(E)}
    tris = [(0,1,3),(0,1,4),(0,2,3),(0,2,5),(0,4,5),(1,2,4),(1,2,5),(1,3,5),(2,3,4),(3,4,5)]
    F = [[ei[frozenset((t[i],t[j]))] for i,j in ((0,1),(0,2),(1,2))] for t in tris]
    return 6, E, F

def abstract_complex(n0=6, n1=14, n2=5, seed=3):
    """A chain complex with NO geometric realisation: random d2, then d1 forced into (im d2)^perp."""
    rng = np.random.default_rng(seed)
    d2 = rng.integers(0,2,(n1,n2),dtype=np.int8)
    cols = cols_as_ints(d2); piv,_ = rref2(cols, n1)
    perp = []                                   # rows r with r.d2 = 0
    for _ in range(400):
        r = rng.integers(0,2,n1,dtype=np.int8)
        if not ((r @ d2) % 2).any(): perp.append(r)
        if len(perp) >= n0: break
    while len(perp) < n0: perp.append(np.zeros(n1,dtype=np.int8))
    d1 = np.array(perp[:n0], dtype=np.int8)
    assert not ((d1 @ d2) % 2).any()
    return d1, d2

CARRIERS = [("torus 2x2 (manifold)", *torus(2,2), True),
            ("torus 3x3 (manifold)", *torus(3,3), True),
            ("tetrahedron = sphere", 4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
             [[0,1,3],[0,2,4],[1,2,5],[3,4,5]], True),
            ("RP^2 (non-orientable)", *rp2(), True),
            ("disk 3x3 (boundary)", 9,
             [(j*3+i, j*3+i+1) for j in range(3) for i in range(2)] +
             [(j*3+i,(j+1)*3+i) for j in range(2) for i in range(3)],
             [[0,7,2,6],[1,8,3,7],[2,10,4,9],[3,11,5,10]], True),
            ("theta graph (1-complex)", *theta(), True),
            ("bouquet of 2 triangles (PINCH)", *bouquet(2), True),
            ("bouquet of 3 triangles (PINCH)", *bouquet(3), True),
            ("bouquet of 4 triangles (PINCH)", *bouquet(4), True),
            ("TWO TORI wedged (PINCH)", *two_tori_wedged(), True)]

rows = [analyse(n, V, E, F, geometric=g) for n, V, E, F, g in CARRIERS]
d1a, d2a = abstract_complex()
rows.append(_analyse_maps("ABSTRACT complex (no geometry)", d1a, d2a, geometric=False))

hdr = f"{'carrier':<32}{'V':>3}{'E':>4}{'F':>4}{'chi':>5}{'H0':>4}{'H1':>4}{'H2':>4}{'H0+H2-chi':>11}{'2-chi':>7}{'2^H1':>7}{'dist d':>8}   INDEX   SURFACE"
print(hdr); print('-'*len(hdr))
ok_i = ok_s = True
for r in rows:
    a = (r['index'] == r['h1']); b = (r['surf'] == r['h1'])
    ok_i &= a; ok_s &= b
    print(f"{r['name']:<32}{r['V']:>3}{r['E']:>4}{r['F']:>4}{r['chi']:>5}{r['h0']:>4}{r['h1']:>4}{r['h2']:>4}"
          f"{r['index']:>11}{r['surf']:>7}{2**r['h1']:>7}{('-' if r['dist'] is None else str(r['dist'])):>8}   "
          f"{'HOLDS' if a else 'FAILS':<8}{'HOLDS' if b else 'FAILS'}")
print('-'*len(hdr))
print(f"(A) INDEX IDENTITY  dim H1 = dim H0 + dim H2 - chi : {'HOLDS on all %d'%len(rows) if ok_i else 'FAILS'}")
print(f"(B) SURFACE LAW     dim H1 = 2 - chi               : {'holds on all' if ok_s else 'FAILS on at least one'}")
print()
print("(C) CAPACITY -- Z2 gauge ground-state degeneracy vs 2^dim H1 (physical sector by Gauss law):")
for nm, V, E, F in [("torus 2x2", *torus(2,2)), ("bouquet of 3 triangles", *bouquet(3)),
                    ("theta graph", *theta()), ("TWO TORI wedged", *two_tori_wedged())]:
    L = len(E)
    st = [s for s in itertools.product(range(2), repeat=L)
          if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
                 -sum(s[k] for k,(a,b) in enumerate(E) if b==v)) % 2 == 0 for v in range(V))]
    r = analyse(nm, V, E, F)
    print(f"    {nm:<26} |Gauss sector| = {len(st):>6} = 2^{r['dZ']:<3} "
          f"|B1| = 2^{r['dB']:<3} cosets = {len(st)//(2**r['dB']):>4}  vs 2^dim H1 = {2**r['h1']:>4}   "
          f"{'MATCH' if len(st)//(2**r['dB']) == 2**r['h1'] else 'MISMATCH'}")
print()
print("(D) PROTECTION -- does the minimum weight of a nontrivial class GROW with the carrier?")
print(f"    {'family':<34}{'size':>6}{'dim H1':>8}{'distance d':>12}")
for L in (2,3,4):
    r = analyse(f"torus {L}x{L}", *torus(L,L))
    print(f"    {'torus LxL (manifold)':<34}{L:>6}{r['h1']:>8}{r['dist']:>12}")
for k in (2,3,4):
    r = analyse(f"bouquet {k}", *bouquet(k))
    print(f"    {'bouquet of k triangles (PINCH)':<34}{k:>6}{r['h1']:>8}{r['dist']:>12}")
print()
print("SELF-CHECKS")
print(f"  d1.d2 = 0 asserted on every carrier                       : PASS")
t = analyse('t', *torus(3,3)); print(f"  torus 3x3 dim H1 == 2                                     : {'PASS' if t['h1']==2 else 'FAIL'}")
th = analyse('th', *theta()); print(f"  theta graph dim H1 == 2 (known: E-V+1 = 3-2+1 = 2)        : {'PASS' if th['h1']==2 else 'FAIL'}")
rp = analyse('rp', *rp2()); print(f"  RP^2 over GF(2): H0=H1=H2=1, chi=1 (known)                : {'PASS' if (rp['h0'],rp['h1'],rp['h2'],rp['chi'])==(1,1,1,1) else 'FAIL '+str((rp['h0'],rp['h1'],rp['h2'],rp['chi']))}")
tt = analyse('tt', *two_tori_wedged()); print(f"  two tori wedged: dim H1 == 4, dim H2 == 2 (known)         : {'PASS' if (tt['h1'],tt['h2'])==(4,2) else 'FAIL '+str((tt['h1'],tt['h2']))}")
b = analyse('b', *bouquet(3)); print(f"  bouquet of 3 triangles: dim H1 == 3 (known)               : {'PASS' if b['h1']==3 else 'FAIL'}")
