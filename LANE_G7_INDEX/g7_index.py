"""G-7: is the capacity law an INDEX? Test dim H1 = dim H0 + dim H2 - chi over GF(2)
by computing every rank explicitly, on the carriers the program already used."""
import numpy as np

def rank2(M):
    if M.size == 0: return 0
    A = M.copy() % 2; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]: A[i] ^= A[r]
        r += 1
        if r == rows: break
    return r

def homology(V, E, F):
    """V = #vertices; E = list of (a,b); F = list of edge-index lists."""
    nV, nE, nF = V, len(E), len(F)
    d1 = np.zeros((nV, nE), dtype=np.int8)          # edges -> vertices
    for k, (a, b) in enumerate(E): d1[a, k] ^= 1; d1[b, k] ^= 1
    d2 = np.zeros((nE, nF), dtype=np.int8)          # faces -> edges
    for k, f in enumerate(F):
        for e in f: d2[e, k] ^= 1
    r1, r2 = rank2(d1), rank2(d2)
    assert rank2((d1 @ d2) % 2) == 0, "d1 d2 != 0 — not a chain complex"
    h0 = nV - r1
    h1 = (nE - r1) - r2
    h2 = nF - r2
    chi = nV - nE + nF
    return dict(V=nV, E=nE, F=nF, chi=chi, h0=h0, h1=h1, h2=h2)

def torus(nx, ny):
    vid = {(i, j): j * nx + i for j in range(ny) for i in range(nx)}
    E, ind = [], {}
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = len(E); E.append((vid[(i, j)], vid[((i + 1) % nx, j)]))
            ind[('v', i, j)] = len(E); E.append((vid[(i, j)], vid[(i, (j + 1) % ny)]))
    F = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)],
          ind[('h', i, (j + 1) % ny)], ind[('v', i, j)]]
         for j in range(ny) for i in range(nx)]
    return nx * ny, E, F

def disk(n):
    vid = {(i, j): j * n + i for j in range(n) for i in range(n)}
    E, ind = [], {}
    for j in range(n):
        for i in range(n - 1): ind[('h', i, j)] = len(E); E.append((vid[(i, j)], vid[(i + 1, j)]))
    for j in range(n - 1):
        for i in range(n): ind[('v', i, j)] = len(E); E.append((vid[(i, j)], vid[(i, j + 1)]))
    F = [[ind[('h', i, j)], ind[('v', i + 1, j)], ind[('h', i, j + 1)], ind[('v', i, j)]]
         for j in range(n - 1) for i in range(n - 1)]
    return n * n, E, F

def tetra():
    E = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    ei = {frozenset(e): k for k, e in enumerate(E)}
    F = [[ei[frozenset((a,b))] for a,b in ((t[0],t[1]),(t[0],t[2]),(t[1],t[2]))]
         for t in ((0,1,2),(0,1,3),(0,2,3),(1,2,3))]
    return 4, E, F

CARRIERS = [("torus 2x2", *torus(2,2)), ("torus 2x3", *torus(2,3)),
            ("torus 3x3", *torus(3,3)), ("tetrahedron (sphere)", *tetra()),
            ("disk 3x3", *disk(3))]

print(f"{'carrier':<22}{'V':>3}{'E':>4}{'F':>4}{'chi':>5}{'H0':>4}{'H1':>4}{'H2':>4}"
      f"{'H0+H2-chi':>11}{'2-chi':>7}{'2^H1':>6}  index law   surface law")
print("-"*106)
ok_index = ok_surface = True
for name, V, E, F in CARRIERS:
    r = homology(V, E, F)
    gen = r['h0'] + r['h2'] - r['chi']
    surf = 2 - r['chi']
    a = (gen == r['h1']); b = (surf == r['h1'])
    ok_index &= a; ok_surface &= b
    print(f"{name:<22}{r['V']:>3}{r['E']:>4}{r['F']:>4}{r['chi']:>5}{r['h0']:>4}{r['h1']:>4}"
          f"{r['h2']:>4}{gen:>11}{surf:>7}{2**r['h1']:>6}  {'HOLDS ' if a else 'FAILS ':<11}"
          f"{'HOLDS' if b else 'FAILS'}")
print("-"*106)
print(f"SELF-CHECK  d1.d2 = 0 on every carrier: PASS (asserted in homology())")
print(f"SELF-CHECK  torus H1 = 2 at three different areas: "
      f"{'PASS' if all(homology(*torus(*s)[0:1]+torus(*s)[1:])['h1']==2 for s in ((2,2),(2,3),(3,3))) else 'FAIL'}")
print()
print(f"INDEX LAW   dim H1 = dim H0 + dim H2 - chi : {'HOLDS ON ALL 5 CARRIERS' if ok_index else 'FAILS'}")
print(f"SURFACE LAW dim H1 = 2 - chi              : {'holds on all 5' if ok_surface else 'FAILS on at least one'}")
