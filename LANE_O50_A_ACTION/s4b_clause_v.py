"""O-50-A  step 4b.  CLAUSE (v) ON THE TORUS through HOMOLOGY, with NO convention (D-23).

   CONTRACTIBLE REGION := the closed edge set of an a x b block of PLAQUETTES with
   a,b <= L-1.  Such a block does not wrap either cycle of T^2, so it IS a topological disk;
   nothing here refers to a diameter, to a fraction of the system, or to "the whole chain".
   a = L (or b = L) makes the block wrap: an ANNULUS, and that is the positive control.

   TEST, exact and enumeration-free.  V(Rg) := {Paulis supported in Rg} cap N(S).  Clause (v)
   holds for Rg iff V(Rg) subset S.  Both are F_2 subspaces, so the test is
        rank(S-basis) == rank(S-basis + V-basis),
   and the number of non-trivial logical classes available inside Rg is 2^(dimV - dim(V cap S)) - 1,
   computed by ranks alone.  No 2^d enumeration anywhere."""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from f2lib import Toric, sp, rank, in_span, nullspace

def region_analysis(T, edges):
    n = T.n; E = sorted(edges)
    coords = [e for e in E] + [n + e for e in E]
    rows = [[(s[n + c] if c < n else s[c - n]) % 2 for c in coords] for s in T.stab]
    ns = nullspace(rows, len(coords))
    V = []
    for b in ns:
        v = [0] * (2 * n)
        for c, bit in zip(coords, b):
            if bit: v[c] = 1
        V.append(v)
    rS = rank(T.stab, 2 * n)
    rSV = rank(T.stab + V, 2 * n)
    dimV = len(V)
    dim_quotient = rSV - rS                     # dim of image of V in N(S)/S
    n_logical_classes = 2 ** dim_quotient - 1
    wit = None
    for v in V:
        if not in_span(v, T.stab, 2 * n): wit = v; break
    return dict(qubits=len(E), dimV=dimV, dim_V_mod_S=dim_quotient,
                n_nontrivial_logical_classes=n_logical_classes,
                clause_v_holds=(dim_quotient == 0), witness=wit)

def plaquette_block(T, a, b, i0, j0):
    L = T.L; E = set()
    for i in range(a + 1):
        for j in range(b): E.add(T.h(i0 + i, j0 + j))
    for i in range(a):
        for j in range(b + 1): E.add(T.v(i0 + i, j0 + j))
    return E

for L in (2, 3, 4):
    T = Toric(L)
    print("=" * 78)
    print(f"L = {L}   n = {T.n}   code distance d = {L} (exact, s1)")
    print("  CONTRACTIBLE (disk) regions -- a x b plaquette blocks, a,b <= L-1")
    worst = None; nreg = 0; fails = 0
    for a in range(1, L):
        for b in range(1, L):
            biggest = None
            for i0 in range(L):
                for j0 in range(L):
                    r = region_analysis(T, plaquette_block(T, a, b, i0, j0))
                    nreg += 1
                    if not r['clause_v_holds']: fails += 1
                    biggest = r
            print(f"    {a}x{b} plaquettes: {biggest['qubits']:2d} qubits, dim V = {biggest['dimV']:2d}, "
                  f"dim V/S = {biggest['dim_V_mod_S']}, non-trivial logicals inside = "
                  f"{biggest['n_nontrivial_logical_classes']}, clause (v) holds = {biggest['clause_v_holds']}")
    print(f"  ALL {nreg} disk regions tested; failures = {fails}; "
          f"CLAUSE (v) HOLDS ON THE TORUS BY HOMOLOGY = {fails == 0}")
    print("  POSITIVE CONTROL -- NON-contractible (wrapping) blocks, a = L or b = L:")
    for (a, b) in [(1, L), (L, 1), (L - 1, L), (L, L)]:
        if a < 1 or b < 1: continue
        r = region_analysis(T, plaquette_block(T, a, b, 0, 0))
        w = r['witness']
        supp = [i for i in range(T.n) if w[i] or w[T.n + i]] if w else None
        cls = None
        if w: cls = (T.x_class(w[:T.n]) if any(w[:T.n]) else None,
                     T.z_class(w[T.n:]) if any(w[T.n:]) else None)
        print(f"    {a}x{b} (wraps): {r['qubits']:2d} qubits, dim V/S = {r['dim_V_mod_S']}, "
              f"non-trivial logicals inside = {r['n_nontrivial_logical_classes']}, "
              f"clause (v) holds = {r['clause_v_holds']}, witness class (x,z) = {cls}")
