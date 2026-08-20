"""O-50-A  step 4c.  CLAUSE (v) ON THE TORUS IS EXACTLY A HOMOLOGY STATEMENT (D-23).

   THE POINT.  Do not argue from region diameter at all.  Show instead that the map
        chi : N(S)/S  ->  H_1(T^2;F2) (+) H^1(T^2;F2) = F_2^4,
        chi(v) = ( z_class(v_z) , x_class(v_x) )
   is an ISOMORPHISM (verified: 4x4 F_2 matrix on the computed logical basis, det = 1), so an
   operator is a non-trivial logical IFF ITS SUPPORT CARRIES A NON-TRIVIAL HOMOLOGY CLASS OF
   THE TORUS.  A region contained in a topological disk has H_1 = 0 and therefore carries no
   such class -- clause (v) then holds for free, with no diameter convention, no "proper
   sub-region" rule, and no appeal to what counts as the whole system.

   VERIFIED NUMERICALLY, region by region:  dim( V(Rg) mod S )  ==  rank chi(V(Rg))
   for EVERY region tested, contractible and wrapping alike.  If those two ever disagreed the
   homological reading would be wrong; they never do."""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
from record_model import symplectic_logicals
from f2lib import Toric, sp, rank, in_span, nullspace

def chi(T, v):
    n = T.n
    return list(T.z_class(v[n:])) + list(T.x_class(v[:n]))

def region_V(T, edges):
    n = T.n; E = sorted(edges)
    coords = [e for e in E] + [n + e for e in E]
    rows = [[(s[n + c] if c < n else s[c - n]) % 2 for c in coords] for s in T.stab]
    V = []
    for b in nullspace(rows, len(coords)):
        v = [0] * (2 * n)
        for c, bit in zip(coords, b):
            if bit: v[c] = 1
        V.append(v)
    return V

def closed_block(T, a, b, i0, j0):
    E = set()
    for i in range(a + 1):
        for j in range(b): E.add(T.h(i0 + i, j0 + j))
    for i in range(a):
        for j in range(b + 1): E.add(T.v(i0 + i, j0 + j))
    return E

def block_wraps(T, a, b):
    """the CLOSED edge set spans a+1 vertex rows and b+1 vertex columns; it fails to embed in a
       disk of T^2 exactly when it covers every row or every column."""
    return (a + 1 >= T.L) or (b + 1 >= T.L)

print("REGION TABLE.  'wraps' is decided by the TOPOLOGY of the closed edge set, not by size.")
for L in (3, 4, 5, 6):
    T = Toric(L); n = T.n
    pairs = symplectic_logicals(T.stab, n)
    flat = [v for pr in pairs for v in pr]
    M = [chi(T, v) for v in flat]
    from s1_toric_f2 import det_f2
    print("=" * 78)
    print(f"L = {L}   n = {n}   d = {L}")
    print(f"  chi on the computed logical basis = {M}   det_F2 = {det_f2(M)}  "
          f"-> chi : N(S)/S -> H_1 (+) H^1 is an ISOMORPHISM: {det_f2(M)==1}")
    rS = rank(T.stab, 2 * n)
    agree = True; disk_fail = 0; disk_n = 0; wrap_n = 0; wrap_carry = 0
    summary = {}
    for a in range(1, L):
        for b in range(1, L):
            for i0 in range(L):
                for j0 in range(L):
                    E = closed_block(T, a, b, i0, j0)
                    V = region_V(T, E)
                    dq = rank(T.stab + V, 2 * n) - rS
                    rc = rank([chi(T, v) for v in V], 4)
                    if dq != rc: agree = False
                    w = block_wraps(T, a, b)
                    if w:
                        wrap_n += 1; wrap_carry += (dq > 0)
                    else:
                        disk_n += 1; disk_fail += (dq > 0)
                    summary[(a, b)] = (len(E), dq, rc, w)
    for (a, b), (q, dq, rc, w) in sorted(summary.items()):
        print(f"    {a}x{b} plaquettes: {q:3d} qubits   dim V/S = {dq}   rank chi(V) = {rc}   "
              f"closed set wraps = {w}   {'CONTRACTIBLE -> clause (v) holds' if not w and dq==0 else ('NON-CONTRACTIBLE -> carries a logical' if w and dq>0 else 'MISMATCH')}")
    print(f"  dim(V/S) == rank chi(V) for every one of {disk_n+wrap_n} regions: {agree}")
    print(f"  CONTRACTIBLE regions tested: {disk_n}, of which carry a logical: {disk_fail}  "
          f"-> CLAUSE (v) HOLDS BY HOMOLOGY, NO CONVENTION: {disk_fail == 0}")
    print(f"  CONTROL, NON-contractible regions tested: {wrap_n}, of which carry a logical: {wrap_carry}  "
          f"-> the instrument registers a (v)-failure when one exists: {wrap_carry == wrap_n}")
