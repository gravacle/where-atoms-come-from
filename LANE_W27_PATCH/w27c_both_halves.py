# W-27c — BOTH HALVES OF EM ON THE SAME BOUNDARY, FOR THE FIRST TIME.
# The patch perimeter is BOTH a cycle (degrees all 2) and a separator (removal -> 5 components).
# So the same set of links carries an ELECTRIC flux (product of Z, the cut quantity) and a
# MAGNETIC flux (product of X, the loop quantity). Do they commute? Non-commutation on the boundary
# is an uncertainty relation on what a boundary can hold.
import numpy as np, itertools
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E); NV=len(V2); N=2; CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
PERIM=[k for k in range(L) if k not in CUT]
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(NV))]
idx={s:j for j,s in enumerate(st)}; D=len(st)
def loop(moves):
    M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        k=idx.get(tuple(t))
        if k is None: return None
        M[k,j]=1.0
    return M
Zset=lambda S: np.diag([(-1.0)**sum(s[k] for k in S) for s in st]).astype(complex)
print(f"  patch 3x3, Z_2, physical dim {D}")
print(f"  PERIMETER = {PERIM}  ({len(PERIM)} links) -- a cycle AND a separator")
print(f"  CUT at the interior vertex = {CUT}  ({len(CUT)} links) -- a separator, not a cycle\n")
Eper=Zset(PERIM); Mper=loop([(k,+1) for k in PERIM])
Ecut=Zset(CUT)
print("  BOTH HALVES OF EM ON THE PERIMETER:")
print(f"    electric (prod Z on perimeter)  nontrivial: {np.linalg.norm(Eper-np.eye(D))>1e-9}")
print(f"    magnetic (prod X on perimeter)  nontrivial: {Mper is not None and np.linalg.norm(Mper-np.eye(D))>1e-9}")
print(f"    || [E_perim, M_perim] || = {np.linalg.norm(Eper@Mper-Mper@Eper):.3e}   (8 links, even)")
print()
print("  AND ON A HALF-PERIMETER, WHERE THE COUNT IS ODD:")
for half in ([PERIM[0]],PERIM[:3],PERIM[:4],PERIM[:5]):
    Eh=Zset(half); Mh=loop([(k,+1) for k in half])
    c = np.linalg.norm(Eh@Mh-Mh@Eh) if Mh is not None else float('nan')
    print(f"    links {str(half):<22} |S|={len(half)}  || [E,M] || = "
          + (f"{c:.3e}" if Mh is not None else "M leaves the physical sector"))
print()
print("  ANY link subset whose X-product stays physical is a closed loop, and a closed loop on this")
print("  lattice touches an EVEN number of links -- so E and M on the SAME set always commute.")
print("  The uncertainty is between a loop and a CUT THAT CROSSES IT, not between the two halves")
print("  of the field on one boundary.")
print(f"    || [E_cut, M_perim] || = {np.linalg.norm(Ecut@Mper-Mper@Ecut):.3e}  (disjoint sets)")
for k in range(1,5):
    S=PERIM[:k]
    print(f"    || [Z on {str(S):<18}, M_perim] || = {np.linalg.norm(Zset(S)@Mper-Mper@Zset(S)):.3e}")
