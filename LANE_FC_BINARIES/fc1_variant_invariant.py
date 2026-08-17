# FRAME CHALLENGE on the corpus's two inherited binaries.
# T2 (embed) clause 1 is EXCLUSION: exhibit a point lying in NEITHER named arm, with a count.
# Binary under test: VARIANT / INVARIANT. Carrier: the OPEN graph from AA3 -- bulk vertices carry
# the Gauss constraint, one dangling link whose far end is BOUNDARY and carries none.
import numpy as np, itertools
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
OPEN=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4),(0,6)]   # link 9 dangles at vertex 6
n=len(OPEN); BULK=list(range(6)); BOUNDARY=[6]
Z=[op(i,Zp,n) for i in range(n)]; X=[op(i,Xp,n) for i in range(n)]
def gauss(v):
    M=np.eye(2**n,dtype=complex)
    for i,(a,b) in enumerate(OPEN):
        if a==v or b==v: M=M@Z[i]
    return M
Gb=[gauss(v) for v in BULK]
print(f"  OPEN carrier: {len(BULK)} bulk vertices, {len(BOUNDARY)} boundary, L={n}")

# the boundary charge: Gauss operator at the UNCONSTRAINED vertex
Q=gauss(6)
print(f"\n== IS THE BOUNDARY CHARGE CENTRAL? (commutes with every bulk constraint) ==")
print(f"   max || [Q, G_v] || over bulk v  =  {max(np.linalg.norm(Q@g-g@Q) for g in Gb):.3e}")
print(f"   Q^2 = I ?  {np.allclose(Q@Q,np.eye(2**n))}   -> two sectors, Q = +1 and Q = -1")

# classify every single-link operator three ways against the BULK gauge group
def cls(A):
    inv = max(np.linalg.norm(A@g-g@A) for g in Gb) < 1e-9
    shifts = np.linalg.norm(A@Q+Q@A) < 1e-9        # anticommutes with Q -> maps sector to sector
    return inv, shifts
print(f"\n== THE THREE-WAY CLASSIFICATION, ON EVERY SINGLE-LINK OPERATOR ==")
print(f"   {'operator':<12}{'commutes w/ bulk Gauss':>24}{'changes the sector':>21}{'category':>26}")
rows={}
for i in range(n):
    for nm,P in (("Z",Zp),("X",Xp)):
        A=op(i,P,n); inv,sh=cls(A)
        cat = "INVARIANT (bulk obs)" if inv and not sh else \
              "SECTOR-CHANGING" if sh else "VARIANT (neither)"
        rows.setdefault(cat,[]).append(f"{nm}_{i}")
        if i in (0,9): print(f"   {nm+'_'+str(i):<12}{str(inv):>24}{str(sh):>21}{cat:>26}")
print()
for k,v in rows.items(): print(f"   {k:<26} {len(v):>3} operators   e.g. {v[:4]}")
print()
print("   >>> EXCLUSION CLAUSE: the SECTOR-CHANGING operators lie in NEITHER named arm.")
print("       They are not bulk-invariant, so the binary calls them 'variant' and discards them.")
print("       But they are exactly the operators that move between superselection sectors --")
print("       i.e. the ones that would CREATE or DESTROY boundary charge. They are not junk.")
print("       VARIANT/INVARIANT IS A TWO-POINT SAMPLE OF A THREE-WAY STRUCTURE.")
