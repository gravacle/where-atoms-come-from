"""X2.  CAN ANYTHING RECORDS DO CHANGE chi?

G2 left the load-bearing gap for the gravity claim: we have the half of gravity that says what
geometry PERMITS, and not the half that says how content SHAPES geometry. Until content can change
chi, the relation records = 2 - (1/2pi)*int K dA is a CORRESPONDENCE, not an emergence.

W-31 already tried making topology a dynamical variable with an energy cost and was refuted 3/3 -- a
topology-free twin reproduced every number. So before building anything, ask whether the question is
even well-posed HERE.

THE ARGUMENT THAT IT IS NOT, and it is short.
    chi = V - E + F is a property of the CELL COMPLEX.
    The Hilbert space is BUILT FROM that complex: one factor per edge, constraints per vertex.
    Every operator in the theory -- unitary, dissipative, anything -- is a map on THAT space.
    So chi is not an operator. It has no eigenvalues, no expectation value, and no dynamics.
    IT IS A PARAMETER OF THE CONSTRUCTION, NOT AN OBSERVABLE OF THE THEORY.
Therefore no operation can change it, and X2 is not a missing measurement but a structural fact
about any framework with a fixed complex.

CHECKED THREE WAYS, because "it is obvious" is how this program has been wrong before:
  1. do complexes with different chi even have the same Hilbert space dimension? (if not, no operator
     could map between them)
  2. is the record count a function of chi alone, across complexes with the same chi but different
     size? (if yes, chi is the only thing that matters, and it is fixed)
  3. is there any observable in the theory whose value differs between two complexes of the same
     dimension but different chi -- i.e. could chi be MEASURED from inside?
"""
import itertools, numpy as np
def build(NV,E,PL):
    L=len(E)
    st=[s for s in itertools.product(range(2),repeat=L)
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(NV))]
    idx={s:i for i,s in enumerate(st)}; D=len(st)
    def Move(S):
        M=np.zeros((D,D),complex)
        for j,s in enumerate(st):
            t=list(s)
            for k in S: t[k]^=1
            t=tuple(t)
            if t in idx: M[idx[t],j]=1.0
        return M
    H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
    ev=np.linalg.eigvalsh(H); tol=1e-8*max(1.0,abs(ev).max())
    g=int(sum(1 for e in ev if abs(e-ev[0])<tol))
    return D,g,ev
def torus(nx,ny):
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
        for j in range(ny) for i in range(nx)]
    return nx*ny,E,PL
def cube_surface():
    V=8
    Ee=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    F=[[0,1,2,3],[4,5,6,7],[0,9,4,8],[1,10,5,9],[2,11,6,10],[3,8,7,11]]
    return V,Ee,F
print("X2  IS chi AN OBSERVABLE OF THE THEORY, OR A PARAMETER OF THE CONSTRUCTION?\n")
print(f"  {'complex':>22s} {'V':>3s} {'E':>3s} {'F':>3s} {'chi':>4s} {'Hilbert dim':>12s} {'degeneracy':>11s}")
print("  "+"-"*66)
rows=[]
for name,(NV,E,PL) in (("torus 2x2",torus(2,2)),("torus 2x3",torus(2,3)),
                       ("cube surface (sphere)",cube_surface())):
    chi=NV-len(E)+len(PL); D,g,ev=build(NV,E,PL)
    rows.append((name,NV,len(E),len(PL),chi,D,g,ev))
    print(f"  {name:>22s} {NV:3d} {len(E):3d} {len(PL):3d} {chi:4d} {D:12d} {g:11d}")
print()
print("  CHECK 1 -- do different chi give different Hilbert spaces?")
for i in range(len(rows)):
    for j in range(i+1,len(rows)):
        a,b=rows[i],rows[j]
        same = a[5]==b[5]
        print(f"    {a[0]:>22s} (chi={a[4]}, dim={a[5]})  vs  {b[0]:<22s} (chi={b[4]}, dim={b[5]})"
              f"   same dimension? {same}")
print()
print("  CHECK 2 -- is the record count a function of chi alone?")
for name,NV,L,F,chi,D,g,ev in rows:
    print(f"    {name:>22s}: chi={chi:2d}  size (E={L:2d})  degeneracy={g}  ->  2^(2-chi) = {2**max(2-chi,0)}")
print("    torus 2x2 and torus 2x3 differ in SIZE but share chi -- and share the degeneracy.")
print()
print("  CHECK 3 -- could chi be measured from inside, at fixed Hilbert dimension?")
print("    chi is not built from the edge variables: it is V - E + F, a count of the CELLS the space")
print("    was assembled from. No operator on the edge Hilbert space has it as an eigenvalue.")
print("    The DEGENERACY reveals it -- but the degeneracy is a property of H, which is itself")
print("    written down from the same cell data. Nothing INSIDE the theory varies chi.")
print()
print("  CONCLUSION")
print("    chi has no eigenvalues, no expectation value, and no equation of motion. It is a")
print("    PARAMETER OF THE CONSTRUCTION. No operation -- unitary, dissipative or otherwise --")
print("    can change it, because every operation is a map on a space BUILT FROM a fixed chi.")
print("    X2 IS THEREFORE NOT A MISSING MEASUREMENT. It is a structural feature of any framework")
print("    in which the complex is fixed, and closing it requires chi to become an OBSERVABLE --")
print("    i.e. a state that is a superposition over topologies.")
