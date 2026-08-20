"""AUDIT 2: ARE THE S_3 'RECORDS' RECORDS BY THIS PROGRAM'S OWN DEFINITION?

O-29 and O-30 called the conjugacy classes and gauge classes of S_3 'records' and measured transport
between them. That identification was IMPORTED from lattice gauge theory and never checked against
the five clauses. A class is a LABEL, not an operator -- it does not obviously satisfy R = R^dagger,
R^2 = I, or any of the rest.

So build the S_3 quantum double as an actual Hamiltonian and put it through model.records(), which
CONSTRUCTS every R satisfying (i)-(iv) rather than searching for one. If the answer is zero, then
O-29/O-30 measured transport between things that are not records by this program's definition, and
they must be rescoped.

CONTROL (D-15): the same pipeline on the toric code, where records are known to exist."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
say("="*104); say("AUDIT 2   ARE THE S_3 'RECORDS' RECORDS BY THE FIVE CLAUSES?"); say("="*104)
G=list(itertools.permutations(range(3)))
def mul(a,b): return tuple(a[b[i]] for i in range(3))
def inv(a):
    r=[0]*3
    for i,x in enumerate(a): r[x]=i
    return tuple(r)
e=(0,1,2); idx={g:i for i,g in enumerate(G)}; n=len(G)
D=n*n
def ket(g1,g2): return idx[g1]*n+idx[g2]
say(f"  D(S_3) on the minimal torus: 1 vertex, 2 edges, 1 face.  Hilbert space |G|^2 = {D}")
# gauge projector A = (1/|G|) sum_h A_h,  A_h |g1,g2> = |h g1 h^-1, h g2 h^-1>
A=np.zeros((D,D))
for h in G:
    for g1 in G:
        for g2 in G:
            A[ket(mul(mul(h,g1),inv(h)),mul(mul(h,g2),inv(h))), ket(g1,g2)] += 1.0/n
# flatness projector B: diagonal, 1 where the plaquette holonomy is trivial
B=np.zeros((D,D))
for g1 in G:
    for g2 in G:
        if mul(mul(g1,g2),mul(inv(g1),inv(g2)))==e: B[ket(g1,g2),ket(g1,g2)]=1.0
for nm,P in (("A (gauge)",A),("B (flat)",B)):
    say(f"    {nm:<12} projector? ||P^2-P|| = {np.linalg.norm(P@P-P):.2e}   rank {int(round(np.trace(P)))}")
H=-(A+B)
M=RecordModel(H.astype(complex))
Pg,kdim=M.ground_space()
say(f"  H = -(A+B).  ground space dimension {kdim}")
say(f"  distinct eigenvalues of H: {[round(float(v),6) for v,_,_ in M.es]}   multiplicities {[m for _,_,m in M.es]}")
say(f"  minimal projections in the commutant: {len(M.projs)}")
try:
    recs=M.records()
    say(f"  model.records() -> {len(recs)} operators satisfying clauses (i)-(iv)")
    if recs:
        R=recs[0]
        say(f"    check on the first: ||R-R^dag||={np.linalg.norm(R-R.conj().T):.2e}  "
            f"||R^2-I||={np.linalg.norm(R@R-np.eye(D)):.2e}  ||[H,R]||={np.linalg.norm(H@R-R@H):.2e}")
except RuntimeError as ex:
    say(f"  model.records() REFUSED: {ex}")
    recs=None
say("")
say("  CONTROL -- the same pipeline on a carrier where records are known to exist")
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def pauli(word):
    Mm=np.array([[1]],dtype=complex)
    for c in word: Mm=np.kron(Mm,{'I':I2,'X':X,'Z':Z}[c])
    return Mm
Hc=-(pauli('XXXX')+pauli('ZZZZ'))
Mc=RecordModel(Hc)
say(f"    [[4,2,2]] dim 16, ground space {Mc.ground_space()[1]}, minimal projections {len(Mc.projs)}")
try:
    rc=Mc.records(); say(f"    model.records() -> {len(rc)} operators satisfying (i)-(iv)")
except RuntimeError as ex:
    say(f"    model.records() REFUSED: {ex}")
say("")
say("  ENUMERATION IS BLOCKED BY O-28, BUT CLAUSE (iv) CAN BE DECIDED WITHOUT IT.")
say("  The framework states: given (i)+(ii), clause (iv) <=> Tr(P_E R) = 0 on EVERY eigenspace of H.")
say("  And (i) forces R^2 = I, so every eigenvalue of R is +-1 and Tr(P_E R) is a sum of dim(E)")
say("  terms each +-1 -- hence Tr(P_E R) = dim(E) (mod 2). AN EIGENSPACE OF ODD DIMENSION CAN NEVER")
say("  GIVE ZERO.")
say("")
say(f"  {'carrier':<22}{'eigenspace multiplicities':>32}{'any ODD?':>11}{'clause (iv) satisfiable?':>26}")
for nm,MM in (("D(S_3) minimal torus",M),("[[4,2,2]] (control)",Mc)):
    mult=[m for _,_,m in MM.es]
    odd=[m for m in mult if m%2==1]
    say(f"  {nm:<22}{str(mult):>32}{str(bool(odd)):>11}  {('NO -- (iv) is UNSATISFIABLE' if odd else 'yes')}")
say("")
say("  DIRECT CHECK, not argument. Clause (iv) needs Tr(P_E R) = 0 on EVERY eigenspace, so the")
say("  binding case is the WORST one, not the easiest. For each eigenspace report the smallest")
say("  |Tr(P_E R)| any sign pattern can achieve, then take the MAXIMUM over eigenspaces.")
say(f"  {'carrier':<22}{'per-eigenspace best |Tr|':>28}{'WORST (binding)':>18}{'(iv) possible?':>17}")
for nm,MM in (("D(S_3) minimal torus",M),("[[4,2,2]] (control)",Mc)):
    per=[min(abs(m-2*j) for j in range(m+1)) for _,_,m in MM.es]
    say(f"  {nm:<22}{str(per):>28}{max(per):>18}{str(max(per)==0):>17}")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  THE GROUND SPACE IS DIMENSION 8, matching the 8 gauge classes counted group-theoretically in")
say("  O-29. That identification cross-checks.")
say("")
say("  BUT NO RECORD EXISTS ON THIS CARRIER. H = -(A+B) on the minimal torus has eigenspaces of")
say("  dimension 8, 13 and 15. Clauses (i)+(ii) force Tr(P_E R) = dim(E) mod 2, so on a 13- or")
say("  15-dimensional eigenspace it can never vanish, and CLAUSE (iv) -- WRITABILITY -- IS")
say("  UNSATISFIABLE. The control has multiplicities that are all even and yields 1260 records.")
say("")
say("  SO O-29 AND O-30 MEASURED TRANSPORT BETWEEN OBJECTS THAT ARE NOT RECORDS BY THIS PROGRAM'S")
say("  DEFINITION ON THAT CARRIER. The transport itself is real group theory and the numbers stand,")
say("  but the word 'record' was imported, not earned. Those rows must be rescoped.")
say("")
say("  IT IS A CARRIER DEFECT, NOT A VERDICT ON NON-ABELIAN CARRIERS. The minimal torus has ONE")
say("  vertex and ONE face, so the gauge and flatness projectors have nowhere to act -- H-4's")
say("  small-carrier warning again. A larger S_3 complex is the next thing to build, and whether its")
say("  multiplicities are even is a question about that carrier, not about S_3.")
