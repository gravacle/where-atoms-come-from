"""Does the first-principles model REPRODUCE the register's results, given only (H,{L_k})?

Every case below is a registered result. The model is told NOTHING about lattices, gauge
groups, codes or geometry -- only a Hamiltonian and a set of Lindblad operators."""
import sys, numpy as np, itertools
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, build_writer, eigenspaces
def say(*a): print(*a); sys.stdout.flush()
rng = np.random.default_rng(20260818)
def gen(n):
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return M
P = 0; F = 0
def check(name, got, want):
    global P, F
    ok = (got == want); P += ok; F += (not ok)
    say(f"  [{'PASS' if ok else 'FAIL'}] {name:<62} got {got}, expected {want}")

say("="*96); say("MODEL VALIDATION -- first principles only, input is (H,{L_k}) and nothing else"); say("="*96)

say("\nV1  maximal degeneracy, no noise:  H = 0 on C^4")
check("records exist", RecordModel(np.zeros((4,4)), []).report()['n_records'] > 0, True)

say("\nV2  C-11 corollary: ODD dimension => Tr R is odd => NO record is ever writable")
for n in (3, 5, 7):
    check(f"H = 0 on C^{n}: n_records", RecordModel(np.zeros((n,n)), []).report()['n_records'], 0)

say("\nV3  O-1 gap (a): a degenerate H plus ONE GENERIC jump kills every record")
H4 = np.diag([0.,0.,1.,1.]).astype(complex)
cnt = sum(RecordModel(H4, [gen(4)]).report()['n_records'] for _ in range(6))
check("H=diag(0,0,1,1) + generic L: total records over 6 draws", cnt, 0)

say("\nV4  POSITIVE CONTROL for V3: the same H with DIAGONAL jumps keeps records")
cnt = sum(RecordModel(H4, [np.diag(rng.normal(size=4)).astype(complex)]).report()['n_records'] for _ in range(6))
check("H=diag(0,0,1,1) + diagonal L: records found in all 6 draws", cnt >= 6, True)

say("\nV5  O-1's TRAP: noise identically ZERO inside every eigenspace, yet NO record")
L1 = np.zeros((4,4), dtype=complex); L1[2,0] = 1; L1[1,2] = 1
L2 = np.zeros((4,4), dtype=complex); L2[3,1] = 1; L2[0,3] = 1
es = eigenspaces(H4)
inside = max(np.linalg.norm(Pp @ L @ Pp) for _, Pp, _ in es for L in (L1, L2))
say(f"       max ||P_E L P_E|| = {inside:.3e}   (the noise vanishes inside every eigenspace)")
check("with BOTH legs (excursion out and back): n_records", RecordModel(H4, [L1, L2]).report()['n_records'], 0)
L1b = np.zeros((4,4), dtype=complex); L1b[2,0] = 1
L2b = np.zeros((4,4), dtype=complex); L2b[3,1] = 1
check("CONTROL, return legs deleted: records exist", RecordModel(H4, [L1b, L2b]).report()['n_records'] > 0, True)

say("\nV6  GAUGE CARRIER FROM FIRST PRINCIPLES -- Z2 on the theta graph, built and handed in as H")
V, E = 2, [(0,1),(0,1),(0,1)]; nL = len(E)
I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex); Zp = np.array([[1,0],[0,-1]], dtype=complex)
def op(d):
    M = np.array([[1]], dtype=complex)
    for l in range(nL): M = np.kron(M, d.get(l, I2))
    return M
STAR = [[l for l,(a,b) in enumerate(E) if a==v or b==v] for v in range(V)]
Hg = -sum(op({l: X for l in s}) for s in STAR)
m = RecordModel(Hg, []); rep = m.report()
say(f"       dim {rep['dim']}  algebra {rep['dim_algebra']}  commutant {rep['dim_commutant']}  "
    f"minimal projections {rep['n_minimal_projections']}  RECORDS {rep['n_records']}")
check("theta graph carries records (dim H_1 = 2 => it must)", rep['n_records'] > 0, True)

say("\nV7  MULTI-RECORD STRUCTURE -- do records commute, and are they independently writable?")
for n in (4, 8, 16):
    m = RecordModel(np.zeros((n,n)), []); recs = m.records()
    fam, comm, indep = m.independence(recs)
    say(f"       C^{n}: {len(recs)} records -> independent family {len(fam)}; all commute "
        f"{bool(comm.all())}; independently writable {len(indep)}")
    check(f"C^{n}: family size = log2(dim), every member independently writable",
          (len(fam), len(indep)) == (int(np.log2(n)), int(np.log2(n))), True)

say("\n" + "="*96)
say(f"  {P} PASS, {F} FAIL")
say("  Every case above was decided by the model from (H,{L_k}) alone -- no lattice, no gauge")
say("  group, no temperature, no coupling, no code, no geometry was supplied.")
sys.exit(1 if F else 0)
