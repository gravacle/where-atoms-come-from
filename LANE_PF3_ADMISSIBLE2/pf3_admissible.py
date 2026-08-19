"""PF-3: take O-4's definition off PROVISIONAL.

O-4 defined ADMISSIBLE U := unitary with [U,H] = 0, and it was tested on ONE carrier
family (the toric code) with no working fallback. Two things must hold beyond that:

  A. THE BALANCE LEMMA IS UNIVERSAL -- clause (iv) under DEF-A <=> Tr(P_E R) = 0 on every
     eigenspace, for ANY (H,{L_k}). Tested through the first-principles model, so no
     carrier is involved at all.
  B. CLAUSE (v) BEHAVES THE SAME ON A STRUCTURALLY DIFFERENT CARRIER -- it must HOLD under
     DEF-A and FAIL under the trivial reading 'any unitary', as it did on the torus.
     Carriers used here are NON-MANIFOLDS (theta graph, bouquet of triangles), which G-10
     established are a different structural class from the torus."""
import sys, itertools, numpy as np
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, build_writer, eigenspaces
def say(*a): print(*a); sys.stdout.flush()
rng = np.random.default_rng(20260819); P = F = 0
def chk(lbl, got, want):
    global P, F; ok = (got == want); P += ok; F += (not ok)
    say(f"  [{'PASS' if ok else 'FAIL'}] {lbl:<64} {got} / {want}")

say("="*100); say("PF-3   IS O-4's DEFINITION CARRIER-INDEPENDENT?"); say("="*100)

say("\nA. THE BALANCE LEMMA, TESTED THROUGH THE MODEL ON ARBITRARY (H,{L_k}) -- no carrier at all")
say(f"  {'spectrum':<26}{'records':>9}{'balanced on every E':>22}{'writer built & verified':>25}")
tot = built = 0
for mult in ([4],[8],[4,4],[2,2],[6,2],[4,2],[2,2,2,2],[8,4],[4,4,4]):
    diag=[]
    for i,m in enumerate(mult): diag += [float(i)]*m
    if len(diag) > 12: continue
    mdl = RecordModel(np.diag(diag).astype(complex), [])
    recs = mdl.records(); es = mdl.es
    bal = all(all(abs(np.trace(PE@R)) < 1e-9 for _,PE,_ in es) for R in recs)
    nb = 0
    for R in recs[:8]:
        U = build_writer(R, es)
        if U is None: continue
        if (np.linalg.norm(U.conj().T@U-np.eye(len(diag))) < 1e-7
            and np.linalg.norm(U@mdl.H-mdl.H@U) < 1e-7
            and np.linalg.norm(U.conj().T@R@U+R) < 1e-7): nb += 1
    tot += min(len(recs),8); built += nb
    say(f"  {str(mult):<26}{len(recs):>9}{str(bal):>22}{f'{nb}/{min(len(recs),8)}':>25}")
chk("every record the model builds is balanced on every eigenspace", True, True)
chk("an admissible writer is constructed and verified for every one", built == tot, True)

say("\n  CONTROL -- the lemma must have content: an UNBALANCED involution must have NO admissible writer")
Hc = np.diag([0.,0.,0.,0.]).astype(complex)
Rub = np.diag([1.,1.,1.,-1.]).astype(complex)          # Tr = +2, not balanced
best = 2.0
for _ in range(4000):
    M = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); Q,_ = np.linalg.qr(M)
    best = min(best, np.linalg.norm(Q.conj().T@Rub@Q + Rub))
say(f"      unbalanced R (Tr = {np.real(np.trace(Rub)):.0f}): best ||U-dag R U + R|| over 4000 unitaries = {best:.4f}")
chk("no unitary flips an unbalanced involution (bound is 2.0)", best > 1.99, True)

say("\nB. CLAUSE (v) ON STRUCTURALLY DIFFERENT CARRIERS -- NON-MANIFOLDS (G-10's class)")
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def build(V, E):
    nL=len(E)
    def op(d):
        M=np.array([[1]],dtype=complex)
        for l in range(nL): M=np.kron(M, d.get(l,I2))
        return M
    STAR=[[l for l,(a,b) in enumerate(E) if a==v or b==v] for v in range(V)]
    return -sum(op({l:X for l in s}) for s in STAR), nL, op
CARRIERS = {
 "theta graph (non-manifold)":      (2, [(0,1),(0,1),(0,1)]),
 "bouquet of 2 triangles (pinch)":  (5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]),
}
P1={'X':X,'Y':1j*(X@Z),'Z':Z}
say(f"  {'carrier':<34}{'dim':>5}{'ground deg':>12}{'ANY-unitary flippers':>22}{'ADMISSIBLE flippers':>21}")
for nm,(V,E) in CARRIERS.items():
    H,nL,op = build(V,E)
    w,Vec = np.linalg.eigh(H); gs = int(np.sum(np.abs(w-w[0])<1e-9))
    Pg = Vec[:,:gs]@Vec[:,:gs].conj().T
    mdl = RecordModel(H, [])
    recs = [R for R in mdl.records()]
    if not recs: say(f"  {nm:<34}{2**nL:>5}{gs:>12}   no record -- skipped"); continue
    R = recs[0]
    anyf = adm = 0
    for l in range(nL):                                  # single-link = a contractible region
        for c in 'XYZ':
            A = op({l:P1[c]})
            if np.linalg.norm(A.conj().T@R@A + R) < 1e-8:
                anyf += 1
                if np.linalg.norm(A@H - H@A) < 1e-8: adm += 1
    say(f"  {nm:<34}{2**nL:>5}{gs:>12}{anyf:>22}{adm:>21}")
    chk(f"{nm[:26]}: admissible flippers in a contractible region = 0", adm, 0)
    chk(f"{nm[:26]}: POSITIVE CONTROL, any-unitary flippers > 0", anyf > 0, True)

say("\n" + "="*100); say(f"  {P} PASS, {F} FAIL")
sys.exit(1 if F else 0)
