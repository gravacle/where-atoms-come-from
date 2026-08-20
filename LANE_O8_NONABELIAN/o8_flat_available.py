"""ROUTE 1: IS FLAT EVER UNAVAILABLE?

The claim needs a carrier where writers CANNOT commute. Before building a non-abelian model, ask
whether the five clauses permit that at all.

THE ARGUMENT TO TEST:
  clause (iii) puts the record inside a DEGENERATE eigenspace of H.
  admissibility is [U,H] = 0, which inside one eigenspace constrains NOTHING beyond preserving it.
  so a writer may be ANY unitary permuting the joint record blocks.
  blocks of EQUAL dimension can be identified with each other, and identifications chosen
  consistently give COMMUTING writers.
  and the blocks DO have equal dimension -- that is exactly clause (iv)'s trace balance (C-11).
  => flat is always available, and clause (iv) is what guarantees it.

If that holds, route 1 is closed from inside: curvature can never be FORCED by any carrier meeting
the clauses. Tested on carriers with degenerate joint record blocks, where the non-abelian freedom
is real and a naive argument would expect an obstruction."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
rng=np.random.default_rng(3)
say("="*104); say("ROUTE 1   IS A FLAT CONNECTION EVER UNAVAILABLE?"); say("="*104)

def build(dim, R):
    """joint record blocks, then a FLAT writer set built by consistent identification"""
    k=len(R)
    Mx=sum((2**i)*Rr for i,Rr in enumerate(R)); w,V=np.linalg.eigh(Mx)
    lab={}
    for c in range(dim):
        s=tuple(int(round(np.real(V[:,c].conj()@Rr@V[:,c]))) for Rr in R)
        lab.setdefault(s,[]).append(c)
    return V, lab

def flat_writers(dim, R, V, lab):
    """Identify every block with the all-plus block by a FIXED basis map, then define the writer
       for record i as 'change the i-th sign, keeping the position within the block'. Consistency
       is automatic, so these commute BY CONSTRUCTION -- if the blocks have equal size."""
    k=len(R); sizes={s:len(c) for s,c in lab.items()}
    if len(set(sizes.values()))!=1: return None, sizes
    Us=[]
    for i in range(k):
        U=np.zeros((dim,dim),dtype=complex)
        for s,cols in lab.items():
            t=list(s); t[i]=-t[i]; t=tuple(t)
            for a,b in zip(cols, lab[t]):          # position-preserving identification
                U += np.outer(V[:,b], V[:,a].conj())
        Us.append(U)
    return Us, sizes

I2=np.eye(2); Zm=np.array([[1,0],[0,-1]],dtype=complex)
def zq(i,nq):
    M=np.array([[1]],dtype=complex)
    for k in range(nq): M=np.kron(M, Zm if k==i else I2)
    return M
# records supplied DIRECTLY rather than derived, so the model's enumeration cap (O-28) does not
# bound the block size we can probe. Three records on nq qubits leaves blocks of dim 2^(nq-3).
CASES=[]
for nq,blab in ((3,"dim 1"),(4,"dim 2"),(5,"dim 4"),(6,"dim 8")):
    CASES.append((f"C^{2**nq}, H=0, blocks of {blab}", np.zeros((2**nq,2**nq),dtype=complex),
                  [zq(i,nq) for i in range(3)]))
CASES.append(("C^16, two 8-fold shells", np.diag(np.repeat([0.0,1.0],8)).astype(complex),
              [zq(i,4) for i in range(3)]))

say(f"  {'carrier':<32}{'records':>9}{'block sizes':>14}{'flat exists':>13}{'max ||[U_i,U_j]||':>20}")
for lbl,H,R in CASES:
    dim=H.shape[0]
    if not all(np.linalg.norm(Rr@H-H@Rr)<1e-9 for Rr in R):
        say(f"  {lbl:<32}   records do not commute with H -- skipped"); continue
    V,lab=build(dim,R)
    Us,sizes=flat_writers(dim,R,V,lab)
    bs=sorted(set(sizes.values()))
    if Us is None:
        say(f"  {lbl:<32}{len(R):>9}{str(bs):>14}{'NO':>13}{'blocks unequal':>20}"); continue
    # verify each is a writer, then measure the holonomy
    good=all(np.linalg.norm(Us[i].conj().T@Us[i]-np.eye(dim))<1e-9 and
             all(np.linalg.norm(Us[i].conj().T@R[j]@Us[i]-((-1 if j==i else 1)*R[j]))<1e-9
                 for j in range(len(R))) for i in range(len(R)))
    mc=max(float(np.linalg.norm(Us[i]@Us[j]-Us[j]@Us[i]))
           for i,j in itertools.combinations(range(len(R)),2))
    say(f"  {lbl:<32}{len(R):>9}{str(bs):>14}{('YES' if good else 'writers FAILED'):>13}{mc:>20.3e}")
say("")
say("  READ")
say("    Every carrier meeting the clauses has blocks of EQUAL dimension -- that is clause (iv)'s")
say("    trace balance (C-11) -- and equal blocks can be identified position-by-position, which makes")
say("    the writers commute BY CONSTRUCTION. Inside one degenerate eigenspace admissibility")
say("    constrains nothing further, so nothing can obstruct the identification.")
say("")
say("    If the right-hand column is zero everywhere, a FLAT connection is ALWAYS AVAILABLE and")
say("    curvature can never be FORCED by any carrier meeting the five clauses -- closing route 1")
say("    from inside the definition rather than by exhausting carriers.")
